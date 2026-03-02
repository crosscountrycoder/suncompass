import { degToRad, earthERadius, earthPRadius } from "./constants.ts";
import {subsolarPoint} from "./suncalc.ts";
import * as mf from "./mathfuncs.ts";
import * as svg from "./svg.ts";
import { generateLODProfile } from "./lookup-tables.ts";
import { moonDistance, sublunarPoint } from "./mooncalc.ts";
import {baseWorldMap, moonIcon, sunIcon} from "./svg-lib.ts";
import { DateTime } from "luxon";

const NUM_POINTS = 180; 

/** Given a lat/long coordinate, angular radius, and number of points, returns a closed regular polygon (approximation to 
 * a circle) with numPoints points centered at [lat, long]. 
 * @param lat Center latitude in degrees
 * @param long Center longitude in degrees
 * @param radius Angular radius in degrees
 * @param numPoints Number of points (optional, defaults to 360 if not specified)
*/
function regularPolygon(lat: number, long: number, radius: number, numPoints: number = NUM_POINTS): mf.PolygonBearing {
    const points: mf.PolygonBearing = [];
    for (let i = 0; i <= numPoints; i++) {
        const bearing = (i * 360) / numPoints;
        points.push(mf.destPoint(lat, long, bearing, radius));
    }
    return points;
}

/** Finds the longitude halfway between two longitudes, restricting to the range [-180, 180) */
function middleLong(long1: number, long2: number): number {
    if (Math.abs(long1 - long2) < 180) {return (long1 + long2)/2;}
    else {return mf.mod((long1 + long2 + 360)/2 + 180, 360) - 180;}
}

export function regularPolygonErrors(lat: number, long: number, radius: number, numPoints: number = NUM_POINTS, refined = false) {
    const polygon = refined ? refinedRegularPolygon(lat, long, radius, numPoints) :
        regularPolygon(lat, long, radius, numPoints);
    console.log(`Total number of points: ${polygon.length}`);
    for (let i = 0; i < polygon.length - 1; i++) {
        const midpointLat = (polygon[i][0] + polygon[i+1][0]) / 2;
        const midpointLong = middleLong(polygon[i][1], polygon[i+1][1]);
        const midpointBearing = (polygon[i][2] + polygon[i+1][2]) / 2;
        const angularDist = mf.distBearing(lat, long, midpointLat, midpointLong)[0];
        const error = Math.abs(angularDist - radius);
        console.log(`${mf.toFixedS(midpointBearing, 4)}: ${mf.toFixedS(error, 4)}`);
    }
}

/** Refines the regular polygon until it is a circle with accuracy 0.01 degrees (i.e. all points are within 0.01 degrees of
 * the actual circle.) */
function refinedRegularPolygon(lat: number, long: number, radius: number, numPoints: number = NUM_POINTS): mf.PolygonBearing {
    const tolerance = 0.01;
    const polygon = regularPolygon(lat, long, radius, numPoints);
    const newPolygon: mf.PolygonBearing = [];
    for (let i = 0; i < polygon.length - 1; i++) {
        const midpointLat = (polygon[i][0] + polygon[i+1][0])/2;
        const midpointLong = middleLong(polygon[i][1], polygon[i+1][1]);
        const error = Math.abs(mf.distBearing(lat, long, midpointLat, midpointLong)[0] - radius);
        if (error < tolerance) {
            if (i === polygon.length - 2) {newPolygon.push(polygon[i], polygon[i+1]);}
            else {newPolygon.push(polygon[i]);}
        }
        else {
            const newPointSet = [polygon[i], polygon[i+1]];
            let numIterations = 0;
            while (numIterations < 20) {
                let numErrors = 0;
                for (let j = newPointSet.length - 1; j >= 1; j--) {
                    const mLat = (newPointSet[j-1][0] + newPointSet[j][0])/2;
                    const mLong = middleLong(newPointSet[j-1][1], newPointSet[j][1]);
                    const mBearing = (newPointSet[j-1][2] + newPointSet[j][2])/2;
                    const err = Math.abs(mf.distBearing(lat, long, mLat, mLong)[0] - radius);
                    if (err >= tolerance) {
                        numErrors++;
                        const newPoint = mf.destPoint(lat, long, mBearing, radius);
                        newPointSet.splice(j, 0, newPoint);
                    }
                }
                if (numErrors === 0) {break;}
                numIterations++;
            }
            if (i === polygon.length - 2) {newPolygon.push(...newPointSet);}
            else {newPolygon.push(...newPointSet.slice(0, -1));}
        }
    }
    return newPolygon;
}

/** Given a latitude/longitude polygon, splits it into two parts if the polygon crosses the 180th meridian.
 * @returns An array of polygons, which contains one polygon if it doesn't cross the 180th meridian, or two polygons if it does.
 */
function splitPolygon(polygon: mf.Polygon): mf.Polygon[] {
    const newPolygon: mf.Polygon[] = [[], [], []]; // [west, east]
    newPolygon[0].push(polygon[0]);
    let curPolygon = 0;
    for (let i=1; i<polygon.length; i++) {
        const prev = polygon[i-1], cur = polygon[i];
        const eastToWest = prev[1] > 90 && cur[1] < -90;
        const westToEast = prev[1] < -90 && cur[1] > 90;
        if (eastToWest) {
            const frac = (180 - prev[1]) / (cur[1] + 360 - prev[1]);
            const latCrossing = prev[0] + frac * (cur[0] - prev[0]);
            newPolygon[curPolygon].push([latCrossing, 180]);
            curPolygon++;
            newPolygon[curPolygon].push([latCrossing, -180]);
        } else if (westToEast) {
            const frac = (prev[1] + 180) / (prev[1] + 360 - cur[1]);
            const latCrossing = prev[0] + frac * (cur[0] - prev[0]);
            newPolygon[curPolygon].push([latCrossing, -180]);
            curPolygon++;
            newPolygon[curPolygon].push([latCrossing, 180]);
        }
        if (i !== polygon.length) {newPolygon[curPolygon].push(cur);} // bearing 360° = bearing 0°
    }
    if (newPolygon[1].length === 0) {return [newPolygon[0]];}
    else {return [[...newPolygon[2], ...newPolygon[0]], newPolygon[1]];}
}

/** For polygons that include either the north or the south pole (but not both), this function sorts the points from
 * west to east, with additional points for longitudes -180 and +180. This is because when a polygon includes one pole, its
 * boundary crosses all lines of longitude and the area within the polygon is simply anything north or south of this boundary.
 */
function sortPolygon(polygon: mf.Polygon): mf.Polygon {
    const newPolygon = polygon.slice(0, -1).sort((a, b) => a[1] - b[1]);
    const [firstLat, firstLong] = newPolygon[0];
    const [lastLat, lastLong] = newPolygon.at(-1)!;
    const frac = (180 - lastLong) / (firstLong + 360 - lastLong);
    const latCrossing = lastLat + frac * (firstLat - lastLat);
    return [[latCrossing, -180], ...newPolygon, [latCrossing, 180]];
}

/** For a polygon stored in [lat, long] order, reverse the points to the SVG compatible [long, lat]. */
function reverseLatLong(polygon: mf.Polygon): mf.Polygon {
    return polygon.map(([lat, lon]) => [lon, lat] as [number, number]);
}

/**
 * Returns the SVG encoding of a polygon approximating a spherical cap with the given angular radius centered at the given
 * latitude and longitude.
 * @param lat Latitude at center
 * @param long Longitude at center
 * @param radius Angular radius
 * @param numPoints Minimum number of points (defaults to 180)
 * @param r RGBA fill color (red) - default is 0
 * @param g RGBA fill color (green) - default is 0
 * @param b RGBA fill color (blue) - default is 0
 * @param a RGBA fill color (alpha/opacity) - default is 0.3
 * @returns The SVG encoding of 
 */
function svgCircleEquirectangular(lat: number, long: number, radius: number, numPoints: number = NUM_POINTS,
    r: number = 0, g: number = 0, b: number = 0, a: number = 0.3): string {
    if (Math.abs(Math.abs(lat) + radius - 90) < 1e-5) {radius -= 2e-5;}
    const containsNorthPole = lat + radius > 90;
    const containsSouthPole = lat - radius < -90;
    const polygon: mf.Polygon = refinedRegularPolygon(lat, long, radius, numPoints).map(([x, y]) => [x, y]);
    const color = `rgba(${r}, ${g}, ${b}, ${a})`;
    if (containsNorthPole) {
        const sortedPolygon: mf.Polygon = [[90, 180], [90, -180], ...sortPolygon(polygon)];
        return svg.pathFromArray(reverseLatLong(sortedPolygon), true, color);
    } else if (containsSouthPole) {
        const sortedPolygon: mf.Polygon = [[-90, 180], [-90, -180], ...sortPolygon(polygon)];
        return svg.pathFromArray(reverseLatLong(sortedPolygon), true, color);
    } else {
        let s: string = "";
        const splitPolygons = splitPolygon(polygon);
        for (const p of splitPolygons) {s += svg.pathFromArray(reverseLatLong(p), true, color);}
        return s;
    }
}

function angularRadiusVisibility(r: number, dist: number): number {
    const earthRadius = (earthERadius + earthPRadius) / 2;
    const eMinRad = (90 - r) * degToRad;
    
    function elevationRad(psiRad: number): number {
        const c = Math.cos(psiRad);
        const denom = Math.sqrt(dist ** 2 + earthRadius ** 2 - 2 * dist * earthRadius * c);
        const sinE = (dist * c - earthRadius) / denom;
        return Math.asin(mf.clamp(sinE));
    }

    function f(psiRad: number): number {
        return elevationRad(psiRad) - eMinRad;
    }

    // Bracket the root.
    let lo = 0;
    let hi = Math.PI - 1e-12;

    // We expect f(0) > 0, but handle edge cases.
    const fLo = f(lo);
    if (fLo <= 0) {return 0;}

    // Find a tighter hi near 90° first, then expand if needed.
    hi = Math.min((Math.PI / 2) + Math.abs(eMinRad) + (2 * Math.PI / 180), Math.PI - 1e-12);
    while (f(hi) > 0 && hi < Math.PI - 1e-12) {
        hi = Math.min(hi + 5 * degToRad, Math.PI - 1e-12); // step 5°
        if (hi >= Math.PI - 1e-12) break;
    }
    if (f(hi) > 0) {
        // If even near π it's still above cutoff (shouldn't happen for realistic inputs),
        // return almost-π in degrees.
        return (Math.PI - 1e-12) / degToRad;
    }

    // Bisection
    for (let i = 0; i < 80; i++) {
        const mid = 0.5 * (lo + hi);
        const fm = f(mid);
        if (fm > 0) lo = mid;
        else hi = mid;
    }

    return (lo + hi) / (2 * degToRad);
}

export function plotSvg(unix: number, numPoints: number = NUM_POINTS, width = 1000, height = 500): 
string {
    const header = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="-180 -90 360 180">\n`;
    const title = `<title>Global Sun and Moon Map - ${DateTime.fromMillis(unix, {zone: "utc"}).toISO()} - SunCompass</title>\n`;
    const [sLat, sLong] = subsolarPoint(generateLODProfile(unix), true);
    const antisolarLong = mf.clamp(sLong + (sLong >= 0 ? -180 : 180), -179.9999, 179.9999);

    let [mLat, mLong] = sublunarPoint(unix, true);
    mLong = mf.clamp(mLong, -179.9999, 179.9999);
    const moonAngularRadius = angularRadiusVisibility(90.8333, moonDistance(unix, true));

    return header + title + baseWorldMap() +
    svgCircleEquirectangular(mLat, mLong, Math.min(89.9999, moonAngularRadius), numPoints, 255, 255, 0, 0.25) +
    svgCircleEquirectangular(-sLat, antisolarLong, 89.1667, numPoints) + 
    svgCircleEquirectangular(-sLat, antisolarLong, 84, numPoints) + 
    svgCircleEquirectangular(-sLat, antisolarLong, 78, numPoints) + 
    svgCircleEquirectangular(-sLat, antisolarLong, 72, numPoints) + 
    moonIcon(mLong, mLat, 7.5) + sunIcon(sLong, sLat, 7.5) + `</g>\n</svg>\n`;
}