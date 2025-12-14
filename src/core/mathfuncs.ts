import { DateTime } from 'luxon';
import {DAY_LENGTH, earthERadius, flattening, J2000UTC, degToRad} from "./constants.ts";

export type Point = [number, number];
export type Polygon = Point[];
export type Polyline = Point[];

/** Divide x by y, rounding the output to the nearest integer with smaller absolute value. */
export function intDiv(x: number, y: number) {
    if (x<0) {return Math.ceil(x/y);}
    else {return Math.floor(x/y);}
}

/** Same as the toMillis() method in Luxon, but truncated to the nearest integer. */
export function ms(date: DateTime) {return Math.trunc(date.toMillis());}

/** Converts Unix milliseconds to approximate "fractional year" for Delta T calculation. */
function fractionalYear(t: number) {
    return t / (365.2425 * DAY_LENGTH) + 1970;
}

/**
 * Adjusts the given elevation angle (elev) of a celestial object to account for atmospheric refraction.
 */
export function refract(elev: number): number {
    const refraction = (elev <= -5/6) ? -0.0089931/Math.tan(elev*degToRad) : (1.02/Math.tan((elev+10.3/(elev+5.11))*degToRad)+0.001927)/60;
    return elev + refraction;
}

/** This function finds the approximate value of ΔT, which is calculated using the formula ΔT = TT - UT1. 
 * 
 * TT is terrestrial time (based on atomic clocks) and UT1 is mean solar time at 0° longitude.
 * 
 * This function's margin of error is 4.8 seconds in 2024, based on the value this function returns (73.8 seconds) versus the
 * real value (69 seconds). The margin of error increases for years before 1800 and after 2100, as the Earth's rotation 
 * varies unpredictably.
 * 
 * Source: https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
 * @param t Unix time in milliseconds. */
export function approxDeltaT(t: number) {
    const y = fractionalYear(t);
    if (y < 500) {
        const u = y/100;
        return 10583.6 - 1014.41 * u + 33.78311 * u**2 - 5.952053 * u**3 - 0.1798452 * u**4 + 0.022174192 * u**5 + 0.0090316521 * u**6;
    }
    else if (y < 1600) {
        const u = (y-1000)/100;
        return 1574.2 - 556.01*u + 71.23472*u**2 + 0.319781*u**3 - 0.853463*u**4 - 0.005050998*u**5 + 0.0083572073*u**6;
    }
    else if (y < 1700) {
        const t = y - 1600;
        return 120 - 0.9808*t - 0.01532*t**2 + t**3/7129;
    }
    else if (y < 1800) {
        const t = y - 1700;
        return 8.83 + 0.1603*t - 0.0059285*t**2 + 0.00013336*t**3 - t**4/1174000;
    }
    else if (y < 1860) {
        const t = y - 1800;
        return 13.72 - 0.332447 * t + 0.0068612 * t**2 + 0.0041116 * t**3 - 0.00037436 * t**4 + 0.0000121272 * t**5 - 0.0000001699 * t**6 + 0.000000000875 * t**7;
    }
    else if (y < 1900) {
        const t = y - 1860;
        return 7.62 + 0.5737 * t - 0.251754 * t**2 + 0.01680668 * t**3 - 0.0004473624 * t**4 + t**5 / 233174;
    }
    else if (y < 1920) {
        const t = y - 1900;
        return -2.79 + 1.494119 * t - 0.0598939 * t**2 + 0.0061966 * t**3 - 0.000197 * t**4;
    }
    else if (y < 1941) {
        const t = y - 1920;
        return 21.20 + 0.84493*t - 0.076100 * t**2 + 0.0020936 * t**3;
    }
    else if (y < 1961) {
        const t = y - 1950;
        return 29.07 + 0.407*t - t**2/233 + t**3 / 2547;
    }
    else if (y < 1986) {
        const t = y - 1975;
        return 45.45 + 1.067*t - t**2/260 - t**3 / 718;
    }
    else if (y < 2005) {
        const t = y - 2000;
        return 63.86 + 0.3345 * t - 0.060374 * t**2 + 0.0017275 * t**3 + 0.000651814 * t**4 + 0.00002373599 * t**5;
    }
    else if (y < 2050) {
        const t = y - 2000;
        return 62.92 + 0.32217 * t + 0.005589 * t**2;
    }
    else if (y < 2150) {return 32 * ((y-1820)/100)**2 - 0.5628 * (2150 - y) - 20;}
    else {
        const u = (y-1820)/100;
        return 32*u**2-20;
    }
}

/** Clamps a number to the range [min, max]. 
 * If min and max are not specified, they default to -1 and 1 respectively.*/
export function clamp(x: number, min=-1, max=1) {
    return Math.max(min, Math.min(x, max));
}

/** Calculates x modulo y, where the output is in the range [0, y). */
export function mod(x: number, y: number) {return ((x % y) + y) % y;}

/** Calculates the Julian century given the Unix timestamp in milliseconds, corrected for delta T. 
 * Note that there is a maximum error of 1 second due to the difference between UT1 and UTC, known as DUT1.
*/
export function jCentury(unix: number) {
    const deltaT0 = Math.round(approxDeltaT(J2000UTC) * 1000);
    let epoch = J2000UTC - deltaT0;

    const deltaT = Math.round(approxDeltaT(unix)*1000) - deltaT0;
    const millis = unix - epoch + deltaT;
    return millis / 3.15576e12; // There are 3.15576e12 milliseconds in a Julian century.
}

/** Calculates Julian day from Unix milliseconds. */
export function jdUTC(unix: number) {return unix / DAY_LENGTH + 2440587.5;}

/**
 * Returns the compass point (ex: NE, SSW) given a compass bearing in degrees.
 * @param bearing Compass bearing, in degrees clockwise from north.
 * @returns Compass point (either N, NNE, NE, ENE, E, ESE, SE, SSE, S, SSW, SW, WSW, W, WNW, NW, or NNW)
 */
export function direction(bearing: number) {
    if (bearing < 0 || bearing >= 360) {bearing = mod(bearing, 360);}
    if (bearing < 11.25) {return "N";}
    else if (bearing < 33.75) {return "NNE";}
    else if (bearing < 56.25) {return "NE";}
    else if (bearing < 78.75) {return "ENE";}
    else if (bearing < 101.25) {return "E";}
    else if (bearing < 123.75) {return "ESE";}
    else if (bearing < 146.25) {return "SE";}
    else if (bearing < 168.75) {return "SSE";}
    else if (bearing < 191.25) {return "S";}
    else if (bearing < 213.75) {return "SSW";}
    else if (bearing < 236.25) {return "SW";}
    else if (bearing < 258.75) {return "WSW";}
    else if (bearing < 281.25) {return "W";}
    else if (bearing < 303.75) {return "WNW";}
    else if (bearing < 326.25) {return "NW";}
    else if (bearing < 348.75) {return "NNW";}
    else {return "N";}
}

/** Returns the time of day in the DateTime as a number of milliseconds, from 0 (00:00:00.000) to 86399999 (23:59:59.999). */
export function convertToMS(date: DateTime) {
    return 1000 * (date.hour * 3600 + date.minute * 60 + date.second) + date.millisecond;
}

/**
 * Calculates one of the roots of a quadratic function, given two points on the quadratic curve and the derivative at one point.
 * @param x0 Initial x-coordinate
 * @param x1 Final x-coordinate (x1 > x0)
 * @param y0 Value of y at x=x0
 * @param y1 Value of y at x=x1
 * @param d0 Value of dy/dx at x=x0
 * @returns The root of the quadratic function between x0 and x1. Note that y0 and y1 should have opposite signs, and there
 * should be exactly one real root between x0 and x1.
 */
export function quadraticZero(x0: number, x1: number, y0: number, y1: number, d0: number) {
    const dx = x1 - x0;
    // y = a*u**2 + b*u + c where u = x-x0
    const a = (y1 - y0 - d0 * dx) / (dx ** 2);
    const b = d0;
    const c = y0;
    const d1 = 2*a*dx + b; // dy/dx at x=x1

    if (y0 === y1) { // constant (y0 = y1 = 0)
        return x0;
    }
    else if (Math.abs(Math.abs(d1 / d0) - 1) <= 1e-6) { // linear
        const frac = y0 / (y0 - y1);
        return clamp(x0 + frac * (x1 - x0), x0, x1);
    }
    else { // quadratic
        const disc = Math.max(0, b**2 - 4*a*c);
        const sign = (y1 > y0) ? 1 : -1;
        const u = (sign * Math.sqrt(disc) - b) / (2 * a);
        return clamp(x0 + u, x0, x1);
    }
}

/** Returns a Luxon DateTime corresponding to the beginning of the given day. */
export function startOfDay(date: DateTime) {
    return date.set({hour: 0, minute: 0, second: 0, millisecond: 0});
}

/** Returns a Luxon DateTime corresponding to the beginning of the next day. */
export function startNextDay(date: DateTime) {
    return date.plus({days: 1}).set({hour: 0, minute: 0, second: 0, millisecond: 0});
}

/** Signed "twice area" of triangle p0,p1,p2 (zero -> collinear). Helper function for isCollinear */
function cross(p0: number[], p1: number[], p2: number[]): number {
    const ax = p1[0] - p0[0], ay = p1[1] - p0[1];
    const bx = p2[0] - p0[0], by = p2[1] - p0[1];
    return ax * by - ay * bx;
}

/** Distance from p1 to the infinite line through p0–p2. Helper function for isCollinear. */
function pointLineDistance(p0: number[], p1: number[], p2: number[]): number {
    const area2 = Math.abs(cross(p0, p1, p2));
    const len = Math.hypot(p2[0] - p0[0], p2[1] - p0[1]);
    return len === 0 ? Math.hypot(p1[0] - p0[0], p1[1] - p0[1]) : area2 / len;
}

/** True if p1 is collinear with p0–p2 within epsilon (absolute distance in same units as coords) */
export function isCollinear(p0: number[], p1: number[], p2: number[], epsilon = 1e-6): boolean {
    return pointLineDistance(p0, p1, p2) <= epsilon;
}

/** Like toFixed() function in JavaScript/TypeScript, but removes trailing zeroes. */
export function toFixedS(n: number, precision: number) {
    if (precision === 0) {return n.toFixed(0);}
    else {return n.toFixed(precision).replace(/\.?0+$/, "");}
}

/** Rotate x, y, z around z-axis by theta degrees using right hand rule.
 * Used to convert ECI coordinates to ECEF. Returns an array [x, y, z]
 */
export function rotateZ(x: number, y: number, z: number, theta: number) {
    const [cosT, sinT] = [Math.cos(theta*degToRad), Math.sin(theta*degToRad)];
    const x2 = x * cosT - y * sinT;
    const y2 = x * sinT + y * cosT;
    return [x2, y2, z];
}

/** Converts geodetic latitude and longitude to rectangular ECEF coordinates in km: [x, y, z] */
export function latLongEcef(lat: number, long: number): number[] {
    const e2 = 2*flattening - flattening**2;
    lat *= degToRad; long *= degToRad;
    const [sinLat, cosLat, sinLong, cosLong] = [Math.sin(lat),Math.cos(lat),Math.sin(long),Math.cos(long)];
    const N = earthERadius / Math.sqrt(1 - e2 * sinLat**2); // radius of curvature in prime vertical
    const [X, Y, Z] = [N*cosLat*cosLong, N*cosLat*sinLong, N*(1-e2)*sinLat]; // observer's ECEF coords
    return [X, Y, Z];
}

/** Converts geocentric latitude, longitude, and distance (in kilometers) to rectangular ECEF coordinates.
 * @param lat Geocentric latitude. (To convert from geodetic, use geodetic2geocentric)
 * @param long Longitude.
 * @param dist Distance from Earth's center in kilometers.
 * @returns ECEF coordinate array: [x, y, z]
 */
export function toEcef(lat: number, long: number, dist: number): number[] {
    lat *= degToRad; long *= degToRad;
    const [cosLat, sinLat, cosLong, sinLong] = [Math.cos(lat), Math.sin(lat), Math.cos(long), Math.sin(long)];
    const [x, y, z] = [dist*cosLat*cosLong, dist*cosLat*sinLong, dist*sinLat];
    return [x, y, z];
}

/** Given the geodetic latitude, longitude, and ECEF of an observer, and the ECEF coordinates of a celestial object, find the
 * elevation and azimuth of the object.
 * @param lat Geodetic latitude of observer.
 * @param long Longitude of observer.
 * @param ecefO ECEF coordinates of observer.
 * @param ecefC ECEF coordinates of celestial object (planet, moon, star).
 * @returns [elevation, azimuth] of celestial object as seen from observer. Both are given in degrees. Elevation is in degrees
 * above the horizon, and azimuth is in degrees clockwise from north (range 0 <= a < 360).
 */
export function elevAzimuth(lat: number, long: number, ecefO: number[], ecefC: number[]): number[] {
    const [xe, ye, ze] = ecefC; // celestial body's ECEF
    const [xo, yo, zo] = ecefO; // observer's ECEF
    const [dx, dy, dz] = [xe-xo, ye-yo, ze-zo];

    // rotate ECEF coordinates to local ENU (east, north, up) at observer
    const [latR, longR] = [lat*degToRad, long*degToRad];
    const [sinLat, cosLat, sinLong, cosLong] = [Math.sin(latR), Math.cos(latR), Math.sin(longR), Math.cos(longR)];
    const E = -sinLong * dx + cosLong * dy;
    const N = -sinLat * cosLong * dx - sinLat * sinLong * dy + cosLat * dz;
    const U =  cosLat * cosLong * dx + cosLat * sinLong * dy + sinLat * dz;

    // convert ENU to elevation and azimuth
    const R = Math.hypot(E, N, U);
    const elev = Math.asin(clamp(U / R)) / degToRad; // elevation above horizon
    const az = mod(Math.atan2(E, N) / degToRad, 360); // degrees clockwise from north
    return [elev, az];
}

/** Given a start date and an end date, both with the same IANA time zone identifier, return an array of Luxon DateTimes with
 * the start of each day within the interval. */
export function dayStarts(start: DateTime, end: DateTime): DateTime[] {
    if (start.zoneName !== end.zoneName) {
        console.log("Start and end must have same time zone");
        return [];
    }
    const dayStarts = [];
    let cur = start.startOf("day");
    while (ms(cur) <= ms(end)) {
        dayStarts.push(cur);
        cur = cur.plus({days: 1}).startOf("day");
    }
    return dayStarts;
}

type Seg = { a: number[]; b: number[] };
// EPS for key-stability with fractional coords
const SNAP = 1e-6;
const snap = (v: number) => Math.round(v / SNAP) * SNAP;
/** Convert a time of day in milliseconds to hh:mm:ss in either 12 or 24-hour format. */
export function convertToHMS(timeOfDay: number, twentyFourHours: boolean) {
    const timeOfDayS = Math.floor(timeOfDay / 1000);
    const second = mod(timeOfDayS, 60);
    const minute = Math.floor(mod(timeOfDayS/60, 60));
    const hour24 = Math.floor(timeOfDayS/3600);
    const hour12 = mod(hour24 - 1, 12) + 1;
    const minString = String(minute).padStart(2, "0");
    const secString = String(second).padStart(2, "0");
    const hourString24 = String(hour24).padStart(2, "0");
    if (twentyFourHours) {return `${hourString24}:${minString}:${secString}`;} // ex: 09:47:29
    else if (hour24 <= 11) {return `${hour12}:${minString}:${secString} am`;} // ex: 9:47:29 am
    else {return `${hour12}:${minString}:${secString} pm`;} // ex: 9:47:29 pm
}

/** Convert a series of intervals into sets of points representing polygons. */
export function intervalsToPolygon(intervals: number[][][]): Polygon[] {
    const normalizeSpans = (spans: number[][]): number[][] => {
        if (!spans || spans.length === 0) return [];
        const s = spans.map(([a, b]) => [Math.min(a, b), Math.max(a, b)] as [number, number]).sort((u, v) => (u[0]-v[0]) || (u[1]-v[1]));
        const out: number[][] = [];
        for (const [a, b] of s) {
            if (out.length === 0 || a > out[out.length-1][1]) {out.push([a, b]);} 
            else {out[out.length - 1][1] = Math.max(out[out.length-1][1], b);}
        }
        return out;
    };
    // symmetric difference of two disjoint, sorted span lists
    const xorSpans = (A: number[][], B: number[][]): number[][] => {
        const evts: {y: number; d: number}[] = [];
        for (const [a, b] of A) {evts.push({ y: a, d: +1 }, { y: b, d: -1 });}
        for (const [a, b] of B) {evts.push({ y: a, d: +1 }, { y: b, d: -1 });}
        evts.sort((u, v) => (u.y - v.y) || (v.d - u.d)); // starts before ends at same y
        const out: number[][] = [];
        let inside = 0;
        let y0 = 0;
        for (const { y, d } of evts) {
            if (inside === 1) out.push([y0, y]);
            inside = (inside + d) & 1;
            y0 = y;
        }
        return out;
    };
    const keyOf = (p: number[]) => `${snap(p[0])}:${snap(p[1])}`;
    const segKey = (u: number[], v: number[]) => {
        const ux = snap(u[0]), uy = snap(u[1]);
        const vx = snap(v[0]), vy = snap(v[1]);
        return (ux < vx || (ux === vx && uy <= vy)) ? `${ux}:${uy}->${vx}:${vy}` : `${vx}:${vy}->${ux}:${uy}`;
    };
    const addAdj = (adj: Map<string, number[][]>, u: number[], v: number[]) => {
        const ku = keyOf(u), kv = keyOf(v);
        if (!adj.has(ku)) adj.set(ku, []);
        if (!adj.has(kv)) adj.set(kv, []);
        adj.get(ku)!.push([snap(v[0]), snap(v[1])]);
        adj.get(kv)!.push([snap(u[0]), snap(u[1])]);
    };

    // build segments
    const W = intervals.length;
    const cols = Array.from({ length: W }, (_, x) => normalizeSpans(intervals[x] || []));
    const horizontal: Seg[] = [];
    for (let x = 0; x < W; x++) {
        for (const [a, b] of cols[x]) {
            horizontal.push({ a: [x, a], b: [x + 1, a] }); // bottom cap
            horizontal.push({ a: [x, b], b: [x + 1, b] }); // top cap
        }
    }

    const vertical: Seg[] = [];
    for (let x = 0; x <= W; x++) {
        const L = x > 0 ? cols[x - 1] : [];
        const R = x < W ? cols[x] : [];
        const diff = xorSpans(L, R);
        for (const [a, b] of diff) {vertical.push({ a: [x, a], b: [x, b] });}
    }
    const segs: Seg[] = vertical.concat(horizontal).map(({ a, b }) => ({a: [snap(a[0]), snap(a[1])], b: [snap(b[0]), snap(b[1])]}));

    // stitch into rings
    const adj = new Map<string, Polygon>();
    for (const { a, b } of segs) {addAdj(adj, a, b);}
    const used = new Set<string>();
    const polygons: Polygon[] = [];

    for (const { a, b } of segs) {
        const startEdge = segKey(a, b);
        if (used.has(startEdge)) continue;

        // seed walk with the exact edge (a -> b)
        used.add(startEdge);
        const polygon: Polygon = [];
        polygon.push([snap(a[0]), snap(a[1])]);
        let prev: Point = [snap(a[0]), snap(a[1])];
        let curr: Point = [snap(b[0]), snap(b[1])];

        while (true) {
            polygon.push(curr);
            const nbrs: Polygon = adj.get(keyOf(curr)) || [];
            // choose neighbor that's NOT prev, prefer the one whose edge isn't used yet
            let next: Point | undefined;
            for (const cand of nbrs) {
                if (cand[0] === prev[0] && cand[1] === prev[1]) continue;
                const k = segKey(curr, cand);
                if (!used.has(k)) { next = cand; break; }
            }
            if (!next) break; // should not happen if all loops are closed

            used.add(segKey(curr, next));
            prev = curr;
            curr = next;
            if (curr[0] === polygon[0][0] && curr[1] === polygon[0][1]) break; // closed
        }
        if (polygon.length >= 4) {polygons.push(polygon);}
    }
    return polygons;
}