import { DateTime } from 'luxon';
import {DAY_LENGTH, earthERadius, flattening, J2000UTC, deltaT, earthPRadius, TAU, HORIZON} from "./constants.ts";

export type Point = [number, number];
export type Polygon = Point[];
export type Polyline = Point[];
export type HourAngleFn = (long: number, unix: number) => number; 

/** Same as the toMillis() method in Luxon, but truncated to the nearest integer. */
export function ms(date: DateTime) {return Math.trunc(date.toMillis());}

/** Converts Unix milliseconds to approximate "fractional year" for Delta T calculation. */
function fractionalYear(t: number) {
    return t / (365.2425 * DAY_LENGTH) + 1970;
}

/** Adjusts the given elevation angle (elev) of a celestial object to account for atmospheric refraction. */
export function refract(elev: number): number {
    const refraction = (elev <= HORIZON) ? 
    -1.569584402829295e-4/Math.tan(elev) : 
    2.96706e-4/Math.tan(elev + 0.00313756/(0.0891863 + elev)) + 5.608094753118464e-7;
    return elev + refraction;
}

/** This function finds the approximate value of ΔT in seconds, which is calculated using the formula ΔT = TT - UT1. TT is
 * terrestrial time (based on atomic clocks) and UT1 is mean solar time at 0° longitude.
 * 
 * @param t Unix time in milliseconds. */
export function approxDeltaT(t: number) {
    const y = fractionalYear(t);
    if (y <= deltaT[0][0]) {return deltaT[0][1];}
    else {
        for (let i=0; i<deltaT.length-1; i++) {
            if (deltaT[i][0] <= y && deltaT[i+1][0] >= y) {
                const frac = (y - deltaT[i][0]) / (deltaT[i+1][0] - deltaT[i][0]);
                const diff = deltaT[i+1][1] - deltaT[i][1];
                return deltaT[i][1] + diff * frac;
            }
        }
        return deltaT[deltaT.length-1][1];
    }
}

/** Clamps a number to the range [min, max]. 
 * If min and max are not specified, they default to -1 and 1 respectively.*/
export function clamp(x: number, min=-1, max=1) {
    return Math.max(min, Math.min(x, max));
}

/** Calculates x modulo y, where the output is in the range [0, y). */
export function mod(x: number, y: number) {return ((x % y) + y) % y;}

/** Evaluates the polynomial with given coefficients at the given x-value, using Horner's method.
 * Example: if coefficients are [2, -3.4, 5.67, -8.9], the polynomial evaluated is 2 - 3.4\*x + 5.67\*x^2 - 8.9\*x^3.
 * The degree of the polynomial is equal to coefficients.length - 1.
 */
export function polynomial(x: number, coefficients: number[]): number {
    let result = coefficients[coefficients.length - 1];
    for (let i = coefficients.length - 2; i >= 0; i--) {
        result = result * x + coefficients[i];
    }
    return result;
}

/** Evaluates the polynomial modulo a specific number. */
export function polymod(x: number, coefficients: number[], modulus: number): number {
    return mod(polynomial(x, coefficients), modulus);
}

/** Calculates the Julian century given the Unix timestamp in milliseconds, corrected for delta T. 
 * Note that there is a maximum error of 1 second due to the difference between UT1 and UTC, known as DUT1.
*/
export function jCentury(unix: number) {
    const deltaT0 = Math.round(approxDeltaT(J2000UTC) * 1000);
    let epoch = J2000UTC - deltaT0;

    const deltaT = Math.round(approxDeltaT(unix) * 1000) - deltaT0;
    const millis = unix - epoch + deltaT;
    return millis / 3.15576e12; // There are 3.15576e12 milliseconds in a Julian century.
}

/** Calculates Julian day from Unix milliseconds. */
export function jdUTC(unix: number) {return unix / DAY_LENGTH + 2440587.5;}

/**
 * Returns the compass point (ex: NE, SSW) given a compass bearing in radians.
 * @param bearing Compass bearing, in radians clockwise from north.
 * @returns Compass point (either N, NNE, NE, ENE, E, ESE, SE, SSE, S, SSW, SW, WSW, W, WNW, NW, or NNW)
 */
export function direction(bearing: number) {
    const directions = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"];
    const adjBearing = mod(bearing + TAU/32, TAU);
    const index = clamp(Math.floor(adjBearing / TAU * 16), 0, 15);
    return directions[index];
}

/**
 * Calculates one of the roots of a quadratic function, given three points on the curve.
 * @param x0 Initial x-coordinate
 * @param y0 y-coordinate at x0
 * @param x1 Middle x-coordinate (x1 > x0)
 * @param y1 y-coordinate at x1
 * @param x2 Final x-coordinate (x2 > x1)
 * @param y2 y-coordinate at x2
 * @returns The root of the quadratic function between x0 and x2. If it is undefined, returns (x0 + x2) / 2 as a placeholder. 
 * If there are two roots, it returns the root which is closer to the center of the range. For SunCompass, there should always 
 * be one solution in the range.
 */
export function quadraticZero(x0: number, y0: number, x1: number, y1: number, x2: number, y2: number): number {
    const u1 = (x1 - x0) / (x2 - x0); // normalize to [0, 1] for numerical stability, u0 = 0 and u2 = 1
    // y = a*u^2 + b*u + c where u = (x - x0) / (x2 - x0)
    const a = y0/u1 - y1/(u1*(1-u1)) + y2/(1-u1);
    const b = -y0*(u1+1)/u1 + y1/(u1*(1-u1)) - y2*u1/(1-u1);
    const c = y0;
    if (Math.abs(2 * a / b) < 1e-4) { // if almost linear
        if (y0 === y2) {return (x0 + x2) / 2;} // placeholder when undefined
        else {return x0 + (y0 / (y0 - y2)) * (x2 - x0);}
    }
    else {
        const disc = b**2 - 4*a*c;
        if (disc < 0) {return (x0 + x2) / 2;} // placeholder when undefined
        else {
            const sol1 = (Math.sqrt(disc) - b) / (2 * a);
            const sol2 = (-Math.sqrt(disc) - b) / (2 * a);
            const solution = (Math.abs(sol1 - 0.5) <= Math.abs(sol2 - 0.5)) ? sol1 : sol2;
            return x0 + solution * (x2 - x0);
        }
    }
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

/** Rotate x, y, z around z-axis by theta radians using right hand rule.
 * Used to convert ECI coordinates to ECEF, by rotating by -sc.gast(). Returns an array [x, y, z]
 */
export function rotateZ(x: number, y: number, z: number, theta: number) {
    const [cosT, sinT] = [Math.cos(theta), Math.sin(theta)];
    const x2 = x * cosT - y * sinT;
    const y2 = x * sinT + y * cosT;
    return [x2, y2, z];
}

/** Converts latitude and longitude in radians to rectangular ECEF coordinates in km: [x, y, z]. */
export function latLongEcef(lat: number, long: number): number[] {
    const e2 = 2*flattening - flattening**2;
    const [sinLat, cosLat, sinLong, cosLong] = [Math.sin(lat),Math.cos(lat),Math.sin(long),Math.cos(long)];
    const N = earthERadius / Math.sqrt(1 - e2 * sinLat**2); // radius of curvature in prime vertical
    const [X, Y, Z] = [N*cosLat*cosLong, N*cosLat*sinLong, N*(1-e2)*sinLat]; // observer's ECEF coords
    return [X, Y, Z];
}

/** Converts latitude, longitude, and distance (in kilometers) to rectangular ECEF coordinates.
 * @param lat Latitude, in radians
 * @param long Longitude, in radians
 * @param dist Distance from Earth's center in kilometers.
 * @returns ECEF coordinate array: [x, y, z]
 */
export function toEcef(lat: number, long: number, dist: number): number[] {
    const [cosLat, sinLat, cosLong, sinLong] = [Math.cos(lat), Math.sin(lat), Math.cos(long), Math.sin(long)];
    const [x, y, z] = [dist*cosLat*cosLong, dist*cosLat*sinLong, dist*sinLat];
    return [x, y, z];
}

/** Given the latitude, longitude, and ECEF of an observer, and the ECEF coordinates of a celestial object, find the
 * elevation and azimuth of the object.
 * @param lat Latitude of observer, in radians.
 * @param long Longitude of observer, in radians.
 * @param ecefO ECEF coordinates of observer.
 * @param ecefC ECEF coordinates of celestial object (planet, moon, star).
 * @returns [elevation, azimuth] of celestial object as seen from observer. Both are given in radians, and should be
 * multiplied by (360/TAU) to convert to degrees.
 */
export function elevAzimuth(lat: number, long: number, ecefO: number[], ecefC: number[]): number[] {
    const [xe, ye, ze] = ecefC; // celestial body's ECEF
    const [xo, yo, zo] = ecefO; // observer's ECEF
    const [dx, dy, dz] = [xe-xo, ye-yo, ze-zo];

    // rotate ECEF coordinates to local ENU (east, north, up) at observer
    const [sinLat, cosLat, sinLong, cosLong] = [Math.sin(lat), Math.cos(lat), Math.sin(long), Math.cos(long)];
    const E = -sinLong * dx + cosLong * dy;
    const N = -sinLat * cosLong * dx - sinLat * sinLong * dy + cosLat * dz;
    const U =  cosLat * cosLong * dx + cosLat * sinLong * dy + sinLat * dz;

    // convert ENU to elevation and azimuth
    const R = Math.hypot(E, N, U);
    const elev = Math.asin(clamp(U / R)); // elevation above horizon
    const az = mod(Math.atan2(E, N), TAU); // radians clockwise from north
    return [elev, az];
}

/** Finds the subpoint of a celestial body - the [latitude, longitude] at which the body is directly overhead - given the
 * body's ECEF coordinates.
 * Latitude and longitude are in radians and are in the range [-TAU/4, TAU/4] and [-TAU/2, TAU/2) respectively.
 */
export function subpoint(ecef: number[]): number[] {
    const a = earthERadius;
    const b = earthPRadius;
    const e2 = 1 - b**2 / a**2;
    
    function F(k: number): number {
        const dx = 1 + k / a**2;
        const dz = 1 + k / b**2;
        const x = ecef[0] / dx, y = ecef[1] / dx, z = ecef[2] / dz;
        return (x**2 + y**2) / a**2 + z**2 / b**2 - 1;
    }

    // Bracket root. k=0 => point at moon vector itself (outside ellipsoid) => F(0) > 0
    // For large k, x,y,z shrink toward 0 => F(k) -> -1. So root exists.
    let lo = 0;
    let hi = b ** 2;
    for (let i=0; i<200 && F(hi) > 0; i++) {hi *= 2;}

    // Bisection
    for (let i = 0; i < 80; i++) {
        const mid = (lo + hi) / 2;
        if (F(mid) > 0) {lo = mid;} 
        else {hi = mid;}
    }
    
    const k = (lo + hi) / 2;
    const dx = 1 + k / a**2;
    const dz = 1 + k / b**2;
    const x = ecef[0] / dx, y = ecef[1] / dx, z = ecef[2] / dz;

    const p = Math.hypot(x, y);
    if (p === 0) { // poles - longitude arbitrary so we use 0
        return [(z >= 0) ? TAU/4 : -TAU/4, 0];
    }

    const lon = Math.atan2(y, x);
    let lat = Math.atan2(z, p*(1-e2)); // initial guess
    let h = 0;

    for (let i=0; i<10; i++) {
        const sinLat = Math.sin(lat);
        const N = a / Math.sqrt(1 - e2 * sinLat ** 2);
        h = p / Math.cos(lat) - N;
        const latNext = Math.atan2(z, p * (1 - e2 * (N / (N + h))));
        if (Math.abs(latNext - lat) < 1e-14) {lat = latNext; break;}
        lat = latNext;
    }

    return [clamp(lat, -TAU/4, TAU/4), mod(lon + TAU/2, TAU) - TAU/2];
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

export function meridianPassings(
    long: number, // longitude in radians
    start: number, // start of day, in Unix milliseconds
    end: number, // end of day, in Unix milliseconds
    targetHourAngle: number, // in radians (ex. 0 for solar noon, or +-TAU/2 for solar midnight)
    hAngle: HourAngleFn // function returning hour angle in radians, range [-TAU/2, TAU/2)
): number[] 
{
    const ha0 = hAngle(long, start), ha1 = hAngle(long, end);
    const haDiff = mod(ha1 - ha0 - TAU/2, TAU) + TAU/2; // change in hour angle (rad) during the day, range [0.5*tau, 1.5*tau)
    const targetAdjusted = mod(targetHourAngle - ha0, TAU); // diff between target hour angle and hour angle at 00:00

    const haRate = haDiff / (end - start); // rate at which the hour angle is changing, in radians per millisecond.
    const times = [];
    for (let i = targetAdjusted; i < haDiff; i += TAU) {
        let t = start + i / haRate;
        const ha = hAngle(long, t);
        const diff = mod(targetHourAngle - ha + TAU/2, TAU) - TAU/2;
        t += diff / haRate;
        times.push(Math.floor(clamp(t, start, end-1)));
    }
    return times;
}

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

type Seg = { a: number[]; b: number[] };
// EPS for key-stability with fractional coords
const SNAP = 1e-6;
const snap = (v: number) => Math.round(v / SNAP) * SNAP;
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