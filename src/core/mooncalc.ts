/** Formulas derived from "Astronomical Algorithms" by Jean Meeus. */

import * as mf from "./mathfuncs.ts";
import * as sc from "./suncalc.ts";
import {BSEARCH_GAP, DAY_LENGTH, degToRad, HORIZON, moonPtl, moonPtld} from "./constants.ts";
import { generateLODProfile, getTimeOfDay, timeZoneLookupTable, type SEvent, type TimeChange } from "./lookup-tables.ts";
import { DateTime } from "luxon";

export type MoonTable = {
    lunarEvents: SEvent[][];
    intervals: number[][][];
    phases: {
        unix: number;
        date: number;
        type: number;
    }[];
    perigees: number[];
    apogees: number[];
    timeZoneTable: TimeChange[];
}

export function moonMeanLongitude(JC: number): number {
    return mf.polynomial(JC, [3.810341023, 8399.7091135216, -2.3157e-5, 3.23904e-8, -2.67713e-10], 2*Math.PI);
}

export function moonMeanElongation(JC: number): number {
    return mf.polynomial(JC, [5.198466741, 7771.3771468129, -2.84489e-5, 3.19735e-8, -1.54365e-10], 2*Math.PI);
}

export function moonMeanAnomaly(JC: number): number {
    return mf.polynomial(JC, [2.355555899, 8328.6914269548, 1.57027e-4, 2.50410e-7, -1.18633e-9], 2*Math.PI);
}

/** Moon argument of latitude */
export function moonArgLat(JC: number): number {
    return mf.polynomial(JC, [1.627905233, 8433.4661581307, -5.93918e-5, -4.94988e-9, 2.02167e-11], 2*Math.PI);
}

/** Sum of all longitude terms in moonPtld (periodic terms for longitude and distance) */
function l(JC: number): number {
    let l = 0;
    const D = moonMeanElongation(JC);
    const M = sc.meanSunAnomaly(JC);
    const Mp = moonMeanAnomaly(JC);
    const F = moonArgLat(JC);
    const E = 1 - 0.002516*JC - 7.4e-6*JC**2;
    for (const curRow of moonPtld) {
        // E ** Math.abs(curRow[1]) multiplies terms which contain -M or M by E, and terms which contain 2M or -2M by E^2.
        l += curRow[4] * Math.sin(curRow[0]*D + curRow[1]*M + curRow[2]*Mp + curRow[3]*F) * E ** Math.abs(curRow[1]);
    }
    return l;
}

/** Sum of all distance terms in moonPtld (periodic terms for longitude and distance) */
function r(JC: number): number {
    let r = 0;
    const D = moonMeanElongation(JC);
    const M = sc.meanSunAnomaly(JC);
    const Mp = moonMeanAnomaly(JC);
    const F = moonArgLat(JC);
    const E = 1 - 0.002516*JC - 7.4e-6*JC**2;
    for (const curRow of moonPtld) {
        r += curRow[5] * Math.cos(curRow[0]*D + curRow[1]*M + curRow[2]*Mp + curRow[3]*F) * E ** Math.abs(curRow[1]);
    }
    return r;
}

/** Derivative of moon-earth distance with respect to time, in kilometers per Julian century. */
function rDeriv(JC: number) {
    let rDeriv = 0;
    const D = moonMeanElongation(JC);
    const M = sc.meanSunAnomaly(JC);
    const Mp = moonMeanAnomaly(JC);
    const F = moonArgLat(JC);
    const E = 1 - 0.002516*JC - 7.4e-6*JC**2;
    const dD = mf.polynomial(JC, [7771.3771468129, 2*-2.84489e-5, 3*3.19735e-8, 4*-1.54365e-10]);
    const dM = mf.polynomial(JC, [628.3019551672, 2*-2.68083e-6, 3*7.1267e-10]);
    const dMp = mf.polynomial(JC, [8328.6914269548, 2*1.57027e-4, 3*2.50410e-7, 4*-1.18633e-9]);
    const dF = mf.polynomial(JC, [8433.4661581307, 2*-5.93918e-5, 3*-4.94988e-9, 4*2.02167e-11]);
    const dE = -0.002516 - 1.48e-5*JC;
    for (const curRow of moonPtld) {
        // Chain rule: if h(x) = f(g(x)), then h'(x) = g'(x) * f'(g(x)) 
        const dCosTerm = (curRow[0]*dD + curRow[1]*dM + curRow[2]*dMp + curRow[3]*dF) * 
            -Math.sin(curRow[0]*D + curRow[1]*M + curRow[2]*Mp + curRow[3]*F);
        const exp = Math.abs(curRow[1]);
        const dETerm = exp * E ** (exp - 1) * dE;
        // Product rule for derivatives: (x * y)' = x' * y + x * y'
        rDeriv += curRow[5] * (dCosTerm * E ** exp + Math.cos(curRow[0]*D + curRow[1]*M + curRow[2]*Mp + curRow[3]*F) * dETerm);
    }
    return rDeriv;
}

/** Sum of all latitude terms in moonPtl (periodic terms for latitude) */
function b(JC: number): number {
    let b = 0;
    const D = moonMeanElongation(JC);
    const M = sc.meanSunAnomaly(JC);
    const Mp = moonMeanAnomaly(JC);
    const F = moonArgLat(JC);
    const E = 1 - 0.002516*JC - 7.4e-6*JC**2;
    for (const curRow of moonPtl) {
        b += curRow[4] * Math.sin(curRow[0]*D + curRow[1]*M + curRow[2]*Mp + curRow[3]*F) * E ** Math.abs(curRow[1]);
    }
    return b;
}

function a(JC: number): number[] {
    return [2.09003 + 2.30120*JC, 0.926595 + 8364.739848*JC, 5.47073 + 8399.684725*JC];
}

/** Variations in longitude due to the actions of Venus, Jupiter, and the flattening of Earth. */
function deltaL(JC: number): number {
    const [a1, a2, a3] = a(JC);
    return 6.90801e-5*Math.sin(a1) + 3.42434e-5*Math.sin(moonMeanLongitude(JC)-moonArgLat(JC)) + 5.55015e-6*Math.sin(a2);
}

/** Variations in latitude due to the actions of Venus, Jupiter, and the flattening of Earth. */
function deltaB(JC: number): number {
    const [a1, a2, a3] = a(JC);
    const meanLong = moonMeanLongitude(JC);
    const meanAnomaly = moonMeanAnomaly(JC);
    const argLat = moonArgLat(JC);
    return -3.90081e-5*Math.sin(meanLong) + 6.66716e-6*Math.sin(a3) + 3.05433e-6*Math.sin(a1-argLat)
    + 3.05433e-6*Math.sin(a1+argLat) + 2.21657e-6*Math.sin(meanLong-meanAnomaly) - 2.00713e-6*Math.sin(meanLong+meanAnomaly);
}

/** Returns the ecliptic latitude and longitude of the moon. Return value is an array: [latitude, longitude], 
 * both measured in radians. */
export function moonLatLong(date: number, unix = false): number[] {
    if (!unix) {
        const lat = mf.clamp(b(date) + deltaB(date), -Math.PI/2, Math.PI/2);
        const long = mf.mod(moonMeanLongitude(date) + l(date) + deltaL(date) + sc.longNutation(date), 2*Math.PI);
        return [lat, long];
    }
    else {return moonLatLong(mf.jCentury(date));}
}

/** Distance from center of earth to center of moon, in kilometers. */
export function moonDistance(date: number, unix = false) {
    if (!unix) {return 385000.56 + r(date);}
    else {return moonDistance(mf.jCentury(date));}
}

/** Derivative of moon-earth distance, in kilometers per second. Used to calculate perigee and apogee. */
export function moonDistanceDeriv(date: number, unix = false) {
    if (!unix) {return rDeriv(date) / (36525 * 86400);}
    else {return moonDistanceDeriv(mf.jCentury(date));}
}

/** Returns the rectangular coordinates [x, y, z] in Earth-centered, Earth-fixed coordinates (ECEF) in kilometers. */
export function moonEcef(unix: number): number[] {
    let [eLat, eLong] = moonLatLong(unix, true);
    const ob = sc.obliquity(unix, true);
    const [sinB,cosB,sinL,cosL,sinE,cosE] = [Math.sin(eLat),Math.cos(eLat),Math.sin(eLong),Math.cos(eLong),Math.sin(ob),Math.cos(ob)];
    
    // Geocentric equatorial (Earth-centered inertial) coordinates
    const dist = moonDistance(unix, true);
    const xeci = dist * cosB * cosL;
    const yeci = dist * (cosB * sinL * cosE - sinB * sinE);
    const zeci = dist * (cosB * sinL * sinE + sinB * cosE);

    // Convert to ECEF coordinates
    const rectCoords = mf.rotateZ(xeci, yeci, zeci, -sc.gast(unix));
    return rectCoords;
}

/** Returns the sublunar point [latitude, longitude].
 * @param unix Unix timestamp in milliseconds.
 * @param degrees Whether to return coordinates in degrees (true) or radians (false). Defaults to false.
 */
export function sublunarPoint(unix: number, degrees = false): number[] {
    const [lat, long] = mf.subpoint(moonEcef(unix));
    return degrees ? [lat / degToRad, long / degToRad] : [lat, long];
}

/** The hour angle of the moon at the given longitude and Unix timestamp, in radians between -pi and pi. 
 * This is equal to the observer's longitude minus the longitude of the sublunar point. */
export function moonHourAngle(longitude: number, unix: number): number {
    return mf.mod(longitude - sublunarPoint(unix)[1] + Math.PI, 2*Math.PI) - Math.PI;
}

/** Returns the moon's position: [elevation, azimuth] in radians. Optionally, the observer's ECEF can be specified in order
 * to avoid repeatedly computing it. */
export function moonPosition(lat: number, long: number, unix: number): number[] {
    return mf.elevAzimuth(lat, long, moonEcef(unix));
}

/**
 * The time at which the moon crosses the observer's meridian, analogous to solar noon for the sun.
 * @param lat Observer's latitude in radians
 * @param long Observer's longitude in radians
 * @param startUnix Start of day in Unix milliseconds
 * @param endUnix End of day in Unix milliseconds
 * @returns An SEvent[] object, containing the Unix time of the meridian passing, the type ("Meridian Passing"), and the
 * elevation and azimuth of the moon as seen from the observer.
 */
export function moonMeridianPassing(lat: number, long: number, startUnix: number, endUnix: number):
SEvent[] {
    const lunarNoonTimes = mf.meridianPassings(long, startUnix, endUnix, 0, moonHourAngle);
    const events: SEvent[] = [];
    for (const time of lunarNoonTimes) {
        const [e, a] = moonPosition(lat, long, time);
        events.push({unix: time, type: "Meridian Passing", elev: e, azimuth: a});
    }
    return events;
}

/** Approximate derivative of the moon's position in radians per second. */
export function moonDerivative(lat: number, long: number, unix: number): number {
    const t0 = unix - 500, t1 = unix + 500;
    const e0 = moonPosition(lat, long, t0)[0], e1 = moonPosition(lat, long, t1)[0];
    return e1 - e0;
}

/**
 * Calculates the extrema of the moon's elevation at a given location - i.e. the points at which the moon's elevation reaches
 * a relative maximum or minimum. Lunar equivalent of maxAndMin() for the sun.
 * @param lat Observer's latitude, radians
 * @param long Observer's longitude, radians
 * @param start Start of day, Unix milliseconds
 * @param end End of day, Unix milliseconds
 * @returns An array of Unix times, including the start of the day, all extrema in lunar elevation, and the start of the next day.
 */
function moonMaxMin(lat: number, long: number, start: number, end: number): number[] {
    const times = [start];
    const interval = (end - start) / 6;
    for (let t = start; t < end; t += interval) {
        let t0 = t, t1 = Math.min(t + interval, end);
        let d0 = moonDerivative(lat, long, t0), d1 = moonDerivative(lat, long, t1);
        if (d0 === 0) {times.push(t0);}
        else if (d0 * d1 < 0) { // derivative changes sign
            while (t1 - t0 > BSEARCH_GAP) {
                const tAvg = (t0+t1)/2;
                const dAvg = moonDerivative(lat, long, tAvg);
                if (d1 * dAvg <= 0) {t0 = tAvg; d0 = dAvg;}
                else {t1 = tAvg; d1 = dAvg;}
            }
            const t = t0 + (d0 / (d0 - d1)) * (t1 - t0); // linear interpolation on 5-minute window
            times.push(t);
        }
    }
    times.push(end);
    return times;
}

/**
 * Returns the time(s) of moonrise on a particular day.
 * @param lat Latitude of observer, radians
 * @param long Longitude of observer, radians
 * @param maxMin Results of moonMaxMin() for the given day and location
 * @param angle The threshold elevation angle, in radians (defaults to -5/1080*pi = -5/6°)
 * @returns An array of SEvent objects, with the time, elevation, and azimuth of the moon, and a type tag.
 */
export function moonrise(lat: number, long: number, maxMin: number[], angle: number = HORIZON): SEvent[] {
    const riseTimes = [];
    for (let i=0; i<maxMin.length-1; i++) {
        let t0 = maxMin[i], t1 = maxMin[i+1];
        let e0 = moonPosition(lat, long, t0)[0], e1 = moonPosition(lat, long, t1)[0];
        if (e0 <= angle && e1 >= angle) {
            while (t1 - t0 > BSEARCH_GAP) {
                const tAvg = (t0+t1)/2;
                const eAvg = moonPosition(lat, long, tAvg)[0];
                if (eAvg <= angle) {t0 = tAvg; e0 = eAvg;}
                else {t1 = tAvg; e1 = eAvg;}
            }
            const tAvg = (t0 + t1) / 2;
            const eAvg = moonPosition(lat, long, tAvg)[0];
            const t = Math.floor(mf.quadraticZero(t0, e0-angle, tAvg, eAvg-angle, t1, e1-angle));
            const [e, a] = moonPosition(lat, long, t);
            riseTimes.push({unix: t, type: "Moonrise", elev: e, azimuth: a});
        }
    }
    return riseTimes;
}

/**
 * Returns the time(s) of moonset on a particular day.
 * @param lat Latitude of observer, radians
 * @param long Longitude of observer, radians
 * @param maxMin Results of moonMaxMin() for the given day and location
 * @param angle The threshold elevation angle, in radians (defaults to -5/1080*pi = -5/6°)
 * @returns An array of SEvent objects, with the time, elevation, and azimuth of the moon, and a type tag.
 */
export function moonset(lat: number, long: number, maxMin: number[], angle: number = HORIZON): SEvent[] {
    const setTimes = [];
    for (let i=0; i<maxMin.length-1; i++) {
        let t0 = maxMin[i], t1 = maxMin[i+1];
        let e0 = moonPosition(lat, long, t0)[0], e1 = moonPosition(lat, long, t1)[0];
        if (e0 >= angle && e1 <= angle) {
            while (t1 - t0 > BSEARCH_GAP) {
                const tAvg = (t0+t1)/2;
                const eAvg = moonPosition(lat, long, tAvg)[0];
                if (eAvg >= angle) {t0 = tAvg; e0 = eAvg;}
                else {t1 = tAvg; e1 = eAvg;}
            }
            const tAvg = (t0 + t1) / 2;
            const eAvg = moonPosition(lat, long, tAvg)[0];
            const t = Math.floor(mf.quadraticZero(t0, e0-angle, tAvg, eAvg-angle, t1, e1-angle));
            const [e, a] = moonPosition(lat, long, t);
            setTimes.push({unix: t, type: "Moonset", elev: e, azimuth: a});
        }
    }
    return setTimes;
}

/** Returns the phase angle of the moon given Unix time, in radians.
 * The value returned is the same as the phase angle i in Astronomical Algorithms, chapter 46. It is the angle between the
 * sun-moon and earth-moon vectors. */
function phaseAngle(unix: number): number {
    const ecefSun = sc.sunEcef(generateLODProfile(unix));
    const ecefMoon = moonEcef(unix);
    const sunVector = [ecefSun[0]-ecefMoon[0], ecefSun[1]-ecefMoon[1], ecefSun[2]-ecefMoon[2]];
    const moonVector = [-ecefMoon[0], -ecefMoon[1], -ecefMoon[2]];

    const dot = sunVector[0] * moonVector[0] + sunVector[1] * moonVector[1] + sunVector[2] * moonVector[2];
    const normSun = Math.hypot(sunVector[0], sunVector[1], sunVector[2]);
    const normMoon = Math.hypot(moonVector[0], moonVector[1], moonVector[2]);
    const i = Math.acos(mf.clamp(dot / (normSun * normMoon)));
    return i;
}

/** Returns the fraction of the Moon's disk that is illuminated. This ranges from 0 to 1.
 * Note that illumination at full moon isn't exactly 1 (except during central lunar eclipses) because the moon's orbit is 
 * tilted relative to the ecliptic. Similarly, illumination at new moon isn't exactly 0 except during a total or annular 
 * solar eclipse. */
export function illumination(unix: number): number {
    return (1 + Math.cos(phaseAngle(unix))) / 2;
}

/** The difference between the ecliptic longitudes of the sun and moon, in radians, normalized to the range [0, 2*pi).
 * A value of 0 indicates a new moon, pi/2 is a first quarter, pi is a half moon and 3*pi/2 is a last quarter. */
export function moonSunLongDiff(unix: number): number {
    return mf.mod(moonLatLong(unix, true)[1] - sc.sunTrueLong(unix, true), 2*Math.PI);
}

/**
 * Returns the moon phases on a given day.
 * @param lat Latitude of observer, radians
 * @param long Longitude of observer, radians
 * @param start Start of the day in Unix milliseconds
 * @param end End of the day in Unix milliseconds
 * @returns If there is a new moon, full moon, or first/last quarter that day, returns an SEvent object with the Unix time of the
 * event, along with a type tag ("New Moon", "First Quarter", "Full Moon", or "Last Quarter") and the moon's elevation and azimuth.
 * Otherwise, returns null.
 */
export function moonPhase(lat: number, long: number, start: number, end: number): SEvent | null {
    const phaseNames = ["New Moon", "First Quarter", "Full Moon", "Last Quarter"];
    let p0 = moonSunLongDiff(start), p1 = moonSunLongDiff(end); // longitude differences betweeen sun and moon
    let t0 = start, t1 = end;
    if (p0 > p1) {p0 -= 2*Math.PI;}
    const i0 = Math.floor(2 * p0 / Math.PI), i1 = Math.floor(2 * p1 / Math.PI);
    if (i0 !== i1) { // if the longitude-difference angle crosses a multiple of pi/2 during the day
        const thresh = i1 * Math.PI / 2;
        while (t1 - t0 > BSEARCH_GAP) {
            const tAvg = (t0 + t1) / 2;
            let pAvg = mf.mod(moonSunLongDiff(tAvg) - p0, 2*Math.PI) + p0;
            if (pAvg <= thresh) {t0 = tAvg; p0 = pAvg;}
            else {t1 = tAvg; p1 = pAvg;}
        }
        const frac = (thresh - p0) / (p1 - p0);
        const t = Math.floor(t0 + frac * (t1 - t0));
        const [e, a] = moonPosition(lat, long, t);
        return {unix: t, type: phaseNames[i1], elev: e, azimuth: a}; 
    }
    else {return null;}
}

/**
 * Returns the moon phase on a given day.
 * @param start Start of the day in Unix milliseconds
 * @param end End of the day in Unix milliseconds
 * @returns A string. Either "New Moon", "Waxing Crescent", "First Quarter", "Waxing Gibbous", "Full Moon",
 * "Waning Gibbous", "Last Quarter", or "Waning Crescent".
 */
export function moonPhaseDay(start: number, end: number): string {
    const diff0 = moonSunLongDiff(start), diff1 = moonSunLongDiff(end);
    if (diff0 <= 0.5*Math.PI) {return (diff1 >= 0.5*Math.PI) ? "First Quarter" : "Waxing Crescent";}
    else if (diff0 <= Math.PI) {return (diff1 >= Math.PI) ? "Full Moon" : "Waxing Gibbous";}
    else if (diff0 <= 1.5*Math.PI) {return (diff1 >= 1.5*Math.PI) ? "Last Quarter" : "Waning Gibbous";}
    else {return (diff1 < diff0) ? "New Moon" : "Waning Crescent";}
}

/** Returns the moon's apsides (perigees and apogees) during a multi-day period, starting with start and ending with
 * end. Both start and end are given in Unix milliseconds. */
export function moonApsides(start: number, end: number) {
    const perigees: number[] = [], apogees: number[] = [];
    for (let t = start; t < end; t += DAY_LENGTH) {
        let t0 = t, t1 = Math.min(t + DAY_LENGTH, end);
        let d0 = moonDistanceDeriv(t0, true), d1 = moonDistanceDeriv(t1, true);
        if (d0 === 0) {(d1 > 0 ? perigees : apogees).push(Math.floor(t0));}
        else if (d0 * d1 < 0) { // sign change in derivative
            while (t1 - t0 > 6e5) { // 6e5 refers to 10-minute window
                const tAvg = (t0 + t1) / 2;
                const dAvg = moonDistanceDeriv(tAvg, true);
                if (d1 * dAvg <= 0) {t0 = tAvg; d0 = dAvg;}
                else {t1 = tAvg; d1 = dAvg;}
            }
            const apsis = Math.floor(t0 + (d0 / (d0 - d1)) * (t1 - t0));
            (d1 > 0 ? perigees : apogees).push(apsis);
        }
    }
    return {apogees: apogees, perigees: perigees};
}

/** Returns all moonrises, moonsets, meridian passings, full moons, new moons, and first/last quarters on a given day. */
export function allMoonEvents(lat: number, long: number, start: number, end: number): SEvent[] {
    const maxMin = moonMaxMin(lat, long, start, end);
    const rise = moonrise(lat, long, maxMin);
    const set = moonset(lat, long, maxMin);
    const meridian = moonMeridianPassing(lat, long, start, end);
    const phase = moonPhase(lat, long, start, end);
    const events = [...rise, ...set, ...meridian];
    if (phase !== null) {events.push(phase);}
    events.sort((a, b) => a.unix - b.unix);
    return events;
}

/** Variant of allMoonEvents, but with a DateTime object as input rather than start/end Unix times. 
 * Note that here, latitude and longitude are given in degrees, not radians, as they are treated as geographical coordinates
 * rather than angles. They are converted to radians for calculation. */
export function moonEventsDay(lat: number, long: number, date: DateTime): SEvent[] {
    const dayStart = mf.ms(date.startOf("day"));
    const dayEnd = mf.ms(date.startOf("day").plus({days: 1}));
    return allMoonEvents(lat * degToRad, long * degToRad, dayStart, dayEnd);
}

/** Returns the intervals during the day in which the moon is above the horizon. */
export function moonIntervals(lat: number, long: number, dayStart: number, events: SEvent[], timeZone: TimeChange[]): number[][] {
    function timeOfDayS(event: SEvent, zone: TimeChange[]) {return Math.floor(getTimeOfDay(event.unix,zone)/1000);}

    const newMoonEvents: SEvent[] = [];
    for (const event of events) {
        if (event.type === "Moonrise" || event.type === "Moonset") {newMoonEvents.push(event);}
    }
    if (newMoonEvents.length === 0) { // if no moonrise or moonset
        if (moonPosition(lat, long, dayStart)[0] >= HORIZON) {return [[0, 86400]];} // up all day
        else {return [];} // down all day
    }
    else {
        const intervals = (newMoonEvents[0].type === "Moonrise") ? [] : [[0, 86400]];
        for (let i=0; i<newMoonEvents.length; i++) {
            if (newMoonEvents[i].type === "Moonrise") {intervals.push([timeOfDayS(newMoonEvents[i], timeZone), 86400]);} 
            else {intervals[intervals.length-1][1] = timeOfDayS(newMoonEvents[i], timeZone);}
        }
        return intervals;
    }
}

/**
 * Generates table of lunar events for an entire year, including moonrise, moonset, meridian passing, phases, perigees and
 * apogees, as well as a time zone lookup table for the given year.
 * @param lat Latitude in degrees
 * @param long Longitude in degrees
 * @param year Year
 * @param zone Time zone in IANA format (ex. "America/Los_Angeles")
 * @returns MoonTable object
 */
export function generateMoonTable(lat: number, long: number, year: number, zone: string): MoonTable {
    lat *= degToRad; long *= degToRad;
    const yearStart = DateTime.fromObject({year: year}, {zone: zone});
    const yearEnd = DateTime.fromObject({year: year + 1}, {zone: zone});
    const dayStarts = mf.dayStarts(yearStart, yearEnd);
    const timeZoneTable = timeZoneLookupTable(dayStarts); // time zone lookup table

    const moonEventsYear: SEvent[][] = []; // note: this includes phases, but not perigee/apogee
    const intervals: number[][][] = [];
    for (let i=0; i<dayStarts.length-1; i++) {
        const start = mf.ms(dayStarts[i]), end = mf.ms(dayStarts[i+1]);
        const evts = allMoonEvents(lat, long, start, end);
        const ints = moonIntervals(lat, long, start, evts, timeZoneTable);
        moonEventsYear.push(evts);
        intervals.push(ints);
    }

    // In the phase table, the type is 0 for new moon, 1 for first quarter, 2 for full moon, 3 for last quarter
    const phases: {unix: number, date: number, type: number}[] = [];
    for (const events of moonEventsYear) {
        for (const event of events) {
            const type = (event.type === "New Moon" ? 0 : event.type === "First Quarter" ? 1 :
                event.type === "Full Moon" ? 2 : event.type === "Last Quarter" ? 3 : -1);
            if (type >= 0) {
                const date = DateTime.fromMillis(event.unix, {zone: zone}).diff(yearStart, ["days", "hours"]).toObject().days!;
                phases.push({unix: event.unix, date: date, type: type});
            }
        }
    }

    const apsides = moonApsides(mf.ms(yearStart), mf.ms(yearEnd));

    return {
        lunarEvents: moonEventsYear,
        intervals: intervals,
        phases: phases,
        perigees: apsides.perigees,
        apogees: apsides.apogees,
        timeZoneTable: timeZoneTable
    }
}