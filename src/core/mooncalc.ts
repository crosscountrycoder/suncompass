/** Formulas derived from "Astronomical Algorithms" by Jean Meeus. */

import * as mf from "./mathfuncs.ts";
import * as sc from "./suncalc.ts";
import {degToRad, moonPtl, moonPtld} from "./constants.ts";
import { generateLODProfile, getTimeOfDay, type SEvent, type TimeChange } from "./lookup-tables.ts";
import type { DateTime } from "luxon";

export function moonMeanLongitude(JC: number): number {
    return mf.mod(218.3164591 + 481267.88134236*JC - 0.0013268*JC**2 + JC**3/538841 - JC**4/65194000, 360);
}

export function moonMeanElongation(JC: number): number {
    return mf.mod(297.8502042 + 445267.1115168*JC - 0.00163*JC**2 + JC**3/545868 - JC**4/113065000, 360);
}

export function moonMeanAnomaly(JC: number): number {
    return mf.mod(134.9634114 + 477198.8676313*JC + 0.008997*JC**2 + JC**3/69699 - JC**4/14712000, 360);
}

/** Moon argument of latitude */
export function moonArgLat(JC: number): number {
    return mf.mod(93.2720993 + 483202.0175273*JC - 0.0034029*JC**2 - JC**3/3526000 + JC**4/863310000, 360);
}

/** Sum of all longitude terms in moonPtld (periodic terms for longitude and distance) */
function l(JC: number): number {
    let l = 0;
    const D = moonMeanElongation(JC);
    const M = sc.meanSunAnomaly(JC);
    const Mp = moonMeanAnomaly(JC);
    const F = moonArgLat(JC);
    const E = 1 - 0.002516*JC - 7.4e-6*JC**2;
    for (let i=0; i<moonPtld.length; i++) {
        const curRow = moonPtld[i];
        let curSum = curRow[4] * Math.sin((curRow[0]*D + curRow[1]*M + curRow[2]*Mp + curRow[3]*F)*degToRad);
        curSum *= (E ** Math.abs(curRow[1])); // Multiply terms which contain -M or M by E, and terms which contain 2M or -2M by E^2.
        l += curSum;
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
    for (let i=0; i<moonPtld.length; i++) {
        const curRow = moonPtld[i];
        let curSum = curRow[5] * Math.cos((curRow[0]*D + curRow[1]*M + curRow[2]*Mp + curRow[3]*F)*degToRad);
        curSum *= (E ** Math.abs(curRow[1])); // Multiply terms which contain -M or M by E, and terms which contain 2M or -2M by E^2.
        r += curSum;
    }
    return r;
}

/** Sum of all latitude terms in moonPtl (periodic terms for latitude) */
function b(JC: number): number {
    let b = 0;
    const D = moonMeanElongation(JC);
    const M = sc.meanSunAnomaly(JC);
    const Mp = moonMeanAnomaly(JC);
    const F = moonArgLat(JC);
    const E = 1 - 0.002516*JC - 7.4e-6*JC**2;
    for (let i=0; i<moonPtl.length; i++) {
        const curRow = moonPtl[i];
        let curSum = curRow[4] * Math.sin((curRow[0]*D + curRow[1]*M + curRow[2]*Mp + curRow[3]*F)*degToRad);
        curSum *= (E ** Math.abs(curRow[1])); // Multiply terms which contain -M or M by E, and terms which contain 2M or -2M by E^2.
        b += curSum;
    }
    return b;
}

function a(JC: number): number[] {
    return [119.75 + 131.849*JC, 53.09 + 479264.29*JC, 313.45 + 481266.484*JC];
}

/** Variations in longitude due to the actions of Venus, Jupiter, and the flattening of Earth. */
function deltaL(JC: number): number {
    const [a1, a2, a3] = a(JC);
    return 3958*Math.sin(a1*degToRad) + 1962*Math.sin((moonMeanLongitude(JC)-moonArgLat(JC))*degToRad) + 318*Math.sin(a2*degToRad);
}

/** Variations in latitude due to the actions of Venus, Jupiter, and the flattening of Earth. */
function deltaB(JC: number): number {
    let [a1, a2, a3] = a(JC);
    a1 *= degToRad; a2 *= degToRad; a3 *= degToRad;
    const meanLong = moonMeanLongitude(JC) * degToRad;
    const meanAnomaly = moonMeanAnomaly(JC) * degToRad;
    const argLat = moonArgLat(JC) * degToRad;
    return -2235*Math.sin(meanLong) + 382*Math.sin(a3) + 175*Math.sin(a1-argLat) +
    175*Math.sin(a1+argLat) + 127*Math.sin(meanLong-meanAnomaly) - 115*Math.sin(meanLong+meanAnomaly);
}

/** Returns the ecliptic latitude and longitude of the moon. Return value is an array: [latitude, longitude], 
 * both measured in degrees. */
export function moonLatLong(date: number, unix = false): number[] {
    if (!unix) {
        let long = moonMeanLongitude(date) + (l(date)+deltaL(date))/1e6 + sc.longNutation(date);
        let lat = (b(date) + deltaB(date))/1e6;
        lat = mf.clamp(lat, -90, 90);
        long = mf.mod(long, 360);
        return [lat, long];
    }
    else {return moonLatLong(mf.jCentury(date));}
}

/** Distance from center of earth to center of moon, in kilometers. */
export function moonDistance(date: number, unix = false) {
    if (!unix) {return 385000.56 + r(date)/1000;}
    else {return moonDistance(mf.jCentury(date));}
}

/** Returns the rectangular coordinates [x, y, z] in Earth-centered, Earth-fixed coordinates (ECEF) in kilometers. */
export function moonEcef(unix: number): number[] {
    let [eLat, eLong] = moonLatLong(unix, true);
    const ob = sc.obliquity(unix, true) * degToRad;
    eLat *= degToRad; eLong *= degToRad;
    const [sinB,cosB,sinL,cosL,sinE,cosE] = [Math.sin(eLat),Math.cos(eLat),Math.sin(eLong),Math.cos(eLong),Math.sin(ob),Math.cos(ob)];
    
    // xq, yq, zq are geocentric equatorial (Earth-centered inertial) coordinates
    const dist = moonDistance(unix, true);
    const xeci = dist * cosB * cosL;
    const yeci = dist * (cosB * sinL * cosE - sinB * sinE);
    const zeci = dist * (cosB * sinL * sinE + sinB * cosE);

    // convert to ECEF coordinates
    const rectCoords = mf.rotateZ(xeci, yeci, zeci, -sc.gast(unix));
    return rectCoords;
}

/** Returns the sublunar point [latitude, longitude], where the moon is directly overhead. */
export function sublunarPoint(unix: number): number[] {
    const [xecef, yecef, zecef] = moonEcef(unix);
    const r = Math.hypot(xecef, yecef, zecef);
    const [ux, uy, uz] = [xecef / r, yecef / r, zecef / r];
    const [lat, lon] = [Math.asin(mf.clamp(uz)), Math.atan2(uy, ux)]; // radians
    return [lat / degToRad, mf.mod(lon / degToRad + 180, 360) - 180];
}

/** Returns the moon's position: [elevation, azimuth] in degrees. Optionally, the observer's ECEF can be specified in order
 * to avoid repeatedly computing it.
*/
export function moonPosition(lat: number, long: number, unix: number): number[] {
    return mf.elevAzimuth(lat, long, mf.latLongEcef(lat, long), moonEcef(unix));
}

/**
 * The time at which the moon crosses the observer's meridian, analogous to solar noon for the sun.
 * @param lat Observer's latitude in degrees
 * @param long Observer's longitude in degrees
 * @param startUnix Start of day in Unix milliseconds
 * @param endUnix End of day in Unix milliseconds
 * @returns An SEvent[] object, containing the Unix time of the meridian passing, the type ("Meridian Passing"), and the
 * elevation and azimuth of the moon as seen from the observer.
 */
export function moonMeridianPassing(lat: number, long: number, startUnix: number, endUnix: number):
SEvent[] {
    const long0 = sublunarPoint(startUnix)[1];
    const long1 = sublunarPoint(endUnix)[1];
    const long0Adjusted = mf.mod(long0 - long, 360);
    const long1Adjusted = mf.mod(long1 - long, 360) - 360;
    const diff = long0Adjusted - long1Adjusted;
    if (diff > 540) { // no meridian passing
        return [];
    } else if (diff < 180) { 
        // 2 meridian passings (possible on DST fall-back days)
        // Example is at 116.4°W in California on 2033-11-04, or 149°W in Alaska on 2063-11-04.
        const lunarPassageRate = (diff+360)/(endUnix-startUnix); // change in sublunar point longitude per millisecond
        let lunarDiff0 = long0Adjusted; // difference between sublunar point & observer longitude at time startUnix
        let lunarNoon0 = startUnix + lunarDiff0 / lunarPassageRate;
        lunarDiff0 = mf.mod(sublunarPoint(lunarNoon0)[1] - long + 180, 360) - 180;
        lunarNoon0 += lunarDiff0 / lunarPassageRate;
        lunarNoon0 = Math.floor(mf.clamp(lunarNoon0, startUnix, endUnix-1));
        const [e0, a0] = moonPosition(lat, long, lunarNoon0);
        const lunarNoon0Obj = {unix: lunarNoon0, type: "Meridian Passing", elev: e0, azimuth: a0};

        let lunarDiff1 = long1Adjusted + 360;
        let lunarNoon1 = endUnix - lunarDiff1 / lunarPassageRate;
        lunarDiff1 = mf.mod(sublunarPoint(lunarNoon1)[1] - long + 180, 360) - 180;
        lunarNoon1 += lunarDiff1 / lunarPassageRate;
        lunarNoon1 = Math.floor(mf.clamp(lunarNoon1, startUnix, endUnix-1));
        const [e1, a1] = moonPosition(lat, long, lunarNoon1);
        const lunarNoon1Obj = {unix: lunarNoon1, type: "Meridian Passing", elev: e1, azimuth: a1};
        return [lunarNoon0Obj, lunarNoon1Obj];
    } else { // 1 meridian passing
        const lunarPassageRate = diff / (endUnix - startUnix);
        let lunarDiff = long0Adjusted;
        let lunarNoon = startUnix + lunarDiff / lunarPassageRate;
        lunarDiff = mf.mod(sublunarPoint(lunarNoon)[1] - long + 180, 360) - 180;
        lunarNoon += lunarDiff / lunarPassageRate;
        lunarNoon = Math.floor(mf.clamp(lunarNoon, startUnix, endUnix-1));
        const [e1, a1] = moonPosition(lat, long, lunarNoon);
        const lunarNoonObj = {unix: lunarNoon, type: "Meridian Passing", elev: e1, azimuth: a1};
        return [lunarNoonObj];
    }
}

/** Derivative of the moon's position in degrees per second. */
export function moonDerivative(lat: number, long: number, unix: number): number {
    const t0 = unix - 500, t1 = unix + 500;
    const e0 = moonPosition(lat, long, t0)[0], e1 = moonPosition(lat, long, t1)[0];
    return e1 - e0;
}

/**
 * Calculates the extrema of the moon's elevation at a given location - i.e. the points at which the moon's elevation reaches
 * a relative maximum or minimum. Lunar equivalent of maxAndMin() for the sun.
 * @param lat Observer's latitude, degrees
 * @param long Observer's longitude, degrees
 * @param start Start of day, Unix milliseconds
 * @param end End of day, Unix milliseconds
 * @returns An array of Unix times, including the start of the day, all extrema in lunar elevation, and the start of the next day.
 */
export function moonMaxMin(lat: number, long: number, start: number, end: number): number[] {
    const times = [start];
    const intervals = [start,start+4*3.6e6,start+8*3.6e6,start+12*3.6e6,start+16*3.6e6,start+20*3.6e6,end];
    for (let i=0; i<intervals.length; i++) {
        let t0 = intervals[i], t1 = intervals[i+1];
        const d0 = moonDerivative(lat, long, t0), d1 = moonDerivative(lat, long, t1);
        if (d0 >= 0 && d1 < 0) { // maximum (i.e. solar noon, or summer solstice at pole)
            while (t1 - t0 > 1) {
                const tAvg = Math.floor((t0+t1)/2);
                const dAvg = moonDerivative(lat, long, tAvg);
                if (dAvg >= 0) {t0 = tAvg;}
                else {t1 = tAvg;}
            }
            times.push(t0);
        }
        else if (d0 <= 0 && d1 > 0) { // minimum (i.e. solar midnight, or winter solstice at pole)
            while (t1 - t0 > 1) {
                const tAvg = Math.floor((t0+t1)/2);
                const dAvg = moonDerivative(lat, long, tAvg);
                if (dAvg <= 0) {t0 = tAvg;}
                else {t1 = tAvg;}
            }
            times.push(t0);
        }
    }
    times.push(end);
    return times;
}

/**
 * Returns the time(s) of moonrise on a particular day.
 * @param lat Latitude of observer
 * @param long Longitude of observer
 * @param maxMin Results of moonMaxMin() for the given day and location
 * @param angle The threshold elevation angle (defaults to -5/6°)
 * @returns An array of SEvent objects, with the time, elevation, and azimuth of the moon, and a type tag.
 */
export function moonrise(lat: number, long: number, maxMin: number[], angle: number = -5/6): SEvent[] {
    const riseTimes = [];
    for (let i=0; i<maxMin.length-1; i++) {
        let t0 = maxMin[i], t1 = maxMin[i+1];
        let [e0, a0] = moonPosition(lat, long, t0);
        let [e1, a1] = moonPosition(lat, long, t1);
        if (e0 <= angle && e1 >= angle) {
            while (t1 - t0 > 1) {
                const tAvg = Math.floor((t0+t1)/2);
                const [eAvg, aAvg] = moonPosition(lat, long, tAvg);
                if (eAvg <= angle) {t0 = tAvg; e0 = eAvg; a0 = aAvg;}
                else {t1 = tAvg; e1 = eAvg; a1 = aAvg;}
            }
            riseTimes.push({unix: t0, type: "Moonrise", elev: e0, azimuth: a0});
        }
    }
    return riseTimes;
}

/**
 * Returns the time(s) of moonset on a particular day.
 * @param lat Latitude of observer
 * @param long Longitude of observer
 * @param maxMin Results of moonMaxMin() for the given day and location
 * @param angle The threshold elevation angle (defaults to -5/6°)
 * @returns An array of SEvent objects, with the time, elevation, and azimuth of the moon, and a type tag.
 */
export function moonset(lat: number, long: number, maxMin: number[], angle: number = -5/6): SEvent[] {
    const setTimes = [];
    for (let i=0; i<maxMin.length-1; i++) {
        let t0 = maxMin[i], t1 = maxMin[i+1];
        let [e0, a0] = moonPosition(lat, long, t0);
        let [e1, a1] = moonPosition(lat, long, t1);
        if (e0 >= angle && e1 <= angle) {
            while (t1 - t0 > 1) {
                const tAvg = Math.floor((t0+t1)/2);
                const [eAvg, aAvg] = moonPosition(lat, long, tAvg);
                if (eAvg >= angle) {t0 = tAvg; e0 = eAvg; a0 = aAvg;}
                else {t1 = tAvg; e1 = eAvg; a1 = aAvg;}
            }
            setTimes.push({unix: t0, type: "Moonset", elev: e0, azimuth: a0});
        }
    }
    return setTimes;
}

/** Returns the phase angle of the moon given Unix time, in RADIANS (not degrees).
 * The value returned is the same as the phase angle i in Astronomical Algorithms, chapter 46. t is the angle between the
 * sun-moon and earth-moon vectors.
*/
function phaseAngle(unix: number): number {
    const lodProfile = generateLODProfile(unix);
    const [subsolarLat, subsolarLong] = sc.subsolarPoint(lodProfile);
    const ecefSun = mf.toEcef(subsolarLat, subsolarLong, lodProfile.distance);
    const ecefMoon = moonEcef(unix);
    const sunVector = [ecefSun[0]-ecefMoon[0], ecefSun[1]-ecefMoon[1], ecefSun[2]-ecefMoon[2]];
    const moonVector = [-ecefMoon[0], -ecefMoon[1], -ecefMoon[2]];

    const dot = sunVector[0] * moonVector[0] + sunVector[1] * moonVector[1] + sunVector[2] * moonVector[2];
    const normSun = Math.hypot(sunVector[0], sunVector[1], sunVector[2]);
    const normMoon = Math.hypot(moonVector[0], moonVector[1], moonVector[2]);
    const i = Math.acos(mf.clamp(dot / (normSun * normMoon)));
    return i;
}

/** Returns the fraction of the Moon's disk that is illuminated. This ranges from 0 (new moon) to 1 (full moon), and is 0.5 at both
 * first and last quarters. */
export function illumination(unix: number): number {
    return (1 + Math.cos(phaseAngle(unix))) / 2;
}

/** The difference between the ecliptic longitudes of the sun and moon, in degrees, normalized to the range [0, 360).
 * A value of 0 indicates a new moon, 90 is a first quarter, 180 is a half moon and 270 is a last quarter.
 */
export function moonSunLongDiff(unix: number): number {
    return mf.mod(moonLatLong(unix, true)[1] - sc.sunTrueLong(unix, true), 360);
}

/**
 * Returns the moon phases on a given day.
 * @param lat Latitude of observer
 * @param long Longitude of observer
 * @param start Start of the day in Unix milliseconds
 * @param end End of the day in Unix milliseconds
 * @returns If there is a new moon, full moon, or first/last quarter that day, returns an SEvent object with the Unix time of the
 * event, along with a type tag ("New Moon", "First Quarter", "Full Moon", or "Last Quarter") and the moon's elevation and azimuth.
 * Otherwise, returns null.
 */
export function moonPhase(lat: number, long: number, start: number, end: number): SEvent | null {
    const diff0 = moonSunLongDiff(start);
    const diff1 = moonSunLongDiff(end);
    if (diff1 < diff0) { // new moon
        let t0 = start, t1 = end;
        while (t1 - t0 > 1) {
            const avg = Math.floor((t0+t1)/2);
            (moonSunLongDiff(avg) > 180) ? (t0 = avg) : (t1 = avg);
        }
        const [e, a] = moonPosition(lat, long, t0);
        return {unix: t0, type: "New Moon", elev: e, azimuth: a};
    }
    else if (diff0 <= 90 && diff1 >= 90) { // first quarter
        let t0 = start, t1 = end;
        while (t1 - t0 > 1) {
            const avg = Math.floor((t0+t1)/2);
            (moonSunLongDiff(avg) <= 90) ? (t0 = avg) : (t1 = avg);
        }
        const [e, a] = moonPosition(lat, long, t0);
        return {unix: t0, type: "First Quarter", elev: e, azimuth: a};
    }
    else if (diff0 <= 180 && diff1 >= 180) { // full moon
        let t0 = start, t1 = end;
        while (t1 - t0 > 1) {
            const avg = Math.floor((t0+t1)/2);
            (moonSunLongDiff(avg) <= 180) ? (t0 = avg) : (t1 = avg);
        }
        const [e, a] = moonPosition(lat, long, t0);
        return {unix: t0, type: "Full Moon", elev: e, azimuth: a};
    }
    else if (diff0 <= 270 && diff1 >= 270) { // last quarter
        let t0 = start, t1 = end;
        while (t1 - t0 > 1) {
            const avg = Math.floor((t0+t1)/2);
            (moonSunLongDiff(avg) <= 270) ? (t0 = avg) : (t1 = avg);
        }
        const [e, a] = moonPosition(lat, long, t0);
        return {unix: t0, type: "Last Quarter", elev: e, azimuth: a};
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
    const diff0 = moonSunLongDiff(start);
    const diff1 = moonSunLongDiff(end);
    if (diff0 <= 90) {return (diff1 >= 90) ? "First Quarter" : "Waxing Crescent";}
    else if (diff0 <= 180) {return (diff1 >= 180) ? "Full Moon" : "Waxing Gibbous";}
    else if (diff0 <= 270) {return (diff1 >= 270) ? "Last Quarter" : "Waning Gibbous";}
    else {return (diff1 < diff0) ? "New Moon" : "Waning Crescent";}
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

/** Variant of allMoonEvents, but with a DateTime object as input rather than start/end Unix times. */
export function moonEventsDay(lat: number, long: number, date: DateTime): SEvent[] {
    const dayStart = mf.ms(date.startOf("day"));
    const dayEnd = mf.ms(date.startOf("day").plus({days: 1}));
    return allMoonEvents(lat, long, dayStart, dayEnd);
}

/** Returns the intervals during the day in which the moon is above the horizon. */
export function moonIntervals(lat: number, long: number, dayStart: number, events: SEvent[], timeZone: TimeChange[]): number[][] {
    const newMoonEvents: SEvent[] = [];
    for (const event of events) {
        if (event.type == "Moonrise" || event.type == "Moonset") {newMoonEvents.push(event);}
    }
    if (newMoonEvents.length == 0) { // if no moonrise or moonset
        if (moonPosition(lat, long, dayStart)[0] >= -5/6) {return [[0, 86400]];} // up all day
        else {return [];} // down all day
    }
    else if (newMoonEvents[0].type == "Moonrise") {
        if (newMoonEvents.length == 1) {
            return [[mf.intDiv(getTimeOfDay(newMoonEvents[0].unix, timeZone),1000), 86400]];
        } else if (newMoonEvents.length == 2) {
            return [[mf.intDiv(getTimeOfDay(newMoonEvents[0].unix, timeZone),1000), 
            mf.intDiv(getTimeOfDay(newMoonEvents[1].unix, timeZone),1000)]];
        } else if (newMoonEvents.length == 3) {
            return [[mf.intDiv(getTimeOfDay(newMoonEvents[0].unix, timeZone),1000), mf.intDiv(getTimeOfDay(newMoonEvents[1].unix, timeZone),1000)],
                [mf.intDiv(getTimeOfDay(newMoonEvents[2].unix, timeZone),1000), 86400]];
        } else {
            return [[mf.intDiv(getTimeOfDay(newMoonEvents[0].unix, timeZone),1000), mf.intDiv(getTimeOfDay(newMoonEvents[1].unix, timeZone),1000)],
                [mf.intDiv(getTimeOfDay(newMoonEvents[2].unix, timeZone),1000), mf.intDiv(getTimeOfDay(newMoonEvents[3].unix, timeZone),1000)]];
        }
    }
    else { // if (newMoonEvents[0].type == "Moonset")
        if (newMoonEvents.length == 1) {
            return [[0, mf.intDiv(getTimeOfDay(newMoonEvents[0].unix, timeZone),1000)]];
        } else if (newMoonEvents.length == 2) {
            return [[0, mf.intDiv(getTimeOfDay(newMoonEvents[0].unix, timeZone),1000)], 
                [mf.intDiv(getTimeOfDay(newMoonEvents[1].unix, timeZone),1000), 86400]];
        } else if (newMoonEvents.length == 3) {
            return [[0, mf.intDiv(getTimeOfDay(newMoonEvents[0].unix, timeZone),1000)], 
                [mf.intDiv(getTimeOfDay(newMoonEvents[1].unix, timeZone),1000), mf.intDiv(getTimeOfDay(newMoonEvents[2].unix, timeZone),1000)]];
        } else {
            return [[0, getTimeOfDay(newMoonEvents[0].unix, timeZone)], 
                [mf.intDiv(getTimeOfDay(newMoonEvents[1].unix, timeZone),1000), mf.intDiv(getTimeOfDay(newMoonEvents[2].unix, timeZone),1000)],
                [mf.intDiv(getTimeOfDay(newMoonEvents[3].unix, timeZone),1000), 86400]];
        }
    }
}