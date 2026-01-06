/* 
The formulas for earth-sun distance, axial tilt, declination, equation of time, and atmospheric refraction of sunlight are 
borrowed from the book "Astronomical Algorithms" by Jean Meeus. The refraction formula is modified slightly to ensure continuity 
when the sun is below the horizon. The formula for solar ecliptic longitude is from the book "Planetary Programs and Tables from 
-4000 to +2800" by Pierre Bretagnon and Jean-Louis Simon.

The subsolar point is calculated using the declination and equation of time, treating UTC as equivalent to UT1 (mean solar time 
at 0 degrees longitude), the two are always within 0.9 seconds of each other. Solar noon and midnight are calculated using the 
equation of time. The position of the sun is calculated from the subsolar point through spherical trigonometry.

For the purposes of this calculator, daytime is defined as the period in which the solar elevation angle is greater than -50 
arcminutes (-5/6 degrees, or -0.1454 radians). This accounts for the sun's angular radius in the sky and refraction of sunlight.
It corresponds to the definition used by most sources, including NOAA's solar calculator.

This site uses the Luxon library to deal with date/time computations. Luxon is used to simplify computation when dealing with
durations, conversion between different time zones, and complexities such as daylight saving time. The geo-tz library is used 
to find the time zone of a geographic coordinate.
*/

import * as mf from "./mathfuncs.ts";
import {degToRad, sunPeriodicTerms, DAY_LENGTH, BSEARCH_GAP, TAU, HORIZON, CIVIL_TWILIGHT, NAUTICAL_TWILIGHT, ASTRO_TWILIGHT} from "./constants.ts";
import {generateLODProfile, estimateLOD, getTimeOfDay, timeZoneLookupTable, longDistLookupTable} from "./lookup-tables.ts";
import type {TimeChange, LODProfile, SEvent} from "./lookup-tables.ts";
import {DateTime} from "luxon";

export type SeasonEvents = {marEquinox: DateTime; junSolstice: DateTime; sepEquinox: DateTime; decSolstice: DateTime;};

export type SunTable = {
    solarEvents: SEvent[][];
    dayLengths: number[];
    solstices: number[];
    solsticeDates: number[];
    perihelion: number[];
    aphelion: number[];
    timeZoneTable: TimeChange[];
}

/** Sun's mean longitude in radians according to formula 27.2 in Astronomical Algorithms. 
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
*/
export function sunMeanLong(date: number, unix = false): number {
    if (!unix) { // if date is specified as a Julian century
        return mf.polymod(date, [4.895063111, 628.3319667476, 5.291887e-6, 3.495482270e-10, -1.14081e-10, -8.77932e-14], TAU);
    }
    else {return sunMeanLong(mf.jCentury(date));}
}

/** Sun's geometric longitude in radians, i.e. longitude excluding aberration and nutation. 
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
*/
export function sunGeomLong(date: number, unix = false): number {
    if (!unix) {
        let long = 4.9353929 + 628.33196168*date;
        for (const curRow of sunPeriodicTerms) {
            long += curRow[0] * Math.sin(curRow[2]+curRow[3]*date);
        }
        return mf.mod(long, TAU);
    }
    else {return sunGeomLong(mf.jCentury(date));}
}

/** Formula 45.3, in page 308 of Astronomical Algorithms. The value returned is in radians */
export function meanSunAnomaly(JC: number): number {
    return mf.polymod(JC, [6.240060127, 628.3019551672, -2.68083e-6, 7.1267e-10], TAU);
}

/** Distance from sun to earth in kilometers. 
 * @param date The timestamp. Can be specified as a Unix timestamp, or a Julian century.
 * @param unix If date is a number, this parameter specifies whether it's measured in Unix milliseconds or Julian centuries.
*/
export function sunDistance(date: number, unix = false): number {
    if (!unix) {
        let dist = 1.0001026; // in astronomical units
        for (const curRow of sunPeriodicTerms) {
            dist += curRow[1] * Math.cos(curRow[2]+curRow[3]*date);
        }
        return 149597870.7 * dist; // convert to kilometers
    }
    else {return sunDistance(mf.jCentury(date));}
}

/** Derivative of sun-earth distance with respect to time, in km per second. Used to calculate perihelion and aphelion. */
export function sunDistanceDerivative(unix: number): number {
    const date = mf.jCentury(unix);
    let deriv = 0;
    for (const curRow of sunPeriodicTerms) { // derivative of cosine is negative sine
        deriv -= curRow[1] * curRow[3] * Math.sin(curRow[2] + curRow[3]*date);
    }
    return 149597870.7 * deriv / (36525 * 86400); // convert AU/century to km/second
}

/**
 * Calculates the sun's apparent ecliptic longitude in radians to within 1.6e-5 rad for years 0-3000. This value is 0 at the
 * March equinox, 0.25*TAU at the June solstice, 0.5*TAU at the September equinox and 0.75*TAU at the December solstice.
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
 */
export function sunTrueLong(date: number, unix = false): number {
    if (!unix) {
        const geoLong = sunGeomLong(date);
        const aberration = 1.7e-6 * Math.cos(3.1 + 628.3014 * date) - 9.93e-5;
        const nutation = -8.34e-5*Math.sin(2.18-33.757*date+3.6e-5*date**2)-6.4e-6*Math.sin(3.51+1256.6639*date+1e-5*date**2);
        return mf.mod(geoLong + aberration + nutation, TAU);
    }
    else {return sunTrueLong(mf.jCentury(date));}
}

/** Nutation in obliquity in radians 
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
*/
export function obNutation(date: number, unix = false) {
    if (!unix) {
        const L = mf.mod(4.895064 + 628.3319663*date, TAU);
        const Lprime = mf.mod(3.810342 + 8399.709113*date, TAU);
        const omega = mf.polymod(date, [2.1824386, -33.75704594, 3.61423e-5, 3.87851e-8], TAU);
        return 4.5e-5*Math.cos(omega) + 2.8e-6*Math.cos(2*L) + 4.8e-7*Math.cos(2*Lprime) - 4.4e-7*Math.cos(2*omega);
    }
    else {return obNutation(mf.jCentury(date));}
}

/** Nutation in longitude in radians 
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
*/
export function longNutation(date: number, unix = false) {
    if (!unix) {
        const L = mf.mod(4.895064 + 628.3319663*date, TAU);
        const Lprime = mf.mod(3.810342 + 8399.709113*date, TAU);
        const omega = mf.polymod(date, [2.1824386, -33.75704594, 3.61423e-5, 3.87851e-8], TAU);
        return -8.34e-5*Math.sin(omega) - 6.40e-6*Math.sin(2*L) - 1.12e-6*Math.sin(2*Lprime) + 1.02e-6*Math.sin(2*omega);
    }
    else {return longNutation(mf.jCentury(date));}
}

/**
 * Returns the obliquity of the ecliptic, or equivalently Earth's axial tilt, in radians.
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
 */
export function obliquity(date: number, unix = false): number {
    if (!unix) {
        const meanObliquity = mf.polynomial(date, [0.409092804, -2.26966e-4, -2.8604e-9, 8.7897e-9]);
        return meanObliquity + obNutation(date);
    }
    else {return obliquity(mf.jCentury(date));}
}

/**
 * Calculates the sun's right ascension in radians. To convert to hours, multiply by (24/TAU).
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
 */
export function sunRA(date: number, unix = false): number {
    if (!unix) {
        const long = sunTrueLong(date);
        const ob = obliquity(date);
        const ra = Math.atan2(Math.sin(long)*Math.cos(ob), Math.cos(long));
        return mf.mod(ra, TAU);
    }
    else {return sunRA(mf.jCentury(date));}
}

/**
 * Equation of time in radians. To convert to minutes, multiply by (1440 / TAU).
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
 */
export function equationOfTime(date: number, unix = false): number { 
    if (!unix) {
        const eot = sunMeanLong(date) - 9.98032e-5 - sunRA(date) + longNutation(date) * Math.cos(obliquity(date));
        return mf.mod(eot + TAU/2, TAU) - TAU/2; // reduce to range [-tau/2, tau/2)
    }
    else {return equationOfTime(mf.jCentury(date));}
}

/** Gives the value of Greenwich apparent sidereal time (GAST) in radians, from equation 11.4 in Astronomical Algorithms.
 * The value returned is in the range 0 <= x < TAU.
 * @param unix Unix timestamp in milliseconds.
 */
export function gast(unix: number): number {
    const JC = (mf.jdUTC(unix) - 2451545) / 36525; // Julian century, but with UTC rather than TT
    const gmst = mf.polynomial(JC, [4.8949612127, 230121.67531543, 6.77071e-6, -4.50873e-10]);
    const correction = longNutation(unix, true) * Math.cos(obliquity(unix, true));
    return mf.mod(gmst + correction, TAU);
}

/**
 * Returns the sun's hour angle in radians, or the observer's longitude minus the longitude of the subsolar point. A value
 * of ±TAU/2 indicates solar midnight and 0 indicates solar noon. To convert to apparent solar time in hours, multiply by 
 * (24/TAU) and add 12.
 * @param long Observer's longitude in radians.
 * @param unix Unix timestamp in milliseconds.
 */
export function sunHourAngle(long: number, unix: number): number {
    const timeEq = equationOfTime(unix, true);
    const mst = TAU * (unix + DAY_LENGTH/2) / DAY_LENGTH + long; // mean solar time
    return mf.mod(timeEq + mst + TAU/2, TAU) - TAU/2; // apparent solar time
}

/**
 * Returns the sun's declination in radians. This is the latitude of the subsolar point.
 * @param long Sun's true ecliptic longitude, in radians
 * @param obliquity Obliquity of the ecliptic (i.e. Earth's axial tilt)
 */
export function declination(long: number, obliquity: number): number {
    return Math.asin(Math.sin(obliquity)*Math.sin(long));
}

/**
 * Returns the subsolar point, or location on Earth at which the sun is directly overhead.
 * @param date Longitude, obliquity, distance (LOD) profile for the given moment.
 * @param degrees Whether to return values in degrees (true) or radians (false). Defaults to false.
 * @returns [latitude, longitude] of subsolar point.
 */
export function subsolarPoint(lod?: LODProfile, degrees = false): number[] {
    if (lod === undefined) {lod = generateLODProfile(Date.now());}
    const lat = declination(lod.longitude, lod.obliquity), long = -sunHourAngle(0, lod.unix);
    if (degrees) {return [mf.clamp(lat/degToRad,-90,90), mf.mod(long/degToRad+180,360)-180];}
    else {return [lat, long];}
}

/** Returns the sun's ECEF coordinates, in kilometers, at the given Unix timestamp in milliseconds. */
export function sunEcef(unix: number): number[] {
    const lod = generateLODProfile(unix);
    const [sunLat, sunLong] = subsolarPoint(lod);
    return mf.toEcef(sunLat, sunLong, lod.distance);
}

/**
 * Returns sun position given latitude, longitude, and DateTime.
 * @param lat Latitude in radians
 * @param long Longitude in radians
 * @param lod Longitude, obliquity, distance (LOD) profile.
 * @returns Array: [elevation, azimuth]. Elevation is in radians above horizon, azimuth is radians clockwise from north
 * Solar elevation is not refracted. To adjust elevation angle for atmospheric refraction, use mf.refract(sunPosition[0])
 */
export function sunPosition(lat: number, long: number, lod: LODProfile): number[] {
    const [sunLat, sunLong] = subsolarPoint(lod);
    return mf.elevAzimuth(lat, long, mf.latLongEcef(lat, long), mf.toEcef(sunLat, sunLong, lod.distance));
}

/**
 * Returns the time(s) of solar noon, along with the sun's position at solar noon.
 * @param lat Latitude in radians.
 * @param long Longitude in radians.
 * @param start Start (as LOD profile)
 * @param end End (as LOD profile)
 * @returns Time(s) of solar noon, along with the sun's position at solar noon.
 */
export function solarNoon(lat: number, long: number, start: LODProfile, end: LODProfile): SEvent[] {
    const solarNoonTimes = mf.meridianPassings(long, start.unix, end.unix, 0, sunHourAngle);
    const events: SEvent[] = [];
    for (const time of solarNoonTimes) {
        const [e, a] = sunPosition(lat, long, estimateLOD(time, start, end));
        events.push({unix: time, type: "Solar Noon", elev: e, azimuth: a});
    }
    return events;
}

/**
 * Returns the time(s) of solar midnight, along with the sun's position at solar midnight.
 * @param lat Latitude in radians.
 * @param long Longitude in radians.
 * @param start Start (as LOD profile)
 * @param end End (as LOD profile)
 * @returns Time(s) of solar midnight, along with the sun's position at solar midnight.
 */
export function solarMidnight(lat: number, long: number, start: LODProfile, end: LODProfile): SEvent[] {
    const solarMidnightTimes = mf.meridianPassings(long, start.unix, end.unix, TAU/2, sunHourAngle);
    const events: SEvent[] = [];
    for (const time of solarMidnightTimes) {
        const [e, a] = sunPosition(lat, long, estimateLOD(time, start, end));
        events.push({unix: time, type: "Solar Midnight", elev: e, azimuth: a});
    }
    return events;
}

/** Returns the approximate derivative of the solar elevation angle at a particular time, in radians per second. */
export function derivative(lat: number, long: number, unix: number, startLOD: LODProfile, endLOD: LODProfile) {
    const t0 = unix - 500, t1 = unix + 500;
    const LOD0 = estimateLOD(t0, startLOD, endLOD), LOD1 = estimateLOD(t1, startLOD, endLOD);
    return sunPosition(lat, long, LOD1)[0] - sunPosition(lat, long, LOD0)[0];
}

/**
 * Returns an array of Unix times, representing the start of the current day, the times at which the derivative of solar
 * elevation angle is 0, and the start of the next day. This is a helper function for "dawn" and "dusk" functions below.
 * @param lat Latitude in radians
 * @param long Longitude in radians
 * @param start Start of day (LOD profile)
 * @param end End of day (LOD profile)
 * @returns Unix timestamps at which sun reaches relative min or max altitude, along with both the start and end of the day.
 */
export function maxAndMin(lat: number, long: number, start: LODProfile, end: LODProfile): number[] {
    const startTime = start.unix, endTime = end.unix;
    const times = [startTime];
    const interval = (endTime - startTime) / 6;
    for (let t = startTime; t < endTime; t += interval) {
        // use binary search to find the time closest to zero derivative
        let t0 = t, t1 = Math.min(t + interval, endTime);
        let d0 = derivative(lat, long, t0, start, end), d1 = derivative(lat, long, t1, start, end);
        if (d0 === 0) {times.push(t0);}
        else if (d0 * d1 < 0) { // derivative changes sign
            while (t1 - t0 > BSEARCH_GAP) {
                const tAvg = (t0+t1)/2;
                const dAvg = derivative(lat, long, tAvg, start, end);
                if (d1 * dAvg <= 0) {t0 = tAvg; d0 = dAvg;}
                else {t1 = tAvg; d1 = dAvg;}
            }
            const t = t0 + (d0 / (d0 - d1)) * (t1 - t0); // linear interpolation on 5-minute window
            times.push(t);
        }
    }
    times.push(endTime);
    return times;
}

/**
 * Calculates the time in the morning at which the sun's elevation reaches the specified angle in radians.
 * @param lat Latitude in radians
 * @param long Longitude in radians
 * @param angle Solar elevation angle in radians
 * @param type "Sunrise", "Civil Dawn", "Nautical Dawn" or "Astro Dawn"
 * @param maxMin Results of maxAndMin() for given day
 * @param startLOD LOD profile for beginning of day
 * @param endLOD LOD profile for beginning of next day
 * @returns SEvent object, which includes a DateTime, the sun's elevation & azimuth and a tag for the type of dawn/sunrise.
 */
export function dawn(lat: number, long: number, angle: number, type: string, maxMin: number[], 
    startLOD: LODProfile, endLOD: LODProfile
): SEvent[] {
    const dawnTimes = [];
    for (let i=0; i<maxMin.length-1; i++) {
        let t0 = maxMin[i], t1 = maxMin[i+1];
        const lod0 = estimateLOD(t0, startLOD, endLOD), lod1 = estimateLOD(t1, startLOD, endLOD);
        let e0 = sunPosition(lat, long, lod0)[0], e1 = sunPosition(lat, long, lod1)[0];
        if (e0 <= angle && e1 >= angle) {
            while (t1 - t0 > BSEARCH_GAP) {
                const avgLOD = estimateLOD((t0+t1)/2, startLOD, endLOD);
                const eAvg = sunPosition(lat, long, avgLOD)[0];
                if (eAvg <= angle) {t0 = avgLOD.unix; e0 = eAvg;}
                else {t1 = avgLOD.unix; e1 = eAvg;}
            }
            // after reducing to sub-5-minute interval, solve using quadratic interpolation
            const tAvg = (t0 + t1) / 2;
            const avgLOD = estimateLOD(tAvg, startLOD, endLOD);
            const eAvg = sunPosition(lat, long, avgLOD)[0];
            const t = Math.floor(mf.quadraticZero(t0, e0-angle, tAvg, eAvg-angle, t1, e1-angle));
            const [e, a] = sunPosition(lat, long, estimateLOD(t, startLOD, endLOD));
            dawnTimes.push({unix: t, type: type, elev: e, azimuth: a});
        }
    }
    return dawnTimes;
}

/**
 * Calculates the time in the evening at which the sun's elevation reaches the specified angle. Angle should be -5/6 for sunset,
 * -6 for civil twilight, -12 for nautical twilight, and -18 for astronomical twilight.
 * @param lat Latitude in radians
 * @param long Longitude in radians
 * @param angle Solar elevation angle in radians
 * @param type "Sunset", "Civil Dusk", "Nautical Dusk", "Astro Dusk"
 * @param maxMin Results of maxAndMin() for given day
 * @returns SEvent object, which includes a DateTime, the sun's elevation & azimuth and a tag for the type of dusk/sunset.
 */
export function dusk(lat: number, long: number, angle: number, type: string, maxMin: number[],
    startLOD: LODProfile, endLOD: LODProfile
): SEvent[] {
    const duskTimes = [];
    for (let i=0; i<maxMin.length-1; i++) {
        let t0 = maxMin[i], t1 = maxMin[i+1];
        const lod0 = estimateLOD(t0, startLOD, endLOD), lod1 = estimateLOD(t1, startLOD, endLOD);
        let e0 = sunPosition(lat, long, lod0)[0], e1 = sunPosition(lat, long, lod1)[0];
        if (e0 >= angle && e1 <= angle) {
            while (t1 - t0 > BSEARCH_GAP) {
                const avgLOD = estimateLOD((t0+t1)/2, startLOD, endLOD);
                const eAvg = sunPosition(lat, long, avgLOD)[0];
                if (eAvg >= angle) {t0 = avgLOD.unix; e0 = eAvg;}
                else {t1 = avgLOD.unix; e1 = eAvg;}
            }
            // after reducing to sub-5-minute interval, solve using quadratic interpolation
            const tAvg = (t0 + t1) / 2;
            const avgLOD = estimateLOD(tAvg, startLOD, endLOD);
            const eAvg = sunPosition(lat, long, avgLOD)[0];
            const t = Math.floor(mf.quadraticZero(t0, e0-angle, tAvg, eAvg-angle, t1, e1-angle));
            const [e, a] = sunPosition(lat, long, estimateLOD(t, startLOD, endLOD));
            duskTimes.push({unix: t, type: type, elev: e, azimuth: a});
        }
    }
    return duskTimes;
}

export function sunrise(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dawn(lat, long, HORIZON, "Sunrise", maxMin, startLOD, endLOD);
} 
export function sunset(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dusk(lat, long, HORIZON, "Sunset", maxMin, startLOD, endLOD);
}
export function civilDawn(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dawn(lat, long, CIVIL_TWILIGHT, "Civil Dawn", maxMin, startLOD, endLOD);
}
export function civilDusk(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dusk(lat, long, CIVIL_TWILIGHT, "Civil Dusk", maxMin, startLOD, endLOD);
}
export function nauticalDawn(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dawn(lat, long, NAUTICAL_TWILIGHT, "Nautical Dawn", maxMin, startLOD, endLOD);
}
export function nauticalDusk(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dusk(lat, long, NAUTICAL_TWILIGHT, "Nautical Dusk", maxMin, startLOD, endLOD);
}
export function astroDawn(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dawn(lat, long, ASTRO_TWILIGHT, "Astro Dawn", maxMin, startLOD, endLOD);
}
export function astroDusk(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dusk(lat, long, ASTRO_TWILIGHT, "Astro Dusk", maxMin, startLOD, endLOD);
}

/**
 * Returns day length in milliseconds (time from sunrise to sunset). If sunset is after midnight or sunrise is before midnight 
 * (due to time zone complexities and DST), it returns the length of the continuous period of daylight that includes noon local time.
 * @param dayStart The Unix time at the start of the day.
 * @param sunEventsYesterday The value of allSunEvents() for yesterday.
 * @param sunEventsToday The value of allSunEvents() for today.
 * @param sunEventsTomorrow The value of allSunEvents() for tomorrow.
 */
export function dayLength(dayStart: number, sunEventsYesterday: SEvent[], sunEventsToday: SEvent[], sunEventsTomorrow: SEvent[]) {
    const rise: SEvent[] = [], set: SEvent[] = [];
    for (const event of sunEventsToday) {
        if (event.type === "Sunrise") {rise.push(event);}
        else if (event.type === "Sunset") {set.push(event);}
    }
    if (rise.length === 0 && set.length === 0) {
        return (sunEventsToday.length === 0 || sunEventsToday[0].elev >= HORIZON) ? DAY_LENGTH : 0;
    }
    else if (rise.length >= 1 && set.length >= 1 && set.at(-1)!.unix > rise[0].unix) {return set.at(-1)!.unix - rise[0].unix;}
    else {
        const riseY: SEvent[] = [], setT: SEvent[] = [];
        for (const e of sunEventsYesterday) {if (e.type === "Sunrise") {riseY.push(e)};}
        for (const e of sunEventsTomorrow) {if (e.type === "Sunset") {setT.push(e)};}
        const MIDDAY = dayStart + DAY_LENGTH/2; // 12 hours after midnight
        if (setT.length >= 1 && rise.length >= 1 && rise[0].unix < MIDDAY) {return setT[0].unix - rise[0].unix;}
        else if (riseY.length >= 1 && set.length >= 1 && set[0].unix >= MIDDAY) {return set[0].unix - riseY.at(-1)!.unix;}
        else {return -1;} // undefined
    }
}

/** Generate sun events for a given day. Latitude and longitude are in radians. */
export function allSunEvents(lat: number, long: number, start: LODProfile, end: LODProfile): SEvent[] {
    const maxMin = maxAndMin(lat, long, start, end);
    const midnight = solarMidnight(lat, long, start, end);
    const adawn = astroDawn(lat, long, maxMin, start, end);
    const ndawn = nauticalDawn(lat, long, maxMin, start, end);
    const cdawn = civilDawn(lat, long, maxMin, start, end);
    const rise = sunrise(lat, long, maxMin, start, end);
    const noon = solarNoon(lat, long, start, end);
    const set = sunset(lat, long, maxMin, start, end);
    const cdusk = civilDusk(lat, long, maxMin, start, end);
    const ndusk = nauticalDusk(lat, long, maxMin, start, end);
    const adusk = astroDusk(lat, long, maxMin, start, end);
    const events = [...midnight, ...adawn, ...ndawn, ...cdawn, ...rise, ...noon, ...set, ...cdusk, ...ndusk, ...adusk];
    events.sort((a, b) => a.unix - b.unix);
    return events;
}

/** Given a date, returns an array with sunrise, sunset, dawn, dusk, solar noon, and solar midnight times for the given date. 
 * Intended for use for a single day only (for example with the sun.ts script). All other uses, including website uses, should
 * use generateSunTable() instead!
 * @param lat Latitude of observer, degrees
 * @param long Longitude of observer, degrees
 * @param date Luxon DateTime (can be any point within the given day)
*/
export function sunEventsDay(lat: number, long: number, date: DateTime): SEvent[] {
    const start = generateLODProfile(mf.ms(date.startOf("day")));
    const end = generateLODProfile(mf.ms(date.plus({days:1}).startOf("day")));
    return allSunEvents(lat * degToRad, long * degToRad, start, end);
}

/** Calculates the Unix millisecond timestamp of the solstice or equinox in the given month and year.
 * Month must be 3, 6, 9 or 12. */
export function calcSolstEq(year: number, month: 3 | 6 | 9 | 12) {
    let t0 = (year + (month - 1) / 12 - 1970) * 31556952000;
    let t1 = t0 + 2629746000; // approx. 1 month after start
    let l0 = sunTrueLong(t0, true), l1 = sunTrueLong(t1, true);
    if (l0 > l1) {l0 -= TAU;}
    const thresh = (month - 3) * TAU/12;
    while (t1 - t0 > 9e5) { // 15 minutes
        const tAvg = (t0 + t1) / 2;
        let lAvg = sunTrueLong(tAvg, true);
        if (lAvg >= l1) {lAvg -= TAU;}
        if (lAvg <= thresh) {t0 = tAvg; l0 = lAvg;}
        else {t1 = tAvg; l1 = lAvg;}
    }
    const frac = (thresh - l0) / (l1 - l0); // linearly interpolate within 7.5-15 minute window
    const t = Math.floor(t0 + frac * (t1 - t0));
    return t;
}

/** Calculates the sun's apsides (aphelion and perihelion) for a particular year bounded by start and end Unix timestamps.
 * Aphelion is when the earth is furthest from the sun, and perihelion is when it is closest. */
export function sunApsides(start: number, end: number) {
    const interval = (end - start) / 12;
    const window = 3.6e6; // 1 hour
    const aphelion: number[] = [];
    const perihelion: number[] = [];
    for (let t = start; t < end; t += interval) {
        let t0 = t, t1 = Math.min(t + interval, end);
        let d0 = sunDistanceDerivative(t0), d1 = sunDistanceDerivative(t1);
        if (d0 === 0) {((d1 > 0) ? perihelion : aphelion).push(t0);}
        else if (d0 * d1 < 0) {
            while (t1 - t0 > window) {
                const tAvg = (t0 + t1) / 2;
                const dAvg = sunDistanceDerivative(tAvg);
                if (d1 * dAvg <= 0) {t0 = tAvg; d0 = dAvg;}
                else {t1 = tAvg; d1 = dAvg;}
            }
            const tApsis = Math.floor(t0 + (d0 / (d0 - d1)) * (t1 - t0)); // linear interpolation of derivative
            ((d1 > 0) ? perihelion : aphelion).push(tApsis);
        }
    }
    return {aphelion: aphelion, perihelion: perihelion}
}

/**
 * Returns intervals of day, civil twilight, nautical twilight, astronomical twilight, and night during a particular day.
 * Time is measured in milliseconds since midnight. It is not adjusted for DST, so "4 pm" is always represented as 16*60*60*1000 =
 * 57600000.
 * @param sunEvents Array returned by allSunEvents
 * @param timeZone Time zone of the given location (ex. "America/Los_Angeles") or a time zone lookup table
 * @returns Array with intervals of [day, civil twilight, nautical twilight, astronomical twilight, night].
 */
export function intervals(sunEvents: SEvent[], timeZone: TimeChange[]) {
    const newSunEvents = []; // sunEvents without solar noon or midnight
    const ints : number[][][] = [[], [], [], [], []]; // intervals of day, civil twilight, nautical twilight, astronomical twilight, and night

    for (const event of sunEvents) {
        if (event.type !== "Solar Noon" && event.type !== "Solar Midnight") {newSunEvents.push(event);}
    }

    if (newSunEvents.length === 0) { // no sunrise, sunset, dawn, or dusk
        const sunAngle = sunEvents[0].elev;
        if (sunAngle < ASTRO_TWILIGHT) {return [[], [], [], [], [[0, DAY_LENGTH]]];}
        else if (sunAngle < NAUTICAL_TWILIGHT) {return [[], [], [], [[0, DAY_LENGTH]], []];}
        else if (sunAngle < CIVIL_TWILIGHT) {return [[], [], [[0, DAY_LENGTH]], [], []];}
        else if (sunAngle < HORIZON) {return [[], [[0, DAY_LENGTH]], [], [], []];}
        else {return [[[0, DAY_LENGTH]], [], [], [], []];}
    }
    
    let etype = newSunEvents[0].type;
    if (etype === "Sunset") {ints[0].push([0, getTimeOfDay(newSunEvents[0].unix, timeZone)]);}
    else if (etype === "Sunrise" || etype === "Civil Dusk") {ints[1].push([0, getTimeOfDay(newSunEvents[0].unix, timeZone)]);}
    else if (etype === "Civil Dawn" || etype === "Nautical Dusk") {ints[2].push([0, getTimeOfDay(newSunEvents[0].unix, timeZone)]);}
    else if (etype === "Nautical Dawn" || etype === "Astro Dusk") {ints[3].push([0, getTimeOfDay(newSunEvents[0].unix, timeZone)]);}
    else if (etype === "Astro Dawn") {ints[4].push([0, getTimeOfDay(newSunEvents[0].unix, timeZone)]);}

    for (let i=0; i<newSunEvents.length-1; i++) {
        etype = newSunEvents[i+1].type;
        const t0 = getTimeOfDay(newSunEvents[i].unix, timeZone);
        const t1 = getTimeOfDay(newSunEvents[i+1].unix, timeZone);
        if (etype === "Sunset") {ints[0].push([t0, t1]);}
        else if (etype === "Sunrise" || etype === "Civil Dusk") {ints[1].push([t0, t1]);}
        else if (etype === "Civil Dawn" || etype === "Nautical Dusk") {ints[2].push([t0, t1]);}
        else if (etype === "Nautical Dawn" || etype === "Astro Dusk") {ints[3].push([t0, t1]);}
        else if (etype === "Astro Dawn") {ints[4].push([t0, t1]);}
    }

    const lastTime = getTimeOfDay(newSunEvents[newSunEvents.length-1].unix, timeZone);
    if (etype === "Sunrise") {ints[0].push([lastTime, DAY_LENGTH]);}
    else if (etype === "Civil Dawn" || etype === "Sunset") {ints[1].push([lastTime, DAY_LENGTH]);}
    else if (etype === "Nautical Dawn" || etype === "Civil Dusk") {ints[2].push([lastTime, DAY_LENGTH]);}
    else if (etype === "Astro Dawn" || etype === "Nautical Dusk") {ints[3].push([lastTime, DAY_LENGTH]);}
    else if (etype === "Astro Dusk") {ints[4].push([lastTime, DAY_LENGTH]);}
    return ints;
}

/**
 * Intervals of daylight and each stage of twilight for use in SVG diagram generation.
 * @param sunEvents Value returned from the "allSunEvents" function on the given day.
 * @param timeZone Time zone, either as IANA identifier or lookup table.
 * @returns Array of arrays of arrays of numbers: [dIntervals, cIntervals, nIntervals, aIntervals].
 * 
 * dIntervals: Intervals of daylight, where the sun's unrefracted elevation angle >= -5/6°.
 * 
 * cIntervals: Intervals of civil twilight or brighter (sun angle >= -6°).
 * 
 * nIntervals: Intervals of nautical twilight or brighter (sun angle >= -12°).
 * 
 * aIntervals: Intervals of astronomical twilight or brighter (sun angle >= -18°).
 */
export function intervalsSvg(sunEvents: SEvent[], timeZone: TimeChange[]): number[][][] {
    const newSunEvents = []; // sunEvents without solar noon or midnight

    for (const event of sunEvents) {
        if (event.type !== "Solar Noon" && event.type !== "Solar Midnight") {newSunEvents.push(event);}
    }

    if (newSunEvents.length === 0) { // no sunrise, sunset, dawn, or dusk
        const s = sunEvents[0].elev;
        if (s < ASTRO_TWILIGHT) {return [[], [], [], []];}
        else if (s < NAUTICAL_TWILIGHT) {return [[], [], [], [[0, 86400]]];}
        else if (s < CIVIL_TWILIGHT) {return [[], [], [[0, 86400]], [[0, 86400]]];}
        else if (s < HORIZON) {return [[], [[0, 86400]], [[0, 86400]], [[0, 86400]]];}
        else {return [[[0, 86400]], [[0, 86400]], [[0, 86400]], [[0, 86400]]];}
    }

    const dIntervals: number[][] = []; // intervals of daylight
    const cIntervals: number[][] = []; // daylight + civil twilight
    const nIntervals: number[][] = []; // daylight + civil twilight + nautical twilight
    const aIntervals: number[][] = []; // daylight + civil twilight + nautical twilight + astronomical twilight
    
    let etype = newSunEvents[0].type;
    let s = Math.floor(getTimeOfDay(newSunEvents[0].unix, timeZone) / 1000);

    // push the first interval
    if (etype === "Astro Dawn") {
        aIntervals.push([s, 86400]);
    } else if (etype === "Nautical Dawn") {
        aIntervals.push([0, 86400]);
        nIntervals.push([s, 86400]);
    } else if (etype === "Civil Dawn") {
        aIntervals.push([0, 86400]);
        nIntervals.push([0, 86400]);
        cIntervals.push([s, 86400]);
    } else if (etype === "Sunrise") {
        aIntervals.push([0, 86400]);
        nIntervals.push([0, 86400]);
        cIntervals.push([0, 86400]);
        dIntervals.push([s, 86400]);
    } else if (etype === "Sunset") {
        aIntervals.push([0, 86400]);
        nIntervals.push([0, 86400]);
        cIntervals.push([0, 86400]);
        dIntervals.push([0, s]);
    } else if (etype === "Civil Dusk") {
        aIntervals.push([0, 86400]);
        nIntervals.push([0, 86400]);
        cIntervals.push([0, s]);
    } else if (etype === "Nautical Dusk") {
        aIntervals.push([0, 86400]);
        nIntervals.push([0, s]);
    } else if (etype === "Astro Dusk") {
        aIntervals.push([0, s]);
    }

    for (let i=1; i<newSunEvents.length; i++) {
        etype = newSunEvents[i].type;
        s = Math.floor(getTimeOfDay(newSunEvents[i].unix, timeZone)/1000);
        if (etype === "Astro Dawn") {aIntervals.push([s, 86400]);}
        else if (etype === "Nautical Dawn") {nIntervals.push([s, 86400]);}
        else if (etype === "Civil Dawn") {cIntervals.push([s, 86400]);}
        else if (etype === "Sunrise") {dIntervals.push([s, 86400]);}
        else if (etype === "Sunset") {dIntervals[dIntervals.length-1][1] = s;}
        else if (etype === "Civil Dusk") {cIntervals[cIntervals.length-1][1] = s;}
        else if (etype === "Nautical Dusk") {nIntervals[nIntervals.length-1][1] = s;}
        else if (etype === "Astro Dusk") {aIntervals[aIntervals.length-1][1] = s;}
    }
    return [dIntervals, cIntervals, nIntervals, aIntervals];
}

export function intervalsNightCivilTwilight(sunEvents: SEvent[], timeZone: TimeChange[]): number[][][] {
    const newSunEvents = []; // sunEvents without solar noon or midnight
    for (const event of sunEvents) {
        if (event.type === "Sunrise" || event.type === "Sunset" || event.type.includes("Civil")) {newSunEvents.push(event);}
    }
    if (newSunEvents.length === 0) {
        const s = sunEvents[0].elev;
        if (s < CIVIL_TWILIGHT) {return [[[0, 86400]], []];}
        else if (s < HORIZON) {return [[], [[0, 86400]]];}
        else {return [[], []];}
    }

    const nIntervals: number[][] = []; // night intervals
    const cIntervals: number[][] = []; // civil twilight intervals

    let etype = newSunEvents[0].type;
    let s = Math.floor(getTimeOfDay(newSunEvents[0].unix, timeZone)/1000);
    if (etype === "Civil Dawn") {
        nIntervals.push([0, s]);
        cIntervals.push([s, 86400]);
    }
    else if (etype === "Sunrise") {
        cIntervals.push([0, s]);
    }
    else if (etype === "Sunset") {
        cIntervals.push([s, 86400]);
    }
    else if (etype === "Civil Dusk") {
        cIntervals.push([0, s]);
        nIntervals.push([s, 86400]);
    }

    for (let i=1; i<newSunEvents.length; i++) {
        etype = newSunEvents[i].type;
        s = Math.floor(getTimeOfDay(newSunEvents[i].unix, timeZone)/1000);
        if (etype === "Civil Dawn") {
            nIntervals[nIntervals.length-1][1] = s;
            cIntervals.push([s, 86400]);
        }
        else if (etype === "Sunrise") {
            cIntervals[cIntervals.length-1][1] = s;
        }
        else if (etype === "Sunset") {
            cIntervals.push([s, 86400]);
        }
        else if (etype === "Civil Dusk") {
            cIntervals[cIntervals.length-1][1] = s;
            nIntervals.push([s, 86400]);
        }
    }
    return [nIntervals, cIntervals];
}

/**
 * Returns the lengths of day combined with different stages of twilight.
 * @param sunEvents The return value of the allSunEvents command at a particular place and date.
 * @param timeZone Time zone, either as IANA identifier or lookup table.
 * @returns An array [t0, t1, t2, t3]. The values are as follows:
 * @t0: Day length
 * @t1: Day + civil twilight
 * @t2: Day + civil twilight + nautical twilight
 * @t3: Day + civil twilight + nautical twilight + astronomical twilight
 */
export function lengths(sunEvents: SEvent[], timeZone: TimeChange[]): number[] {
    const newSunEvents = []; // sunEvents without solar noon or midnight
    const durations = [0, 0, 0, 0]; // durations
    for (const event of sunEvents) {
        if (event.type !== "Solar Noon" && event.type !== "Solar Midnight") {newSunEvents.push(event);}
    }

    if (newSunEvents.length === 0) { // no sunrise, sunset, dawn, or dusk
        const s = sunEvents[0].elev;
        if (s < ASTRO_TWILIGHT) {return [0, 0, 0, 0];} // night all day
        else if (s < NAUTICAL_TWILIGHT) {return [0, 0, 0, 86400];} // astronomical twilight all day
        else if (s < CIVIL_TWILIGHT) {return [0, 0, 86400, 86400];} // nautical twilight all day
        else if (s < HORIZON) {return [0, 86400, 86400, 86400];} // civil twilight all day
        else {return [86400, 86400, 86400, 86400];} // daylight all day
    }

    let etype = newSunEvents[0].type;
    let ms = getTimeOfDay(newSunEvents[0].unix, timeZone);
    if (etype === "Nautical Dawn" || etype === "Astro Dusk") {durations[3] += ms;}
    else if (etype === "Civil Dawn" || etype === "Nautical Dusk") {durations[2] += ms;}
    else if (etype === "Sunrise" || etype === "Civil Dusk") {durations[1] += ms;}
    else if (etype === "Sunset") {durations[0] += ms;}

    for (let i=0; i<newSunEvents.length-1; i++) {
        etype = newSunEvents[i+1].type;
        ms = getTimeOfDay(newSunEvents[i+1].unix, timeZone) - getTimeOfDay(newSunEvents[i].unix, timeZone);
        if (etype === "Nautical Dawn" || etype === "Astro Dusk") {durations[3] += ms;}
        else if (etype === "Civil Dawn" || etype === "Nautical Dusk") {durations[2] += ms;}
        else if (etype === "Sunrise" || etype === "Civil Dusk") {durations[1] += ms;}
        else if (etype === "Sunset") {durations[0] += ms;}
    }

    ms = DAY_LENGTH - getTimeOfDay(newSunEvents[newSunEvents.length-1].unix, timeZone);
    if (etype === "Astro Dawn" || etype === "Nautical Dusk") {durations[3] += ms;}
    else if (etype === "Nautical Dawn" || etype === "Civil Dusk") {durations[2] += ms;}
    else if (etype === "Civil Dawn" || etype === "Sunset") {durations[1] += ms;}
    else if (etype === "Sunrise") {durations[0] += ms;}

    return [Math.floor(durations[0]/1000), 
    Math.floor((durations[0]+durations[1])/1000), 
    Math.floor((durations[0]+durations[1]+durations[2])/1000),
    Math.floor((durations[0]+durations[1]+durations[2]+durations[3])/1000)];
}

/** Generate an object containing sunrise, sunset, twilight, solar noon and solar midnight times, day lengths,
 * solstices and equinoxes (grouped under "solstices"), perihelion, aphelion, and the time zone lookup table for
 * the given zone in the given year.
 * Note that latitude and longitude are given in degrees here, not radians, because they are treated as geographic coordinates
 * rather than angles.
 */
export function generateSunTable(lat: number, long: number, year: number, zone: string): SunTable {
    lat *= degToRad; long *= degToRad;
    const yearStart = DateTime.fromObject({year: year}, {zone: zone});
    const yearEnd = DateTime.fromObject({year: year + 1}, {zone: zone});

    // dayStarts has 1-day buffer due to nature of the dayLengths function
    const dayStarts = mf.dayStarts(yearStart.minus({days:1}), yearEnd.plus({days:1}));
    const timeZoneTable = timeZoneLookupTable(dayStarts.slice(1, -1)); // time zone lookup table
    const lodLookupTable = longDistLookupTable(dayStarts);

    const sunEventsYear: SEvent[][] = [];
    for (let i=0; i<lodLookupTable.length-1; i++) {
        const start = lodLookupTable[i], end = lodLookupTable[i+1];
        sunEventsYear.push(allSunEvents(lat, long, start, end));
    }
    const sunEventsYearM = sunEventsYear.slice(1, -1);

    const dayLengths: number[] = [];
    for (let i=1; i<sunEventsYear.length-1; i++) {
        const l = dayLength(mf.ms(dayStarts[i]), sunEventsYear[i-1], sunEventsYear[i], sunEventsYear[i+1]);
        dayLengths.push(l);
    }

    const solstices = [calcSolstEq(year, 3), calcSolstEq(year, 6), calcSolstEq(year, 9), calcSolstEq(year, 12)];
    const solsticeDates = [
        DateTime.fromMillis(solstices[0], {zone: zone}).diff(yearStart, ["days", "hours"]).toObject().days!,
        DateTime.fromMillis(solstices[1], {zone: zone}).diff(yearStart, ["days", "hours"]).toObject().days!,
        DateTime.fromMillis(solstices[2], {zone: zone}).diff(yearStart, ["days", "hours"]).toObject().days!,
        DateTime.fromMillis(solstices[3], {zone: zone}).diff(yearStart, ["days", "hours"]).toObject().days!
    ]; // days of the solstice/equinox, where January 1 is "day 0"

    const apsides = sunApsides(mf.ms(yearStart), mf.ms(yearEnd));

    return {
        solarEvents: sunEventsYearM,
        dayLengths: dayLengths,
        solstices: solstices,
        solsticeDates: solsticeDates,
        perihelion: apsides.perihelion,
        aphelion: apsides.aphelion,
        timeZoneTable: timeZoneTable
    }
}