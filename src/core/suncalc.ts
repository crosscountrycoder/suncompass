/* 
The formulas for eccentricity, earth-sun distance, axial tilt, declination, equation of time, and atmospheric refraction of sunlight 
are borrowed from the book "Astronomical Algorithms" by Jean Meeus. The refraction formula is modified slightly to ensure continuity 
when the sun is below the horizon. The formula for solar ecliptic longitude is from the book "Planetary Programs and Tables from 
-4000 to +2800" by Pierre Bretagnon and Jean-Louis Simon.

The subsolar point is calculated using the declination and equation of time, treating UTC as equivalent to UT1 (mean solar time at 0
degrees longitude), the two are always within 0.9 seconds of each other. Solar noon and midnight are calculated using the equation of 
time. The position of the sun is calculated from the subsolar point through spherical trigonometry.

For the purposes of this calculator, daytime is defined as the period in which the solar elevation angle is greater than -50 
arcminutes, or -5/6 of a degree. This accounts for the sun's angular radius in the sky and refraction of sunlight. It corresponds to 
the definition used by most sources, including NOAA's solar calculator.

This site uses the Luxon library to deal with date/time computations. Luxon is used to simplify computation when dealing with
durations, conversion between different time zones, and complexities such as daylight saving time. The geo-tz library is used to find 
the time zone of a geographic coordinate.
*/

import * as mf from "./mathfuncs.ts";
import {degToRad, sunPeriodicTerms, DAY_LENGTH, BSEARCH_GAP} from "./constants.ts";
import {generateLODProfile, estimateLOD, getTimeOfDay} from "./lookup-tables.ts";
import type {TimeChange, LODProfile, SEvent} from "./lookup-tables.ts";
import {DateTime} from "luxon";

export type SeasonEvents = {marEquinox: DateTime; junSolstice: DateTime; sepEquinox: DateTime; decSolstice: DateTime;};

/** Sun's mean longitude according to formula 27.2 in Astronomical Algorithms. 
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
*/
export function sunMeanLong(date: number, unix = false): number {
    if (!unix) { // if date is specified as a Julian century
        const m = date / 10; // Julian millennium
        return mf.mod(280.4664567 + 360007.6982779*m + 0.03032028*m**2 + m**3/49931 - m**4/15299 - m**5/1988000, 360);
    }
    else {return sunMeanLong(mf.jCentury(date));}
}

/** Sun's geometric longitude, i.e. longitude excluding aberration and nutation. 
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
*/
export function sunGeomLong(date: number, unix = false): number {
    if (!unix) {
        const U = date / 100;
        let long = 4.9353929 + 62833.196168*U;
        for (let i=0; i<sunPeriodicTerms.length; i++) {
            const curRow = sunPeriodicTerms[i];
            long += (1e-7 * curRow[0] * Math.sin(curRow[2]+curRow[3]*U));
        }
        return mf.mod(long / degToRad, 360);
    }
    else {return sunGeomLong(mf.jCentury(date));}
}

/** Formula 45.3, in page 308 of Astronomical Algorithms */
export function meanSunAnomaly(JC: number): number {return mf.mod(357.5291092 + 35999.0502909*JC - 1.536e-4*JC**2 + JC**3/24490000, 360);}
export function eccentricity(JC: number): number {return 0.016708617 - 4.2037e-5*JC - 1.236e-7*JC**2;}

/* Unused function
export function equationOfCenter(JC: number): number {
    const anom = meanSunAnomaly(JC) * degToRad;
    return (1.9146-0.004817*JC-1.4e-5*JC**2)*Math.sin(anom) + (0.019993-1.01e-4*JC)*Math.sin(2*anom) + 2.9e-4*Math.sin(3*anom);
} */

/** Distance from sun to earth in kilometers. 
 * @param date The timestamp. Can be specified as a Luxon DateTime, a Unix timestamp, or a Julian century.
 * @param unix If date is a number, this parameter specifies whether it's measured in Unix milliseconds or Julian centuries.
*/
export function sunDistance(date: number | DateTime, unix = false): number {
    if (typeof(date) === "number") {
        if (!unix) {
            const U = date / 100;
            let dist = 1.0001026; // in astronomical units
            for (let i=0; i<sunPeriodicTerms.length; i++) {
                const curRow = sunPeriodicTerms[i];
                dist += (1e-7 * curRow[1] * Math.cos(curRow[2]+curRow[3]*U));
            }
            return 149597870.7 * dist; // convert to kilometers
        }
        else {return sunDistance(mf.jCentury(date));}
    }
    else {return sunDistance(mf.jCentury(mf.ms(date)));}
}

/**
 * Calculates the sun's apparent ecliptic longitude to within 0.0009 degrees for years 0-3000. This value is 0 at the March equinox,
 * 90 at the June solstice, 180 at the September equinox, and 270 at the December solstice.
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
 */
export function sunTrueLong(date: number, unix = false): number {
    if (!unix) {
        const U = date / 100;
        const geoLong = sunGeomLong(date);
        const aberration = 1e-7*(17*Math.cos(3.1+62830.14*U)-993)/degToRad;
        const nutation = 1e-7*(-834*Math.sin(2.18-3375.7*U+0.36*U**2)-64*Math.sin(3.51+125666.39*U+0.1*U**2))/degToRad;
        return mf.mod(geoLong + aberration + nutation, 360);
    }
    else {return sunTrueLong(mf.jCentury(date));}
}

/** Nutation in obliquity in degrees 
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
*/
export function obNutation(date: number, unix = false) {
    if (!unix) {
        const L = mf.mod(280.4665 + 36000.7698*date, 360)*degToRad;
        const Lprime = mf.mod(218.3165 + 481267.8813*date, 360)*degToRad;
        const omega = mf.mod(125.04452 - 1934.136261*date + 0.0020708*date**2 + date**3/450000, 360)*degToRad;
        return (9.2*Math.cos(omega) + 0.57*Math.cos(2*L) + 0.1*Math.cos(2*Lprime) - 0.09*Math.cos(2*omega)) / 3600;
    }
    else {return obNutation(mf.jCentury(date));}
}

/** Nutation in longitude in degrees 
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
*/
export function longNutation(date: number, unix = false) {
    if (!unix) {
        const L = mf.mod(280.4665 + 36000.7698*date, 360)*degToRad;
        const Lprime = mf.mod(218.3165 + 481267.8813*date, 360)*degToRad;
        const omega = mf.mod(125.04452 - 1934.136261*date + 0.0020708*date**2 + date**3/450000, 360)*degToRad;
        return (-17.2*Math.sin(omega) - 1.32*Math.sin(2*L) - 0.23*Math.sin(2*Lprime) + 0.21*Math.sin(2*omega)) / 3600;
    }
    else {return longNutation(mf.jCentury(date));}
}

/**
 * Returns the obliquity of the ecliptic, or equivalently Earth's axial tilt.
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
 */
export function obliquity(date: number, unix = false): number {
    if (!unix) {
        const meanObliquity = 23.4392911 + (-46.815*date - 5.9e-4*date**2 + 1.813e-3*date**3) / 3600;
        return meanObliquity + obNutation(date);
    }
    else {return obliquity(mf.jCentury(date));}
}

/**
 * Calculates the sun's right ascension in degrees. To convert to hours, divide by 15.
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
 */
export function sunRA(date: number, unix = false): number {
    if (!unix) {
        const long = sunTrueLong(date) * degToRad;
        const ob = obliquity(date) * degToRad;
        const ra = Math.atan2(Math.sin(long)*Math.cos(ob), Math.cos(long));
        return mf.mod(ra / degToRad, 360);
    }
    else {return sunRA(mf.jCentury(date));}
}

/**
 * Equation of time in minutes (apparent solar time - mean solar time). Based on equation 27.1 in Astronomical Algorithms.
 * @param date The timestamp.
 * @param unix If true, date is measured in Unix milliseconds. If false, date is Julian centuries since J2000 epoch.
 */
export function equationOfTime(date: number, unix = false): number { 
    if (!unix) {
        const eot = sunMeanLong(date) - 0.0057183 - sunRA(date) + longNutation(date) * Math.cos(obliquity(date)*degToRad); // in degrees
        return mf.mod(4 * eot + 720, 1440) - 720; // reduce to range [-720, 720] if absolute value too large
    }
    else {return equationOfTime(mf.jCentury(date));}
}

/** Gives the value of Greenwich apparent sidereal time (GAST) in degrees, from equation 11.4 in Astronomical Algorithms.
 * The value returned is in the range 0 <= x < 360.
 * @param unix Unix timestamp in milliseconds.
 */
export function gast(unix: number): number {
    const JD = mf.jdUTC(unix); // Julian day but using UTC (approximation to UT1) rather than TT
    const JC = (JD - 2451545) / 36525; // Julian century, but with UTC rather than TT
    const gmst = 280.46061837 + 360.98564736629*(JD-2451545) + 3.87933e-4*JC**2 - JC**3/38710000;
    const correction = longNutation(unix, true) * Math.cos(obliquity(unix, true) * degToRad);
    return mf.mod(gmst + correction, 360);
}

/**
 * Returns apparent solar time given longitude and time, in minutes after solar midnight. 0 is solar midnight, 720 is solar noon.
 * @param longitude Longitude in degrees.
 * @param unix Unix timestamp in milliseconds.
 */
export function solarTime(longitude: number, unix: number): number {
    const timeEq = equationOfTime(unix, true);
    const mst = unix / 60000 + longitude * 4; // mean solar time in minutes
    const ast = mf.mod(timeEq + mst, 1440); // apparent solar time
    return ast;
}

/**
 * Returns the sun's declination in degrees. This is the latitude of the subsolar point.
 * @param longitude Sun's true ecliptic longitude, in degrees
 * @param obliquity Obliquity of the ecliptic (i.e. Earth's axial tilt)
 */
export function declination(longitude: number, obliquity: number): number {
    return Math.asin(Math.sin(obliquity*degToRad)*Math.sin(longitude*degToRad)) / degToRad;
}

/**
 * Returns the subsolar point, or location on Earth at which the sun is directly overhead.
 * @param date Longitude, obliquity, distance (LOD) profile for the given moment.
 * @returns [latitude, longitude] of subsolar point
 */
export function subsolarPoint(lod?: LODProfile): number[] {
    if (lod === undefined) {lod = generateLODProfile(Date.now());}
    const subsolarLat = declination(lod.longitude, lod.obliquity);
    const soltime0 = solarTime(0, lod.unix); // solar time at Greenwich meridian (longitude 0)
    const subsolarLong = mf.mod(-soltime0/4, 360) - 180;
    return [subsolarLat, subsolarLong];
}

/**
 * Returns sun position given latitude, longitude, and DateTime.
 * @param lat Latitude in degrees
 * @param long Longitude in degrees
 * @param lod Longitude, obliquity, distance (LOD) profile.
 * @returns Array: [elevation, azimuth]. Elevation is in degrees above horizon, azimuth is degrees clockwise from north
 * Solar elevation is not refracted. To find the solar elevation angle adjusted for atmospheric refraction, use mf.refract(sunPosition[0])
 */
export function sunPosition(lat: number, long: number, lod: LODProfile): number[] {
    const [sunLat, sunLong] = subsolarPoint(lod);
    const sunEcef = mf.toEcef(sunLat, sunLong, lod.distance);
    /* Note: Geodetic latitude of subsolar point is the same as geocentric latitude at which the sun-Earth center line intersects 
    the ellipsoid. Subsolar point is where the surface normal intersects the sun. */
    return mf.elevAzimuth(lat, long, mf.latLongEcef(lat, long), sunEcef);
}

/**
 * Returns the time(s) of solar noon, along with the sun's position at solar noon.
 * @param latitude Latitude in degrees.
 * @param longitude Longitude in degrees.
 * @param start Start (as LOD profile)
 * @param end End (as LOD profile)
 * @returns Time(s) of solar noon, along with the sun's position at solar noon.
 */
export function solarNoon(lat: number, long: number, start: LODProfile, end: LODProfile): SEvent[] {
    const startT = start.unix, endT = end.unix;
    const st00 = solarTime(long, startT), st24 = solarTime(long, endT);
    const stDiff = mf.mod(st24 - st00 - 720, 1440) + 720;
    const solarTimeRate = stDiff / (endT - startT); // rate of solar time change relative to actual time
    if (st00 > 600 && st00 <= 720 && st24 > 720 && st24 < 840) { // 2 solar noons in a day
        let solarNoon0 = startT + (720 - st00) / solarTimeRate;
        solarNoon0 += ((720 - solarTime(long, solarNoon0)) * 60000); // refine to 1 ms precision
        solarNoon0 = Math.floor(mf.clamp(solarNoon0, startT, endT-1));
        const [e0, a0] = sunPosition(lat, long, estimateLOD(solarNoon0,start,end)); // solar elevation/azimuth at solarNoon0

        let solarNoon1 = endT - 60000 * (st24 - 720) / solarTimeRate;
        solarNoon1 += ((720 - solarTime(long, solarNoon1)) * 60000);
        solarNoon1 = Math.floor(mf.clamp(solarNoon1, startT, endT-1));
        const [e1, a1] = sunPosition(lat, long, estimateLOD(solarNoon1,start,end)); // solar elevation/azimuth at solarNoon1
        return [{unix:solarNoon0,type:"Solar Noon",elev:e0,azimuth:a0},{unix:solarNoon1,type:"Solar Noon",elev:e1,azimuth:a1}];
    }
    else if (st00 > 720 && st00 < 840 && st24 > 600 && st24 <= 720) { // 0 solar noons in a day
        return [];
    }
    else { // 1 solar noon in a day
        let solarNoon = startT + mf.mod(720 - st00, 1440) / solarTimeRate;
        solarNoon += ((720 - solarTime(long, solarNoon)) * 60000);
        solarNoon = Math.floor(mf.clamp(solarNoon, startT, endT-1));
        const [e, a] = sunPosition(lat, long, estimateLOD(solarNoon,start,end));
        return [{unix: solarNoon, type: "Solar Noon", elev: e, azimuth: a}];
    }
}

/**
 * Returns the time(s) of solar midnight, along with the sun's position at solar midnight.
 * @param latitude Latitude in degrees.
 * @param longitude Longitude in degrees.
 * @param start Start (as LOD profile)
 * @param end End (as LOD profile)
 * @returns Time(s) of solar midnight, along with the sun's position at solar midnight.
 */
export function solarMidnight(lat: number, long: number, start: LODProfile, end: LODProfile): SEvent[] {
    const startT = start.unix, endT = end.unix;
    const st00 = solarTime(long, startT), st24 = solarTime(long, endT);
    const stDiff = mf.mod(st24 - st00 - 720, 1440) + 720;
    const solarTimeRate = stDiff / (endT - startT); // rate of solar time change relative to actual time

    if (st00 > 1320 && st24 < 120) { // 2 solar midnights in a day
        let solarMidnight0 = startT + (1440 - st00) / solarTimeRate;
        solarMidnight0 += ((720-mf.mod(solarTime(long,solarMidnight0)+720,1440))*60000);
        solarMidnight0 = Math.floor(mf.clamp(solarMidnight0, startT, endT-1));
        const [e0, a0] = sunPosition(lat, long, estimateLOD(solarMidnight0,start,end)); // solar elevation/azimuth at solarMidnight0

        let solarMidnight1 = endT - st24 / solarTimeRate;
        solarMidnight1 += ((720-mf.mod(solarTime(long,solarMidnight1)+720,1440))*60000);
        solarMidnight1 = Math.floor(mf.clamp(solarMidnight1, startT, endT-1));
        const [e1, a1] = sunPosition(lat, long, estimateLOD(solarMidnight1,start,end)); // solar elevation/azimuth at solarMidnight1
        return [{unix: solarMidnight0, type: "Solar Midnight", elev: e0, azimuth: a0},
            {unix: solarMidnight1, type: "Solar Midnight", elev: e1, azimuth: a1}];
    }
    else if (st00 < 120 && st24 > 1320) { // 0 solar midnights in a day
        return [];
    }
    else { // 1 solar midnight in a day
        let solarMidnight = endT - st24 / solarTimeRate;
        solarMidnight += ((720-mf.mod(solarTime(long,solarMidnight)+720,1440))*60000);
        solarMidnight = Math.floor(mf.clamp(solarMidnight, startT, endT-1));
        const [e, a] = sunPosition(lat, long, estimateLOD(solarMidnight,start,end));
        return [{unix: solarMidnight, type: "Solar Midnight", elev: e, azimuth: a}];
    }
}

/** Returns the approximate derivative of the solar elevation angle at a particular time, in degrees per second. */
export function derivative(lat: number, long: number, unix: number, startLOD: LODProfile, endLOD: LODProfile) {
    const t0 = unix - 500, t1 = unix + 500;
    const LOD0 = estimateLOD(t0, startLOD, endLOD), LOD1 = estimateLOD(t1, startLOD, endLOD);
    return sunPosition(lat, long, LOD1)[0] - sunPosition(lat, long, LOD0)[0];
}

/**
 * Returns an array of Unix times, representing the start of the current day, the times at which the derivative of solar
 * elevation angle is 0, and the start of the next day. This is a helper function for "dawn" and "dusk" functions below.
 * @param lat Latitude in degrees
 * @param long Longitude in degrees
 * @param start Start of day (LOD profile)
 * @param end End of day (LOD profile)
 * @returns Unix timestamps at which sun reaches relative min or max altitude, along with both the start and end of the day.
 */
export function maxAndMin(lat: number, long: number, start: LODProfile, end: LODProfile): number[] {
    const startTime = start.unix, endTime = end.unix;
    const times = [startTime];
    const intervals = [startTime,startTime+4*3.6e6,startTime+8*3.6e6,startTime+12*3.6e6,startTime+16*3.6e6,startTime+20*3.6e6,endTime];
    for (let i=0; i<intervals.length-1; i++) {
        // use binary search to find the time closest to zero derivative
        let t0 = intervals[i], t1 = intervals[i+1];
        let d0 = derivative(lat, long, t0, start, end), d1 = derivative(lat, long, t1, start, end);
        if (d0 >= 0 && d1 < 0) { // maximum (i.e. solar noon, or summer solstice at pole)
            while (t1 - t0 > BSEARCH_GAP) {
                const tAvg = Math.floor((t0+t1)/2);
                const dAvg = derivative(lat, long, tAvg, start, end);
                if (dAvg >= 0) {t0 = tAvg; d0 = dAvg;}
                else {t1 = tAvg; d1 = dAvg;}
            }
            // use 5-minute window, then linear interpolation to make calculation faster
            const t = t0 + (d0 / (d0 - d1)) * (t1 - t0);
            times.push(t);
        }
        else if (d0 <= 0 && d1 > 0) { // minimum (i.e. solar midnight, or winter solstice at pole)
            while (t1 - t0 > BSEARCH_GAP) {
                const tAvg = Math.floor((t0+t1)/2);
                const dAvg = derivative(lat, long, tAvg, start, end);
                if (dAvg <= 0) {t0 = tAvg; d0 = dAvg;}
                else {t1 = tAvg; d1 = dAvg;}
            }
            const t = t0 + (d0 / (d0 - d1)) * (t1 - t0);
            times.push(t);
        }
    }
    times.push(endTime);
    return times;
}

/**
 * Calculates the time in the morning at which the sun's elevation reaches the specified angle. Angle should be -5/6 for sunrise,
 * -6 for civil twilight, -12 for nautical twilight, and -18 for astronomical twilight.
 * @param lat Latitude in degrees
 * @param long Longitude in degrees
 * @param angle Solar elevation angle in degrees
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
        let e0 = sunPosition(lat, long, lod0)[0];
        let e1 = sunPosition(lat, long, lod1)[0];
        if (e0 <= angle && e1 >= angle) {
            while (t1 - t0 > BSEARCH_GAP) {
                const avgLOD = estimateLOD(Math.floor((t0+t1)/2), startLOD, endLOD);
                const eAvg = sunPosition(lat, long, avgLOD)[0];
                if (eAvg <= angle) {t0 = avgLOD.unix; e0 = eAvg;}
                else {t1 = avgLOD.unix; e1 = eAvg;}
            }
            // after reducing to 1-minute interval, use quadratic interpolation
            const d0 = derivative(lat, long, t0, startLOD, endLOD) / 1000;
            const t = Math.floor(mf.quadraticZero(t0, t1, e0-angle, e1-angle, d0));
            const [e, a] = sunPosition(lat, long, estimateLOD(t, startLOD, endLOD));
            dawnTimes.push({unix: t, type: type, elev: e, azimuth: a});
        }
    }
    return dawnTimes;
}

/**
 * Calculates the time in the evening at which the sun's elevation reaches the specified angle. Angle should be -5/6 for sunset,
 * -6 for civil twilight, -12 for nautical twilight, and -18 for astronomical twilight.
 * @param lat Latitude in degrees
 * @param long Longitude in degrees
 * @param angle Solar elevation angle in degrees
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
        let e0 = sunPosition(lat, long, lod0)[0];
        let e1 = sunPosition(lat, long, lod1)[0];
        if (e0 >= angle && e1 <= angle) {
            while (t1 - t0 > BSEARCH_GAP) {
                const avgLOD = estimateLOD(Math.floor((t0+t1)/2), startLOD, endLOD);
                const eAvg = sunPosition(lat, long, avgLOD)[0];
                if (eAvg >= angle) {t0 = avgLOD.unix; e0 = eAvg;}
                else {t1 = avgLOD.unix; e1 = eAvg;}
            }
            // after reducing to 1-minute interval, use quadratic interpolation
            const d0 = derivative(lat, long, t0, startLOD, endLOD) / 1000;
            const t = Math.floor(mf.quadraticZero(t0, t1, e0-angle, e1-angle, d0));
            const [e, a] = sunPosition(lat, long, estimateLOD(t, startLOD, endLOD));
            duskTimes.push({unix: t, type: type, elev: e, azimuth: a});
        }
    }
    return duskTimes;
}

export function sunrise(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dawn(lat, long, -5/6, "Sunrise", maxMin, startLOD, endLOD);
} 
export function sunset(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dusk(lat, long, -5/6, "Sunset", maxMin, startLOD, endLOD);
}
export function civilDawn(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dawn(lat, long, -6, "Civil Dawn", maxMin, startLOD, endLOD);
}
export function civilDusk(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dusk(lat, long, -6, "Civil Dusk", maxMin, startLOD, endLOD);
}
export function nauticalDawn(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dawn(lat, long, -12, "Nautical Dawn", maxMin, startLOD, endLOD);
}
export function nauticalDusk(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dusk(lat, long, -12, "Nautical Dusk", maxMin, startLOD, endLOD);
}
export function astroDawn(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dawn(lat, long, -18, "Astro Dawn", maxMin, startLOD, endLOD);
}
export function astroDusk(lat: number, long: number, maxMin: number[], startLOD: LODProfile, endLOD: LODProfile) {
    return dusk(lat, long, -18, "Astro Dusk", maxMin, startLOD, endLOD);
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
        return (sunEventsToday.length === 0 || sunEventsToday[0].elev >= -5/6) ? DAY_LENGTH : 0;
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
 * @param lat Latitude of observer
 * @param long Longitude of observer
 * @param date Luxon DateTime (can be any point within the given day)
*/
export function sunEventsDay(lat: number, long: number, date: DateTime): SEvent[] {
    const start = generateLODProfile(mf.ms(date.startOf("day")));
    const end = generateLODProfile(mf.ms(date.plus({days:1}).startOf("day")));
    return allSunEvents(lat, long, start, end);
}

/** Calculates the Unix millisecond timestamp of the solstice or equinox in the given month and year.
 * Month must be 3, 6, 9 or 12.
 */
export function calcSolstEq(year = DateTime.now().toUTC().year, month: number) {
    let t0 = mf.ms(DateTime.fromObject({year:year, month:month, day:10}, {zone: "utc"}));
    let t1 = t0 + 18*DAY_LENGTH; // 18 days after start
    while (t1 - t0 > 1) {
        const avg = Math.floor((t0+t1)/2);
        if (month === 3) {(sunTrueLong(avg, true) >= 180) ? t0 = avg : t1 = avg;}
        else {(sunTrueLong(avg, true) <= 30*(month-3)) ? t0 = avg : t1 = avg;}
    }
    return t0;
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
        if (sunAngle < -18) {return [[], [], [], [], [[0, DAY_LENGTH]]];}
        else if (sunAngle < -12) {return [[], [], [], [[0, DAY_LENGTH]], []];}
        else if (sunAngle < -6) {return [[], [], [[0, DAY_LENGTH]], [], []];}
        else if (sunAngle < -5/6) {return [[], [[0, DAY_LENGTH]], [], [], []];}
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
        if (s < -18) {return [[], [], [], []];}
        else if (s < -12) {return [[], [], [], [[0, 86400]]];}
        else if (s < -6) {return [[], [], [[0, 86400]], [[0, 86400]]];}
        else if (s < -5/6) {return [[], [[0, 86400]], [[0, 86400]], [[0, 86400]]];}
        else {return [[[0, 86400]], [[0, 86400]], [[0, 86400]], [[0, 86400]]];}
    }

    const dIntervals: number[][] = []; // intervals of daylight
    const cIntervals: number[][] = []; // daylight + civil twilight
    const nIntervals: number[][] = []; // daylight + civil twilight + nautical twilight
    const aIntervals: number[][] = []; // daylight + civil twilight + nautical twilight + astronomical twilight
    
    let etype = newSunEvents[0].type;
    let s = mf.intDiv(getTimeOfDay(newSunEvents[0].unix, timeZone),1000);

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
        s = mf.intDiv(getTimeOfDay(newSunEvents[i].unix, timeZone),1000);
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
        if (s < -6) {return [[[0, 86400]], []];}
        else if (s < -5/6) {return [[], [[0, 86400]]];}
        else {return [[], []];}
    }

    const nIntervals: number[][] = []; // night intervals
    const cIntervals: number[][] = []; // civil twilight intervals

    let etype = newSunEvents[0].type;
    let s = mf.intDiv(getTimeOfDay(newSunEvents[0].unix, timeZone),1000);
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
        s = mf.intDiv(getTimeOfDay(newSunEvents[i].unix, timeZone),1000);
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
        if (s < -18) {return [0, 0, 0, 0];} // night all day
        else if (s < -12) {return [0, 0, 0, 86400];} // astronomical twilight all day
        else if (s < -6) {return [0, 0, 86400, 86400];} // nautical twilight all day
        else if (s < -5/6) {return [0, 86400, 86400, 86400];} // civil twilight all day
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

    return [mf.intDiv(durations[0],1000), 
    mf.intDiv(durations[0]+durations[1],1000), 
    mf.intDiv(durations[0]+durations[1]+durations[2],1000),
    mf.intDiv(durations[0]+durations[1]+durations[2]+durations[3],1000)];
}