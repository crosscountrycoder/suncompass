import {DateTime} from "luxon";
import * as mf from "./mathfuncs.ts";
import {sunTrueLong, sunDistance, obliquity} from "./suncalc.ts";
import {illumination} from "./mooncalc.ts";
import { DAY_LENGTH, degToRad, HORIZON, TAU } from "./constants.ts";

/** Object representing the change in a location's time zone. 
 * @param unix the Unix timestamp at which the change occurs.
 * @param offset the time zone's new offset, in minutes from UTC
 * @param change True if the time zone changes, false if it's just the initial time zone for the year.
*/
export type TimeChange = {unix: number; offset: number; change: boolean};

/** Object representing solar ecliptic longitude, obliquity, and distance (LOD) at a particular Unix moment.
 * @param unix Unix timestamp in milliseconds.
 * @param longitude Solar ecliptic longitude in radians, range [0, TAU).
 * @param obliquity Obliquity of ecliptic in radians.
 * @param distance Sun-earth distance in kilometers.
 */
export type LODProfile = {unix: number; longitude: number; obliquity: number; distance: number};

/** Object represents an event involving the sun or moon. 
 * @param unix Unix timestamp in milliseconds.
 * @param type Event type. Can be "Solar Midnight", "Astro Dawn", "Nautical Dawn", "Civil Dawn", "Sunrise", "Solar Noon", "Sunset",
 * "Civil Dusk", "Nautical Dusk", "Astro Dusk", "Moonrise", "Moonset", or "Moon Meridian" (i.e. when the moon passes the observer's)
 * meridian.
 * @param elev Sun's elevation above the horizon in radians. Not refracted.
 * @param azimuth Sun's azimuth or compass bearing, in radians clockwise from north.
*/
export type SEvent = {unix: number, type: string, elev: number, azimuth: number};

/** Given the value returned by dayStarts, create a "lookup table" showing when the time offsets change during the given period. */
export function timeZoneLookupTable(dayStarts: DateTime[]): TimeChange[] {
    const changeAtStart = (dayStarts[0].offset !== dayStarts[0].minus(1).offset);
    const firstChange: TimeChange = {unix: mf.ms(dayStarts[0]), offset: dayStarts[0].offset, change: changeAtStart};
    const table: TimeChange[] = [firstChange];

    for (let i=1; i<dayStarts.length; i++) {
        const prevDay = dayStarts[i-1], curDay = dayStarts[i];
        if (prevDay.offset !== curDay.offset) {
            // If the time changes during this day, use binary search to find where it changes.
            let t0 = mf.ms(prevDay), t1 = mf.ms(curDay);
            while (t1 - t0 > 1) {
                const avg = Math.floor((t0 + t1)/2);
                const avgTime = DateTime.fromMillis(avg, {zone: curDay.zone});
                if (avgTime.offset === prevDay.offset) {t0 = avg;}
                else {t1 = avg;}
            }
            table.push({unix: t1, offset: curDay.offset, change: true});
        }
    }
    return table;
}

/** Get UTC offset (minutes) from a Unix timestamp and a time zone lookup table (see mathfuncs.timeZoneLookupTable()) */
export function getOffsetFromTable(unix: number, table: TimeChange[]): number {
    let offset = 0;
    for (const change of table) {
        if (unix >= change.unix) {offset = change.offset;}
        else {break;}
    }
    return offset;
}

export function generateLODProfile(unix: number): LODProfile {
    return {unix: unix, longitude: sunTrueLong(unix, true), obliquity: obliquity(unix, true), distance: sunDistance(unix, true)};
}

/** Create a lookup table for estimating longitude, obliquity, and distance throughout the year using linear interpolation. */
export function longDistLookupTable(dayStarts: DateTime[]): LODProfile[] {
    const table = [];
    for (const date of dayStarts) {table.push(generateLODProfile(mf.ms(date)));}
    return table;
}

/** Estimates solar ecliptic longitude, obliquity, and distance at the given Unix time using linear interpolation 
 * from the start of the day to the end of the day.
 */
export function estimateLOD(unix: number, start: LODProfile, end: LODProfile): LODProfile {
    const diffLong = mf.mod(end.longitude - start.longitude, TAU);
    const diffObliquity = end.obliquity - start.obliquity;
    const diffDistance = end.distance - start.distance;

    const fraction = (unix - start.unix) / (end.unix - start.unix);
    const estLongitude = mf.mod(start.longitude + fraction * diffLong, TAU);
    const estObliquity = start.obliquity + fraction * diffObliquity;
    const estDistance = start.distance + fraction * diffDistance;
    return {unix: unix, longitude: estLongitude, obliquity: estObliquity, distance: estDistance};
}

/** Given a Unix timestamp and the time zone (as a lookup table), return the time of day in milliseconds, as a number
 * between 0 and 86399999 inclusive. If the time is expressed in 24-hour format as h:m:s.ms (where ms is milliseconds),
 * then the return value is equal to ms + 1000*s + 60000*m + 3600000*h.
 * 
 * The return value does not depend on DST changes. For example, 10 am will always return 36000000, even if there was a
 * DST time change earlier in the day. In these cases, getTimeOfDay() is not the same as "milliseconds since midnight".
 */
export function getTimeOfDay(unix: number, zone: TimeChange[]): number {
    const offset = getOffsetFromTable(unix, zone);
    return mf.mod(unix + offset * 60000, DAY_LENGTH); // Unix time is in UTC, adding the offset converts to local time zone
}

/**
 * Converts the sun event to a string (printable using console.log or process.stdout.write).
 * @param event The event as an SEvent object, with the Unix timestamp and event type.
 * @param zoneTable Time zone lookup table.
 * @param twentyFourHours Whether to print time in 12-hour or 24-hour format.
 * @returns A string representing the solar event, printable using console.log or process.stdout.write.
 */
export function sunEventString(event: SEvent, zoneTable: TimeChange[], twentyFourHours = false): string {
    const eventType = event.type.padStart(14);
    const timeOfDay = Math.floor(getTimeOfDay(event.unix, zoneTable));
    const timeString = mf.convertToHMS(timeOfDay, twentyFourHours);
    const elevStr = (mf.refract(event.elev) / degToRad).toFixed(4) + "°";
    const azStr = (event.azimuth / degToRad).toFixed(4) + "° " + mf.direction(event.azimuth).padStart(3);
    const eventStr = `${eventType} | ${timeString.padStart(11)} | ${elevStr.padStart(9)} | ${azStr.padStart(13)}`;

    const bold = event.type === "Sunrise" || event.type === "Sunset" || event.type === "Solar Noon";
    let [r, g, b] = [128, 128, 128];
    if (event.type === "Sunrise" || event.type === "Sunset") {[r, g, b] = [255, 255, 0];}
    else if (event.elev >= HORIZON) {[r, g, b] = [255, 255, 255];}
    const colorStr = `\x1b[38;2;${r};${g};${b}m`;
    const boldStr = "\x1b[1m";
    const resetStr = "\x1b[0m";
    if (bold) {return colorStr + boldStr + eventStr + resetStr;}
    else {return colorStr + eventStr + resetStr;}
}

/**
 * Converts the moon event to a string (printable using console.log or process.stdout.write).
 * @param event The event as an SEvent object, with the Unix timestamp and event type.
 * @param zoneTable Time zone lookup table.
 * @param twentyFourHours Whether to print time in 12-hour or 24-hour format.
 * @returns A string representing the solar event, printable using console.log or process.stdout.write.
 */
export function moonEventString(event: SEvent, zoneTable: TimeChange[], twentyFourHours = false): string {
    const eventType = event.type.padStart(16);
    const timeOfDay = Math.floor(getTimeOfDay(event.unix, zoneTable));
    const timeString = mf.convertToHMS(timeOfDay, twentyFourHours);
    const elevStr = (mf.refract(event.elev) / degToRad).toFixed(4) + "°";
    const azStr = (event.azimuth / degToRad).toFixed(4) + "° " + mf.direction(event.azimuth).padStart(3);
    const illum = illumination(event.unix);
    const illumStr = (100*illum).toFixed(2) + "%";
    const eventStr = `${eventType} | ${timeString.padStart(11)} | ${elevStr.padStart(9)} | ${azStr.padStart(13)} | ${illumStr.padStart(12)}`;

    const bold = event.type === "Moonrise" || event.type === "Moonset" || event.type === "Meridian Passing";
    let [r, g, b] = [128, 128, 128];
    if (event.type === "Moonrise" || event.type === "Moonset") {[r, g, b] = [255, 255, 0];}
    else if (event.elev >= HORIZON) {[r, g, b] = [255, 255, 255];}
    const colorStr = `\x1b[38;2;${r};${g};${b}m`;
    const boldStr = "\x1b[1m";
    const resetStr = "\x1b[0m";
    if (bold) {return colorStr + boldStr + eventStr + resetStr;}
    else {return colorStr + eventStr + resetStr;}
}