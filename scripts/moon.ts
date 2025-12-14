import * as mc from "../src/core/mooncalc.ts";
import {DateTime} from "luxon";
import * as mf from "../src/core/mathfuncs.ts";
import {find} from "geo-tz";
import {timeZoneLookupTable, moonEventString} from "../src/core/lookup-tables.ts";

const args = process.argv;
let lat: number, long: number, zone: string, date: DateTime | undefined;
if (args.length === 4) { 
    /* Accepts coordinates. Example: "npx ts-node scripts/moon.ts 40.75 -73.99" gives moonrise and moonset times for 
    Manhattan, New York City in Eastern Time. */
    [lat, long] = [Number(args[2]), Number(args[3])];
    zone = find(lat, long)[0];
    date = DateTime.now().setZone(zone);
}
else if (args.length === 5) {
    /* Coordinates and date. Example: "npx ts-node scripts/moon.ts 40.75 -73.99 2025-06-20" gives times for June 20, 2025
    in Manhattan in EDT.
    The specific time can also be specified in the date. Example:
    "2025-06-21T02:04:07" -> June 21, 2025, 2:42:07 am
    "2025-12-21T15:02:45" -> December 21, 2025, 3:02:45 pm
    */
    [lat, long] = [Number(args[2]), Number(args[3])];
    zone = find(lat, long)[0];
    date = DateTime.fromISO(args[4], {zone: zone});
}
else {
    console.log("Syntax:")
    console.log("\x1b[1mnpx ts-node scripts/moon.ts <latitude> <longitude> [time]\x1b[0m");
    console.log("Latitude and longitude are expressed as decimals, in the ranges -90 <= lat <= 90 and -180 <= long <= 180.");
    console.log("Time is expressed as an ISO timestamp, ex: 2025-06-20 (start of day on June 20, 2025) or 2025-12-21T15:02:45 "
        + "(Dec. 21, 2025 at 15:02:45 local time)");
    process.exit(1);
}

if (Math.abs(lat) >= 90) {console.log("Latitude must be between -90 and 90, exclusive (use ±89.9999 for poles)");}
else if (Math.abs(long) > 180) {console.log("Longitude must be between -180 and 180");}
else {
    const dayStarts = mf.dayStarts(DateTime.fromObject({year: date.year}, {zone: date.zone}), 
        DateTime.fromObject({year: date.year+1}, {zone: date.zone}));
    const timeZoneTable = timeZoneLookupTable(dayStarts);
    const [elev, az] = mc.moonPosition(lat, long, mf.ms(date));
    const [sublunarLat, sublunarLong] = mc.sublunarPoint(mf.ms(date));
    const apparentElev = mf.refract(elev);
    const distance = mc.moonDistance(mf.ms(date), true);
    const illumination = mc.illumination(mf.ms(date));
    const phase = mc.moonPhaseDay(mf.ms(date.startOf("day")), mf.ms(date.startOf("day").plus({days: 1})));
    
    console.log(zone);
    console.log(`${date.toLocaleString(DateTime.DATETIME_FULL_WITH_SECONDS)}`);
    console.log(`Current moon elevation: ${elev.toFixed(4)}° (After refraction: ${apparentElev.toFixed(4)}°)`);
    console.log(`Current moon bearing: ${az.toFixed(4)}° (${mf.direction(az)})`);
    console.log(`Sublunar point: ${sublunarLat.toFixed(4)}, ${sublunarLong.toFixed(4)}`);
    console.log(`Moon-Earth distance: ${Math.round(distance)} km (${Math.round(distance/1.609344)} mi)`);
    console.log(`Illumination: ${(100*illumination).toFixed(2)}%`);
    console.log(`Phase: ${phase}`);
    console.log();

    const lunarEvents = mc.moonEventsDay(lat, long, date);
    console.log("           Event |        Time | Elevation |       Bearing | Illumination"); // header
    for (const event of lunarEvents) {console.log(moonEventString(event, timeZoneTable));}
}