import { DateTime, Duration } from "luxon";
import {find} from "geo-tz";
import * as sc from "../src/core/suncalc.ts";
import * as mf from "../src/core/mathfuncs.ts";
import { generateLODProfile, timeZoneLookupTable, sunEventString } from "../src/core/lookup-tables.ts";
import fs from "fs";
import path from "path";

type RawSeasonRecord = {year: number; marEquinox: string; junSolstice: string; sepEquinox: string; decSolstice: string;};
function solsticesEquinoxes(): RawSeasonRecord[] {
  const jsonPath = path.join("public", "data", "solstices_equinoxes.json");
  const raw = fs.readFileSync(jsonPath, "utf8");
  return JSON.parse(raw) as RawSeasonRecord[];
}

const args = process.argv;
let lat: number, long: number, zone: string, date: DateTime | undefined;
if (args.length === 2) {
    [lat, long] = [34.42,-119.85]; // the location around the University of California, Santa Barbara
    zone = find(lat, long)[0];
    date = DateTime.now().setZone(zone);   
}
else if (args.length === 4) { 
    /* Accepts coordinates. Example: "npx ts-node scripts/rise-set.ts 40.75 -73.99" gives rise-set times for Manhattan, 
    New York City in Eastern Time. */
    [lat, long] = [Number(args[2]), Number(args[3])];
    zone = find(lat, long)[0];
    date = DateTime.now().setZone(zone);
}
else if (args.length === 5) {
    /* Coordinates and date. Example: "npx ts-node scripts/rise-set.ts 40.75 -73.99 2025-06-20" gives times for June 20, 2025 
    in Manhattan in EDT.
    The date argument can be replaced with "me" for March equinox, "js" for June solstice, "se" for Sep equinox or "ds" 
    for Dec solstice.
    The specific time can also be specified in the date. Example:
    "2025-06-21T02:04:07" -> June 21, 2025, 2:42:07 am
    "2025-12-21T15:02:45" -> December 21, 2025, 3:02:45 pm
    */
    [lat, long] = [Number(args[2]), Number(args[3])];
    zone = find(lat, long)[0];
    const curYear = DateTime.now().setZone(zone).year;
    const solstices = solsticesEquinoxes();
    if (args[4] === "me") {date = DateTime.fromISO(solstices[curYear].marEquinox, {zone: zone});}
    else if (args[4] === "js") {date = DateTime.fromISO(solstices[curYear].junSolstice, {zone: zone});}
    else if (args[4] === "se") {date = DateTime.fromISO(solstices[curYear].sepEquinox, {zone: zone});}
    else if (args[4] === "ds") {date = DateTime.fromISO(solstices[curYear].decSolstice, {zone: zone});}
    else {date = DateTime.fromISO(args[4], {zone: zone});}
}
else {
    console.log("Invalid argument");
    process.exit(1);
}

if (Math.abs(lat) >= 90) {console.log("Latitude must be between -90 and 90, exclusive (use ±89.9999 for poles)");}
else if (Math.abs(long) > 180) {console.log("Longitude must be between -180 and 180");}
else {
    const lod = generateLODProfile(mf.ms(date));
    const subsolarPoint = sc.subsolarPoint(lod);
    const [elev, az] = sc.sunPosition(lat, long, lod);
    const apparentElev = mf.refract(elev);
    const dist = lod.distance;

    const dayStarts = mf.dayStarts(DateTime.fromObject({year: date.year}, {zone: date.zone}), 
    DateTime.fromObject({year: date.year+1}, {zone: date.zone}));
    const timeZoneTable = timeZoneLookupTable(dayStarts);

    console.log(zone);
    console.log(`${date.toLocaleString(DateTime.DATETIME_FULL_WITH_SECONDS)}`);
    console.log(`Current sun elevation: ${elev.toFixed(4)}° (After refraction: ${apparentElev.toFixed(4)}°)`);
    console.log(`Current sun bearing: ${az.toFixed(4)}° (${mf.direction(az)})`);
    console.log(`Subsolar point: ${subsolarPoint[0].toFixed(4)}, ${subsolarPoint[1].toFixed(4)}`);
    console.log(`Sun-earth distance: ${dist.toFixed(0)} km (${(dist/1.609344).toFixed(0)} mi)`);

    // Print day length
    const solarEvents = sc.sunEventsDay(lat, long, date);
    const solarEventsY = sc.sunEventsDay(lat, long, date.minus({days: 1}));
    const solarEventsT = sc.sunEventsDay(lat, long, date.plus({days: 1}));
    const curDayStart = mf.ms(date.startOf("day"));
    const dayLength = sc.dayLength(curDayStart, solarEventsY, solarEvents, solarEventsT);
    if (dayLength === -1) {console.log("Day length: undefined");}
    else {console.log(`Day length: ${Duration.fromMillis(1000*Math.round(dayLength/1000)).toFormat("h:mm:ss")}`);}
    console.log();

    // Print sunrise, sunset, twilight, solar noon, and solar midnight times
    console.log("         Event |        Time | Elevation |       Bearing"); // header
    for (const event of solarEvents) {console.log(sunEventString(event, timeZoneTable));}
}