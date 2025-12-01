import {sunEventsDay, allSunEvents, intervals, lengths} from "../src/core/suncalc.ts";
import {find} from "geo-tz";
import {generateMoonSvg, generateSunSvg} from "../src/core/gen-svg.ts";
import { DateTime } from "luxon";
import fs from "fs";
import * as mf from "../src/core/mathfuncs.ts";
import { timeZoneLookupTable, longDistLookupTable } from "../src/core/lookup-tables.ts";
import path from "path";
import { allMoonEvents, moonIntervals, moonPhase } from "../src/core/mooncalc.ts";

type RawSeasonRecord = {year: number; marEquinox: number; junSolstice: number; sepEquinox: number; decSolstice: number;};
function solsticesEquinoxes(): RawSeasonRecord[] {
    const jsonPath = path.join("public", "data", "solstices_equinoxes.json");
    const raw = fs.readFileSync(jsonPath, "utf8");
    return JSON.parse(raw) as RawSeasonRecord[];
}

const start = performance.now();

const args = process.argv;
if (args.length < 4 || args.length > 6) {
    console.log("Syntax: npx ts-node scripts/generate-svg.ts <lat> <long> [year] [location-name]");
    process.exit(1);
}

const [lat, long] = [Number(args[2]), Number(args[3])];
if (Math.abs(lat) >= 90) {
    console.log("Latitude must be between -90 and 90, exclusive (use ±89.9999 for poles)");
    process.exit(1);
}
else if (Math.abs(long) > 180) {
    console.log("Longitude must be between -180 and 180");
    process.exit(1);
}
const timeZone = find(lat, long)[0];

/** If year is specified, set the year to what is given. Otherwise, set to current year in given location. */
const year = (args.length >= 5) ? Number(args[4]) : DateTime.now().setZone(timeZone).year;

const daylengthFileName = 
(args.length == 6) ? `./diagrams/${args[5]}-day-lengths-${year}.svg` : `./diagrams/day-lengths.svg`;
const risesetFileName = 
(args.length == 6) ? `./diagrams/${args[5]}-sunrise-sunset-${year}.svg` : `./diagrams/sunrise-sunset.svg`;
const moonFileName =
(args.length == 6) ? `./diagrams/${args[5]}-moon-${year}.svg` : `./diagrams/moon.svg`;

const sunEvents = [];
const startDate = DateTime.fromObject({year: year}, {zone: timeZone});
const endDate = DateTime.fromObject({year: year+1}, {zone: timeZone});
const dateList = mf.dayStarts(startDate, endDate);
const tzLookupTable = timeZoneLookupTable(dateList);
const lodLookupTable = longDistLookupTable(dateList);

for (let i=0; i<dateList.length-1; i++) {
    const startLOD = lodLookupTable[i], endLOD = lodLookupTable[i+1];
    const curDaySunEvents = allSunEvents(lat, long, startLOD, endLOD);
    sunEvents.push(curDaySunEvents);
}

const solstEq = solsticesEquinoxes();
const solstEqDT = [
    DateTime.fromMillis(solstEq[year].marEquinox, {zone: timeZone}),
    DateTime.fromMillis(solstEq[year].junSolstice, {zone: timeZone}),
    DateTime.fromMillis(solstEq[year].sepEquinox, {zone: timeZone}),
    DateTime.fromMillis(solstEq[year].decSolstice, {zone: timeZone}),
]

// Generate day length SVG
const daylengthSvg = generateSunSvg({events: sunEvents, type: "length", timeZone: tzLookupTable, solsticesEquinoxes: solstEqDT});
fs.writeFileSync(daylengthFileName, daylengthSvg, "utf8");
console.log(`File written to ${daylengthFileName}`);

// Generate sunrise-sunset SVG
const risesetSvg = generateSunSvg({events: sunEvents, type: "rise-set", timeZone: tzLookupTable, solsticesEquinoxes: solstEqDT});
fs.writeFileSync(risesetFileName, risesetSvg, "utf8");
console.log(`File written to ${risesetFileName}`);

// Do moon calculations
const moonEvents = [];
const moonInts = [];
const newMoons = [];
const fullMoons = [];
for (let i=0; i<dateList.length-1; i++) {
    const startUnix = mf.ms(dateList[i]), endUnix = mf.ms(dateList[i+1]);
    const curDayMoonEvents = allMoonEvents(lat, long, startUnix, endUnix);
    const curDayIntervals = moonIntervals(lat, long, startUnix, curDayMoonEvents, tzLookupTable);
    const phase = moonPhase(lat, long, startUnix, endUnix);
    if (phase !== null) {
        if (phase.type == "New Moon") {newMoons.push(DateTime.fromMillis(phase.unix, {zone: timeZone}));}
        else if (phase.type == "Full Moon") {fullMoons.push(DateTime.fromMillis(phase.unix, {zone: timeZone}));}
    }
    moonEvents.push(curDayMoonEvents);
    moonInts.push(curDayIntervals);
}

const moonSvg = generateMoonSvg({sunEvents: sunEvents, moonIntervals: moonInts, timeZone: tzLookupTable, newMoons: newMoons,
    fullMoons: fullMoons
});
fs.writeFileSync(moonFileName, moonSvg, "utf8");
console.log(`File written to ${moonFileName}`);

const end = performance.now();
console.log(`Took ${((end-start)/1000).toFixed(3)} seconds`);