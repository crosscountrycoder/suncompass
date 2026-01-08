import {sunEventsDay, allSunEvents, intervals, lengths, generateSunTable} from "../src/core/suncalc.ts";
import {find} from "geo-tz";
import {generateMoonSvg, generateSunSvg} from "../src/core/gen-svg.ts";
import { DateTime } from "luxon";
import fs from "fs";
import * as mf from "../src/core/mathfuncs.ts";
import { generateMoonTable } from "../src/core/mooncalc.ts";

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
const year = args.length >= 5 ? Number(args[4]) : DateTime.now().setZone(timeZone).year;

const daylengthFileName = 
args.length === 6 ? `./diagrams/${args[5]}-day-lengths-${year}.svg` : `./diagrams/day-lengths.svg`;
const risesetFileName = 
args.length === 6 ? `./diagrams/${args[5]}-sunrise-sunset-${year}.svg` : `./diagrams/sunrise-sunset.svg`;
const moonFileName =
args.length === 6 ? `./diagrams/${args[5]}-moon-${year}.svg` : `./diagrams/moon.svg`;

const sunTable = generateSunTable(lat, long, year, timeZone);

// Generate day length SVG
const daylengthSvg = generateSunSvg({sunTable: sunTable, type: "length"});
fs.writeFileSync(daylengthFileName, daylengthSvg, "utf8");
console.log(`File written to ${daylengthFileName}`);

// Generate sunrise-sunset SVG
const risesetSvg = generateSunSvg({sunTable: sunTable, type: "rise-set"});
fs.writeFileSync(risesetFileName, risesetSvg, "utf8");
console.log(`File written to ${risesetFileName}`);

// Do moon calculations
const moonTable = generateMoonTable(lat, long, year, timeZone);
const moonSvg = generateMoonSvg({sunTable: sunTable, moonTable: moonTable});
fs.writeFileSync(moonFileName, moonSvg, "utf8");
console.log(`File written to ${moonFileName}`);

const end = performance.now();
console.log(`Took ${((end-start)/1000).toFixed(3)} seconds`);