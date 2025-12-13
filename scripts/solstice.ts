// Returns the time of solstices and equinoxes for a given year and time zone. If the year and/or time zone are unspecified, it defaults to the current year and local time zone of the computer.
// Example queries:
// "node solstice.js" returns solstices and equinoxes for current year and local time zone.
// "node solstice.js 2025" returns solstices and equinoxes for 2025 in the local time zone.
// "node solstice.js 2026 utc" returns solstices and equinoxes for 2026 in UTC.
// "node solstice.js 2025 America/Los_Angeles" returns solstices and equinoxes for 2025 in Pacific Time.

import {DateTime} from "luxon";
import {subsolarPoint} from "../src/core/suncalc.ts";
import { generateLODProfile } from "../src/core/lookup-tables.ts";
import {ms} from "../src/core/mathfuncs.ts";
import fs from "fs";
import path from "path";

let year: number;
let zone: string;

type RawSeasonRecord = {year: number; marEquinox: string; junSolstice: string; sepEquinox: string; decSolstice: string;};
function solsticesEquinoxes(): RawSeasonRecord[] {
  const jsonPath = path.join("public", "data", "solstices_equinoxes.json");
  const raw = fs.readFileSync(jsonPath, "utf8");
  return JSON.parse(raw) as RawSeasonRecord[];
}

const args = process.argv;
if (args.length === 2) {year = DateTime.now().year;}
else {year = Number(args[2]);}
if (args.length <= 3) {zone = "local";}
else {zone = args[3];}

const obj = solsticesEquinoxes()[year];
const [mar, jun, sep, dec] = [
    DateTime.fromISO(obj.marEquinox, {zone: zone}), 
    DateTime.fromISO(obj.junSolstice, {zone: zone}), 
    DateTime.fromISO(obj.sepEquinox, {zone: zone}), 
    DateTime.fromISO(obj.decSolstice, {zone: zone})
];
const [marSSP, junSSP, sepSSP, decSSP] = [
    subsolarPoint(generateLODProfile(ms(mar))),
    subsolarPoint(generateLODProfile(ms(jun))),
    subsolarPoint(generateLODProfile(ms(sep))),
    subsolarPoint(generateLODProfile(ms(dec)))
];

console.log("March equinox: " + mar.toFormat("MMM d, y HH:mm:ss ZZZZ"));
console.log("Subsolar point: " + marSSP[0].toFixed(4) + ", " + marSSP[1].toFixed(4));
console.log();

console.log("June solstice: " + jun.toFormat("MMM d, y HH:mm:ss ZZZZ"));
console.log("Subsolar point: " + junSSP[0].toFixed(4) + ", " + junSSP[1].toFixed(4));
console.log();

console.log("September equinox: " + sep.toFormat("MMM d, y HH:mm:ss ZZZZ"));
console.log("Subsolar point: " + sepSSP[0].toFixed(4) + ", " + sepSSP[1].toFixed(4));
console.log();

console.log("December solstice: " + dec.toFormat("MMM d, y HH:mm:ss ZZZZ"));
console.log("Subsolar point: " + decSSP[0].toFixed(4) + ", " + decSSP[1].toFixed(4));