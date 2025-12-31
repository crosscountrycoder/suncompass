import {DateTime} from "luxon";
import {subsolarPoint, calcSolstEq} from "../src/core/suncalc.ts";
import { generateLODProfile } from "../src/core/lookup-tables.ts";

let year: number;
let zone: string;

const args = process.argv;
if (args.length === 2) {year = DateTime.now().year;}
else {year = Number(args[2]);}
if (args.length <= 3) {zone = "local";}
else {zone = args[3];}

const mar = calcSolstEq(year, 3);
const jun = calcSolstEq(year, 6);
const sep = calcSolstEq(year, 9);
const dec = calcSolstEq(year, 12);

const [marSSP, junSSP, sepSSP, decSSP] = [
    subsolarPoint(generateLODProfile(mar)),
    subsolarPoint(generateLODProfile(jun)),
    subsolarPoint(generateLODProfile(sep)),
    subsolarPoint(generateLODProfile(dec))
];

const marDT = DateTime.fromMillis(mar, {zone: zone});
const junDT = DateTime.fromMillis(jun, {zone: zone});
const sepDT = DateTime.fromMillis(sep, {zone: zone});
const decDT = DateTime.fromMillis(dec, {zone: zone});

console.log("March equinox: " + marDT.toFormat("MMM d, y HH:mm:ss ZZZZ"));
console.log("Subsolar point: " + marSSP[0].toFixed(4) + ", " + marSSP[1].toFixed(4));
console.log();

console.log("June solstice: " + junDT.toFormat("MMM d, y HH:mm:ss ZZZZ"));
console.log("Subsolar point: " + junSSP[0].toFixed(4) + ", " + junSSP[1].toFixed(4));
console.log();

console.log("September equinox: " + sepDT.toFormat("MMM d, y HH:mm:ss ZZZZ"));
console.log("Subsolar point: " + sepSSP[0].toFixed(4) + ", " + sepSSP[1].toFixed(4));
console.log();

console.log("December solstice: " + decDT.toFormat("MMM d, y HH:mm:ss ZZZZ"));
console.log("Subsolar point: " + decSSP[0].toFixed(4) + ", " + decSSP[1].toFixed(4));