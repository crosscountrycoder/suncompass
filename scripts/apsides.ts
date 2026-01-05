import {DateTime} from "luxon";
import {sunApsides, sunDistance} from "../src/core/suncalc.ts";
import {moonApsides, moonDistance} from "../src/core/mooncalc.ts";
import {ms} from "../src/core/mathfuncs.ts";

let year: number;
let zone: string;

const args = process.argv;
if (args.length === 2) {year = DateTime.now().year;}
else {year = Number(args[2]);}
if (args.length <= 3) {zone = "local";}
else {zone = args[3];}

const startYear = ms(DateTime.fromObject({year: year}, {zone: zone}));
const endYear = ms(DateTime.fromObject({year: year + 1}, {zone: zone}));
const apsidesSun = sunApsides(startYear, endYear);
const apsidesMoon = moonApsides(startYear, endYear);

console.log("Sun apsides:");

let perihelionStr = "Perihelion: ";
for (const x of apsidesSun.perihelion) {
    const distance = sunDistance(x, true);
    perihelionStr += (DateTime.fromMillis(x, {zone: zone}).toFormat("MMM d, y HH:mm:ss ZZZZ") + 
    `\x1b[38;2;128;128;128m` +
    ` (${Math.round(distance)} km / ${Math.round(distance/1.609344)} mi)\x1b[0m, `);
}
perihelionStr = perihelionStr.slice(0, -2);
console.log(perihelionStr);

let aphelionStr = "Aphelion: ";
for (const x of apsidesSun.aphelion) {
    const distance = sunDistance(x, true);
    aphelionStr += (DateTime.fromMillis(x, {zone: zone}).toFormat("MMM d, y HH:mm:ss ZZZZ") + 
    `\x1b[38;2;128;128;128m` +
    ` (${Math.round(distance)} km / ${Math.round(distance/1.609344)} mi)\x1b[0m, `);
}
aphelionStr = aphelionStr.slice(0, -2);
console.log(aphelionStr);

console.log();
console.log("Moon apsides:")

console.log("Perigees:");
for (const x of apsidesMoon.perigees) {
    const distance = moonDistance(x, true);
    console.log(DateTime.fromMillis(x, {zone: zone}).toFormat("MMM d, y HH:mm:ss ZZZZ") +
    ` \x1b[38;2;128;128;128m(${Math.round(distance)} km / ${Math.round(distance/1.609344)} mi)\x1b[0m`);
}
console.log("Apogees:");
for (const x of apsidesMoon.apogees) {
    const distance = moonDistance(x, true);
    console.log(DateTime.fromMillis(x, {zone: zone}).toFormat("MMM d, y HH:mm:ss ZZZZ") +
    ` \x1b[38;2;128;128;128m(${Math.round(distance)} km / ${Math.round(distance/1.609344)} mi)\x1b[0m`);
}