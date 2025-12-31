import {DateTime} from "luxon";
import {sunApsides} from "../src/core/suncalc.ts";
import {moonApsides} from "../src/core/mooncalc.ts";
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
    perihelionStr += (DateTime.fromMillis(x, {zone: zone}).toFormat("MMM d, y HH:mm:ss ZZZZ") + ", ");
}
perihelionStr = perihelionStr.slice(0, -2);
console.log(perihelionStr);

let aphelionStr = "Aphelion: ";
for (const x of apsidesSun.aphelion) {
    aphelionStr += (DateTime.fromMillis(x, {zone: zone}).toFormat("MMM d, y HH:mm:ss ZZZZ") + ", ");
}
aphelionStr = aphelionStr.slice(0, -2);
console.log(aphelionStr);

console.log();
console.log("Moon apsides:")

console.log("Perigees:");
for (const x of apsidesMoon.perigees) {
    console.log(DateTime.fromMillis(x, {zone: zone}).toFormat("MMM d, y HH:mm:ss ZZZZ"));
}
console.log("Apogees:");
for (const x of apsidesMoon.apogees) {
    console.log(DateTime.fromMillis(x, {zone: zone}).toFormat("MMM d, y HH:mm:ss ZZZZ"));
}