import * as mf from "../src/core/mathfuncs.ts";
import * as sc from "../src/core/suncalc.ts";
import * as mc from "../src/core/mooncalc.ts";
import * as gwm from "../src/core/world-map.ts";
import { writeFileSync } from "node:fs";
import { DateTime } from "luxon";
import { generateLODProfile } from "../src/core/lookup-tables.ts";

const start = performance.now();
const args = process.argv;

if (args.length < 2 || args.length > 4) {
    console.log("Syntax: npx ts-node scripts/test.ts [time] [timeZone]");
    process.exit(1);
}

const unix: number = (() => {
    if (args.length === 2) return Date.now();
    else if (args.length === 3) {
        const date = new Date();
        const solstEqMonth: Record<string, 3 | 6 | 9 | 12> = { me: 3, js: 6, se: 9, ds: 12 };
        const month: 3 | 6 | 9 | 12 = solstEqMonth[args[2]];
        return month ? sc.calcSolstEq(date.getFullYear(), month) : mf.ms(DateTime.fromISO(args[2]));
    }
    else {return mf.ms(DateTime.fromISO(args[2], {zone: args[3]}));}
})();

const svg = gwm.plotSvg(unix); // example args
writeFileSync("./diagrams/world-map.svg", svg, "utf8");

console.log(DateTime.fromMillis(unix, {zone: "utc"}).toFormat("MMM d, y HH:mm:ss.SSS ZZZZ"));
const [sLat, sLong] = sc.subsolarPoint(generateLODProfile(unix), true);
const [mLat, mLong] = mc.sublunarPoint(unix, true);
console.log(`Subsolar point: ${sLat.toFixed(4)}, ${sLong.toFixed(4)}`);
console.log(`Sublunar point: ${mLat.toFixed(4)}, ${mLong.toFixed(4)}`);
console.log("Wrote world-map.svg");
const end = performance.now();
console.log(`Took ${((end-start)/1000).toFixed(3)} seconds`);