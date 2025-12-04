import * as fs from "node:fs";
import {calcSolstEq} from "../src/core/suncalc.ts";
import {DateTime} from "luxon";

const data = [];

const [startYear, endYear] = [0, 2500];
for (let year=startYear; year<=endYear; year++) {
    const [me, js, se, ds] = [
        DateTime.fromMillis(calcSolstEq(year, 3), {zone: "utc"}).toISO(), 
        DateTime.fromMillis(calcSolstEq(year, 6), {zone: "utc"}).toISO(), 
        DateTime.fromMillis(calcSolstEq(year, 9), {zone: "utc"}).toISO(), 
        DateTime.fromMillis(calcSolstEq(year, 12), {zone: "utc"}).toISO()
    ];
    const obj = {year: year, marEquinox: me, junSolstice: js, sepEquinox: se, decSolstice: ds};
    data.push(obj);
}

const jsonContent = JSON.stringify(data, null, 2);
fs.writeFile("public/data/solstices_equinoxes.json", jsonContent, "utf8", function (err) {
    if (err) {
        console.log("Error occurred - file not saved");
        throw err;
    }
    else {console.log("File saved successfully");}
});