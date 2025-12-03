/**
 * This TypeScript file contains code to generate SVG chart of day length and sunrise/sunset/twilight times for an entire year.
 * "generateSvg" generates an SVG file showing (from top to bottom) the lengths of day, civil twilight, nautical twilight, astronomical
 * twilight and night for an entire year
 */

import {getTimeOfDay} from "./lookup-tables.ts";
import * as mf from "./mathfuncs.ts";
import {intervalsSvg, lengths, intervalsNightCivilTwilight} from "./suncalc.ts";
import {DAY_LENGTH} from "./constants.ts";
import {DateTime} from "luxon";
import type {SEvent, TimeChange} from "./lookup-tables.ts";

const svgClose = "</svg>";
const sunColors = ["#80c0ff", "#0060c0", "#004080", "#002040", "#000000"];

type sunSvgOptions = {
    events: SEvent[][];
    type: string;
    timeZone: TimeChange[];
    solsticesEquinoxes: DateTime[];
    svgWidth?: number;
    svgHeight?: number;
    leftPadding?: number;
    rightPadding?: number;
    topPadding?: number;
    bottomPadding?: number;
    textSize?: number;
    font?: string;
    language?: string;
    gridInterval?: number;
    gridlineWidth?: number
};

type moonSvgOptions = {
    sunEvents: SEvent[][];
    moonIntervals: number[][][];
    timeZone: TimeChange[];
    newMoons: DateTime[];
    fullMoons: DateTime[];
    svgWidth?: number;
    svgHeight?: number;
    leftPadding?: number;
    rightPadding?: number;
    topPadding?: number;
    bottomPadding?: number;
    textSize?: number;
    font?: string;
    language?: string;
    gridInterval?: number;
    gridlineWidth?: number
};

/** Generates the opening of an SVG */
function svgOpen(width: number, height: number): string {
    return `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">\n`
}

/**
 * Generates SVG code for a polygon from an array of points, with the specified fill color, stroke color and stroke width
 * @param points An array of [x, y] points, ex: [[0, 0], [0, 100], [100, 0]]
 * @param fillColor Fill color of polygon.
 * @param strokeColor Stroke color of polygon.
 * @param strokeWidth Stroke width of polygon.
 * @param precision Number of digits after the decimal point to round pixel coordinates.
 * @returns SVG string for the given polygon.
 */
function polygonFromArray(
    points: mf.Polygon,
    fillColor: string = "none",
    strokeColor: string = "none",
    strokeWidth: number = 0,
    precision: number = 2,
): string {
    const simplifiedPoints = simplifyCollinear(points); 
    const ptsAttr = simplifiedPoints.map(([x,y]) => `${mf.toFixedS(x,precision)},${mf.toFixedS(y,precision)}`).join(" "); // format the "x,y x,y ..." string
    return `<polygon points="${ptsAttr}" fill="${fillColor}" stroke="${strokeColor}" stroke-width="${strokeWidth}"/>\n`;
}

/** Generates SVG code for a polyline from an array of points with the specified stroke color, width, and precision (digits after the
 * decimal point in the coordinates). */
function polylineFromArray(points: mf.Polyline, color: string = "#000000", width: number = 1, precision: number = 2): string {
    const ptsAttr = points.map(([x,y]) => `${mf.toFixedS(x,precision)},${mf.toFixedS(y,precision)}`).join(" "); // format the "x,y x,y ..." string
    return `<polyline points="${ptsAttr}" fill="none" stroke="${color}" stroke-width="${width}"/>\n`;
}

/** Generates SVG code for a rectangle with the top-left corner at the given x and y cordinates, and the given width, height,
 * fill and stroke colors. */
function rectangleSvg(x: number, y: number, width: number, height: number, fillColor: string = "none", strokeColor: string = "none",
    strokeWidth: number = 0) {
    return `<rect x="${x}" y="${y}" width="${width}" height="${height}" fill="${fillColor}" stroke="${strokeColor}" stroke-width="${strokeWidth}"/>\n`
}

/** Generates SVG code for a text box with the given text, centered at the given x and y coordinate, with the given font and font size. 
 * Text anchor can be "start" (left-aligned), "middle" (centered), or "end" (right-aligned). Alignment baseline can be either 
 * "text-before-edge" (top-aligned), "middle" (centered), or "text-after-edge" (bottom-aligned). */
function textSvg(
    text: string, 
    x: number, 
    y: number, 
    fontSize: number, 
    font: string, 
    textColor: string, 
    textAnchor: string,
    alignmentBaseline: string,
    precision: number = 2
): string {
    return `<text x="${mf.toFixedS(x,precision)}" y="${mf.toFixedS(y,precision)}" font-family="${font}" font-size="${fontSize}pt"`
    + ` text-anchor="${textAnchor}" alignment-baseline="${alignmentBaseline}" fill="${textColor}">${text}</text>\n`;
}

/** Generates an SVG line from (x1, y1) to (x2, y2) with the given color and width. */
function lineSvg(p1: mf.Point, p2: mf.Point, color: string, width: number, precision: number = 2): string {
    return `<line x1="${mf.toFixedS(p1[0],precision)}" y1="${mf.toFixedS(p1[1],precision)}" x2="${mf.toFixedS(p2[0],precision)}"`
    + ` y2="${mf.toFixedS(p2[1],precision)}" stroke="${color}" stroke-width="${width}"/>\n`;
}

/** Returns an array of month abbreviations in the given language, represented by a language code, such as "en" (English), "es"
 * (Spanish), "fr" (French), "zh" (Chinese). So far there is only English - I plan to add more when I localize the site. */
function months(language: string = "en"): string[] {
    return ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
}

/** Edges of the months, used for drawing gridlines. */
function monthEdges(leapYear: boolean = false): number[] {
    if (leapYear) {return [0,31,60,91,121,152,182,213,244,274,305,335,366];}
    else {return [0,31,59,90,120,151,181,212,243,273,304,334,365];}
}

/** Simplifies a polygon or polyline (represented as points) to remove collinear points. */
function simplifyCollinear(points: mf.Polygon | mf.Polyline) {
    if (points.length <= 2) {return points;}
    const newPoints = [points[0], points[1]];
    for (let i=2; i<points.length; i++) {
        const prev = newPoints[newPoints.length-2];
        const cur = newPoints[newPoints.length-1];
        const next = points[i];
        if (mf.isCollinear(prev, cur, next)) {newPoints[newPoints.length-1] = next;}
        else {newPoints.push(next);}
    }
    return newPoints;
}

/** Transforms times of day and year into SVG diagram coordinates based on the image parameters. 
 * The point is transformed in place, and this function doesn't return anything. */
function coordinateTransform(point: mf.Point, options: sunSvgOptions | moonSvgOptions, addOneHalf: boolean = false) {
    const days = ("events" in options) ? options.events.length : options.sunEvents.length;
    if (addOneHalf) {point[0] += 0.5;}
    point[0] = options.leftPadding! + options.svgWidth! * (point[0] / days);
    point[1] = options.topPadding! + options.svgHeight! * (1 - point[1] / DAY_LENGTH);
}

/**
 * Returns a string containing an SVG diagram for either day/twilight/night lengths throughout the year, or the times of day in which
 * day, night, and each stage of twilight occur.
 * Parameters should be passed in an object. All parameters except events, type, timeZone, and solsticesEquinoxes are optional.
 * @param events Values of "allSunEvents" for each day of the year.
 * @param type Set to "length" to generate a day/night/twilight length chart, or "rise-set" to generate a chart with times of day.
 * @param timeZone Time zone, either as an IANA string or a lookup table (see mathfuncs.timeZoneLookupTable)
 * @param solsticesEquinoxes Solstices and equinoxes for the given year, as an array of four DateTimes. They will appear as green lines on the diagram.
 * @param svgWidth Width of the chart (not the entire SVG file) in pixels. Defaults to 1000.
 * @param svgHeight Height of the chart (not the entire SVG file) in pixels. Defaults to 500.
 * @param leftPadding Padding to the left of the carpet plot in pixels. Defaults to 25.
 * @param rightPadding Padding to the right of the carpet plot in pixels. Defaults to 10.
 * @param topPadding Padding above the carpet plot in pixels. Defaults to 10.
 * @param bottomPadding Padding below the carpet plot in pixels. Defaults to 25.
 * @param textSize The font size to use for the axis labels, in points. Defaults to 11.
 * @param font The font family to use for the axis labels. Defaults to Arial.
 * @param language The language used for month abbreviations, represented as a 2-letter code.
 * @param gridInterval Y axis interval. Defaults to 2 (i.e. 2 hours between gridlines).
 * @param gridlineWidth Width of gridlines. Defaults to 0.5 (pixels).
 * @returns A string for the carpet plot, with gridlines and axis labels, that can be saved into an SVG file.
 * The total width of the SVG file is equal to svgWidth + leftPadding + rightPadding. The height is equal to svgHeight + topPadding +
 * bottomPadding.
 */
export function generateSunSvg(options: sunSvgOptions): string {
    const {events,type,timeZone,solsticesEquinoxes,svgWidth=1000,svgHeight=500,leftPadding=25,rightPadding=10,
        topPadding=10,bottomPadding=25,textSize=11,font="Arial",language="en",gridInterval=2,gridlineWidth=0.5} = options;
    const days = events.length; // 365 days for common years, 366 for leap years
    const isLeapYear = (days == 366);
    const nOptions = // normalized options
    {...options,svgWidth,svgHeight,leftPadding,rightPadding,topPadding,bottomPadding,textSize,font,language,gridInterval,gridlineWidth};

    /** Converts a set of durations into an array of points representing a polygon. This includes coordinate transforms. */
    function durationsToArray(durations: number[]): mf.Polygon[] {
        const p: mf.Polygon[] = [[]]; // p is short for polygons
        for (let i=0; i<days; i++) {
            if (durations[i] > 0) {
                if (i == 0 || durations[i-1] == 0) {p[0].push([i, 0], [i, durations[i]], [i+1, durations[i]]);}
                else {p[p.length-1].push([i, durations[i]], [i+1, durations[i]]);}
                if (i == days-1) {p[p.length-1].push([days, 0]);}
            }
            else if (i != 0 && durations[i-1] > 0) {p[p.length-1].push([i, 0]);}
        }
        for (const polygon of p) {
            for (const point of polygon) {coordinateTransform(point, nOptions);}
        }
        return p;
    }

    /** Used to draw lines representing solar noon on the graph. */
    function solarNoonLines() {
        const solarNoons: number[][] = [];
        for (const evts of events) {
            const curDay: number[] = [];
            for (const event of evts) {
                if (event.type == "Solar Noon") {curDay.push(getTimeOfDay(event.unix, timeZone));}
            }
            solarNoons.push(curDay);
        }
        
        const groups: mf.Polyline[] = []; // a group of multiple lines, each representing a cluster of solar noons
        for (const solarNoon of solarNoons[0]) {groups.push([[0, solarNoon]]);}
        for (let i=1; i<days; i++) { // for each day of the year
            for (const noon of solarNoons[i]) { // for each solar noon of the day (may be more than 1)
                let flag: boolean = false;
                for (const group of groups) {
                    if (Math.abs(noon - group[group.length-1][1]) < DAY_LENGTH/48 && group[group.length-1][0] == i-1) {
                        flag = true;
                        group.push([i, noon]);
                        break;
                    }
                }
                if (!flag) {groups.push([[i, noon]]);}
            }
        }
        for (const line of groups) { // convert to SVG coordinates
            for (const point of line) {coordinateTransform(point, nOptions, true);}
        }
        const lines: string[] = [];
        for (const line of groups) {lines.push(polylineFromArray(line, "#ff0000"));}
        return lines;
    }

    /** Used to draw lines representing solar midnight on the graph. */
    function solarMidnightLines() {
        const solarMidnights: number[][] = [];
        for (const evts of events) {
            const curDay: number[] = [];
            for (const event of evts) {
                if (event.type == "Solar Midnight") {curDay.push(getTimeOfDay(event.unix, timeZone));}
            }
            solarMidnights.push(curDay);
        }
        
        const groups: mf.Polyline[] = []; // a group of multiple lines (number[][]), each representing a cluster of solar midnights
        for (const solarMidnight of solarMidnights[0]) {
            groups.push([[0, solarMidnight]]);
        }
        for (let i=1; i<days; i++) { // for each day of the year
            for (const midnight of solarMidnights[i]) { // for each solar midnight of the day (may be more than 1)
                let flag: boolean = false;
                for (const group of groups) {
                    if (Math.abs(midnight - group[group.length-1][1]) < DAY_LENGTH/48 && group[group.length-1][0] == i-1) {
                        flag = true;
                        group.push([i, midnight]);
                        break;
                    }
                }
                if (!flag) {groups.push([[i, midnight]]);}
            }
        }
        for (const line of groups) { // convert to SVG coordinates
            for (const point of line) {coordinateTransform(point, nOptions, true);}
        }
        const lines: string[] = [];
        for (const line of groups) {lines.push(polylineFromArray(line, "#0000ff"));}
        return lines;
    }

    function toPolygons(intervals: number[][][], color: string): string[] {
        const polygons = mf.intervalsToPolygon(intervals);
        for (const polygon of polygons) {
            for (const point of polygon) {coordinateTransform(point, nOptions);}
        }
        const strings: string[] = [];
        for (const polygon of polygons) {strings.push(polygonFromArray(polygon, color));}
        return strings;
    }

    // generate SVG opening and background
    const imageWidth = svgWidth + leftPadding + rightPadding;
    const imageHeight = svgHeight + topPadding + bottomPadding;
    let svgString = svgOpen(imageWidth, imageHeight);
    svgString += rectangleSvg(0, 0, imageWidth, imageHeight, "#ffffff"); // white background

    if (type == "length") { // day/twilight/night length plot
        const dLengths: number[] = []; // day lengths
        const cLengths: number[] = []; // day + civil twilight lengths
        const nLengths: number[] = []; // day + civil + nautical twilight lengths
        const aLengths: number[] = []; // day + civil + nautical + astronomical twilight lengths
        for (const e of events) {
            const dur = lengths(e, timeZone);
            dLengths.push(dur[0]);
            cLengths.push(dur[1]);
            nLengths.push(dur[2]);
            aLengths.push(dur[3]);
        }
        const dp = durationsToArray(dLengths); // day polygons
        const cp = durationsToArray(cLengths); // civil twilight polygons
        const np = durationsToArray(nLengths); // nautical twilight polygons
        const ap = durationsToArray(aLengths); // astronomical twilight polygons
        
        // construct SVG day length diagram
        svgString += rectangleSvg(leftPadding, topPadding, svgWidth, svgHeight, sunColors[4]); // night
        for (const polygon of ap) {svgString += polygonFromArray(polygon, sunColors[3]);} // astronomical twilight
        for (const polygon of np) {svgString += polygonFromArray(polygon, sunColors[2]);} // nautical twilight
        for (const polygon of cp) {svgString += polygonFromArray(polygon, sunColors[1]);} // civil twilight
        for (const polygon of dp) {svgString += polygonFromArray(polygon, sunColors[0]);} // daylight
    }
    else if (type == "rise-set") { // sunrise, sunset, dusk, dawn plot
        const aIntervals: number[][][] = []; // intervals of astronomical twilight or brighter
        const nIntervals: number[][][] = []; // intervals of nautical twilight or brighter
        const cIntervals: number[][][] = []; // intervals of civil twilight or brighter
        const dIntervals: number[][][] = []; // intervals of daylight

        for (const event of events) {
            const int = intervalsSvg(event, timeZone);
            aIntervals.push(int[3]);
            nIntervals.push(int[2]);
            cIntervals.push(int[1]);
            dIntervals.push(int[0]);
        }
        
        const aPolygons = toPolygons(aIntervals, sunColors[3]); // astronomical twilight
        const nPolygons = toPolygons(nIntervals, sunColors[2]); // nautical twilight
        const cPolygons = toPolygons(cIntervals, sunColors[1]); // civil twilight
        const dPolygons = toPolygons(dIntervals, sunColors[0]); // daylight

        const allPolygons = [...aPolygons, ...nPolygons, ...cPolygons, ...dPolygons];

        svgString += rectangleSvg(leftPadding, topPadding, svgWidth, svgHeight, sunColors[4]); // night
        for (const polygon of allPolygons) {svgString += polygon;} // twilight + daylight

        const noonMidnightLines = [...solarMidnightLines(), ...solarNoonLines()];
        for (const line of noonMidnightLines) {svgString += line;}
    }
    
    // draw solstices and equinoxes as green lines
    for (const date of Object.values(solsticesEquinoxes)) {
        const newYear = DateTime.fromISO(`${date.year}-01-01`, {zone: date.zone});
        const x = date.diff(newYear, ['days', 'hours']).days;
        const p1: mf.Point = [x, DAY_LENGTH], p2: mf.Point = [x, 0];
        coordinateTransform(p1, nOptions, true); coordinateTransform(p2, nOptions, true);
        svgString += lineSvg(p1, p2, "#00c000", 1);
    }

    // draw y-axis and gridlines
    for (let i=0; i<=24; i+=gridInterval) {
        const y = topPadding + (1-i/24) * svgHeight;
        svgString += textSvg(String(i), leftPadding-5, y, textSize, font, "#000000", "end", "middle");
        svgString += lineSvg([leftPadding, y], [leftPadding+svgWidth, y], "#808080", gridlineWidth);
    }

    // draw x-axis and gridlines
    for (let i=0; i<12; i++) {
        const x = monthEdges(isLeapYear)[i];
        const xText = (x + monthEdges(isLeapYear)[i+1])/2;
        const p1: mf.Point = [x, DAY_LENGTH], p2: mf.Point = [x, 0], pText: mf.Point = [xText, 0];
        coordinateTransform(p1, nOptions);
        coordinateTransform(p2, nOptions);
        coordinateTransform(pText, nOptions);
        svgString += textSvg(months(language)[i], pText[0], pText[1]+12, textSize, font, "#000000", "middle", "middle");
        svgString += lineSvg(p1, p2, "#808080", gridlineWidth);
    }
    svgString += lineSvg([leftPadding+svgWidth, topPadding], [leftPadding+svgWidth, topPadding+svgHeight], 
        "#808080", gridlineWidth); // right boundary of diagram
    svgString += svgClose; // complete SVG
    return svgString;
}

export function generateMoonSvg(options: moonSvgOptions) {
    const {sunEvents,moonIntervals,timeZone,newMoons,fullMoons,svgWidth=1000,svgHeight=500,leftPadding=25,rightPadding=10,
    topPadding=10,bottomPadding=25,textSize=11,font="Arial",language="en",gridInterval=2,gridlineWidth=0.5} = options;
    const days = sunEvents.length; // 365 days for common years, 366 for leap years
    const isLeapYear = (days == 366);
    const nOptions = // normalized options
    {...options,svgWidth,svgHeight,leftPadding,rightPadding,topPadding,bottomPadding,textSize,font,language,gridInterval,gridlineWidth};
    
    function toPolygons(intervals: number[][][], color: string): string[] {
        const polygons = mf.intervalsToPolygon(intervals);
        for (const polygon of polygons) {
            for (const point of polygon) {coordinateTransform(point, nOptions);}
        }
        const strings: string[] = [];
        for (const polygon of polygons) {strings.push(polygonFromArray(polygon, color));}
        return strings;
    }

    // generate SVG opening and background
    const imageWidth = svgWidth + leftPadding + rightPadding;
    const imageHeight = svgHeight + topPadding + bottomPadding;
    let svgString = svgOpen(imageWidth, imageHeight);
    svgString += rectangleSvg(0, 0, imageWidth, imageHeight, "#ffffff"); // white background

    const nIntervals: number[][][] = [];
    const cIntervals: number[][][] = [];
    for (let i=0; i<sunEvents.length; i++) {
        const [nIntervalsI, cIntervalsI] = intervalsNightCivilTwilight(sunEvents[i], timeZone);
        nIntervals.push(nIntervalsI);
        cIntervals.push(cIntervalsI);
    }
    const nPolygons = toPolygons(nIntervals, "rgba(0, 0, 0, 0.5)");
    const cPolygons = toPolygons(cIntervals, "rgba(0, 0, 0, 0.25)");
    
    const moonPolygons = toPolygons(moonIntervals, "#80c0ff");
    const allPolygons = [...moonPolygons, ...nPolygons, ...cPolygons];
    for (const polygon of allPolygons) {svgString += polygon;}

    // draw lines for new and full moons
    for (const date of newMoons) {
        const newYear = DateTime.fromISO(`${date.year}-01-01`, {zone: date.zone});
        const x = date.diff(newYear, ['days', 'hours']).days;
        const lines: mf.Point[][] = [];
        for (const interval of moonIntervals[x]) {
            lines.push([[x, interval[0]], [x, interval[1]]]);
        }
        for (const line of lines) {
            coordinateTransform(line[0], nOptions, true);
            coordinateTransform(line[1], nOptions, true);
            svgString += lineSvg(line[0], line[1], "#ff0000", 1);
        }
    }
    for (const date of fullMoons) {
        const newYear = DateTime.fromISO(`${date.year}-01-01`, {zone: date.zone});
        const x = date.diff(newYear, ['days', 'hours']).days;
        const lines: mf.Point[][] = [];
        for (const interval of moonIntervals[x]) {
            lines.push([[x, interval[0]], [x, interval[1]]]);
        }
        for (const line of lines) {
            coordinateTransform(line[0], nOptions, true);
            coordinateTransform(line[1], nOptions, true);
            svgString += lineSvg(line[0], line[1], "#000080", 1);
        }
    }

    // draw y-axis and gridlines
    for (let i=0; i<=24; i+=gridInterval) {
        const y = topPadding + (1-i/24) * svgHeight;
        svgString += textSvg(String(i), leftPadding-5, y, textSize, font, "#000000", "end", "middle");
        svgString += lineSvg([leftPadding, y], [leftPadding+svgWidth, y], "#000000", gridlineWidth);
    }

    // draw x-axis and gridlines
    for (let i=0; i<12; i++) {
        const x = monthEdges(isLeapYear)[i];
        const xText = (x + monthEdges(isLeapYear)[i+1])/2;
        const p1: mf.Point = [x, DAY_LENGTH], p2: mf.Point = [x, 0], pText: mf.Point = [xText, 0];
        coordinateTransform(p1, nOptions);
        coordinateTransform(p2, nOptions);
        coordinateTransform(pText, nOptions);
        svgString += textSvg(months(language)[i], pText[0], pText[1]+12, textSize, font, "#000000", "middle", "middle");
        svgString += lineSvg(p1, p2, "#000000", gridlineWidth);
    }
    svgString += lineSvg([leftPadding+svgWidth, topPadding], [leftPadding+svgWidth, topPadding+svgHeight], 
        "#000000", gridlineWidth); // right boundary of diagram
    svgString += svgClose; // complete SVG
    return svgString;
}