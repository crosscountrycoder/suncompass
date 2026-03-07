/**
 * This TypeScript file contains code to generate SVG chart of day length and sunrise/sunset/twilight times for an entire year.
 * "generateSvg" generates an SVG file showing (from top to bottom) the lengths of day, civil twilight, nautical twilight, astronomical
 * twilight and night for an entire year
 */

import {getTimeOfDay} from "./lookup-tables.ts";
import * as mf from "./mathfuncs.ts";
import {intervalsSvg, lengths, intervalsNightCivilTwilight, type SunTable} from "./suncalc.ts";
import type {MoonTable} from "./mooncalc.ts";

const sunColors = ["#80c0ff", "#0060c0", "#004080", "#002040", "#000000"];

type sunSvgOptions = {
    sunTable: SunTable;
    type: string;
    title: string;
    svgWidth?: number;
    svgHeight?: number;
    diagramWidth?: number;
    diagramHeight?: number;
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
    sunTable: SunTable;
    moonTable: MoonTable;
    title: string;
    svgWidth?: number;
    svgHeight?: number;
    diagramWidth?: number;
    diagramHeight?: number;
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

export function indent(n: number): string {return "    ".repeat(n);}

/** Generates the opening of an SVG */
export function svgOpen(width: number, height: number, viewBoxWidth: number, viewBoxHeight: number, title: string): string {
    return `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${viewBoxWidth} ${viewBoxHeight}">\n`
    + indent(1) + `<title>${title}</title>\n`
    + indent(1) + `<style>.nss line, .nss path, .nss rect {vector-effect: non-scaling-stroke;}</style>\n`
    + indent(1) + `<g class="nss">\n`;
}

export function simplifyCollinear(points: mf.Polygon): mf.Polygon {
    if (points.length <= 2) return [...points];
    const simplified: mf.Polygon = [points[0]];
    for (let i = 1; i < points.length - 1; i++) {
        const [x0, y0] = simplified[simplified.length - 1];
        const [x1, y1] = points[i];
        const [x2, y2] = points[i + 1];
        const sameX = x0 === x1 && x1 === x2;
        const sameY = y0 === y1 && y1 === y2;
        if (!sameX && !sameY) {simplified.push(points[i]);}
    }
    simplified.push(points[points.length - 1]);
    return simplified;
}

export function pathFromArray(
    points: mf.Polygon[],
    closed: boolean,
    fillColor?: string,
    strokeColor?: string,
    strokeWidth?: number,
    id?: string,
    orthogonal: boolean = false,
    precision: number = 2,
): string {
    const subpaths: string[] = [];
    for (const poly0 of points) {
        const poly = simplifyCollinear(poly0);
        if (poly.length === 0) continue;
        const [firstX, firstY] = poly[0];
        const moveTo = `M${mf.toFixedS(firstX, precision)},${mf.toFixedS(firstY, precision)}`;
        let lineStr = "";
        if (orthogonal) {
            const cmds: string[] = [];
            for (let i = 1; i < poly.length; i++) {
                const [x0, y0] = poly[i - 1];
                const [x1, y1] = poly[i];
                if (x1 === x0) {cmds.push(`V${mf.toFixedS(y1, precision)}`);}
                else if (y1 === y0) {cmds.push(`H${mf.toFixedS(x1, precision)}`);}
                else {cmds.push(`L${mf.toFixedS(x1, precision)},${mf.toFixedS(y1, precision)}`);}
            }
            lineStr = cmds.join("");
        }
        else {
            lineStr = poly.slice(1).map(([x, y]) => `L${mf.toFixedS(x, precision)},${mf.toFixedS(y, precision)}`).join("");
        }
        const closeCmd = closed ? "Z" : "";
        subpaths.push(lineStr ? `${moveTo}${lineStr}${closeCmd}` : `${moveTo}${closeCmd}`);
    }

    if (subpaths.length === 0) return "";
    const d = subpaths.join("");
    const attrs = [
        id === undefined ? `` : `id="${id}"`,
        fillColor === undefined ? `` : `fill="${fillColor}"`,
        strokeColor === undefined ? `` : `stroke="${strokeColor}"`,
        strokeWidth === undefined ? `` : `stroke-width="${strokeWidth}"`,
        `d="${d}"`
    ].filter(Boolean).join(" ");

    return `<path ${attrs}/>\n`;
}

/** Generates SVG code for a rectangle with the top-left corner at the given x and y cordinates, and the given width, height,
 * fill and stroke colors. */
export function rectangleSvg(
    x: number, 
    y: number, 
    width: number,
    height: number, 
    fillColor?: string, 
    strokeColor?: string,
    strokeWidth?: number,
    id?: string
) {
    const idString = (id === undefined) ? `` : `id="${id}"`;
    return `<rect ${idString} x="${x}" y="${y}" width="${width}" height="${height}"`
    + (fillColor === undefined ? `` : ` fill="${fillColor}"`) 
    + (strokeColor === undefined ? `` : ` stroke="${strokeColor}"`)
    + (strokeWidth === undefined ? `` : ` stroke-width="${strokeWidth}"`) + `/>\n`;
}

/** Generates SVG code for a text box with the given text, centered at the given x and y coordinate, with the given font and font size. 
 * Text anchor can be "start" (left-aligned), "middle" (centered), or "end" (right-aligned). Alignment baseline can be either 
 * "text-before-edge" (top-aligned), "middle" (centered), or "text-after-edge" (bottom-aligned). */
export function textSvg(text: string, x: number, y: number, precision: number = 2): string {
    return `<text x="${mf.toFixedS(x,precision)}" y="${mf.toFixedS(y,precision)}">${text}</text>\n`;
}

/** Generates an SVG line from (x1, y1) to (x2, y2) with the given color and width. */
export function lineSvg(p1: mf.Point, p2: mf.Point, color?: string, width?: number, precision: number = 2): string {
    return `<line x1="${mf.toFixedS(p1[0],precision)}" y1="${mf.toFixedS(p1[1],precision)}" x2="${mf.toFixedS(p2[0],precision)}"`
    + ` y2="${mf.toFixedS(p2[1],precision)}"` + (color === undefined ? `` : ` stroke="${color}"`)
    + (width === undefined ? `` : ` stroke-width="${width}"`) + `/>\n`;
}

/** Returns an array of month abbreviations in the given language, represented by a language code, such as "en" (English), "es"
 * (Spanish), "fr" (French), "zh" (Chinese). So far there is only English - I plan to add more when I localize the site. */
export function months(language: string = "en"): string[] {
    return ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
}

/** Edges of the months, used for drawing gridlines. */
export function monthEdges(leapYear: boolean = false): number[] {
    if (leapYear) {return [0,31,60,91,121,152,182,213,244,274,305,335,366];}
    else {return [0,31,59,90,120,151,181,212,243,273,304,334,365];}
}

/** This returns an array of two strings: [lineString, textString].
 * lineString is intended to be enclosed within the <g id="graph"> block. textString should be outside the block.
 */
export function generateGrid(
    options: sunSvgOptions | moonSvgOptions, 
    gridlineColor: string, 
    lineIndent: number,
    textIndent: number
    ): [string, string] {
    // draw y-axis and gridlines
    const gridInterval = options.gridInterval!;
    const topPadding = options.topPadding!;
    const diagramWidth = options.diagramWidth!;
    const diagramHeight = options.diagramHeight!;
    const leftPadding = options.leftPadding!;
    const textSize = options.textSize!;
    const font = options.font!;
    const gridlineWidth = options.gridlineWidth!;
    const days = options.sunTable.solarEvents.length;
    const isLeapYear = (days === 366);
    const scaleX = diagramWidth / days;
    const language = options.language!;

    let textString =
    `${indent(textIndent)}<g id="gridtext" font-family="${font}" font-size="${textSize}pt" dominant-baseline="middle" fill="#000000">\n` + 
    `${indent(textIndent+1)}<g text-anchor="end">\n`; // string containing x and y axis labels

    for (let i=0; i<=24; i+=gridInterval) {
        const y = topPadding + (1-i/24) * diagramHeight;
        textString += indent(textIndent+2) + textSvg(String(i), leftPadding-5, y);
    }
    textString += `${indent(textIndent+1)}</g>\n${indent(textIndent+1)}<g text-anchor="middle">\n`;
    
    for (let i=0; i<12; i++) {
        const xText = (monthEdges(isLeapYear)[i] + monthEdges(isLeapYear)[i+1])/2 * scaleX + leftPadding;
        const y = topPadding+diagramHeight;
        textString += indent(textIndent+2) + textSvg(months(language)[i], xText, y+12);
    }
    textString += `${indent(textIndent+1)}</g>\n${indent(textIndent)}</g>\n`;

    let verticalLines = ``;
    for (let i=0; i<13; i++) {verticalLines += `M${monthEdges(isLeapYear)[i]},0V86400`;}
    let horizontalLines = ``;
    for (let i=0; i<=24; i+=gridInterval) {horizontalLines += `M0,${i*3600}H${days}`;}

    const lineString = indent(lineIndent) +
    `<path id="gridlines" stroke="${gridlineColor}" stroke-width="${gridlineWidth}" fill="none" d="${horizontalLines}${verticalLines}"/>\n`;
    return [lineString, textString];
}

/**
 * Returns a string containing an SVG diagram for either day/twilight/night lengths throughout the year, or the times of day in which
 * day, night, and each stage of twilight occur.
 * Parameters should be passed in an object. All parameters except sunTable, type, timeZone, and solsticesEquinoxes are optional.
 * @param sunTable Values of "generateSunTable" for each day of the year.
 * @param type Set to "length" to generate a day/night/twilight length chart, or "rise-set" to generate a chart with times of day.
 * @param title Title of the SVG.
 * @param svgWidth Width of the SVG file.
 * @param svgHeight Height of the SVG file.
 * @param diagramWidth Width of the SVG's view box.
 * @param diagramHeight Height of the SVG's view box.
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
 * The total width of the SVG file is equal to diagramWidth + leftPadding + rightPadding. The height is equal to diagramHeight + topPadding +
 * bottomPadding.
 */
export function generateSunSvg(options: sunSvgOptions): string {
    const {sunTable,type,title,svgWidth=1035,svgHeight=535,diagramWidth=1000,diagramHeight=500,leftPadding=25,rightPadding=10,
        topPadding=10,bottomPadding=25,textSize=11,font="Arial,sans-serif",language="en",gridInterval=2,gridlineWidth=0.5} = options;
    const days = sunTable.solarEvents.length; // 365 days for common years, 366 for leap years
    const nOptions = // normalized options
    {...options,diagramWidth,diagramHeight,leftPadding,rightPadding,topPadding,bottomPadding,textSize,font,language,gridInterval,gridlineWidth};

    /** Converts a set of durations into an array of points representing a polygon. This includes coordinate transforms. */
    function durationsToArray(durations: number[]): mf.Polygon[] {
        const p: mf.Polygon[] = [[]]; // p is short for polygons
        for (let i=0; i<days; i++) {
            if (durations[i] > 0) {
                if (i === 0 || durations[i-1] === 0) {p[0].push([i, 0], [i, durations[i]], [i+1, durations[i]]);}
                else {p[p.length-1].push([i, durations[i]], [i+1, durations[i]]);}
                if (i === days-1) {p[p.length-1].push([days, 0]);}
            }
            else if (i !== 0 && durations[i-1] > 0) {p[p.length-1].push([i, 0]);}
        }
        return p;
    }

    /** Used to draw lines representing solar noon on the graph. */
    function solarNoonLines() {
        const solarNoons: number[][] = [];
        for (const evts of sunTable.solarEvents) {
            const curDay: number[] = [];
            for (const event of evts) {
                if (event.type === "Solar Noon") {curDay.push(Math.floor(getTimeOfDay(event.unix, sunTable.timeZoneTable)/1000));}
            }
            solarNoons.push(curDay);
        }
        
        const groups: mf.Polygon[] = []; // a group of multiple lines, each representing a cluster of solar noons
        for (const solarNoon of solarNoons[0]) {groups.push([[0, solarNoon], [1, solarNoon]]);}
        for (let i=1; i<days; i++) { // for each day of the year
            for (const noon of solarNoons[i]) { // for each solar noon of the day (may be more than 1)
                let flag: boolean = false;
                for (const group of groups) {
                    if (Math.abs(noon - group[group.length-1][1]) < 1800 && group[group.length-1][0] === i) {
                        flag = true;
                        group.push([i, noon], [i+1, noon]);
                        break;
                    }
                }
                if (!flag) {groups.push([[i, noon], [i+1, noon]]);}
            }
        }
        return indent(3) + pathFromArray(groups, false, undefined, "#ff0000", undefined, "solar-noon", true);
    }

    /** Used to draw lines representing solar midnight on the graph. */
    function solarMidnightLines() {
        const solarMidnights: number[][] = [];
        for (const evts of sunTable.solarEvents) {
            const curDay: number[] = [];
            for (const event of evts) {
                if (event.type === "Solar Midnight") {curDay.push(Math.floor(getTimeOfDay(event.unix, sunTable.timeZoneTable)/1000));}
            }
            solarMidnights.push(curDay);
        }
        
        const groups: mf.Polygon[] = []; // a group of multiple lines (number[][]), each representing a cluster of solar midnights
        for (const solarMidnight of solarMidnights[0]) {groups.push([[0, solarMidnight], [1, solarMidnight]]);}
        for (let i=1; i<days; i++) { // for each day of the year
            for (const midnight of solarMidnights[i]) { // for each solar midnight of the day (may be more than 1)
                let flag: boolean = false;
                for (const group of groups) {
                    if (Math.abs(midnight - group[group.length-1][1]) < 1800 && group[group.length-1][0] === i) {
                        flag = true;
                        group.push([i, midnight], [i+1, midnight]);
                        break;
                    }
                }
                if (!flag) {groups.push([[i, midnight], [i+1, midnight]]);}
            }
        }
        return indent(3) + pathFromArray(groups, false, undefined, "#0000ff", undefined, "solar-midnight", true);
    }

    function toPolygons(intervals: number[][][], color: string, id?: string): string {
        const polygons = mf.intervalsToPolygon(intervals);
        return pathFromArray(polygons, true, color, undefined, undefined, id, true);
    }

    // generate SVG opening and background
    const vbWidth = diagramWidth + leftPadding + rightPadding;
    const vbHeight = diagramHeight + topPadding + bottomPadding;
    let svgString = svgOpen(svgWidth, svgHeight, vbWidth, vbHeight, title);
    svgString += indent(2) + rectangleSvg(0, 0, vbWidth, vbHeight, "#ffffff", undefined, undefined, "background");

    // special coordinate system for graph (x is in days, y is in seconds)
    const scaleX = diagramWidth / days;
    const scaleY = -diagramHeight / 86400;
    svgString += indent(2) + `<g id="graph" transform="translate(${leftPadding},${topPadding+diagramHeight})` + 
    ` scale(${scaleX.toPrecision(8)},${scaleY.toPrecision(8)})" stroke-width="1" fill="none">\n`;

    if (type === "length") { // day/twilight/night length plot
        const dLengths: number[] = []; // day lengths
        const cLengths: number[] = []; // day + civil twilight lengths
        const nLengths: number[] = []; // day + civil + nautical twilight lengths
        const aLengths: number[] = []; // day + civil + nautical + astronomical twilight lengths
        for (const e of sunTable.solarEvents) {
            const dur = lengths(e, sunTable.timeZoneTable);
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
        svgString += indent(3) + rectangleSvg(0, 0, days, 86400, "#000000", undefined, undefined, "night"); // night
        svgString += indent(3) + pathFromArray(ap, true, sunColors[3], undefined, undefined, "astro-twilight", true);
        svgString += indent(3) + pathFromArray(np, true, sunColors[2], undefined, undefined, "nautical-twilight", true);
        svgString += indent(3) + pathFromArray(cp, true, sunColors[1], undefined, undefined, "civil-twilight", true);
        svgString += indent(3) + pathFromArray(dp, true, sunColors[0], undefined, undefined, "daylight", true);
    }
    else if (type === "rise-set") { // sunrise, sunset, dusk, dawn plot
        const aIntervals: number[][][] = []; // intervals of astronomical twilight or brighter
        const nIntervals: number[][][] = []; // intervals of nautical twilight or brighter
        const cIntervals: number[][][] = []; // intervals of civil twilight or brighter
        const dIntervals: number[][][] = []; // intervals of daylight

        for (const event of sunTable.solarEvents) {
            const int = intervalsSvg(event, sunTable.timeZoneTable);
            aIntervals.push(int[3]);
            nIntervals.push(int[2]);
            cIntervals.push(int[1]);
            dIntervals.push(int[0]);
        }
        
        svgString += indent(3) + rectangleSvg(0, 0, days, 86400, "#000000", undefined, undefined, "night");
        const aPolygons = toPolygons(aIntervals, sunColors[3], "astro-twilight");
        const nPolygons = toPolygons(nIntervals, sunColors[2], "nautical-twilight");
        const cPolygons = toPolygons(cIntervals, sunColors[1], "civil-twilight");
        const dPolygons = toPolygons(dIntervals, sunColors[0], "daylight");
        const allPolygons = [aPolygons, nPolygons, cPolygons, dPolygons];

        for (const polygon of allPolygons) {svgString += indent(3) + polygon;}
        svgString += solarMidnightLines() + solarNoonLines();
    }
    
    // draw solstices and equinoxes as green lines
    const [me, js, se, ds] = Object.values(sunTable.solsticeDates);
    svgString += indent(3) + `<path id="solstices-equinoxes" stroke="#00c000" d="M${me+0.5},0V86400M${js+0.5},0V86400M${se+0.5},0V86400M${ds+0.5},0V86400"/>\n`;

    // add gridlines
    const [gridLines, gridText] = generateGrid(nOptions, "#808080", 3, 2);
    svgString += (gridLines + indent(2) + "</g>\n" + gridText + indent(1) + "</g>\n</svg>\n");
    return svgString;
}

/**
 * Returns a string containing an SVG diagram showing moonrise and moonset times for every day of the year, along with new and
 * full moons and overlays for night and civil twilight.
 * Parameters should be passed in an object. All parameters except events, type, timeZone, and solsticesEquinoxes are optional.
 * @param sunTable Values of "generateSunTable" for each day of the year.
 * @param moonTable Values of "generateMoonTable" for each day of the year.
 * @param title Title of the SVG.
 * @param svgWidth Width of the SVG file.
 * @param svgHeight Height of the SVG file.
 * @param diagramWidth Width of the SVG's view box.
 * @param diagramHeight Height of the SVG's view box.
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
 * The total width of the SVG file is equal to diagramWidth + leftPadding + rightPadding. The height is equal to diagramHeight + topPadding +
 * bottomPadding.
 */
export function generateMoonSvg(options: moonSvgOptions) {
    const {sunTable,moonTable,title,svgWidth=1035,svgHeight=535,diagramWidth=1000,diagramHeight=500,leftPadding=25,rightPadding=10,
        topPadding=10,bottomPadding=25,textSize=11,font="Arial,sans-serif",language="en",gridInterval=2,gridlineWidth=0.5} = options;
    const days = sunTable.solarEvents.length; // 365 days for common years, 366 for leap years
    const nOptions = // normalized options
    {...options,diagramWidth,diagramHeight,leftPadding,rightPadding,topPadding,bottomPadding,textSize,font,language,gridInterval,gridlineWidth};
    
    function toPolygons(intervals: number[][][], color?: string, id?: string): string {
        const polygons = mf.intervalsToPolygon(intervals);
        return pathFromArray(polygons, true, color, undefined, undefined, id, true);
    }

    // generate SVG opening and background
    const viewBoxWidth = diagramWidth + leftPadding + rightPadding;
    const viewBoxHeight = diagramHeight + topPadding + bottomPadding;
    let svgString = svgOpen(svgWidth, svgHeight, viewBoxWidth, viewBoxHeight, title);
    svgString += indent(2) + rectangleSvg(0, 0, viewBoxWidth, viewBoxHeight, "#ffffff", undefined, undefined, "background");

    // special coordinate system for graph (x is in days, y is in seconds)
    const scaleX = diagramWidth / days;
    const scaleY = -diagramHeight / 86400;
    svgString += indent(2) + `<g id="graph" transform="translate(${leftPadding},${topPadding+diagramHeight})` + 
    ` scale(${scaleX.toPrecision(8)},${scaleY.toPrecision(8)})" stroke-width="1">\n`;
    
    // add light blue polygons (when moon above horizon)
    svgString += indent(3) + toPolygons(moonTable.intervals, "#80c0ff", "moonlight");

    // draw lines for new and full moons
    let newMoonString = indent(3) + `<path id="new-moons" stroke="#ff0000" d="`;
    let fullMoonString = indent(3) + `<path id="full-moons" stroke="#0000ff" d="`;
    for (const phase of moonTable.phases) {
        const x = phase.date;
        for (const interval of moonTable.intervals[x]) {
            if (phase.type === 0) {newMoonString += `M${x+0.5},${interval[0]}V${interval[1]}`;}
            else if (phase.type === 2) {fullMoonString += `M${x+0.5},${interval[0]}V${interval[1]}`;}
        }
    }
    newMoonString += `"/>\n`;
    fullMoonString += `"/>\n`;
    svgString += newMoonString + fullMoonString;

    // add night and civil twilight overlays
    const nIntervals: number[][][] = [];
    const cIntervals: number[][][] = [];
    for (let i=0; i<days; i++) {
        const [nIntervalsI, cIntervalsI] = intervalsNightCivilTwilight(sunTable.solarEvents[i], sunTable.timeZoneTable);
        nIntervals.push(nIntervalsI);
        cIntervals.push(cIntervalsI);
    }
    svgString += indent(3) + toPolygons(nIntervals, "rgba(0, 0, 0, 0.5)", "night"); // night overlay
    svgString += indent(3) + toPolygons(cIntervals, "rgba(0, 0, 0, 0.25)", "civil-twilight"); // civil twilight overlay

    const [gridLines, gridText] = generateGrid(nOptions, "#000000", 3, 2);
    svgString += (gridLines + indent(2) + "</g>\n" + gridText + indent(1) + "</g>\n</svg>\n");
    return svgString;
}