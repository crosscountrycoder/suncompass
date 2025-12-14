Contains the astronomical formulas used in Sun Compass. There are five executable TypeScript files (in the scripts folder), which rely on functions contained within other files.

Note: All files MUST be run from the project's root directory, not from the "scripts" directory!

## sun.ts

For any place and time on Earth, gives the position of the sun, the distance between the sun and earth, and the times of sunrise, sunset, solar noon, solar midnight, and civil, nautical and astronomical twilight. Usage is as follows:

`npx ts-node scripts/sun.ts <latitude> <longitude> [time]`: Prints the sun position and distance at the given time, subsolar point (where the sun is directly overhead), and sunrise, sunset, twilight, solar noon, and solar midnight times for the day. 

The latitude and longitude are given in signed decimal degrees, so the coordinates for downtown Los Angeles would be specified as `34.05 -118.25` (not "34.05°N 118.25°W" or "34°3'N 118°15'W"). Time is specified in ISO format, which is the concatenation of the date in yyyy-mm-dd format, the letter "T", and the time in 24-hour hh:mm:ss format. For example `2025-12-21T15:02:45` refers to 3:02:45 pm on December 21, 2025 in the local time zone of the given coordinates.

Time may also be specified as `me`, `js`, `se`, or `ds`, which refer to the current year's March equinox, June solstice, September equinox, and December solstice respectively. If time is omitted, it defaults to the current time at the given location.

Note: to find data for the north or south poles, set the latitude to ±89.9999, not ±90. Solar noon and midnight are undefined at  latitude ±90.

## solstice.ts

`npx ts-node scripts/solstice.ts [year] [timezone]`: Prints the times of solstices and equinoxes in the given year and time zone. The time zones are either in IANA time zone format (ex. "America/Los_Angeles" for Pacific Time) or a fixed UTC offset (ex.  "utc", "utc-8", "utc+5:30"). If the year and time zone are omitted, they default to the current year and the device's time zone respectively.

## gen-solsteq-json.ts

`npx ts-node scripts/gen-solsteq-json.ts`: Generates a JSON file containing solstice and equinox times in UTC for all years from 0 to 2500.

## generate-svg.ts

`npx ts-node scripts/generate-svg.ts <latitude> <longitude> [year] [location]`: Creates three SVG diagrams:

**day-lengths.svg**: Shows the lengths of day, night, and each stage of twilight for every day of the year. There are also green  lines, which represent solstices and equinoxes.

**sunrise-sunset.svg**: Shows sunrise, sunset, and twilight times for every day of the year. Green lines represent solstices and  equinoxes, the red line is solar noon, and the blue line is solar midnight. In some locations, there are abrupt "jumps" in the  graph. This represents daylight saving time, also known as summer time.

**moon.svg**: Shows moonrise, moonset, and moon phases for the entire year. Times at which the moon is above the horizon are shaded in light blue. New and full moons are marked as red and blue lines, respectively. The dark shaded overlays represent night and civil twilight.

If the year is not specified, it defaults to the current year. If location is specified, the file names will state the given location name and year. For example, `npx ts-node scripts/generate-svg.ts 40.75 -73.99 2025 nyc` will create files titled `nyc-day-lengths-2025.svg`, `nyc-sunrise-sunset-2025.svg` and `nyc-moon-2025.svg`. (The coordinates 40.75, -73.99 are in midtown Manhattan.)

The SVG files `ucsb-day-lengths-2025`, `ucsb-moon-2025` and `ucsb-sunrise-sunset-2025` are for coordinates 34.42, -119.85, on the campus of UC Santa Barbara, where I go to college.

To create charts for the north pole or south pole, set the coordinates to `89.9999 0` or `-89.9999 0`. Do not set the latitude to 90 or -90.

## moon.ts
`npx ts-node scripts/moon.ts <latitude> <longitude> [time]`: Prints the position (elevation and compass bearing), sublunar point, distance, phase and illumination of the moon at the given latitude, longitude, and time. Additionally, moonrise, moonset, and meridian transit times are shown. If the phase is "new moon", "first quarter", "full moon", or "last quarter", it shows the exact times of these events.

Unlike for sun.ts, time cannot be specified as "me", "js", "se" or "ds". If time is omitted, it defaults to the current time at the given location.