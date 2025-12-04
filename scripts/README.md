Contains the astronomical formulas used in Sun Compass. There are five executable TypeScript files (in the scripts folder), which 
rely on functions contained within other files.

Note: All files MUST be run from the project's root directory, not from the "scripts" directory!

## rise-set.ts

For any place and time on Earth, gives the position of the sun, the distance between the sun and earth, and the times of sunrise, 
sunset, solar noon, solar midnight, and civil, nautical and astronomical twilight. Usage is as follows:

`npx ts-node scripts/rise-set.ts`: Prints a help menu with the syntax and usage of this command.

`npx ts-node scripts/rise-set.ts <latitude> <longitude>`: Prints the current sun position and distance, and sunrise, sunset and 
twilight times for today. Coordinates are given in decimal degrees, for example `npx ts-node rise-set.ts 34.05 -118.25` gives 
data for Los Angeles, California.

`npx ts-node scripts/rise-set.ts <latitude> <longitude> <time>`: Prints the sun's position and distance at the given latitude and 
longitude coordinates. For example, `npx ts-node rise-set.ts 40.75 -73.99 2025-12-31T16:15:00` gives the sun position and 
distance for Manhattan, New York City on December 31, 2025 at 16:15 (4:15 pm) eastern standard time, along with sunrise, sunset, 
and twilight times for December 31 in Manhattan.

`npx ts-node scripts/rise-set.ts <latitude> <longitude> <me/js/se/ds>` gives sun position and sunrise/sunset/twilight data for 
the given latitude and longitude at the March equinox (me), June solstice (js), September equinox (se), or December solstice (ds).

Note: to find data for the north or south poles, set the latitude to ±89.9999, not ±90. Solar noon and midnight are undefined at 
latitude ±90.

## solstice.ts

`npx ts-node scripts/solstice.ts`: Prints the solstices and equinoxes in the current year (in the local time zone of the device) 
along with the subsolar point (the position on Earth at which the sun is directly overhead) at each solstice and equinox.

`npx ts-node scripts/solstice.ts <year>`: Prints the times of solstices and equinoxes in the given year, in the device's local 
time zone, along with subsolar points.

`npx ts-node scripts/solstice.ts <year> <time zone>`: Prints the times of solstices and equinoxes in the given year and time 
zone. 
The time zones are either in IANA time zone format (ex. "America/Los_Angeles" for Pacific Time) or a fixed UTC offset (ex. "utc", 
"utc-8", "utc+5:30").

## gen-solsteq-json.ts

`npx ts-node scripts/gen-solsteq-json.ts`: Generates a JSON file containing solstice and equinox times in UTC for all years from 
0 to 2500.

## generate-svg.ts

`npx ts-node scripts/generate-svg.ts <latitude> <longitude>`: Creates two SVG diagrams:

**day-lengths.svg**: Shows the lengths of day, night, and each stage of twilight for every day of the year. There are also green 
lines, which represent solstices and equinoxes.

**sunrise-sunset.svg**: Shows sunrise, sunset, and twilight times for every day of the year. Green lines represent solstices and 
equinoxes, the red line is solar noon, and the blue line is solar midnight. In some locations, there are abrupt "jumps" in the 
graph. 
This represents daylight saving time, also known as summer time.

The files "ucsb-day-lengths-2025.svg" and "ucsb-sunrise-sunset-2025.svg" show day lengths and sunrise/sunset times in 2025 for 
coordinates 34.42, -119.85 (on the campus of UC Santa Barbara, where I go to college).

To create charts for the north pole or south pole, set the coordinates to `89.9999 0` or `-89.9999 0`. Do not set the latitude to 
90 or -90.

## moonrise-set.ts
`npx ts-node scripts/moonrise-set.ts <latitude> <longitude>`: Prints the times of moonrise and moonset along with other lunar 
data.