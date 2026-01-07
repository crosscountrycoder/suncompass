// This file contains constant values such as the conversion factor between degrees and radians, as well as periodic terms for 
// calculation of the sun and moon positions.

/** The number of radians in a degree, or the factor to multiply by when converting an angle from degrees to radians. 
 * To convert from radians to degrees, divide by degToRad.
*/
export const degToRad = Math.PI/180;

/** Solar elevation thresholds for sunrise/sunset (HORIZON), civil, nautical and astronomical twilight. */
export const HORIZON = (-5/6) * degToRad;
export const CIVIL_TWILIGHT = -6 * degToRad;
export const NAUTICAL_TWILIGHT = -12 * degToRad;
export const ASTRO_TWILIGHT = -18 * degToRad;

/** Earth's equatorial radius in kilometers (WGS84 ellipsoid) */
export const earthERadius = 6378.137;

/** Earth's polar radius in kilometers (WGS84 ellipsoid) */
export const earthPRadius = 6356.7523142;

/** Flattening of the WGS84 ellipsoid */
export const flattening = (earthERadius - earthPRadius) / earthERadius;

/** The Unix timestamp corresponding to 2000-01-01 at 12:00:00 UTC. */
export const J2000UTC = 946728000000;

/** Number of milliseconds in 24 hours. */
export const DAY_LENGTH = 86400000;

/** Gap to use for binary search - after this use linear or quadratic interpolation. */
export const BSEARCH_GAP = 300000; // 5 minutes

/** These periodic terms are adopted from the book "Planetary Programs and Tables from -4000 to +2800".
 * Column 0: Coefficient of sine argument (for longitude, in radians)
 * Column 1: Coefficient of cosine argument (for distance, in astronomical units)
 * Column 2: Constant added inside the sine and cosine functions
 * Column 3: Coefficient inside the sine and cosine functions
*/
export const sunPeriodicTerms = [
    [0.0403406, 0, 4.721964, 0.01621043],
    [0.0195207, -0.0097597, 5.937458, 628.30348067],
    [0.0119433, -0.0059715, 1.115589, 628.30821524],
    [0.0112392, -0.0056188, 5.781616, 628.29634302],
    [3.891e-4, -1.556e-4, 5.5474, 1256.605691],
    [2.819e-4, -1.126e-4, 1.512, 1256.609845],
    [1.721e-4, -8.61e-5, 4.1897, 628.324766],
    [0, 9.41e-5, 1.163, 0.00813],
    [6.6e-5, -2.64e-5, 5.415, 1256.5931],
    [3.5e-5, -1.63e-5, 4.315, 575.3385],
    [3.34e-5, 0, 4.553, -0.33931],
    [3.14e-5, 3.09e-5, 5.198, 7771.37715],
    [2.68e-5, -1.58e-5, 5.989, 786.04191],
    [2.42e-5, 0, 2.911, 0.05412],
    [2.34e-5, -5.4e-6, 1.423, 393.02098],
    [1.58e-5, 0, 0.061, -0.34861],
    [1.32e-5, -9.3e-6, 2.317, 1150.67698],
    [1.29e-5, -2e-6, 3.193, 157.74337],
    [1.14e-5, 0, 2.828, 52.9667],
    [9.9e-6, -4.7e-6, 0.52, 588.4927],
    [9.3e-6, 0, 4.65, 52.9611],
    [8.6e-6, 0, 4.35, -39.807],
    [7.8e-6, -3.3e-6, 2.75, 522.3769],
    [7.2e-6, -3.2e-6, 4.5, 550.7647],
    [6.8e-6, 0, 3.23, 2.6108],
    [6.4e-6, -1e-6, 1.22, 157.7385],
    [4.6e-6, -1.6e-6, 0.14, 1884.9103],
    [3.8e-6, 0, 3.44, -77.5655],
    [3.7e-6, 0, 4.37, 2.6489],
    [3.2e-6, -2.4e-6, 1.14, 1179.0627],
    [2.9e-6, -1.3e-6, 2.84, 550.7575],
    [2.8e-6, 0, 5.96, -79.6139],
    [2.7e-6, -9e-7, 5.09, 1884.8981],
    [2.7e-6, 0, 1.72, 21.3219],
    [2.5e-6, -1.7e-6, 2.56, 1097.7103],
    [2.4e-6, -1.1e-6, 1.92, 548.6856],
    [2.1e-6, 0, 0.09, 254.4393],
    [2.1e-6, 3.1e-6, 5.98, -557.3143],
    [2e-6, -1e-6, 4.03, 606.9774],
    [1.8e-6, 0, 4.27, 21.3279],
    [1.7e-6, -1.2e-6, 0.79, 1097.7163],
    [1.4e-6, 0, 4.24, -77.5282],
    [1.3e-6, -5e-7, 2.01, 1884.9191],
    [1.3e-6, 0, 2.65, 2.0781],
    [1.3e-6, 0, 4.98, 294.2463],
    [1.2e-6, 0, 0.93, -0.0799],
    [1e-6, 0, 2.21, 469.4114],
    [1e-6, 0, 3.59, -0.6829],
    [1e-6, 0, 1.5, 214.6325],
    [1e-6, -9e-7, 2.55, 1572.084],
];

/** Periodic terms for moon's ecliptic longitude and distance. Each row contains 6 values, in the following order:
 * D: multiple of mean elongation of the moon
 * M: multiple of sun's mean anomaly
 * M': multiple of moon's mean anomaly
 * F: multiple of moon's argument of latitude
 * Σl: coefficient of sine argument (for longitude, in radians)
 * Σr: coefficient of cosine argument (for distance, in kilometers)
 */
export const moonPtld =
[
    [0, 0, 1, 0, 0.109759812, -20905.355],
    [2, 0, -1, 0, 0.022235966, -3699.111],
    [2, 0, 0, 0, 0.011489747, -2955.968],
    [0, 0, 2, 0, 0.003728337, -569.925],
    [0, 1, 0, 0, -0.003230884, 48.888],
    [0, 0, 0, 2, -0.00199547, -3.149],
    [2, 0, -2, 0, 0.001026131, 246.158],
    [2, -1, -1, 0, 9.9599e-4, -152.138],
    [2, 0, 1, 0, 9.30644e-4, -170.733],
    [2, -1, 0, 0, 7.98628e-4, -204.586],
    [0, 1, -1, 0, -7.14241e-4, -129.62],
    [1, 0, 0, 0, -6.05978e-4, 108.743],
    [0, 1, 1, 0, -5.30283e-4, 104.755],
    [2, 0, 0, -2, 2.67507e-4, 10.321],
    [0, 0, 1, 2, -2.18655e-4, 0],
    [0, 0, 1, -2, 1.91637e-4, 79.661],
    [4, 0, -1, 0, 1.86314e-4, -34.782],
    [0, 0, 3, 0, 1.75126e-4, -23.21],
    [4, 0, -2, 0, 1.49191e-4, -21.636],
    [2, 1, -1, 0, -1.37672e-4, 24.208],
    [2, 1, 0, 0, -1.18089e-4, 30.824],
    [1, 0, -1, 0, -9.0111e-5, -8.379],
    [1, 1, 0, 0, 8.704e-5, -16.675],
    [2, -1, 1, 0, 7.0441e-5, -12.831],
    [2, 0, 2, 0, 6.9708e-5, -10.445],
    [4, 0, 0, 0, 6.7387e-5, -11.65],
    [2, 0, -3, 0, 6.3966e-5, 14.403],
    [0, 1, -2, 0, -4.6932e-5, -7.003],
    [2, 0, -1, 2, -4.5413e-5, 0],
    [2, -1, -2, 0, 4.1713e-5, 10.056],
    [1, 0, 1, 0, -4.098e-5, 6.322],
    [2, -2, 0, 0, 3.9026e-5, -9.884],
    [0, 1, 2, 0, -3.7001e-5, 5.751],
    [0, 2, 0, 0, -3.6111e-5, 0],
    [2, -2, -1, 0, 3.5744e-5, -4.95],
    [2, 0, 1, -2, -3.0945e-5, 4.13],
    [2, 0, 0, 2, -2.7838e-5, 0],
    [4, -1, -1, 0, 2.1206e-5, -3.958],
    [0, 0, 2, 2, -1.9373e-5, 0],
    [3, 0, -1, 0, -1.5568e-5, 3.258],
    [2, 1, 1, 0, -1.4137e-5, 2.616],
    [4, -1, -2, 0, 1.3247e-5, -1.897],
    [0, 2, -1, 0, -1.2444e-5, -2.117],
    [2, 2, -1, 0, -1.2217e-5, 2.354],
    [2, 1, -2, 0, 1.206e-5, 0],
    [2, -1, 0, -2, 1.0402e-5, 0],
    [4, 0, 1, 0, 9.582e-6, -1.423],
    [0, 0, 4, 0, 9.372e-6, -1.117],
    [4, -1, 0, 0, 9.076e-6, -1.571],
    [1, 0, -2, 0, -8.5e-6, -1.739],
    [2, 1, 0, -2, -6.964e-6, 0],
    [0, 0, 2, -2, -6.65e-6, -4.421],
    [1, 1, 1, 0, 6.126e-6, 0],
    [3, 0, -2, 0, -5.934e-6, 0],
    [4, 0, -3, 0, 5.76e-6, 0],
    [2, -1, 2, 0, 5.707e-6, 0],
    [0, 2, 1, 0, -5.637e-6, 1.165],
    [1, 1, -1, 0, 5.219e-6, 0],
    [2, 0, 3, 0, 5.131e-6, 0],
    [2, 0, -1, -2, 0, 8.752],
];

/** Periodic terms for moon's ecliptic latitude. Each row contains 6 values, in the following order:
 * D: multiple of mean elongation of the moon
 * M: multiple of sun's mean anomaly
 * M': multiple of moon's mean anomaly
 * F: multiple of moon's argument of latitude
 * Σl: coefficient of sine argument, in radians
 */
export const moonPtl = 
[
    [0, 0, 0, 1, 0.089502613],
    [0, 0, 1, 1, 0.004897429],
    [0, 0, 1, -1, 0.004846657],
    [2, 0, 0, -1, 0.003023556],
    [2, 0, -1, 1, 9.67139e-4],
    [2, 0, -1, -1, 8.07581e-4],
    [2, 0, 0, 1, 5.68506e-4],
    [0, 0, 2, 1, 3.00162e-4],
    [2, 0, 1, -1, 1.61722e-4],
    [0, 0, 2, -1, 1.53973e-4],
    [2, -1, 0, -1, 1.43396e-4],
    [2, 0, -2, -1, 7.5468e-5],
    [2, 0, 1, 1, 7.3304e-5],
    [2, 1, 0, -1, -5.8626e-5],
    [2, -1, -1, 1, 4.2987e-5],
    [2, -1, 0, 1, 3.8589e-5],
    [2, -1, -1, -1, 3.6041e-5],
    [0, 1, -1, -1, -3.2638e-5],
    [4, 0, -1, -1, 3.1905e-5],
    [0, 1, 0, 1, -3.1311e-5],
    [0, 0, 0, 3, -3.0526e-5],
    [0, 1, -1, 1, -2.7314e-5],
    [1, 0, 0, 1, -2.6023e-5],
    [0, 1, 1, 1, -2.5744e-5],
    [0, 1, 1, -1, -2.4609e-5],
    [0, 1, 0, -1, -2.3457e-5],
    [1, 0, 0, -1, -2.33e-5],
    [0, 0, 3, 1, 1.9321e-5],
    [4, 0, 0, -1, 1.782e-5],
    [4, 0, -1, 1, 1.4539e-5],
    [0, 0, 1, -3, 1.3561e-5],
    [4, 0, -2, 1, 1.1711e-5],
    [2, 0, 0, -3, 1.0594e-5],
    [2, 0, 2, -1, 1.0402e-5],
    [2, -1, 1, -1, 8.57e-6],
    [2, 0, -2, 1, -7.871e-6],
    [0, 0, 3, -1, 7.662e-6],
    [2, 0, 2, 1, 7.365e-6],
    [2, 0, -3, -1, 7.348e-6],
    [2, 1, -1, 1, -6.388e-6],
    [2, 1, 0, 1, -6.126e-6],
    [4, 0, 0, 1, 5.777e-6],
    [2, -1, 1, 1, 5.498e-6],
    [2, -2, 0, -1, 5.271e-6],
    [0, 0, 1, 3, -4.939e-6],
    [2, 1, 1, -1, -3.997e-6],
    [1, 1, 0, -1, 3.892e-6],
    [1, 1, 0, 1, 3.892e-6],
    [0, 1, -2, -1, -3.84e-6],
    [2, 1, -1, -1, -3.84e-6],
    [1, 0, 1, 1, -3.229e-6],
    [2, -1, -2, -1, 3.159e-6],
    [0, 1, 2, 1, -3.089e-6],
    [4, 0, -2, -1, 3.072e-6],
    [4, -1, -1, -1, 2.897e-6],
    [1, 0, 1, -1, -2.862e-6],
    [4, 0, 1, -1, 2.304e-6],
    [1, 0, -1, -1, -2.077e-6],
    [4, -1, 0, -1, 2.007e-6],
    [2, -2, 0, 1, 1.868e-6],
];

/** Approximate values of delta T for selected years from 0 through 2500. Delta T is linearly interpolated between years.
 * 
 * Sources:
 *  
 * https://eclipse.gsfc.nasa.gov/SEhelp/deltat2004.html#tab1 (before 1700)
 * 
 * https://maia.usno.navy.mil/ser7/deltat.data (1700 - 2025)
 * 
 * https://maia.usno.navy.mil/ser7/finals.all (projection for 2027)
 * 
 * https://eclipse.gsfc.nasa.gov/LEcat5/deltatpoly.html (projections for 2050 - 2500)
 * 
 */
export const deltaT = [
    [0, 10580],
    [100, 9600],
    [200, 8640],
    [300, 7680],
    [400, 6700],
    [500, 5710],
    [600, 4740],
    [700, 3810],
    [800, 2960],
    [900, 2200],
    [1000, 1570],
    [1100, 1090],
    [1200, 740],
    [1300, 490],
    [1400, 320],
    [1500, 200],
    [1600, 120],
    [1700, 21],
    [1720, 21.1],
    [1740, 13.5],
    [1760, 14.8],
    [1780, 15.6],
    [1800, 12.6],
    [1820, 11.1],
    [1840, 6.2],
    [1860, 7.4],
    [1880, -5.4],
    [1900, -2.7],
    [1910, 10.4],
    [1920, 21.4],
    [1930, 24],
    [1940, 24.4],
    [1950, 29.2],
    [1960, 33.2],
    [1970, 40.2],
    [1980, 50.5],
    [1985, 54.3],
    [1990, 56.9],
    [1995, 60.8],
    [2000, 63.8],
    [2005, 64.7],
    [2010, 66.1],
    [2015, 67.6],
    [2020, 69.4],
    [2025, 69.1],
    [2027, 69.1],
    [2050, 93],
    [2100, 203],
    [2150, 328],
    [2200, 442],
    [2250, 572],
    [2300, 717],
    [2400, 1056],
    [2500, 1460]
];