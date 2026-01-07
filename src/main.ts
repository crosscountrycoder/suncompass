// src/main.ts
import { generateLODProfile } from "./core/lookup-tables";
import { subsolarPoint } from "./core/suncalc";

function updateSubsolarPoint() {
    const [lat, lon] = subsolarPoint(generateLODProfile(Date.now()), true); // subsolar lat/lon in degrees
    const ssp = document.getElementById("position");
    if (!ssp) {
        console.warn("Element with id='position' not found");
        return;
    }
    ssp.textContent = `${lat.toFixed(4)}, ${lon.toFixed(4)}`;
}

function startSubsolarUpdates() { // Align to the next integer second
    const now = new Date();
    const delay = 1000 - now.getMilliseconds();
    setTimeout(() => {
        updateSubsolarPoint();
        setInterval(updateSubsolarPoint, 1000); // Then update every 1000 ms
    }, delay);
}

window.addEventListener("DOMContentLoaded", () => {
    startSubsolarUpdates();
});