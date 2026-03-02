// src/main.ts
import "./style.css";
import { generateLODProfile } from "./core/lookup-tables";
import { subsolarPoint } from "./core/suncalc";

// English-only for now: force / -> /en and any /xx -> /en
function enforceEnglishPath(): boolean {
    const { pathname, search, hash } = window.location;

    // Root -> /en
    if (pathname === "/") {
        window.location.replace("/en" + search + hash);
        return true;
    }

    // If someone visits /es, /fr, etc. before supported, force /en
    const firstSegment = pathname.split("/").filter(Boolean)[0] ?? "";
    const base = firstSegment.toLowerCase().split("-")[0];

    const looksLikeLang = /^[a-z]{2}$/.test(base);
    if (looksLikeLang && base !== "en") {
        window.location.replace("/en" + search + hash);
        return true;
    }

    // Optional: if user visits /en-US, normalize to /en for now
    if (base === "en" && firstSegment.toLowerCase() !== "en") {
        window.location.replace("/en" + search + hash);
        return true;
    }

    return false;
}

function updateSubsolarPoint() {
    const [lat, lon] = subsolarPoint(generateLODProfile(Date.now()), true); // subsolar lat/lon in degrees
    const ssp = document.getElementById("position");
    if (!ssp) {
        console.warn("Element with id='position' not found");
        return;
    }
    ssp.textContent = `${lat.toFixed(4)}, ${lon.toFixed(4)}`;
}

function startSubsolarUpdates() {
// Align to the next integer second
const now = new Date();
const delay = 1000 - now.getMilliseconds();
setTimeout(() => {
    updateSubsolarPoint();
    setInterval(updateSubsolarPoint, 1000);
}, delay);
}

window.addEventListener("DOMContentLoaded", () => {
    if (enforceEnglishPath()) return;
    startSubsolarUpdates();
});