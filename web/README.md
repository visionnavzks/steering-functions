# Steering Functions — Web UI

An interactive browser-based demo for the **steering-functions** library.
The path algorithms (Dubins and Reeds-Shepp) are implemented directly in JavaScript
from the same mathematical foundations as the C++ library.

## Quick Start

Open `index.html` in any modern browser — no build step or server required.

```bash
# Optional: serve with a local HTTP server
python3 -m http.server 8080
# then open http://localhost:8080
```

## Features

| Feature | Description |
|---------|-------------|
| **Dubins paths** | 6 path types: LSL, RSR, RSL, LSR, RLR, LRL (forward motion only) |
| **Reeds-Shepp paths** | All 48 candidate paths across 18 type families (forward + backward motion) |
| **Interactive canvas** | Click to place start/goal, drag to set heading |
| **Pan & zoom** | Mouse wheel to zoom; Alt+drag or middle-button drag to pan |
| **All-paths view** | Toggle between shortest path only or all candidates |
| **Path info panel** | Shows path type, total length, and per-segment breakdown |
| **Color coding** | Green segments = forward motion · Red/orange segments = backward motion |
| **Cusp markers** | Yellow dots mark direction-reversal points in Reeds-Shepp paths |

## Controls

### Keyboard shortcuts

| Key | Action |
|-----|--------|
| `S` / `1` | Switch to **Start** placement mode |
| `G` / `2` | Switch to **Goal** placement mode |
| `D` | Switch to **Dubins** algorithm |
| `E` | Switch to **Reeds-Shepp** algorithm |
| `A` | Toggle **All paths** / **Shortest path** display |
| `R` | **Reset** everything to defaults |

### Mouse

- **Click** on the canvas → set the selected state's position
- **Click + drag** → set orientation (angle from click point to mouse)
- **Scroll wheel** → zoom in/out (anchored at cursor)
- **Alt + drag** or **Middle-button drag** → pan the viewport

## Algorithm Notes

Both implementations are direct JavaScript ports of the C++ source in
`src/steering_functions/src/`, using the same normalization, formula indices,
and path-type tables:

- **Dubins** — formulas from Dubins (1957); normalised distance + angle coordinates
- **Reeds-Shepp** — formulas 8.1–8.11 from Reeds & Shepp (1990); CSC, CCC, CCCC,
  CCSC, and CCSCC path families with time-flip and reflect symmetries

The turning radius slider (20–250 px) directly corresponds to `1/κ_max` in the C++
library. Negative `delta_s` values in the control sequence indicate backward motion.
