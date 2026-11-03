# Boundary Layer Lab

An interactive, browser-based simulation of the **flat-plate boundary layer** — a teaching and exploration tool for fluid mechanics students.

**Live demo:** [boundary-layer-lab.vercel.app](https://boundary-layer-lab.vercel.app)

![Fluids covered](https://img.shields.io/badge/fluids-air%20%2F%20water%20%2F%20glycerin-blue)

---

## About

When a fluid flows over a flat plate at zero pressure gradient, viscosity slows the fluid near the wall, forming a thin *boundary layer* that grows along the plate. Boundary Layer Lab visualizes this phenomenon in real time and computes the standard engineering solutions:

- **Laminar region** — the exact **Blasius** similarity solution (tabulated profile, interpolated in code)
- **Turbulent region** — the **1/7-power-law** velocity profile
- **Transition** at a user-adjustable critical Reynolds number (default `Re_x = 5×10⁵`)

The simulation shows the layer growth along the plate (δ envelope), the velocity profile at any inspect station, tracers advected with the flow, and standard boundary-layer readouts (δ, δ*, θ, shape factor H, skin friction C_f, local and total Re).

## Features

- **Three fluids** — air, water, and glycerin, each with its own kinematic viscosity ν and density ρ
- **Adjustable parameters** — free-stream speed `U∞`, plate length `L`, critical Reynolds number, and inspect station `x/L`
- **Three views**
  - *Simulation* — animated flow bench with tracers and layer growth
  - *Equations & Calculations* — step-by-step working for the current parameters
  - *Theory & Background* — the experiment, scientific principles, and the scientists behind the solutions
- **Toggleable overlays** — tracers, profiles, δ envelope, and the δ_Pohlhausen cubic approximation
- **Demo presets** — laminar only, transition on plate, water/slow
- Click anywhere on the plate to inspect that station

## The equations used

| Quantity | Laminar (exact) | Turbulent (approx.) |
|---|---|---|
| δ(x) | `5.0x / √Re_x` | `0.37x / Re_x^{1/5}` |
| δ* / θ | from Blasius profile (0.344δ / 0.133δ) | 1/7-power law |
| C_f | `0.664 / √Re_x` | `0.059 / Re_x^{1/5}` |

## Running locally

No build step or dependencies — it's plain HTML/CSS/JS.

```bash
git clone https://github.com/sahilchoudhary2006/Boundary_layer-Fluid-Mechanics-.git
cd Boundary_layer-Fluid-Mechanics-
# open index.html in your browser, or serve it:
python -m http.server 8000
```

Then visit `http://localhost:8000`.

## Project structure

```
index.html   — markup, controls, and static sections
style.css    — layout and styling
script.js    — simulation: Blasius interpolation, profiles, canvas rendering
```

## Tech

- Vanilla JavaScript + Canvas (no frameworks, no build tools)
- Google Fonts (Chivo / Chivo Mono)
- Deployed on Vercel

## Credits

Built as a fluid mechanics learning project. Boundary-layer theory per Schlichting,
*Boundary-Layer Theory*.
