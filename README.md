# Solid-motor kick stage for a 3U CubeSat — design portfolio

A design portfolio built around one vehicle: a free-flying solid-motor kick
stage that rides a rideshare port and moves a 3U CubeSat to a different
orbital plane in a single burn.

| Project | Topic | Folder |
|---|---|---|
| P5a | System requirements and preliminary sizing | `sizing/` |
| P2 | Nozzle with ablative liner | `nozzle/` |
| P1 | CFRP filament-wound motor case | `motor-case/` |
| P4 | Topology-optimized thrust structure | `thrust-structure/` |
| P3 | 3U CubeSat structure with sandwich panels | `cubesat-structure/` |
| P5b | System integration | `integration/` |

## P5a — preliminary sizing

```bash
cd sizing
python kick_stage_sizing.py
```

Every assumption sits at the top of the script in the `Assumptions` class,
marked with the project that will replace it. The delta-v trade is written to
`sizing/out/dv_trade.csv`.

```bash
python orbit_lifetime.py
```

Screens orbital lifetime against the FCC 5-year post-mission disposal rule.
It is why the mission is a plane change and not an apogee raise: a single
solid burn that raises apogee from 500 km leaves the spent stage in orbit for
decades, while a plane change at 500 km circular lets it re-enter in about
1–4 years.

Baseline (TBR): 500 km circular injection, single 500 m/s burn at a node
(525 m/s with margin, ~3.9° inclination change), Isp 280 s, motor
propellant fraction 0.72 (STAR 6A heritage; the CFRP case in P1 targets
≥ 0.80), 100 N nominal thrust, 3-axis control through the burn with a small
nitrogen cold-gas RCS
→ 3.0 kg propellant, ~83 s burn, 17.3 kg gross.
