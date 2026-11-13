# Solid-motor kick stage for a 3U CubeSat — design portfolio

A design portfolio built around one vehicle: a free-flying solid-motor kick
stage that rides a rideshare port and delivers a 3U CubeSat to its final
orbit in a single burn.

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

Baseline (TBR): 500 km circular injection, 500 m/s single burn
(525 m/s with margin), Isp 280 s → 2.6 kg propellant, 14.9 kg gross.
