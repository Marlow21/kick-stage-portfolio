"""
P5a - Kick stage preliminary sizing (mass budget + delta-v trade)

Solid-motor kick stage on a rideshare port that carries a 3U CubeSat
deployer and releases the CubeSat after ONE burn at perigee of the
injection orbit.

Every number in ASSUMPTIONS is a placeholder to be refined:
  - Isp            -> P2 (NASA CEA / RocketCEA)
  - motor fraction -> P1 (CFRP case) + P2 (nozzle)
  - thrust struct. -> P4
  - deployer mass  -> datasheet of the chosen commercial deployer

Run:  python kick_stage_sizing.py
"""

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

G0 = 9.80665            # m/s^2, standard gravity
MU_EARTH = 3.986004418e14  # m^3/s^2
R_EARTH = 6378.137e3    # m, equatorial radius


# ----------------------------------------------------------------------------
# Assumptions (edit here)
# ----------------------------------------------------------------------------
@dataclass
class Assumptions:
    # Mission
    h_injection_km: float = 500.0   # circular rideshare injection orbit [TBR]
    dv_baseline: float = 500.0      # m/s, single burn [TBR]
    dv_margin: float = 0.05         # 5 % delta-v margin [proposed]

    # Propulsion
    isp_vac: float = 280.0          # s, AP/HTPB/Al small motor [assumed -> P2]
    motor_prop_fraction: float = 0.80  # m_prop / m_motor_loaded [assumed -> P1/P2]
    a_max_g: float = 5.0            # max axial accel on CubeSat during burn [TBD]

    # Masses [kg]
    m_cubesat: float = 6.00         # CDS Rev 14.1, Table 1, 3U max
    fixed_dry: dict = field(default_factory=lambda: {
        "deployer_3U": 2.0,         # [placeholder -> datasheet]
        "thrust_structure": 0.8,    # [placeholder -> P4]
        "avionics_power_acs": 1.0,  # [placeholder]
        "separation_ring_flight": 0.5,  # [placeholder -> rideshare UG]
        "harness_misc": 0.3,        # [placeholder]
    })
    dry_margin: float = 0.20        # system margin on stage dry mass [proposed]


# ----------------------------------------------------------------------------
# Physics
# ----------------------------------------------------------------------------
def mass_ratio(dv: float, isp: float) -> float:
    """Tsiolkovsky: m0/mf."""
    return math.exp(dv / (G0 * isp))


def size_stage(dv: float, a: Assumptions) -> dict:
    """
    Closed-form propellant mass with motor inert mass proportional to
    propellant:  m_inert_motor = k * m_p,  k = (1 - zeta) / zeta.

        m_f = m_cs + (1 + mu) * (m_fixed + k * m_p)
        m_p = (R - 1) * m_f
    =>  m_p = (R - 1)(m_cs + (1 + mu) m_fixed) / (1 - (R - 1)(1 + mu) k)
    """
    R = mass_ratio(dv, a.isp_vac)
    k = (1.0 - a.motor_prop_fraction) / a.motor_prop_fraction
    mu = a.dry_margin
    m_fixed = sum(a.fixed_dry.values())

    denom = 1.0 - (R - 1.0) * (1.0 + mu) * k
    if denom <= 0:
        raise ValueError(f"No closure at dv={dv} m/s: motor too heavy for this Isp.")

    m_p = (R - 1.0) * (a.m_cubesat + (1.0 + mu) * m_fixed) / denom
    m_motor_inert = k * m_p
    m_stage_dry = (1.0 + mu) * (m_fixed + m_motor_inert)
    m_f = a.m_cubesat + m_stage_dry
    m_0 = m_f + m_p
    f_max = m_f * a.a_max_g * G0   # thrust ceiling set at burnout (lightest)

    return {
        "dv_m_s": dv,
        "mass_ratio": R,
        "m_prop_kg": m_p,
        "m_motor_inert_kg": m_motor_inert,
        "m_motor_loaded_kg": m_p + m_motor_inert,
        "m_stage_dry_kg": m_stage_dry,
        "m_burnout_kg": m_f,
        "m_gross_kg": m_0,
        "total_impulse_kN_s": m_p * a.isp_vac * G0 / 1e3,
        "thrust_max_N": f_max,
    }


def apogee_after_perigee_burn(h_circ_km: float, dv: float) -> float:
    """Apogee altitude [km] after a tangential burn from a circular orbit."""
    r_p = R_EARTH + h_circ_km * 1e3
    v_p = math.sqrt(MU_EARTH / r_p) + dv
    inv_a = 2.0 / r_p - v_p**2 / MU_EARTH
    if inv_a <= 0:
        return float("inf")  # escape
    a_sma = 1.0 / inv_a
    r_a = 2.0 * a_sma - r_p
    return (r_a - R_EARTH) / 1e3


def plane_change_deg(h_circ_km: float, dv: float) -> float:
    """Pure inclination change [deg] the same dv buys at circular speed."""
    v = math.sqrt(MU_EARTH / (R_EARTH + h_circ_km * 1e3))
    return math.degrees(2.0 * math.asin(min(1.0, dv / (2.0 * v))))


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
def main() -> None:
    a = Assumptions()
    out_dir = Path(__file__).parent / "out"
    out_dir.mkdir(exist_ok=True)

    dv_design = a.dv_baseline * (1.0 + a.dv_margin)
    base = size_stage(dv_design, a)

    print("=== Kick stage preliminary sizing ===")
    print(f"Injection orbit        : {a.h_injection_km:.0f} km circular")
    print(f"Baseline dv            : {a.dv_baseline:.0f} m/s  (+{a.dv_margin:.0%} -> {dv_design:.0f} m/s design)")
    print(f"Isp (vac) assumed      : {a.isp_vac:.0f} s")
    print(f"Motor prop. fraction   : {a.motor_prop_fraction:.2f}")
    print(f"Dry-mass margin        : {a.dry_margin:.0%}\n")

    print("Mass budget at design dv:")
    print(f"  CubeSat 3U            {a.m_cubesat:6.2f} kg")
    for name, m in a.fixed_dry.items():
        print(f"  {name:<22}{m:6.2f} kg")
    print(f"  motor inert           {base['m_motor_inert_kg']:6.2f} kg")
    print(f"  margin ({a.dry_margin:.0%} of stage dry){base['m_stage_dry_kg'] - sum(a.fixed_dry.values()) - base['m_motor_inert_kg']:6.2f} kg")
    print(f"  ---------------------------------")
    print(f"  burnout mass          {base['m_burnout_kg']:6.2f} kg")
    print(f"  propellant            {base['m_prop_kg']:6.2f} kg")
    print(f"  GROSS (liftoff)       {base['m_gross_kg']:6.2f} kg")
    print(f"  total impulse         {base['total_impulse_kN_s']:6.2f} kN*s")
    print(f"  thrust ceiling ({a.a_max_g:.0f} g) {base['thrust_max_N']:6.0f} N\n")

    rows = []
    print(f"{'dv':>6} {'m_p':>7} {'m_motor':>8} {'m_gross':>8} {'I_tot':>8} {'F_max':>7} {'apogee':>8} {'di':>6}")
    print(f"{'m/s':>6} {'kg':>7} {'kg':>8} {'kg':>8} {'kN*s':>8} {'N':>7} {'km':>8} {'deg':>6}")
    for dv in range(300, 801, 100):
        r = size_stage(dv, a)
        r["apogee_alt_km"] = apogee_after_perigee_burn(a.h_injection_km, dv)
        r["equiv_plane_change_deg"] = plane_change_deg(a.h_injection_km, dv)
        rows.append(r)
        print(f"{dv:6.0f} {r['m_prop_kg']:7.2f} {r['m_motor_loaded_kg']:8.2f} {r['m_gross_kg']:8.2f} "
              f"{r['total_impulse_kN_s']:8.2f} {r['thrust_max_N']:7.0f} {r['apogee_alt_km']:8.0f} "
              f"{r['equiv_plane_change_deg']:6.2f}")

    csv_path = out_dir / "dv_trade.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nTrade table written to {csv_path}")


if __name__ == "__main__":
    main()
