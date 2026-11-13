"""
P5a - Kick stage preliminary sizing (mass budget + delta-v trade)

Solid-motor kick stage on a rideshare port that carries a 3U CubeSat
deployer and releases the CubeSat after ONE burn at perigee of the
injection orbit. 3-axis attitude control through the burn with a small
cold-gas RCS (decided 12 Nov 2026).

Every number in ASSUMPTIONS is a placeholder to be refined:
  - Isp, thrust    -> P2 (NASA CEA / RocketCEA, grain design)
  - motor fraction -> P1 (CFRP case) + P2 (nozzle)
  - thrust struct. -> P4
  - deployer mass  -> ISISPACE ISIPOD 3U datasheet (2.0 kg)
  - RCS hardware   -> component datasheets (TBD)

Run:  python kick_stage_sizing.py
"""

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

G0 = 9.80665               # m/s^2, standard gravity
MU_EARTH = 3.986004418e14  # m^3/s^2
R_EARTH = 6378.137e3       # m, equatorial radius


# ----------------------------------------------------------------------------
# Assumptions (edit here)
# ----------------------------------------------------------------------------
@dataclass
class Assumptions:
    # Mission
    h_injection_km: float = 500.0   # circular rideshare injection orbit [TBR]
    dv_baseline: float = 500.0      # m/s, single burn [TBR]
    dv_margin: float = 0.05         # 5 % delta-v margin [proposed]

    # Propulsion (solid motor)
    isp_vac: float = 280.0          # s, AP/HTPB/Al small motor [assumed -> P2]
    motor_prop_fraction: float = 0.80  # m_prop / m_motor_loaded [assumed -> P1/P2]
    thrust_nominal: float = 100.0   # N, low thrust to keep RCS small [decided, refine in P2]

    # Burn attitude control (cold-gas RCS)
    cg_offset_m: float = 0.001      # lateral CG offset from thrust axis [assumed]
    misalign_deg: float = 0.25      # nozzle thrust misalignment [assumed]
    nozzle_to_cg_m: float = 0.30    # moment arm of misalignment [assumed]
    rcs_arm_m: float = 0.15         # RCS thruster moment arm [assumed]
    rcs_isp: float = 65.0           # s, nitrogen cold gas [assumed]
    rcs_margin: float = 2.0         # factor on burn-control impulse [proposed]

    # Masses [kg]
    m_cubesat: float = 6.00         # CDS Rev 14.1, Table 1, 3U max
    fixed_dry: dict = field(default_factory=lambda: {
        "deployer_ISIPOD_3U": 2.0,  # ISIPOD datasheet
        "thrust_structure": 0.8,    # [placeholder -> P4]
        "avionics_power_acs": 1.0,  # [placeholder]
        "rcs_hardware": 0.8,        # tank, valves, 4+ thrusters [placeholder]
        "separation_ring_8in": 0.5, # [placeholder -> rideshare UG]
        "harness_misc": 0.3,        # [placeholder]
    })
    dry_margin: float = 0.20        # system margin on stage dry mass [proposed]


# ----------------------------------------------------------------------------
# Physics
# ----------------------------------------------------------------------------
def mass_ratio(dv: float, isp: float) -> float:
    """Tsiolkovsky: m0/mf."""
    return math.exp(dv / (G0 * isp))


def effective_offset(a: Assumptions) -> float:
    """Lever arm [m] of the thrust line about the CG."""
    return a.cg_offset_m + a.nozzle_to_cg_m * math.sin(math.radians(a.misalign_deg))


def rcs_gas_mass(total_impulse: float, a: Assumptions) -> float:
    """
    Cold gas needed to hold attitude through the burn.
    Angular momentum to absorb H = e * I_total (independent of thrust).
    RCS impulse = H / arm; gas = impulse / (Isp * g0), times margin.
    """
    h = effective_offset(a) * total_impulse
    return a.rcs_margin * (h / a.rcs_arm_m) / (a.rcs_isp * G0)


def size_stage(dv: float, a: Assumptions) -> dict:
    """
    Motor inert mass proportional to propellant: m_inert = k * m_p,
    k = (1 - zeta) / zeta. RCS gas depends on total impulse, so iterate.

        m_f = m_cs + (1 + mu) * (m_fixed + k * m_p) + m_gas
        m_p = (R - 1) * m_f
    RCS gas is carried as non-propulsive mass through the burn (conservative).
    """
    R = mass_ratio(dv, a.isp_vac)
    k = (1.0 - a.motor_prop_fraction) / a.motor_prop_fraction
    mu = a.dry_margin
    m_fixed = sum(a.fixed_dry.values())

    m_gas = 0.0
    for _ in range(50):
        denom = 1.0 - (R - 1.0) * (1.0 + mu) * k
        if denom <= 0:
            raise ValueError(f"No closure at dv={dv} m/s: motor too heavy for this Isp.")
        m_p = (R - 1.0) * (a.m_cubesat + (1.0 + mu) * m_fixed + m_gas) / denom
        i_tot = m_p * a.isp_vac * G0
        m_gas_new = rcs_gas_mass(i_tot, a)
        if abs(m_gas_new - m_gas) < 1e-6:
            break
        m_gas = m_gas_new

    m_motor_inert = k * m_p
    m_stage_dry = (1.0 + mu) * (m_fixed + m_motor_inert)
    m_f = a.m_cubesat + m_stage_dry + m_gas
    m_0 = m_f + m_p
    e = effective_offset(a)
    t_dist = a.thrust_nominal * e

    return {
        "dv_m_s": dv,
        "mass_ratio": R,
        "m_prop_kg": m_p,
        "m_motor_inert_kg": m_motor_inert,
        "m_motor_loaded_kg": m_p + m_motor_inert,
        "m_rcs_gas_kg": m_gas,
        "m_stage_dry_kg": m_stage_dry,
        "m_burnout_kg": m_f,
        "m_gross_kg": m_0,
        "total_impulse_kN_s": i_tot / 1e3,
        "burn_time_s": i_tot / a.thrust_nominal,
        "accel_burnout_g": a.thrust_nominal / (m_f * G0),
        "disturbance_torque_Nm": t_dist,
        "momentum_to_absorb_Nms": e * i_tot,
        "rcs_thrust_needed_N": t_dist / a.rcs_arm_m,
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
    m_fixed = sum(a.fixed_dry.values())

    print("=== Kick stage preliminary sizing ===")
    print(f"Injection orbit        : {a.h_injection_km:.0f} km circular")
    print(f"Baseline dv            : {a.dv_baseline:.0f} m/s  (+{a.dv_margin:.0%} -> {dv_design:.0f} m/s design)")
    print(f"Isp (vac) assumed      : {a.isp_vac:.0f} s")
    print(f"Motor prop. fraction   : {a.motor_prop_fraction:.2f}")
    print(f"Nominal thrust         : {a.thrust_nominal:.0f} N")
    print(f"Dry-mass margin        : {a.dry_margin:.0%}\n")

    print("Mass budget at design dv:")
    print(f"  CubeSat 3U              {a.m_cubesat:6.2f} kg")
    for name, m in a.fixed_dry.items():
        print(f"  {name:<24}{m:6.2f} kg")
    print(f"  motor inert             {base['m_motor_inert_kg']:6.2f} kg")
    margin_label = f"margin ({a.dry_margin:.0%} stage dry)"
    print(f"  {margin_label:<24}{base['m_stage_dry_kg'] - m_fixed - base['m_motor_inert_kg']:6.2f} kg")
    print(f"  RCS gas (N2)            {base['m_rcs_gas_kg']:6.2f} kg")
    print(f"  -----------------------------------")
    print(f"  burnout mass            {base['m_burnout_kg']:6.2f} kg")
    print(f"  propellant              {base['m_prop_kg']:6.2f} kg")
    print(f"  GROSS (separation)      {base['m_gross_kg']:6.2f} kg\n")

    print("Burn at design dv:")
    print(f"  total impulse           {base['total_impulse_kN_s']:6.2f} kN*s")
    print(f"  burn time               {base['burn_time_s']:6.0f} s")
    print(f"  accel at burnout        {base['accel_burnout_g']:6.2f} g")
    print(f"  disturbance torque      {base['disturbance_torque_Nm']:6.3f} N*m")
    print(f"  momentum to absorb      {base['momentum_to_absorb_Nms']:6.1f} N*m*s")
    print(f"  RCS thrust per couple   {base['rcs_thrust_needed_N']:6.2f} N (arm {a.rcs_arm_m} m, no margin)\n")

    rows = []
    print(f"{'dv':>6} {'m_p':>7} {'m_motor':>8} {'m_gas':>6} {'m_gross':>8} {'I_tot':>7} {'t_burn':>7} {'apogee':>8} {'di':>6}")
    print(f"{'m/s':>6} {'kg':>7} {'kg':>8} {'kg':>6} {'kg':>8} {'kN*s':>7} {'s':>7} {'km':>8} {'deg':>6}")
    for dv in range(300, 801, 100):
        r = size_stage(dv, a)
        r["apogee_alt_km"] = apogee_after_perigee_burn(a.h_injection_km, dv)
        r["equiv_plane_change_deg"] = plane_change_deg(a.h_injection_km, dv)
        rows.append(r)
        print(f"{dv:6.0f} {r['m_prop_kg']:7.2f} {r['m_motor_loaded_kg']:8.2f} {r['m_rcs_gas_kg']:6.2f} "
              f"{r['m_gross_kg']:8.2f} {r['total_impulse_kN_s']:7.2f} {r['burn_time_s']:7.0f} "
              f"{r['apogee_alt_km']:8.0f} {r['equiv_plane_change_deg']:6.2f}")

    csv_path = out_dir / "dv_trade.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nTrade table written to {csv_path}")


if __name__ == "__main__":
    main()
