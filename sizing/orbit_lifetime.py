"""
P5a - Orbital lifetime estimate (DEB-01)

Orbit-averaged drag decay of an elliptical orbit using Gauss's equations
with a purely tangential drag acceleration, a non-rotating atmosphere and
a piecewise-exponential density model (CIRA-72 based table as tabulated in
Vallado, "Fundamentals of Astrodynamics and Applications", Table 8-4,
mean solar activity).

Solar activity changes upper-atmosphere density by large factors, so the
result is reported for density multipliers that bracket low, mean and high
activity. This is an order-of-magnitude screening tool, not a replacement
for NASA DAS / ESA DRAMA / STELA.

Result (13 Nov 2026): the apogee-raise mission (500 x 2,660 km) leaves the
spent stage in orbit for decades to centuries, far beyond the FCC 5-year
post-mission disposal rule, so the mission was changed to a plane change
that keeps the stage at 500 km circular.

Run:  python orbit_lifetime.py
"""

import math

MU = 3.986004418e14   # m^3/s^2
RE = 6378.137e3       # m

# base altitude [km], base density [kg/m^3], scale height [km]
DENSITY_TABLE = [
    (100, 5.297e-7, 5.877), (110, 9.661e-8, 7.263), (120, 2.438e-8, 9.473),
    (130, 8.484e-9, 12.636), (140, 3.845e-9, 16.149), (150, 2.070e-9, 22.523),
    (180, 5.464e-10, 29.740), (200, 2.789e-10, 37.105), (250, 7.248e-11, 45.546),
    (300, 2.418e-11, 53.298), (350, 9.518e-12, 53.628), (400, 3.725e-12, 58.515),
    (450, 1.585e-12, 60.828), (500, 6.967e-13, 63.822), (600, 1.454e-13, 71.835),
    (700, 3.614e-14, 88.667), (800, 1.170e-14, 124.64), (900, 5.245e-15, 181.05),
    (1000, 3.019e-15, 268.00),
]


def density(h_km: float) -> float:
    if h_km < 100:
        h0, rho0, H = DENSITY_TABLE[0]
    else:
        h0, rho0, H = DENSITY_TABLE[0]
        for row in DENSITY_TABLE:
            if h_km >= row[0]:
                h0, rho0, H = row
    return rho0 * math.exp(-(h_km - h0) / H)


def averaged_rates(a: float, e: float, B: float, k_rho: float, n: int = 720):
    """Orbit-averaged da/dt, de/dt [m/s, 1/s] for tangential drag."""
    p = a * (1 - e * e)
    h = math.sqrt(MU * p)
    T = 2 * math.pi * math.sqrt(a**3 / MU)
    da = de = 0.0
    dnu = 2 * math.pi / n
    for i in range(n):
        nu = (i + 0.5) * dnu
        r = p / (1 + e * math.cos(nu))
        v = math.sqrt(MU * (2 / r - 1 / a))
        rho = k_rho * density((r - RE) / 1e3)
        aT = -0.5 * B * rho * v * v            # tangential acceleration
        dt = r * r / h * dnu
        da += 2 * a * a * v * aT / MU * dt
        de += 2 * (e + math.cos(nu)) / v * aT * dt
    return da / T, de / T


def lifetime_years(hp_km: float, ha_km: float, m: float, cd: float, area: float,
                   k_rho: float = 1.0, max_years: float = 200.0) -> float:
    B = cd * area / m
    rp, ra = RE + hp_km * 1e3, RE + ha_km * 1e3
    a, e = (rp + ra) / 2, (ra - rp) / (ra + rp)
    t = 0.0
    while t < max_years * 365.25 * 86400:
        hp = (a * (1 - e) - RE) / 1e3
        if hp < 120:
            return t / (365.25 * 86400)
        da, de = averaged_rates(a, e, B, k_rho)
        # step so that perigee altitude drops by at most ~1 km (or 30 days)
        dhp = abs(da * (1 - e) - a * de)
        dt = min(30 * 86400, 1e3 / dhp if dhp > 0 else 30 * 86400)
        a += da * dt
        e = max(0.0, e + de * dt)
        t += dt
    return float("inf")


def main() -> None:
    cases = {
        # name: (mass kg, Cd, mean tumbling area m^2)
        "Spent stage (8.3 kg, 0.10 m^2)": (8.3, 2.2, 0.10),
        "3U CubeSat (6.0 kg, 0.045 m^2)": (6.0, 2.2, 0.045),
    }
    orbits = {
        "500 km circular (baseline: plane change)": (500, 500),
        "500 x 2,660 km (rejected: apogee raise)": (500, 2660),
        "400 x 2,500 km (rejected: lower perigee)": (400, 2500),
    }
    activity = {"low (x0.5)": 0.5, "mean (x1)": 1.0, "high (x2.5)": 2.5}

    print("Orbital lifetime [years], perigee decay to 120 km")
    for oname, (hp, ha) in orbits.items():
        print(f"\n{oname}")
        for cname, (m, cd, A) in cases.items():
            vals = [lifetime_years(hp, ha, m, cd, A, k) for k in activity.values()]
            txt = "  ".join(f"{k}: {v:6.1f}" for k, v in zip(activity, vals))
            print(f"  {cname:<32} {txt}")


if __name__ == "__main__":
    main()
