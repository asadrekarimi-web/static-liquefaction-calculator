"""
Static liquefaction strength ratios: Sadrekarimi (2016) compared with other methods.

Undrained yield, su(yield)/s'vo, and liquefied, su(liq)/s'vo, strength ratios from
Sadrekarimi (2016) for compression, simple shear and extension modes of shearing,
accounting for slope-induced anisotropy Kc, compared with Olson & Stark (2002),
Idriss & Boulanger (2015) and Robertson (2010). Clean sand (FC = 0 %) is assumed.

When the strength-softening potential IB < 0 for a mode, the soil is taken as
not liquefiable and results for that mode are reported as N/A.

Edit the INPUTS block and run:  python strength_ratios.py

References
  Sadrekarimi, A. (2016) Geotech. Geol. Eng. 34:1135-1154
  Olson, S.M. & Stark, T.D. (2002) Can. Geotech. J. 39:629-647
  Olson, S.M. & Stark, T.D. (2003) J. Geotech. Geoenviron. Eng. 129(8):727-737  (yield, reference only)
  Idriss, I.M. & Boulanger, R.W. (2015) Soil Dyn. Earthq. Eng. 68:57-68
  Robertson, P.K. (2010) J. Geotech. Geoenviron. Eng. 136(6):842-853
"""
import math

# ============================================================================
# INPUTS
# ============================================================================
INPUTS = dict(
    test="SPT",           # "SPT" or "CPT" - only this test type is used in all methods
    N1_60=6.0,            # (N1)60, blows/ft   (used when test = "SPT")
    qc1_MPa=3.0,          # qc1, MPa           (used when test = "CPT")
    phi_deg=33.0,         # tan(phi') cap for Idriss & Boulanger and Robertson
    # --- Sadrekarimi (2016) ---
    slope_deg=30.0,       # slope angle a, degrees -> Kc = (1 - sin a)/(1 + sin a)
    Kc_override=None,     # give a value to bypass the slope-angle estimate
    mode="simple shear",  # "compression", "simple shear" or "extension" - reference for comparison
    b_comp=0.3, b_ss=0.3, b_ext=1.0,                    # post-liquefaction b
    phi_cs_comp=35.0, phi_cs_ss=35.0, phi_cs_ext=30.0,  # critical-state friction angles, deg
)

PA_MPA = 0.1013  # atmospheric pressure, MPa
MODES = ("compression", "simple shear", "extension")


def tan_d(deg):
    return math.tan(math.radians(deg))


# ============================================================================
# Other methods (FC = 0 %)
# ============================================================================
def olson_stark_liq(test, x):
    """Olson & Stark (2002) liquefied strength ratio, +/- 0.03."""
    if test == "SPT":
        return 0.03 + 0.0075 * x, ("(N1)60 > 12: beyond data range" if x > 12 else "")
    return 0.03 + 0.0143 * x, ("qc1 > 6.5 MPa: beyond data range" if x > 6.5 else "")


def olson_stark_yield(test, x):
    """Olson & Stark (2003) yield strength ratio, +/- 0.04 (reference only)."""
    return 0.205 + (0.0075 * x if test == "SPT" else 0.0143 * x)


def ib_spt(N):          # void redistribution significant
    return math.exp(N / 16 + ((N - 16) / 21.2) ** 3 - 3.0)


def ib_spt_neg(N):      # void redistribution negligible
    return ib_spt(N) * (1 + math.exp(N / 2.4 - 6.6))


def ib_cpt(q):          # void redistribution significant
    return math.exp(q / 24.5 - (q / 61.7) ** 2 + (q / 106) ** 3 - 4.42)


def ib_cpt_neg(q):      # void redistribution negligible
    return ib_cpt(q) * (1 + math.exp(q / 11.1 - 9.82))


def robertson_liq(Q, phi_deg):
    """Robertson (2010); FC = 0 -> Kc = 1, Qtn,cs = Qtn ~ qc1/pa. Valid for Qtn,cs <= 70."""
    if Q > 70:
        return None, "Qtn,cs > 70: likely dilative; correlation not applicable"
    r = (0.02199 - 0.0003124 * Q) / (1 - 0.02676 * Q + 0.0001783 * Q ** 2)
    return min(max(r, 0.03), tan_d(phi_deg)), "Qtn taken as qc1/pa"


def other_methods(test, x, phi_deg):
    cap = tan_d(phi_deg)
    rows = []
    v, note = olson_stark_liq(test, x)
    rows.append(("Olson & Stark (2002)", f"{test} correlation", v, note))
    if test == "SPT":
        f_sig, f_neg, arg = ib_spt, ib_spt_neg, x
        rows_b = ("void redistribution significant", "void redistribution negligible")
    else:
        f_sig, f_neg, arg = ib_cpt, ib_cpt_neg, x / PA_MPA
        rows_b = (f"qc1Ncs-Sr = {arg:.1f}, void redistribution significant",
                  f"qc1Ncs-Sr = {arg:.1f}, void redistribution negligible")
    for f, lab in zip((f_sig, f_neg), rows_b):
        r = f(arg)
        rows.append(("Idriss & Boulanger (2015)", lab, min(r, cap), "capped at tan(phi')" if r > cap else ""))
    if test == "SPT":
        rows.append(("Robertson (2010)", "CPT-only method", None, "not applicable to SPT"))
    else:
        v, note = robertson_liq(x / PA_MPA, phi_deg)
        rows.append(("Robertson (2010)", f"Qtn,cs = {x / PA_MPA:.1f}", v, note))
    return rows


# ============================================================================
# Sadrekarimi (2016)
# ============================================================================
# su(liq)/s'1c = a IB^2 + b IB + c for each mode (laboratory relationships)
MODE_CURVE = {
    "compression": (0.045, -0.320, 0.275),
    "simple shear": (0.053, -0.309, 0.251),
    "extension": (0.120, -0.345, 0.228),
}


def kc_from_slope(alpha_deg):
    s = math.sin(math.radians(alpha_deg))
    return (1 - s) / (1 + s)


def ib_from_su(su, mode):
    a, b, c = MODE_CURVE[mode]
    disc = b * b - 4 * a * (c - su)
    return float("nan") if disc < 0 else (-b - math.sqrt(disc)) / (2 * a)


def su_from_ib(IB, mode):
    a, b, c = MODE_CURVE[mode]
    return a * IB * IB + b * IB + c


def sadrekarimi_mode(su, Kc, mode, b, phi_cs):
    """Returns dict of ru, r*u, I*B, IB, su(liq), su(yield); na=True when IB < 0."""
    t = tan_d(45 + phi_cs / 2) ** 2
    root = math.sqrt(b * b - b + 1)
    if mode == "compression":
        ru = Kc - 2 / ((t - 1) * root) * su
    else:
        ru = Kc - 2 * (1 + b * t - b) / ((t - 1) * root) * su
    k = (1 - Kc) / math.sqrt(Kc)
    if mode == "compression":
        rs = ru / (1 - 0.60 * (1 - Kc))
        IBs = 0.553 * rs ** 2 + 1.621 * rs - 1.214
        IB = IBs + k
    elif mode == "simple shear":
        aH = 0.330 + 9.089 * Kc - 11.376 * b - 5.833 * Kc ** 2 + 15.333 * b ** 2 - 5.848 * b * Kc
        bH = -0.577 * b ** 2 + 1.117 * b
        rs = ru / (1 - aH * (1 - Kc)) + 0.85 * bH * Kc * (1 - ru)
        IBs = 1.003 * rs ** 2 + 0.908 * rs - 0.963
        IB = IBs + 0.80 * k
    else:
        rs = ru / Kc
        IBs = 0.959 * rs ** 2 + 0.076 * rs - 0.038
        IB = IBs - 0.1 * k
    if IB < 0:
        return dict(na=True, flag="Not liquefiable (IB < 0)")
    flags = []
    yld = su / (1 - IB) if IB < 1 else float("inf")
    if IB >= 1:
        flags.append("IB >= 1: outside correlation")
    if ru < 0:
        flags.append("ru < 0: su(liq) high for this Kc")
    return dict(na=False, ru=ru, rs=rs, IBs=IBs, IB=IB, liq=su, yld=yld, flag="; ".join(flags))


def sadrekarimi_2016(test, x, slope_deg=30.0, Kc_override=None, b_comp=0.3, b_ss=0.3, b_ext=1.0,
                     phi_cs_comp=35.0, phi_cs_ss=35.0, phi_cs_ext=30.0):
    Kc = Kc_override if Kc_override is not None else kc_from_slope(slope_deg)
    if test == "SPT":       # simple-shear field correlation
        ss, tol, valid = 0.046 + 0.006 * x, 0.017, 3 < x < 12
    else:
        ss, tol, valid = 0.022 + 0.017 * x, 0.021, 0.5 < x < 8
    pars = {"compression": (b_comp, phi_cs_comp), "simple shear": (b_ss, phi_cs_ss),
            "extension": (b_ext, phi_cs_ext)}
    modes = {}
    for case, s0 in (("best", ss), ("lo", ss - tol), ("hi", ss + tol)):
        IB_ss = ib_from_su(s0, "simple shear")   # same strain-softening potential for all modes
        for m, (b, phi) in pars.items():
            su = s0 if m == "simple shear" else su_from_ib(IB_ss, m)
            r = sadrekarimi_mode(su, Kc, m, b, phi)
            r["IB_map"] = IB_ss
            modes.setdefault(m, {})[case] = r
    return dict(Kc=Kc, valid=valid, modes=modes, pars=pars)


# ============================================================================
# Output
# ============================================================================
def fmt(x, nd=3):
    if x is None:
        return "-"
    if isinstance(x, float) and (math.isinf(x) or math.isnan(x)):
        return "N/A"
    return f"{x:.{nd}f}"


def rng(M, key):
    v = [M[c][key] for c in ("lo", "best", "hi") if not M[c]["na"] and math.isfinite(M[c][key])]
    return f"{min(v):.3f}-{max(v):.3f}" if v else "N/A"


def run(p):
    test = p["test"].upper()
    x = p["N1_60"] if test == "SPT" else p["qc1_MPa"]
    sk = sadrekarimi_2016(test, x, p["slope_deg"], p["Kc_override"], p["b_comp"], p["b_ss"],
                          p["b_ext"], p["phi_cs_comp"], p["phi_cs_ss"], p["phi_cs_ext"])
    xl = f"(N1)60 = {x}" if test == "SPT" else f"qc1 = {x} MPa"

    print(f"\nSADREKARIMI (2016)   {xl}   Kc = {sk['Kc']:.3f}"
          + ("" if sk["valid"] else "   ** outside correlation range **"))
    print(f"{'Mode':<14}{'b':>5}{'phi_cs':>8}{'su(liq)/svo':>13}{'range':>14}{'ru':>8}{'r*u':>8}"
          f"{'I*B':>8}{'IB':>8}{'su(yield)/svo':>15}{'range':>14}  Notes")
    for m in MODES:
        M, b = sk["modes"][m], sk["modes"][m]["best"]
        bb, ph = sk["pars"][m]
        tag = "  <== reference" if m == p["mode"] else ""
        if b["na"]:
            print(f"{m:<14}{bb:>5}{ph:>8}" + "".join(f"{'N/A':>{w}}" for w in (13, 14, 8, 8, 8, 8, 15, 14))
                  + f"  {b['flag']}{tag}")
            continue
        print(f"{m:<14}{bb:>5}{ph:>8}{fmt(b['liq']):>13}{rng(M, 'liq'):>14}{fmt(b['ru']):>8}"
              f"{fmt(b['rs']):>8}{fmt(b['IBs']):>8}{fmt(b['IB']):>8}{fmt(b['yld']):>15}"
              f"{rng(M, 'yld'):>14}  {b['flag']}{tag}")
    print(f"IB used to map simple-shear su(liq) to compression/extension: "
          f"{sk['modes']['simple shear']['best']['IB_map']:.3f}")

    ref = sk["modes"][p["mode"]]["best"]
    r_liq = None if ref["na"] else ref["liq"]
    r_yld = None if ref["na"] else ref["yld"]

    def d(v, r):
        if v is None or r is None or not math.isfinite(r):
            return "-", "-"
        if abs(v - r) < 5e-4:
            return "0.000", "0%"
        return f"{v - r:+.3f}", f"{100 * (v - r) / r:+.0f}%"

    print(f"\nCOMPARISON vs Sadrekarimi (2016), {p['mode']}   "
          f"[su(liq)/svo = {fmt(r_liq) if r_liq is not None else 'N/A'}, "
          f"su(yield)/svo = {fmt(r_yld) if r_yld is not None else 'N/A'}]"
          + ("   -> not liquefiable (IB < 0)" if ref["na"] else ""))
    print(f"{'Method':<27}{'Basis':<52}{'ratio':>8}{'delta':>9}{'%diff':>8}  Notes")
    print("-" * 116)
    for meth, basis, v, note in other_methods(test, x, p["phi_deg"]):
        dd, pc = d(v, r_liq)
        print(f"{meth:<27}{('liq: ' + basis)[:51]:<52}{fmt(v):>8}{dd:>9}{pc:>8}  {note}")
    y = olson_stark_yield(test, x)
    dd, pc = d(y, r_yld)
    print(f"{'Olson & Stark (2003)*':<27}{'yield (+/-0.04, reference only)':<52}{fmt(y):>8}{dd:>9}{pc:>8}")
    print("delta = other - Sadrekarimi;  % diff relative to Sadrekarimi")
    return sk


if __name__ == "__main__":
    run(INPUTS)
