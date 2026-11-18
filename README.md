# Static liquefaction strength ratio calculator

Browser-based calculator for undrained **yield** strength ratio, su(yield)/σ′vo, and **liquefied** (post-liquefaction) strength ratio, su(liq)/σ′vo, of cohesionless soils from SPT or CPT penetration resistance.

**Live tool:** open `index.html` from the GitHub Pages site of this repository.

## What it does
- Implements the procedure of **Sadrekarimi (2016)**, which accounts for the mode of shearing (compression, simple shear, extension) and slope-induced anisotropy, Kc = (1 − sin α)/(1 + sin α), for an infinite slope of angle α.
- Reports both strength ratios for all three modes, with ranges reflecting the scatter of the field correlation. A mode is reported as N/A (not liquefiable) when the strain-softening potential IB < 0.
- For comparison, gives liquefied strength ratios from Olson & Stark (2002), Idriss & Boulanger (2015) and Robertson (2010, CPT only), and the yield strength ratio of Olson & Stark (2003).
- Clean sand (FC = 0%) is assumed throughout.

## Files
- `index.html` — the web calculator (self-contained; open in any browser).
- `strength_ratios.py` — Python implementation of the same calculations. Edit the `INPUTS` block and run `python strength_ratios.py`.

## References
- Sadrekarimi, A. (2016). Static liquefaction analysis considering principal stress directions and anisotropy. *Geotechnical and Geological Engineering*, 34, 1135–1154.
- Olson, S.M., and Stark, T.D. (2002). Liquefied strength ratio from liquefaction flow failure case histories. *Canadian Geotechnical Journal*, 39, 629–647.
- Olson, S.M., and Stark, T.D. (2003). Yield strength ratio and liquefaction analysis of slopes and embankments. *Journal of Geotechnical and Geoenvironmental Engineering*, 129(8), 727–737.
- Idriss, I.M., and Boulanger, R.W. (2015). 2nd Ishihara Lecture: SPT- and CPT-based relationships for the residual shear strength of liquefied soils. *Soil Dynamics and Earthquake Engineering*, 68, 57–68.
- Robertson, P.K. (2010). Evaluation of flow liquefaction and liquefied strength using the cone penetration test. *Journal of Geotechnical and Geoenvironmental Engineering*, 136(6), 842–853.

## Disclaimer
Intended for screening, teaching and comparison of methods. It is not a substitute for site-specific investigation, laboratory testing and engineering judgment.
