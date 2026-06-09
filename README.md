# 6N-wing-finestructure

**Part XXXV — Fine Structure of the Wing-Transition Memory on the 6N Skeleton: the Lag-Two Echo and the 1/log Decay**

Ruqing Chen · GUT Geoservice Inc., Montreal · June 2026

Companion code and data for Part XXXV of *Arithmetic Geodynamics on the 6N Skeleton*. A short refinement
of Part XXXIV: it settles two fine questions about the prime wing-transition memory. **Everything here
is a measured count over the primes — no fitted parameters, no fabricated numbers.**

## Background

Part XXXIV showed the prime wing string (R = 6N+1, L = 6N−1) is not first-order Markov but carries an
effective second-order memory: the switch sequence `σ_i = 1[w_{i+1} ≠ w_i]` has lag-1 autocorrelation
≈ −0.033. Two fine questions were left open; we answer both with all 50,847,531 steps over primes
3 < p ≤ 10⁹.

## Findings (measured)

1. **The lag-2 echo is real, not an AR(1) shadow.** For a first-order (AR(1)) switch process the lag-2
   autocorrelation would be exactly `ρ₁² = +0.00098`. Measured: `ρ₂ = +0.00417` — **4.3× larger**. The
   partial autocorrelation `PACF(2) = +0.00319` at **23σ**, with `ρ₃₊` vanishing. So the switch process
   has memory length **exactly two**: a negative nearest-neighbour term and a *positive* lag-2 echo a
   tenth its size, then nothing.

2. **The bias fades as a leading 1/log N.** Across 13 quarter-decade windows, `δ₁ = ½ − P(same)` falls
   0.074 → 0.053 with `δ₁·log N ≈ 1.05` (nearly flat → leading 1/log N); a single power fit gives an
   effective exponent 0.85, the shortfall a slowly-varying (plausibly log log) correction the range
   can't resolve. The second-order signal `|ρ₁|` fades in lock-step: `δ₁/|ρ₁| ≈ 1.7` across five
   decades — first and second order are one fading effect at two orders.

**Scope (honest).** This is a refinement of Part XXXIV, not a new direction. The structure is contained
in the Hardy–Littlewood r-tuple heuristic; the vanishing of the bias is the Oliver–Soundararajan
prediction. **We prove no theorem and propose no mechanism.** What this adds are the specific measured
numbers: the lag-2 partial autocorrelation, the exact-two memory length, the leading 1/log N exponent,
and the locked first-to-second-order ratio.

## Reproducing

```bash
pip install -r requirements.txt
cd code
python3 explore_decay.py   # quick look: lag-2 + decay windows (console)
python3 final_fine.py      # full measurement to N=1e9 -> data/fine_*.csv   (~17 s, ~0.6 GB RAM)
python3 makefigs_fine.py   # reads the CSVs -> figures/p35_fig1.pdf, p35_fig2.pdf
```

The sieve frees its 1 GB boolean array before forming the prime list; paths resolve relative to the
script. NumPy 2.x compliant; single-threaded.

## Files

```
code/    explore_decay.py   final_fine.py   makefigs_fine.py
data/    fine_autocorr.csv    lag, rho, sigma, AR1_prediction        [lags 1-8]
         fine_lag2_test.csv   rho1, rho2, AR1 pred, PACF2, sigma
         fine_decay.csv       lnN, P_same, delta1, abs_rho1, delta1*lnN, ratio
         fine_decay_fit.csv   signal, power_alpha, C
figures/ p35_fig1.pdf  p35_fig2.pdf
paper/   paper35.tex   paper35.pdf
```

All data files are plain CSV.

## Citation

See `CITATION.cff`. The paper is archived on Zenodo (DOI in the citation file once minted).

## License

MIT (see `LICENSE`).
