# Fluidic Learning Channel (FLC) - 2x2 Network Re-Implementation

Re-implementation from the author's reference settings of the `2x2` FLC network simulations in R. Wu, L. Y. Lee, S. Terrile, and H. Hauser (*IEEE RoboSoft 2026*, [DOI: 10.1109/RoboSoft67810.2026.11522860](https://doi.org/10.1109/RoboSoft67810.2026.11522860)), matching the author's Python port ([wurui1991/flc-reproduction-notes](https://github.com/wurui1991/flc-reproduction-notes)).

## Quick Start (Tested with Python 3.13, Standard Library Only)

```bash
python reproduce_bristol_flc_clean.py
```

Uses only the Python standard library (`json`, `math`, `pathlib`), executes built-in tolerance self-checks, writes `results_bristol_flc_clean.json`, and overwrites this auto-generated `README.md` from the computed JSON values (manual edits to `README.md` will be lost on re-run). All episode numbers are **1-based** (`ep 1 .. ep N`).

## 1. Author Reference Settings (Agreement of Implementations)
- Base FLC resistance: `R_base = 7000 Pa*s/ml` (`D = 1.43e-4 ml/(s*Pa)`), offsets `delta = (0, -0.2, +0.2, 0)` for Fig. 4D/E and `(0, 0, 0, 0)` for Fig. 5G-I.
- Pump & valve conductances: single `4.0 ml/s` ideal flow source, open input valves modeled with conductance `1e6 ml/(s*Pa)`, closed output valves modeled with `1e-12 ml/(s*Pa)`, open `Out1 = 1/700 ml/(s*Pa)` (`700 Pa*s/ml`), open `Out2 = 1/7 ml/(s*Pa)` (`7 Pa*s/ml`).
- Episode-mean update & normalized deficit: conductances update once per episode from time-weighted mean flow `Q_bar`, and `Eq. (4)` uses normalized deficit `((Qc - Q_bar)+ / Qc)^2`.

| Simulation Run (1-based episodes) | P0 (Pa) | Ep 1 (Pa) | Ep 20 (Pa) | Min Error (Pa) | Ep 140 (Pa) | Ep 400 (Pa) | Long-Horizon & Reference Values |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| Irreversible, Fig. 4D/E (140 ep) | 1245.9 | 776.0 | 133.0 | 105.0 (ep 44) | **771.9** | - | Author's Python `772` (`3522 s`), MATLAB `787`, ours `771.9` (`3521.7 s`) |
| Reversible, Fig. 5G-I (1,400 ep; extended to 10,000) | 1274.1 | 764.5 | 46.3 | 6.1 (ep 259) | 23.2 | **27.1** (Py: `27`) | Python port: `101 Pa` (`ep 1,400`, `17,047 s`; ours `101.4 Pa`, `17047.1 s`), `107.6 Pa` (`ep 10,000`, vs `~108 Pa` in author's `reproduce.py --long`). Original MATLAB (`reproduce.py` header): `108 Pa` (`ep 1,000`), `127 Pa` (`ep 1,400`); author's reference time for the experimental protocol: `18,830 s`. |

*Simulation vs. physical experiment:* These runs reproduce the lumped resistive **simulation** (`772 Pa` in Python / `787 Pa` in MATLAB at ep 140), whereas in the physical PVA hardware experiment (`Fig. 4C`) total error reaches `1192 Pa` at episode 137.

*Channel resistances & conductances on the plateau:* At episode 1,400, resistances in `(in1->out1, in2->out1, in1->out2, in2->out2)` order are `(6996.7, 6958.5, 3409.3, 6999.2) Pa*s/ml` (matching the author's `reproduce.py` printout `(6997, 6958, 3409, 6999)`). Channel `in1->out2` is `2.0532x D_ini` at episode 1,400 and `2.0322x D_ini` at episode 10,000, while the other three channels relax near `D_ini` (at ep 1,400: `in2->out1 (D1)` = 1.0060x, `in2->out2 (D3)` = 1.0001x, `in1->out1 (D4)` = 1.0005x; at ep 10,000: `in2->out1 (D1)` = 1.0000x, `in2->out2 (D3)` = 1.0000x, `in1->out1 (D4)` = 1.0000x). In the irreversible Fig. 4D/E run, all four channels approach the `10x` ceiling by episode 140 (`8.9x-9.9x D_ini`), consistent with late forgetting.

## 2. Initial Blind Attempt (Prior Parameter Assumptions)
Before the reference parameters were shared, our initial blind run assumed `D_base = 1.9e-4 ml/(s*Pa)`, `2.5 mm x 20 cm` tubes on inputs and `Out2`, `q_pump = 3.40-4.08 ml/s`, continuous `dt = 0.2 s` integration, and unnormalized `(Qc - Q)+^2`. That gave `729-875 Pa` at ep 1, `136-162 Pa` min near ep 40-45, and `413-554 Pa` at ep 140. Under continuous unnormalized `Eq. (4)` at `q_pump = 3.40 ml/s`, error rose to `~713 Pa` at episode 400 (`~8,811 s` due to longer `~22 s` episodes under those assumptions), consistent with relaxation toward `D_ini` reported by the author for unnormalized continuous `Eq. (4)`.

## Acknowledgments
Many thanks to Rui Wu (University of Bristol) for sharing the reference simulation settings, Python port ([wurui1991/flc-reproduction-notes](https://github.com/wurui1991/flc-reproduction-notes)), and raw experimental data.
