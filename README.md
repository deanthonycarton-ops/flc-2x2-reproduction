# Fluidic Learning Channel (FLC) - 2x2 Reproduction & Clog-Recovery Study

Independent numerical reproduction (with explicitly stated parameter assumptions) of the `2x2` Fluidic Learning Channel experiment and reversible simulation in:

> R. Wu, C. Lee, S. Terrile, and H. Hauser, *"Physical Learning in Soft Fluidic Channels through Experience,"* IEEE RoboSoft 2026.

## Quick Start

```bash
pip install -r requirements.txt
python reproduce_bristol_flc_clean.py
```

Outputs `results_bristol_flc_clean.json`.

## 1. Assumed Parameters (Not Explicitly Stated in the Paper)
- Initial base FLC conductance: `D_base = 1.9e-4 ml/(s*Pa)` (estimated from Fig. 2C at `t = 0`), with initial-conductance offsets `delta = (0, -0.2, +0.2, 0)` from Sec. IV-C.
- Pump flow rate: tested at `3.40 ml/s` and `4.08 ml/s` (single common pump header per Fig. 4A).
- Input and Out2 tubing: `2.5 mm ID x 20 cm` silicone tubing (`Out1` uses the published `1.8 mm ID x 20 cm` constrictor tube).
- Validation valve state: both `Out1` and `Out2` open during the `2 s` validation windows.

With these assumptions, total error drops from `729-875 Pa` at episode 0 to `152-184 Pa` by episode 20, reaches a minimum of `136-162 Pa` around episodes 39-44, and rises to `413-554 Pa` by episode 140.

## 2. Reversible `Eq. (4)` Sensitivity & Unknown 2-Channel Clog Benchmark (28 Pairs)
In a redundant 8-channel `2x2` topology at equal total conductance budget (`sum D = 0.003683 ml/(s*Pa)`), all **28 possible dual-channel clog pairs** (24 route-preserving pairs + 4 total route-severing pairs) were evaluated against an open-loop minimax static design optimized across failure states:

| Arm (`survivable 24 pairs` / `all 28 pairs`) | Clean Error (Pa) | Peak Spike at Clog (Pa) | Final Mean Error (Pa) | Final Worst-Pair Error (Pa) |
| :--- | :---: | :---: | :---: | :---: |
| `minimax_static_open_loop` | `185.2` / `195.8` | `868.8` / `1691.1` | `568.3` / `715.0` | `868.8` / `1691.1` |
| `taxis_rev_ideal_ode` | `82.0` / `82.0` | `2517.7` / `2605.9` | `120.7` / `266.6` | `515.8` / `1542.3` |
| `rev_eq4_paper` | `110.9` / `110.9` | `3145.6` / `3145.6` | `466.5` / `658.5` | `970.8` / `2660.5` |

### Caveats & Limitations
1. **Idealized ODE vs `Eq. (4)`:** Worst-case final error improves over `minimax_static_open_loop` only for the idealized linear-decay ODE (`taxis_rev_ideal_ode`, not yet demonstrated as a passive material) on the 24 route-preserving pairs (`515.8` vs `868.8 Pa`; nearly equal across all 28 pairs at `1542.3` vs `1691.1 Pa`). Under cumulative integration of `Eq. (4)` (`rev_eq4_paper`), worst-case error exceeds minimax static (`970.8` vs `868.8 Pa` on 24 pairs; `2660.5` vs `1691.1 Pa` on all 28), improving only on the 24-pair mean (`466.5` vs `568.3 Pa`).
2. **Transient failure spike:** Online retraining incurs a `2.9x-3.6x` higher immediate error spike right at channel loss (`2517.7-3145.6 Pa` vs `868.8 Pa` on survivable pairs).
3. **External feedback:** Online adaptation requires pressure sensing and teacher valve actuation; a closed-loop adjustable valve baseline is not modeled here.
