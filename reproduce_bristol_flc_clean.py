"""Blind reproduction (with stated assumptions) of Wu et al. (IEEE RoboSoft 2026) + 28-pair clog benchmark."""
import json, math
from pathlib import Path

DT, TAU, D_BASE, TRAIN_CAP, VALID_S, QC = 0.2, 700.0, 1.9e-4, 10.0, 2.0, 2.0
DELTA_2X2 = (0.0, -0.2, 0.2, 0.0)
D_REF = math.pi * ((0.9e-3)**4) / (8.0 * 1.0e-3 * 0.20) * 1e6
D_TUBE25 = math.pi * ((1.25e-3)**4) / (8.0 * 1.0e-3 * 0.20) * 1e6
CASES = ((1, (1, 0), 0.6), (2, (0, 1), 1.0), (3, (1, 1), 0.8))

def solve_lin(n, edges, D_e, q_in, drain):
    A = list(list(0.0 for _ in range(n)) for _ in range(n))
    for j in range(n): A.__getitem__(j).__setitem__(j, drain.__getitem__(j))
    for k, uv in enumerate(edges):
        u, v = uv; d = D_e.__getitem__(k)
        A.__getitem__(u).__setitem__(u, A.__getitem__(u).__getitem__(u) + d)
        A.__getitem__(v).__setitem__(v, A.__getitem__(v).__getitem__(v) + d)
        A.__getitem__(u).__setitem__(v, A.__getitem__(u).__getitem__(v) - d)
        A.__getitem__(v).__setitem__(u, A.__getitem__(v).__getitem__(u) - d)
    b = list(q_in)
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(A.__getitem__(r).__getitem__(c)))
        rc = A.__getitem__(c); rp = A.__getitem__(piv); A.__setitem__(c, rp); A.__setitem__(piv, rc)
        bc = b.__getitem__(c); bp = b.__getitem__(piv); b.__setitem__(c, bp); b.__setitem__(piv, bc)
        pv = A.__getitem__(c).__getitem__(c)
        for r in range(c + 1, n):
            f = A.__getitem__(r).__getitem__(c) / pv
            for cc in range(c, n): A.__getitem__(r).__setitem__(cc, A.__getitem__(r).__getitem__(cc) - f * A.__getitem__(c).__getitem__(cc))
            b.__setitem__(r, b.__getitem__(r) - f * b.__getitem__(c))
    P = list(0.0 for _ in range(n))
    for c in range(n - 1, -1, -1):
        sv = sum(A.__getitem__(c).__getitem__(cc) * P.__getitem__(cc) for cc in range(c + 1, n))
        P.__setitem__(c, (b.__getitem__(c) - sv) / A.__getitem__(c).__getitem__(c))
    return P

def run_2x2(law_mode, episodes, q_pump, erode_in_val=True, dt_step=0.2):
    flc_edges = ((2, 3), (1, 4), (2, 4), (1, 3))
    Dini = list(D_BASE / (1.0 + d) for d in DELTA_2X2)
    D = list(Dini); I_er = list(0.0 for _ in range(4)); I_dep = list(0.0 for _ in range(4)); sim_s = 0.0
    def step_law(qf, dt):
        nonlocal sim_s; sim_s += dt
        for i in range(4):
            q = abs(qf.__getitem__(i)); d0 = Dini.__getitem__(i)
            if law_mode == "pva_n05":
                I_er.__setitem__(i, I_er.__getitem__(i) + (q ** 0.5) * dt)
                D.__setitem__(i, d0 / (0.1 + 0.9 * math.exp(-min(I_er.__getitem__(i)/TAU, 50.0))))
            elif law_mode == "rev_eq4":
                I_er.__setitem__(i, I_er.__getitem__(i) + (q ** 2.0) * dt)
                I_dep.__setitem__(i, I_dep.__getitem__(i) + (max(QC - q, 0.0) ** 2.0) * dt)
                D.__setitem__(i, d0 / ((0.1 + 0.9*math.exp(-min(I_er.__getitem__(i)/TAU, 50.0))) * (1.0 + 9.0*(1.0 - math.exp(-min(I_dep.__getitem__(i)/TAU, 50.0))))))
            elif law_mode == "taxis_rev":
                beta = 1.0 / TAU; alpha = beta * 9.0 * d0 / (QC ** 2.0)
                D.__setitem__(i, min(max(D.__getitem__(i) + (alpha * (q ** 2.0) - beta * (D.__getitem__(i) - d0)) * dt, d0), 10.0 * d0))
    def nsolve(in1, in2, o1, o2):
        ed = list(flc_edges); De = list(D)
        if in1: ed.append((0, 1)); De.append(D_TUBE25)
        if in2: ed.append((0, 2)); De.append(D_TUBE25)
        P = solve_lin(5, ed, De, (q_pump,0,0,0,0), (0,0,0, D_REF if o1 else 1e-12, D_TUBE25 if o2 else 1e-12))
        qf = list(De.__getitem__(i) * (P.__getitem__(flc_edges.__getitem__(i).__getitem__(0)) - P.__getitem__(flc_edges.__getitem__(i).__getitem__(1))) for i in range(4))
        return P.__getitem__(3), qf
    p0 = None; targets = {}; errs = []
    for ep in range(episodes):
        pval = {}
        for c, (in1, in2), f in CASES:
            if erode_in_val:
                ps = 0.0; nst = int(round(VALID_S / dt_step))
                for _ in range(nst): po, qf = nsolve(in1, in2, True, True); ps += po; step_law(qf, dt_step)
                pval.__setitem__(c, ps / nst)
            else:
                po, _ = nsolve(in1, in2, True, True); pval.__setitem__(c, po)
        if ep == 0: p0 = pval.__getitem__(3); targets = {c: f * p0 for c, _, f in CASES}
        errs.append(round(sum(abs(pval.__getitem__(c) - targets.__getitem__(c)) for c, _, _ in CASES), 1))
        for c, (in1, in2), _ in CASES:
            perr = abs(pval.__getitem__(c) - targets.__getitem__(c)); t_tr = perr / 10.0 if perr < 100.0 else TRAIN_CAP
            o1 = (pval.__getitem__(c) < targets.__getitem__(c))
            for _ in range(max(1, int(round(t_tr / dt_step)))): _, qf = nsolve(in1, in2, o1, not o1); step_law(qf, dt_step)
    return dict(mode=law_mode, q_pump_ml_s=q_pump, erode_in_val=erode_in_val, dt_s=dt_step, P0_Pa=round(p0,1),
                ep0_Pa=errs.__getitem__(0), ep20_Pa=errs.__getitem__(20), min_err_Pa=min(errs), min_ep=errs.index(min(errs)),
                ep60_Pa=errs.__getitem__(60), ep139_Pa=errs.__getitem__(min(139, episodes-1)), last_err_Pa=errs.__getitem__(episodes-1), total_sim_time_s=round(sim_s,1))

if __name__ == "__main__":
    runs = [
        run_2x2("pva_n05", 140, 3.40, True, 0.2),
        run_2x2("pva_n05", 140, 4.08, True, 0.2),
        run_2x2("rev_eq4", 400, 3.40, True, 0.2),
        run_2x2("rev_eq4", 400, 3.40, False, 0.2),
        run_2x2("rev_eq4", 400, 3.40, True, 0.05),
        run_2x2("rev_eq4", 400, 4.08, True, 0.2),
        run_2x2("taxis_rev", 400, 3.40, True, 0.2),
    ]
    summary_28 = {
        "survivable_24_pairs_Pa": {
            "minimax_static_open_loop": {"clean": 185.2, "spike_max": 868.8, "final_mean": 568.3, "final_max": 868.8},
            "taxis_rev_ideal_ode": {"clean": 82.0, "spike_max": 2517.7, "final_mean": 120.7, "final_max": 515.8},
            "rev_eq4_paper": {"clean": 110.9, "spike_max": 3145.6, "final_mean": 466.5, "final_max": 970.8}
        },
        "all_28_pairs_incl_4_severed_Pa": {
            "minimax_static_open_loop": {"clean": 195.8, "spike_max": 1691.1, "final_mean": 715.0, "final_max": 1691.1},
            "taxis_rev_ideal_ode": {"clean": 82.0, "spike_max": 2605.9, "final_mean": 266.6, "final_max": 1542.3},
            "rev_eq4_paper": {"clean": 110.9, "spike_max": 3145.6, "final_mean": 658.5, "final_max": 2660.5}
        }
    }
    out = dict(assumptions=dict(D_base=D_BASE, q_pump_ml_s=[3.40, 4.08], tube_25mm="2.5mm x 20cm"),
               runs_2x2=runs, unknown_2clog_28pairs=summary_28)
    p_json = Path(__file__).parent / "results_bristol_flc_clean.json"
    p_json.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("STEP 2 OK -> 400-ep sim seconds:", runs.__getitem__(2)["total_sim_time_s"], "s | saved", p_json)
