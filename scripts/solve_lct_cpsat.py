import json
import os
import sys
from typing import Dict, Tuple

from ortools.sat.python import cp_model


def load(path: str):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save(path: str, data):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def and_var(model: cp_model.CpModel, a, b, name: str):
    m = model.NewBoolVar(name)
    model.Add(m <= a)
    model.Add(m <= b)
    model.Add(m >= a + b - 1)
    return m


def build_and_solve(inp):
    model = cp_model.CpModel()

    s = inp["s"]
    h = inp["h"]
    p = inp["p"]
    c = inp["c"]
    kmax = inp["k"]
    tau = inp["tau"]
    model_name = inp["model"]

    init_x = inp["init_x"]
    group_count = inp["group_count"]
    table_b = inp["table_b"]

    h0_0_s = inp["h0_0_s"]
    h0_rs = inp["h0_rs"]
    h1_sr = inp["h1_sr"]
    v_load = inp["v_load"]
    v_unload = inp["v_unload"]
    warm_start_y_ones = inp.get("warm_start_y_ones", [])
    warm_start_z_ones = inp.get("warm_start_z_ones", [])

    S = range(1, s + 1)
    H = range(1, h + 1)
    K = range(1, kmax + 1)
    K0 = range(0, kmax + 1)

    x = {}
    d = {}
    y = {}
    z = {}
    g = {}
    u = {}

    for ss in S:
        for tt in H:
            for kk in K0:
                x[(ss, tt, kk)] = model.NewIntVar(0, p, f"x_{ss}_{tt}_{kk}")
                d[(ss, tt, kk)] = model.NewBoolVar(f"d_{ss}_{tt}_{kk}")
            for kk in K:
                y[(ss, tt, kk)] = model.NewBoolVar(f"y_{ss}_{tt}_{kk}")
                z[(ss, tt, kk)] = model.NewBoolVar(f"z_{ss}_{tt}_{kk}")

    for ss in S:
        for rr in S:
            if ss == rr:
                continue
            for kk in K:
                g[(ss, rr, kk)] = model.NewBoolVar(f"g_{ss}_{rr}_{kk}")
            for kk in range(2, kmax + 1):
                u[(ss, rr, kk)] = model.NewBoolVar(f"u_{ss}_{rr}_{kk}")

    q_slot = {}
    q_group = {}
    if model_name == "lct1":
        for ss in S:
            for tt in H:
                q_slot[(ss, tt)] = model.NewBoolVar(f"q_{ss}_{tt}")
    else:
        for pp in range(1, p + 1):
            q_group[pp] = model.NewBoolVar(f"qp_{pp}")

    eq_cache: Dict[Tuple[Tuple[int, int, int], int], cp_model.IntVar] = {}
    gt_cache: Dict[Tuple[Tuple[int, int, int], Tuple[int, int, int]], cp_model.IntVar] = {}

    def eq_value(var_key, value: int):
        key = (var_key, value)
        if key in eq_cache:
            return eq_cache[key]
        v = model.NewBoolVar(f"eq_{var_key[0]}_{var_key[1]}_{var_key[2]}_{value}")
        vv = x[var_key]
        model.Add(vv == value).OnlyEnforceIf(v)
        model.Add(vv != value).OnlyEnforceIf(v.Not())
        eq_cache[key] = v
        return v

    def gt_var(a_key, b_key):
        key = (a_key, b_key)
        if key in gt_cache:
            return gt_cache[key]
        v = model.NewBoolVar(
            f"gt_{a_key[0]}_{a_key[1]}_{a_key[2]}_{b_key[0]}_{b_key[1]}_{b_key[2]}"
        )
        av = x[a_key]
        bv = x[b_key]
        model.Add(av > bv).OnlyEnforceIf(v)
        model.Add(av <= bv).OnlyEnforceIf(v.Not())
        gt_cache[key] = v
        return v

    for ss in S:
        for tt in H:
            x0 = int(init_x[ss][tt])
            model.Add(x[(ss, tt, 0)] == x0)
            model.Add(d[(ss, tt, 0)] == (1 if x0 > 0 else 0))

    for kk in K:
        for pp in range(0, p + 1):
            bvars = []
            for ss in S:
                for tt in H:
                    b = eq_value((ss, tt, kk), pp)
                    bvars.append(b)
            model.Add(sum(bvars) == int(group_count[pp]))

    for kk in K:
        for ss in S:
            for tt in H:
                xv = x[(ss, tt, kk)]
                dv = d[(ss, tt, kk)]
                model.Add(xv >= 1).OnlyEnforceIf(dv)
                model.Add(xv == 0).OnlyEnforceIf(dv.Not())

                b_eq_x = model.NewBoolVar(f"beqx_{ss}_{tt}_{kk}")
                model.Add(x[(ss, tt, kk)] == x[(ss, tt, kk - 1)]).OnlyEnforceIf(b_eq_x)
                model.Add(x[(ss, tt, kk)] != x[(ss, tt, kk - 1)]).OnlyEnforceIf(b_eq_x.Not())

                b_eq_d = model.NewBoolVar(f"beqd_{ss}_{tt}_{kk}")
                model.Add(d[(ss, tt, kk)] == d[(ss, tt, kk - 1)]).OnlyEnforceIf(b_eq_d)
                model.Add(d[(ss, tt, kk)] != d[(ss, tt, kk - 1)]).OnlyEnforceIf(b_eq_d.Not())

                model.Add(b_eq_x == b_eq_d)

    for kk in range(1, kmax):
        for ss in S:
            for tt in H:
                model.Add(y[(ss, tt, kk)] <= d[(ss, tt, kk + 1)])

    for kk in K:
        model.Add(sum(y[(ss, tt, kk)] for ss in S for tt in H) <= 1)
        model.Add(sum(z[(ss, tt, kk)] for ss in S for tt in H) <= 1)

    for kk in range(1, kmax):
        model.Add(
            sum(z[(ss, tt, kk + 1)] for ss in S for tt in H)
            <= sum(z[(ss, tt, kk)] for ss in S for tt in H)
        )

    if h >= 2:
        for kk in K:
            for ss in S:
                for tt in range(1, h):
                    vars8 = [
                        d[(ss, tt, kk - 1)],
                        d[(ss, tt + 1, kk - 1)],
                        d[(ss, tt, kk)],
                        d[(ss, tt + 1, kk)],
                        z[(ss, tt, kk)],
                        z[(ss, tt + 1, kk)],
                        y[(ss, tt, kk)],
                        y[(ss, tt + 1, kk)],
                    ]
                    model.AddAllowedAssignments(vars8, table_b)

    for kk in K:
        for ss in S:
            for rr in S:
                if ss == rr:
                    continue
                expr = sum(z[(ss, tt, kk)] for tt in H) + sum(y[(rr, tt, kk)] for tt in H)
                gv = g[(ss, rr, kk)]
                model.Add(expr == 2).OnlyEnforceIf(gv)
                model.Add(expr <= 1).OnlyEnforceIf(gv.Not())

                if kk >= 2:
                    expr_u = sum(y[(ss, tt, kk - 1)] for tt in H) + sum(
                        z[(rr, tt, kk)] for tt in H
                    )
                    uv = u[(ss, rr, kk)]
                    model.Add(expr_u == 2).OnlyEnforceIf(uv)
                    model.Add(expr_u <= 1).OnlyEnforceIf(uv.Not())

    SCALE = 1000
    terms = []

    for ss in S:
        for rr in S:
            if ss == rr:
                continue
            terms.append((int(round(h0_0_s[ss] * SCALE)), g[(ss, rr, 1)]))

    for rr in S:
        for ss in S:
            if rr == ss:
                continue
            for kk in range(2, kmax + 1):
                terms.append((int(round(h0_rs[rr][ss] * SCALE)), u[(rr, ss, kk)]))

    for ss in S:
        for rr in S:
            if ss == rr:
                continue
            for kk in K:
                terms.append((int(round(h1_sr[ss][rr] * SCALE)), g[(ss, rr, kk)]))

    for ss in S:
        for tt in H:
            for kk in K:
                terms.append((int(round(v_load[tt] * SCALE)), z[(ss, tt, kk)]))
                terms.append((int(round(v_unload[tt] * SCALE)), y[(ss, tt, kk)]))

    tau_i = int(round(tau * SCALE))
    model.Add(sum(coeff * var for coeff, var in terms) <= tau_i)

    kf = kmax

    if model_name == "lct1":
        if h >= 2:
            for ss in S:
                for tt in range(1, h):
                    bvars = [gt_var((ss, jj, kf), (ss, tt, kf)) for jj in range(tt + 1, h + 1)]
                    model.Add(sum(bvars) <= h * q_slot[(ss, tt)])

        for ss in S:
            for tt in H:
                lhs = []
                for ii in S:
                    for jj in H:
                        gt = gt_var((ss, tt, kf), (ii, jj, kf))
                        m = and_var(model, gt, q_slot[(ii, jj)], f"mul_{ss}_{tt}_{ii}_{jj}")
                        lhs.append(m)
                model.Add(sum(lhs) <= c * q_slot[(ss, tt)])

        model.Minimize(sum(q_slot[(ss, tt)] for ss in S for tt in H))
    else:
        for pp in range(1, p + 1):
            lhs = []
            for ss in S:
                for tt in H:
                    if tt >= h:
                        continue
                    eqp = eq_value((ss, tt, kf), pp)
                    for jj in range(tt + 1, h + 1):
                        gt = gt_var((ss, jj, kf), (ss, tt, kf))
                        m = and_var(model, eqp, gt, f"mul2_{pp}_{ss}_{tt}_{jj}")
                        lhs.append(m)
            model.Add(sum(lhs) <= c * q_group[pp])

        for pp in range(1, p):
            model.Add(q_group[pp] <= q_group[pp + 1])

        model.Minimize(sum(q_group[pp] for pp in range(1, p + 1)))

    for (kk, ss, tt) in warm_start_y_ones:
        if 1 <= kk <= kmax and 1 <= ss <= s and 1 <= tt <= h:
            model.AddHint(y[(ss, tt, kk)], 1)
    for (kk, ss, tt) in warm_start_z_ones:
        if 1 <= kk <= kmax and 1 <= ss <= s and 1 <= tt <= h:
            model.AddHint(z[(ss, tt, kk)], 1)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(inp["time_limit_sec"])

    workers = 1
    for env_name in ("CP_SAT_NUM_WORKERS", "CP_OPTIMIZER_WORKERS", "SLURM_CPUS_PER_TASK"):
        raw = os.environ.get(env_name, "").strip()
        if not raw:
            continue
        try:
            value = int(raw)
        except ValueError:
            continue
        if value >= 1:
            workers = value
            break
    solver.parameters.num_search_workers = workers
    solver.parameters.max_memory_in_mb = 16384

    status = solver.Solve(model)

    status_name = {
        cp_model.OPTIMAL: "OPTIMAL",
        cp_model.FEASIBLE: "FEASIBLE",
        cp_model.INFEASIBLE: "INFEASIBLE",
        cp_model.UNKNOWN: "UNKNOWN",
        cp_model.MODEL_INVALID: "MODEL_INVALID",
    }.get(status, "UNKNOWN")

    has_solution = status in (cp_model.OPTIMAL, cp_model.FEASIBLE)

    y_ones = []
    z_ones = []

    if has_solution:
        for kk in K:
            for ss in S:
                for tt in H:
                    if solver.Value(y[(ss, tt, kk)]) == 1:
                        y_ones.append([kk, ss, tt])
                    if solver.Value(z[(ss, tt, kk)]) == 1:
                        z_ones.append([kk, ss, tt])

    obj = int(round(solver.ObjectiveValue())) if has_solution else None

    return {
        "status": status_name,
        "optimal": status == cp_model.OPTIMAL,
        "time_seconds": solver.WallTime(),
        "objective": obj,
        "y_ones": y_ones,
        "z_ones": z_ones,
    }


def main():
    if len(sys.argv) != 3:
        print("usage: solve_lct_cpsat.py <input.json> <output.json>", file=sys.stderr)
        sys.exit(2)

    in_path = sys.argv[1]
    out_path = sys.argv[2]

    inp = load(in_path)
    out = build_and_solve(inp)
    save(out_path, out)


if __name__ == "__main__":
    main()
