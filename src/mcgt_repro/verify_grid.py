"""Compute the twelve-vertex witness and compare it with the expected file."""

import json
import sys
from pathlib import Path

from mcgt_repro.build_gauge_tangent import gauge_tangent_columns
from mcgt_repro.build_internal_jacobian import internal_jacobian_rows
from mcgt_repro.build_reference_designs import design_a_rows, design_b_rows
from mcgt_repro.exact_rank import (
    apply_matrix,
    is_zero_columns,
    nullspace,
    rank_columns,
    rank_rows,
)
from mcgt_repro.graph_fixture import (
    EDGES,
    N_VERTICES,
    cycle_counts,
    declared_cycle_rows,
    tree_indices,
)


SECTOR = 2
ROOT = Path(__file__).resolve().parents[2]
EXPECTED_PATH = ROOT / "expected" / "twelve_vertex_witness.json"


def compute():
    counts = cycle_counts()
    declared = declared_cycle_rows()
    rows, domain = internal_jacobian_rows(EDGES, declared, SECTOR, N_VERTICES)
    jacobian_rank = rank_rows(rows)
    kernel = nullspace(rows, domain)
    gauge = gauge_tangent_columns(EDGES, SECTOR, N_VERTICES)
    gauge_rank = rank_columns(gauge)
    image_is_zero = is_zero_columns(apply_matrix(rows, gauge))
    tree = tree_indices()
    rows_a = design_a_rows(N_VERTICES, len(EDGES), tree, SECTOR)
    rows_b = design_b_rows(N_VERTICES, len(EDGES), tree, SECTOR, root=0)
    nullity = domain - jacobian_rank
    missing = nullity - gauge_rank if image_is_zero else None
    return {
        "vertices": N_VERTICES,
        "free_edges": len(EDGES),
        "components": counts["component_count"],
        "r": SECTOR,
        "q_max": counts["q_max"],
        "q": counts["q"],
        "declared_rows_lie_in_kernel": counts["declared_rows_lie_in_kernel"],
        "cycle_basis_lies_in_kernel": counts["cycle_basis_lies_in_kernel"],
        "tree_is_spanning_forest": counts["tree_is_spanning_forest"],
        "jacobian_rows": len(rows),
        "jacobian_columns": domain,
        "rank_J_int": jacobian_rank,
        "dim_ker_J_int": len(kernel),
        "rank_T_gauge": gauge_rank,
        "gauge_image_is_zero": image_is_zero,
        "dim_H_miss": missing,
        "rank_A_A_on_T_gauge": rank_columns(apply_matrix(rows_a, gauge)),
        "rank_A_B_on_T_gauge": rank_columns(apply_matrix(rows_b, gauge)),
        "rank_A_A_on_ker": rank_columns(apply_matrix(rows_a, kernel)),
        "rank_A_B_on_ker": rank_columns(apply_matrix(rows_b, kernel)),
        "arithmetic": "exact rational",
        "floating_point_tolerance": None,
    }


def _matches(observed, expected):
    if not observed["gauge_image_is_zero"]:
        return False
    if observed["dim_ker_J_int"] != observed["jacobian_columns"] - observed["rank_J_int"]:
        return False
    keys = (
        "vertices",
        "free_edges",
        "components",
        "r",
        "q_max",
        "q",
        "jacobian_rows",
        "jacobian_columns",
        "rank_J_int",
        "dim_ker_J_int",
        "rank_T_gauge",
        "dim_H_miss",
        "rank_A_A_on_T_gauge",
        "rank_A_B_on_T_gauge",
        "rank_A_A_on_ker",
        "rank_A_B_on_ker",
    )
    return all(observed[key] == expected[key] for key in keys)


def format_report(observed, passed):
    zero = "exact zero" if observed["gauge_image_is_zero"] else "not zero"
    result = "PASS" if passed else "FAIL"
    return "\n".join(
        [
            "Twelve-vertex witness",
            "---------------------",
            "",
            "|V|                     = %s" % observed["vertices"],
            "|E_free|                = %s" % observed["free_edges"],
            "r                       = %s" % observed["r"],
            "q_max                   = %s" % observed["q_max"],
            "q                       = %s" % observed["q"],
            "",
            "rank(J_int)             = %s" % observed["rank_J_int"],
            "dim ker(J_int)          = %s" % observed["dim_ker_J_int"],
            "rank(T_gauge)           = %s" % observed["rank_T_gauge"],
            "dim H_miss              = %s" % observed["dim_H_miss"],
            "",
            "rank(A_A | T_gauge)     = %s" % observed["rank_A_A_on_T_gauge"],
            "rank(A_B | T_gauge)     = %s" % observed["rank_A_B_on_T_gauge"],
            "rank(A_A | ker J_int)   = %s" % observed["rank_A_A_on_ker"],
            "rank(A_B | ker J_int)   = %s" % observed["rank_A_B_on_ker"],
            "",
            "J_int(T_gauge)          = %s" % zero,
            "",
            "arithmetic              = exact rational",
            "floating-point tolerance = none",
            "",
            "RESULT                  = %s" % result,
            "",
        ]
    )


def main():
    expected = json.loads(EXPECTED_PATH.read_text())
    observed = compute()
    passed = _matches(observed, expected)
    report = format_report(observed, passed)
    sys.stdout.write(report)
    payload = dict(observed)
    payload["result"] = "PASS" if passed else "FAIL"
    (Path.cwd() / "reproduced_results.json").write_text(
        json.dumps(payload, indent=2) + "\n"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
