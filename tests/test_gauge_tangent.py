import json
import unittest
from pathlib import Path

from mcgt_repro.build_gauge_tangent import gauge_tangent_columns
from mcgt_repro.build_internal_jacobian import internal_jacobian_rows
from mcgt_repro.exact_rank import apply_matrix, is_zero_columns, nullspace, rank_columns, rank_rows
from mcgt_repro.graph_fixture import EDGES, N_VERTICES, declared_cycle_rows


EXPECTED = json.loads(
    (Path(__file__).resolve().parents[1] / "expected" / "twelve_vertex_witness.json").read_text()
)


class GaugeTangentTest(unittest.TestCase):
    def test_rank_inclusion_and_quotient(self):
        rows, domain = internal_jacobian_rows(
            EDGES, declared_cycle_rows(), EXPECTED["r"], N_VERTICES
        )
        gauge = gauge_tangent_columns(EDGES, EXPECTED["r"], N_VERTICES)
        gauge_rank = rank_columns(gauge)
        self.assertEqual(gauge_rank, EXPECTED["rank_T_gauge"])
        self.assertTrue(is_zero_columns(apply_matrix(rows, gauge)))
        nullity = len(nullspace(rows, domain))
        self.assertEqual(nullity, domain - rank_rows(rows))
        self.assertEqual(nullity, EXPECTED["dim_ker_J_int"])
        self.assertEqual(nullity - gauge_rank, EXPECTED["dim_H_miss"])


if __name__ == "__main__":
    unittest.main()
