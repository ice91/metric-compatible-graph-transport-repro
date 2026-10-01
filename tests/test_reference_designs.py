import json
import unittest
from pathlib import Path

from mcgt_repro.build_gauge_tangent import gauge_tangent_columns
from mcgt_repro.build_internal_jacobian import internal_jacobian_rows
from mcgt_repro.build_reference_designs import design_a_rows, design_b_rows
from mcgt_repro.exact_rank import apply_matrix, nullspace, rank_columns
from mcgt_repro.graph_fixture import EDGES, N_VERTICES, declared_cycle_rows, tree_indices


EXPECTED = json.loads(
    (Path(__file__).resolve().parents[1] / "expected" / "twelve_vertex_witness.json").read_text()
)


class ReferenceDesignTest(unittest.TestCase):
    def test_restricted_ranks(self):
        sector = EXPECTED["r"]
        rows, domain = internal_jacobian_rows(EDGES, declared_cycle_rows(), sector, N_VERTICES)
        gauge = gauge_tangent_columns(EDGES, sector, N_VERTICES)
        kernel = nullspace(rows, domain)
        tree = tree_indices()
        rows_a = design_a_rows(N_VERTICES, len(EDGES), tree, sector)
        rows_b = design_b_rows(N_VERTICES, len(EDGES), tree, sector, root=0)
        self.assertEqual(
            rank_columns(apply_matrix(rows_a, gauge)), EXPECTED["rank_A_A_on_T_gauge"]
        )
        self.assertEqual(
            rank_columns(apply_matrix(rows_b, gauge)), EXPECTED["rank_A_B_on_T_gauge"]
        )
        self.assertEqual(
            rank_columns(apply_matrix(rows_a, kernel)), EXPECTED["rank_A_A_on_ker"]
        )
        self.assertEqual(
            rank_columns(apply_matrix(rows_b, kernel)), EXPECTED["rank_A_B_on_ker"]
        )


if __name__ == "__main__":
    unittest.main()
