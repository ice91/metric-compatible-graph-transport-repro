import json
import unittest
from pathlib import Path

from mcgt_repro.build_internal_jacobian import internal_jacobian_rows
from mcgt_repro.exact_rank import nullspace, rank_rows
from mcgt_repro.graph_fixture import EDGES, N_VERTICES, declared_cycle_rows
from mcgt_repro.real_matrix_basis import gl_basis, gl_coords


EXPECTED = json.loads(
    (Path(__file__).resolve().parents[1] / "expected" / "twelve_vertex_witness.json").read_text()
)


class InternalJacobianTest(unittest.TestCase):
    def setUp(self):
        self.rows, self.domain = internal_jacobian_rows(
            EDGES, declared_cycle_rows(), EXPECTED["r"], N_VERTICES
        )

    def test_basis_is_dual(self):
        for index, matrix in enumerate(gl_basis(EXPECTED["r"])):
            coordinates = gl_coords(matrix)
            self.assertEqual(coordinates[index], 1)
            self.assertTrue(
                all(coordinates[other] == 0 for other in range(len(coordinates)) if other != index)
            )

    def test_shape_and_rank(self):
        self.assertEqual(len(self.rows), EXPECTED["jacobian_rows"])
        self.assertEqual(self.domain, EXPECTED["jacobian_columns"])
        self.assertTrue(all(len(row) == self.domain for row in self.rows))
        self.assertEqual(rank_rows(self.rows), EXPECTED["rank_J_int"])

    def test_nullity(self):
        kernel = nullspace(self.rows, self.domain)
        self.assertEqual(len(kernel), EXPECTED["dim_ker_J_int"])
        self.assertEqual(self.domain - rank_rows(self.rows), EXPECTED["dim_ker_J_int"])


if __name__ == "__main__":
    unittest.main()
