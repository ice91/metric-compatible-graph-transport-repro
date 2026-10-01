import json
import unittest
from pathlib import Path

from mcgt_repro.exact_rank import rank_rows
from mcgt_repro.graph_fixture import (
    boundary_matrix,
    cycle_space_basis,
    declared_cycle_rows,
    rows_in_kernel,
)


EXPECTED = json.loads(
    (Path(__file__).resolve().parents[1] / "expected" / "twelve_vertex_witness.json").read_text()
)


class CycleSpaceTest(unittest.TestCase):
    def test_declared_and_basis_lie_in_kernel(self):
        boundary = boundary_matrix()
        self.assertTrue(rows_in_kernel(declared_cycle_rows(), boundary))
        self.assertTrue(rows_in_kernel(cycle_space_basis(), boundary))

    def test_cycle_space_rank(self):
        self.assertEqual(rank_rows(cycle_space_basis()), EXPECTED["q_max"])

    def test_declared_row_rank(self):
        self.assertEqual(rank_rows(declared_cycle_rows()), EXPECTED["q"])


if __name__ == "__main__":
    unittest.main()
