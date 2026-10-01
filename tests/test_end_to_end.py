import json
import unittest
from pathlib import Path

from mcgt_repro.verify_grid import EXPECTED_PATH, _matches, compute


class EndToEndTest(unittest.TestCase):
    def test_computed_witness_matches_expected_file(self):
        expected = json.loads(EXPECTED_PATH.read_text())
        observed = compute()
        self.assertTrue(observed["gauge_image_is_zero"])
        self.assertTrue(observed["declared_rows_lie_in_kernel"])
        self.assertTrue(observed["tree_is_spanning_forest"])
        self.assertEqual(
            observed["dim_ker_J_int"] - observed["rank_T_gauge"],
            expected["dim_H_miss"],
        )
        self.assertTrue(_matches(observed, expected))


if __name__ == "__main__":
    unittest.main()
