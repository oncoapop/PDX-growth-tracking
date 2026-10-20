import unittest

import pandas as pd

from plot_tumour_growth_by_cage import aggregate_for_plot, matplotlib, normalize_cage


class CageGroupingTests(unittest.TestCase):
    def test_plotting_uses_headless_backend(self) -> None:
        self.assertEqual(matplotlib.get_backend().lower(), "agg")

    def test_parenthesised_cage_numbers_remain_distinct(self) -> None:
        self.assertEqual(normalize_cage("251220(1)"), "251220(1)")
        self.assertEqual(normalize_cage("251220(2)"), "251220(2)")
        self.assertEqual(normalize_cage("251220(3)"), "251220(3)")

    def test_parenthesised_cage_label_is_preserved(self) -> None:
        self.assertEqual(normalize_cage("251220(Cage 3)"), "251220(Cage 3)")

    def test_trailing_mouse_suffix_stays_in_parent_cage(self) -> None:
        self.assertEqual(normalize_cage("260409(2)-1"), "260409(2)")

    def test_plot_groups_do_not_merge_parenthesised_cages(self) -> None:
        rows = []
        for mouse_id in ("251220(1)", "251220(3)"):
            rows.append(
                {
                    "date": pd.Timestamp("2026-06-22"),
                    "cage": normalize_cage(mouse_id),
                    "sheet": "Xenograft",
                    "tumour_model": "SA501 Breast",
                    "mouse_id": mouse_id,
                    "tumour_site": "B",
                    "tumour_volume_mm3": 100.0,
                }
            )

        result = aggregate_for_plot(pd.DataFrame(rows), "sum")

        self.assertEqual(
            set(result["plot_group"]),
            {
                "251220(1) SA501 Breast (Xenograft)",
                "251220(3) SA501 Breast (Xenograft)",
            },
        )


if __name__ == "__main__":
    unittest.main()
