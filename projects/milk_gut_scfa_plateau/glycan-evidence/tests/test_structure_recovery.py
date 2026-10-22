import copy
import hashlib
import tempfile
import unittest
from pathlib import Path

from glycan_evidence.structure_recovery import DEFAULT_CONFIG, recover_structure, validate_file_binding


def fixture():
    binding = {"path": "bound.tsv", "bytes": 12, "sha256": "a" * 64}
    panel = [{"candidate_id": "ref", "inner_length": 824, "class_id": "short_824",
              "full_sequence_sha256": "b" * 64, "inner_sequence_sha256": "c" * 64}]
    candidate = {"case_id": "sample", "source_kind": "own_library_ordinary_assembly",
        "structure_id": "s1", "contig": "contig", "contig_length": 5000,
        "contig_sequence_sha256": "d" * 64, "interval_start": 100,
        "interval_end": 3000, "left_boundary": 1000, "right_boundary": 1824,
        "orientation": "+", "inner_length": 824, "anchor_screen_pass": True,
        "inner_edges_complete": True, "outer_anchors_complete": True,
        "interval_sequence_sha256": "b" * 64, "inner_sequence_sha256": "c" * 64,
        "coverage": {"unit": "primary_mate_alignment_read_depth", "filter": "quality_filtered",
                     "segment": "full_anchored_interval", "start": 100, "end": 3000,
                     "breadth_ge_1": 1, "min_read_depth": 5},
        "boundaries": {side: {"position": pos, "evidence": "direct_read_spanning",
                               "qnames": 5, "distinct_pair_starts": 3}
                       for side, pos in [("left", 1000), ("right", 1824)]},
        "source_bindings": {key: copy.deepcopy(binding) for key in ["receipt", "structures", "coverage", "boundaries"]},
        "base_concordance": None}
    return panel, candidate


class RecoveryTests(unittest.TestCase):
    def run_case(self, candidate=None, refs=None):
        panel, default = fixture()
        return recover_structure("sample", refs or [], panel, [default if candidate is None else candidate])

    def test_primary_and_strict_and_fair_comparator(self):
        result = self.run_case(refs=["short_824"])
        self.assertEqual(result["assembly_supported_classes"], ["short_824"])
        self.assertEqual(result["assembly_supported_classes"], result["assembly_only_supported_classes"])
        self.assertEqual(result["strict_assembly_supported_classes"], ["short_824"])
        self.assertFalse(result["discordance"])

    def test_single_fragment_primary_only(self):
        _, row = fixture()
        row["coverage"]["min_read_depth"] = 1
        for bound in row["boundaries"].values(): bound.update(qnames=1, distinct_pair_starts=1)
        result = self.run_case(row)["candidates"][0]
        self.assertTrue(result["assembly_supported"])
        self.assertFalse(result["strict_supported"])

    def test_reference_only_cannot_create_assembly(self):
        panel, _ = fixture()
        result = recover_structure("sample", ["short_824"], panel, [])
        self.assertEqual(result["assembly_supported_classes"], [])
        self.assertEqual(result["reference_only_classes"], ["short_824"])

    def test_partial_inner_is_not_exact_sequence_or_length(self):
        _, row = fixture()
        row.update(inner_edges_complete=False, outer_anchors_complete=False,
                   inner_length=None, inner_sequence_sha256=None)
        candidate = self.run_case(row)["candidates"][0]
        self.assertFalse(candidate["assembly_supported"])
        self.assertIsNone(candidate["length_class"])
        self.assertIsNone(candidate["new_sequence_candidate"])

    def test_partial_outer_can_support_inner_and_has_unknown_full_match(self):
        _, row = fixture(); row["outer_anchors_complete"] = False
        result = self.run_case(row)["candidates"][0]
        self.assertTrue(result["assembly_supported"])
        self.assertIsNone(result["exact_full_sequence_in_panel"])

    def test_novel_length_and_sequence_separate(self):
        _, row = fixture(); row.update(inner_length=825, right_boundary=1825,
                                       inner_sequence_sha256="e"*64, interval_sequence_sha256="f"*64)
        row["boundaries"]["right"]["position"] = 1825
        result = self.run_case(row)["candidates"][0]
        self.assertTrue(result["outside_panel_length_class"])
        self.assertTrue(result["new_sequence_candidate"])
        self.assertEqual(result["length_class"], "inner_length_825")

    def test_known_length_can_have_new_sequence(self):
        _, row = fixture(); row["inner_sequence_sha256"] = "e"*64
        result = self.run_case(row)["candidates"][0]
        self.assertFalse(result["outside_panel_length_class"])
        self.assertTrue(result["new_sequence_candidate"])

    def test_unobserved_boundary_and_incomplete_coverage(self):
        _, row = fixture(); row["boundaries"]["right"].update(qnames=0, distinct_pair_starts=0)
        row["coverage"].update(breadth_ge_1=.99, min_read_depth=0)
        result = self.run_case(row)["candidates"][0]
        self.assertEqual(result["primary_failures"], ["incomplete_quality_filtered_coverage", "right_boundary_unobserved"])

    def test_rejects_cross_library_and_source(self):
        for key, value in [("case_id", "other"), ("source_kind", "public_truth")]:
            _, row = fixture(); row[key] = value
            with self.assertRaises(ValueError): self.run_case(row)

    def test_rejects_duplicates_missing_fields_invalid_coordinates(self):
        panel, row = fixture()
        with self.assertRaises(ValueError): recover_structure("sample", [], panel, [row, row])
        del row["source_bindings"]["coverage"]
        with self.assertRaises(ValueError): self.run_case(row)
        _, row = fixture(); row["right_boundary"] = 6000
        with self.assertRaises(ValueError): self.run_case(row)

    def test_rejects_nonfinite_and_impossible_counts(self):
        for value in [float("nan"), float("inf"), True, -1, 1.01]:
            _, row = fixture(); row["coverage"]["breadth_ge_1"] = value
            with self.assertRaises(ValueError): self.run_case(row)
        _, row = fixture(); row["boundaries"]["left"]["distinct_pair_starts"] = 6
        with self.assertRaises(ValueError): self.run_case(row)

    def test_config_is_fixed(self):
        panel, row = fixture(); cfg = copy.deepcopy(DEFAULT_CONFIG); cfg["strict"]["min_boundary_qnames"] = 1
        with self.assertRaises(ValueError): recover_structure("sample", [], panel, [row], cfg)

    def test_does_not_mutate_inputs(self):
        panel, row = fixture(); before=copy.deepcopy((panel,row))
        recover_structure("sample", [], panel, [row])
        self.assertEqual((panel,row),before)

    def test_base_conflict_is_preserved_without_sequence_correction(self):
        _, row = fixture(); row["base_concordance"] = {"positions":2900,"zero_acgt_support_positions":0,
            "assembly_base_no_qname_support_positions":1,"alternative_exceeds_assembly_positions":2}
        result=self.run_case(row)["candidates"][0]
        self.assertTrue(result["assembly_supported"])
        self.assertEqual(result["base_concordance"],row["base_concordance"])
        self.assertEqual(result["interval_sequence_sha256"],row["interval_sequence_sha256"])

    def test_file_binding(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"source.tsv";path.write_bytes(b"original")
            binding={"path":str(path),"bytes":8,"sha256":hashlib.sha256(b"original").hexdigest()}
            validate_file_binding(binding)
            path.write_bytes(b"changed!")
            with self.assertRaises(ValueError): validate_file_binding(binding)

    def test_non_h5_panel_does_not_invent_h5_alias(self):
        panel, row = fixture()
        panel[0].update(inner_length=4898,class_id="custom_4898")
        row.update(contig_length=10000,interval_end=8000,right_boundary=5898,inner_length=4898)
        row["coverage"]["end"]=8000
        row["boundaries"]["right"]["position"]=5898
        result=recover_structure("sample",[],panel,[row])
        self.assertEqual(result["assembly_supported_classes"],["custom_4898"])

    def test_extra_base_audit_counts_are_validated(self):
        for value in ["1", -1, True, float("nan")]:
            _, row=fixture()
            row["base_concordance"]={"positions":2900,"zero_acgt_support_positions":0,
                "assembly_base_no_qname_support_positions":0,"alternative_exceeds_assembly_positions":0,
                "alternate_ge2_qnames_both_orientations_positions":value}
            with self.assertRaises(ValueError): self.run_case(row)


if __name__ == "__main__": unittest.main()
