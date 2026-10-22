"""Synthetic tests only; no real sample data or phenotype files are read."""
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import profile_sequence_audit as audit

RUNTIME = all(importlib.util.find_spec(x) is not None for x in ("Bio", "h5py", "numpy"))
RULES = audit.read_json(Path(audit.__file__).with_name("PROFILE_RULES.json"))


def variant(ref="A", con="T", counts=None, alleles=1, cls="SNS", cryptic=False):
    counts = counts or {"A": 0, "C": 0, "G": 0, "T": 10}
    total = sum(counts.values())
    secondary = max((b for b in "ACGT" if b != con), key=lambda b: counts[b])
    return {"ref_base": ref, "con_base": con, "var_base": secondary, "A": str(counts["A"]), "C": str(counts["C"]),
            "G": str(counts["G"]), "T": str(counts["T"]), "allele_count": str(alleles), "class": cls,
            "cryptic": str(cryptic), "position_coverage": str(total), "mm": "0",
            "con_freq": str(counts[con]/total), "var_freq": str(counts[secondary]/total), "ref_freq": str(counts.get(ref, 0)/total)}


class ConsensusTests(unittest.TestCase):
    def test_positive_strand_substitution(self):
        sequence, _, mixed, _ = audit.make_consensus("ATGAAATAA", [10]*9, {3: variant()}, 0, 1)
        self.assertEqual(sequence, "ATGTAATAA")
        self.assertFalse(mixed)

    def test_negative_strand_coordinate_and_complement(self):
        # Reverse complement of ATGAAATAA is TTATTTCAT. T->A at forward 105
        # changes coding position 4 from A to T when the locus starts at 100.
        ref = "TTATTTCAT"
        row = variant(ref="T", con="A", counts={"A": 10, "C": 0, "G": 0, "T": 0})
        sequence, _, _, _ = audit.make_consensus(ref, [10]*9, {105: row}, 100, -1)
        self.assertEqual(sequence, "ATGTAATAA")

    def test_depth4_masks_even_without_snv_row(self):
        sequence, reasons, _, _ = audit.make_consensus("ATGAAATAA", [10, 4, 0, 10, 10, 10, 10, 10, 10], {})
        self.assertEqual(sequence, "ANNAAATAA")
        self.assertEqual(reasons[:3], ["", "depth_below_5", "zero_depth"])

    def test_missing_scaffold_is_all_unknown(self):
        sequence, reasons, _, _ = audit.make_consensus("ATGAAATAA", None, {})
        self.assertEqual(sequence, "N"*9)
        self.assertEqual(set(reasons), {"coverage_unavailable"})

    def test_ambiguous_reference_masked(self):
        sequence, reasons, _, _ = audit.make_consensus("ATGNAATAA", [10]*9, {})
        self.assertEqual(sequence[3], "N")
        self.assertEqual(reasons[3], "ambiguous_reference")

    def test_mixed_non_tied_keeps_original_major(self):
        row = variant(counts={"A": 3, "C": 0, "G": 0, "T": 7}, alleles=2, cls="con_SNV")
        sequence, _, mixed, _ = audit.make_consensus("ATGAAATAA", [10]*9, {3: row})
        self.assertEqual(sequence, "ATGTAATAA")
        self.assertEqual(mixed, {3})
        _, _, _, stops = audit.protein_effects("ATGAAATAA", sequence, mixed)
        self.assertEqual(stops, [{"aa_position1": 2, "codon_mixed": True}])

    def test_tied_major_is_masked_not_arbitrarily_picked(self):
        row = variant(counts={"A": 5, "C": 0, "G": 0, "T": 5}, alleles=2, cls="SNV")
        sequence, reasons, mixed, _ = audit.make_consensus("ATGAAATAA", [10]*9, {3: row})
        self.assertEqual(sequence[3], "N")
        self.assertEqual(reasons[3], "tied_major_alleles")
        self.assertEqual(mixed, {3})

    def test_cryptic_original_call_masks(self):
        sequence, reasons, _, _ = audit.make_consensus("ATGAAATAA", [10]*9, {3: variant(cryptic=True)})
        self.assertEqual(sequence[3], "N")
        self.assertEqual(reasons[3], "original_call_cryptic")

    def test_unresolved_original_caller_masks(self):
        row = variant(alleles=0, cls="DivergentSite")
        sequence, reasons, _, _ = audit.make_consensus("ATGAAATAA", [10]*9, {3: row})
        self.assertEqual(sequence[3], "N")
        self.assertEqual(reasons[3], "original_caller_unresolved")

    def test_depth_or_reference_mismatch_fails_closed(self):
        with self.assertRaises(ValueError):
            audit.make_consensus("ATGAAATAA", [9]*9, {3: variant()})
        with self.assertRaises(ValueError):
            audit.make_consensus("ATGCAATAA", [10]*9, {3: variant()})

    def test_unrecognized_caller_class_fails_closed(self):
        with self.assertRaises(ValueError):
            audit.call_existing(variant(cls="future_unknown"), "A", 10)

    def test_premature_stop_and_missing_codon(self):
        ref, aa, changes, stops = audit.protein_effects("ATGAAATAA", "ATGTAATAA", set())
        self.assertEqual((ref, aa), ("MK*", "M**"))
        self.assertEqual(stops, [{"aa_position1": 2, "codon_mixed": False}])
        self.assertEqual(audit.protein_effects("ATGAAATAA", "ATGNAATAA", set())[3], [])
        self.assertEqual(audit.protein_effects("ATGAAATAA", "ATGAAATAA", set())[3], [])

    def test_alternative_table11_start(self):
        self.assertEqual(audit.translate("GTGAAATAA"), "MK*")
        self.assertEqual(audit.translate("NTGAAATAA"), "XK*")

    def test_negative_mixed_offset(self):
        row = variant(ref="T", con="A", counts={"A": 7, "C": 0, "G": 0, "T": 3}, alleles=2, cls="con_SNV")
        seq, _, mixed, _ = audit.make_consensus("TTATTTCAT", [10]*9, {105: row}, 100, -1)
        self.assertEqual(mixed, {3})
        self.assertEqual(audit.protein_effects("ATGAAATAA", seq, mixed)[3][0]["codon_mixed"], True)

    def test_mask_runs_half_open(self):
        self.assertEqual(audit.mask_runs(["", "zero_depth", "zero_depth", "depth_below_5", ""], 100),
                         [{"start0": 101, "end0": 103, "length": 2, "reason": "zero_depth"},
                          {"start0": 103, "end0": 104, "length": 1, "reason": "depth_below_5"}])


class SourceTests(unittest.TestCase):
    def test_source_hash_fail_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/"source"
            path.write_text("one")
            binding = audit.Bindings()
            original = audit.record(path)
            binding.add(path, original)
            path.write_text("two")
            with self.assertRaises(ValueError):
                binding.stable()
            with self.assertRaises(ValueError):
                audit.Bindings().add(path, original)

    def test_raw_cumulative_disagreement_fails(self):
        row = variant()
        raw = {("c", 0): row}
        audit.crosscheck_variant_tables(raw, raw)
        changed = dict(row, con_base="C")
        with self.assertRaises(ValueError):
            audit.crosscheck_variant_tables(raw, {("c", 0): changed})
        with self.assertRaises(ValueError):
            audit.crosscheck_variant_tables(raw, {})

    def test_no_existing_directory_overwrite(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(FileExistsError):
                audit.run_audit(Path(temp), "SYNTHETIC", Path(temp))

    def test_incomplete_collection_fails_and_preserves_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root/"samples.tsv").write_text("sample_id\nSYNTHETIC\n")
            with self.assertRaises(ValueError):
                audit.collect(root/"samples.tsv", root, root/"failed")
            self.assertTrue((root/"failed/FAILED.json").exists())
            self.assertFalse((root/"failed/COMPLETE").exists())


class LocusTests(unittest.TestCase):
    def fixture(self, historical=False):
        target = {"locus_id": "historic" if historical else "gene", "trait_id": "NanH2" if historical else "AfcA", "role": "positive_enzyme_seed",
                  "assembly_accession": "genome", "contig_id": "contig", "start0": 0, "end0": 9, "strand": 1,
                  "evidence_type": "curated_historical_CDS" if historical else "ordinary_reference_CDS",
                  "annotation_conflict": historical, "reference_cds_sha256": "synthetic"}
        ev = {("" if historical else "gene_")+k: v for k,v in {"mean_depth": 10., "breadth_depth1": 1., "breadth_depth5": 1.}.items()}
        ev["region_status" if historical else "gene_status"] = "historical_region_supported_annotation_conflict" if historical else "supported"
        ev["annotation_conflict"] = historical
        carrier = {"breadth_depth5": 1., "background_median_depth": 10., "background_supported_cds_count": 600, "background_eligible": True}
        upstream = {"scaffold": "contig", "start": "0", "end": "8", "coverage": "10", "breadth": "1", "breadth_minCov": "1"}
        return target, ev, carrier, upstream

    def test_historical_conflict_retained_not_intact_orf(self):
        t,e,c,u = self.fixture(True)
        result = audit.locus_result("SYNTHETIC", t, "ATGAAATAA", [10]*9, {}, e, c, None, "MK", RULES)
        summary = result[0]
        self.assertTrue(summary["annotation_conflict"])
        self.assertTrue(summary["complete_substitution_consensus"])
        self.assertFalse(summary["intact_sample_orf_asserted"])
        self.assertFalse(summary["indels_assessed"])
        self.assertIn("historical_annotation_conflict", summary["flags"])

    def test_ordinary_complete_consensus_still_not_orf(self):
        t,e,c,u = self.fixture()
        summary = audit.locus_result("SYNTHETIC", t, "ATGAAATAA", [10]*9, {}, e,c,u,"MK",RULES)[0]
        self.assertEqual(set(summary), set(audit.LOCUS_FIELDS))
        self.assertTrue(summary["full_depth5_callable"])
        self.assertFalse(summary["intact_sample_orf_asserted"])
        self.assertFalse(summary["full_length_phasing_assessed"])

    def test_comparator_denominator_separate(self):
        t,e,c,u = self.fixture()
        one = audit.locus_result("SYNTHETIC", t, "ATGAAATAA", [10]*9, {},e,c,u,"MK",RULES)[0]
        two = dict(one, role="specificity_comparator", trait_id="Bga42B")
        ds = audit.denominators("SYNTHETIC", [one,two])
        self.assertEqual({(r["role"], r["fixed_locus_denominator"]) for r in ds}, {("positive_enzyme_seed",1),("specificity_comparator",1)})

    def test_identical_reference_seed_needs_no_alignment(self):
        self.assertEqual(audit.seed_positions("MK*", "MK", RULES), ("identical_reference_seed", {0:0,1:1}))

    @unittest.skipUnless(RUNTIME, "installed HPC scientific runtime required")
    def test_seed_alignment_ambiguity_suppresses_residue_claim(self):
        status, mapping = audit.seed_positions("MAAA", "MAA", RULES)
        self.assertEqual(status, "ambiguous_optimal_reference_seed_alignment")
        self.assertEqual(mapping, {})


@unittest.skipUnless(RUNTIME, "installed HPC scientific runtime required")
class RuntimeSyntheticTests(unittest.TestCase):
    def test_hdf5_subset_and_absent_scaffold(self):
        import h5py
        import numpy as np
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/"toy.hd5"
            with h5py.File(path,"w") as handle:
                handle["c::0"] = np.array([[2,5,9],[1,3,5]], dtype="int64")
            with h5py.File(path,"r") as handle:
                self.assertEqual(audit.dense_locus_cov(handle,"c",2,5,6).tolist(), [0,5,0])
                self.assertIsNone(audit.dense_locus_cov(handle,"missing",0,3,6))

    def test_synthetic_complete_16_grid_and_failure_marker(self):
        traits=["AfcA","AfcB","LnbB","Bga42A","Bga2A","NanH2"]
        identities=[(t,"positive_enzyme_seed") for t in traits for _ in range(2)] + [("Bga42B","specificity_comparator")]*2 + [("Bga42C","specificity_comparator"),("NanH1","specificity_comparator")]
        roster=[]; genes={}; regions={}; carriers={"g":{"breadth_depth5":"","background_median_depth":"","background_supported_cds_count":"","background_eligible":False}}
        for i,(trait,role) in enumerate(identities):
            historical=i==11
            t={"locus_id":"toy"+str(i),"trait_id":trait,"role":role,"assembly_accession":"g","contig_id":"c","start0":0,"end0":9,"strand":1,"annotation_conflict":historical,"evidence_type":"curated_historical_CDS" if historical else "ordinary_reference_CDS","reference_cds_sha256":"synthetic"}
            roster.append(t)
            ev={"trait_id":trait,"assembly_accession":"g","contig_id":"c","annotation_conflict":historical}
            prefix="" if historical else "gene_"
            ev.update({prefix+k:"" for k in ("mean_depth","breadth_depth1","breadth_depth5")})
            ev["region_status" if historical else "gene_status"]="unresolved"
            (regions if historical else genes)[t["locus_id"]]=ev
        context=(roster,{"c":"ATGAAATAA"},{t:"MK" for t,_ in identities},{})
        profile=({"profile_status":"NO_ANALYSABLE_PAIRS_NOT_BIOLOGICAL_ABSENCE"},genes,regions,carriers,None,None,None,{},None)
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            with patch.object(audit,"load_roster",return_value=context), patch.object(audit,"load_profile",return_value=profile):
                result=audit.run_audit(root,"SYNTHETIC",root/"ok")
            self.assertEqual(result["loci"],16)
            summary=audit.read_json(root/"ok/SUMMARY.json")
            self.assertEqual(summary["row_counts"]["trait_denominators.tsv"],9)
            self.assertEqual((root/"ok/COMPLETE").read_text().strip(), audit.sha(root/"ok/SUMMARY.json"))
            self.assertTrue(all(r["consensus_N_bases"]=="9" for r in audit.table(root/"ok/locus_summary.tsv")))
            # Construct a fully synthetic 50-sample sealed grid to check collection.
            sample_ids = ["SYNTHETIC_"+str(i) for i in range(50)]
            for sample in sample_ids:
                folder = root/"samples"/sample
                folder.mkdir(parents=True)
                cloned = dict(summary, sample_id=sample)
                outputs=[]
                for entry in summary["outputs"]:
                    name=entry["path"]
                    content=(root/"ok"/name).read_text()
                    (folder/name).write_text(content.replace("SYNTHETIC",sample))
                    outputs.append({"path":name,"bytes":(folder/name).stat().st_size,"sha256":audit.sha(folder/name)})
                cloned["outputs"]=outputs
                cloned["roster_sha256"]=audit.sha(folder/"ROSTER.tsv")
                (folder/"SUMMARY.json").write_text(json.dumps(cloned))
                (folder/"COMPLETE").write_text(audit.sha(folder/"SUMMARY.json")+"\n")
            manifest=root/"manifest.tsv"
            manifest.write_text("sample_id\n"+"\n".join(sample_ids)+"\n")
            result=audit.collect(manifest,root/"samples",root/"merged")
            self.assertEqual(result["locus_rows"],800)
            self.assertEqual(len(audit.table(root/"merged/trait_denominators.tsv")),450)
            tampered=root/"samples"/sample_ids[0]/"variants.tsv"
            tampered.write_text(tampered.read_text()+"corruption\n")
            with self.assertRaises(ValueError):
                audit.collect(manifest,root/"samples",root/"corrupt_collection")
            self.assertFalse((root/"corrupt_collection/COMPLETE").exists())
            with patch.object(audit,"load_roster",side_effect=ValueError("synthetic source failure")):
                with self.assertRaises(ValueError):
                    audit.run_audit(root,"SYNTHETIC",root/"bad")
            self.assertTrue((root/"bad/FAILED.json").exists())
            self.assertFalse((root/"bad/COMPLETE").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
