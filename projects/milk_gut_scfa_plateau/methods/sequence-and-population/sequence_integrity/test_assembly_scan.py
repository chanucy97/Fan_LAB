"""Small synthetic checks; no cohort execution or scheduler access."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import assembly_scan as a


class AssemblyScanTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write(self, name, text):
        path = self.root/name
        path.write_text(text, encoding="ascii")
        return path

    def trait_fixture(self):
        return {**{f"p{i}": {"role": "positive_enzyme_seed", "annotation_conflict": "", "claim_ceiling": "DNA candidate"} for i in range(6)},
                **{f"c{i}": {"role": "specificity_comparator", "annotation_conflict": "", "claim_ceiling": "comparator"} for i in range(3)}}

    def module_fixture(self):
        return [{"module_id": f"m{i}", "required_all": "p0;p1"} for i in range(4)]

    def test_complete_and_alternative_start_translation(self):
        for start in ("ATG", "GTG", "TTG"):
            row = a.inspect_orf(start+"AAATAA", "MK*", "00", start)
            self.assertTrue(row["complete_predicted_orf"])
        self.assertFalse(a.inspect_orf("ATGTAATAA", "M**", "00", "ATG")["complete_predicted_orf"])
        self.assertFalse(a.inspect_orf("ATGNNNTAA", "MX*", "00", "ATG")["complete_predicted_orf"])

    def test_partial_edge_not_forced_methionine(self):
        for start, amino in (("TTG", "L"), ("GTG", "V")):
            result = a.inspect_orf(start+"AAATAA", amino+"K*", "10", "Edge")
            self.assertTrue(result["translation_matches"])
            self.assertFalse(result["complete_predicted_orf"])
            self.assertFalse(a.inspect_orf(start+"AAATAA", "MK*", "10", "Edge")["translation_matches"])

    def test_translation_mismatch_and_frame_recorded(self):
        result = a.inspect_orf("ATGAAATAA", "MA*", "00", "ATG")
        self.assertFalse(result["translation_matches"])
        self.assertFalse(result["complete_predicted_orf"])
        self.assertFalse(a.inspect_orf("ATGAAATA", "MK", "11", "Edge")["in_frame"])
        with self.assertRaises(ValueError):
            a.inspect_orf("ATGAAATAA", "MK*", "0", "ATG")

    def prediction_fixture(self, strand):
        sequence = "ATGAAATAA"
        assembly = {"contig_with_underscores": sequence if strand == 1 else a.revcomp(sequence)}
        attrs = "ID=1_1;partial=00;start_type=ATG;rbs_motif=None"
        header = f"contig_with_underscores_1 # 1 # 9 # {strand} # {attrs}"
        faa = self.write("test.faa", f">{header}\nMK*\n")
        fna = self.write("test.fna", f">{header}\n{sequence}\n")
        gff = self.write("test.gff", f"##gff-version 3\ncontig_with_underscores\tProdigal_v2.6.3\tCDS\t1\t9\t20\t{'+' if strand == 1 else '-'}\t0\t{attrs};conf=99.0\n")
        return assembly, faa, fna, gff

    def test_actual_negative_strand_coordinates_and_hash(self):
        assembly, faa, fna, gff = self.prediction_fixture(-1)
        orfs, genes, proteins = a.parse_predictions(assembly, faa, fna, gff, "synthetic")
        self.assertEqual(orfs[0]["strand"], -1)
        self.assertEqual(orfs[0]["cds_sha256"], a.seqsha("ATGAAATAA"))
        self.assertEqual(orfs[0]["protein_sha256"], a.seqsha("MK"))
        self.assertEqual(genes[orfs[0]["gene_id"]]["assembly_accession"], "synthetic")
        self.assertTrue(orfs[0]["complete_predicted_orf"])
        self.assertEqual(list(proteins.values()), ["MK"])

    def test_prediction_coordinate_or_cds_mismatch_fails(self):
        assembly, faa, fna, gff = self.prediction_fixture(1)
        assembly["contig_with_underscores"] = "ATGCCCTAA"
        with self.assertRaisesRegex(ValueError, "oriented assembly"):
            a.parse_predictions(assembly, faa, fna, gff, "synthetic")

    def test_fasta_duplicate_and_empty_fail(self):
        with self.assertRaises(ValueError):
            a.fasta(self.write("dup.fa", ">a\nAAA\n>a\nAAA\n"))
        with self.assertRaises(ValueError):
            a.fasta(self.write("empty_record.fa", ">a\n"))
        self.assertEqual(a.fasta(self.write("empty.fa", "")), ({}, {}))

    def test_tblastn_negative_and_far_apart_hsps_never_become_orf(self):
        seeds, assembly = {"p0": "M"*10}, {"contig": "A"*1000}
        line1 = "p0\tcontig\t100\t5\t10\t1000\t1\t5\t990\t976\t-1\t1e-20\t80\t0\t0\t5\n"
        line2 = "p0\tcontig\t100\t5\t10\t1000\t6\t10\t30\t16\t-1\t1e-20\t80\t0\t0\t5\n"
        hits, summary = a.parse_tblastn(self.write("tb.tsv", line1+line2), seeds, assembly, self.trait_fixture(), {"minimum_identity_percent": 90, "maximum_evalue": 1e-10})
        self.assertEqual(len(hits), 2)
        self.assertEqual(summary[0]["query_union_aa"], 10)
        self.assertEqual(summary[0]["high_identity_query_union_aa"], 10)
        self.assertEqual(summary[0]["query_union_percent"], 100)
        self.assertEqual(summary[0]["high_identity_query_union_percent"], 100)
        self.assertEqual(set(summary[0]), set(a.TB_SUMMARY_FIELDS))
        self.assertEqual(summary[0]["largest_uncovered_subject_gap_nt"], 945)
        self.assertFalse(summary[0]["single_orf_asserted"])
        self.assertEqual(hits[0]["btop"], "5")
        with self.assertRaisesRegex(ValueError, "frame/strand"):
            a.parse_tblastn(self.write("bad.tsv", line1.replace("\t-1\t", "\t1\t")), seeds, assembly, self.trait_fixture(), {"minimum_identity_percent": 90, "maximum_evalue": 1e-10})

    def test_union_counts_overlapping_query_once(self):
        self.assertEqual(a.union_length([(1, 50), (40, 90), (100, 100)]), 91)

    def test_empty_assembly_retains_all_denominators(self):
        traits, modules, contigs = a.evidence_tables("empty", "empty_assembly", self.trait_fixture(), self.module_fixture(), [], [], [], [], [])
        self.assertEqual(len(traits), 6)
        self.assertEqual(len(modules), 4)
        self.assertEqual(contigs, [])
        self.assertEqual({r["status"] for r in traits}, {"assembly_not_evaluable"})
        self.assertTrue(all(not r["negative_is_absence"] for r in traits))

    def test_comparator_not_a_positive_and_ambiguity_not_module(self):
        assignments = [{"trait_id": "p0", "trait_role": "positive_enzyme_seed", "gene_id": "g0", "mapping_status": "exact_seed_peptide", "competitor_ambiguity": True},
                       {"trait_id": "c0", "trait_role": "specificity_comparator", "gene_id": "gc", "mapping_status": "exact_seed_peptide", "competitor_ambiguity": False}]
        orfs = [{"gene_id": "g0", "complete_predicted_orf": True, "contig_id": "ctg"}, {"gene_id": "gc", "complete_predicted_orf": True, "contig_id": "ctg"}]
        traits, modules, contigs = a.evidence_tables("x", "assembly_scanned", self.trait_fixture(), self.module_fixture(), orfs, assignments, [], [], [])
        self.assertEqual(traits[0]["status"], "homology_ambiguous_candidate")
        self.assertEqual(len(traits), 6)
        self.assertEqual(contigs, [])
        self.assertTrue(all(r["status"] == "candidate_set_not_resolved" for r in modules))

    def test_same_contig_vs_dispersed(self):
        assignments = [{"trait_id": f"p{i}", "trait_role": "positive_enzyme_seed", "gene_id": f"g{i}", "mapping_status": "high_similarity_candidate", "competitor_ambiguity": False} for i in range(2)]
        orfs = [{"gene_id": f"g{i}", "complete_predicted_orf": True, "contig_id": f"ctg{i}"} for i in range(2)]
        _, modules, _ = a.evidence_tables("x", "assembly_scanned", self.trait_fixture(), self.module_fixture(), orfs, assignments, [], [], [])
        self.assertEqual(modules[0]["status"], "complete_candidates_dispersed_in_sample")
        orfs[1]["contig_id"] = "ctg0"
        _, modules, contigs = a.evidence_tables("x", "assembly_scanned", self.trait_fixture(), self.module_fixture(), orfs, assignments, [], [], [])
        self.assertEqual(modules[0]["status"], "complete_candidates_same_assembly_contig")
        self.assertFalse(modules[0]["same_strain_asserted"])
        self.assertTrue(contigs[0]["assembly_colocation_only"])

    def test_original_classifier_comparator_and_exact_ambiguity(self):
        mapper = a.load_mapper(a.HERE.parent/"trait_mapping/map_traits.py")
        traits = self.trait_fixture()
        seeds = {"p0": "MK", "c0": "MA"}
        proteins = {"g0": "MK"}
        genes = {"g0": {"assembly_accession": "synthetic", "contig_id": "ctg", "start0": 0, "end0": 9, "strand": 1, "locus_tag": "", "old_locus_tag": "", "protein_id": "g0", "instrain_gene_exported": "False"}}
        def hit(query, subject, score):
            return {"qseqid": query, "sseqid": subject, "pident": 100, "qcov": 100, "scov": 100, "bitscore": score, "evalue": 1e-20, "qstart": 1, "sstart": 1}
        rules = a.read_json(a.HERE.parent/"trait_mapping/RULES.json")
        data = mapper.classify(seeds, traits, genes, proteins, [hit("g0", "p0", 100), hit("g0", "c0", 95)], [hit("p0", "g0", 100)], rules)
        self.assertEqual(data[0]["mapping_status"], "exact_seed_peptide")
        self.assertTrue(data[0]["competitor_ambiguity"])
        self.assertIn("c0", data[0]["seed_near_best_ids"])

    def test_fresh_output_and_completion_binding(self):
        out = a.fresh(self.root/"run")
        with self.assertRaises(FileExistsError):
            a.fresh(out)
        (out/"result.txt").write_text("synthetic", encoding="ascii")
        a.complete_run(out, {"schema": "synthetic.v1"})
        self.assertEqual(a.verify_complete(out)["status"], "COMPLETE")
        (out/"result.txt").write_text("changed", encoding="ascii")
        with self.assertRaisesRegex(ValueError, "changed"):
            a.verify_complete(out)

    def test_unbound_output_or_failed_marker_rejected(self):
        out = a.fresh(self.root/"run")
        a.complete_run(out, {"schema": "synthetic.v1"})
        (out/"extra").write_text("x")
        with self.assertRaisesRegex(ValueError, "unbound"):
            a.verify_complete(out)
        (out/"extra").unlink()
        (out/"FAILED.json").write_text("{}")
        with self.assertRaisesRegex(ValueError, "completion"):
            a.verify_complete(out)
        with self.assertRaisesRegex(ValueError, "completed or failed"):
            a.complete_run(out, {"schema": "synthetic.v1"})

    def test_empty_run_end_to_end_with_version_stubs(self):
        assembly = self.write("empty.fa", "")
        config_path = self.write("config.json", "{}")
        binding = a.record(config_path)
        sample = {"sample_id": "synthetic_empty", "dyad": "synthetic", "window": "toy", "historical_contigs": 0,
                  "historical_bases": 0, "assembly": a.record(assembly)}
        prepared = a.fresh(self.root/"prepared")
        a.write_json(prepared/"PREPARED.json", {"schema": "sequence-integrity.assembly-prepared.v1", "config_path": str(config_path),
                                              "bindings": [binding], "samples": [sample], "rules": {}})
        a.complete_run(prepared, {"schema": "sequence-integrity.assembly-prepare-run.v1", "bindings": [binding]})
        config = {"tools": {k: "synthetic_version_only" for k in ("prodigal", "diamond", "tblastn", "makeblastdb")},
                  "mapping_source": str(a.HERE.parent/"trait_mapping/map_traits.py"), "mapping_rules": str(a.HERE.parent/"trait_mapping/RULES.json")}
        def stub_version(runner, name, argv):
            self.assertTrue(name.endswith("_version"))
            runner.records.append({"name": name, "argv": list(argv), "returncode": 0, "synthetic": True})
            a.write_json(runner.out/"COMMANDS.json", runner.records)
        fixture = (config, {}, [binding], [sample], {"unused": "M"}, self.trait_fixture(), self.module_fixture())
        with patch.object(a, "source_inputs", return_value=fixture), patch.object(a.Commands, "run", new=stub_version):
            result = a.run_sample(prepared, 0, self.root/"sample")
        self.assertEqual(result["assembly_state"], "empty_assembly")
        self.assertEqual(len(a.rows(self.root/"sample/TRAIT_STATUS.tsv")), 6)
        self.assertEqual(len(a.rows(self.root/"sample/MODULE_SAMPLE_STATUS.tsv")), 4)
        self.assertEqual(len(a.read_json(self.root/"sample/COMMANDS.json")), 4)
        self.assertEqual(a.verify_complete(self.root/"sample")["status"], "COMPLETE")
        self.assertFalse((self.root/"sample/predicted.gff").exists())

    def test_failed_command_keeps_logs_and_returncode(self):
        out = a.fresh(self.root/"commands")
        runner = a.Commands(out)
        fake = type("Result", (), {"returncode": 9})()
        with patch.object(a.subprocess, "run", return_value=fake):
            with self.assertRaises(RuntimeError):
                runner.run("synthetic", ["fake_binary", "arg"])
        self.assertEqual(a.read_json(out/"COMMANDS.json")[0]["returncode"], 9)
        self.assertTrue((out/"logs/synthetic.stderr.txt").exists())


if __name__ == "__main__":
    unittest.main()
