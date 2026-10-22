"""Synthetic serialization/coverage and authorization tests; no cohort data IO."""
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace as NS
import hashlib
import json
import pickle
import tempfile
import unittest
import numpy as np
import driver
import core
from test_core import ASSAY, RULES, call


class DriverTests(unittest.TestCase):
    def test_public_adapter_requires_one_named_candidate(self):
        base={"schema":"single-candidate-readback-rules.public-adapter.v1","status":"FROZEN_FOR_SINGLE_CANDIDATE_READBACK","candidate_count":1,"sample_id":"SYNTHETIC_SAMPLE"}
        for changed in ({"sample_id":""},{"sample_id":None},{"candidate_count":2},{"schema":"other"}):
            with self.assertRaises(ValueError):
                driver.verify_authorization({**base,**changed},"a","a",True)

    def test_draft_cannot_run_or_create_outputs(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"rules.json";p.write_text(json.dumps({"status":"DRAFT_NOT_AUTHORIZED"}))
            args=NS(rules=str(p),rules_sha256=driver.digest(p),out=str(Path(d)/"out"),root_authorized_run=True)
            with self.assertRaises(ValueError):driver.run(args)
            self.assertFalse(Path(args.out).exists())

    def test_exact_sha_and_explicit_authorization(self):
        rules={"schema":"single-candidate-readback-rules.public-adapter.v1","status":"FROZEN_FOR_SINGLE_CANDIDATE_READBACK","candidate_count":1,"sample_id":"SYNTHETIC_SAMPLE"}
        driver.verify_authorization(rules,"sha","sha",True)
        for actual,expected,authorized in [("a","b",True),("a","a",False)]:
            with self.assertRaises(ValueError):driver.verify_authorization(rules,actual,expected,authorized)

    def test_rdic_single_stream_hash_and_shape(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"Rdic.pickle";p.write_bytes(pickle.dumps({"r1":{"q1","q2"},"unused":{"x"}},protocol=5))
            self.assertEqual(driver.load_rdic(p,driver.digest(p),["r1"]),{"r1":{"q1","q2"}})
            with self.assertRaises(ValueError):driver.load_rdic(p,"wrong",["r1"])
            with self.assertRaises(ValueError):driver.load_rdic(p,driver.digest(p),["absent"])

    def test_pickle_globals_forbidden(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"Rdic.pickle";p.write_bytes(pickle.dumps({"r1":Path("not primitive")}))
            with self.assertRaises(pickle.UnpicklingError):driver.load_rdic(p,driver.digest(p),["r1"])

    def test_covt_original_values_then_positions_axis(self):
        # Deliberately distinct scales, so swapping values/coordinates cannot pass.
        h={"r1::0":np.array([[2,4,7,8],[99,101,102,1000]])}
        self.assertEqual(driver.read_cov_window(h,"r1",100,104),{101:4,102:7})
        self.assertIsNone(driver.read_cov_window(h,"absent",100,104))
        with self.assertRaises(ValueError):driver.read_cov_window({"r1::1":h["r1::0"]},"r1",100,104)
        with self.assertRaises(ValueError):driver.read_cov_window({"r1::0":np.array([[2,4],[101,101]])},"r1",100,104)

    def test_parity_failure_blocks_positive_label_and_counts(self):
        summary,main,diag=driver.interpret_calls([call()], [ASSAY],asdict(RULES),{"all_windows_match":False})
        self.assertEqual(summary["status"],"BLOCKED_COVT_PARITY")
        self.assertIsNone(summary["technical_label"])
        self.assertFalse(summary["event_counts_generated"])
        self.assertEqual((main,diag),([],[]))

    def test_partial_parity_failure(self):
        class Bam:
            def pileup(self,ref,**kw):return []
        windows=[{"reference":"r1","comparison_window0":[100,104],"original_splits0_inclusive":[[0,4999]]},
                 {"reference":"r2","comparison_window0":[200,204],"original_splits0_inclusive":[[0,4999]]}]
        cov={"r1::0":np.empty((2,0)),"r2::0":np.array([[1],[201]])}
        p=driver.parity_report(Bam(),cov,{"r1":set(),"r2":set()},windows)
        self.assertFalse(p["all_windows_match"])
        self.assertEqual([w["status"] for w in p["windows"]],["MATCH","MISMATCH_BLOCK_EVENT_INTERPRETATION"])

    def test_stat_and_fresh_output(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"source";p.write_bytes(b"synthetic")
            rec={"path":str(p),"historic_bytes":9,"historic_mtime_ns":p.stat().st_mtime_ns,"sha256":driver.digest(p)}
            self.assertIsNone(driver.stat_bound(rec)["current_full_sha256"])
            p.write_bytes(b"changed")
            with self.assertRaises(ValueError):driver.stat_bound(rec)


if __name__=="__main__":unittest.main()
