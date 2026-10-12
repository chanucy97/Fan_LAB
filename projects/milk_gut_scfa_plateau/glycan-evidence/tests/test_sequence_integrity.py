"""Sequence outputs are part of the independently verifiable artifact set."""
import tempfile
import unittest
from pathlib import Path

from glycan_evidence.integrity import seal_outputs, verify_outputs


class SequenceIntegrityTests(unittest.TestCase):
    def test_each_coding_sequence_output_is_bound_and_tampering_is_rejected(self):
        for suffix, content in (("fna", ">fixture\nATGAAATAA\n"), ("faa", ">fixture\nMK\n")):
            with self.subTest(suffix=suffix), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                (root / "COMPLETE").write_text("done\n", encoding="utf-8")
                path = root / ("localized." + suffix)
                path.write_text(content, encoding="utf-8")
                result = seal_outputs(root)
                self.assertEqual(result["report_files"], 2)
                replacement = content.replace("K", "A") if suffix == "faa" else content.replace("AAA", "CCC")
                self.assertEqual(len(content), len(replacement))
                path.write_text(replacement, encoding="utf-8")
                with self.assertRaisesRegex(ValueError, "Artifact failed integrity check"):
                    verify_outputs(root)


if __name__ == "__main__":
    unittest.main()
