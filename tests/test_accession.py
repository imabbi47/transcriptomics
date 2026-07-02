"""Unit tests for accession classification (pure logic, no network)."""
from __future__ import annotations

import unittest

from transcriptomics.accession import AccessionType, classify, parse_many


class TestClassify(unittest.TestCase):
    def test_geo(self):
        self.assertEqual(classify("GSE12345").type, AccessionType.GEO_SERIES)
        self.assertEqual(classify("GSM98765").type, AccessionType.GEO_SAMPLE)

    def test_insdc_runs(self):
        for acc in ("SRR000001", "ERR1234567", "DRR009999"):
            self.assertEqual(classify(acc).type, AccessionType.RUN, acc)

    def test_experiment_sample_study(self):
        self.assertEqual(classify("SRX123").type, AccessionType.EXPERIMENT)
        self.assertEqual(classify("ERS999").type, AccessionType.SAMPLE)
        self.assertEqual(classify("SRP000123").type, AccessionType.STUDY)
        self.assertEqual(classify("DRP000001").type, AccessionType.STUDY)

    def test_bioproject_biosample(self):
        self.assertEqual(classify("PRJNA63443").type, AccessionType.BIOPROJECT)
        self.assertEqual(classify("PRJEB12345").type, AccessionType.BIOPROJECT)
        self.assertEqual(classify("SAMN0000123").type, AccessionType.BIOSAMPLE)
        self.assertEqual(classify("SAMEA104").type, AccessionType.BIOSAMPLE)

    def test_case_insensitive_and_whitespace(self):
        acc = classify("  srr000001 ")
        self.assertEqual(acc.type, AccessionType.RUN)
        self.assertEqual(acc.value, "SRR000001")  # normalised upper-case

    def test_unknown(self):
        acc = classify("not_an_accession")
        self.assertEqual(acc.type, AccessionType.UNKNOWN)
        self.assertFalse(acc.is_resolvable)

    def test_helpers(self):
        self.assertTrue(classify("GSE1").is_geo)
        self.assertTrue(classify("SRR1").is_run)
        self.assertFalse(classify("SRP1").is_run)

    def test_parse_many_skips_blanks(self):
        parsed = parse_many(["SRR1", "", "  ", "GSE2"])
        self.assertEqual([p.value for p in parsed], ["SRR1", "GSE2"])


if __name__ == "__main__":
    unittest.main()
