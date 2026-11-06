import csv
import gzip
import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from dbsnp_liftover.core import (
    LiftoverError,
    build_index,
    convert_bigbed_manifest,
    convert_manifest,
    input_effect_is_source_alt,
    source_alt_maps_to_target_alt,
)


class LiftoverPipelineTests(unittest.TestCase):
    def test_complementary_alleles_are_matched(self):
        self.assertTrue(input_effect_is_source_alt("C", "T", "A", "G"))
        self.assertTrue(source_alt_maps_to_target_alt("A", "G", "T", "C"))
        self.assertFalse(source_alt_maps_to_target_alt("A", "G", "C", "T"))

    def _write_vcf(self, path: Path, rows):
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "wt", encoding="utf-8") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for row in rows:
                handle.write("\t".join(row) + "\n")

    def test_build_index_and_flip_effects_with_duplicate_input_variant(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_vcf = root / "source.vcf.gz"
            target_vcf = root / "target.vcf.gz"
            self._write_vcf(source_vcf, [["chr1", "100", "rs10", "A", "G", ".", ".", "."]])
            self._write_vcf(target_vcf, [["1", "80", "rs10", "G", "A", ".", ".", "."]])
            source_db, target_db = root / "source.sqlite", root / "target.sqlite"
            build_index(source_vcf, source_db, "GRCh38", "155")
            build_index(target_vcf, target_db, "GRCh37", "155")

            input_path = root / "study.tsv"
            input_path.write_text(
                "CHR\tPOS\tEA\tOA\tBETA\tEAF\n"
                "chr1\t100\tG\tA\t0.25\t0.2\n"
                "chr01\t0100\tg\ta\t-0.5\t0.7\n",
                encoding="utf-8",
            )
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "sample_id\tpath\tchr\tpos\tea\toa\tbeta\teaf\n"
                "study\tstudy.tsv\tCHR\tPOS\tEA\tOA\tBETA\tEAF\n",
                encoding="utf-8",
            )
            output = root / "out"
            summary = convert_manifest(manifest, source_db, target_db, output)
            self.assertEqual(summary, [("study", 2, 1)])

            with (output / "study.lifted.tsv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["LIFTOVER_STATUS"], "mapped_flip")
            self.assertEqual(rows[0]["LIFTOVER_POS"], "80")
            self.assertEqual(rows[0]["LIFTOVER_BETA"], "-0.25")
            self.assertEqual(rows[0]["LIFTOVER_EAF"], "0.8")
            self.assertEqual(rows[0]["LIFTOVER_SOURCE_ASSEMBLY"], "GRCh38")
            self.assertEqual(rows[0]["LIFTOVER_TARGET_DBSNP_RELEASE"], "155")
            self.assertEqual(rows[0]["LIFTOVER_VARIANT_ID"], rows[1]["LIFTOVER_VARIANT_ID"])
            self.assertEqual(rows[1]["LIFTOVER_BETA"], "0.5")

            with (output / "study.variants.tsv").open(encoding="utf-8") as handle:
                variants = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(variants), 1)

    def test_compressed_outputs_are_gzip_tsv_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_vcf, target_vcf = root / "source.vcf", root / "target.vcf"
            row = ["1", "5", "rs1", "A", "C", ".", ".", "."]
            self._write_vcf(source_vcf, [row])
            self._write_vcf(target_vcf, [["1", "8", "rs1", "A", "C", ".", ".", "."]])
            source_db, target_db = root / "s.sqlite", root / "t.sqlite"
            build_index(source_vcf, source_db, "buildA", "1")
            build_index(target_vcf, target_db, "buildB", "1")

            (root / "input.tsv").write_text(
                "CHR\tPOS\tEA\tOA\tBETA\n1\t5\tC\tA\t0.2\n", encoding="utf-8"
            )
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "sample_id\tpath\tchr\tpos\tea\toa\tbeta\n"
                "study\tinput.tsv\tCHR\tPOS\tEA\tOA\tBETA\n",
                encoding="utf-8",
            )
            output = root / "out"

            summary = convert_manifest(
                manifest, source_db, target_db, output, compressed_output=True
            )

            self.assertEqual(summary, [("study", 1, 1)])
            lifted = output / "study.lifted.tsv.gz"
            variants = output / "study.variants.tsv.gz"
            self.assertTrue(lifted.is_file())
            self.assertTrue(variants.is_file())
            self.assertFalse((output / "study.lifted.tsv").exists())
            self.assertFalse((output / "study.variants.tsv").exists())
            with gzip.open(lifted, "rt", encoding="utf-8") as handle:
                lifted_rows = list(csv.DictReader(handle, delimiter="\t"))
            with gzip.open(variants, "rt", encoding="utf-8") as handle:
                variant_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(lifted_rows[0]["LIFTOVER_STATUS"], "mapped")
            self.assertEqual(len(variant_rows), 1)

    def test_target_missing_and_ambiguous_records_are_not_guessed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_vcf, target_vcf = root / "source.vcf", root / "target.vcf"
            self._write_vcf(source_vcf, [["1", "5", "rs1", "A", "C", ".", ".", "."]])
            self._write_vcf(
                target_vcf,
                [
                    ["1", "7", "rs1", "A", "C", ".", ".", "."],
                    ["1", "8", "rs1", "A", "C", ".", ".", "."],
                ],
            )
            source_db, target_db = root / "s.sqlite", root / "t.sqlite"
            build_index(source_vcf, source_db, "buildA", "1")
            build_index(target_vcf, target_db, "buildB", "2")
            with sqlite3.connect(target_db) as connection:
                self.assertEqual(
                    connection.execute("SELECT COUNT(*) FROM variants").fetchone()[0], 2
                )

            input_path, manifest = root / "input.tsv", root / "manifest.tsv"
            input_path.write_text("C\tP\tE\tO\tB\n1\t5\tC\tA\t1\n", encoding="utf-8")
            manifest.write_text(
                "sample_id\tpath\tchr\tpos\tea\toa\tbeta\n"
                "x\tinput.tsv\tC\tP\tE\tO\tB\n",
                encoding="utf-8",
            )
            convert_manifest(manifest, source_db, target_db, root / "out")
            with (root / "out" / "x.lifted.tsv").open(encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["LIFTOVER_STATUS"], "target_ambiguous")

    def test_input_rsid_is_not_used_to_identify_source_variant(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_vcf, target_vcf = root / "source.vcf", root / "target.vcf"
            row = ["1", "5", "rs1", "A", "C", ".", ".", "."]
            self._write_vcf(source_vcf, [row])
            self._write_vcf(target_vcf, [row])
            source_db, target_db = root / "s.sqlite", root / "t.sqlite"
            build_index(source_vcf, source_db, "buildA", "1")
            build_index(target_vcf, target_db, "buildB", "1")

            (root / "input.tsv").write_text(
                "C\tP\tE\tO\tB\tR\n1\t5\tC\tA\t0.2\trs999\n", encoding="utf-8"
            )
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "sample_id\tpath\tchr\tpos\tea\toa\tbeta\trsid\n"
                "x\tinput.tsv\tC\tP\tE\tO\tB\tR\n",
                encoding="utf-8",
            )
            convert_manifest(manifest, source_db, target_db, root / "out")
            with (root / "out" / "x.lifted.tsv").open(encoding="utf-8") as handle:
                result = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(result["LIFTOVER_STATUS"], "mapped")
            self.assertEqual(result["LIFTOVER_RSID"], "rs1")

    def test_each_manifest_sample_can_choose_its_own_output_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_vcf, target_vcf = root / "source.vcf", root / "target.vcf"
            row = ["1", "5", "rs1", "A", "C", ".", ".", "."]
            self._write_vcf(source_vcf, [row])
            self._write_vcf(target_vcf, [["1", "8", "rs1", "A", "C", ".", ".", "."]])
            source_db, target_db = root / "s.sqlite", root / "t.sqlite"
            build_index(source_vcf, source_db, "buildA", "1")
            build_index(target_vcf, target_db, "buildB", "1")

            (root / "input.tsv").write_text(
                "CHR\tPOS\tEA\tOA\tBETA\n1\t5\tC\tA\t0.2\n", encoding="utf-8"
            )
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "sample_id\tpath\toutput_dir\tchr\tpos\tea\toa\tbeta\n"
                "sample_a\tinput.tsv\tresults/a\tCHR\tPOS\tEA\tOA\tBETA\n"
                "sample_b\tinput.tsv\tresults/b\tCHR\tPOS\tEA\tOA\tBETA\n",
                encoding="utf-8",
            )

            summary = convert_manifest(manifest, source_db, target_db, None)

            self.assertEqual(summary, [("sample_a", 1, 1), ("sample_b", 1, 1)])
            for sample_id, relative_dir in (("sample_a", "results/a"), ("sample_b", "results/b")):
                self.assertTrue((root / relative_dir / (sample_id + ".lifted.tsv")).is_file())
                self.assertTrue((root / relative_dir / (sample_id + ".variants.tsv")).is_file())

    def test_source_ref_alt_are_used_to_select_the_source_rsid(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_vcf, target_vcf = root / "source.vcf", root / "target.vcf"
            self._write_vcf(
                source_vcf,
                [
                    ["1", "5", "rs1", "A", "C", ".", ".", "."],
                    ["1", "5", "rs2", "A", "G", ".", ".", "."],
                ],
            )
            self._write_vcf(target_vcf, [["1", "8", "rs1", "A", "C", ".", ".", "."]])
            source_db, target_db = root / "s.sqlite", root / "t.sqlite"
            build_index(source_vcf, source_db, "buildA", "1")
            build_index(target_vcf, target_db, "buildB", "1")

            (root / "input.tsv").write_text(
                "C\tP\tE\tO\tB\tR\tREF\tALT\n1\t5\tC\tA\t0.2\trs2\tA\tC\n",
                encoding="utf-8",
            )
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "sample_id\tpath\tchr\tpos\tea\toa\tbeta\trsid\tsource_ref\tsource_alt\n"
                "x\tinput.tsv\tC\tP\tE\tO\tB\tR\tREF\tALT\n",
                encoding="utf-8",
            )
            convert_manifest(manifest, source_db, target_db, root / "out")
            with (root / "out" / "x.lifted.tsv").open(encoding="utf-8") as handle:
                result = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(result["LIFTOVER_STATUS"], "mapped")
            self.assertEqual(result["LIFTOVER_RSID"], "rs1")

    def test_bigbed_indexes_can_be_persisted_and_use_a_scratch_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_bigbed, target_bigbed = root / "source.bb", root / "target.bb"
            source_bigbed.touch()
            target_bigbed.touch()
            input_path = root / "input.tsv"
            input_path.write_text(
                "CHR\tPOS\tEA\tOA\tBETA\n1\t5\tC\tA\t0.2\n", encoding="utf-8"
            )
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "sample_id\tpath\toutput_dir\tdatabase_dir\tchr\tpos\tea\toa\tbeta\n"
                "sample\tinput.tsv\tresults\tindexes\tCHR\tPOS\tEA\tOA\tBETA\n",
                encoding="utf-8",
            )
            scratch = root / "scratch"

            def fake_run_tool(command, operation):
                result_path = Path(command[2]) if "source-coordinate" in operation else Path(command[-1])
                row = "chr1\t4\t5\trs1\tA\t1\tC\n"
                if "source-coordinate" in operation:
                    result_path.write_text(
                        "chrom\tchromStart\tchromEnd\tname\tref\taltCount\talts\n" + row,
                        encoding="utf-8",
                    )
                else:
                    result_path.write_text(row, encoding="utf-8")

            with patch("dbsnp_liftover.core._run_tool", side_effect=fake_run_tool):
                summary = convert_bigbed_manifest(
                    manifest,
                    source_bigbed,
                    target_bigbed,
                    None,
                    "GRCh38",
                    "155",
                    "GRCh37",
                    "155",
                    database_dir=None,
                    work_dir=scratch,
                    sqlite_mmap_size=0,
                )

                self.assertEqual(summary, [("sample", 1, 1)])
                self.assertTrue((root / "indexes" / "source.sqlite").is_file())
                self.assertTrue((root / "indexes" / "target.sqlite").is_file())
                self.assertEqual(list(scratch.iterdir()), [])
                with (root / "results" / "sample.lifted.tsv").open(encoding="utf-8") as handle:
                    result = next(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(result["LIFTOVER_STATUS"], "mapped")
                with sqlite3.connect(root / "indexes" / "source.sqlite") as connection:
                    self.assertEqual(
                        connection.execute("SELECT rsid FROM variants").fetchone()[0], "rs1"
                    )

                (root / "results" / "sample.lifted.tsv").unlink()
                (root / "results" / "sample.variants.tsv").unlink()
                with self.assertRaisesRegex(LiftoverError, "existing database files"):
                    convert_bigbed_manifest(
                        manifest,
                        source_bigbed,
                        target_bigbed,
                        None,
                        "GRCh38",
                        "155",
                        "GRCh37",
                        "155",
                        database_dir=None,
                        work_dir=scratch,
                        sqlite_mmap_size=0,
                    )


if __name__ == "__main__":
    unittest.main()
