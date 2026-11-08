from __future__ import annotations

import csv
import gzip
import os
import re
import shutil
import sqlite3
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


PRIMARY_CHROMOSOMES = frozenset(
    [str(number) for number in range(1, 23)] + ["X", "Y", "MT"]
)


class LiftoverError(ValueError):
    """Raised when an input, index, or mapping is invalid."""


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    rsid: str
    ref: str
    alt: str


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return path.open("r", encoding="utf-8", errors="replace", newline="")


def open_output_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "wt", compresslevel=1, encoding="utf-8", newline="")
    return path.open("w", encoding="utf-8", newline="")


def normalize_chrom(value: str) -> str:
    chrom = value.strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    if not chrom:
        raise LiftoverError("chromosome is empty")
    if chrom.upper() in {"M", "MT"}:
        return "MT"
    if chrom.isdigit():
        return str(int(chrom))
    return chrom.upper()


def is_primary_chromosome(chrom: str) -> bool:
    return chrom in PRIMARY_CHROMOSOMES


def normalize_allele(value: str) -> str:
    allele = value.strip().upper()
    if not allele or not re.fullmatch(r"[ACGTN]+", allele):
        raise LiftoverError("alleles must contain only A, C, G, T, or N")
    return allele


def reverse_complement(allele: str) -> str:
    return allele.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def input_effect_is_source_alt(ea: str, oa: str, ref: str, alt: str) -> Optional[bool]:
    """Return whether EA is source ALT, accepting a strand-complemented pair."""
    if ea == alt and oa == ref:
        return True
    if ea == ref and oa == alt:
        return False
    comp_ea, comp_oa = reverse_complement(ea), reverse_complement(oa)
    if comp_ea == alt and comp_oa == ref:
        return True
    if comp_ea == ref and comp_oa == alt:
        return False
    return None


def source_alt_maps_to_target_alt(
    source_ref: str, source_alt: str, target_ref: str, target_alt: str
) -> Optional[bool]:
    """Return whether the source ALT is target ALT (rather than target REF)."""
    if source_ref == target_ref and source_alt == target_alt:
        return True
    if source_ref == target_alt and source_alt == target_ref:
        return False
    comp_ref = reverse_complement(source_ref)
    comp_alt = reverse_complement(source_alt)
    if comp_ref == target_ref and comp_alt == target_alt:
        return True
    if comp_ref == target_alt and comp_alt == target_ref:
        return False
    return None


def _rsids(value: str) -> List[str]:
    return sorted(
        {
            token.lower()
            for token in value.split(";")
            if re.fullmatch(r"rs\d+", token, flags=re.IGNORECASE)
        }
    )


def build_index(vcf_path: Path, output_path: Path, assembly: str, release: str) -> int:
    vcf_path = vcf_path.resolve()
    output_path = output_path.resolve()
    if not vcf_path.is_file():
        raise LiftoverError("VCF does not exist: {}".format(vcf_path))
    if output_path.exists():
        raise LiftoverError("refusing to overwrite existing index: {}".format(output_path))
    if not assembly.strip() or not release.strip():
        raise LiftoverError("assembly and release must be non-empty")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(str(output_path))
    inserted = 0
    try:
        connection.executescript(
            """
            CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE variants (
                chrom TEXT NOT NULL,
                pos INTEGER NOT NULL,
                rsid TEXT NOT NULL,
                ref TEXT NOT NULL,
                alt TEXT NOT NULL,
                PRIMARY KEY (chrom, pos, rsid, ref, alt)
            ) WITHOUT ROWID;
            CREATE INDEX variants_by_rsid ON variants (rsid);
            """
        )
        connection.executemany(
            "INSERT INTO metadata(key, value) VALUES (?, ?)",
            [("assembly", assembly.strip()), ("release", release.strip()), ("vcf", str(vcf_path))],
        )
        with open_text(vcf_path) as source:
            saw_header = False
            for line_number, line in enumerate(source, start=1):
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    saw_header = True
                    continue
                if line.startswith("#") or not line.strip():
                    continue
                if not saw_header:
                    raise LiftoverError("VCF data found before #CHROM header at line {}".format(line_number))
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 8:
                    raise LiftoverError("malformed VCF row at line {}".format(line_number))
                rsids = _rsids(fields[2])
                if not rsids:
                    continue
                try:
                    chrom = normalize_chrom(fields[0])
                    pos = int(fields[1])
                    ref = normalize_allele(fields[3])
                except (ValueError, LiftoverError) as exc:
                    raise LiftoverError("invalid VCF variant at line {}: {}".format(line_number, exc))
                if pos < 1:
                    raise LiftoverError("VCF position must be one-based at line {}".format(line_number))
                for alt_value in fields[4].split(","):
                    try:
                        alt = normalize_allele(alt_value)
                    except LiftoverError:
                        continue
                    if alt == ref:
                        continue
                    for rsid in rsids:
                        cursor = connection.execute(
                            "INSERT OR IGNORE INTO variants VALUES (?, ?, ?, ?, ?)",
                            (chrom, pos, rsid, ref, alt),
                        )
                        inserted += cursor.rowcount
        if not saw_header:
            raise LiftoverError("VCF is missing its #CHROM header")
        if inserted == 0:
            raise LiftoverError("no biallelic rsID variants were indexed")
        connection.commit()
    except Exception:
        connection.close()
        output_path.unlink(missing_ok=True)
        raise
    connection.close()
    return inserted


def _create_index(path: Path, assembly: str, release: str, source_path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(str(path))
    connection.executescript(
        """
        CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE variants (
            chrom TEXT NOT NULL,
            pos INTEGER NOT NULL,
            rsid TEXT NOT NULL,
            ref TEXT NOT NULL,
            alt TEXT NOT NULL,
            PRIMARY KEY (chrom, pos, rsid, ref, alt)
        ) WITHOUT ROWID;
        CREATE INDEX variants_by_rsid ON variants (rsid);
        """
    )
    connection.executemany(
        "INSERT INTO metadata(key, value) VALUES (?, ?)",
        [
            ("assembly", assembly),
            ("release", release),
            ("bigbed", str(source_path.resolve())),
        ],
    )
    return connection


def _configure_sqlite_mmap(connection: sqlite3.Connection, mmap_size: int) -> None:
    if mmap_size < 0:
        raise LiftoverError("SQLite mmap size must be zero or a positive byte count")
    connection.execute("PRAGMA mmap_size = {}".format(mmap_size))


def _bigbed_rows(path: Path, has_header: bool):
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if has_header:
            header = next(reader, None)
            expected = ["chrom", "chromStart", "chromEnd", "name", "ref", "altCount", "alts"]
            if header is None or header[: len(expected)] != expected:
                raise LiftoverError("unexpected bigBedToBed TSV header")
        for line_number, fields in enumerate(reader, start=2 if has_header else 1):
            if not fields or len(fields) < 7:
                raise LiftoverError("malformed BigBed row at line {}".format(line_number))
            yield fields


def _insert_bigbed_rows(connection: sqlite3.Connection, path: Path, has_header: bool) -> int:
    inserted = 0
    for fields in _bigbed_rows(path, has_header):
        try:
            chrom = normalize_chrom(fields[0])
            pos = int(fields[1]) + 1
            ref = normalize_allele(fields[4])
        except (ValueError, LiftoverError):
            continue
        if pos < 1:
            continue
        rsids = _rsids(fields[3])
        if not rsids:
            continue
        alts = fields[6].rstrip(",").split(",")
        for alt_value in alts:
            try:
                alt = normalize_allele(alt_value)
            except LiftoverError:
                continue
            if alt == ref:
                continue
            for rsid in rsids:
                cursor = connection.execute(
                    "INSERT OR IGNORE INTO variants VALUES (?, ?, ?, ?, ?)",
                    (chrom, pos, rsid, ref, alt),
                )
                inserted += cursor.rowcount
    return inserted


def _manifest_source_queries(
    manifest_path: Path, manifest_rows: List[Dict[str, str]]
) -> List[Tuple[str, int]]:
    coordinates = set()
    for row in manifest_rows:
        input_path = Path(row["path"])
        if not input_path.is_absolute():
            input_path = manifest_path.parent / input_path
        input_path = input_path.resolve()
        if not input_path.is_file():
            raise LiftoverError("input does not exist: {}".format(input_path))
        handle, delimiter = _open_table(input_path, (row.get("delimiter") or "").strip())
        try:
            reader = csv.DictReader(handle, delimiter=delimiter)
            columns = reader.fieldnames or []
            if len(columns) != len(set(columns)):
                raise LiftoverError("duplicate column names in {}".format(input_path))
            mappings = {}
            for name in (
                "chr",
                "pos",
                "ea",
                "oa",
                "beta",
                "rsid",
                "eaf",
                "se",
                "p",
                "source_ref",
                "source_alt",
            ):
                value = (row.get(name) or "").strip()
                if value:
                    if value not in columns:
                        raise LiftoverError(
                            "column '{}' configured for {} is absent from {}".format(
                                value, row["sample_id"], input_path
                            )
                        )
                    mappings[name] = value
            for input_row in reader:
                key = _canonical_variant_key(
                    *(_value(input_row, mappings[field]) for field in ("chr", "pos", "ea", "oa"))
                )
                if not key[0].startswith("INVALID:"):
                    coordinates.add((key[0], int(key[1])))
        finally:
            handle.close()
    return sorted(coordinates)


def _run_tool(command: List[str], operation: str) -> None:
    try:
        result = subprocess.run(command, check=False, capture_output=True, text=True)
    except OSError as exc:
        raise LiftoverError("could not run {}: {}".format(operation, exc))
    if result.returncode:
        detail = result.stderr.strip() or result.stdout.strip() or "exit code {}".format(result.returncode)
        raise LiftoverError("{} failed: {}".format(operation, detail))


def convert_bigbed_manifest(
    manifest_path: Path,
    source_bigbed: Path,
    target_bigbed: Path,
    output_dir: Optional[Path],
    source_assembly: str,
    source_release: str,
    target_assembly: str,
    target_release: str,
    bigbed_to_bed: str = "bigBedToBed",
    bigbed_named_items: str = "bigBedNamedItems",
    database_dir: Optional[Path] = None,
    work_dir: Optional[Path] = None,
    sqlite_mmap_size: int = 256 * 1024 * 1024,
    compressed_output: bool = False,
    prefer_unique_primary: bool = False,
    allow_partial_source_rsids: bool = False,
) -> List[Tuple[str, int, int]]:
    """Query BigBed by input coordinates and by rsID without materializing either file."""
    manifest_path = manifest_path.resolve()
    source_bigbed = source_bigbed.resolve()
    target_bigbed = target_bigbed.resolve()
    if not manifest_path.is_file():
        raise LiftoverError("manifest does not exist: {}".format(manifest_path))
    if not source_bigbed.is_file() or not target_bigbed.is_file():
        raise LiftoverError("both source and target BigBed files must exist")
    if source_bigbed == target_bigbed:
        raise LiftoverError("source and target BigBed files must be different")
    labels = (source_assembly, source_release, target_assembly, target_release)
    if not all(value.strip() for value in labels):
        raise LiftoverError("source/target assembly and dbSNP release labels are required")

    manifest_rows = _read_manifest(manifest_path)
    coordinates = _manifest_source_queries(manifest_path, manifest_rows)
    sample_output_paths = _sample_output_paths(
        manifest_path, manifest_rows, output_dir, compressed_output
    )
    database_output_paths = _database_output_paths(
        manifest_path, manifest_rows, database_dir
    )
    if database_output_paths:
        sample_paths = {path for paths in sample_output_paths for path in paths}
        if sample_paths.intersection(database_output_paths):
            raise LiftoverError("database output paths conflict with sample output files")
        existing = [str(path) for path in database_output_paths if path.exists()]
        if existing:
            raise LiftoverError(
                "refusing to overwrite existing database files: {}".format(", ".join(existing))
            )
    if sqlite_mmap_size < 0:
        raise LiftoverError("SQLite mmap size must be zero or a positive byte count")

    resolved_work_dir = None
    if work_dir is not None:
        resolved_work_dir = work_dir.resolve()
        resolved_work_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(
        prefix="dbsnp_liftover_", dir=str(resolved_work_dir) if resolved_work_dir else None
    ) as temporary:
        temp = Path(temporary)
        source_db_path, target_db_path = temp / "source.sqlite", temp / "target.sqlite"
        source_db = _create_index(
            source_db_path, source_assembly.strip(), source_release.strip(), source_bigbed
        )
        try:
            for offset in range(0, len(coordinates), 25000):
                chunk = coordinates[offset : offset + 25000]
                bed_path, output_path = temp / "source.bed", temp / "source.tsv"
                with bed_path.open("w", encoding="utf-8") as bed:
                    for chrom, pos in chunk:
                        ucsc_chrom = "chrM" if chrom == "MT" else "chr" + chrom
                        bed.write("{}\t{}\t{}\n".format(ucsc_chrom, pos - 1, pos))
                _run_tool(
                    [
                        bigbed_to_bed,
                        str(source_bigbed),
                        str(output_path),
                        "-bed={}".format(bed_path),
                        "-tsv",
                    ],
                    "bigBedToBed source-coordinate lookup",
                )
                _insert_bigbed_rows(source_db, output_path, has_header=True)
                print(
                    "Source BigBed coordinate lookup: {} / {}".format(
                        min(offset + len(chunk), len(coordinates)), len(coordinates)
                    ),
                    flush=True,
                )
            source_db.commit()
        finally:
            source_db.close()

        with sqlite3.connect(str(source_db_path)) as lookup:
            _configure_sqlite_mmap(lookup, sqlite_mmap_size)
            rsids = [
                item[0]
                for item in lookup.execute("SELECT DISTINCT rsid FROM variants ORDER BY rsid")
            ]

        target_db = _create_index(
            target_db_path, target_assembly.strip(), target_release.strip(), target_bigbed
        )
        try:
            for offset in range(0, len(rsids), 25000):
                chunk = rsids[offset : offset + 25000]
                names_path, output_path = temp / "rsids.txt", temp / "target.bed"
                names_path.write_text("\n".join(chunk) + "\n", encoding="utf-8")
                _run_tool(
                    [
                        bigbed_named_items,
                        "-nameFile",
                        str(target_bigbed),
                        str(names_path),
                        str(output_path),
                    ],
                    "bigBedNamedItems target-rsID lookup",
                )
                _insert_bigbed_rows(target_db, output_path, has_header=False)
                print(
                    "Target BigBed rsID lookup: {} / {}".format(
                        min(offset + len(chunk), len(rsids)), len(rsids)
                    ),
                    flush=True,
                )
            target_db.commit()
        finally:
            target_db.close()

        if database_output_paths is not None:
            _persist_bigbed_indexes(source_db_path, target_db_path, database_output_paths)

        return convert_manifest(
            manifest_path,
            source_db_path,
            target_db_path,
            output_dir,
            sqlite_mmap_size=sqlite_mmap_size,
            compressed_output=compressed_output,
            prefer_unique_primary=prefer_unique_primary,
            allow_partial_source_rsids=allow_partial_source_rsids,
        )


def _read_metadata(connection: sqlite3.Connection) -> Dict[str, str]:
    try:
        return dict(connection.execute("SELECT key, value FROM metadata"))
    except sqlite3.Error as exc:
        raise LiftoverError("not a dbSNP liftover index: {}".format(exc))


def _read_manifest(path: Path) -> List[Dict[str, str]]:
    with open_text(path) as source:
        reader = csv.DictReader(source, delimiter="\t")
        required = {"sample_id", "path", "chr", "pos", "ea", "oa", "beta"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise LiftoverError("manifest is missing columns: {}".format(", ".join(sorted(missing))))
        rows = list(reader)
    if not rows:
        raise LiftoverError("manifest contains no input files")
    seen = set()
    for row in rows:
        sample_id = (row.get("sample_id") or "").strip()
        if not sample_id or not re.fullmatch(r"[A-Za-z0-9._-]+", sample_id):
            raise LiftoverError("sample_id must use only letters, numbers, '.', '_' or '-'")
        if sample_id in seen:
            raise LiftoverError("duplicate sample_id in manifest: {}".format(sample_id))
        seen.add(sample_id)
        for field in ("path", "chr", "pos", "ea", "oa", "beta"):
            if not (row.get(field) or "").strip():
                raise LiftoverError("manifest field '{}' is empty for {}".format(field, sample_id))
    return rows


def _sample_output_paths(
    manifest_path: Path,
    manifest_rows: List[Dict[str, str]],
    default_output_dir: Optional[Path],
    compressed_output: bool = False,
) -> List[Tuple[Path, Path, Path]]:
    resolved_default = default_output_dir.resolve() if default_output_dir is not None else None
    suffix = ".tsv.gz" if compressed_output else ".tsv"
    paths = []
    for row in manifest_rows:
        configured_dir = (row.get("output_dir") or "").strip()
        if configured_dir:
            sample_dir = Path(configured_dir)
            if not sample_dir.is_absolute():
                sample_dir = manifest_path.parent / sample_dir
            sample_dir = sample_dir.resolve()
        elif resolved_default is not None:
            sample_dir = resolved_default
        else:
            raise LiftoverError(
                "output_dir must be set in the manifest for sample '{}' or supplied with --output-dir".format(
                    row["sample_id"]
                )
            )
        paths.append(
            (
                sample_dir / "{}.lifted{}".format(row["sample_id"], suffix),
                sample_dir / "{}.variants{}".format(row["sample_id"], suffix),
                sample_dir / "{}.unmapped{}".format(row["sample_id"], suffix),
            )
        )

    all_paths = [path for sample_paths in paths for path in sample_paths]
    normalized = [str(path) for path in all_paths]
    if len(normalized) != len(set(normalized)):
        raise LiftoverError("multiple manifest rows resolve to the same output file")
    existing = [str(path) for path in all_paths if path.exists()]
    if existing:
        raise LiftoverError("refusing to overwrite output files: {}".format(", ".join(existing)))
    return paths


def _database_output_paths(
    manifest_path: Path,
    manifest_rows: List[Dict[str, str]],
    default_database_dir: Optional[Path],
) -> Optional[Tuple[Path, Path]]:
    configured_values = [(row.get("database_dir") or "").strip() for row in manifest_rows]
    if not any(configured_values) and default_database_dir is None:
        return None

    configured_dirs = []
    for value in configured_values:
        if value:
            path = Path(value)
            if not path.is_absolute():
                path = manifest_path.parent / path
            configured_dirs.append(path.resolve())
        elif default_database_dir is not None:
            configured_dirs.append(default_database_dir.resolve())
        else:
            raise LiftoverError(
                "database_dir must be set for every sample or supplied with --database-dir"
            )

    if not configured_dirs:
        return None
    unique_dirs = set(configured_dirs)
    if len(unique_dirs) != 1:
        raise LiftoverError(
            "database_dir must resolve to one shared directory for all samples in this manifest"
        )
    database_dir = configured_dirs[0]
    return database_dir / "source.sqlite", database_dir / "target.sqlite"


def _persist_bigbed_indexes(
    source_db_path: Path, target_db_path: Path, output_paths: Tuple[Path, Path]
) -> None:
    source_output, target_output = output_paths
    database_dir = source_output.parent
    database_dir.mkdir(parents=True, exist_ok=True)
    if source_output.exists() or target_output.exists():
        existing = [str(path) for path in output_paths if path.exists()]
        raise LiftoverError(
            "refusing to overwrite existing database files: {}".format(", ".join(existing))
        )

    with tempfile.TemporaryDirectory(prefix=".dbsnp_liftover_index_", dir=str(database_dir)) as staging:
        staged_source, staged_target = Path(staging) / "source.sqlite", Path(staging) / "target.sqlite"
        shutil.copy2(source_db_path, staged_source)
        shutil.copy2(target_db_path, staged_target)
        os.link(str(staged_source), str(source_output))
        try:
            os.link(str(staged_target), str(target_output))
        except OSError:
            source_output.unlink(missing_ok=True)
            raise


def _open_table(path: Path, delimiter_override: str = ""):
    handle = open_text(path)
    first_line = handle.readline()
    if not first_line:
        handle.close()
        raise LiftoverError("input table is empty: {}".format(path))
    if delimiter_override:
        delimiter_names = {"tab": "\t", "comma": ",", "semicolon": ";"}
        delimiter = delimiter_names.get(delimiter_override.lower(), delimiter_override)
        if delimiter not in {"\t", ",", ";"}:
            handle.close()
            raise LiftoverError("delimiter must be tab, comma, semicolon, or one of those characters")
    else:
        try:
            delimiter = csv.Sniffer().sniff(first_line, delimiters="\t,;").delimiter
        except csv.Error:
            delimiter = "\t"
    handle.seek(0)
    return handle, delimiter


def _value(row: Dict[str, Optional[str]], field: str) -> str:
    return (row.get(field) or "").strip()


def _canonical_variant_key(
    chrom_value: str, pos_value: str, ea_value: str, oa_value: str
) -> Tuple[str, str, str, str]:
    try:
        chrom = normalize_chrom(chrom_value)
        pos = int(pos_value)
        ea, oa = normalize_allele(ea_value), normalize_allele(oa_value)
        if pos < 1:
            raise LiftoverError("POS must be one-based and positive")
        return chrom, str(pos), ea, oa
    except (ValueError, LiftoverError):
        return "INVALID:" + chrom_value.strip(), pos_value.strip(), ea_value.strip(), oa_value.strip()


def _find_source_variants(
    connection: sqlite3.Connection,
    chrom: str,
    pos: int,
    ea: str,
    oa: str,
    source_ref: str,
    source_alt: str,
) -> List[Tuple[Variant, bool]]:
    records = connection.execute(
        "SELECT rsid, ref, alt FROM variants WHERE chrom = ? AND pos = ?",
        (chrom, pos),
    )
    matches = []
    for rsid, ref, alt in records:
        if source_ref and source_alt:
            if not (
                (source_ref == ref and source_alt == alt)
                or (
                    reverse_complement(source_ref) == ref
                    and reverse_complement(source_alt) == alt
                )
            ):
                continue
        effect_is_alt = input_effect_is_source_alt(ea, oa, ref, alt)
        if effect_is_alt is not None:
            matches.append((rsid, ref, alt, effect_is_alt))
    unique = sorted(set(matches))
    return [
        (Variant(chrom, pos, rsid, ref, alt), effect_is_alt)
        for rsid, ref, alt, effect_is_alt in unique
    ]


def _find_target_variant(
    connection: sqlite3.Connection,
    source: Variant,
    prefer_unique_primary: bool = False,
) -> Tuple[str, Optional[Variant], Optional[bool]]:
    rows = connection.execute(
        "SELECT chrom, pos, ref, alt FROM variants WHERE rsid = ?",
        (source.rsid,),
    )
    candidates = set()
    for chrom, pos, ref, alt in rows:
        same_direction = source_alt_maps_to_target_alt(source.ref, source.alt, ref, alt)
        if same_direction is not None:
            candidates.add((chrom, pos, ref, alt, same_direction))
    if not candidates:
        return "target_variant_not_found", None, None
    if prefer_unique_primary and len(candidates) > 1:
        primary_candidates = {
            candidate for candidate in candidates if is_primary_chromosome(candidate[0])
        }
        if len(primary_candidates) == 1:
            candidates = primary_candidates
    if len(candidates) != 1:
        return "target_ambiguous", None, None
    chrom, pos, ref, alt, same_direction = next(iter(candidates))
    return "mapped", Variant(chrom, pos, source.rsid, ref, alt), same_direction


OUTPUT_FIELDS = [
    "LIFTOVER_STATUS",
    "LIFTOVER_VARIANT_ID",
    "LIFTOVER_RSID",
    "LIFTOVER_CHR",
    "LIFTOVER_POS",
    "LIFTOVER_REF",
    "LIFTOVER_ALT",
    "LIFTOVER_EA",
    "LIFTOVER_OA",
    "LIFTOVER_BETA",
    "LIFTOVER_EAF",
    "LIFTOVER_NOTE",
    "LIFTOVER_SOURCE_ASSEMBLY",
    "LIFTOVER_SOURCE_DBSNP_RELEASE",
    "LIFTOVER_TARGET_ASSEMBLY",
    "LIFTOVER_TARGET_DBSNP_RELEASE",
]

VARIANT_FIELDS = [
    "LIFTOVER_VARIANT_ID",
    "CHR",
    "POS",
    "EA",
    "OA",
    "STATUS",
    "RSID",
    "TARGET_CHR",
    "TARGET_POS",
    "TARGET_REF",
    "TARGET_ALT",
    "BETA_FLIPPED",
    "SOURCE_ASSEMBLY",
    "SOURCE_DBSNP_RELEASE",
    "TARGET_ASSEMBLY",
    "TARGET_DBSNP_RELEASE",
]


def _format_number(value: float) -> str:
    return format(value, ".15g")


def _map_unique_variant(
    source_db: sqlite3.Connection,
    target_db: sqlite3.Connection,
    chrom_value: str,
    pos_value: str,
    ea_value: str,
    oa_value: str,
    source_ref_value: str = "",
    source_alt_value: str = "",
    prefer_unique_primary: bool = False,
    allow_partial_source_rsids: bool = False,
) -> Dict[str, str]:
    result = {field: "" for field in OUTPUT_FIELDS}
    try:
        chrom = normalize_chrom(chrom_value)
        pos = int(pos_value)
        ea, oa = normalize_allele(ea_value), normalize_allele(oa_value)
        source_ref = normalize_allele(source_ref_value) if source_ref_value else ""
        source_alt = normalize_allele(source_alt_value) if source_alt_value else ""
        if pos < 1:
            raise LiftoverError("POS must be one-based and positive")
        if bool(source_ref) != bool(source_alt):
            raise LiftoverError("source REF and ALT must either both be provided or both be empty")
    except (ValueError, LiftoverError) as exc:
        result["LIFTOVER_STATUS"] = "invalid_input"
        result["LIFTOVER_NOTE"] = str(exc)
        return result
    result["LIFTOVER_EA"], result["LIFTOVER_OA"] = ea, oa
    source_candidates = _find_source_variants(
        source_db, chrom, pos, ea, oa, source_ref, source_alt
    )
    if not source_candidates:
        result["LIFTOVER_STATUS"] = "source_variant_not_found"
        return result
    if len(source_candidates) == 1:
        source, effect_is_source_alt = source_candidates[0]
        result["LIFTOVER_RSID"] = source.rsid
        target_status, target, same_direction = _find_target_variant(
            target_db, source, prefer_unique_primary
        )
        result["LIFTOVER_STATUS"] = target_status
        if target is None or same_direction is None:
            return result
    else:
        candidate_rsids = sorted({source.rsid for source, _ in source_candidates})
        if not allow_partial_source_rsids:
            result["LIFTOVER_STATUS"] = "source_ambiguous"
            result["LIFTOVER_NOTE"] = (
                "multiple source rsIDs matched; enable partial source rsID resolution "
                "to compare their target mappings: {}".format(";".join(candidate_rsids))
            )
            return result

        resolved = []
        unmapped_source_rsids = []
        ambiguous_target_rsids = []
        for source, effect_is_alt in source_candidates:
            target_status, target, same_direction = _find_target_variant(
                target_db, source, prefer_unique_primary
            )
            if target_status == "target_variant_not_found":
                unmapped_source_rsids.append(source.rsid)
                continue
            if target_status != "mapped" or target is None or same_direction is None:
                ambiguous_target_rsids.append(source.rsid)
                continue
            flip = effect_is_alt != same_direction
            resolved.append(
                (
                    source.rsid,
                    target,
                    effect_is_alt,
                    same_direction,
                    flip,
                )
            )
        signatures = {
            (
                target.chrom,
                target.pos,
                target.ref,
                target.alt,
                same_direction,
                flip,
            )
            for _, target, _, same_direction, flip in resolved
        }
        unresolved_candidates = len(resolved) != len(source_candidates)
        can_ignore_missing = (
            bool(resolved)
            and not ambiguous_target_rsids
            and len(resolved) + len(unmapped_source_rsids) == len(source_candidates)
        )
        if (
            not resolved
            or len(signatures) != 1
            or (unresolved_candidates and not can_ignore_missing)
        ):
            result["LIFTOVER_STATUS"] = "source_ambiguous"
            result["LIFTOVER_NOTE"] = (
                "candidate source rsIDs do not resolve to one agreed target variant "
                "and effect direction: {}".format(";".join(candidate_rsids))
            )
            return result
        _, target, effect_is_source_alt, same_direction, _ = resolved[0]
        result["LIFTOVER_RSID"] = ";".join(candidate_rsids)
        note = "candidate source rsIDs resolve consistently: {}".format(
            ";".join(candidate_rsids)
        )
        if unmapped_source_rsids:
            note += "; no compatible target mapping for: {}".format(
                ";".join(sorted(set(unmapped_source_rsids)))
            )
        result["LIFTOVER_NOTE"] = note

    flip = effect_is_source_alt != same_direction
    result.update(
        {
            "LIFTOVER_CHR": target.chrom,
            "LIFTOVER_POS": str(target.pos),
            "LIFTOVER_REF": target.ref,
            "LIFTOVER_ALT": target.alt,
            "LIFTOVER_EA": target.alt,
            "LIFTOVER_OA": target.ref,
        }
    )
    result["LIFTOVER_BETA"] = "FLIP_REQUIRED" if flip else "UNCHANGED"
    result["LIFTOVER_EAF"] = "FLIP_REQUIRED" if flip else "UNCHANGED"
    result["LIFTOVER_STATUS"] = "mapped_flip" if flip else "mapped"
    return result


def convert_manifest(
    manifest_path: Path,
    source_index_path: Path,
    target_index_path: Path,
    output_dir: Optional[Path],
    sqlite_mmap_size: int = 256 * 1024 * 1024,
    compressed_output: bool = False,
    prefer_unique_primary: bool = False,
    allow_partial_source_rsids: bool = False,
) -> List[Tuple[str, int, int]]:
    manifest_path = manifest_path.resolve()
    source_index_path = source_index_path.resolve()
    target_index_path = target_index_path.resolve()
    if not manifest_path.is_file():
        raise LiftoverError("manifest does not exist: {}".format(manifest_path))
    if not source_index_path.is_file() or not target_index_path.is_file():
        raise LiftoverError("both source and target index files must exist")
    if source_index_path == target_index_path:
        raise LiftoverError("source and target indexes must be different files")
    if sqlite_mmap_size < 0:
        raise LiftoverError("SQLite mmap size must be zero or a positive byte count")

    manifest_rows = _read_manifest(manifest_path)
    sample_output_paths = _sample_output_paths(
        manifest_path, manifest_rows, output_dir, compressed_output
    )
    input_paths = []
    for row in manifest_rows:
        input_path = Path(row["path"])
        if not input_path.is_absolute():
            input_path = manifest_path.parent / input_path
        input_path = input_path.resolve()
        if not input_path.is_file():
            raise LiftoverError("input does not exist: {}".format(input_path))
        input_paths.append(input_path)

    summaries = []
    with sqlite3.connect(str(source_index_path)) as source_db, sqlite3.connect(
        str(target_index_path)
    ) as target_db:
        _configure_sqlite_mmap(source_db, sqlite_mmap_size)
        _configure_sqlite_mmap(target_db, sqlite_mmap_size)
        source_meta, target_meta = _read_metadata(source_db), _read_metadata(target_db)
        if not all(
            metadata.get(field)
            for metadata in (source_meta, target_meta)
            for field in ("assembly", "release")
        ):
            raise LiftoverError("source and target indexes must declare their assembly and release")
        for row, input_path, output_paths in zip(
            manifest_rows, input_paths, sample_output_paths
        ):
            lifted_path, variants_path, unmapped_path = output_paths
            lifted_path.parent.mkdir(parents=True, exist_ok=True)
            variants_path.parent.mkdir(parents=True, exist_ok=True)
            unmapped_path.parent.mkdir(parents=True, exist_ok=True)
            handle, delimiter = _open_table(input_path, (row.get("delimiter") or "").strip())
            try:
                reader = csv.DictReader(handle, delimiter=delimiter)
                columns = reader.fieldnames or []
                if len(columns) != len(set(columns)):
                    raise LiftoverError("duplicate column names in {}".format(input_path))
                reserved_columns = set(columns) & set(OUTPUT_FIELDS)
                if reserved_columns:
                    raise LiftoverError(
                        "input contains reserved output columns: {}".format(
                            ", ".join(sorted(reserved_columns))
                        )
                    )
                mappings = {}
                for name in (
                    "chr",
                    "pos",
                    "ea",
                    "oa",
                    "beta",
                    "rsid",
                    "eaf",
                    "se",
                    "p",
                    "source_ref",
                    "source_alt",
                ):
                    value = (row.get(name) or "").strip()
                    if value:
                        if value not in columns:
                            raise LiftoverError(
                                "column '{}' configured for {} is absent from {}".format(
                                    value, row["sample_id"], input_path
                                )
                            )
                        mappings[name] = value
                source_rows: Dict[Tuple[str, str, str, str], Dict[str, str]] = {}
                source_alleles: Dict[Tuple[str, str, str, str], set] = {}
                row_count = 0
                for input_row in reader:
                    row_count += 1
                    key = _canonical_variant_key(
                        *(_value(input_row, mappings[field]) for field in ("chr", "pos", "ea", "oa"))
                    )
                    source_ref = _value(input_row, mappings["source_ref"]).upper() if "source_ref" in mappings else ""
                    source_alt = _value(input_row, mappings["source_alt"]).upper() if "source_alt" in mappings else ""
                    if source_ref or source_alt:
                        source_alleles.setdefault(key, set()).add((source_ref, source_alt))
                    if key not in source_rows:
                        source_rows[key] = {}
                for key in source_rows:
                    allele_values = source_alleles.get(key, set())
                    if len(allele_values) > 1:
                        source_rows[key] = {field: "" for field in OUTPUT_FIELDS}
                        source_rows[key]["LIFTOVER_STATUS"] = "source_allele_conflict"
                    else:
                        source_ref, source_alt = next(iter(allele_values), ("", ""))
                        source_rows[key] = _map_unique_variant(
                            source_db,
                            target_db,
                            key[0],
                            key[1],
                            key[2],
                            key[3],
                            source_ref,
                            source_alt,
                            prefer_unique_primary,
                            allow_partial_source_rsids,
                        )
                    source_rows[key].update(
                        {
                            "LIFTOVER_SOURCE_ASSEMBLY": source_meta["assembly"],
                            "LIFTOVER_SOURCE_DBSNP_RELEASE": source_meta["release"],
                            "LIFTOVER_TARGET_ASSEMBLY": target_meta["assembly"],
                            "LIFTOVER_TARGET_DBSNP_RELEASE": target_meta["release"],
                        }
                    )
                if len(source_rows) == 0 and row_count == 0:
                    raise LiftoverError("input contains a header but no data rows: {}".format(input_path))

                output_columns = list(columns) + [field for field in OUTPUT_FIELDS if field not in columns]
                variant_index = {
                    key: "V{:08d}".format(number)
                    for number, key in enumerate(source_rows, start=1)
                }
                with open_output_text(lifted_path) as lifted_handle, open_output_text(
                    unmapped_path
                ) as unmapped_handle:
                    writer = csv.DictWriter(
                        lifted_handle, fieldnames=output_columns, delimiter="\t", lineterminator="\n"
                    )
                    unmapped_writer = csv.DictWriter(
                        unmapped_handle,
                        fieldnames=output_columns,
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    unmapped_writer.writeheader()
                    handle.seek(0)
                    reader = csv.DictReader(handle, delimiter=delimiter)
                    for input_row in reader:
                        key = _canonical_variant_key(
                            *(_value(input_row, mappings[field]) for field in ("chr", "pos", "ea", "oa"))
                        )
                        mapped = dict(source_rows[key])
                        mapped["LIFTOVER_VARIANT_ID"] = variant_index[key]
                        status = mapped["LIFTOVER_STATUS"]
                        if status.startswith("mapped"):
                            flip = status == "mapped_flip"
                            beta = _value(input_row, mappings["beta"])
                            errors = []
                            if flip:
                                try:
                                    mapped["LIFTOVER_BETA"] = _format_number(-float(beta))
                                except ValueError:
                                    mapped["LIFTOVER_BETA"] = ""
                                    errors.append("beta")
                            else:
                                mapped["LIFTOVER_BETA"] = beta
                            if "eaf" in mappings:
                                eaf = _value(input_row, mappings["eaf"])
                                if flip and eaf:
                                    try:
                                        frequency = float(eaf)
                                        if not 0 <= frequency <= 1:
                                            raise ValueError
                                        mapped["LIFTOVER_EAF"] = _format_number(1 - frequency)
                                    except ValueError:
                                        mapped["LIFTOVER_EAF"] = ""
                                        errors.append("eaf")
                                else:
                                    mapped["LIFTOVER_EAF"] = eaf
                            else:
                                mapped["LIFTOVER_EAF"] = ""
                            if errors:
                                mapped["LIFTOVER_STATUS"] = "mapped_invalid_{}".format("_".join(errors))
                                mapped["LIFTOVER_NOTE"] = "could not safely flip {}".format(", ".join(errors))
                        else:
                            mapped["LIFTOVER_BETA"] = ""
                            mapped["LIFTOVER_EAF"] = ""
                        output_row = {**input_row, **mapped}
                        writer.writerow(output_row)
                        if not mapped["LIFTOVER_STATUS"].startswith("mapped"):
                            unmapped_writer.writerow(output_row)

                with open_output_text(variants_path) as variants_handle:
                    writer = csv.DictWriter(
                        variants_handle, fieldnames=VARIANT_FIELDS, delimiter="\t", lineterminator="\n"
                    )
                    writer.writeheader()
                    for key, mapped in source_rows.items():
                        status = mapped["LIFTOVER_STATUS"]
                        flip = status == "mapped_flip"
                        writer.writerow(
                            {
                                "LIFTOVER_VARIANT_ID": variant_index[key],
                                "CHR": key[0].removeprefix("INVALID:"),
                                "POS": key[1],
                                "EA": key[2].upper(),
                                "OA": key[3].upper(),
                                "STATUS": status,
                                "RSID": mapped["LIFTOVER_RSID"],
                                "TARGET_CHR": mapped["LIFTOVER_CHR"],
                                "TARGET_POS": mapped["LIFTOVER_POS"],
                                "TARGET_REF": mapped["LIFTOVER_REF"],
                                "TARGET_ALT": mapped["LIFTOVER_ALT"],
                                "BETA_FLIPPED": str(flip).lower(),
                                "SOURCE_ASSEMBLY": source_meta["assembly"],
                                "SOURCE_DBSNP_RELEASE": source_meta["release"],
                                "TARGET_ASSEMBLY": target_meta["assembly"],
                                "TARGET_DBSNP_RELEASE": target_meta["release"],
                            }
                        )
            finally:
                handle.close()
            mapped_count = sum(
                item["LIFTOVER_STATUS"].startswith("mapped") for item in source_rows.values()
            )
            summaries.append((row["sample_id"], row_count, mapped_count))
    return summaries
