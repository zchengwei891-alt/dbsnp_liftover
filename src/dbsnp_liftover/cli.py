from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

from .core import LiftoverError, build_index, convert_bigbed_manifest, convert_manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dbsnp-liftover",
        description="Map GWAS variants between dbSNP assembly snapshots.",
    )
    commands = parser.add_subparsers(dest="command", required=True)

    index = commands.add_parser("index-vcf", help="build a local SQLite index from a dbSNP VCF")
    index.add_argument("--vcf", required=True, type=Path)
    index.add_argument("--assembly", required=True, help="assembly label, e.g. GRCh38")
    index.add_argument("--release", required=True, help="dbSNP release, e.g. 155")
    index.add_argument("--output", required=True, type=Path)

    convert = commands.add_parser("convert", help="convert every input listed in a TSV manifest")
    convert.add_argument("--manifest", required=True, type=Path)
    convert.add_argument("--source-index", required=True, type=Path)
    convert.add_argument("--target-index", required=True, type=Path)
    convert.add_argument(
        "--output-dir",
        type=Path,
        help="default output directory for manifest rows without output_dir",
    )
    convert.add_argument(
        "--sqlite-mmap-mib",
        type=int,
        default=256,
        help="SQLite memory-mapped read window in MiB; use 0 to disable (default: 256)",
    )
    convert.add_argument(
        "--compress-output",
        action="store_true",
        help="write .tsv.gz outputs instead of plain TSV files",
    )
    convert.add_argument(
        "--prefer-unique-primary",
        action="store_true",
        help="resolve multiple target matches when exactly one compatible candidate is on a primary chromosome",
    )
    convert.add_argument(
        "--allow-partial-source-rsids",
        action="store_true",
        help="allow consistent mapped source rsIDs when other candidates have no compatible target mapping",
    )

    bigbed = commands.add_parser(
        "convert-bigbed",
        help="query source/target BigBed files directly without full-file conversion",
    )
    bigbed.add_argument("--manifest", required=True, type=Path)
    bigbed.add_argument("--source-bigbed", required=True, type=Path)
    bigbed.add_argument("--target-bigbed", required=True, type=Path)
    bigbed.add_argument("--source-assembly", required=True)
    bigbed.add_argument("--source-release", required=True)
    bigbed.add_argument("--target-assembly", required=True)
    bigbed.add_argument("--target-release", required=True)
    bigbed.add_argument(
        "--output-dir",
        type=Path,
        help="default output directory for manifest rows without output_dir",
    )
    bigbed.add_argument(
        "--database-dir",
        type=Path,
        help="directory for persistent source.sqlite and target.sqlite lookup indexes",
    )
    bigbed.add_argument(
        "--work-dir",
        type=Path,
        help="directory for temporary query files and working SQLite indexes (e.g. /dev/shm)",
    )
    bigbed.add_argument(
        "--sqlite-mmap-mib",
        type=int,
        default=256,
        help="SQLite memory-mapped read window in MiB; use 0 to disable (default: 256)",
    )
    bigbed.add_argument(
        "--compress-output",
        action="store_true",
        help="write .tsv.gz outputs instead of plain TSV files",
    )
    bigbed.add_argument(
        "--prefer-unique-primary",
        action="store_true",
        help="resolve multiple target matches when exactly one compatible candidate is on a primary chromosome",
    )
    bigbed.add_argument(
        "--allow-partial-source-rsids",
        action="store_true",
        help="allow consistent mapped source rsIDs when other candidates have no compatible target mapping",
    )
    bigbed.add_argument("--bigbed-to-bed", default="bigBedToBed")
    bigbed.add_argument("--bigbed-named-items", default="bigBedNamedItems")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "index-vcf":
            count = build_index(args.vcf, args.output, args.assembly, args.release)
            print(
                "Indexed {} variants (dbSNP {}; {}) -> {}".format(
                    count, args.release, args.assembly, args.output
                )
            )
        elif args.command == "convert":
            summaries = convert_manifest(
                args.manifest,
                args.source_index,
                args.target_index,
                args.output_dir,
                sqlite_mmap_size=args.sqlite_mmap_mib * 1024 * 1024,
                compressed_output=args.compress_output,
                prefer_unique_primary=args.prefer_unique_primary,
                allow_partial_source_rsids=args.allow_partial_source_rsids,
            )
        elif args.command == "convert-bigbed":
            summaries = convert_bigbed_manifest(
                args.manifest,
                args.source_bigbed,
                args.target_bigbed,
                args.output_dir,
                args.source_assembly,
                args.source_release,
                args.target_assembly,
                args.target_release,
                args.bigbed_to_bed,
                args.bigbed_named_items,
                args.database_dir,
                args.work_dir,
                args.sqlite_mmap_mib * 1024 * 1024,
                args.compress_output,
                args.prefer_unique_primary,
                args.allow_partial_source_rsids,
            )
        else:
            raise LiftoverError("unknown command: {}".format(args.command))
        if args.command != "index-vcf":
            for sample_id, rows, mapped in summaries:
                print("{}: rows={}, mapped_unique_variants={}".format(sample_id, rows, mapped))
    except (LiftoverError, OSError) as exc:
        parser.exit(2, "error: {}\n".format(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
