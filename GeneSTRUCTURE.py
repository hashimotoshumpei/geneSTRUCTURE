import csv
import argparse
import os
from parse_utils import (
    parse_gff_for_transcript,
    parse_gff_for_region,
    parse_deletions,
    parse_insertions,
    parse_snps,
    parse_domains,
    print_features_as_gff
)
from draw_utils import draw_gene_structure, draw_region_gene_structures
from welcome_message import print_welcome_message


def parse_args():
    parser = argparse.ArgumentParser(
        description="Draw gene structure SVG from GFF/GTF and CSV input"
    )

    parser.add_argument(
        "--gff", "--gtf",
        dest="gff_file",
        required=True,
        help="Path to GFF or GTF file"
    )

    parser.add_argument(
        "--input", "-i",
        dest="input_csv",
        default=None,
        help="Input CSV file describing gene features (transcript mode)"
    )

    parser.add_argument(
        "--output", "-o",
        dest="output_prefix",
        required=True,
        help="Output file prefix or directory (path allowed)"
    )

    # Region mode arguments
    parser.add_argument(
        "--chr",
        dest="chromosome",
        default=None,
        help="Chromosome/SeqID for region mode"
    )

    parser.add_argument(
        "--start",
        type=int,
        default=None,
        help="Region start position (1-based)"
    )

    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help="Region end position (1-based)"
    )

    args = parser.parse_args()
    validate_args(parser, args)
    return args


def validate_args(parser, args):
    """Validate argument combinations"""
    region_args = [args.chromosome, args.start, args.end]
    has_region = any(region_args)
    has_input = args.input_csv is not None

    # Both modes specified
    if has_region and has_input:
        parser.error("Cannot use --input with region mode (--chr, --start, --end)")

    # Region mode requires all three arguments
    if has_region and not all(region_args):
        parser.error("Region mode requires all of --chr, --start, and --end")

    # Neither mode specified
    if not has_region and not has_input:
        parser.error("Either --input or region mode (--chr, --start, --end) is required")


def main():
    args = parse_args()

    gff_file = args.gff_file
    output_prefix = args.output_prefix

    print_welcome_message()

    # Ensure output directory exists
    os.makedirs(output_prefix, exist_ok=True)

    # Region mode
    if args.chromosome and args.start and args.end:
        genes = parse_gff_for_region(gff_file, args.chromosome, args.start, args.end)

        if not genes:
            print(f"No transcripts found in region {args.chromosome}:{args.start}-{args.end}")
            return

        labels = [g.gene_id for g in genes]

        for gene in genes:
            gene.add_introns()

        output_svg = f"{output_prefix}/{args.chromosome}_{args.start}-{args.end}.svg"

        draw_region_gene_structures(
            genes,
            labels,
            args.start,
            args.end,
            output_svg
        )
        print(f"Finished! : {output_svg}")
        print(f"  Transcripts: {len(genes)}")

    # Transcript mode (CSV input)
    else:
        input_csv = args.input_csv
        with open(input_csv, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                transcript_id = row["transcript_id"]

                snps_raw = row.get("snp", "").replace(";", "_")
                snps = f"_SNP_{snps_raw}" if snps_raw else ""

                deletions_raw = row.get("deletions", "").replace(";", "_")
                deletions = f"_DEL_{deletions_raw}" if deletions_raw else ""

                insertions_raw = row.get("insertions", "").replace(";", "_")
                insertions = f"_INS_{insertions_raw}" if insertions_raw else ""

                domains_raw = row.get("domains", "").replace(":", "_").replace(";", "_")
                domains = f"_DOM_{domains_raw}" if domains_raw else ""

                snp_positions = parse_snps(row.get("snp", ""))
                deletion_regions_relative = parse_deletions(row.get("deletions", ""))
                insertion_positions = parse_insertions(row.get("insertions", ""))
                domain_defs = parse_domains(row.get("domains", ""))

                gene = parse_gff_for_transcript(gff_file, transcript_id)
                if not gene:
                    print(f"Skip: {transcript_id} not found")
                    continue

                gene.add_introns()
                gene.to_relative()

                # domain（AA座標 → ゲノム）
                for start_aa, end_aa, name in domain_defs:
                    gene.add_domain_from_protein_coords(start_aa, end_aa, name)

                # variants
                gene.update_features_with_deletions(deletion_regions_relative)
                gene.add_insertions(insertion_positions)
                gene.add_snps(snp_positions)

                output_svg = f"{output_prefix}/{transcript_id}{snps}{deletions}{insertions}{domains}.svg"

                draw_gene_structure(gene, output_svg)
                print(f"Finished! : {output_svg}")


if __name__ == "__main__":
    main()
