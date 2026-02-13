import csv
import argparse
from parse_utils import (
    parse_gff_for_transcript,
    parse_deletions,
    parse_insertions,
    parse_snps,
    parse_domains,
    print_features_as_gff
)
from draw_utils import draw_gene_structure
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
        required=True,
        help="Input CSV file describing gene features"
    )

    parser.add_argument(
        "--output", "-o",
        dest="output_prefix",
        required=True,
        help="Output file prefix or directory (path allowed)"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    gff_file = args.gff_file
    input_csv = args.input_csv
    output_prefix = args.output_prefix

    print_welcome_message()

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

            # 出力
            #print_features_as_gff(gene)


            output_svg = f"{output_prefix}/{transcript_id}{snps}{deletions}{insertions}{domains}.svg"

            draw_gene_structure(gene, output_svg)
            print(f"✏️ Finished! 🎉 : {output_svg}")


if __name__ == "__main__":
    main()
