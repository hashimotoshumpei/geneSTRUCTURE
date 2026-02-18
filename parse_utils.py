from gene_classes import GeneFeature, GeneStructure

# =====================
# INPUTパーサ
# ====================

def parse_insertions(insertion_str):
    """
    "100;200;1000" -> [100, 200, 1000]
    """
    if not insertion_str:
        return []
    return [int(x) for x in insertion_str.split(";") if x.strip()]

def parse_deletions(deletion_str):
    """
    "12-2000;3000-3500" -> [(12,2000),(3000,3500)]
    """
    if not deletion_str:
        return []
    regions = []
    for part in deletion_str.split(";"):
        s, e = part.split("-")
        regions.append((int(s), int(e)))
    return regions

def parse_snps(snp_str):
    """
    "100;200;350" -> [100, 200, 350]
    """
    if not snp_str:
        return []
    return [int(x) for x in snp_str.split(";") if x.strip()]

def parse_domains(domain_str):
    """
    "1-200:domain1;201-500:domain2"
    -> [(1,200,'domain1'), (201,500,'domain2')]
    """
    if not domain_str:
        return []
    domains = []
    for part in domain_str.split(";"):
        coord, name = part.split(":")
        s, e = coord.split("-")
        domains.append((int(s), int(e), name))
    return domains


# =====================
# GFFパーサ
# ====================


def parse_gff_for_transcript(gff_file, transcript_id):
    gene_structure = None
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
            if f"Parent={transcript_id}" not in attributes and f"ID={transcript_id}" not in attributes:
                continue
            if gene_structure is None:
                gene_structure = GeneStructure(transcript_id, seqid, strand)
            if strand == '+':
                feature = GeneFeature(seqid, int(start), int(end), feature_type, strand)
            elif strand == '-':
                feature = GeneFeature(seqid, int(end)*-1, int(start)*-1, feature_type, strand)
            gene_structure.add_feature(feature)
    return gene_structure


def parse_gff_for_region(gff_file, seqid, region_start, region_end):
    """
    Extract all transcripts within the specified genomic region.

    Args:
        gff_file: Path to GFF/GTF file
        seqid: Chromosome/SeqID
        region_start: Region start position (1-based)
        region_end: Region end position (1-based)

    Returns:
        List[GeneStructure]: List of GeneStructure objects for transcripts
                             that overlap with the specified region
    """
    # First pass: find all mRNA/transcript IDs in the region
    transcript_ids = set()
    transcript_info = {}  # transcript_id -> (seqid, strand)

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue

            line_seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
            start, end = int(start), int(end)

            # Check if this is on the target chromosome
            if line_seqid != seqid:
                continue

            # Look for mRNA or transcript features
            if feature_type in ('mRNA', 'transcript'):
                # Check if overlaps with region
                if end >= region_start and start <= region_end:
                    # Extract ID from attributes
                    attr_dict = parse_attributes(attributes)
                    transcript_id = attr_dict.get('ID')
                    if transcript_id:
                        transcript_ids.add(transcript_id)
                        transcript_info[transcript_id] = (line_seqid, strand)

    # Second pass: collect features for each transcript
    gene_structures = {}

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) != 9:
                continue

            line_seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
            start, end = int(start), int(end)

            # Parse attributes
            attr_dict = parse_attributes(attributes)

            # Check if this feature belongs to one of our transcripts
            parent_id = attr_dict.get('Parent')
            feature_id = attr_dict.get('ID')

            # Match by Parent or ID
            matched_transcript = None
            if parent_id in transcript_ids:
                matched_transcript = parent_id
            elif feature_id in transcript_ids:
                matched_transcript = feature_id

            if matched_transcript is None:
                continue

            # Create GeneStructure if not exists
            if matched_transcript not in gene_structures:
                t_seqid, t_strand = transcript_info[matched_transcript]
                gene_structures[matched_transcript] = GeneStructure(
                    matched_transcript, t_seqid, t_strand
                )

            # Add feature (keep original genomic coordinates for region mode)
            gene = gene_structures[matched_transcript]
            feature = GeneFeature(line_seqid, start, end, feature_type, strand)
            gene.add_feature(feature)

    # Sort by start position and return as list
    result = list(gene_structures.values())
    result.sort(key=lambda g: min(f.start for f in g.features) if g.features else 0)

    return result


def parse_attributes(attr_string):
    """
    Parse GFF3 attribute string into a dictionary.

    Args:
        attr_string: Attribute string (e.g., "ID=mRNA1;Parent=gene1;Name=foo")

    Returns:
        dict: Attribute key-value pairs
    """
    attr_dict = {}
    for item in attr_string.split(';'):
        item = item.strip()
        if '=' in item:
            key, value = item.split('=', 1)
            attr_dict[key] = value
    return attr_dict



def get_terminal_feature(features):
    """
    右端（最大 end）にある feature を返す。
    優先順位: three_prime_UTR > CDS > exon
    """
    priority = ['three_prime_UTR', 'CDS', 'exon']

    for ftype in priority:
        candidates = [f for f in features if f.feature_type == ftype]
        if candidates:
            return max(candidates, key=lambda f: f.end)

    return None

def print_features_as_gff(gene: GeneStructure):
    for f in gene.get_sorted_features():
        attr_str = ';'.join(f"{k}={v}" for k, v in f.attributes.items()) if f.attributes else '.'
        print("\t".join([
            f.seqid,
            "manual",                # source
            f.feature_type,
            str(f.start),
            str(f.end),
            ".",                     # score
            f.strand,
            ".",                     # phase
            attr_str
        ]))
