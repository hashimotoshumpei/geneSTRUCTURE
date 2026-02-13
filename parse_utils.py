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
