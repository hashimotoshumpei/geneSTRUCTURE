from color_utils import get_domain_color
from config import DOMAIN_COLOR_PALETTE

# =====================
# クラス定義
# =====================

class GeneFeature:
    def __init__(self, seqid, start, end, feature_type, strand, attributes=None):
        self.seqid = seqid
        self.start = start
        self.end = end
        self.feature_type = feature_type
        self.strand = strand
        self.attributes = attributes or {}

class GeneStructure:
    def __init__(self, gene_id, seqid, strand):
        self.gene_id = gene_id
        self.seqid = seqid
        self.strand = strand
        self.features = []
        self.insertions = []
        self.snps = []
        self.domain_color_map = {}

    def add_feature(self, feature: GeneFeature):
        self.features.append(feature)

    def get_sorted_features(self):
        #reverse = True if self.strand == '-' else False
        return sorted(self.features, key=lambda f: f.start, reverse=False)

    def add_insertions(self, insertion_positions):
        self.insertions = insertion_positions

    def add_snps(self, snp_positions):
        self.snps = snp_positions

    def add_introns(self):
        # exon / CDS / UTR をまとめて処理
        exon_like_list = sorted(
            [f for f in self.features if f.feature_type in ('exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR')],
            key=lambda x: x.start
        )

        for i in range(len(exon_like_list) - 1):
            intron_start = exon_like_list[i].end + 1
            intron_end = exon_like_list[i + 1].start - 1
            if intron_start <= intron_end:
                intron = GeneFeature(self.seqid, intron_start, intron_end, 'intron', self.strand, {})
                self.features.append(intron)

    def add_domains(self, domain_regions):
        for domain in domain_regions:
            start = domain['start']
            end = domain['end']
            name = domain.get('name', '')
            color = domain.get('color', '')
            domain_feature = GeneFeature(
                self.seqid,
                start,
                end,
                'domain',
                self.strand,
                attributes={'name': name, 'color': color}
            )
            self.features.append(domain_feature)

    def update_features_with_deletions(self, deletion_regions):

        new_features = []

        for i, feature in enumerate(self.features):
            f_start, f_end = feature.start, feature.end
            segments = [(f_start, f_end)]  # featureの元の範囲

            for del_start, del_end in deletion_regions:
                updated_segments = []

                if i == 0:
                    new_features.append(GeneFeature(
                        self.seqid, del_start, del_end,
                        'deletion', self.strand, {}
                    ))

                for seg_start, seg_end in segments:
                    # 削除領域と重なっていなければそのまま残す
                    if seg_end < del_start or seg_start > del_end:
                        updated_segments.append((seg_start, seg_end))
                    else:
                        # 左端が削除領域より前
                        if seg_start < del_start:
                            updated_segments.append((seg_start, del_start - 1))
                        # 右端が削除領域より後
                        if seg_end > del_end:
                            updated_segments.append((del_end + 1, seg_end))
                segments = updated_segments

            # 分割後の有効セグメントが残っていれば追加（完全削除されたらスキップ）
            for start, end in segments:
                if start <= end:
                    new_features.append(GeneFeature(
                        seqid=feature.seqid,
                        start=start,
                        end=end,
                        feature_type=feature.feature_type,
                        strand=feature.strand,
                        attributes=feature.attributes
                    ))

        # 結果を更新
        self.features = new_features


    def to_relative(self):
        cds_list = [f for f in self.features if f.feature_type in ('exon', 'CDS')]
        if not cds_list:
            return 0
        #if self.strand == '+':
        anchor = min(cds_list, key=lambda f: f.start).start
        #else:
            #anchor = max(cds_list, key=lambda f: f.start).end
        for f in self.features:
            f.start = f.start - anchor + 1
            f.end = f.end - anchor + 1
        min_start = min(f.start for f in self.features)
        return min_start

    def add_domain_from_protein_coords(self, start_aa: int, end_aa: int, domain_name: str):
        """
        アミノ酸座標（1-based）を基に、CDSからcDNA、そしてゲノム座標へと変換して
        ドメイン領域をfeaturesに追加する。

        例: start_aa=10, end_aa=100 は、cDNA上で28-300に相当
        """

        # アミノ酸座標 → cDNA 座標（1-based）
        cdna_start = (start_aa - 1) * 3 + 1
        cdna_end = end_aa * 3

        # CDS features を取得してストランド順に並べ替え
        cds_features = [f for f in self.features if f.feature_type == 'CDS']

        cds_sorted = sorted(cds_features, key=lambda f: f.start)
        # if self.strand == '+':
        #     cds_sorted = sorted(cds_features, key=lambda f: f.start)
        # else:
        #     cds_sorted = sorted(cds_features, key=lambda f: f.start, reverse=True)

        gdna_segments = []
        current_cdna_pos = 1

        for cds in cds_sorted:
            cds_len = cds.end - cds.start + 1
            next_cdna_pos = current_cdna_pos + cds_len - 1

            # このCDSにドメインが含まれているか？
            if next_cdna_pos < cdna_start:
                current_cdna_pos = next_cdna_pos + 1
                continue
            if current_cdna_pos > cdna_end:
                break

            # オーバーラップする部分だけを計算
            overlap_start = max(cdna_start, current_cdna_pos)
            overlap_end = min(cdna_end, next_cdna_pos)
            offset_start = overlap_start - current_cdna_pos
            offset_end = overlap_end - current_cdna_pos

            # ゲノム座標に変換

            g_start = cds.start + offset_start
            g_end = cds.start + offset_end
            # if self.strand == '+':
            #     g_start = cds.start + offset_start
            #     g_end = cds.start + offset_end
            # else:
            #     g_end = cds.end - offset_start
            #     g_start = cds.end - offset_end

            # ドメイン feature を追加
            color = get_domain_color(domain_name, self.domain_color_map, DOMAIN_COLOR_PALETTE)

            domain_feature = GeneFeature(
                seqid=self.seqid,
                start=g_start,
                end=g_end,
                feature_type='domain',
                strand=self.strand,
                attributes={'name': domain_name, 'color': color}
            )
            self.features.append(domain_feature)

            current_cdna_pos = next_cdna_pos + 1
