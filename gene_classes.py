from enum import Enum
from color_utils import get_domain_color
from config import DOMAIN_COLOR_PALETTE


# =====================
# 座標モード
# =====================

class CoordinateMode(str, Enum):
    """座標モードの列挙型"""
    RELATIVE = "relative"  # 相対座標（デフォルト）
    ABSOLUTE = "absolute"  # 絶対座標（ゲノム座標）


# =====================
# バリアントクラス
# =====================

class Insertion:
    def __init__(self, position, length):
        self.position = position
        self.length = length

class Snp:
    def __init__(self, position):
        self.position = position

class Deletion:
    def __init__(self, start, end):
        self.start = start
        self.end = end


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
        self.deletion_regions = []
        self.domain_color_map = {}
        self.anchor = 0  # 基準となるゲノム座標を保存

    def add_feature(self, feature: GeneFeature):
        self.features.append(feature)

    def get_sorted_features(self):
        return sorted(self.features, key=lambda f: f.start, reverse=False)

    def add_insertions(self, insertions):
        """
        List[Insertion] または List[int] または List[list] を受け取り、内部で Insertion オブジェクトとして保持
        """
        self.insertions = []
        for ins in insertions:
            if isinstance(ins, Insertion):
                self.insertions.append(ins)
            elif isinstance(ins, (list, tuple)) and len(ins) >= 2:
                pos = ins[0]
                length = ins[1]
                self.insertions.append(Insertion(pos, length))
            else:
                # 単なる位置 (int)
                self.insertions.append(Insertion(ins, 1))

    def add_snps(self, snps):
        """
        List[Snp] または List[int] または List[list] を受け取り保持
        """
        self.snps = []
        for s in snps:
            if isinstance(s, Snp):
                self.snps.append(s)
            elif isinstance(s, (list, tuple)) and len(s) >= 1:
                self.snps.append(Snp(s[0]))
            else:
                self.snps.append(Snp(s))

    def normalize_features(self):
        """
        Feature の正規化処理（イントロン追加を含む）
        1. exon + CDS + UTR → exon を削除
        2. exon + CDS (UTRなし) → exon と CDS の差分から UTR を計算し、exon を削除
        3. exon のみ → そのまま維持
        4. イントロンを追加
        """
        exons = [f for f in self.features if f.feature_type == 'exon']
        cds_list = [f for f in self.features if f.feature_type == 'CDS']
        utrs = [f for f in self.features if f.feature_type in ('five_prime_UTR', 'three_prime_UTR')]

        # Case 1 & 2: CDS がある場合
        if cds_list:
            # UTR がない場合、exon と CDS の差分から UTR を計算
            if not utrs and exons:
                self._compute_utrs_from_exon_cds(exons, cds_list)

            # exon を削除（CDS + UTR で表現するため）
            self.features = [f for f in self.features if f.feature_type != 'exon']

        # Case 3: exon のみの場合はそのまま

        # イントロンを追加
        self.add_introns()

    def _compute_utrs_from_exon_cds(self, exons, cds_list):
        """
        exon と CDS の差分から UTR を計算して追加
        """
        # CDS の全体範囲を取得
        cds_start = min(c.start for c in cds_list)
        cds_end = max(c.end for c in cds_list)

        for exon in exons:
            if self.strand == '-':
                # マイナスストランド: 5' UTR は CDS より大きな座標、3' UTR は CDS より小さな座標
                # 5' UTR: CDS の終了から exon の終了まで
                if exon.end > cds_end and exon.start <= cds_end + 1:
                    utr_start = max(exon.start, cds_end + 1)
                    if utr_start <= exon.end:
                        self.features.append(GeneFeature(
                            self.seqid, utr_start, exon.end,
                            'five_prime_UTR', self.strand, {}
                        ))

                # 3' UTR: exon の開始から CDS の開始まで
                if exon.start < cds_start and exon.end >= cds_start - 1:
                    utr_end = min(exon.end, cds_start - 1)
                    if exon.start <= utr_end:
                        self.features.append(GeneFeature(
                            self.seqid, exon.start, utr_end,
                            'three_prime_UTR', self.strand, {}
                        ))
            else:
                # プラスストランド: 5' UTR は CDS より小さな座標、3' UTR は CDS より大きな座標
                # 5' UTR: exon の開始から CDS の開始まで
                if exon.start < cds_start and exon.end >= cds_start - 1:
                    utr_end = min(exon.end, cds_start - 1)
                    if exon.start <= utr_end:
                        self.features.append(GeneFeature(
                            self.seqid, exon.start, utr_end,
                            'five_prime_UTR', self.strand, {}
                        ))

                # 3' UTR: CDS の終了から exon の終了まで
                if exon.end > cds_end and exon.start <= cds_end + 1:
                    utr_start = max(exon.start, cds_end + 1)
                    if utr_start <= exon.end:
                        self.features.append(GeneFeature(
                            self.seqid, utr_start, exon.end,
                            'three_prime_UTR', self.strand, {}
                        ))

    def add_introns(self):
        # 以前のイントロンがあれば削除
        self.features = [f for f in self.features if f.feature_type != 'intron']
        
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
            # domain_regions can be list of dict or list of tuple
            if isinstance(domain, dict):
                start = domain['start']
                end = domain['end']
                name = domain.get('name', '')
                color = domain.get('color', '')
            else:
                # tuple (start, end, name) or (start, end, name, color)
                start = domain[0]
                end = domain[1]
                name = domain[2] if len(domain) > 2 else ''
                color = domain[3] if len(domain) > 3 else ''

            if not color:
                color = get_domain_color(name, self.domain_color_map, DOMAIN_COLOR_PALETTE)
            else:
                self.domain_color_map[name] = color

            domain_feature = GeneFeature(
                self.seqid,
                start,
                end,
                'domain',
                self.strand,
                attributes={'name': name, 'color': color}
            )
            self.features.append(domain_feature)

    def get_full_extent(self):
        """SNPや挿入を含めた、遺伝子構造の真の開始・終了座標を返す"""
        starts = [f.start for f in self.features]
        ends = [f.end for f in self.features]

        for snp in self.snps:
            pos = getattr(snp, 'position', snp)
            starts.append(pos)
            ends.append(pos)

        for ins in self.insertions:
            pos = getattr(ins, 'position', ins)
            length = getattr(ins, 'length', 1)
            starts.append(pos)
            ends.append(pos + length - 1)

        if not starts:
            return 1, 1
        return min(starts), max(ends)

    def update_features_with_deletions(self, deletion_regions):
        """
        List[Deletion] または List[tuple] または List[dict] を受け取り処理
        """
        self.deletion_regions = []
        for d in deletion_regions:
            if isinstance(d, Deletion):
                self.deletion_regions.append(d)
            elif isinstance(d, (list, tuple)) and len(d) >= 2:
                self.deletion_regions.append(Deletion(d[0], d[1]))
            elif isinstance(d, dict):
                self.deletion_regions.append(Deletion(d['start'], d['end']))

        new_features = []
        structural_types = {'exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'intron', 'mRNA', 'transcript'}

        # まずデリーション自体をフィーチャーとして追加
        for d in self.deletion_regions:
            new_features.append(GeneFeature(
                self.seqid, d.start, d.end,
                'deletion', self.strand, {}
            ))

        for feature in self.features:
            if feature.feature_type == 'deletion':
                continue

            # 非構造的要素（ドメイン等）の場合、デリーションと重なれば削除する
            if feature.feature_type not in structural_types:
                overlaps = False
                for d in self.deletion_regions:
                    if not (feature.end < d.start or feature.start > d.end):
                        overlaps = True
                        break
                if overlaps:
                    continue

            f_start, f_end = feature.start, feature.end
            segments = [(f_start, f_end)]  # featureの元の範囲

            for d in self.deletion_regions:
                del_start, del_end = d.start, d.end
                updated_segments = []

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

            # 分割後の有効セグメントが残っていれば追加
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

        # SNPと挿入のフィルタリング
        if self.deletion_regions:
            self.snps = [
                s for s in self.snps 
                if not any(d.start <= s.position <= d.end for d in self.deletion_regions)
            ]
            self.insertions = [
                i for i in self.insertions 
                if not any(d.start <= i.position <= d.end for d in self.deletion_regions)
            ]

    def to_relative(self):
        exon_like = [f for f in self.features if f.feature_type in ('exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'mRNA', 'transcript')]
        if not exon_like:
            return 0

        # プラス鎖: 最小値が基準 (anchor)
        # マイナス鎖: 最大値が基準 (anchor)
        all_coords = []
        for f in exon_like:
            all_coords.append(f.start)
            all_coords.append(f.end)
        
        if self.strand == '-':
            self.anchor = max(all_coords)
        else:
            self.anchor = min(all_coords)

        # すべてのフィーチャーを変換
        for f in self.features:
            if self.strand == '-':
                s = self.anchor - f.start + 1
                e = self.anchor - f.end + 1
                f.start = min(s, e)
                f.end = max(s, e)
            else:
                f.start = f.start - self.anchor + 1
                f.end = f.end - self.anchor + 1

        # SNPと挿入も変換
        if hasattr(self, 'snps') and self.snps:
            for s in self.snps:
                pos = s.position
                if self.strand == '-':
                    s.position = self.anchor - pos + 1
                else:
                    s.position = pos - self.anchor + 1
        
        if hasattr(self, 'insertions') and self.insertions:
            for ins in self.insertions:
                pos = ins.position
                if self.strand == '-':
                    ins.position = self.anchor - pos + 1
                else:
                    ins.position = pos - self.anchor + 1

        # デリーション領域も変換
        if hasattr(self, 'deletion_regions') and self.deletion_regions:
            for d in self.deletion_regions:
                s_orig, e_orig = d.start, d.end
                if self.strand == '-':
                    s = self.anchor - s_orig + 1
                    e = self.anchor - e_orig + 1
                    d.start, d.end = min(s, e), max(s, e)
                else:
                    d.start, d.end = s_orig - self.anchor + 1, e_orig - self.anchor + 1

        return 1

    def add_domain_from_protein_coords(self, start_aa: int, end_aa: int, domain_name: str):
        """
        アミノ酸座標（1-based）を基に、CDSからcDNA、そしてゲノム座標へと変換して
        ドメイン領域をfeaturesに追加する。
        """

        # アミノ酸座標 → cDNA 座標（1-based）
        cdna_start = (start_aa - 1) * 3 + 1
        cdna_end = end_aa * 3

        # CDS features を取得してストランド順に並べ替え
        cds_features = [f for f in self.features if f.feature_type == 'CDS']
        if self.strand == '-':
            cds_sorted = sorted(cds_features, key=lambda f: f.start, reverse=True)
        else:
            cds_sorted = sorted(cds_features, key=lambda f: f.start)

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
            if self.strand == '-':
                g_end = cds.end - offset_start
                g_start = cds.end - offset_end
            else:
                g_start = cds.start + offset_start
                g_end = cds.start + offset_end

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
