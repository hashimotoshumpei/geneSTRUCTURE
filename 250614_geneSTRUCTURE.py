import svgwrite

# =====================
# ★描画の色やスタイル設定
# =====================

FEATURE_COLORS = {
    'exon': 'lightblue',
    'CDS': 'lightblue',
    'five_prime_UTR': 'orange',
    'three_prime_UTR': 'lightgreen',
    'intron': 'black',
    'domain': 'green',
    'deletion': 'none',  # 塗りつぶしなし（点線表示用）
    'highlight_intron': 'blue',
}

FEATURE_OUTLINES = {
    'exon': 'black',
    'CDS': 'black',
    'five_prime_UTR': 'black',
    'three_prime_UTR': 'black',
    'domain': 'black',
}

FEATURE_OUTLINE_ENABLED = {
    'exon': True,
    'CDS': True,
    'five_prime_UTR': True,
    'three_prime_UTR': True,
    'domain': True,
}

FEATURE_OUTLINE_WIDTHS = {
    'exon': 1,
    'CDS': 1,
    'five_prime_UTR': 1,
    'three_prime_UTR': 1,
    'domain': 1,
    'intron': 1,
}

LEFT_MARGIN = 50  # 左側マージン

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

    def add_feature(self, feature: GeneFeature):
        self.features.append(feature)

    def get_sorted_features(self):
        #reverse = True if self.strand == '-' else False
        return sorted(self.features, key=lambda f: f.start, reverse=False)

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
        if self.strand == '-':
            cds_sorted = sorted(cds_features, key=lambda f: f.start)
        else:
            cds_sorted = sorted(cds_features, key=lambda f: f.start, reverse=True)

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
            if self.strand == '-':
                g_start = cds.start + offset_start
                g_end = cds.start + offset_end
            else:
                g_end = cds.end - offset_start
                g_start = cds.end - offset_end

            # ドメイン feature を追加
            domain_feature = GeneFeature(
                seqid=self.seqid,
                start=g_start,
                end=g_end,
                feature_type='domain',
                strand=self.strand,
                attributes={'name': domain_name}
            )
            self.features.append(domain_feature)

            current_cdna_pos = next_cdna_pos + 1



# =====================
# GFFパーサ
# =====================

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

# =====================
# 描画関数
# =====================

def draw_gene_structure(gene: GeneStructure, output_svg: str, scale=2, extra_padding=100, shrink_factor=30.0):

    min_start = gene.to_relative()
    all_features = gene.get_sorted_features()
    max_end = max(f.end / shrink_factor for f in all_features)

    shift = -min_start if min_start < 0 else 0
    
    
    canvas_width = LEFT_MARGIN + (max_end + shift / shrink_factor) * scale + extra_padding + 300
    canvas_height = 300  # 凡例分のスペースを確保

    dwg = svgwrite.Drawing(output_svg, size=(canvas_width, canvas_height))
    y_pos = 50
    height_feature = 15
    max_x_coord = LEFT_MARGIN + (max_end + shift / shrink_factor) * scale

    for feat in all_features:
        x_start = LEFT_MARGIN + (feat.start / shrink_factor + shift / shrink_factor) * scale
        x_end = LEFT_MARGIN + (feat.end / shrink_factor + shift / shrink_factor) * scale
        width = x_end - x_start

        if feat.feature_type == 'domain':
            continue 

        if feat.feature_type == 'deletion':
            dwg.add(
                dwg.rect(
                    insert=(x_start, y_pos),
                    size=(width, height_feature),
                    fill='none',
                    stroke='red',
                    stroke_dasharray="5,5",
                    stroke_width=2
                )
            )
        elif feat.feature_type in ('exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'):
            fill_color = FEATURE_COLORS.get(feat.feature_type, 'gray')
            stroke_color = FEATURE_OUTLINES.get(feat.feature_type, 'black')
            stroke_width = FEATURE_OUTLINE_WIDTHS.get(feat.feature_type, 1)
            outline_enabled = FEATURE_OUTLINE_ENABLED.get(feat.feature_type, True)

            dwg.add(
                dwg.rect(
                    insert=(x_start, y_pos),
                    size=(width, height_feature),
                    fill=fill_color,
                    stroke=stroke_color if outline_enabled else 'none',
                    stroke_width=stroke_width
                )
            )
        elif feat.feature_type == 'intron':
            if x_start < x_end:
                y_line = y_pos + height_feature // 2
                dwg.add(
                    dwg.line(
                        start=(x_start, y_line),
                        end=(x_end, y_line),
                        stroke=FEATURE_COLORS.get('intron', 'black'),
                        stroke_width=FEATURE_OUTLINE_WIDTHS.get('intron', 1)
                    )
                )
    for feat in all_features:
        if feat.feature_type == 'domain':
            x_start = LEFT_MARGIN + (feat.start / shrink_factor + shift / shrink_factor) * scale
            x_end = LEFT_MARGIN + (feat.end / shrink_factor + shift / shrink_factor) * scale
            width = x_end - x_start

            dwg.add(
                dwg.rect(
                    insert=(x_start, y_pos),
                    size=(width, height_feature),
                    fill=FEATURE_COLORS.get('domain', 'green'),
                    stroke=FEATURE_OUTLINES.get('domain', 'black'),
                    stroke_width=FEATURE_OUTLINE_WIDTHS.get('domain', 1)
                )
            )

    # === 凡例 ===
    legend_x = max_x_coord + 100
    legend_y = 30
    box_size = 12
    spacing = 20
    legend_items = [
        ('domain', 'Domain'),
        ('CDS', 'Exon/CDS'),
        ('five_prime_UTR', "5' UTR"),
        ('three_prime_UTR', "3' UTR"),
        ('intron', 'Intron'),
        ('deletion', 'Deletion')
    ]
    for i, (feat_key, label) in enumerate(legend_items):
        y_legend = legend_y + i * spacing
        if feat_key == 'deletion':
            dwg.add(dwg.rect(
                insert=(legend_x, y_legend),
                size=(box_size, box_size),
                fill='none',
                stroke='red',
                stroke_dasharray="5,5",
                stroke_width=2
            ))
        elif feat_key == 'intron':
            y_line = y_legend + box_size // 2
            dwg.add(dwg.line(
                start=(legend_x, y_line),
                end=(legend_x + box_size, y_line),
                stroke=FEATURE_COLORS.get('intron', 'black'),
                stroke_width=FEATURE_OUTLINE_WIDTHS.get('intron', 1)
            ))
        else:
            color = FEATURE_COLORS.get(feat_key, 'gray')
            dwg.add(dwg.rect(
                insert=(legend_x, y_legend),
                size=(box_size, box_size),
                fill=color,
                stroke='black'
            ))
        dwg.add(dwg.text(
            label,
            insert=(legend_x + box_size + 5, y_legend + box_size - 2),
            font_size='12px',
            fill='black'
        ))

    dwg.save()

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


# =====================
# インプット例と実行
# =====================

gff_file = './gff3/IRGSP-1.0_representative/transcripts.gff'
#transcript_id = 'Os01t0100100-01'
#transcript_id = 'Os04t0648800-01'
transcript_id = 'Os06t0160700-01'
### 追加して表示したい情報 ###
#deletion_regions_relative = [(12,2000)]
deletion_regions_relative = []

domains = [
    {'start':   200, 'end': 500, 'name': 'Kinase', 'color': 'red'},
    {'start': 600, 'end': 800, 'name': 'ATPase', 'color': 'blue'}
]
#domains = []


# =====================
# メイン処理
# =====================

gene = parse_gff_for_transcript(gff_file, transcript_id)

if gene:
    #print_features_as_gff(gene)
    gene.add_introns()
    #print_features_as_gff(gene)
    gene.to_relative()
    print_features_as_gff(gene)
    gene.add_domain_from_protein_coords(1, 120, 'domain1')
    #gene.add_domains(domains)

    # デリーション処理
    gene.update_features_with_deletions(deletion_regions_relative)

    # 出力
    print_features_as_gff(gene)
    output_svg = f'{transcript_id}_with_relative_deletions.svg'
    draw_gene_structure(gene, output_svg)
    print(f'描画しました: {output_svg}')
