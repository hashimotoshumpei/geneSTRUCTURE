import svgwrite
from gene_classes import GeneStructure
from parse_utils import get_terminal_feature
from color_utils import get_or_create_gradient
from config import (
    utr_gradation, exon_gradation, domain_gradation,
    FEATURE_COLORS, LEFT_MARGIN
)

# 描画関数（修正）
def draw_gene_structure(gene, output_svg, scale=2, extra_padding=100, shrink_factor=30.0):
    

    min_start = gene.to_relative()
    all_features = gene.get_sorted_features()
    terminal_feature = get_terminal_feature(all_features)


    print(f"terminal_feature: {terminal_feature.feature_type}")
    max_end = max(f.end / shrink_factor for f in all_features)
    shift = -min_start if min_start < 0 else 0
    canvas_width = LEFT_MARGIN + (max_end + shift / shrink_factor) * scale + extra_padding + 300
    canvas_height = 300

    dwg = svgwrite.Drawing(output_svg, size=(canvas_width, canvas_height))
    grad_dict = {}
    y_pos = 50
    height_feature = 15
    max_x_coord = LEFT_MARGIN + (max_end + shift / shrink_factor) * scale

    for i, feat in enumerate(all_features):
        x_start = LEFT_MARGIN + (feat.start / shrink_factor + shift / shrink_factor) * scale
        x_end = LEFT_MARGIN + (feat.end / shrink_factor + shift / shrink_factor) * scale
        width = x_end - x_start

        if feat.feature_type == 'domain':
            continue

        if feat.feature_type == 'deletion':

            print("直前のfeature_type", all_features[i-1].feature_type)

            # くの字型の折れ線（xの半分の長さ）
            y_line = y_pos + height_feature // 2
            mid_x = x_start + (x_end - x_start) / 2
            offset = 10  # くの字の高さ
            dwg.add(
                dwg.polyline(
                    points=[
                        (x_start, y_line),
                        (mid_x, y_line - offset),
                        (x_end, y_line)
                    ],
                    fill='none',
                    stroke='black',
                    stroke_width=1,
		    stroke_dasharray="2,2"
                )
            )
        elif feat.feature_type in ('exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'):
            base_color = FEATURE_COLORS.get(feat.feature_type, 'gray')
            fill_color = base_color

            # グラデーション設定
            if feat.feature_type in ('exon', 'CDS') and exon_gradation == "on":
                fill_color = f'url(#{get_or_create_gradient(dwg, base_color, grad_dict)})'
            elif feat.feature_type in ('five_prime_UTR', 'three_prime_UTR') and utr_gradation == "on":
                fill_color = f'url(#{get_or_create_gradient(dwg, base_color, grad_dict)})'

            if feat is terminal_feature:
                # 右端が尖った polygon
                tip = height_feature // 2

                dwg.add(
                    dwg.polygon(
                        points=[
                            (x_start, y_pos),
                            (x_end - tip, y_pos),
                            (x_end, y_pos + height_feature / 2),
                            (x_end - tip, y_pos + height_feature),
                            (x_start, y_pos + height_feature)
                        ],
                        fill=fill_color,
                        stroke='black',
                        stroke_width=1
                    )
                )
            else:
                # 通常の四角
                dwg.add(
                    dwg.rect(
                        insert=(x_start, y_pos),
                        size=(width, height_feature),
                        fill=fill_color,
                        stroke='black',
                        stroke_width=1
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
                        stroke_width=1
                    )
                )


    # === Insertions ===
    triangle_width = 8
    triangle_height = 6
    y_triangle = y_pos - 8  # exon の少し上

    for ins_pos in getattr(gene, "insertions", []):
        x = LEFT_MARGIN + (ins_pos / shrink_factor + shift / shrink_factor) * scale

        dwg.add(
            dwg.polygon(
                points=[
                    (x - triangle_width / 2, y_triangle),
                    (x + triangle_width / 2, y_triangle),
                    (x, y_triangle + triangle_height)
                ],
                fill="black",
                stroke="black",
                stroke_width=1.5
            )
        )
    
    # === SNPs ===
    snp_extend_up = 8     # 上にどれだけ伸ばすか
    snp_extend_down = 8   # 下にどれだけ伸ばすか

    y_snp_top = y_pos - snp_extend_up
    y_snp_bottom = y_pos + height_feature + snp_extend_down

    for snp_pos in getattr(gene, "snps", []):
        x = LEFT_MARGIN + (snp_pos / shrink_factor + shift / shrink_factor) * scale

        dwg.add(
            dwg.line(
                start=(x, y_snp_top),
                end=(x, y_snp_bottom),
                stroke="black",
                stroke_width=1.2
            )
        )

    # domainは最後に描画
    for feat in all_features:
        if feat.feature_type == 'domain':
            x_start = LEFT_MARGIN + (feat.start / shrink_factor + shift / shrink_factor) * scale
            x_end = LEFT_MARGIN + (feat.end / shrink_factor + shift / shrink_factor) * scale
            width = x_end - x_start

            domain_color = feat.attributes.get('color',
            FEATURE_COLORS.get('domain', 'green')
            )
            if domain_gradation == "on":
                domain_color = f'url(#{get_or_create_gradient(dwg, domain_color, grad_dict)})'

            dwg.add(
                dwg.rect(
                    insert=(x_start, y_pos),
                    size=(width, height_feature),
                    fill=domain_color,
                    stroke='black',
                    stroke_width=1
                )
            )

    # === 凡例に出現する feature type を収集 ===
    present_feature_types = set(f.feature_type for f in all_features)

    # domain 名（色付き）
    present_domains = gene.domain_color_map  # 空なら {} のまま

    # === 凡例 ===
    legend_x = max_x_coord + 100
    legend_y = 30
    box_size = 12
    spacing = 20

    legend_items = []

    if 'CDS' in present_feature_types or 'exon' in present_feature_types:
        legend_items.append(('CDS', 'Exon/CDS'))

    if 'five_prime_UTR' in present_feature_types:
        legend_items.append(('five_prime_UTR', "5' UTR"))

    if 'three_prime_UTR' in present_feature_types:
        legend_items.append(('three_prime_UTR', "3' UTR"))

    if 'intron' in present_feature_types:
        legend_items.append(('intron', 'Intron'))

    if 'deletion' in present_feature_types:
        legend_items.append(('deletion', 'Deletion'))

    has_insertion = bool(getattr(gene, "insertions", []))
    if has_insertion:
        legend_items.append(('insertion', 'Insertion'))

    has_snp = bool(getattr(gene, "snps", []))
    if has_snp:
        legend_items.append(('snp', 'SNP'))



    # domain は domain_color_map に基づいて追加
    for domain_name, color in gene.domain_color_map.items():
        legend_items.append(('domain', domain_name))

    for i, (feat_key, label) in enumerate(legend_items):
        y_legend = legend_y + i * spacing

        # === Deletion ===
        if feat_key == 'deletion':
            y_mid = y_legend + box_size // 2
            x0 = legend_x
            x1 = legend_x + box_size // 2
            x2 = legend_x + box_size

            dwg.add(
                dwg.polyline(
                    points=[
                        (x0, y_mid),
                        (x1, y_mid - 6),
                        (x2, y_mid)
                    ],
                    fill="none",
                    stroke="black",
                    stroke_width=1.5,
		    stroke_dasharray="2,2"
                )
            )

        # === Insertion ===
        elif feat_key == 'insertion':
            y_mid = y_legend + box_size // 2
            x0 = legend_x
            x1 = legend_x + box_size // 2
            x2 = legend_x + box_size

            dwg.add(
                dwg.polygon(
                    points=[
                        (x0, y_mid - 4),
                        (x2, y_mid - 4),
                        (x1, y_mid + 4)
                    ],
                    fill="black",
                    stroke="black",
                    stroke_width=1.5
                )
            )

        # === SNP ===
        elif feat_key == 'snp':
            dwg.add(
                dwg.line(
                    start=(legend_x + box_size // 2, y_legend),
                    end=(legend_x + box_size // 2, y_legend + box_size),
                    stroke="black",
                    stroke_width=1.2
                )
            )


        # === Intron ===
        elif feat_key == 'intron':
            y_line = y_legend + box_size // 2
            dwg.add(
                dwg.line(
                    start=(legend_x, y_line),
                    end=(legend_x + box_size, y_line),
                    stroke=FEATURE_COLORS.get('intron', 'black'),
                    stroke_width=1
                )
            )

        # === Exon / UTR / Domain ===
        else:
            if feat_key == 'domain':
                base_color = gene.domain_color_map[label]
                use_grad = (domain_gradation == "on")
            else:
                base_color = FEATURE_COLORS.get(feat_key, 'gray')
                use_grad = (
                    (feat_key in ('CDS', 'exon') and exon_gradation == "on") or
                    (feat_key in ('five_prime_UTR', 'three_prime_UTR') and utr_gradation == "on")
                )

            fill_color = base_color
            if use_grad:
                fill_color = f'url(#{get_or_create_gradient(dwg, base_color, grad_dict)})'

            dwg.add(
                dwg.rect(
                    insert=(legend_x, y_legend),
                    size=(box_size, box_size),
                    fill=fill_color,
                    stroke="black",
                    stroke_width=1
                )
            )

        # ラベル
        dwg.add(
            dwg.text(
                label,
                insert=(legend_x + box_size + 6, y_legend + box_size - 2),
                font_size="11px",
                fill="black"
            )
    )

    dwg.save()
