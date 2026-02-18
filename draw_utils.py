import svgwrite
from typing import List
from gene_classes import GeneStructure
from parse_utils import get_terminal_feature
from color_utils import get_or_create_gradient
from config import (
    utr_gradation, exon_gradation, domain_gradation,
    FEATURE_COLORS, LEFT_MARGIN, FEATURE_OUTLINES, FEATURE_OUTLINE_WIDTHS,
    FEATURE_OUTLINE_ENABLED
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

    # === canvas height を legend に合わせて再計算 ===
    gene_bottom = y_pos + height_feature + 20
    legend_bottom = legend_y + len(legend_items) * spacing

    bottom_padding = 20
    final_canvas_height = max(gene_bottom, legend_bottom) + bottom_padding

    # SVG の高さを更新
    dwg['height'] = final_canvas_height


    dwg.save()


# =====================
# Helper functions for region mode
# =====================

def get_tick_params(range_size: int) -> tuple:
    """
    Return appropriate tick interval and unit based on range size.

    Args:
        range_size: Size of the display range (bp)

    Returns:
        (tick_interval, unit_label, divisor)
        e.g., (1000, "kb", 1000) -> tick every 1kb, labels "1 kb", "2 kb"...
    """
    if range_size >= 10_000_000:  # 10Mb or more
        return 1_000_000, "Mb", 1_000_000
    elif range_size >= 1_000_000:  # 1Mb or more
        return 100_000, "kb", 1000
    elif range_size >= 100_000:   # 100kb or more
        return 10_000, "kb", 1000
    elif range_size >= 5_000:     # 5kb or more
        return 1_000, "kb", 1000
    elif range_size >= 1_000:     # 1kb or more
        return 500, "bp", 1
    else:
        return 100, "bp", 1


def draw_region_gene_structures(
    genes: List[GeneStructure],
    labels: List[str],
    region_start: int,
    region_end: int,
    output_svg: str,
    show_labels: bool = True,
    gene_spacing: int = 50,
    label_spacing: int = 10,
    scale: float = 2,
    shrink_factor: float = 30.0
):
    """
    Draw multiple gene structures on a common coordinate axis.
    Non-overlapping genes are placed on the same track (row).

    Args:
        genes: List of GeneStructure objects to draw
        labels: Labels for each gene
        region_start: Start coordinate of display region (genomic coordinate)
        region_end: End coordinate of display region (genomic coordinate)
        output_svg: Output SVG file path
        show_labels: Whether to show labels
        gene_spacing: Spacing between tracks (pixels)
        label_spacing: Spacing between label and gene structure (pixels)
        scale: Scale factor
        shrink_factor: Coordinate shrink factor
    """
    height_feature = 15

    # Calculate coordinate range for each gene
    gene_ranges = []
    for idx, gene in enumerate(genes):
        features = gene.get_sorted_features()
        if features:
            gene_start = min(f.start for f in features)
            gene_end = max(f.end for f in features)
        else:
            gene_start = 0
            gene_end = 0
        gene_ranges.append({
            'idx': idx,
            'gene': gene,
            'label': labels[idx],
            'start': gene_start,
            'end': gene_end
        })

    # Sort by start coordinate
    gene_ranges.sort(key=lambda x: x['start'])

    # Track assignment algorithm (non-overlapping genes on same track)
    tracks = []  # Track end coordinates
    gene_track_assignments = []  # (gene_info, track_idx)

    for gene_info in gene_ranges:
        gene_start = gene_info['start']
        gene_end = gene_info['end']

        # Find available track
        track_found = False
        min_gap = 500  # Minimum gap between genes (bp)

        for track_idx, track_end in enumerate(tracks):
            if gene_start > track_end + min_gap:
                # Can place on this track
                tracks[track_idx] = gene_end
                gene_track_assignments.append((gene_info, track_idx))
                track_found = True
                break

        if not track_found:
            # Create new track
            tracks.append(gene_end)
            gene_track_assignments.append((gene_info, len(tracks) - 1))

    num_tracks = len(tracks)

    # Calculate coordinate range (including overflow)
    all_starts = [g['start'] for g in gene_ranges if g['start'] > 0]
    all_ends = [g['end'] for g in gene_ranges if g['end'] > 0]

    # Determine drawing range
    if all_starts and all_ends:
        draw_start = min(region_start, min(all_starts))
        draw_end = max(region_end, max(all_ends))
    else:
        draw_start = region_start
        draw_end = region_end

    # Calculate axis width
    axis_width = (draw_end - draw_start) / shrink_factor * scale

    # Canvas dimensions
    extra_padding = 100
    canvas_width = LEFT_MARGIN + axis_width + extra_padding + 300

    label_height = 15 if show_labels else 0
    track_height = height_feature + label_height + label_spacing
    top_margin = 50  # Space for coordinate axis
    canvas_height = top_margin + num_tracks * (track_height + gene_spacing) + 150

    # Create SVG
    dwg = svgwrite.Drawing(output_svg, size=(canvas_width, canvas_height))
    grad_dict = {}

    # Draw coordinate axis (top)
    axis_y = top_margin - 20
    dwg.add(dwg.line(
        start=(LEFT_MARGIN, axis_y),
        end=(LEFT_MARGIN + axis_width, axis_y),
        stroke='black',
        stroke_width=1
    ))

    # Draw tick marks
    tick_interval, unit_label, divisor = get_tick_params(draw_end - draw_start)
    first_tick = ((draw_start // tick_interval) + 1) * tick_interval

    for tick_pos in range(first_tick, draw_end + 1, tick_interval):
        x = LEFT_MARGIN + (tick_pos - draw_start) / shrink_factor * scale

        # Tick line
        dwg.add(dwg.line(
            start=(x, axis_y),
            end=(x, axis_y + 5),
            stroke='black',
            stroke_width=1
        ))

        # Label
        if divisor == 1:
            tick_label = f"{tick_pos} {unit_label}"
        else:
            tick_label = f"{tick_pos // divisor} {unit_label}"
        dwg.add(dwg.text(
            tick_label,
            insert=(x, axis_y - 3),
            font_size='9px',
            fill='black',
            text_anchor='middle'
        ))

    # Draw each gene
    for gene_info, track_idx in gene_track_assignments:
        gene = gene_info['gene']
        label = gene_info['label']
        all_features = gene.get_sorted_features()
        y_pos = top_margin + track_idx * (track_height + gene_spacing)

        # Calculate gene center X (for label placement)
        gene_center_x = None
        if all_features:
            gene_start = min(f.start for f in all_features)
            gene_end = max(f.end for f in all_features)
            gene_center_x = LEFT_MARGIN + ((gene_start + gene_end) / 2 - draw_start) / shrink_factor * scale

        # Collect feature types for terminal feature detection
        terminal_feature = get_terminal_feature(all_features)

        # Draw features (except domain)
        for feat in all_features:
            x_start = LEFT_MARGIN + (feat.start - draw_start) / shrink_factor * scale
            x_end = LEFT_MARGIN + (feat.end - draw_start) / shrink_factor * scale
            width = x_end - x_start

            if feat.feature_type == 'domain':
                continue

            if feat.feature_type == 'deletion':
                # Zigzag line for deletion
                y_line = y_pos + height_feature // 2
                mid_x = x_start + (x_end - x_start) / 2
                offset = 10
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

                # Gradient settings
                if feat.feature_type in ('exon', 'CDS') and exon_gradation == "on":
                    fill_color = f'url(#{get_or_create_gradient(dwg, base_color, grad_dict)})'
                elif feat.feature_type in ('five_prime_UTR', 'three_prime_UTR') and utr_gradation == "on":
                    fill_color = f'url(#{get_or_create_gradient(dwg, base_color, grad_dict)})'

                stroke_color = FEATURE_OUTLINES.get(feat.feature_type, 'black')
                stroke_width = FEATURE_OUTLINE_WIDTHS.get(feat.feature_type, 1)

                if feat is terminal_feature:
                    # Arrow-shaped polygon for terminal feature
                    tip = height_feature // 2
                    dwg.add(
                        dwg.polygon(
                            points=[
                                (x_start, y_pos),
                                (x_end, y_pos),
                                (x_end + tip, y_pos + height_feature / 2),
                                (x_end, y_pos + height_feature),
                                (x_start, y_pos + height_feature)
                            ],
                            fill=fill_color,
                            stroke=stroke_color,
                            stroke_width=stroke_width
                        )
                    )
                else:
                    dwg.add(
                        dwg.rect(
                            insert=(x_start, y_pos),
                            size=(width, height_feature),
                            fill=fill_color,
                            stroke=stroke_color,
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

        # Draw domains (on top)
        for feat in all_features:
            if feat.feature_type == 'domain':
                x_start = LEFT_MARGIN + (feat.start - draw_start) / shrink_factor * scale
                x_end = LEFT_MARGIN + (feat.end - draw_start) / shrink_factor * scale
                width = x_end - x_start

                domain_color = feat.attributes.get('color', FEATURE_COLORS.get('domain', 'green'))
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

        # Draw label below gene structure (centered)
        if show_labels and gene_center_x is not None:
            dwg.add(dwg.text(
                label,
                insert=(gene_center_x, y_pos + height_feature + label_spacing + 10),
                font_size='10px',
                fill='black',
                font_family='monospace',
                text_anchor='middle'
            ))

    # === Legend (top right) ===
    legend_x = LEFT_MARGIN + axis_width + 50
    legend_y = 30
    box_size = 12
    spacing = 20

    # Collect present feature types
    present_feature_types = set()
    for gene in genes:
        for f in gene.get_sorted_features():
            present_feature_types.add(f.feature_type)

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

    if 'domain' in present_feature_types:
        legend_items.append(('domain', 'Domain'))

    for i, (feat_key, label_text) in enumerate(legend_items):
        y_legend = legend_y + i * spacing

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
        elif feat_key == 'intron':
            y_line = y_legend + box_size // 2
            dwg.add(dwg.line(
                start=(legend_x, y_line),
                end=(legend_x + box_size, y_line),
                stroke=FEATURE_COLORS.get('intron', 'black'),
                stroke_width=1
            ))
        else:
            base_color = FEATURE_COLORS.get(feat_key, 'gray')
            fill_color = base_color

            use_grad = (
                (feat_key in ('CDS', 'exon') and exon_gradation == "on") or
                (feat_key in ('five_prime_UTR', 'three_prime_UTR') and utr_gradation == "on") or
                (feat_key == 'domain' and domain_gradation == "on")
            )
            if use_grad:
                fill_color = f'url(#{get_or_create_gradient(dwg, base_color, grad_dict)})'

            dwg.add(dwg.rect(
                insert=(legend_x, y_legend),
                size=(box_size, box_size),
                fill=fill_color,
                stroke='black'
            ))

        dwg.add(dwg.text(
            label_text,
            insert=(legend_x + box_size + 5, y_legend + box_size - 2),
            font_size='12px',
            fill='black'
        ))

    dwg.save()
