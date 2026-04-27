import math
import svgwrite
from typing import List
from gene_classes import GeneStructure, CoordinateMode
from parse_utils import get_terminal_feature
from color_utils import get_or_create_gradient
from config import (
    utr_gradation, exon_gradation, domain_gradation,
    FEATURE_COLORS, LEFT_MARGIN, FEATURE_OUTLINES, FEATURE_OUTLINE_WIDTHS,
    FEATURE_OUTLINE_ENABLED
)


# =====================
# スケールバー関数
# =====================

def get_tick_params(range_size: int, shrink_factor: float = 30.0, scale: float = 2.0) -> tuple:
    """
    範囲サイズと物理的なスケールに応じて、重なり合わない適切な目盛り間隔と単位を返す
    """
    # 目標とする最小ピクセル間隔（ラベルが重ならないように）
    min_pixel_step = 50 
    # 最小ピクセル間隔を bp に換算
    min_bp_step = min_pixel_step * shrink_factor / scale
    
    # 小さな範囲でも最低限の目盛りが出るように調整
    if range_size <= 0:
        return 1, "bp", 1

    # 1, 2, 5 の倍数の中から、min_bp_step に近い適切な値を探す
    exponent = math.floor(math.log10(min_bp_step))
    magnitude = 10 ** exponent
    
    candidates = [1 * magnitude, 2 * magnitude, 5 * magnitude, 10 * magnitude]
    step = 10 * magnitude
    for c in candidates:
        if c >= min_bp_step:
            step = c
            break
    
    # 範囲に対して目盛りが少なすぎる（3個未満）場合は、一段階細かいステップを検討
    if range_size / step < 2.5 and step / 2 >= min_bp_step / 2:
        if step == 10 * magnitude: step = 5 * magnitude
        elif step == 5 * magnitude: step = 2 * magnitude
        elif step == 2 * magnitude: step = 1 * magnitude
        else: step = magnitude / 2
            
    step = int(step) if step >= 1 else 1
    
    # 単位の決定
    if step >= 1_000_000:
        return step, "Mb", 1_000_000
    elif step >= 1_000:
        return step, "kb", 1000
    else:
        return step, "bp", 1


def get_insertion_base_width(length_bp: int, shrink_factor: float, scale: float) -> float:
    """
    挿入の長さに応じて逆三角形の底辺幅を計算
    """
    # 実際のbp長をスケール変換
    scaled_width = (length_bp / shrink_factor) * scale

    # 最小幅と最大幅を設定
    min_width = 8
    max_width = 40

    return max(min_width, min(scaled_width, max_width))


def get_baseline_segments(actual_min_start: int, actual_max_end: int, deletion_regions: List[any]) -> List[tuple]:
    """
    全体の開始・終了座標とデリーション領域を基に、
    デリーションを避けたベースラインのセグメントリストを返す
    """
    if actual_min_start >= actual_max_end:
        return []
    
    segments = [(actual_min_start, actual_max_end)]
    
    for deletion in deletion_regions:
        if hasattr(deletion, 'start'):
            del_start, del_end = deletion.start, deletion.end
        elif isinstance(deletion, dict):
            del_start, del_end = deletion['start'], deletion['end']
        elif isinstance(deletion, (list, tuple)) and len(deletion) == 2:
            del_start, del_end = deletion
        else:
            continue
            
        new_segments = []
        for seg_start, seg_end in segments:
            if seg_end < del_start or seg_start > del_end:
                new_segments.append((seg_start, seg_end))
            else:
                if seg_start < del_start:
                    new_segments.append((seg_start, del_start - 1))
                if seg_end > del_end:
                    new_segments.append((del_end + 1, seg_end))
        segments = new_segments
        
    return [s for s in segments if s[0] < s[1]]


# 描画関数
def draw_gene_structure(gene, output_svg, scale=2, extra_padding=100, shrink_factor=30.0,
                        coordinate_mode="relative"):
    """
    遺伝子構造をSVGに描画する
    """
    # 描画用に全フィーチャーをソートして取得
    all_features = gene.get_sorted_features()
    terminal_feature = get_terminal_feature(all_features)
    
    # Calculate true extents including SNPs and Insertions
    actual_min_start, actual_max_end = gene.get_full_extent()

    # 描画用にシフト (内部座標を0付近に)
    shift = -actual_min_start
    range_bp = actual_max_end - actual_min_start

    # 座標軸用スペース
    axis_height = 40

    canvas_width = LEFT_MARGIN + (range_bp / shrink_factor) * scale + extra_padding + 300
    canvas_height = 300 + axis_height

    dwg = svgwrite.Drawing(output_svg, size=(canvas_width, canvas_height))
    grad_dict = {}
    y_pos = 50 + axis_height
    height_feature = 15
    max_x_coord = LEFT_MARGIN + (range_bp / shrink_factor) * scale

    # === 座標軸の描画 ===
    axis_y = 30
    x_axis_start = LEFT_MARGIN
    x_axis_end = max_x_coord

    # 座標軸の線
    dwg.add(dwg.line(
        start=(x_axis_start, axis_y),
        end=(x_axis_end, axis_y),
        stroke='black',
        stroke_width=1
    ))

    # 目盛りの計算
    tick_interval, unit_label, divisor = get_tick_params(range_bp, shrink_factor, scale)

    # coordinate_mode に応じて表示用の開始座標を決定
    display_anchor = gene.anchor if coordinate_mode == "absolute" else 1
    
    # 良い感じの目盛り値を計算するために、表示値ベースで最初の目盛りを決定
    if coordinate_mode == "absolute" and gene.strand == '-':
        # マイナスストランドの場合、tick_val が増えると display_tick_val は減る
        max_display_val = display_anchor - actual_min_start + 1
        first_tick_label = math.floor(max_display_val / tick_interval) * tick_interval
        first_tick_val = display_anchor - first_tick_label + 1
    else:
        min_display_val = display_anchor + actual_min_start - 1
        first_tick_label = math.ceil(min_display_val / tick_interval) * tick_interval
        first_tick_val = first_tick_label - display_anchor + 1

    # tick_val は内部相対座標
    for tick_val in range(int(first_tick_val), int(actual_max_end) + 1, int(tick_interval)):
        if tick_val < actual_min_start - 0.1 or tick_val > actual_max_end + 0.1:
            continue
            
        x = LEFT_MARGIN + ((tick_val - actual_min_start) / shrink_factor) * scale

        # 目盛り線
        dwg.add(dwg.line(
            start=(x, axis_y),
            end=(x, axis_y + 5),
            stroke='black',
            stroke_width=1
        ))

        # ラベル
        if coordinate_mode == "absolute":
            if gene.strand == '-':
                display_tick_val = display_anchor - tick_val + 1
            else:
                display_tick_val = display_anchor + tick_val - 1
        else:
            display_tick_val = tick_val
            
        display_tick_val = abs(display_tick_val)

        if divisor == 1:
            tick_label = f"{display_tick_val} {unit_label}"
        else:
            tick_label = f"{display_tick_val // divisor} {unit_label}"

        dwg.add(dwg.text(
            tick_label,
            insert=(x, axis_y - 5),
            font_size='9px',
            fill='black',
            text_anchor='middle'
        ))

    # イントロン風のベースラインを描画 (デリーション領域を避ける)
    deletion_list = [f for f in all_features if f.feature_type == 'deletion']
    baseline_segments = get_baseline_segments(actual_min_start, actual_max_end, deletion_list)
    y_line = y_pos + height_feature // 2
    for seg_start, seg_end in baseline_segments:
        x_base_start = LEFT_MARGIN + (seg_start + shift) / shrink_factor * scale
        x_base_end = LEFT_MARGIN + (seg_end + shift) / shrink_factor * scale
        dwg.add(
            dwg.line(
                start=(x_base_start, y_line),
                end=(x_base_end, y_line),
                stroke=FEATURE_COLORS.get('intron', 'black'),
                stroke_width=FEATURE_OUTLINE_WIDTHS.get('intron', 1)
            )
        )

    for feat in all_features:
        x_start = LEFT_MARGIN + (feat.start + shift) / shrink_factor * scale
        x_end = LEFT_MARGIN + (feat.end + shift) / shrink_factor * scale
        width = x_end - x_start

        if feat.feature_type == 'domain':
            continue

        if feat.feature_type == 'deletion':
            # くの字型の折れ線
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

            stroke_color = FEATURE_OUTLINES.get(feat.feature_type, 'black')
            stroke_width = FEATURE_OUTLINE_WIDTHS.get(feat.feature_type, 1)

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
                        stroke=stroke_color,
                        stroke_width=stroke_width
                    )
                )
            else:
                # 通常の四角
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
            # Skip drawing here as we used baseline_segments for introns
            pass

    # === Insertions ===
    triangle_height = 6
    y_triangle = y_pos - 8  # exon の少し上

    for ins in getattr(gene, "insertions", []):
        if hasattr(ins, 'position'):
            ins_pos = ins.position
            ins_length = getattr(ins, 'length', 1)
            ins_color = getattr(ins, 'color', 'black')
        else:
            ins_pos = ins
            ins_length = 1
            ins_color = "black"

        x = LEFT_MARGIN + (ins_pos + shift) / shrink_factor * scale
        # 挿入の長さに応じて幅を計算
        base_width = get_insertion_base_width(ins_length, shrink_factor, scale)

        dwg.add(
            dwg.polygon(
                points=[
                    (x - base_width / 2, y_triangle),
                    (x + base_width / 2, y_triangle),
                    (x, y_triangle + triangle_height)
                ],
                fill=ins_color,
                stroke=ins_color,
                stroke_width=1.5
            )
        )
    
    # === SNPs ===
    snp_extend_up = 8     # 上にどれだけ伸ばすか
    snp_extend_down = 8   # 下にどれだけ伸ばすか

    y_snp_top = y_pos - snp_extend_up
    y_snp_bottom = y_pos + height_feature + snp_extend_down

    for snp in getattr(gene, "snps", []):
        if hasattr(snp, 'position'):
            snp_pos = snp.position
            snp_color = getattr(snp, 'color', 'black')
        else:
            snp_pos = snp
            snp_color = "black"

        x = LEFT_MARGIN + (snp_pos + shift) / shrink_factor * scale
        dwg.add(
            dwg.line(
                start=(x, y_snp_top),
                end=(x, y_snp_bottom),
                stroke=snp_color,
                stroke_width=1.2
            )
        )

    # domain
    for feat in all_features:
        if feat.feature_type == 'domain':
            x_start = LEFT_MARGIN + (feat.start + shift) / shrink_factor * scale
            x_end = LEFT_MARGIN + (feat.end + shift) / shrink_factor * scale
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

    # === 凡例に出現する feature type を収集 ===
    present_feature_types = set(f.feature_type for f in all_features)

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
    if 'intron' in present_feature_types or baseline_segments:
        legend_items.append(('intron', 'Intron'))
    if 'deletion' in present_feature_types:
        legend_items.append(('deletion', 'Deletion'))
    if getattr(gene, "insertions", []):
        legend_items.append(('insertion', 'Insertion'))
    if getattr(gene, "snps", []):
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
    shrink_factor: float = 30.0,
    coordinate_mode: str = "absolute"
):
    """
    共通座標軸上に複数の遺伝子構造を描画する
    座標が重複しない遺伝子は同じトラック（行）に横並びで配置
    """
    height_feature = 15

    # 各遺伝子の座標範囲を計算
    gene_ranges = []
    for idx, gene in enumerate(genes):
        # Calculate true extents including SNPs and Insertions
        gs, ge = gene.get_full_extent()
            
        gene_ranges.append({
            'idx': idx, 'gene': gene, 'label': labels[idx], 'start': gs, 'end': ge
        })

    # 開始座標でソート
    gene_ranges.sort(key=lambda x: x['start'])

    # トラック配置アルゴリズム（重複しない遺伝子は同じトラックに配置）
    tracks = []  # 各トラックの終了座標を保持
    gene_track_assignments = []  # 各遺伝子のトラック番号

    for gene_info in gene_ranges:
        gs, ge = gene_info['start'], gene_info['end']

        # 配置可能なトラックを探す（余白を考慮）
        track_found = False
        min_gap = 500  # 遺伝子間の最小間隔（bp）

        for track_idx, track_end in enumerate(tracks):
            if gs > track_end + min_gap:
                # このトラックに配置可能
                tracks[track_idx] = ge
                gene_track_assignments.append((gene_info, track_idx))
                track_found = True
                break

        if not track_found:
            # 新しいトラックを作成
            tracks.append(ge)
            gene_track_assignments.append((gene_info, len(tracks) - 1))

    num_tracks = len(tracks)

    # 全遺伝子の座標範囲を計算（はみ出しを含む）
    all_starts = [g['start'] for g in gene_ranges if g['start'] > 0]
    all_ends = [g['end'] for g in gene_ranges if g['end'] > 0]

    # 描画範囲を決定（領域指定とはみ出しを考慮）
    draw_start = min(region_start, min(all_starts)) if all_starts else region_start
    draw_end = max(region_end, max(all_ends)) if all_ends else region_end
    
    # 座標軸の幅を計算
    range_bp = draw_end - draw_start
    axis_width = (range_bp / shrink_factor) * scale

    # Canvas幅
    extra_padding = 100
    canvas_width = LEFT_MARGIN + axis_width + extra_padding + 300

    # Canvas高さ（ラベルは遺伝子構造の下に表示するため、トラックごとに追加スペース）
    label_height = 15 if show_labels else 0
    track_height = height_feature + label_height + label_spacing
    top_margin = 50  # 座標軸用のスペース
    canvas_height = top_margin + num_tracks * (track_height + gene_spacing) + 150

    # メモリ上にSVGを作成
    dwg = svgwrite.Drawing(output_svg, size=(canvas_width, canvas_height))
    grad_dict = {}

    # 座標軸を描画（上部）
    axis_y = top_margin - 20
    dwg.add(dwg.line(start=(LEFT_MARGIN, axis_y), end=(LEFT_MARGIN + axis_width, axis_y), stroke='black', stroke_width=1))

    # 目盛りを描画
    tick_interval, unit_label, divisor = get_tick_params(range_bp, shrink_factor, scale)
    first_tick = math.ceil(draw_start / tick_interval) * tick_interval

    for tick_pos in range(int(first_tick), int(draw_end) + 1, int(tick_interval)):
        # 描画範囲外の tick は描画しない
        if tick_pos < draw_start - 0.1 or tick_pos > draw_end + 0.1:
            continue
        x = LEFT_MARGIN + (tick_pos - draw_start) / shrink_factor * scale
        # 目盛り線
        dwg.add(dwg.line(start=(x, axis_y), end=(x, axis_y + 5), stroke='black', stroke_width=1))

        # ラベル (coordinate_mode に応じて表示値を変える)
        if coordinate_mode == "relative":
            display_tick_val = tick_pos - draw_start + 1
        else:
            # 絶対座標の場合は絶対値（プラス表示）を保証
            display_tick_val = abs(tick_pos)

        if divisor == 1:
            tick_label = f"{display_tick_val} {unit_label}"
        else:
            tick_label = f"{display_tick_val // divisor} {unit_label}"
        dwg.add(dwg.text(tick_label, insert=(x, axis_y - 3), font_size='9px', fill='black', text_anchor='middle'))

    # 各遺伝子を描画
    for gene_info, track_idx in gene_track_assignments:
        gene = gene_info['gene']
        label = gene_info['label']
        all_features = gene.get_sorted_features()
        y_pos = top_margin + track_idx * (track_height + gene_spacing)
        
        # 遺伝子の中心X座標を計算（ラベル配置用）
        gene_center_x = None
        if all_features:
            gs, ge = min(f.start for f in all_features), max(f.end for f in all_features)
            gene_center_x = LEFT_MARGIN + ((gs + ge) / 2 - draw_start) / shrink_factor * scale

        terminal_feature = get_terminal_feature(all_features)
        
        deletion_list = [f for f in all_features if f.feature_type == 'deletion']
        # イントロン風のベースラインを描画 (デリーション領域を避ける)
        baseline_segments = get_baseline_segments(max(draw_start, gene_info['start']), min(draw_end, gene_info['end']), deletion_list)
        y_line = y_pos + height_feature // 2
        for s_start, s_end in baseline_segments:
            xb_start = LEFT_MARGIN + (s_start - draw_start) / shrink_factor * scale
            xb_end = LEFT_MARGIN + (s_end - draw_start) / shrink_factor * scale
            dwg.add(dwg.line(start=(xb_start, y_line), end=(xb_end, y_line), stroke=FEATURE_COLORS.get('intron', 'black'), stroke_width=FEATURE_OUTLINE_WIDTHS.get('intron', 1)))

        # フィーチャーを描画（ドメイン以外）
        for feat in all_features:
            # X座標 = 描画範囲の開始位置からのオフセット
            x_start = LEFT_MARGIN + (feat.start - draw_start) / shrink_factor * scale
            x_end = LEFT_MARGIN + (feat.end - draw_start) / shrink_factor * scale
            width = x_end - x_start

            if feat.feature_type == 'domain': continue
            if feat.feature_type == 'deletion':
                # くの字型の折れ線
                mid_x = x_start + width / 2
                offset = 10
                dwg.add(dwg.polyline(points=[(x_start, y_line), (mid_x, y_line - offset), (x_end, y_line)], fill='none', stroke='black', stroke_width=1, stroke_dasharray="2,2"))
            elif feat.feature_type in ('exon', 'CDS', 'five_prime_UTR', 'three_prime_UTR'):
                base_color = FEATURE_COLORS.get(feat.feature_type, 'gray')
                fill_color = base_color

                # グラデーション設定
                if feat.feature_type in ('exon', 'CDS') and exon_gradation == "on":
                    fill_color = f'url(#{get_or_create_gradient(dwg, base_color, grad_dict)})'
                elif feat.feature_type in ('five_prime_UTR', 'three_prime_UTR') and utr_gradation == "on":
                    fill_color = f'url(#{get_or_create_gradient(dwg, base_color, grad_dict)})'
                
                stroke_color = FEATURE_OUTLINES.get(feat.feature_type, 'black')
                stroke_width = FEATURE_OUTLINE_WIDTHS.get(feat.feature_type, 1)

                if feat is terminal_feature:
                    # 右端が尖った polygon
                    tip = height_feature // 2
                    dwg.add(dwg.polygon(points=[(x_start, y_pos), (x_end - tip, y_pos), (x_end, y_pos + height_feature / 2), (x_end - tip, y_pos + height_feature), (x_start, y_pos + height_feature)], fill=fill_color, stroke=stroke_color, stroke_width=stroke_width))
                else:
                    # 通常の四角
                    dwg.add(dwg.rect(insert=(x_start, y_pos), size=(width, height_feature), fill=fill_color, stroke=stroke_color, stroke_width=stroke_width))

        # === Insertions ===
        triangle_height = 6
        y_triangle = y_pos - 8

        for ins in getattr(gene, "insertions", []):
            if hasattr(ins, 'position'):
                ins_pos = ins.position
                ins_length = getattr(ins, 'length', 1)
                ins_color = getattr(ins, 'color', 'black')
            else:
                ins_pos = ins
                ins_length = 1
                ins_color = "black"

            x = LEFT_MARGIN + (ins_pos - draw_start) / shrink_factor * scale
            base_width = get_insertion_base_width(ins_length, shrink_factor, scale)

            dwg.add(
                dwg.polygon(
                    points=[
                        (x - base_width / 2, y_triangle),
                        (x + base_width / 2, y_triangle),
                        (x, y_triangle + triangle_height)
                    ],
                    fill=ins_color,
                    stroke=ins_color,
                    stroke_width=1.5
                )
            )

        # === SNPs ===
        snp_extend_up = 8
        snp_extend_down = 8
        y_snp_top = y_pos - snp_extend_up
        y_snp_bottom = y_pos + height_feature + snp_extend_down

        for snp in getattr(gene, "snps", []):
            if hasattr(snp, "position"):
                snp_pos = snp.position
                snp_color = getattr(snp, "color", "black")
            else:
                snp_pos = snp
                snp_color = "black"

            x = LEFT_MARGIN + (snp_pos - draw_start) / shrink_factor * scale
            dwg.add(
                dwg.line(
                    start=(x, y_snp_top),
                    end=(x, y_snp_bottom),
                    stroke=snp_color,
                    stroke_width=1.2
                )
            )

        # ドメインを描画（上層）
        for feat in all_features:
            if feat.feature_type == 'domain':
                x_start = LEFT_MARGIN + (feat.start - draw_start) / shrink_factor * scale
                x_end = LEFT_MARGIN + (feat.end - draw_start) / shrink_factor * scale
                domain_color = feat.attributes.get('color', FEATURE_COLORS.get('domain', 'green'))
                if domain_gradation == "on":
                    domain_color = f'url(#{get_or_create_gradient(dwg, domain_color, grad_dict)})'
                dwg.add(dwg.rect(insert=(x_start, y_pos), size=(x_end - x_start, height_feature), fill=domain_color, stroke='black', stroke_width=1))

        # ラベルを遺伝子構造の下に描画（中央揃え）
        if show_labels and gene_center_x is not None:
            dwg.add(dwg.text(label, insert=(gene_center_x, y_pos + height_feature + label_spacing + 10), font_size='10px', fill='black', font_family='monospace', text_anchor='middle'))

    # === 凡例の動的生成 ===
    legend_x = LEFT_MARGIN + axis_width + 50
    legend_y = 30
    box_size, spacing = 12, 20
    present_feature_types = set()
    for gene in genes:
        for f in gene.get_sorted_features(): present_feature_types.add(f.feature_type)
    
    legend_items = []
    if 'CDS' in present_feature_types or 'exon' in present_feature_types: legend_items.append(('CDS', 'Exon/CDS'))
    if 'five_prime_UTR' in present_feature_types: legend_items.append(('five_prime_UTR', "5' UTR"))
    if 'three_prime_UTR' in present_feature_types: legend_items.append(('three_prime_UTR', "3' UTR"))
    if 'intron' in present_feature_types or any(get_baseline_segments(g['start'], g['end'], [f for f in g['gene'].features if f.feature_type == 'deletion']) for g in gene_ranges):
        legend_items.append(('intron', 'Intron'))
    if 'deletion' in present_feature_types: legend_items.append(('deletion', 'Deletion'))
    if any(getattr(g['gene'], "insertions", []) for g in gene_ranges): legend_items.append(('insertion', 'Insertion'))
    if any(getattr(g['gene'], "snps", []) for g in gene_ranges): legend_items.append(('snp', 'SNP'))
    
    all_domain_colors = {}
    for g in genes: all_domain_colors.update(g.domain_color_map)
    for domain_name, color in all_domain_colors.items(): legend_items.append(('domain', domain_name))

    for i, (feat_key, label_text) in enumerate(legend_items):
        y_legend = legend_y + i * spacing
        if feat_key == 'deletion':
            # くの字型
            y_mid = y_legend + box_size // 2
            dwg.add(dwg.polyline(points=[(legend_x, y_mid), (legend_x + box_size // 2, y_mid - 6), (legend_x + box_size, y_mid)], fill="none", stroke="black", stroke_width=1.5, stroke_dasharray="2,2"))
        elif feat_key == 'intron':
            y_line = y_legend + box_size // 2
            dwg.add(dwg.line(start=(legend_x, y_line), end=(legend_x + box_size, y_line), stroke=FEATURE_COLORS.get('intron', 'black'), stroke_width=1))
        else:
            base_color = FEATURE_COLORS.get(feat_key, 'gray')
            if feat_key == 'domain': base_color = all_domain_colors.get(label_text, base_color)
            fill_color = base_color
            dwg.add(dwg.rect(insert=(legend_x, y_legend), size=(box_size, box_size), fill=fill_color, stroke='black'))
        dwg.add(dwg.text(label_text, insert=(legend_x + box_size + 5, y_legend + box_size - 2), font_size='12px', fill='black'))

    dwg.save()
