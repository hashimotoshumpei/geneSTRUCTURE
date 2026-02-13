# =====================
# 色とグラデーション設定フラグ
# =====================

utr_gradation = "off"
exon_gradation = "off"
domain_gradation = "off"

FEATURE_COLORS = {
    'exon': 'dodgerblue',  # lightblue
    'CDS': 'dodgerblue',
    'five_prime_UTR': 'white',  # orange
    'three_prime_UTR': 'white',  # lightgreen
    'intron': 'black',
    'deletion': 'none',
    'highlight_intron': 'blue',
}

# =====================
# ★描画の色やスタイル設定
# =====================

# FEATURE_COLORS = {
#     'exon': 'lightblue',
#     'CDS': 'lightblue',
#     'five_prime_UTR': 'orange',
#     'three_prime_UTR': 'lightgreen',
#     'intron': 'black',
#     'domain': 'green',
#     'deletion': 'none',  # 塗りつぶしなし（点線表示用）
#     'highlight_intron': 'blue',
# }

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


DOMAIN_COLOR_PALETTE = [
    "lightblue",  # blue
    "lightcoral",  # orange
    "lightgreen",  # green
    "#d62728",  # red
    "#9467bd",  # purple
    "#8c564b",  # brown
    "#e377c2",  # pink
    "#7f7f7f",  # gray
    "#bcbd22",  # olive
    "#17becf",  # cyan
]
