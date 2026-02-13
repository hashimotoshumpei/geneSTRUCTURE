![img1](images/gs_logo.png)

![GitHub Release](https://img.shields.io/github/v/release/hashimotoshumpei/geneSTRUCTURE)
![GitHub Release Date](https://img.shields.io/github/release-date/hashimotoshumpei/geneSTRUCTURE)
![GitHub last commit](https://img.shields.io/github/last-commit/hashimotoshumpei/geneSTRUCTURE)
![GitHub License](https://img.shields.io/github/license/hashimotoshumpei/geneSTRUCTURE)

# A high-quality visualization tool for gene structures
## Citation
Hashimoto, yamada and Izawa. in prepareing.

## Getting Started

### Installation

You can git clone the Github repo and install it locally with the following:

```
git clone https://github.com/hashimotoshumpei/geneSTRUCTURE.git
cd geneSTRUCTURE
python geneSTRUCTURE.py --help
```

### Requirement

* Python3
* svgwrite

### Input data

You can use following csv files.

|  transcript_id  | deletions | insertions |  snps  |                    domains                    |
| :-------------: | :-------: | :--------: | :-----: | :-------------------------------------------: |
| Os06t0160700-01 | 4000-5000 |  100;500  | 300;600 | 1-200:domain1;200-300:domain2;300-400:domain3 |
| Os04t0648800-01 |   1-350   |            |        | 1-100:Kinase;101-200:domain1;201-300:domain2 |

### Options
Here is a list of the arguments that can be used with this tool.

| Flag                       | Description                                                            |
| -------------------------- | ---------------------------------------------------------------------- |
| `-h`, `--help`         | Displays this help message and basic documentation.                    |
| `-i`, `--input`        | Specifies the file path of the input CSV file [required].              |
| `-o`, `--output`       | Specifies the dir path of the output image file.                      |

## Run
```
python geneSTRUCTURE.py --input examples/gene_input.csv --gff examples/gff/transcripts.gff --output ./examples
```


## Examples

### 1. Simple
![simple](examples/Os06t0160700-01.svg)

### 2. Deletion(s)
![simple](examples/Os06t0160700-01_DEL_4000-5000.svg)

### 3. Insertion(s)
![simple](examples/Os06t0160700-01_INS_100_500.svg)

### 4. Domain(s)
![simple](examples/Os06t0160700-01_DOM_1-200_domain1_200-300_domain2_300-400_domain3.svg)

## Settings
For custamizing visualization, please edit config.py.

```
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

```
