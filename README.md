![img1](images/gs_logo.png)

![GitHub Release](https://img.shields.io/github/v/release/hashimotoshumpei/GenoSee)
![GitHub Release Date](https://img.shields.io/github/release-date/hashimotoshumpei/GenoSee)
![GitHub last commit](https://img.shields.io/github/last-commit/hashimotoshumpei/GenoSee)
![GitHub License](https://img.shields.io/github/license/hashimotoshumpei/GenoSee)
`<img src="https://img.shields.io/badge/-Python-F9DC3E.svg?logo=python&style=flat">`

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
| `-s`, `--species`      | Name of species to be analyzed [required].                             |
| `-o`, `--output`       | Specifies the file path of the output image file.                      |
| `--fill`                 | Fill between markers [default: on].                                    |
| `--drawing_mode`         | Selects the drawing mode: normal, compare or zoomed [default: normal]. |
| `--coloring_mode`        | Selects the color mode: 2-color or 3-color [default: 2-color]          |
| `--display_marker_names` | Adds marker names at each position [default: on].                      |
| `--chr`                  | Chromosome number of the region to zoom in.                            |
| `--start`                | Start position of the region to zoom in.                               |
| `--end`                  | End position of the region to zoom in.                                 |
| `--dpi`                  | dpi of the output image(s).                                            |

## Examples

### 1. Chromosome-wide graphical genotype

```
python GenoSee.py --input ./examples/Oryza_sativa_phased_100_markers_10_samples.csv --species Oryza_sativa
```

![img3](image/Os_filled_2-color_normal.png)

## Settings

### Color

The lengths of each chromosome for the following species are pre-registered in the database (chromosome_length.json)

* *Oryza sativa*
* *Sorghum bicolor*
* *Zea mays*
* *Triticum aestivum*
* *Hordeum vulgare*
* *Glycine max*
* *Solanum lycopersicum*
* *Arabidopsis thaliana*

### How to add new species ?

You can add chromosome lengths of new species from genome fasta file by using `add_chromosome_lengths.py` in this repo like this;

``python add_chromosome_lengths.py path/to/your/file.fasta.gz path/to/your/chromosome_length_database.json your_species_name``

> [!NOTE]
> BioPython package needs to be installed before use.
>
> `pip install BioPython`

### How to make an input format from a VCF file ?

By using `make_input_file_from_VCF.py`, you can create an input file for GenoSee from a VCF file.

``python make_input_file_from_VCF.py path/to/your/VCF``

> [!NOTE]
> PyVCF package needs to be installed before use.
>
> `pip install PyVCF`

### How to custamize visualization results ?

You can adjust the output diagram by directly modifying the variables within the functions responsible for each drawing mode. Each drawing mode corresponds to the following functions in `plotting.py`.

| Drawing Mode | Responsible Function Name |
| ------------ | ------------------------- |
| normal       | create_normal_plot        |
| compare      | create_comparison_plot    |
| zoomed       | create_zoomed_plot        |

Adjustable variables correspond to the following elements within the diagram.

![img10](image/readme_img_2.png)

### What color sets are used to visualize ?

Following color pallettes are pre-registered in the database (color_set.json)

* normal
  ![img11](image/normal.png)
* grays
  ![img12](image/grays.png)
* reds
  ![img13](image/reds.png)
* blues
  ![img14](image/blues.png)

### How to add new color sets ?

By editing `color_set.json` directly, you can use new color sets for visualization. Please run  `GenoSee.py` with `--color_palette` as an argument like `--color_palette new_color_set` . In this version of GenoSee, The specification of the color set is somewhat redundant. "0|0", "1|1", "0|1" and ".|." are used in the `3-color` mode, whereas "0", "1" and "." are used in the `2-color` mode.

```
"normal": {
        "0|0": "gold",
        "1|1": "dodgerblue",
        "0": "gold",
        "1": "dodgerblue",
        "0|1": "limegreen",
        ".|.": "gray",
        ".": "gray"
}
```

### What is ColabGenoSee ?

![img](image/ColabGenoSee.png)
ColabGenoSee is the Google Colab version of GenoSee.
To use ColabGenoSee, first move the ColabGenoSee folder to your Google Drive, then open ColabGenoSee.ipynb in Google Colab and run it.
