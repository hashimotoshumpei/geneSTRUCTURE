# geneSTRUCTURE

# What is it ?
geneSTRUCTURE is a tool for drawing gene structures from gonome annotation files. Domain information can also be added!
![image](/img/1.png) 

# How to use ?
### Requirments
* Python 3.8 or later
* reportlab 3.6.7
* numpy 1.21.2
### Installation

```bash
pip install reportlab
pip install numpy
```

### Usage
```bash
git clone https://github.com/geneSTRUCTURE/geneSTRUCTURE.git
cd geneSTRUCTURE
python geneSTRUCTURE.py
```

### NOTE
You must run `geneSTRUCTURE.py` with `config.ini` file which includes all settings.  
```bash
[mode_setting]
# Choose "basic" or "domain"
mode = basic

[file_settings]
gff_path = ./example/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.51.gff3
gene_id = SORBI_3007G204600
file_name = bs_220329_h20_my5_domain3_zz.pdf

[domain_settings] 
number_of_domains = 3
# domain positions (in the protein sequence)
domain1_AA_start = 90
domain1_AA_end = 186
domain2_AA_start = 236
domain2_AA_end = 300
domain3_AA_start = 350
domain3_AA_end = 500

[color_settings]
UTR_color = #d3d3d3
Exon_color = #000000
line_color = #000000
domain1_color = #87cefa
domain2_color = #ffb6c1
domain3_color = #98fb98

[drawing_settings]
margin_x = 50
margin_y = 5
gene_h = 20
# zigzag or straight 
intron_shape = zigzag
```
Contents in the config file are shown here. 
![image](/img/2.png) 
