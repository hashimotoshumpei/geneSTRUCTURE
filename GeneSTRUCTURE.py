
from reportlab.pdfgen import canvas
import re
import sys
import numpy as np
import configparser
from logging import getLogger, StreamHandler, FileHandler, DEBUG, INFO, WARNING, Formatter

############################### set logger ##################################
logger = getLogger(__name__)
logger.setLevel(DEBUG)

sh = StreamHandler()
sh.setLevel(INFO)
sh.setFormatter(Formatter("%(asctime)s %(levelname)8s %(message)s"))

fh = FileHandler(filename = 'geneSTRUCTURE.log', mode = 'w')
fh.setLevel(DEBUG)
fh.setFormatter(Formatter("%(asctime)s %(levelname)8s %(message)s"))

logger.addHandler(sh)
logger.addHandler(fh)
#---------------------------------------------------------------------------END#
############################### Functions ##################################
def get_transcript_id(gff_path, gene_id):
    """_summary_

    Args:
        gff_path (_type_): _description_
        id (_type_): _description_

    Returns:
        _type_: _description_
    """
    if gff_path.endswith('.gff3'):
        transcript_id = []
        strand = []
        with open(gff_path, mode = 'r') as inp:
            for i, line in enumerate(inp):
                if line.startswith('#'):
                    pass
                else:
                    line = line.split('\t')
                    id_col = re.split('[:;]', line[8])

                    if len(id_col) < 4: continue
                    # Use "continue", not "pass".

                    elif id_col[3] == gene_id:

                        transcript_id = id_col[1]
                        strand = line[6]

        return transcript_id, strand
    
    else:
        return False

                
def get_structure(transcript_id):
    """_summary_

    Args:
        transcript_id (_type_): _description_

    Returns:
        _type_: _description_
    """

    exon_pos = []
    exon_len = []
    total_length = 0
    five_prime_UTR = 0
    three_prime_UTR = 0

    with open(gff_path, mode = 'r') as inp: 
        for line in inp:
            line = line.rstrip('\r|\n|\r\n')
            if line.startswith('#'):
                pass
            else:
                line = line.split('\t')
                id_col = re.split('[:;]', line[8])
                if id_col[1] == transcript_id:

                    if line[2] == 'mRNA':
                        total_length = (int(line[4]) - int(line[3]))/10
                    elif line[2] == 'exon':
                        exon_pos.append(int(line[3]))
                        exon_pos.append(int(line[4]))
                        exon_len.append((int(line[4]) - int(line[3]))/10)
                    elif line[2] == 'five_prime_UTR':
                        five_prime_UTR = (int(line[4]) - int(line[3]))/10
                    elif line[2] == 'three_prime_UTR':
                        three_prime_UTR = (int(line[4]) - int(line[3]))/10
    
    return total_length, exon_pos, five_prime_UTR, three_prime_UTR

def cDNA_pos2gDNA_pos(cDNA_exon_pos, domain_cDNA_pos):
    """_summary_

    Args:
        cDNA_exon_pos (_type_): _description_
        domain_cDNA_pos (_type_): _description_

    Returns:
        _type_: _description_
    """

    x2 = np.sort(np.append(cDNA_exon_pos, domain_cDNA_pos))
    index = int(np.where(x2 == domain_cDNA_pos)[0])
    gDNA_pos = domain_cDNA_pos + cumsum_intron_len[index-1]

    return gDNA_pos


def color_convert(color16):
    """_summary_

    Args:
        color16 (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    color = color16.lstrip('#')
    color = list(color)

    lib = {'A':10, 'B':11, 'C':12, 'D':13, 'E':14, 'F':15,
            'a':10, 'b':11, 'c':12, 'd':13, 'e':14, 'f':15}

    for i,j in lib.items():
        for k,l in enumerate(color):
            if l == i:
                color[k] = j

    color = np.array([int(i) for i in color]).reshape(-1,2)

    RGB_255 = color @ np.array([16**1, 16**0]).T
    RGB_1 = RGB_255/255

    return RGB_1
    
#---------------------------------------------------------------------------END#
############################### config setting ##################################
inifile = configparser.ConfigParser()
inifile.read('./config.ini')

mode = inifile.get('mode_setting', 'mode')

# file settings
gff_path = inifile.get('file_settings', 'gff_path')
gene_id = inifile.get('file_settings', 'gene_id')
file_name = inifile.get('file_settings', 'file_name')

# color setting
UTR_color = color_convert(inifile.get('color_settings', 'UTR_color'))
Exon_color = color_convert(inifile.get('color_settings', 'Exon_color'))
line_color = color_convert(inifile.get('color_settings', 'line_color'))

# drawing settings
margin_x = int(inifile.get('drawing_settings', 'margin_x'))
margin_y = int(inifile.get('drawing_settings', 'margin_y'))
intron_shape = inifile.get('drawing_settings', 'intron_shape')
gene_h = int(inifile.get('drawing_settings', 'gene_h'))
#---------------------------------------------------------------------------END#

##################################### main script ############################################ 

# Get transcript ID
if not get_transcript_id(gff_path, gene_id):
    logger.info('Invalid input file format. Only .gff3 file is acceptable.')
    sys.exit()

else: transcript_id, strand = get_transcript_id(gff_path, gene_id)

if transcript_id == []:
    logger.info(f'Gene ID "{gene_id}" was not found.')
    sys.exit()

else:
    logger.info(f'Gene ID "{gene_id}" was found.')
    logger.info(f'Gene structure will be made for transcript id "{transcript_id}"')
    total_length, exon_pos, five_prime_UTR, three_prime_UTR = get_structure(transcript_id)
'''NOTE

'''

if five_prime_UTR == 0:
    logger.info("There was no annotation for 5'UTR")

if three_prime_UTR == 0:
    logger.info("There was no annotation for 3'UTR")

#エキソン＆イントロン長を取得
if strand == '-':
    exon_intron_length = np.asarray([(exon_pos[i+1] - exon_pos[i])/10 for i in range(len(exon_pos)-1)])[::-1]
else: exon_intron_length = np.asarray([(exon_pos[i+1] - exon_pos[i])/10 for i in range(len(exon_pos)-1)])

exon_len = np.asarray([exon_intron_length[i] for i in range(len(exon_intron_length)) if i % 2 == 0])
intron_len = np.asarray([exon_intron_length[i] for i in range(len(exon_intron_length)) if i % 2 == 1])

# 累積イントロン長
cumsum_intron_len = np.append(np.append(0, np.cumsum(intron_len)), 0)
#cDNAでの位置
cDNA_exon_pos = np.append(0, np.cumsum(exon_len))

# 基本の作図用
if strand == '-':
    x = np.abs(exon_pos - np.max(exon_pos))[::-1]/10 + 50
else: x  = (exon_pos - np.min(exon_pos))/10 + 50


############################### printing setting ##################################

pagesize_w = total_length + margin_x * 2
pagesize_h = gene_h + margin_y * 2
center_line_y = margin_y + gene_h/2

page = canvas.Canvas(file_name, pagesize=(pagesize_w,pagesize_h))
page.setStrokeColorRGB(line_color[0], line_color[1], line_color[2])
page.setFillColorRGB(Exon_color[0], Exon_color[1], Exon_color[2])
# 線の太さを変更
page.setLineWidth(1)

# 線の描画
for i in range(0,len(x),2):
    page.rect(x[i],margin_y,exon_intron_length[i], gene_h, fill=True)

if intron_shape == 'zigzag':

    mid_intron = []

    for i in range(1,len(x)-1,2):
        mid = (x[i] + x[i+1])/2
        mid_intron.append(mid)

    for i,j in enumerate(range(1,len(x)-1,2)):
        page.line(x[j], center_line_y, mid_intron[i], 0)
        page.line(mid_intron[i], 0, x[j+1], center_line_y)
    
elif intron_shape == 'straight':
    for i in range(1,len(x)-1,2):
        page.line(x[i], center_line_y, x[i+1], center_line_y)

else:
    print('Invalid intron shape')

page.setFillColorRGB(UTR_color[0], UTR_color[1], UTR_color[2])

if five_prime_UTR != 0:
    page.rect(x[0], margin_y, five_prime_UTR, gene_h, fill = True)

if three_prime_UTR != 0:
    page.rect(x[-1]-three_prime_UTR, margin_y, three_prime_UTR, gene_h, fill = True)





            
############################### Domain mode ##################################
if mode == 'domain':

    n_domain = int(inifile.get('domain_settings', 'number_of_domains'))

    for i in range(n_domain):
        AA_start = int(inifile.get('domain_settings', f'domain{i+1}_AA_start'))
        AA_end = int(inifile.get('domain_settings', f'domain{i+1}_AA_end'))
        color = color_convert(inifile.get('color_settings', f'domain{i+1}_color'))

        cDNA_start = (AA_start * 3)/10 + five_prime_UTR
        cDNA_end = (AA_end * 3)/10 + five_prime_UTR

        gDNA_start = cDNA_pos2gDNA_pos(cDNA_exon_pos, cDNA_start) + 50 
        gDNA_end = cDNA_pos2gDNA_pos(cDNA_exon_pos, cDNA_end) + 50

        if gDNA_end > x[-1]:
            logger.warning(f'The end position of domain{i+1} is out of range.') 
        # Even in this case, the end point of the domain is as same as that of codeing region.

        domain_pos = x[(gDNA_start <= x) & (x <= gDNA_end)]
        domain_pos = np.sort(np.append([gDNA_start, gDNA_end], domain_pos))
        domain_len = np.asarray([(domain_pos[i+1] - domain_pos[i]) for i in range(len(domain_pos)-1)])

        page.setFillColorRGB(color[0], color[1], color[2])

        for j in range(0,len(domain_len),2):
            page.rect(domain_pos[j],margin_y,domain_len[j], gene_h, fill=True)
#---------------------------------------------------------------------------END#

page.save() 
logger.info(f'Gene structure was successfully saved as "{file_name}"')

#---------------------------------------------------------------------------END#



