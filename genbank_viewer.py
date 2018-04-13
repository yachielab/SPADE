from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from Bio import Alphabet
from Bio import SeqIO
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import math

#matplotlib.rcParams['pgf.preamble'] = [r'\usepackage{hyperref}']
template = {'legend.numpoints': 1, 'axes.axisbelow': True, 'font.sans-serif': [u'Arial', u'Liberation Sans', u'Bitstream Vera Sans', u'sans-serif'], 'axes.labelcolor': '.15', 'ytick.major.size': 0.05, 'axes.grid': False, 'ytick.minor.size': 0.0, 'legend.scatterpoints': 1, 'axes.edgecolor': "black", 'grid.color': 'white', 'legend.frameon': False, 'ytick.color': '.15', 'xtick.major.size': 2.0, 'figure.facecolor': "#EAEAF2", 'xtick.color': '.15', 'xtick.minor.size': 3.0, 'font.family': [u'sans-serif'], 'xtick.direction': u'out', 'lines.solid_capstyle': u'round', 'grid.linestyle': u'-', 'image.cmap': u'Greys', 'axes.facecolor': "white", 'text.color': '.15', 'ytick.direction': u'out', 'axes.linewidth': 0.5}
sns.set_context("poster", font_scale=1.0, rc={"lines.linewidth": 2.0}) 
#sns.set_style(template) 

class GenbView(object):
    def __init__(self,gb,fig):
        self.gb    = gb
        self.ax    = None
        self.fig   = fig
        self.width = 0.025  
        self.url_list = [] 
    
    def make_parameter_list(self, start_pos, end_pos): 
        self.arrows_parameters = []
        self.start_pos = start_pos
        self.end_pos   = end_pos 
        partial_gb     = self.extract(start_pos, end_pos, make_file = False)
        self.length    = end_pos - start_pos
        self.hl        = (end_pos-start_pos) * 0.01
        ypos  = 1.00  # height of arrow positio
        num   = 0 
        e_num = 0
        box   = [[0] * 10001]
        pre_box_start = 0 
        pre_box_end   = -1 
        for i, afeat in enumerate(partial_gb.features):
            parameter_dict = {} 
            if afeat.location.start < start_pos:
                start = start_pos   
            elif "region" in afeat.qualifiers.keys() and len(afeat.qualifiers["region"]) == 1:
                start = int(afeat.qualifiers["region"][0].split(":")[0])
            else:
                start = afeat.location.start

            if afeat.location.end > end_pos:
                end = end_pos 
            elif "region" in afeat.qualifiers.keys() and len(afeat.qualifiers["region"]) == 1:
                end = int(afeat.qualifiers["region"][0].split(":")[1]) 
            else:
                end = afeat.location.end
            
            if afeat.type != "gene" and afeat.type != "mRNA" and afeat.type != "exon" and afeat.type != "source":
                if afeat.strand == 1:
                    xpos  = start - start_pos
                    tdx   = end - start
                    dx    = end - start - self.hl
                    if afeat.type == "CDS": 
                        edgecolor = "#C481A1"
                        direction = "right"
                    else:
                        edgecolor = "#C9BD74"
                        direction = "right"
                    if dx < 0:
                        dx = self.hl
                else: 
                    xpos = end - start_pos
                    tdx  = start - end
                    dx   = start - end + self.hl
                    if afeat.type == "CDS":
                        edgecolor = "#74AFB9"
                        direction = "left"
                    else:
                        edgecolor = "#C9BD74"
                        direction = "left"
                    if dx > 0:
                        dx = -1.0 * self.hl
                if "note" in afeat.qualifiers.keys() and "SPADE" in afeat.qualifiers["note"][0]:
                    pass 
                else:
                    if "old_locus_tag" in afeat.qualifiers.keys():
                        locus_tag = afeat.qualifiers["old_locus_tag"][0]
                    elif "locus_tag" in afeat.qualifiers.keys(): 
                        locus_tag = afeat.qualifiers["locus_tag"][0]
                    else:
                        locus_tag = "-"

                    if "product" in afeat.qualifiers.keys() and 4.5 * len(afeat.qualifiers["product"][0]) < 346 * (1.0 * abs(dx) / self.length):
                        txt_annotation = afeat.qualifiers["product"][0] 
                    elif "gene" in afeat.qualifiers.keys() and 4.5 * len(afeat.qualifiers["gene"][0]) < 346 * (1.0 * abs(dx) / self.length):
                        txt_annotation = afeat.qualifiers["gene"][0] 
                    if "locus_tag" in afeat.qualifiers.keys() and 4.5 * len(afeat.qualifiers["locus_tag"][0]) < 346 * (1.0 * abs(dx) / self.length): 
                        txt_annotation = afeat.qualifiers["locus_tag"][0]
                    elif "old_locus_tag" in afeat.qualifiers.keys() and 4.5 * len(afeat.qualifiers["old_locus_tag"][0]) < 346 * (1.0 * abs(dx) / self.length):
                        txt_annotation = afeat.qualifiers["old_locus_tag"][0]
                    elif "rpt_family" in afeat.qualifiers.keys() and 4.5 * len(afeat.qualifiers["rpt_family"][0]) < 346 * (1.0 * abs(dx) / self.length):
                        txt_annotation = afeat.qualifiers["rpt_family"][0]
                    else:
                        txt_annotation = afeat.type
                    
                    
                    hl = self.hl
                    #if abs(dx) + abs(self.hl) == end_pos - start_pos:
                    #    dx = 10 * dx
                    #    if afeat.strand == 1:
                    #        xpos = xpos-100
                    #    if afeat.strand == -1: 
                    #        xpos = xpos+100
                    
                    #if afeat.strand == 1 and xpos == 0:   
                    #    xpos = xpos-100
                    #    dx   = dx+100
                    
                    #if afeat.strand == -1 and xpos == end_pos-start_pos: 
                    #    xpos = xpos+100
                    #    dx   = dx-100
                   
                    if afeat.strand == -1 and xpos + dx - self.hl  <= 0.1:
                        dx   = dx-self.hl
                        hl = 0 
                    parameter_dict["start"]      = hl
                    parameter_dict["start"]      = start
                    parameter_dict["end"]        = end 
                    parameter_dict["xpos"]       = xpos
                    parameter_dict["dx"]         = dx
                    parameter_dict["edgecolor"]  = edgecolor
                    parameter_dict["direction"]  = direction 
                    parameter_dict["facecolor"]  = "#FFFFEF"
                    parameter_dict["annotation"] = txt_annotation
                    flag = 0
                    box_start = int(10000 * (xpos) * 1.0/self.length)
                    box_end   = int(math.ceil(10000 * (xpos+tdx) * 1.0/self.length))
                    box_start,box_end = sorted([box_start,box_end]) 
                    for j in range(len(box)):
                        if box[j][box_start] == 0:
                            box[j][box_start:box_end] = [1] * (box_end - box_start)
                            ypos = 1.0 - 0.03 * j
                            flag = 1
                            break 
                        else:
                            pass 
                    if flag == 0:
                        ypos = 1.0 - 0.03 * (j + 1)  
                        box.append([0] * 10001)
                        box[-1][box_end:box_start] = [1] * (box_end - box_start)
                    pre_box_start = box_start
                    pre_box_end   = box_end
                    parameter_dict["ypos"] = ypos
                    self.arrows_parameters.append(parameter_dict)
        self.box_num = len(box)
        if sum(box[0]) == 0:
            return 0 
        return len(box)                 
    
    def view_feature(self, ax):
        self.ax = ax
        url_list = [] 
        for n, parameters in enumerate(self.arrows_parameters):
            arrow = self.ax.arrow(parameters["xpos"], parameters["ypos"], parameters["dx"], 0, width=self.width, head_length=self.hl, head_width=self.width, facecolor="#FFFFEF", edgecolor=parameters["edgecolor"], lw=1.5)
            arrow.set_url("Annotation_" + str(n))
            url_list.append("https://www.ncbi.nlm.nih.gov/nuccore/%s?from=%d&to=%d" % (self.gb.id,parameters["start"],parameters["end"]))
            if 4.5 * len(parameters["annotation"]) < 346 * (1.0 * abs(parameters["dx"]) / self.length):
                pass 
                self.ax.text(parameters["xpos"] + 0.5 * parameters["dx"], parameters["ypos"], "%s" % (parameters["annotation"]), horizontalalignment='center', verticalalignment='center', fontsize=6, family="Courier New", weight="bold")        
        self.ax.set_xlim(0,self.length) 
        self.ax.set_ylim(1.01-0.03 * self.box_num,1.02) 
        self.ax.spines["top"].set_color("none")
        self.ax.spines["bottom"].set_color("none")
        self.ax.spines["right"].set_color("none")
        self.ax.spines["left"].set_color("none")
        self.ax.yaxis.set_ticks([])
        self.ax.yaxis.set_ticklabels([])
        self.ax.xaxis.set_ticks([])
        self.ax.xaxis.set_ticklabels([])
        url_output = open("annotaiton_url_list.txt","w")
        for url in url_list:
            url_output.write(url + "\n")
    
    def extract(self, start_pos, end_pos, make_file=False): 
        range_set = set(range(start_pos, end_pos)) 
        partial_gb = SeqRecord(Seq(str(self.gb.seq[start_pos:end_pos]),Alphabet.DNAAlphabet())) 
        for afeat in self.gb.features:
            afeat_range = set(range(afeat.location.start, afeat.location.end))
            if len(afeat_range & range_set) > 0:
                partial_gb.features.append(afeat) 
        if make_file == True:
            record_handle = open(partial_gb.id + "_" + str(start_pos) + "_" + str(end_pos),"w")
            SeqIO.write(partial_gb, record_handle, "genbank") 
        return partial_gb
              

if __name__ == "__main__":
    fig,ax     = plt.subplots(figsize=(18,2))
    gb  = SeqIO.read(sys.argv[1],"genbank")
    if set(str(gb.seq)) == set(["N"]):
        fasta = SeqIO.read(sys.argv[2],"fasta")
        gb.seq = fasta.seq
    else:
        pass
    SeqIO.read(sys.argv[1],"genbank")
    GV = GenbView(gb,ax,fig) 
    GV.view_feature(133000,136000,feature_name="all")
    fig.savefig("test_gv.pdf")


