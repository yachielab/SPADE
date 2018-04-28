#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys 
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
import shutil
import genbank_viewer as gv
from Bio import SeqIO
from weblogolib import *
import warnings
#warnings.filterwarnings('ignore')

template = {'legend.numpoints': 1, 'axes.axisbelow': True, 'axes.labelcolor': '.15', 'ytick.major.width': 1.0, 'ytick.major.size': 4.0, 'axes.grid': False, 'ytick.minor.size': 0.0, 'legend.scatterpoints': 1, 'axes.edgecolor': "black", 'grid.color': 'white', 'legend.frameon': False, 'ytick.color': '.15', 'xtick.major.size': 4.0, 'xtick.major.width': 1.0, 'figure.facecolor': "#EAEAF2", 'xtick.color': '.15', 'xtick.minor.size': 3.0, 'xtick.direction': u'out', 'lines.solid_capstyle': u'round', 'grid.linestyle': u'-', 'image.cmap': u'Greys', 'axes.facecolor': "white", 'text.color': '.15', 'ytick.direction': u'out', 'axes.linewidth': 1.0}
fonts = [font.split("/")[-1] for font in fm.findSystemFonts()]
if "Helvetica.ttf" in fonts:
    sns.set(font = "Helvetica")
else:
    sns.set(font = "Arial")

sns.set_context("poster", font_scale=1.4, rc={"lines.linewidth": 1.0}) 
sns.set_style(template)

def save_data(matrix, kmer_count_signal, repeat_unit_array, gap_repeat_unit_array, matrix_sumx, periods, candidate, k=10, thresh=20, peak_period_list=[],xlabel_a="Position (bp)", fig_name="hoge", key="", gb="None", dtype="DNA", p_start=0, p_end=0,strand=0):
    peak_matrix_list = []   
    matrix_sumx_list = []

    for i, peak in enumerate(peak_period_list):
        p_min = peak-5 if peak-5 > 0 else 0
        p_max = peak+5 if peak+5 < 1000 else 1000
        peak_matrix_list.append(matrix[p_min:p_max,:]) 
        matrix_sumx_list.append(matrix_sumx[p_min:p_max])
        np.savetxt("peak_matrix_" + str(peak) + ".txt",matrix[p_min:p_max,:].astype(np.int),delimiter="\t", fmt="%.0f")
        np.savetxt("rpt_unit_signal_" + str(peak) + ".txt",repeat_unit_array[i],delimiter="\t", fmt="%.0f")
        np.savetxt("avg_intensity_" + str(peak) + ".txt",matrix_sumx[p_min:p_max],delimiter="\t", fmt="%.0f") 
     
    np.savetxt("kmer_count_signal.txt",kmer_count_signal,delimiter="\t", fmt="%.0f") 
    np.savetxt("periods.txt",periods,delimiter="\t", fmt="%.0f") 

    return peak_matrix_list, matrix_sumx_list

def make_figure(peak_matrix_list, kmer_count_signal, repeat_unit_array, gap_repeat_unit_array, aln_repeat_unit_array, unit_length_list, matrix_sumx_list, periods, candidate, k=10, thresh=20, peak_period_list=[], xlabel_a="Position (bp)", fig_name="hoge", key="", gb="None",  dtype="DNA", p_start=0, p_end=0, Format="pdf", html_path="../../../../script/html", strand=0):
    if dtype == "protein": 
        xlabel_a = "Position (aa)"
    start, end     = candidate[2],candidate[3]
    r_start, r_end = candidate[0], candidate[1] 
    
    #data sorting with max average intensity
    zip_set = list(zip(peak_matrix_list, repeat_unit_array, matrix_sumx_list, peak_period_list))
    zip_set.sort(key = lambda x: np.max(x[2]))
    peak_matrix_list, repeat_unit_array, matrix_sumx_list, peak_period_list = list(zip(*zip_set)) 
    
    #create figure opbject
    fig = plt.figure(figsize=(9,10))

    #make ax1, It is for K-mer cumulative count.
    pos_data = np.array(range(len(kmer_count_signal))) 
    ax1 = plt.axes([0.2, 0.805, 0.55, 0.10])
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    ax1.fill_between(pos_data,0,kmer_count_signal,facecolor="black")
    ax1.set_ylim(0,1.1 * max(kmer_count_signal))
    ax1.tick_params(axis='x', which='major')
    if start >= 0 and end >= 0:
        ticks = [i for i in range(int(math.floor(start/10.0) * 10.0),int(math.ceil(end/10.0) * 10.0),10)]
        hl    = len(kmer_count_signal) / 2.0
        space = round(hl * 1.0 / math.pow(10,int(math.log10(hl)))) * math.pow(10,int(math.log10(hl))) 
        for tick in ticks:
            if (tick-ticks[0]) * 1.0 / (2*hl) >= 0.15 and tick % math.pow(10,int(math.log10(hl))) == 0:
                if end > 1.0e+5:
                    tick_labels = [tick, tick + space]
                else: 
                    tick_labels = [tick-0.5*space, tick, tick + 0.5*space, tick + space, tick + 1.5 * space] 
                    if tick-0.5*space < start:
                        tick_labels.remove(tick-0.5*space) 
                    if tick+1.5*space > end:
                        tick_labels.remove(tick+1.5*space) 
                break
         
        ax1.xaxis.set_ticks(list(map(lambda x: x-start, tick_labels)))
        ax1_xticklocs = ax1.xaxis.get_ticklocs()
        ax1.xaxis.set_ticklabels(list(map(lambda x: format(int(x), ',') if x != "" else "", tick_labels)), fontsize=22) 
    
    ax1.tick_params(top=False,right=False,pad=2)
    ax1_yticklocs = ax1.yaxis.get_ticklocs()
    if len(ax1_yticklocs) > 5:
        ax1_yticklocs = [int(ax1_yticklocs[i]) for i in range(len(ax1_yticklocs)) if i%2 == 0]
        ax1.yaxis.set_ticks(ax1_yticklocs) 
    ax1.set_xlabel(xlabel_a, fontsize=22, labelpad=2) 
    ax1.set_ylabel("Cumulative\n%d-mer count" % (k), fontsize=22) 
    ax1.yaxis.set_ticklabels(list(map(str,ax1_yticklocs)),fontsize=22) 
    ax1.yaxis.set_ticklabels(list(map(lambda x: str(int(x)), ax1_yticklocs)), fontsize=22) 
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    if type(gb) != "None":
        GV = gv.GenbView(gb,fig)
        if dtype == "nucl":
            box_num = GV.make_parameter_list(0,len(kmer_count_signal))
            ax_gv = plt.axes([0.2,0.72-0.0175*(box_num-1),0.55,0.0175*box_num])
            GV.view_feature(ax_gv)
        elif dtype == "prot":
            box_num = 0
    ax1.set_xlim(0,len(pos_data))
    
    start_y = 0.695 - 0.0160 * box_num 
    vmin = 0
    vmax = np.max(list(map(np.max,peak_matrix_list))) 
    new_periods = periods 
    for i, peak in enumerate(peak_period_list):
        matrix      = peak_matrix_list[i]
        matrix_sumx = matrix_sumx_list[i]

        p_min = peak-5 if peak-5 > 0 else 0
        p_max = peak+5 if peak+5 < 1000 else 1000

        axn1 = plt.axes([0.2,start_y-i*0.13-0.09,0.55,0.08])
        img  = axn1.imshow(matrix, extent=[0, matrix.shape[1], p_max, p_min], aspect="auto", interpolation="nearest", cmap="YlGnBu_r", vmin=vmin, vmax=vmax)  
        axn1.set_xlim(0,matrix.shape[1])
        axn1.set_ylabel("Period",fontsize=22)
        axn1_yticklocs = axn1.yaxis.get_ticklocs()
        if len(axn1_yticklocs) > 4:
            axn1_yticklocs = [int(j)  for j in axn1_yticklocs if j%4 == 0]
            axn1.yaxis.set_ticks(axn1_yticklocs) 
        axn1.yaxis.set_ticklabels(list(map(lambda x: str(int(x)), axn1_yticklocs)), fontsize=22) 
        axn1.set_ylim(p_max,p_min)
        axn1.tick_params(top=False)
        axn1.tick_params(right=False) 
        axn1.tick_params(pad=4) 

        if i == 0:
            axn2 = plt.axes([0.2, start_y-i*0.13 + 0.02, 0.55, 0.02])
        else:
            axn2 = plt.axes([0.2, start_y-i*0.13, 0.55, 0.02])
                
        cmap = matplotlib.cm.spring
        norm = matplotlib.colors.Normalize(vmin=0.0,vmax=1.0)
        axn2.patch.set_facecolor('#FFFFF5')
        ul = unit_length_list[i]
        for unitx, unity in repeat_unit_array[i]:
            for x,y in zip(unitx, unity):
                axn2.fill_between([x,x+1],[y,y],0,color=cmap(norm((unitx[-1]-x)*1.0/ul))) 
        
        x_array = []
        y_array = [] 
        for unitx, unity in gap_repeat_unit_array[i]:
            axn2.fill_between(unitx,unity,0,color="#FFFFF5",lw=0.5) 
            x_array.extend(unitx) 
            y_array.extend(unity) 
        
        x_array = []
        y_array = [] 
        for unitx, unity in aln_repeat_unit_array[i]:
            axn2.add_patch(matplotlib.patches.Rectangle((unitx[0], 0.04), unitx[-1]-unitx[0], 0.96, alpha=1, facecolor='none', lw=0.5,edgecolor="k"))
            x_array.extend(unitx) 
            y_array.extend(unity) 
        
        x_array2 = [0] + x_array + [matrix.shape[1]] 
        y_array2 = [0] + y_array + [0]
        axn2.fill_between(x_array2,y_array2,1.05,color="#FFFFF5",lw=0)
        
        axn2.set_ylim(0.025,1.05)
        axn2.set_xlim(0,matrix.shape[1]) 
        axn2.tick_params(top=False)
        axn2.tick_params(right=False)
        axn2.yaxis.set_label_coords(-0.09,0.5)
        axn2.spines['right'].set_visible(False)
        axn2.spines['left'].set_visible(False)
        axn2.spines['top'].set_visible(False)
        axn2.spines['bottom'].set_visible(False)
        axn2.yaxis.set_ticks([])
        axn2.yaxis.set_ticklabels([])
        axn2.xaxis.set_ticks([])
        axn2.xaxis.set_ticklabels([])       
        axn2.tick_params(pad=4) 
       
        axn3 = plt.axes([0.760,start_y-i*0.13-0.09,0.1,0.08])
        axn3.barh(new_periods[p_min:p_max+1],matrix_sumx,facecolor="lightgrey",align="center",linewidth=0.5,edgecolor="k")
        axn3.set_xlim(0,np.max(matrix_sumx) * 1.1)
        axn3.set_ylim(p_max,p_min)
        axn3.tick_params(top=False)
        axn3.tick_params(right=False) 
        axn3.tick_params(left=False)
        axn3.tick_params(bottom=False)        
        axn3.yaxis.set_ticks_position('right')
        axn3.yaxis.set_ticks([peak])
        axn3.yaxis.set_ticklabels(["%d" %(peak)], fontsize=22)
        axn3.tick_params(pad=4)

        if i == len(peak_period_list) - 1:
            axn1.set_xlabel(xlabel_a, fontsize=22)
            #axn1.xaxis.set_ticks(ax1_xticklocs)
            axn1.xaxis.set_ticks([])
            axn1.xaxis.set_ticklabels([]) 
            axn3_xticklocs = axn3.xaxis.get_ticklocs()
            axn3.xaxis.set_ticks([]) 
            axn3.xaxis.set_ticklabels([]) 
            axn3.tick_params(bottom=True,pad=6)        
            axn3.set_xlabel("Average\n Intensity (a.u.)",fontsize=22) 
            cbar_ax = fig.add_axes([0.120, start_y - 0.120  - i * 0.13, 0.10, 0.018])
            fig.colorbar(img, cax=cbar_ax, orientation="horizontal", ticks=[0, int(np.max(matrix))])
            cbar_ax.xaxis.set_ticks([]) 
            cbar_ax.xaxis.set_ticklabels([]) 
            cbar_ax.set_xlabel("Intensity (a.u.)", fontsize=22, labelpad=8)
            cbar_ax.tick_params(pad=6)
        else:
            axn1.xaxis.set_ticks([])
            axn1.xaxis.set_ticklabels([])
            axn3.xaxis.set_ticks([])
            axn3.xaxis.set_ticklabels([])
     
    ax3_title = plt.axes([0.2, start_y + 0.4, 0.55, 0.03])
    ax3_title.text(0.5,0.47, "Periodic mortif area", horizontalalignment='center',verticalalignment='center', fontsize=22)
    ax3_title.xaxis.set_ticks([])
    ax3_title.xaxis.set_ticklabels([])
    ax3_title.yaxis.set_ticks([])
    ax3_title.yaxis.set_ticklabels([])
    ax3_title.spines['right'].set_visible(False)
    ax3_title.spines['top'].set_visible(False)
    ax3_title.spines['left'].set_visible(False)
    ax3_title.spines['bottom'].set_visible(False)

    axn1_title = plt.axes([0.2, start_y, 0.55, 0.02])
    axn1_title.text(0.5,0.15,"Position-Period Matrix", horizontalalignment='center',verticalalignment='center', fontsize=22)
    axn1_title.xaxis.set_ticks([])
    axn1_title.xaxis.set_ticklabels([])
    axn1_title.yaxis.set_ticks([])
    axn1_title.yaxis.set_ticklabels([])
    axn1_title.spines['right'].set_visible(False)
    axn1_title.spines['top'].set_visible(False)
    axn1_title.spines['left'].set_visible(False)
    axn1_title.spines['bottom'].set_visible(False) 

    ax7 = plt.axes([0.2,0.21,0.55,0.73])
    ax7.set_xlim(0, len(pos_data))
    ax7.set_ylim(0.17,0.90)
    ax7.patch.set_alpha(0.0)
    ax7.xaxis.set_ticks([])
    ax7.xaxis.set_ticklabels([])
    ax7.yaxis.set_ticks([])
    ax7.yaxis.set_ticklabels([])
    ax7.spines['right'].set_visible(False)
    ax7.spines['top'].set_visible(False)
    ax7.spines['left'].set_visible(False)
    ax7.spines['bottom'].set_visible(False)
    ax7.tick_params(pad=4)
    fig.patch.set_alpha(0.0)  
    
    rptunit_ax  = fig.add_axes([0.760, 0.7-0.0160 * (box_num-1), 0.1, 0.018])
    cmap        = matplotlib.cm.spring_r
    cmaplist    = [cmap(i) for i in range(256)]
    bounds      = np.linspace(0,255,256)
    norm        = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    cmap        = cmap.from_list('Custom cmap', cmaplist, 256)
    cb          = matplotlib.colorbar.ColorbarBase(rptunit_ax, cmap=cmap, norm=norm, spacing='proportional', ticks=[], boundaries=bounds, orientation='horizontal')
     
    rptunit_ax.tick_params(pad=20,) 
    cb.ax.text(-0.3,-1.32,"Repeat unit",fontsize=22) 
    cb.outline.set_linewidth(0.5)
    
    if Format == "svg":
        path = os.getcwd() 
        shutil.copytree(html_path, path + "/html")
        os.chdir("html") 
        plt.savefig(fig_name + "." + Format, transparent=False)
        svg       = open(fig_name + "." + Format) 
        url_list  = [url.rstrip() for url in open("../annotaiton_url_list.txt","r")]
        new_svg   = set_url(svg,url_list)
        html_temp = open("template.html").read()
        new_html  = html_temp % (fig_name,new_svg)
        open(fig_name + ".html","w").write(new_html)  
    else:
        plt.savefig(fig_name + "." + Format, transparent=False)

    plt.close()

def set_url(svg,url_list):
    url_num = 0
    new_svg = ""
    for line in svg: 
        if "Annotation_" in line:
            new_svg += '<a href="%s" target="_blank" class="livepreview" data-position="bottom">\n' % url_list[url_num]
            url_num += 1
        else:
            new_svg += line
    return new_svg

def motif_logo(dtype):
    fasta_name = "align.unit_seq.fasta"
    fasta = open(fasta_name)
    seqs  = read_seq_data(fasta)
    data  = LogoData.from_seqs(seqs) 
    if dtype == "nucl":
        data.alphabet = std_alphabets["dna"] 
        options = LogoOptions()
        options.logo_font = "Arial_Bold"
        options.number_fontsize    = 34
        options.fontsize           = 34
        options.logo_margin        = 30
        options.stacks_per_line    = 150
        options.scale_width        = 0
        options.stack_width        = 15
        options.stack_aspect_ratio = 5
        options.stroke_width       = 1.5
        options.tic_length         = 7
        options.fineprint          = ""
        options.color_scheme = std_color_schemes["classic"]
        options.yaxis_label        = "Bits"
        options.yaxis_tic_interval = 1
        options.yaxis_minor_tic_ratio = 2   
        format = LogoFormat(data, options)
        fout = open("weblogo.pdf","wb")
        try:
            fout = open("weblogo.pdf","wb")
            fout.write(pdf_formatter(data, format))
            fout.close() 
        except:
            fout.close()
            warnings.warn("Ghostscript was not found. If you want to generate LOGO in pdf format, please install Ghostscript.\n For the present, LOGO was output in EPS format") 
            os.system("rm -rf weblogo.pdf") 
            fout = open("weblogo.eps","wb")
            fout.write(eps_formatter(data, format))
            fout.close() 
        try:
            fout = open("weblogo.txt","wb")
            fout.write(txt_formatter(data, format).decode("utf-8"))    
        except:
            fout = open("weblogo.txt","wb")
            fout.write(txt_formatter(data, format))    
        fout.close()

    else:
        data.alphabet = std_alphabets["protein"]
        options = LogoOptions()
        options.logo_font  = "Arial_Bold"
        options.text_foint = "Helvetica"   
        options.stacks_per_line = 150
        options.number_fontsize    = 34
        options.fontsize           = 34
        options.logo_margin        = 30
        options.stacks_per_line    = 150
        options.scale_width        = 0
        options.stack_width        = 15
        options.stack_aspect_ratio = 5
        options.stroke_width       = 1.5
        options.tic_length         = 7
        options.fineprint          = ""
        options.color_scheme       = std_color_schemes["chemistry"]
        options.yaxis_label        = "Bits"
        options.yaxis_tic_interval = 2
        options.yaxis_minor_tic_ratio = 2   
        format = LogoFormat(data, options)
        try:
            fout = open("weblogo.pdf","wb")
            fout.write(pdf_formatter(data, format))
            fout.close() 
        except:
            fout.close()
            warnings.warn("Ghostscript was not found. If you want to generate LOGO in pdf format, please install Ghostscript.\n For the present, LOGO was output in EPS format") 
            os.system("rm -rf weblogo.pdf") 
            fout = open("weblogo.eps","wb")
            fout.write(eps_formatter(data, format))
            fout.close()
        try:
            fout = open("weblogo.txt","wb")
            fout.write(txt_formatter(data, format).decode("utf-8"))    
        except:
            fout = open("weblogo.txt","wb")
            fout.write(txt_formatter(data, format))    
        fout.close()

def load_data(dtype, strand, ksize ,thresh, Format="pdf"):
    dir_list  = os.listdir("./")
    peak_list = []
    unit_length_list        = [] 
    peak_matrix_list        = [] 
    rpt_unit_array_list     = []
    gap_rpt_unit_array_list = []
    aln_rpt_unit_array_list = []
    avg_intensity_list      = []
    
    for file_name in dir_list:
        kmer_count_signal = []
        kmer_file = open("kmer.tsv")
        kmer_file.readline() 
        for line in kmer_file:
            line = line.rstrip().split("\t") 
            kmer_count_signal.append(int(line[1])) 
        kmer_coount_signal = np.array(kmer_count_signal) 
        gbk = SeqIO.read("repeat.gbk","genbank")

    for feature in gbk.features: 
        if "note" in feature.qualifiers.keys():
            if "SPADE" in feature.qualifiers["note"][0]:
                peak_list.append(int(feature.qualifiers["period"][0]))
                candidate   = list(map(int,feature.qualifiers["hrr_range"][0].split(",")))
                if dtype == "prot":
                    p_start = feature.location.start
                    p_end   = feature.location.end

    for peak in peak_list: 
        for file_name in dir_list:
            if "align.unit_seq.fasta" == file_name:
                se_sets     = []
                gap_se_sets = []
                aln_se_sets = []
                for line in open(file_name).readlines():
                    if line[0] == ">":
                        line = line.rstrip().split("_")
                        s = int(line[1]) 
                        s = s if s > 0 else 0
                    else:
                        unit_length = len(line.rstrip()) 
                        unit_length = peak if unit_length > peak else unit_length
                        seq = line.rstrip() 
                        e = s + unit_length 
                        for m, char in enumerate(seq):
                            if char == "-": 
                                pass
                            else: 
                                break
                        
                        for n, char in enumerate(reversed(seq)):
                            if char == "-": 
                                pass
                            else:
                                break
                        
                        se_sets.append([s,e])
                        aln_se_sets.append([s+m,e-n]) 
                        if m > 0 and len(gap_se_sets) == 0:
                            gap_se_sets.append([s,s+m])
                if n > 0:
                    gap_se_sets.append([e-n,e])
                repeat_unit_array_x = []
                repeat_unit_array_y = []
                for start, end in se_sets:
                    unit_x = [start]  
                    unit_y = [0] 
                    for j in range(start, end+1): 
                        unit_x.append(j)  
                        unit_y.append(1) 
                    unit_x.append(j) 
                    unit_y.append(1)
                    unit_x.append(j) 
                    unit_y.append(0)
                    repeat_unit_array_x.append(unit_x) 
                    repeat_unit_array_y.append(unit_y) 
                rpt_unit_array = zip(repeat_unit_array_x, repeat_unit_array_y) 
                rpt_unit_array_list.append(rpt_unit_array) 
                 
                repeat_unit_array_x = []
                repeat_unit_array_y = []
                for start, end in gap_se_sets:
                    unit_x = [start]  
                    unit_y = [0] 
                    for j in range(start, end+1): 
                        unit_x.append(j)  
                        unit_y.append(1) 
                    unit_x.append(j) 
                    unit_y.append(1)
                    unit_x.append(j) 
                    unit_y.append(0)
                    repeat_unit_array_x.append(unit_x) 
                    repeat_unit_array_y.append(unit_y) 
                rpt_unit_array = zip(repeat_unit_array_x, repeat_unit_array_y) 
                gap_rpt_unit_array_list.append(rpt_unit_array) 
                
                repeat_unit_array_x = []
                repeat_unit_array_y = []
                for start, end in aln_se_sets:
                    unit_x = [start]  
                    unit_y = [0] 
                    for j in range(start, end+1): 
                        unit_x.append(j)  
                        unit_y.append(1) 
                    unit_x.append(j) 
                    unit_y.append(1)
                    unit_x.append(j) 
                    unit_y.append(0)
                    repeat_unit_array_x.append(unit_x) 
                    repeat_unit_array_y.append(unit_y) 
                rpt_unit_array = zip(repeat_unit_array_x, repeat_unit_array_y) 
                aln_rpt_unit_array_list.append(rpt_unit_array) 
                unit_length_list.append(unit_length) 
                
                avg_intensity = []                  
                pdist_file = open("pdist.tsv") 
                pdist_file.readline() 
                for line in pdist_file: 
                    avg_intensity.append(float(line.rstrip().split("\t")[1])) 
                
                avg_intensity = np.array(avg_intensity) 
                if peak > 25:
                    avg_intensity_list.append(avg_intensity[int(len(avg_intensity)/2)-5:int(len(avg_intensity)/2)+5+1]) 
                else:
                    avg_intensity_list.append(avg_intensity) 
                
                matrix_file = open("ppm4vis.tsv")
                matrix_file.readline() 
                matrix_file.readline() 
                peak_matrix = [] 
                periods     = [] 
                for line in matrix_file: 
                    line = list(map(int, line.rstrip().split("\t")))
                    peak_matrix.append(line[1:]) 
                    periods.append(line[0]) 
                periods = list(range(0,periods[-1]+1))  
                peak_matrix = np.array(peak_matrix) 
                if peak > 25:
                    peak_matrix = peak_matrix[int(peak_matrix.shape[0]/2)-5:int(peak_matrix.shape[0]/2)+5+1,:] 
                else:
                    pass 

                if strand == -1:
                    peak_matrix = peak_matrix[:,::-1]
                peak_matrix_list.append(peak_matrix)
        else:
            pass 

    if dtype == "nucl":
        make_figure(peak_matrix_list, kmer_count_signal, rpt_unit_array_list, gap_rpt_unit_array_list, aln_rpt_unit_array_list, unit_length_list, avg_intensity_list, periods, candidate, peak_period_list=peak_list, gb=gbk, fig_name="periodic_repeat", thresh=thresh, dtype=dtype ,Format=Format, k=ksize, strand=strand)
    elif dtype == "prot":
        make_figure(peak_matrix_list, kmer_count_signal, rpt_unit_array_list, gap_rpt_unit_array_list, aln_rpt_unit_array_list, unit_length_list, avg_intensity_list, periods, candidate, peak_period_list=peak_list, gb=gbk, fig_name="periodic_repeat", thresh=thresh, dtype=dtype ,Format=Format, p_start=p_start, p_end=p_end, k=ksize, strand=strand)
    return peak_list

if __name__ == "__main__":
    if sys.argv[1] == "nucl":
        peak_list = load_data("nucl",1,10,0) 
        motif_logo("nucl")
    
    elif sys.argv[1] == "prot":
        peak_list = load_data("prot",1,3,0) 
        motif_logo("prot")
