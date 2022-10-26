
import os.path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from MEClass3.io_functions import eprint
from MEClass3.io_functions import read_gene_file
from MEClass3.io_functions import read_region_file

def extract_raw_data(start_ref_pt, end_ref_pt, chrom, strand, raw_data_dict_chr, raw_data_index, raw_data_index_size, param, args):
    raw_x_data = []
    raw_y_data = []
    start_pos = start_ref_pt - param.region_size - 100 #pull a bit extra to avoid boundary issues
    end_pos = end_ref_pt + param.region_size + 100 #pull a bit extra to avoid boundary issues
    if start_pos < 0:
        start_pos = 0
    start_idx = int(start_pos / raw_data_index_size) - 1
    end_idx = int(end_pos / raw_data_index_size) + 1
    for idx in range(start_idx, end_idx + 1):
        for pos in raw_data_index[(chrom, idx)]:
            if start_pos <= pos and pos <= end_pos:
                if pos < start_ref_pt:
                    rel_pos =  pos - start_ref_pt
                elif start_ref_pt <= pos and pos <= end_ref_pt:
                    rel_pos = 0
                else:
                    rel_pos =  pos - end_ref_pt
                if strand == '+':
                    raw_x_data.append(rel_pos)
                else:
                    raw_x_data.append(-rel_pos)
                raw_y_data.append(raw_data_dict_chr[pos])
    return raw_x_data, raw_y_data

def grab_region_data(anno_file, type):
    region_dict = dict()
    if anno_file and os.path.exists(anno_file):
        if type == 'tss':
            region_list = read_gene_file(anno_file)
        elif type == 'enh':
            region_list = read_region_file(anno_file)
        else:
            eprint(f"Can't recognize anno_type: {type}\n")
            exit()
    elif anno_file:
        eprint("Must specify anno_file param for feature annotation to add raw data.\n")
        exit()
    else:
        eprint(f"Can't find anno_file: {anno_file}\n")
        exit()
    for region in region_list:
        start, end = pull_start_end(region, type)
        region_dict[region.id] = (region.chr, start, end, region.strand)
    return region_dict

def pull_start_end(anno, type):
    if type == 'tss':
        start = anno.tss()
        end = start
    elif type == 'enh':
        start = anno.start
        end = anno.end
    return start, end

def plot_interp(x, y, raw_x_data, raw_y_data, gene, data_type, anno_type, out_file, param, args):
    fig, ax = plt.subplots( nrows=1, ncols=1 )
    x_range = [-param.region_size,param.region_size]
    ax.plot([0,0],[-1,1], linestyle='dashed', color='black')
    ax.plot(x_range,[0,0], color='black')
    ax.plot(x,y, color=args.interp_data_color)
    if len(raw_x_data) > 0:
        ax.plot(raw_x_data, raw_y_data,
                marker="o",
                markersize=2,
                linestyle='None',
                color=args.raw_data_color)
    plt.title(gene, loc='left')
    ax.set_ylabel(r'$\Delta$' + data_type)
    if anno_type == 'tss':
        ax.set_xlabel("Distance to TSS (bp)")
    elif anno_type == 'enh':
        ax.set_xlabel("Distance to Enhancer (bp)")
    else:
        ax.set_xlabel("Relative Position (bp)")
    plt.ylim(-args.max_y_val,args.max_y_val)
    plt.xlim(x_range)
    if anno_type == 'tss':
        feat_label = "TSS"
        x_tick_minorLocator = MultipleLocator(1000)
        ax.xaxis.set_minor_locator(x_tick_minorLocator)
    elif anno_type == 'enh':
        feat_label = "Enhancer"
        x_tick_minorLocator = MultipleLocator(100)
        ax.xaxis.set_minor_locator(x_tick_minorLocator)
    else:
        feat_label = anno_type
        x_tick_minorLocator = MultipleLocator(100)
        ax.xaxis.set_minor_locator(x_tick_minorLocator)
    plt.xticks((-param.region_size,0,param.region_size),(str(-param.region_size),feat_label,str(param.region_size)))
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()
    return
