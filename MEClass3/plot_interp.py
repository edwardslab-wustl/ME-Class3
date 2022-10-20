
import argparse
import os.path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import format_args_to_print
from MEClass3.io_functions import mk_output_dir
from MEClass3.io_functions import eprint
from MEClass3.io_functions import read_param_file_list
from MEClass3.io_functions import read_gene_file
from MEClass3.io_functions import index_raw_data

def exec_plot_interp(args):
    with open(args.logfile, 'w') as log_FH:
        print_to_log(log_FH, format_args_to_print(args))
        mk_output_dir(args.output_path)
        if args.interp_param_file:
            param_file = args.interp_param_file
        else:
            param_file = args.interp_file.strip(".csv") + ".param"
        if not os.path.exists(param_file):
            eprint(f"Can't find suitable param_file: {param_file}\n")
        x_data,param = read_param_file_list(param_file)
        if args.number > 0:
            data = pd.read_csv(args.interp_file, nrows=args.number)
        else:
            data = pd.read_csv(args.interp_file)
        feat_cols = [ x for x in data if x.startswith(args.data_type) ]
        eprint("reading in raw data")
        if args.raw_data:
            raw_data_dict, raw_data_index, raw_data_index_size = index_raw_data(args.raw_data, param.region_size)
            region_dict = grab_region_data(args.anno_file, args.anno_type)
        eprint("looping genes")
        for idx in data.index:
            raw_x_data = []
            raw_y_data = []
            y_data = list(data.loc[idx, feat_cols])
            id = data.iloc[idx,0]
            (gene_id, sample_id) = id.split('-')
            eprint(f"  on gene {gene_id}")
            file_tag = '.'.join([gene_id, sample_id, args.tag, args.anno_type,args.data_type])
            out_file = args.output_path + "/" + file_tag + '.png'
            if args.raw_data:
                if gene_id in region_dict:
                    (region_chr, region_ref_pt, region_strand) = region_dict[gene_id]
                    if region_chr in raw_data_dict:
                        raw_x_data, raw_y_data = extract_raw_data(gene_id, region_ref_pt, region_chr, region_strand,
                                                                  raw_data_dict[region_chr], raw_data_index, raw_data_index_size, param, args)
                    else:
                        print_to_log(log_FH, f"Skipping raw data for {gene_id}. Can't find chrom in raw_data, chrom: {region_chr}\n")
                else:
                    print_to_log(log_FH, f"Skipping raw data for {gene_id}. Can't find gene in anno_data.\n")
                if len(raw_x_data) == 0:
                    print_to_log(log_FH, f"Skipping raw data for {gene_id}. Can't find any raw_data.\n")
            plot_interp(x_data, y_data, raw_x_data, raw_y_data, gene_id,args.data_type, args.anno_type, out_file, args)

def extract_raw_data( id, ref_pt, chrom, strand, raw_data_dict_chr, raw_data_index, raw_data_index_size, param, args):
    raw_x_data = []
    raw_y_data = []
    start_pos = ref_pt - param.region_size - 100 #pull a bit extra to avoid boundary issues
    end_pos = ref_pt + param.region_size + 100 #pull a bit extra to avoid boundary issues
    if start_pos < 0:
        start_pos = 0
    start_idx = int(start_pos / raw_data_index_size) - 1
    end_idx = int(end_pos / raw_data_index_size) + 1
    for idx in range(start_idx, end_idx + 1):
        for pos in raw_data_index[(chrom, idx)]:
            if start_pos <= pos and pos <= end_pos:
                if strand == '+':
                    raw_x_data.append(pos - ref_pt)
                else:
                    raw_x_data.append(-pos + ref_pt)
                raw_y_data.append(raw_data_dict_chr[pos])
    return raw_x_data, raw_y_data

def grab_region_data(anno_file, type):
    region_dict = dict()
    if anno_file and os.path.exists(anno_file):
        if type == 'tss':
            region_list = read_gene_file(anno_file)
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
        region_dict[region.id] = (region.chr, region.tss(), region.strand)
    return region_dict

def plot_interp(x, y, raw_x_data, raw_y_data, gene, data_type, anno_type, out_file, args):
    fig, ax = plt.subplots( nrows=1, ncols=1 )
    ax.plot([0,0],[-1,1], linestyle='dashed', color='black')
    ax.plot([-5000,5000],[0,0], color='black')
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
    else:
        ax.set_xlabel("Position (bp)")
    plt.ylim(-args.max_y_val,args.max_y_val)
    plt.xlim([-5000,5000])
    plt.xticks((-5000,0,5000),('-5kb','TSS','+5kb'))
    x_tick_minorLocator = MultipleLocator(1000)
    ax.xaxis.set_minor_locator(x_tick_minorLocator)
    plt.savefig(out_file, bbox_inches='tight')
    plt.close()
    return

def exec_plot_interp_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-i','--interp_file', default=argparse.SUPPRESS, required=True, help='interpolation file')
    parser_general = parser.add_argument_group('general arguments')
    parser_general.add_argument('-t', '--tag', default='interp', help='Tag for output')
    parser_general.add_argument('--interp_param_file', default=None, 
        help='for interp file <base>.csv assumed to be <base>.param.csv')
    parser_general.add_argument('--raw_data', default=None, 
        help='Raw bedgraph data. Will only plot if supplied.')
    parser_general.add_argument('--anno_file', default=None, 
        help='Must supply this if you want to plot raw data. Ignored otherwise.')
    parser_general.add_argument('--anno_type', default="tss", 
        choices=["tss", "enh"],
        help='region or gene annotation file')
    parser_general.add_argument('--data_type', default="mC", 
        choices=["mC", "hmC", "other"],
        help='type of data')
    parser_general.add_argument('-n', '--number', default=0, type=int,
        help='Number of genes to print. Starts at top of file. Set to 0 to print all.')
    #parser_general.add_argument('--num_proc', type=int, default=0, help='number of processes to run')
    parser_general.add_argument('-o', '--output_path', action='store', dest='output_path',
        default='interp_plots', help='Path to Output')
    parser_general.add_argument('--logfile', dest='logfile', default='plot_interp.log', help='log file')
    parser_plot = parser.add_argument_group('plotting arguments')
    parser_plot.add_argument('--max_y_val', default=0.8, type=float,
        help='Set y-axis limits as [-max_y_val, max_y_val].')
    parser_plot.add_argument('--raw_data_color', default='red', type=str,
        help='Marker color for raw data.')
    parser_plot.add_argument('--interp_data_color', default='blue', type=str,
        help='Interpolation line color.')
    #parser._action_groups.reverse()
    return(parser)