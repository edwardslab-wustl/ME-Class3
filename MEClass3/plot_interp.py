
import argparse
import os.path

import pandas as pd

from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import format_args_to_print
from MEClass3.io_functions import mk_output_dir
from MEClass3.io_functions import eprint
from MEClass3.io_functions import read_param_file_list
from MEClass3.io_functions import index_raw_data
from MEClass3.plot_interp_functions import grab_region_data
from MEClass3.plot_interp_functions import extract_raw_data
from MEClass3.plot_interp_functions import plot_interp

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
        #eprint("reading in raw data")
        if args.raw_data:
            raw_data_dict, raw_data_index, raw_data_index_size = index_raw_data(args.raw_data, param.region_size)
            region_dict = grab_region_data(args.anno_file, args.anno_type)
        #eprint("looping interp lines")
        for idx in data.index:
            raw_x_data = []
            raw_y_data = []
            y_data = list(data.loc[idx, feat_cols])
            id = data.iloc[idx,0]
            if args.anno_type == 'tss':
                (gene_id, sample_id) = id.split('-')
                #eprint(f"  on gene {gene_id}")
                file_tag = '.'.join([gene_id, sample_id, args.tag, args.anno_type, args.data_type])
                anno_info = ",".join([sample_id,gene_id])
                feat_id = gene_id
            elif args.anno_type == 'enh':
                (enh_info, gene_id, sample_id) = id.split('-')
                #eprint(f"  processing {enh_info} {gene_id}")
                file_tag = '.'.join([gene_id, sample_id, enh_info, args.tag, args.anno_type, args.data_type])
                anno_info = ",".join([sample_id,gene_id,enh_info])
                feat_id = enh_info + '-' + gene_id
            if args.raw_data:
                if feat_id in region_dict:
                    (region_chr, region_start_ref_pt, region_end_ref_pt, region_strand) = region_dict[feat_id]
                    if region_chr in raw_data_dict:
                        raw_x_data, raw_y_data = extract_raw_data(region_start_ref_pt, region_end_ref_pt, region_chr, region_strand,
                                                                raw_data_dict[region_chr], raw_data_index, raw_data_index_size, param, args)
                    else:
                        print_to_log(log_FH, f"Skipping raw data for {feat_id}. Can't find chrom in raw_data, chrom: {region_chr}\n")
                else:
                    print_to_log(log_FH, f"Skipping raw data for {feat_id}. Can't find feat in anno_data.\n")
                if len(raw_x_data) == 0:
                    print_to_log(log_FH, f"Skipping raw data for {feat_id} {gene_id}. Can't find any raw_data.\n")
            out_file = args.output_path + "/" + file_tag + '.png'
            plot_interp(x_data, y_data, raw_x_data, raw_y_data, anno_info, args.data_type, args.anno_type, out_file, param, args)

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
    parser_plot.add_argument('--max_y_val', default=1.0, type=float,
        help='Set y-axis limits as [-max_y_val, max_y_val].')
    parser_plot.add_argument('--raw_data_color', default='red', type=str,
        help='Marker color for raw data.')
    parser_plot.add_argument('--interp_data_color', default='blue', type=str,
        help='Interpolation line color.')
    #parser._action_groups.reverse()
    return(parser)