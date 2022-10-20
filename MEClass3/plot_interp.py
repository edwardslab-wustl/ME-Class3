
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
        x_data = read_param_file_list(param_file)
        if args.number > 0:
            data = pd.read_csv(args.interp_file, nrows=args.number)
        else:
            data = pd.read_csv(args.interp_file)
        feat_cols = [ x for x in data if x.startswith(args.data_type) ]
        for idx in data.index:
            y_data = list(data.loc[idx, feat_cols])
            gene_id = data.iloc[idx,0]
            tag = '.'.join([args.tag, args.anno_type,args.data_type])
            out_file = args.output_path + "/" + gene_id + '.' + tag + '.png'
            plot_interp(x_data, y_data, gene_id,args.data_type, args.anno_type, out_file, args)

def plot_interp(x, y, gene, data_type, anno_type, out_file, args):
    fig, ax = plt.subplots( nrows=1, ncols=1 )
    ax.plot([0,0],[-1,1], linestyle='dashed', color='black')
    ax.plot([-5000,5000],[0,0], color='black')
    ax.plot(x,y)
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
    #parser._action_groups.reverse()
    return(parser)