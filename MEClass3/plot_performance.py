
#--------------------------------------------
# -*- coding: utf-8 -*-
"""
Created on June 8th 2017

@author: cschlosberg & jedwards
"""
#-------------------------------------------

import copy

import matplotlib as mpl
import pandas as pd

from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import format_args_to_print
from MEClass3.plot_performance_functions import Sample
from MEClass3.plot_performance_functions import calculate_labels
from MEClass3.plot_performance_functions import calculate_results
from MEClass3.plot_performance_functions import calculate_roc
from MEClass3.plot_performance_functions import print_roc
from MEClass3.plot_performance_functions import plot_roc
from MEClass3.plot_performance_functions import calculate_acc_rejectrate
from MEClass3.plot_performance_functions import print_acc_rejectrate
from MEClass3.plot_performance_functions import plot_acc_rejectrate
from MEClass3.plot_performance_functions import set_matplotlib_params


def exec_plot_performance(args):
    with open(args.logfile, 'w') as log_FH:
        print_to_log(log_FH, format_args_to_print(args))
        samples = list() 
        for csv_file in args.input_files:
            df = pd.read_csv(csv_file).set_index('gene_id-sample_name')
            if args.split_samples:
                for s in df['sample_name'].unique():
                    sample_name = csv_file.split(".")[0] + ',' + s
                    legend = csv_file.split('.')[0]+'.'+csv_file.split('.')[1] + ',' + s
                    df_tmp = df[df['sample_name'] == s]
                    expression = df_tmp['expr_value'].tolist() 
                    gene_ID = df_tmp['gene_id'].tolist()
                    predictions = df_tmp['prob_up'].tolist()
                    sample = Sample(sample_name, legend, predictions, expression, gene_ID)
                    samples.append(sample)
            else:
                sample_name = csv_file.split('.')[0]
                legend = csv_file.split('.')[0]+'.'+csv_file.split('.')[1]
                expression = df['expr_value'].tolist() 
                gene_ID = df['gene_id'].tolist()
                predictions = df['prob_up'].tolist()
                sample = Sample(sample_name, legend, predictions, expression, gene_ID)
                samples.append(sample)
                del df

        calculate_labels(samples, args)
        results = calculate_results(samples,args)
        results_P = copy.deepcopy(results)
        results_N = copy.deepcopy(results)
        calculate_roc(results,args)
        print_roc(results, args)
        calculate_acc_rejectrate(results,'all',args)
        calculate_acc_rejectrate(results_P,'pos',args)
        calculate_acc_rejectrate(results_N,'neg',args)
        print_acc_rejectrate(results, 'all', args)
        print_acc_rejectrate(results_P, 'pos', args)
        print_acc_rejectrate(results_N, 'neg', args)
        if args.plot_results:
            version = int(mpl.__version__.split('.')[0])
            set_matplotlib_params(args, version)
            outFile = args.outFileBase + ".roc.png"
            plot_roc(results, outFile, args, version)
            outFile = args.outFileBase + ".all.acc_rejectrate.png"
            plot_acc_rejectrate(results, outFile, args, version)
            outFile = args.outFileBase + ".pos.acc_rejectrate.png"
            plot_acc_rejectrate(results_P, outFile, args, version)
            outFile = args.outFileBase + ".neg.acc_rejectrate.png"
            plot_acc_rejectrate(results_N, outFile, args, version)
    return

def exec_plot_performance_help(parser):
    ### Required positional arguments
    parser.add_argument('input_files',nargs='+',help="prediction files in csv format") 
    parser.add_argument('-o', '--outFileBase',
        default='results',
        help="base tag for output files, default = \"results\"")
#    parser.add_argument('--label_file',dest='label_file',default='',
#            help="label file if there is a single common label file for all\
#            .pred files. Otherwise (default behaviour) assume for evey\
#            <base>.pred file there is a corresponding <base>.label file.")
    parser.add_argument('--plot_results',
            action='store_true',
            help="Plot results.")
    parser.set_defaults(plot_results=True)
    parser.add_argument('--split_samples', action='store_true',
            default=False,
            help="Split results by sample.")
    parser.add_argument('--logfile', dest='logfile', default='plot_performance.log', help='log file')
    parser.add_argument('--verbose',dest='verbose',
            action='store_true',
            help="Print extra information to stderr.")
    parser.set_defaults(verbose=False)
    ### Advanced Plotting Options
    advanced_plot_group = parser.add_argument_group(
            title='Advanced Plotting Options')
    advanced_plot_group.add_argument('--acc_rejectrate_steps',
            help="Number of points to use in accuracy vs reject rate plot, 11 steps \
                    would be points every 0.1, default=101",
            type=int,
            default=101)
    advanced_plot_group.add_argument('--legendLoc',
            help="location of legend for plots, default=best",
            choices=["best","right"],
            default="best")
    advanced_plot_group.add_argument('--lineWidth',
            help="Line width for plots, default=t.5",
            type=float,
            default=2.5)
    return parser

#if __name__ == "__main__":
#    main()    
    
    
    
