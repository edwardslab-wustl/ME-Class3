
import argparse

import pandas as pd

from MEClass3.io_functions import format_args_to_print
from MEClass3.io_functions import print_to_log

def pull_gene_list(df, accuracy):
    up_df = df[df['prob_up'] >= accuracy]
    dn_df = df[df['prob_dn'] >= accuracy]
    result = pd.concat([up_df,dn_df])
    return result

def exec_pull_gene_list(args):
    with open(args.logfile, 'w') as log_FH:
        print_to_log(log_FH, format_args_to_print(args))
        df = pd.read_csv( args.input_file, index_col=[0] ).reset_index()
        tag = (args.input_file).rstrip('.csv')
        if args.out_file == None:
            gene_file = f"{tag}.acc_{args.accuracy}.geneList.csv"
        else:
            gene_file = args.out_file
        gene_list  = pull_gene_list(df, args.accuracy)
        if not args.print_all_preds:
            gene_list = gene_list[gene_list['expr_flag'] * gene_list['expr_pred'] == 1]
        gene_list.to_csv(gene_file, index=False)
    return

def exec_pull_gene_list_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-i', '--input_file', 
        required=True,
        default=argparse.SUPPRESS,
        help='Dataframe output from classification')
    parser.add_argument('-o', '--out_file', default=None,
        help='output file. If None, then if input is <base>.csv, will store results at <base>.acc_<accuracy>.geneList.csv')
    parser.add_argument('--accuracy', type=float, default=0.85,
        help='Accuracy cutoff for gene list')
    parser.add_argument('--print_all_preds', action='store_true', default=False,
        help='Print all predictions, correct and incorrect')
    parser.add_argument('--logfile', default='pull_gene_list.log', help='log file')
    parser._action_groups.reverse()
    return(parser)
