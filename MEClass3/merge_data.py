
import sys
import argparse

import numpy as np
import pandas as pd

from MEClass3.sample import read_sample_file
from MEClass3.io_functions import print_to_log

def add_tss_interp(df_interp, file):
    df_merged = ''
    tss_interp = pd.read_csv(file)
    if len(df_interp) == 0:
        df_merged = tss_interp
    else:
        df_merged = pd.concat([df_interp, tss_interp], axis=1)
    del tss_interp
    df_merged = df_merged.fillna(0)
    return df_merged

def fix_column_headers(df, tag):
    new_headers = []
    for num, header in enumerate(df.columns):
        if num == 0 or num == len(df.columns) - 1 or num == len(df.columns) - 2:
            new_headers.append(header)
        else:
            new_headers.append(header + "_" + tag)
    df.columns = new_headers
    return df
        
def add_enh_interp(df_interp, file):
    df_merged = ''
    enh_interp_all = pd.read_csv(file)
    enh_interp_all['enh_loc'] = enh_interp_all['enh_loc-gene_id-sample_name'].apply(lambda x: x.split('-',1)[0])
    enh_interp_all['gene_id-sample_name'] = enh_interp_all['enh_loc-gene_id-sample_name'].apply(lambda x: x.split('-',1)[1])
    if len(df_interp) != 0:
        df_merged = df_interp
    for enh_loc in enh_interp_all['enh_loc'].unique():
        df_tmp_interp = enh_interp_all[enh_interp_all['enh_loc'] == enh_loc]
        df_tmp_interp = fix_column_headers(df_tmp_interp, enh_loc)
        if len(df_merged) == 0:
            df_merged = df_tmp_interp
        else:
            tmp = df_tmp_interp.drop('enh_loc-gene_id-sample_name', axis=1)
            tmp.drop('enh_loc', axis=1, inplace=True)
            df_merged = df_merged.merge(tmp, on='gene_id-sample_name')
            del tmp
        del df_tmp_interp
    return df_merged

def exec_merge_data(args):
    expr_file = args.expr
    expr_floor_value = args.expr_floor_val
    pair_list = read_sample_file(args.input_list)
    df_expr_all = ( pd.read_table(expr_file, index_col=False) ).set_index('gene_id')
    with open(args.logfile, 'w') as log_FH:
        for sample_pair in pair_list:
            print_to_log(log_FH, "processing " + sample_pair.name + "\n")
            df_expr = df_expr_all.loc[:, [sample_pair.tag1, sample_pair.tag2]]
            df_interp = ''
            if args.tss:
                interp_file = args.file_path + "/" + sample_pair.name + args.tss_base + ".csv"
                print_to_log(log_FH, "\tadding tss interp data: " + interp_file + "\n")
                df_interp = add_tss_interp(df_interp, interp_file)
            if args.enh:
                interp_file = args.file_path + "/" + sample_pair.name + args.enh_base + ".csv"
                print_to_log(log_FH, "\tadding enh interp data: " + interp_file + "\n")
                df_interp = add_enh_interp(df_interp, interp_file)
                
            df_interp['gene_id'] = df_interp['gene_id-sample_name'].apply(lambda x: x.split('-')[0])
            df_interp['sample_name'] = df_interp['gene_id-sample_name'].apply(lambda x: x.split('-')[1])    
            df_interp = df_interp.set_index('gene_id-sample_name')
            # Remove Duplicates
            # df_interp = df_interp.groupby(df_interp.index).first()
            # Remove Duplicates
            # df_expr = df_expr.groupby(df_expr.index).first()
            if args.floor_expr:    # Floor expression values for selected cell types
                df_expr.loc[df_expr[sample_pair.tag1] < expr_floor_value, sample_pair.tag1] = expr_floor_value
                df_expr.loc[df_expr[sample_pair.tag2] < expr_floor_value, sample_pair.tag2] = expr_floor_value
                
            if args.dexpr_flag == 0:   
                df_expr[sample_pair.name] = np.where( df_expr[sample_pair.tag1] > df_expr[sample_pair.tag2],
                                                     -(df_expr[sample_pair.tag1] / df_expr[sample_pair.tag2]),
                                                      (df_expr[sample_pair.tag2] / df_expr[sample_pair.tag1]) )
                df_expr[sample_pair.name] = np.where( df_expr[sample_pair.tag1 == 1.0, 0.0, df_expr[sample_pair.name]])
            elif args.dexpr_flag == 1:   
                df_expr[sample_pair.name] = np.where( df_expr[sample_pair.tag1] > df_expr[sample_pair.tag2],
                                                     -(df_expr[sample_pair.tag1] / df_expr[sample_pair.tag2]),
                                                      (df_expr[sample_pair.tag2] / df_expr[sample_pair.tag1]) )
            elif args.dexpr_flag == 2:   
                df_expr[sample_pair.name] = np.log2( df_expr[sample_pair.tag2] / df_expr[sample_pair.tag1] )
            elif args.dexpr_flag == 3:   
                df_expr[sample_pair.name] = np.where( df_expr[sample_pair.tag1] > df_expr[sample_pair.tag2],
                                                     -(df_expr[sample_pair.tag1] / df_expr[sample_pair.tag2]),
                                                      (df_expr[sample_pair.tag2] / df_expr[sample_pair.tag1]) )
                
            df_interp['expr_value'] = 0.0  # Float type
            df_interp['expr_flag'] = 0      # int type
            for gene_item in df_interp.index:
            #    print(gene_item)
                try:
                    df_interp.loc[gene_item, 'expr_value'] = df_expr.loc[ gene_item.split('-')[0], gene_item.split('-')[1] ]
                except KeyError:
                    print_to_log(log_FH, gene_item.split('-')[0]+' not found in expression file.\n')
                    df_interp.loc[gene_item, 'expr_value'] = 0
            if args.dexpr_flag == 0: # me-class demo
                df_interp['expr_flag'] = np.where( df_interp['expr_value'] > 0, 1, df_interp['expr_flag'] )
                df_interp['expr_flag'] = np.where( df_interp['expr_value'] < 0, -1, df_interp['expr_flag'] )
                df_interp['expr_flag'] =  df_interp['expr_flag'].astype(int)
            elif args.dexpr_flag == 1: # me-class paper
                df_interp['expr_flag'] = np.where( df_interp['expr_value'] >= 2, 1, df_interp['expr_flag'] )
                df_interp['expr_flag'] = np.where( df_interp['expr_value'] <= -2, -1, df_interp['expr_flag'] )
                df_interp['expr_flag'] =  df_interp['expr_flag'].astype(int)
            elif args.dexpr_flag == 2:
                df_interp['expr_flag'] = np.where( df_interp['expr_value'] > 0.0, 1, df_interp['expr_flag'] )
                df_interp['expr_flag'] = np.where( df_interp['expr_value'] < 0.0, -1, df_interp['expr_flag'] )
                df_interp['expr_flag'] = df_interp['expr_flag'].astype(int)
            elif args.dexpr_flag == 3: # custom
                df_interp['expr_flag'] = np.where( df_interp['expr_value'] >= args.expr_cutoff, 1, df_interp['expr_flag'] )
                df_interp['expr_flag'] = np.where( df_interp['expr_value'] <= -args.expr_cutoff, -1, df_interp['expr_flag'] )
                df_interp['expr_flag'] =  df_interp['expr_flag'].astype(int)
            # Finally write data for classifier
            out_file = args.output_path + '/' + sample_pair.name + '.interp_expr_data.csv' 
            print_to_log(log_FH, "\tprinting to outfile: " + out_file)
            df_interp.to_csv(out_file, sep=',')

def exec_merge_data_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-i', '--input_list', action='store',
        dest='input_list', required=True, 
        default=argparse.SUPPRESS,
        help='Input list of sample names and file locations for pairings.')
    parser_required.add_argument('-e', '--expr',
        default=argparse.SUPPRESS,
        required=True, help='Name of expression file')
    parser.add_argument('-p', '--file_path',
        default='intermediate_files', help='Path to directory with interpolation files')
    parser.add_argument('-o', '--output_path',
        default='intermediate_files', help='Path to directory to store output files')
    parser.add_argument('--tss', action='store_true', default=False, help='Use tss interpolation data')
    parser.add_argument('--tss_base', type=str, default='_gene_interp', help='Base part of tss interpolation file name')
    parser.add_argument('--enh', action='store_true', default=False, help='Use enh interpolation data')
    parser.add_argument('--enh_base', type=str, default='_enh_interp', help='Base part of enh interpolation file name')
    #parser_required.add_argument('-expr', action='store', dest='expr_inp', required=True, help='Expression file')
    #parser_required.add_argument('-intrp', action='store', dest='intrp_inp', required=True, help='Interpolation CSV file')
    parser.add_argument('--floor_expr', type=bool, default=True, help='Floor expression value?')
    parser.add_argument('--expr_floor_val', dest='expr_floor_val', type=float, default=5.0, help='Expression floor value')
    parser.add_argument('--diff_expr_flag', dest='dexpr_flag', type=int, default=1, help='Method for differential expression')
    parser.add_argument('--expr_cutoff', dest='expr_cutoff', type=int, default=2, help='fold change in expression cutoff, must use -def 3')
    parser.add_argument('--logfile', action='store', dest='logfile',
        default='merge_data.log', help='log file')
    parser._action_groups.reverse()
    return(parser)
