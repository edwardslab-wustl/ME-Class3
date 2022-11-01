
import argparse
import multiprocessing as mp
from itertools import repeat

import pandas as pd

from MEClass3.io_functions import mk_output_dir
from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import eprint
from MEClass3.io_functions import format_args_to_print
from MEClass3.sample import read_bed_methyl_data
from MEClass3.sample import read_sample_file
from MEClass3.sample import read_bedms_methyl_data

def exec_preprocess(args):
    with open(args.logfile, 'w') as log_FH:
        print_to_log(log_FH, format_args_to_print(args))
        pair_list = read_sample_file(args.input_list)
        if args.sample:
            sample_set_flag = False
            for sample_pair in pair_list:
                if args.sample == sample_pair.name:
                    pair_list = [ sample_pair ]
                    sample_set_flag = True
            if not sample_set_flag:
                eprint("Can't find sample name " + args.sample + " in pair list. Check pair list and sample name and rerun.")
                exit()
        output_path = args.output_path
        mk_output_dir(output_path)
        preprocess_results = []
        if args.num_proc > 1:
            arg_iterable = zip(pair_list, repeat(output_path), repeat(args.data_format), repeat(args.sig_digits), repeat(args.divide_by_100))
            with mp.Pool(processes=args.num_proc) as pool:
                preprocess_results = pool.starmap(preprocess_sample, arg_iterable)
        else:
            for sample_pair in pair_list:
                (name, success_flag, error) = preprocess_sample(sample_pair, output_path, args.data_format, args.sig_digits, args.divide_by_100)
                preprocess_results.append((name,success_flag,error))
        for (name, success_flag, error) in preprocess_results:
            if success_flag:
                print_to_log(log_FH, f"Sucessfully processed {name}\n")
            else:
                print_to_log(log_FH, f"Failed to process {name}: {error}")
    return

def preprocess_sample(sample, out_path, data_format, sig_digits, divide_by_100_flag ):
    success_flag = False
    error = ''
    out_file = out_path + '/' + sample.name +'.bedgraph'
    if data_format == 'bed':
        df1_bed = read_bed_methyl_data(sample.file1, 1, divide_by_100_flag)
        df2_bed = read_bed_methyl_data(sample.file2, 2, divide_by_100_flag)
        success_flag = True
    elif data_format == 'bedms':
        df1_bed = read_bedms_methyl_data(sample.file1, 1, divide_by_100_flag)
        df2_bed = read_bedms_methyl_data(sample.file2, 2, divide_by_100_flag)
        success_flag = True
    else:
        error = f"Unrecognized data format: {data_format}\n"
        success_flag = False
    if success_flag:
        df_merged = df_merged = pd.concat( [ df1_bed, df2_bed ], axis=1, join='inner')
        df_merged['dcm'] = (df_merged['value2'] - df_merged['value1']).round(sig_digits)
        cols_to_keep = ['chrom1', 'start1', 'end1', 'dcm']
        df_merged[cols_to_keep].to_csv(out_file, sep='\t', na_rep='NA', header=None, index=False)
        df1_bed.drop(df1_bed.index, inplace=True)
        df2_bed.drop(df2_bed.index, inplace=True)
        del df1_bed, df2_bed, df_merged
    return sample.name, success_flag, error
            
def exec_preprocess_old(args):
    pair_list = read_sample_file(args.input_list)
    if args.sample:
        sample_set_flag = False
        for sample_pair in pair_list:
            if args.sample == sample_pair.name:
                pair_list = [ sample_pair ]
                sample_set_flag = True
        if not sample_set_flag:
            eprint("Can't find sample name " + args.sample + " in pair list. Check pair list and sample name and rerun.")
            exit()
    output_path = args.output_path
    mk_output_dir(output_path)
    with open(args.logfile, 'w') as log_FH:
        print_to_log(log_FH, format_args_to_print(args))
        for sample_pair in pair_list:
            print_to_log(log_FH, "processing " + sample_pair.name + "\n")
            if args.data_format == 'bed':
                df1_bed = read_bed_methyl_data(sample_pair.file1, 1, args.divide_by_100)
                df2_bed = read_bed_methyl_data(sample_pair.file2, 2, args.divide_by_100)
            elif args.data_format == 'bedms':
                df1_bed = read_bedms_methyl_data(sample_pair.file1, 1, args.divide_by_100)
                df2_bed = read_bedms_methyl_data(sample_pair.file2, 2, args.divide_by_100)
            else:
                eprint("Unrecognized data format: " + args.data_format + "\n")
                exit
            df_merged = df_merged = pd.concat( [ df1_bed, df2_bed ], axis=1, join='inner')
            df_merged['dcm'] = (df_merged['value2'] - df_merged['value1']).round(args.sig_digits)
            output_file = sample_pair.name +'.bedgraph'
            cols_to_keep = ['chrom1', 'start1', 'end1', 'dcm']
            df_merged[cols_to_keep].to_csv(output_path+'/'+output_file, sep='\t', na_rep='NA', header=None, index=False)
            df1_bed.drop(df1_bed.index, inplace=True)
            df2_bed.drop(df2_bed.index, inplace=True)
            del df1_bed, df2_bed, df_merged
             
def exec_preprocess_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-i', '--input_list', action='store',
        dest='input_list', required=True,
        default=argparse.SUPPRESS,
        help='Input list of sample names and file locations for pairings.')
    parser.add_argument('--sample', default=None, help='Use to select specific sample from pair list and only run that one.')
    parser.add_argument('-o', '--output_path', action='store',
        default='intermediate_files', help='Path to Output')
    parser.add_argument('--data_format', choices=('bed', 'bedms'), default='bed',
        help='format of data files, see README.md for more info on types')
    parser.add_argument('--sig_digits', action='store', type=int,
        default=3, help='Significant digits for methylation difference')
    parser.add_argument('--divide_by_100', action='store_true', default=False,
        help="divide methylation values by 100 to put on range [0,1]")
    parser.add_argument('--num_proc', type=int, default=1, help='number of processes to run')
    parser.add_argument('--logfile', action='store', dest='logfile',
        default='preprocess.log', help='log file')
    parser._action_groups.reverse()
    return(parser)



###########################
# OLD CODE
###########################
      
# def exec_proc_sample_wheel(args):
#     ctp_fofn = args.ctp_inp
#    path_to_output = args.pto_inp
#    expr_file = args.expr_inp
#        
#    df_expr = pd.read_table(expr_file).set_index('gene_id')
#    expr_cell_ids = list(df_expr.columns.values)
#    
#    with open(ctp_fofn, 'r') as ctp_fofn_file:
#        ctp_fofn_file_lines = ctp_fofn_file.readlines()
#    ctp_fofn_file.close()
#    
#    #ctp_fofn_file_lines_fltrd = list( item for item in ctp_fofn_file_lines if item.strip().split()[0] in expr_cell_ids )
#    ctp_fofn_file_lines_fltrd = ctp_fofn_file_lines
#    
#    # Wheel difference
#    i = 0
#    while i < len(ctp_fofn_file_lines_fltrd):
#        first_info = ctp_fofn_file_lines_fltrd[i-1] 
#        second_info = ctp_fofn_file_lines_fltrd[i]
#
#        bg1 = first_info.strip().split()[1] # First .bedgraph
#        bg2 = second_info.strip().split()[1] # Second .bedgraph
#        tag = first_info.strip().split()[0]+'_'+second_info.strip().split()[0]
#
#        df1_bed = pd.read_table(bg1, index_col=False, \
#                na_values = 'NA', names = ['chrom1', 'start1', 'end1', 'value1'])
##        df1_bed.sort_values(by=['chrom1', 'start1']) # sort (doesn't change answer)    
#        df1_bed['idx_value1'] = df1_bed['chrom1']+'_'+df1_bed['start1'].astype(str)
#        df1_bed = df1_bed.set_index('idx_value1')
#
#        df2_bed = pd.read_table(bg2, index_col=False, \
#                na_values = 'NA', names = ['chrom2', 'start2', 'end2', 'value2'])
##        df2_bed.sort_values(by=['chrom2', 'start2']) # sort (doesn't change answer)
#        df2_bed['idx_value2'] = df2_bed['chrom2']+'_'+df2_bed['start2'].astype(str)
#        df2_bed = df2_bed.set_index('idx_value2')
#
#        df_merged = df_merged = pd.concat( [ df1_bed, df2_bed ], axis=1, join='inner')
#        df_merged['dcm'] = (df_merged['value2'] - df_merged['value1']).round(2)
##        df_merged['start1end'] = df_merged['start1'] + 1
#
#        output_file = tag+'.bedgraph'
#        cols_to_keep = ['chrom1', 'start1', 'end1', 'dcm']
#        df_merged[cols_to_keep].to_csv(path_to_output+'/'+output_file, sep='\t', na_rep='NA', header=None, index=False)
#        
#        i += 1
#        
#        df1_bed.drop(df1_bed.index, inplace=True)
#        df2_bed.drop(df2_bed.index, inplace=True)
#        df1_bed.drop(df1_bed.index, inplace=True)
#        del df1_bed, df2_bed, df_merged
#
#
#def exec_proc_sample_wheel_help(parser):
#    parser_required = parser.add_argument_group('required arguments')
#    parser_required.add_argument('-ctp', action='store', dest='ctp_inp', required=True, help='Cell type Fofn')
#    parser_required.add_argument('-expr', action='store', dest='expr_inp', required=True, help='Name of expression file')
#    parser.add_argument('-pto', action='store', dest='pto_inp', default='.', help='Path to Output')
#    parser._action_groups.reverse()
#    return(parser)