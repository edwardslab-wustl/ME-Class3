import pandas as pd

from MEClass3.io_functions import mk_output_dir
from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import eprint
from MEClass3.sample import read_bed_methyl_data
from MEClass3.sample import read_sample_file
from MEClass3.sample import read_bedms_methyl_data

def exec_proc_sample(args):
    pair_list = read_sample_file(args.input_list)
    output_path = args.output_path
    mk_output_dir(output_path)
    with open(args.logfile, 'w') as log_FH:
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
  
def exec_proc_sample_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-i', '--input_list', action='store',
        dest='input_list', required=True, help='Input list of sample names and file locations for pairings.')
    #parser_required.add_argument('-e', '--expr', action='store',
    #    dest='expr_input', required=True, help='Name of expression file')
    parser.add_argument('-o', '--output_path', action='store', dest='output_path',
        default='intermediate_files', help='Path to Output')
    parser.add_argument('--data_format', choices=('bed', 'bedms'), default='bed',
        help='format of data files, see README.md for more info on types')
    parser.add_argument('--sig_digits', action='store', type=int,
        default=3, help='Significant digits for methylation difference')
    parser.add_argument('--divide_by_100', action='store_true', default=False,
        help="divide methylation values by 100 to put on range [0,1]")
    parser.add_argument('--logfile', action='store', dest='logfile',
        default='proc_sample.log', help='log file')
    parser._action_groups.reverse()
    return(parser)



###########################
# OLD CODE
###########################
      
def exec_proc_sample_wheel(args):
    ctp_fofn = args.ctp_inp
    path_to_output = args.pto_inp
    expr_file = args.expr_inp
        
    df_expr = pd.read_table(expr_file).set_index('gene_id')
    expr_cell_ids = list(df_expr.columns.values)
    
    with open(ctp_fofn, 'r') as ctp_fofn_file:
        ctp_fofn_file_lines = ctp_fofn_file.readlines()
    ctp_fofn_file.close()
    
    #ctp_fofn_file_lines_fltrd = list( item for item in ctp_fofn_file_lines if item.strip().split()[0] in expr_cell_ids )
    ctp_fofn_file_lines_fltrd = ctp_fofn_file_lines
    
    # Wheel difference
    i = 0
    while i < len(ctp_fofn_file_lines_fltrd):
        first_info = ctp_fofn_file_lines_fltrd[i-1] 
        second_info = ctp_fofn_file_lines_fltrd[i]

        bg1 = first_info.strip().split()[1] # First .bedgraph
        bg2 = second_info.strip().split()[1] # Second .bedgraph
        tag = first_info.strip().split()[0]+'_'+second_info.strip().split()[0]

        df1_bed = pd.read_table(bg1, index_col=False, \
                na_values = 'NA', names = ['chrom1', 'start1', 'end1', 'value1'])
#        df1_bed.sort_values(by=['chrom1', 'start1']) # sort (doesn't change answer)    
        df1_bed['idx_value1'] = df1_bed['chrom1']+'_'+df1_bed['start1'].astype(str)
        df1_bed = df1_bed.set_index('idx_value1')

        df2_bed = pd.read_table(bg2, index_col=False, \
                na_values = 'NA', names = ['chrom2', 'start2', 'end2', 'value2'])
#        df2_bed.sort_values(by=['chrom2', 'start2']) # sort (doesn't change answer)
        df2_bed['idx_value2'] = df2_bed['chrom2']+'_'+df2_bed['start2'].astype(str)
        df2_bed = df2_bed.set_index('idx_value2')

        df_merged = df_merged = pd.concat( [ df1_bed, df2_bed ], axis=1, join='inner')
        df_merged['dcm'] = (df_merged['value2'] - df_merged['value1']).round(2)
#        df_merged['start1end'] = df_merged['start1'] + 1

        output_file = tag+'.bedgraph'
        cols_to_keep = ['chrom1', 'start1', 'end1', 'dcm']
        df_merged[cols_to_keep].to_csv(path_to_output+'/'+output_file, sep='\t', na_rep='NA', header=None, index=False)
        
        i += 1
        
        df1_bed.drop(df1_bed.index, inplace=True)
        df2_bed.drop(df2_bed.index, inplace=True)
        df1_bed.drop(df1_bed.index, inplace=True)
        del df1_bed, df2_bed, df_merged


def exec_proc_sample_wheel_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-ctp', action='store', dest='ctp_inp', required=True, help='Cell type Fofn')
    parser_required.add_argument('-expr', action='store', dest='expr_inp', required=True, help='Name of expression file')
    parser.add_argument('-pto', action='store', dest='pto_inp', default='.', help='Path to Output')
    parser._action_groups.reverse()
    return(parser)