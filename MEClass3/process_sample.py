import pandas as pd

def exec_proc_sample(args):
    
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
