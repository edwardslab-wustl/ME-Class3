import pandas as pd

def exec_filter(args):

    # Try keeping format fixed to ['gene_id', 'chrom', 'chromStart', 'chromEnd', 'strand', 'code_info', 'name', 'description']
    #
    expr_file = args.expr_inp
    annot_file = args.annot_inp
    output_annot_file_name = args.out_annot_inp 
    output_expr_file_name = args.out_expr_inp 
    annot_type = args.annot_type    
    #
    valid_chrom = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', \
             'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', \
             'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', \
             'chrX', 'chrY']
    #
    # Read expression data and set gene_id as index
    df_expr = ( pd.read_table(expr_file, index_col=False) ).set_index('gene_id')
    # Remove duplicate gene_id entries and keep first
    df_expr = df_expr.groupby(df_expr.index).first()
    # Read annotation data
    if annot_type == 'Gencode':
        df_annot = ( pd.read_table(annot_file, index_col=False, na_values = 'NA', names = \
                ['gene_id', 'chrom', 'chromStart', 'chromEnd', 'strand', 'code_info', 'name', 'description']) ).set_index('gene_id')
    if annot_type == 'RefSeq':
        df_annot = ( pd.read_table(annot_file, index_col=False, na_values = 'NA', names = \
                ['gene_id', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', \
                    'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']) ).set_index('gene_id')
    # Remove duplicate gene_id entries and keep firs
    df_annot = df_annot.groupby(df_annot.index).first()
    # Select only those for which expr daata is available
    df_filter = df_annot[(df_annot.index).isin(df_expr.index)]
    # add chr
    if annot_type == 'Gencode': 
        df_filter.loc[ :, 'chrom'] = 'chr' + df_filter['chrom'].astype(str)
    # chromosome filter
    df_filter = df_filter[(df_filter['chrom']).isin(valid_chrom)]
    # Write it to a file
    df_filter.to_csv(output_annot_file_name, sep='\t', na_rep='NA') #, header=None)
    #
    # filter and Write expression file after removing duplicates
    df_expr = df_expr[(df_expr.index).isin(df_filter.index)] 
    df_expr.to_csv(output_expr_file_name, sep='\t', na_rep='NA', float_format='%.3f') #, header=None)

def exec_filter_help(parser):
    #parser_required = parser.add_argument_group('required arguments') #then use parer_required for required args
    parser.add_argument('-expr', action='store', dest='expr_inp', help='Name of first file')
    parser.add_argument('-annot', action='store', dest='annot_inp', help='Name of second file')
    parser.add_argument('-oex', action='store', dest='out_expr_inp', help='Name of output expr file')
    parser.add_argument('-oan', action='store', dest='out_annot_inp', help='Name of output annot file')
    parser.add_argument('-apt', action='store', dest='annot_type', default='RefSeq', help='Gencode/RefSeq')
    #parser._action_groups.reverse()
    return(parser)
