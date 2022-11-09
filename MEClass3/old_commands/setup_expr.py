
import sys

import numpy as np
import pandas as pd

def exec_setup_expr(args):
    #
    interp_file = args.intrp_inp
    expr_file = args.expr_inp
    expr_floor_value = args.efv_inp
    #
    # Read in interpolation csv file.
    df_interp =  pd.read_csv(interp_file)
    # Add sample and gene_id column
    df_interp['gene_id'] = df_interp['gene_id-sample_name'].apply(lambda x: x.split('-')[0])
    df_interp['sample_name'] = df_interp['gene_id-sample_name'].apply(lambda x: x.split('-')[1])    
    # List of sample names
    sample_names =  list( set( df_interp['sample_name'].tolist() ) )
    cell_type_names = []
    for item in sample_names:
        cell_type_names.extend( item.strip().split('_') )
    #
    # Read expression data and set gene_id as index
    df_expr = ( pd.read_table(expr_file, index_col=False) ).set_index('gene_id')
    # Remove Duplicates
#    df_expr = df_expr.groupby(df_expr.index).first()
    # Floor expression values for selected cell types
    if args.floor_expr:    
        for cell_type in cell_type_names:
            df_expr.loc[ df_expr[cell_type] < expr_floor_value, cell_type ] = expr_floor_value
    #
    # Fold change in expression value [1]. [2].  [3]. log2(floored_Expr_B/floored_Expr_A)
    for sample in sample_names:
        if args.dexpr_flag == 0:
            df_expr[sample] = np.where( df_expr[ sample.split('_')[0] ] > df_expr[ sample.split('_')[1] ], \
                        -( df_expr[ sample.split('_')[0] ]/df_expr[ sample.split('_')[1] ] ), \
                        ( df_expr[ sample.split('_')[1] ]/df_expr[ sample.split('_')[0] ] ) )
            df_expr[sample] = np.where( df_expr[sample] == 1.0, 0.0, df_expr[sample] )
#            df_expr[sample] = np.where( np.logical_or((df_expr[sample] == 1.0), \
#                                (df_expr[sample] == -1.0)), 0.0, df_expr[sample] )
        #
        if args.dexpr_flag == 1: # 
            df_expr[sample] = np.where( df_expr[ sample.split('_')[0] ] > df_expr[ sample.split('_')[1] ], \
                        -( df_expr[ sample.split('_')[0] ]/df_expr[ sample.split('_')[1] ] ), \
                        ( df_expr[ sample.split('_')[1] ]/df_expr[ sample.split('_')[0] ] ) )
        #
        if args.dexpr_flag == 2:
            df_expr[sample] = np.log2( df_expr[ sample.split('_')[1] ]/df_expr[ sample.split('_')[0] ] )
            
        if args.dexpr_flag == 3: # 
            df_expr[sample] = np.where( df_expr[ sample.split('_')[0] ] > df_expr[ sample.split('_')[1] ], \
                        -( df_expr[ sample.split('_')[0] ]/df_expr[ sample.split('_')[1] ] ), \
                        ( df_expr[ sample.split('_')[1] ]/df_expr[ sample.split('_')[0] ] ) )
    #
    # processing of interf df
    #
    df_interp = df_interp.set_index('gene_id-sample_name')
    # Remove Duplicates
#    df_interp = df_interp.groupby(df_interp.index).first()
    #
    df_interp['expr_value'] = 0.0  # Float type
    df_interp['expr_flag'] = 0      # int type
    # 
    for gene_item in df_interp.index:
    #    print(gene_item)
        try:
            df_interp.loc[gene_item, 'expr_value'] = df_expr.loc[ gene_item.split('-')[0], gene_item.split('-')[1] ]
        except KeyError:
            sys.stdout.write(gene_item.split('-')[0]+' not found in expression file.\n')
            df_interp.loc[gene_item, 'expr_value'] = 0
        #
    if args.dexpr_flag == 0: # me-class demo
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] > 0, 1, df_interp['expr_flag'] )
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] < 0, -1, df_interp['expr_flag'] )
        df_interp['expr_flag'] =  df_interp['expr_flag'].astype(int)
    #
    if args.dexpr_flag == 1: # me-class paper
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] >= 2, 1, df_interp['expr_flag'] )
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] <= -2, -1, df_interp['expr_flag'] )
        df_interp['expr_flag'] =  df_interp['expr_flag'].astype(int)
    #
    if args.dexpr_flag == 2:
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] > 0.0, 1, df_interp['expr_flag'] )
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] < 0.0, -1, df_interp['expr_flag'] )
        df_interp['expr_flag'] = df_interp['expr_flag'].astype(int)
        
    if args.dexpr_flag == 3: # custom
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] >= args.expr_cutoff, 1, df_interp['expr_flag'] )
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] <= -args.expr_cutoff, -1, df_interp['expr_flag'] )
        df_interp['expr_flag'] =  df_interp['expr_flag'].astype(int)
    #
    #
    # Finally write data for classifier
    df_interp.to_csv('interp_expr_data.csv', sep=',')

def exec_setup_expr_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-intrp', action='store', dest='intrp_inp', required=True, help='Interpolation CSV file')
    parser_required.add_argument('-expr', action='store', dest='expr_inp', required=True, help='Expression file')
    parser.add_argument('-fef', action='store', dest='floor_expr', type=bool, default=True, help='Floor expression value?')
    parser.add_argument('-efv', action='store', dest='efv_inp', type=float, default=5.0, help='Expression floor value')
    parser.add_argument('-def', action='store', dest='dexpr_flag', type=int, default=1, help='Method for differential expression')
    parser.add_argument('--expr_cutoff', action='store', dest='expr_cutoff', type=int, default=2, help='fold change in expression cutoff, must use -def 3')
    parser._action_groups.reverse()
    return(parser)
