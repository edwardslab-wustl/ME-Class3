
import sys
import argparse

import pandas as pd
import numpy as np
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score

def exec_eval(args):

    df = pd.read_csv( args.dfi_inp, index_col=[0] ).reset_index()
    steps = args.steps_inp
    tag = ((args.dfi_inp).strip().split('.'))[0]
    output_file_name = ((args.dfi_inp).strip('.csv'))+'_acc_rej.txt'
    df['expr_pred'] = df['expr_pred'].astype(int)
    totalGenes = 0
    totalGenes_P = 0
    totalGenes_N = 0
    totalGenes = df['gene_id-sample_name'].count()
    expr_flag_value_count = df['expr_flag'].value_counts()
    try:
            totalGenes_P = pd.Series(expr_flag_value_count)[1]
    except KeyError:
            totalGenes_P = 0
    try:
            totalGenes_N = pd.Series(expr_flag_value_count)[-1]
    except KeyError:
            totalGenes_N = 0
    output_file = open(output_file_name,'w')
    sample_names = list( set( df['sample_name'].tolist() ) )
    init_gene_values = [0] * len( sample_names )
    sample_freq_90 = sample_freq_85 = sample_freq_80 = sample_freq_75  = pd.Series( init_gene_values, index=sample_names ) 
    gene_count_90_per_acc, gene_count_85_per_acc, gene_count_80_per_acc, gene_count_75_per_acc = 0, 0, 0, 0
    for threshold in np.linspace(0,0.5,steps):
        TP_P = 0
        TP_N = 0
        numGenes_P = 0
        numGenes_N = 0
        df_tmp = df[ (df['prob_up'] >= 1-threshold) | (df['prob_dw'] >= 1-threshold) ]
        df_tmp_P = df_tmp[ (df_tmp['prob_up'] >= 1-threshold) ]
        numGenes_P = df_tmp_P['gene_id-sample_name'].count()
        df_tmp_P = df_tmp_P[ (df_tmp_P['expr_flag'] == 1) & (df_tmp_P['expr_pred'] == 1) ]
        TP_P = df_tmp_P['gene_id-sample_name'].count()
        df_tmp_N = df_tmp[ (df_tmp['prob_dw'] >= 1-threshold) ]
        numGenes_N = df_tmp_N['gene_id-sample_name'].count()
        df_tmp_N = df_tmp_N[ ( df_tmp_N['expr_flag'] == -1 ) & (df_tmp_N['expr_pred'] == -1) ]
        TP_N = df_tmp_N['gene_id-sample_name'].count()
        #sys.stdout.write( str(threshold)+'\t'+str(numGenes_P)+'\t'+str(TP_P)+'\t'+str(numGenes_N)+'\t'+str(TP_N)+'\n')
        accuracy = 0
        rejectRate = 0
        if (numGenes_P + numGenes_N) > 0:
                accuracy = float(TP_P + TP_N) / float(numGenes_P + numGenes_N)
        if totalGenes > 0:
                rejectRate = float(numGenes_P + numGenes_N) / float(totalGenes)
        if ( accuracy >= 0.90) :
                gene_count_90_per_acc = ( numGenes_P + numGenes_N )
                sample_freq_90 = df_tmp['sample_name'].value_counts()
        if ( accuracy >= 0.85) :
                gene_count_85_per_acc = ( numGenes_P + numGenes_N )
                sample_freq_85 = df_tmp['sample_name'].value_counts()
        if ( accuracy >= 0.80) :
                gene_count_80_per_acc = ( numGenes_P + numGenes_N )
                sample_freq_80 = df_tmp['sample_name'].value_counts()
        if ( accuracy >= 0.75) :
                gene_count_75_per_acc = ( numGenes_P + numGenes_N )
                sample_freq_75 = df_tmp['sample_name'].value_counts()
        df_tmp = df_tmp.iloc[0:0]
        df_tmp_P = df_tmp_P.iloc[0:0]
        df_tmp_N = df_tmp_N.iloc[0:0]
        output_file.write(str(threshold)+'\t'+str(rejectRate)+'\t'+str(accuracy)+'\t'+str( numGenes_P + numGenes_N )+'\n')

    output_file.close()
    y_true = df['expr_flag'].tolist()
    y_pred = df['expr_pred'].tolist()
    sys.stdout.write( 'F1 Score (macro):\t'+str( (f1_score(y_true, y_pred, average='macro') ))+'\n' )
    sys.stdout.write( 'F1 Score (micro):\t'+str( (f1_score(y_true, y_pred, average='micro') ))+'\n' )
    sys.stdout.write( 'Accuracy score:\t'+str( accuracy_score(y_true, y_pred) )+'\n' )
    sys.stdout.write( '#Genes with >90% accuracy:\t'+str( gene_count_90_per_acc )+'\n' )
    with open( tag+'_genes_with_90_percent_accuracy.txt', 'w') as output_90per:
        for sample in sample_names:
            try:
                output_90per.write( str(pd.Series( sample_freq_90 )[sample]) +'\t'+ sample +'\n' )
            except KeyError:
                output_90per.write( str(0) +'\t'+ sample +'\n' )
    sys.stdout.write( '#Genes with >85% accuracy:\t'+str( gene_count_85_per_acc )+'\n' )
    with open( tag+'_genes_with_85_percent_accuracy.txt', 'w') as output_85per:
        for sample in sample_names:
            try:
                output_85per.write( str(pd.Series( sample_freq_85 )[sample]) +'\t'+ sample +'\n' )
            except KeyError:
                output_85per.write( str(0) +'\t'+ sample +'\n' )
    sys.stdout.write( '#Genes with >80% accuracy:\t'+str( gene_count_80_per_acc )+'\n' )
    with open( tag+'_genes_with_80_percent_accuracy.txt', 'w') as output_80per:
        for sample in sample_names:
            try:
                output_80per.write( str(pd.Series( sample_freq_80 )[sample]) +'\t'+ sample +'\n' )
            except KeyError:
                output_80per.write( str(0) +'\t'+ sample +'\n' )
    sys.stdout.write( '#Genes with >75% accuracy:\t'+str( gene_count_75_per_acc )+'\n' )
    with open( tag+'_genes_with_75_percent_accuracy.txt', 'w') as output_75per:
        for sample in sample_names:
            try:
                output_75per.write( str(pd.Series( sample_freq_75 )[sample]) +'\t'+ sample +'\n' )
            except KeyError:
                output_75per.write( str(0) +'\t'+ sample +'\n' )

    sys.stdout.write( '#Total_Genes:\t'+str( totalGenes )+'\n' )

def exec_eval_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('--dfi', action='store', dest='dfi_inp',
        required=True,
        default=argparse.SUPPRESS,
        help='Dataframe output from classification')
    parser.add_argument('--nstp', action='store', dest='steps_inp',
        type=int, default=101, help='Number of steps in acc_reject rate curve')
    parser._action_groups.reverse()
    return(parser)
