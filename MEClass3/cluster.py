# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 11:30:15 2016

@author: cschlosberg & jredwards
"""
 
#import warnings
import argparse
from os.path import exists

import pandas as pd
import scipy.cluster.hierarchy    

from MEClass3.cluster_functions import read_pred_file
from MEClass3.cluster_functions import cluster_plot_heatmap
from MEClass3.cluster_functions import subset_and_normalize_data
from MEClass3.cluster_functions import print_individual_cluster_averages
from MEClass3.cluster_functions import normalize_expression
from MEClass3.cluster_functions import select_features
from MEClass3.io_functions import format_args_to_print
from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import eprint
from MEClass3.io_functions import read_params_from_interp_2

#warnings.filterwarnings("ignore", category=UserWarning)

def exec_cluster(args):
    with open(args.logfile, 'w') as log_FH:
        print_to_log(log_FH, format_args_to_print(args))
        pred_data = read_pred_file(args.pred_file)
        pred_data_filtered = pred_data[ \
            ((args.lowerPredBound <= pred_data['prob_up']) & (pred_data['prob_up'] <= args.upperPredBound)) | \
            ((args.lowerPredBound <= pred_data['prob_dn']) & (pred_data['prob_dn'] <= args.upperPredBound)) ]
        interp_data = ''
        found_data_flag = False
        param_data_dict = dict()
        param_dict = dict()
        for sample in pred_data_filtered['sample_name'].unique():
            sample_file = args.interp_data_path + '/' + sample + ".interp_expr_data.csv"
            if not exists(sample_file):
                sample_file = args.interp_data_path + '/' + sample + ".interp_expr_data.csv.gz"
            if exists(sample_file):
                print_to_log(log_FH, f"Reading in file: {sample_file}\n")
                found_data_flag = True
                param_data_dict, param_dict = read_params_from_interp_2(sample_file)
                if len(interp_data) == 0:
                    interp_data = pd.read_csv(sample_file, comment='#')
                    interp_data.drop(['gene_id','sample_name','expr_value','expr_flag'], axis=1, inplace=True)
                else:
                    tmp = pd.read_csv(sample_file, comment='#')
                    tmp.drop(['gene_id','sample_name','expr_value','expr_flag'], axis=1, inplace=True)
                    interp_data = pd.concat([interp_data, tmp], axis=0, ignore_index=True)
            else:
                print_to_log(log_FH, f"Skipping {sample}, can't find file: {sample_file}\n")
        if not found_data_flag:
            eprint(f"Can't find any interpolation data: {sample_file}!\nCheck file paths and names.\n")
            exit()
        merged_data = interp_data.merge(pred_data_filtered, on='gene_id-sample_name')
        merged_data.drop('clf_flag', axis=1, inplace=True)
        #merged_data.dropna()
        feat_cols = select_features(merged_data, args.anno_type, args.data_type)
        if len(feat_cols) == 0:
            eprint(f"could not find any features: {args.anno_type} {args.data_type}\nCheck --anno_type and --data_type match features in data files.\n")
            exit()
        merged_data_vals = subset_and_normalize_data(merged_data, feat_cols)
        norm_Y = normalize_expression(merged_data['expr_value']) 
        linkage = scipy.cluster.hierarchy.linkage(merged_data_vals,method=args.linkage_method,metric='euclidean')
        fcluster = scipy.cluster.hierarchy.fcluster(linkage,args.numClusters,criterion='maxclust')
        #cluster_tags = [float(1)/ct for ct in fcluster] 
        cluster_tags = [float(ct % 2) for ct in fcluster] 
        uniq_clusters = list(set(fcluster))
        cluster_plot_heatmap(merged_data_vals,norm_Y,linkage,cluster_tags,args)
        cluster_info = print_individual_cluster_averages(uniq_clusters,fcluster,merged_data,param_data_dict, param_dict, args)
        print_to_log(log_FH, "\n".join(cluster_info))
        return
        
        
def exec_cluster_help(parser):
    parser_required = parser.add_argument_group('required arguments')
#    parser_required.add_argument('interp_files', metavar='interp_files', type=str, nargs='+',
#        default=argparse.SUPPRESS,
#        help='interpolation files for classification')
#    parser_required.add_argument('-i', '--input_list', dest='input_list',
#        default=argparse.SUPPRESS,
#        required=True, help='Input list of sample names and file locations for pairings.')
    parser_required.add_argument('-p', '--pred_file', 
        default=argparse.SUPPRESS,
        required=True, help='Prediction file from classifiier.')
    parser_general = parser.add_argument_group('general arguments')
    parser_general.add_argument('--out_base', default='cluster_results', help='base name for output files and plots')
    parser_general.add_argument('--interp_data_path',
        default='intermediate_files', help='Path to directory with merged interpolation files from merge data step.')
    parser_general.add_argument('--anno_type', default="tss", 
        choices=["tss", "enh","all"],
        help='region or gene annotation file')
    parser_general.add_argument('--data_type', default="mC", 
        choices=["mC", "hmC", "other", "all"],
        help='type of data. all uses all data from all annotations and overrides --anno_type')
    parser_general.add_argument('--logfile', default='cluster.log', help='log file')
    #parser.add_argument('--cluster_data', default="both", choices=["meth_only", "hmC_only", "both"],
    #    help="which signatures/data to cluster based on. (default: both)")
    parser_heatmap = parser.add_argument_group('heatmap arguments')
    parser_heatmap.add_argument('--color_max_meth_diff',default=0.2,type=float,
        help="Sets color bar scale for heatmap.")
    parser_clustering = parser.add_argument_group('clustering arguments')
    parser_clustering.add_argument('--features', default="tss",
        choices=["tss", "enh"],
        help="Features to use for clustering.")
    parser_clustering.add_argument('--numClusters',default=3,type=int,
        help="number of clusters.")
    parser_clustering.add_argument('--linkage_method', default="complete", choices=["single", "complete", "average", "weighted", "ward", "median", "centroid"],
        help="linkage method for clustering. See scipy.cluster.hierarchy.linkage online documentation for more info.")
    parser_clustering.add_argument('--upperBound',default=2500,type=int,
        help="upperBound of window for clustering in number of features, in bp relative to TSS")
    parser_clustering.add_argument('--lowerBound',default=-500,type=int,
        help="lowerBound of window for clustering in number of features, in bp relative to TSS")
    parser_clustering.add_argument('--upperPredBound',default=1.0,type=float,
        help="upper prediction score bound for clustering.")
    parser_clustering.add_argument('--lowerPredBound',default=0.7,type=float,
        help="lower prediction score bound for clustering.")
    parser_clusterplots = parser.add_argument_group('cluster line plot arguments')
    parser_clusterplots.add_argument('--confidence_interval',default=95, type=int,
        help="confidence interval [0,100] for aggregate cluster plots.")
    parser_clusterplots.add_argument('--max_y_val',default=0.4, type=float,
        help="max y-value for plots.")
    parser_clusterplots.add_argument('--min_genes_per_cluster',default=10, type=int,
        help="min number of genes needed in a cluster to print it")
    parser_clusterplots.add_argument('--min_cluster_purity',default=0.75, type=float,
        help="min purity of the cluster to print it")
    parser_clusterplots.add_argument('--tight_layout', action='store_true',
        help="May help with layout of cluster plots if the plot does not fill the figure")
    #parser._action_groups.reverse()
    return parser
    
