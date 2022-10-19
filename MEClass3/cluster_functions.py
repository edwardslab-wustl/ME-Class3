import sys

import sklearn.preprocessing 
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

class recursion_depth:
# from : https://www.codingem.com/python-maximum-recursion-depth/
    def __init__(self, limit):
        self.limit = limit
        self.default_limit = sys.getrecursionlimit()
    def __enter__(self):
        sys.setrecursionlimit(self.limit)
    def __exit__(self, type, value, traceback):
        sys.setrecursionlimit(self.default_limit)

def read_pred_file(file):
    result = pd.read_csv(file)
    return result

def normalize_expression(E):
    cm_norm = mpl.colors.Normalize(vmin=-1.,vmax=1.) 
    return cm_norm(E)

def subset_and_normalize_data( df, columns):
    result = df.loc[:,columns]
    result = sklearn.preprocessing.normalize(result)
    return result

def check_purity(Y):
    direction = 'up-reg'
    if (len(Y) > 0):
        up_cnt = 0
        dn_cnt = 0
        for y in Y:
            if y>0:
                up_cnt += 1
            elif y < 0:
                dn_cnt += 1
        purity = up_cnt / (up_cnt + dn_cnt)
        if purity < 0.5:
            purity = 1-purity
            direction = 'dn-reg'
    else:
        purity = -1
    return purity,direction

def cluster_plot_heatmap(df, norm_Y, linkage, cluster_tags, args):
    filename = args.out_base + ".meth.clustermap.png"
    title = "Delta Meth."
    meth_cmap = sns.diverging_palette(240,10,n=15,as_cmap=True)
    cluster_cmap = plt.get_cmap("cool")
    pred_cmap = plt.get_cmap("RdYlGn_r")
    # idea is to add column color for samples
    #sample_cmap = plt.get_cmap("gnuplot2")
    #sample_tags = [float(1)/(s+1) for s in S]    
    #row_colors = [sample_cmap(sample_tags),pred_cmap(norm_Y),cluster_cmap(cluster_tags)]
    row_colors = [pred_cmap(norm_Y),cluster_cmap(cluster_tags)]
    vmin = -args.color_max_meth_diff
    vmax = args.color_max_meth_diff
    #sns.set(style="white")
    sns.set(font_scale=2)
    with recursion_depth(5000):
        sns.clustermap( df, row_colors=row_colors,col_cluster = False,
                           figsize=(35,25),  method=args.linkage_method, row_linkage=linkage,
                           cmap=meth_cmap, linewidths = 0,
                           xticklabels=False,yticklabels=False,
                           vmin = vmin, vmax = vmax,
                           dendrogram_ratio=(.08, .2))
    plt.title(title)
    plt.savefig(filename)
    plt.close()   
    return

def print_individual_cluster_averages(uniq_clusters,fcluster,df,args):    
    cluster_info = []
    for cluster in uniq_clusters:
        idx = [ i for i,fc in enumerate(fcluster) if cluster == fc ]
        cluster_data = df.iloc[idx]
        purity,expression_direction = check_purity(cluster_data['expr_flag'])
        cluster_info.append("cluster: " + str(cluster) + 
                "; Expr. Dir: " + expression_direction +
                "; purity: " + str(purity) +
                "; num_genes: " + str(len(cluster_data['expr_flag'])) + "\n")
        if cluster_data.shape[0] >= args.min_genes_per_cluster and purity >= args.min_cluster_purity:
            #cluster_labels = [l for i,l in enumerate(cluster_data['gene_id']) if i in idx]
#            average_print_helper_meth_cpg(Xs,Cs,str(cluster),of_base,cluster_labels,purity,expression_direction,data_info,args)
            average_print_helper_meth_cpg(cluster_data,str(cluster),purity,expression_direction,args)
            outFile = (args.out_base+".meth_cpg.cluster_%s"%(cluster)+".csv")
            keep_cols= ['gene_id-sample_name','gene_id','sample_name',
                        'expr_value','expr_flag','prob_dn','prob_up','expr_pred']
            cluster_data[keep_cols].to_csv(outFile, index=False)
    return cluster_info

def average_print_helper_meth_cpg(data,cluster,purity,expression_direction,args):
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    plt.figure(figsize=(8,5)) 
    y_cols = [x for x in data if x.startswith('ftss')]
    y_data = data[y_cols].transpose(copy=True)
    y_data.reset_index(inplace=True)
    y_data_m = pd.melt(y_data, id_vars='index')
    y_data_m['hue'] = pd.Series(['signal' for x in range(len(y_data_m.index))])
    y_data_m_lab = replace_axis_labels(y_data_m, 'index')
    ax = sns.lineplot( data=y_data_m_lab,
                       x='index',
                       y='value',
                       ci=args.confidence_interval,
                       hue='hue',
                       legend=None )
    #lgd=plt.legend(loc=5,bbox_to_anchor=(1.7,0.5),handlelength=1,handletextpad=0.5)
    plt.title("Cluster: %s (n=%d, %s, purity=%5.2f)" % (cluster,data.shape[0],expression_direction,purity))
    plt.xlabel("Position relative to TSS (bp)")
    plt.ylabel(r'$\Delta$mCG/CG')
    md_pt = 0
    plt.plot([md_pt,md_pt],[-1,1],'k-',alpha=0.5)
    plt.ylim([-args.max_y_val,args.max_y_val])
    plt.yticks(np.arange(-args.max_y_val, args.max_y_val*1.05, step=0.2))
    plt.xlim([-5000,5000])
    plt.xticks((-5000,md_pt,5000),('-5kb','TSS','+5kb'))
    x_tick_minorLocator = MultipleLocator(1000)
    ax.xaxis.set_minor_locator(x_tick_minorLocator)
    if args.tight_layout:
        plt.tight_layout()
    plt.savefig(args.out_base+".meth_cpg.cluster_%s"%(cluster)+".png", bbox_inches='tight')
    #plt.savefig(args.out_base+".meth_cpg.cluster_%s"%(cluster)+".png", bbox_extra_artists=(lgd,), bbox_inches='tight')
    plt.close()
    return

def replace_axis_labels(df, column):
    new_df = df.copy(deep=True)
    bin_size = 20
    x_shift = 5000
    val_dict = dict()
    for i,x in enumerate(new_df[column].unique()):
        x_val = (i * bin_size) - x_shift
        new_df.loc[new_df[column] == x, column] = x_val
    return new_df
