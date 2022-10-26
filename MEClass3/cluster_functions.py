import sys

import sklearn.preprocessing 
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
#import matplotlib.gridspec

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
    result = pd.DataFrame(result, columns=columns)
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

def cluster_plot_heatmap_old(df, norm_Y, linkage, cluster_tags, args):
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
                           xticklabels=True,yticklabels=False,
                           vmin = vmin, vmax = vmax,
                           dendrogram_ratio=(.08, .2))
                           #xticklabels=False,yticklabels=False,
    ##plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0) cg = sns.clustermap
    plt.title(title)
    plt.savefig(filename)
    plt.close()   
    return

#def cluster_plot_heatmap(df, norm_Y, linkage, cluster_tags, param_dict, param_data_dict, args):
def cluster_plot_heatmap(df, norm_Y, linkage, cluster_tags, param_dict, args):
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
    fig_width = 35
    fontsize=30
    fontsize2=40
    tick_height = 0.8
    text_height = 0.82
    anno_label_height = 0.86
    with recursion_depth(5000):
        sns.clustermap( df, row_colors=row_colors,col_cluster = False,
                           figsize=(fig_width,25),  method=args.linkage_method, row_linkage=linkage,
                           cmap=meth_cmap, linewidths = 0,
                           xticklabels=False,yticklabels=False,
                           vmin = vmin, vmax = vmax,
                           dendrogram_ratio=(.10, .2))
                           #xticklabels=False,yticklabels=False,
    ##plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0) cg = sns.clustermap
    plt.title(title, fontsize=fontsize)
    left_side = 585 / (fig_width * 100)
    right_side = 3464 / (fig_width * 100)
    plt.figtext(left_side - 140 / (fig_width * 100), tick_height, "Expr.", fontsize=fontsize2, horizontalalignment='center', rotation='vertical')
    plt.figtext(left_side - 60 / (fig_width * 100), tick_height, "Cluster", fontsize=fontsize2, horizontalalignment='center', rotation='vertical')
    data_type = 'mC'
    if args.anno_type == 'tss':
        anno_id = data_type + '-' + args.anno_type
        anno_label = 'TSS'
        region_size = param_dict[anno_id].region_size
        add_tss_labels(plt, anno_label, left_side, right_side, region_size, text_height, anno_label_height, tick_height, fontsize, fontsize2 )
    elif args.anno_type == 'enh':
        anno_id = data_type + '-' + args.anno_type
        enh_labels = pull_enh_labels(df)
        region_size = param_dict[anno_id].region_size
        add_enh_labels(plt, enh_labels, left_side, right_side, region_size, text_height, anno_label_height, tick_height, fontsize, fontsize2 )
    elif args.anno_type == 'all':
        enh_labels = pull_enh_labels(df)
        num_enh_regions = len(enh_labels)
        enh_region_size = param_dict[data_type + '-' + 'enh'].region_size
        tss_region_size = param_dict[data_type + '-' + 'tss'].region_size
        enh_start = left_side + (right_side - left_side) * (2 * tss_region_size / (2 * tss_region_size + enh_region_size * num_enh_regions)) 
        enh_start = enh_start - 27 / (fig_width * 100)
        add_tss_labels(plt, 'TSS', left_side, enh_start, tss_region_size, text_height, anno_label_height, tick_height, fontsize, fontsize2, rotate=False )
        add_enh_labels(plt, enh_labels, enh_start, right_side, enh_region_size, text_height, anno_label_height, tick_height, fontsize, fontsize2, rotate=True )
    plt.savefig(filename)
    plt.close()   
    return

def pull_enh_labels(df):
    result = []
    for feat in df.columns:
        if feat.split('-')[1].startswith('enh'):
            y = feat.split('_', 1)[1]
            if y not in result:
                result.append(y) 
    return result

def add_tss_labels (plt, anno_label, left_side, right_side, region_size, text_height, anno_label_height, tick_height, fontsize, fontsize2, rotate=False):
    label_rotation='horizontal'
    if rotate:
        label_rotation='vertical'
    md_pt = (right_side - left_side ) / 2 + left_side
    plt.figtext(md_pt, text_height, '0', fontsize=fontsize, horizontalalignment='center', rotation='vertical')
    plt.figtext(md_pt, anno_label_height, anno_label, fontsize=fontsize2, horizontalalignment='center', rotation=label_rotation)
    plt.figtext(left_side, text_height, str(-region_size), fontsize=fontsize, horizontalalignment='center', rotation='vertical')
    plt.figtext(right_side, text_height, str(region_size), fontsize=fontsize, horizontalalignment='center', rotation='vertical')
    tick_start = left_side
    tick_end = right_side
    tick_step = (right_side - left_side) / 10
    i = tick_start
    while i <= tick_end:
        plt.figtext(i, tick_height, '|', fontsize=fontsize, horizontalalignment='center')
        i += tick_step
    return

def add_enh_labels (plt, enh_labels, left_side, right_side, region_size, text_height, anno_label_height, tick_height, fontsize, fontsize2, rotate=False):
    label_rotation='horizontal'
    if rotate:
        label_rotation='vertical'
    num_regions = len(enh_labels)
    center_start = left_side + 0.5 * (right_side - left_side) / num_regions
    step_size = (right_side - left_side) / num_regions
    for i in range(0,num_regions):
        label = "0" 
        anno_label = enh_labels[i]
        x_position = center_start + step_size * i
        plt.figtext(x_position, tick_height, '|', fontsize=fontsize, horizontalalignment='center')
        plt.figtext(x_position, text_height, label, fontsize=fontsize, horizontalalignment='center', rotation='vertical')
        plt.figtext(x_position, anno_label_height, anno_label, fontsize=fontsize2, horizontalalignment='center', rotation=label_rotation)
        x_boundary = left_side + step_size * i
        plt.figtext(x_boundary, tick_height, '|', fontsize=fontsize, horizontalalignment='center')
        plt.figtext(x_boundary, text_height, str(region_size), fontsize=fontsize, horizontalalignment='center', rotation='vertical')
    x_boundary = left_side + step_size * num_regions
    plt.figtext(x_boundary, tick_height, '|', fontsize=fontsize, horizontalalignment='center')
    plt.figtext(x_boundary, text_height, str(region_size), fontsize=fontsize, horizontalalignment='center', rotation='vertical')
    return

### OLD STUFF
    #char_width_corr = 10
    #plt.figtext(md_pt - char_width_corr * len(anno_label) / 3500, text_height, anno_label, fontsize=fontsize)
    #plt.figtext(left_side - char_width_corr * len(str(region_size +1)) / 3500, text_height, str(-region_size), fontsize=fontsize)
    #plt.figtext(right_side - char_width_corr * len(str(region_size)) / 3500, text_height, str(region_size), fontsize=fontsize)
        
    #heatmap.gs.update(left=0.05, bottom=0.1)
    #gs2 = matplotlib.gridspec.GridSpec(1,1, left=0.2, top=0.10)
    #ax2 = heatmap.fig.add_subplot(gs2[0])
    #region_size = 5000
    #plt.plot([-region_size,region_size],[1,1],'k-',alpha=0.5)
    #ax2.set_facecolor('white')
    #ax2.get_yaxis().set_visible(False)
    #ax2.get_xaxis().set_visible(True)
    #ax2.set_xlim([-region_size,region_size])
    #feat_label ='TSS'
    #ax2.xaxis.set_tick_params(width=5, size=10)
    #x_tick_minorLocator = MultipleLocator(1000)
    #ax2.xaxis.set_minor_locator(x_tick_minorLocator)
    #ax2.set_xticks((-region_size,0,region_size),(str(region_size),feat_label,str(region_size)))
    #md_pt = 0
    #plt.plot([md_pt,md_pt],[-1,1],'k-',alpha=0.5)
    #ax2.set_ylim([-0.1,0.1])

def print_individual_cluster_averages(uniq_clusters, fcluster, df, param_data_dict, param_dict, args):    
    cluster_info = []
    data_list = []
    anno_list = []
    if args.data_type == 'all':
        data_list = ['mC', 'hmC']
    else:
        data_list = [ args.data_type ]
    if args.anno_type == 'all':
        anno_list = ['tss', 'enh']
    else:
        anno_list = [ args.anno_type ]
    for data_type in data_list:
        for anno_type in anno_list:
            anno_id = data_type + '-' + anno_type
            x_data = pd.DataFrame(param_data_dict[ anno_id ], columns=['info', 'x_data'])
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
                    average_print_helper_meth_cpg(cluster_data,str(cluster),purity,expression_direction,x_data,param_dict,anno_id,args)
                    outFile = (args.out_base+".meth_cpg.cluster_%s"%(cluster)+".csv")
                    keep_cols= ['gene_id-sample_name','gene_id','sample_name',
                                'expr_value','expr_flag','prob_dn','prob_up','expr_pred']
                    cluster_data[keep_cols].to_csv(outFile, index=False)
    return cluster_info

def select_features( df, anno_type, data_type):
    if data_type == 'all':
        select_cols = [ x for x in df ]
    elif anno_type == 'all':
        feat_tag = data_type
        select_cols = [ x for x in df if x.startswith(feat_tag) ]
    else:
        feat_tag = data_type + '-' + anno_type
        select_cols = [ x for x in df if x.startswith(feat_tag) ]
    return select_cols

def average_print_helper_meth_cpg(data,cluster,purity,expression_direction,x_data,param_dict,anno_id,args):
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    plt.figure(figsize=(8,5)) 
    y_cols = select_features(data, args.anno_type, args.data_type)
    y_data = data[y_cols].transpose(copy=True)
    y_data.reset_index(inplace=True)
    param = param_dict[ anno_id ]
    x_range = [-param.region_size,param.region_size]
    [data_type,anno_type] = anno_id.split('-')
    if anno_type == 'enh':
        y_data[['info','hue']] = y_data['index'].str.split('_', n=1, expand=True)
        legend = 'full'
    else:
        y_data['hue'] = pd.Series(['signal' for x in range(len(y_data.index))])
        y_data['info'] = y_data.iloc[:,0]
        legend = None
    all_data = pd.merge(x_data,y_data)
    y_data_m = pd.melt(all_data, id_vars=['x_data','hue', 'index', 'info'])
    ax = sns.lineplot( data=y_data_m,
                       x='x_data',
                       y='value',
                       ci=args.confidence_interval,
                       hue='hue',
                       legend=legend )
    plt.title("Cluster: %s (n=%d, %s, purity=%5.2f)" % (cluster,data.shape[0],expression_direction,purity))
    if anno_type == 'tss':
        ax.set_xlabel("Distance to TSS (bp)")
    elif anno_type == 'enh':
        ax.set_xlabel("Distance to Enhancer (bp)")
        plt.legend(loc='right', title='',bbox_to_anchor=(1.4,0.5),handlelength=1,handletextpad=0.5,frameon=False)
    else:
        ax.set_xlabel("Relative Position (bp)")
    plt.ylabel(r'$\Delta$' + data_type)
    md_pt = 0
    plt.plot([md_pt,md_pt],[-1,1],'k-',alpha=0.5)
    plt.ylim([-args.max_y_val,args.max_y_val])
    plt.yticks(np.arange(-args.max_y_val, args.max_y_val*1.05, step=0.2))
    plt.xlim(x_range)
    if anno_type == 'tss':
        feat_label = "TSS"
        x_tick_minorLocator = MultipleLocator(1000)
        ax.xaxis.set_minor_locator(x_tick_minorLocator)
    elif anno_type == 'enh':
        feat_label = "Enhancer"
        x_tick_minorLocator = MultipleLocator(100)
        ax.xaxis.set_minor_locator(x_tick_minorLocator)
    else:
        feat_label = anno_type
        x_tick_minorLocator = MultipleLocator(100)
        ax.xaxis.set_minor_locator(x_tick_minorLocator)
    plt.xticks((-param.region_size,0,param.region_size),(str(param.region_size),feat_label,str(param.region_size)))
    if args.tight_layout:
        plt.tight_layout()
    plot_file = args.out_base + f".meth_cpg.cluster_{cluster}.{data_type}.{anno_type}.png"
    plt.savefig(plot_file, bbox_inches='tight')
    plt.close()
    return

def replace_axis_labels(df, column):
    new_df = df.copy(deep=True)
    bin_size = 20
    x_shift = 5000
    #val_dict = dict()
    for i,x in enumerate(new_df[column].unique()):
        x_val = (i * bin_size) - x_shift
        new_df.loc[new_df[column] == x, column] = x_val
    return new_df
