import sys
from collections import defaultdict

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
                           figsize=(35,25), row_linkage=linkage,
                           cmap=meth_cmap, linewidths = 0,
                           xticklabels=True,yticklabels=False,
                           vmin = vmin, vmax = vmax,
                           dendrogram_ratio=(.08, .2))
                           ##figsize=(35,25),  method=args.linkage_method, row_linkage=linkage,
    plt.title(title)
    plt.savefig(filename)
    plt.close()   
    return

def cluster_plot_heatmap(df, norm_Y, linkage, cluster_tags, param_dict, data_type, args):
    filename = args.out_base + f".{data_type}.clustermap.png"
    title = r'$\Delta$' + f"{data_type}"
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
    plt.title(title, fontsize=fontsize)
    left_side = 585 / (fig_width * 100)
    right_side = 3464 / (fig_width * 100)
    plt.figtext(left_side - 140 / (fig_width * 100), tick_height, "Expr.", fontsize=fontsize2, horizontalalignment='center', rotation='vertical')
    plt.figtext(left_side - 60 / (fig_width * 100), tick_height, "Cluster", fontsize=fontsize2, horizontalalignment='center', rotation='vertical')
    #data_type = 'mC'
    x_labels = defaultdict(list)
    if args.anno_type == 'tss':
        anno_id = data_type + '-' + args.anno_type
        anno_label = 'TSS'
        left_label = param_dict[anno_id].left_label()
        right_label = param_dict[anno_id].right_label()
        #region_size = param_dict[anno_id].size()
        #region_size = param_dict[anno_id].region_size
        x_labels = add_tss_labels(x_labels, left_side, right_side, left_label, right_label) 
        add_cat_labels(plt, ['TSS'], left_side, right_side, anno_label_height, fontsize2, rotate=False )
    elif args.anno_type == 'enh':
        anno_id = data_type + '-' + args.anno_type
        enh_labels = pull_enh_labels(df)
        left_label = param_dict[anno_id].left_label()
        right_label = param_dict[anno_id].right_label()
        #region_size = param_dict[anno_id].size()
        #region_size = param_dict[anno_id].region_size
        x_labels = add_enh_labels(x_labels, enh_labels, left_side, right_side, left_label, right_label )
        add_cat_labels(plt, enh_labels, left_side, right_side, anno_label_height, fontsize2, rotate=False )
    elif args.anno_type == 'all':
        num_tss_feat = len([ x for x in df.columns if x.startswith(data_type + '-' + 'tss')])
        num_enh_feat = len([ x for x in df.columns if x.startswith(data_type + '-' + 'enh')])
        enh_labels = pull_enh_labels(df)
        #num_enh_regions = len(enh_labels)
        #enh_region_size = param_dict[data_type + '-' + 'enh'].region_size
        #tss_region_size = param_dict[data_type + '-' + 'tss'].region_size
        #enh_region_size = param_dict[data_type + '-' + 'enh'].size()
        enh_left_label = param_dict[data_type + '-' + 'enh'].left_label()
        enh_right_label = param_dict[data_type + '-' + 'enh'].right_label()
        #tss_region_size = param_dict[data_type + '-' + 'tss'].size()
        tss_left_label = param_dict[data_type + '-' + 'tss'].left_label()
        tss_right_label = param_dict[data_type + '-' + 'tss'].right_label()
        #enh_start = left_side + (right_side - left_side) * (tss_region_size / (tss_region_size + enh_region_size * num_enh_regions)) 
        enh_start = left_side + (right_side - left_side) * (num_tss_feat / (num_tss_feat + num_enh_feat)) 
        #enh_start = enh_start - 27 / (fig_width * 100)
        x_labels = add_tss_labels(x_labels, left_side, enh_start, tss_left_label, tss_right_label) 
        x_labels = add_enh_labels(x_labels, enh_labels, enh_start, right_side, enh_left_label, enh_right_label )
        add_cat_labels(plt, enh_labels, enh_start, right_side, anno_label_height, fontsize2, rotate=True )
        add_cat_labels(plt, ['TSS'], left_side, right_side, anno_label_height, fontsize2, rotate=False )
    for x_boundary,label_list in x_labels.items():
        plt.figtext(x_boundary, tick_height, '|', fontsize=fontsize, horizontalalignment='center')
        plt.figtext(x_boundary, text_height, ",".join(label_list), fontsize=fontsize, horizontalalignment='center', rotation='vertical')
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

#def add_tss_cat_label (plt, left_side, right_side, anno_label, anno_label_height, fontsize2, label_rotation=False):
#    md_pt = (right_side - left_side ) / 2 + left_side
#    plt.figtext(md_pt, anno_label_height, anno_label, fontsize=fontsize2, horizontalalignment='center', rotation=label_rotation)
#    return

def add_tss_labels (labels, left_side, right_side, left_label, right_label):
    if int(left_label) < 0 and 0 < int(right_label):
        zero_label_scale_factor = -float(left_label) / (float(right_label) - float(left_label))
        zero_label_coord = left_side + zero_label_scale_factor * (right_side - left_side)
        labels = add_label(labels,"0", zero_label_coord)
    labels = add_label(labels,left_label, left_side)
    labels = add_label(labels,right_label, right_side)
    return labels

def add_label(dict, label, key):
    if key in dict:
        dict[key].insert(0, str(label))
    else:
        dict[key] = [str(label)]
    return dict

def add_cat_labels (plt, enh_labels, left_side, right_side, anno_label_height, fontsize2, rotate=False):
    label_rotation='horizontal'
    if rotate:
        label_rotation='vertical'
    num_regions = len(enh_labels)
    center_start = left_side + 0.5 * (right_side - left_side) / num_regions
    step_size = (right_side - left_side) / num_regions
    for i in range(0,num_regions):
        anno_label = enh_labels[i]
        x_position = center_start + step_size * i
        plt.figtext(x_position, anno_label_height, anno_label, fontsize=fontsize2, horizontalalignment='center', rotation=label_rotation)
    return
        
def add_enh_labels (labels, enh_labels, left_side, right_side, left_label, right_label):
    num_regions = len(enh_labels)
    #center_start = left_side + 0.5 * (right_side - left_side) / num_regions
    step_size = (right_side - left_side) / num_regions
    for i in range(0,num_regions):
        coord = left_side + step_size * i
        labels = add_label(labels,left_label, coord)
        coord2 = coord + step_size
        labels = add_label(labels,right_label, coord2)
        if left_label < 0 and 0 < right_label:
            center =  coord + step_size * -left_label / (right_label - left_label)
            labels = add_label(labels,"0", center)
    return labels

### OLD STUFF
    #char_width_corr = 10
    #plt.figtext(md_pt - char_width_corr * len(anno_label) / 3500, text_height, anno_label, fontsize=fontsize)
    #plt.figtext(left_side - char_width_corr * len(str(region_size +1)) / 3500, text_height, str(-region_size), fontsize=fontsize)
    #plt.figtext(right_side - char_width_corr * len(str(region_size)) / 3500, text_height, str(region_size), fontsize=fontsize)

def print_individual_cluster_averages(uniq_clusters, fcluster, df, param_data_dict, param_dict, args):    
    cluster_info = []
    data_list = []
    anno_list = []
    print_csv_flag = True
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
                    if print_csv_flag:
                        outFile = (args.out_base+f".cluster_{cluster}.csv")
                        keep_cols= ['gene_id-sample_name','gene_id','sample_name',
                                    'expr_value','expr_flag','prob_dn','prob_up','expr_pred']
                        cluster_data[keep_cols].to_csv(outFile, index=False)
                    average_print_helper_meth_cpg(cluster_data,str(cluster),purity,expression_direction,x_data,param_dict,anno_id,args)
            print_csv_flag = False
    return cluster_info

def select_features( df, anno_type, data_features):
    data_feature_list = pull_feature_list(data_features)
    final_col_list = []
    for data_type in data_feature_list:
        if anno_type == 'all':
            feat_tag = data_type
        else:
            feat_tag = data_type + '-' + anno_type
        select_cols = [ x for x in df if x.startswith(feat_tag) ]
        final_col_list.extend(select_cols)
    return final_col_list

def pull_feature_list ( features ):
    data_feature_list = []
    if features == 'all':
        data_feature_list.append('mC')
        data_feature_list.append('hmC')
    else:
        data_feature_list = [features]
    return data_feature_list

def select_cluster_features( df, param_data_dict, param_dict, features, data_features, args ):
    data_feature_list = pull_feature_list(data_features)
    select_cols = []
    for data_type in data_feature_list:
        feat_tag = data_type + '-' + features
        data_dict = dict(param_data_dict[feat_tag])
        param = param_dict[feat_tag]
        dn, up = pull_up_dn_bounds(param, features, args)
        select_cols_tmp = [ x for x in df if x.startswith(feat_tag) ]
        for feat in select_cols_tmp:
            test_feat =''
            if features == 'tss':
                test_feat = feat
            elif features == 'enh':
                [enh_feat, enh_info] = feat.split('_', 1)
                if args.enh_tag == 'all':
                    test_feat = enh_feat
                elif ',' in args.enh_tag and enh_info in args.enh_tag.split(','):
                    test_feat = enh_feat
                elif args.enh_tag == enh_info:
                    test_feat = enh_feat
                else:
                    continue
            if test_feat in data_dict and dn <= data_dict[test_feat] and data_dict[test_feat] <= up:
                select_cols.append(feat)
    return select_cols

def pull_up_dn_bounds(param, features, args):
    dn = -param.region_size
    up = param.region_size
    if features == 'tss':
        if args.tss_upperBound < up:
            up = args.tss_upperBound
        if args.tss_lowerBound > dn:
            dn = args.tss_lowerBound
    elif features == 'enh':
        if args.enh_upperBound < up:
            up = args.enh_upperBound
        if args.enh_lowerBound > dn:
            dn = args.enh_lowerBound
    if param.lowerBound !=0 or param.upperBound != 0:
        if param.upperBound < up:
            up = param.uppperBound
        if param.lowerBound > dn:
            dn = param.lowerBound
    return dn, up

def average_print_helper_meth_cpg(data,cluster,purity,expression_direction,x_data,param_dict,anno_id,args):
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    plt.figure(figsize=(8,5)) 
    y_cols = select_features(data, args.anno_type, args.data_type)
    y_data = data[y_cols].transpose(copy=True)
    y_data.reset_index(inplace=True)
    param = param_dict[ anno_id ]
    #x_range = [-param.region_size,param.region_size]
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
    setup_x_axis(plt, ax, param, anno_type, args)
    if args.tight_layout:
        plt.tight_layout()
    plot_file = args.out_base + f".cluster_{cluster}.{data_type}.{anno_type}.png"
    plt.savefig(plot_file, bbox_inches='tight')
    plt.close()
    return

def setup_x_axis(figure, axis, param, anno_type, args):
    x_range = []
    if param.lowerBound !=0 or param.upperBound !=0:
        x_range = [param.lowerBound,param.upperBound]
    else:
        x_range = [-param.region_size,param.region_size]
    figure.xlim(x_range)
    if anno_type == 'tss':
        feat_label = "TSS"
        x_tick_minorLocator = MultipleLocator(1000)
        axis.xaxis.set_minor_locator(x_tick_minorLocator)
    elif anno_type == 'enh':
        feat_label = "Enhancer"
        plt.legend(loc='right', title='',bbox_to_anchor=(1.4,0.5),handlelength=1,handletextpad=0.5,frameon=False)
        x_tick_minorLocator = MultipleLocator(100)
        axis.xaxis.set_minor_locator(x_tick_minorLocator)
    else:
        feat_label = anno_type
        x_tick_minorLocator = MultipleLocator(100)
        axis.xaxis.set_minor_locator(x_tick_minorLocator)
    axis.set_xlabel(f"Distance to {feat_label} (bp)")
    figure.xticks((x_range[0],0,x_range[1]))
    return
    #return figure, axis

def replace_axis_labels(df, column):
    new_df = df.copy(deep=True)
    bin_size = 20
    x_shift = 5000
    #val_dict = dict()
    for i,x in enumerate(new_df[column].unique()):
        x_val = (i * bin_size) - x_shift
        new_df.loc[new_df[column] == x, column] = x_val
    return new_df

def check_features(anno, feat):
    return_flag = True
    if anno == 'all':
        return_flag = True
    elif anno == feat:
        return_flag = True
    else:
        return_flag = False
    return return_flag
