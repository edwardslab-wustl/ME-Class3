
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt

from MEClass3.cluster_functions import select_features
from MEClass3.plot_functions import setup_x_axis


def read_interp_files(file_list):
    df = ''
    for file in file_list:
        if len(df) == 0:
            df = pd.read_csv(file, comment='#')
        else:
            df_tmp = pd.read_csv(file, comment='#')
            df_new = pd.concat([df, df_tmp], ignore_index=True)
            df = df_new
            del df_tmp
            del df_new
    return df

def normalize_labels(df):
    # Normalize up and down
    expr_flag_value_count = df['expr_flag'].value_counts()
    num_up = pd.Series(expr_flag_value_count)[1]
    num_dn = pd.Series(expr_flag_value_count)[-1]
    if num_up > num_dn:
        idx = df.index[ (df['expr_flag'] == 1) ]
        df.loc[np.random.choice(idx, size=(num_up-num_dn), replace=False), 'expr_flag'] = 0
        del idx
    elif num_dn > num_up:
        idx = df.index[ (df['expr_flag'] == -1) ]
        df.loc[np.random.choice(idx, size=(num_dn-num_up), replace=False), 'expr_flag'] = 0
        del idx
    return df
    
def plot_featureImportance(data, param_dict, param_data_dict, args):
    sns.set(font_scale=1.8)
    sns.set_style("ticks")
    plt.figure(figsize=(8,5)) 
    data_orig = data.copy(deep=True)
    for anno_id, param in param_dict.items():
        data = data_orig.copy(deep=True)
        param = param_dict[ anno_id ]
        [data_type,anno_type] = anno_id.split('-')
        data.set_index('feature')
        data_T = data.transpose(copy=True)
        data_T.columns = data_T.iloc[0]
        data_T = data_T[1:]
        y_cols = select_features(data_T, anno_type, data_type)
        y_data = data_T[y_cols].transpose(copy=True)
        y_data.reset_index(inplace=True)
        x_data = pd.DataFrame(param_data_dict[ anno_id ], columns=['x_info', 'x_data'])
        if anno_type == 'enh':
            y_data[['info','hue']] = y_data['feature'].str.split('_', n=1, expand=True)
            legend = 'full'
        else:
            y_data['hue'] = pd.Series(['signal' for x in range(len(y_data.index))])
            y_data['info'] = y_data.iloc[:,0]
            legend = None
        all_data = pd.merge(x_data,y_data, left_on='x_info', right_on='info')
        all_data.drop(['sum','feature','x_info'], axis=1, inplace=True)
        y_data_m = pd.melt(all_data, id_vars=['x_data','hue', 'info'])
        ax = sns.lineplot(  data=y_data_m,
                            x='x_data',
                            y='value',
                            hue='hue',
                            legend=legend )
        md_pt = 0
        plt.plot([md_pt,md_pt],[-1,1],'k-',alpha=0.5)
        plt.ylabel('Feature Importance')
        if args.featureImportance_max_y > 0:
            max_y = args.featureImportance_max_y
        else:
            max_y = 1.1 * y_data_m['value'].max()
        plt.ylim([0,max_y])
        setup_x_axis(plt, ax, param, anno_type, args)
        #if args.tight_layout:
        #    plt.tight_layout()
        plot_file = args.tag + f".featureImportance.{data_type}.{anno_type}.png"
        plt.savefig(plot_file, bbox_inches='tight')
        plt.close()
    return