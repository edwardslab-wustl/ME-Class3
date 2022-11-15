
import pandas as pd
#import numpy as np

from MEClass3.io_functions import read_params_from_interp_2
from MEClass3.cluster_functions import select_cluster_features
#from MEClass3.io_functions import eprint

def add_interp_header(file, header_list, args):
    with open(file, 'r') as FH:
        for line in FH:
            if line.startswith('#'):
                (comment, info, data) = line.strip().split()
                info2 = info.split(':')[1]
                anno_id = info2.split(',')[0] 
                (data_type, anno_type) = anno_id.split('-')
                if anno_type == 'enh':
                    new_info = ','.join([info,str(args.enh_lowerBound), str(args.enh_upperBound), args.enh_tag])
                elif anno_type == 'tss':
                    new_info = ','.join([info,str(args.tss_lowerBound), str(args.tss_upperBound)])
                new_line = " ".join([comment, new_info, data]) + "\n"
                header_list.append(new_line)
    return header_list

def add_tss_interp(df_interp, file, data_type, args):
    df_merged = ''
    #tss_interp = pd.read_csv(file, comment='#', dtype=np.float16, converters = {'gene_id-sample_name': str})
    tss_interp = pd.read_csv(file, comment='#')
    param_data_dict, param_dict = read_params_from_interp_2(file)
    feat_cols = select_cluster_features(tss_interp, param_data_dict, param_dict, 'tss', data_type, args)
    feat_cols.append('gene_id-sample_name')
    tss_interp = tss_interp.loc[:,feat_cols]
    if len(df_interp) == 0:
        df_merged = tss_interp
    else:
        df_merged = df_interp.merge(tss_interp, on='gene_id-sample_name')
    del tss_interp
    df_merged = df_merged.fillna(0)
    return df_merged

def fix_column_headers(df, tag):
    new_headers = []
    for num, header in enumerate(df.columns):
        if num == 0 or num == len(df.columns) - 1 or num == len(df.columns) - 2:
            new_headers.append(header)
        else:
            new_headers.append(header + "_" + tag)
    df.columns = new_headers
    return df
        
def add_enh_interp(df_interp, file, data_type, args):
    df_merged = ''
    #enh_interp_all = pd.read_csv(file, comment='#', dtype=np.float16, converters = {'enh_loc-gene_id-sample_name': str})
    enh_interp_all = pd.read_csv(file, comment='#') 
    param_data_dict, param_dict = read_params_from_interp_2(file)
    enh_interp_all['enh_loc'] = enh_interp_all['enh_loc-gene_id-sample_name'].apply(lambda x: x.split('-',1)[0])
    enh_interp_all['gene_id-sample_name'] = enh_interp_all['enh_loc-gene_id-sample_name'].apply(lambda x: x.split('-',1)[1])
    if len(df_interp) != 0:
        df_merged = df_interp
    for enh_loc in enh_interp_all['enh_loc'].unique():
        df_tmp_interp = enh_interp_all[enh_interp_all['enh_loc'] == enh_loc]
        df_tmp_interp = fix_column_headers(df_tmp_interp, enh_loc)
        feat_cols = select_cluster_features(df_tmp_interp, param_data_dict, param_dict, 'enh', data_type, args)
        feat_cols.extend(['enh_loc-gene_id-sample_name', 'gene_id-sample_name'])
        df_tmp_interp = df_tmp_interp.loc[:,feat_cols]
        df_tmp_interp.drop_duplicates(subset=None, keep='first', inplace=True)
        #eprint(f"enh_loc: {enh_loc}")
        #eprint(f"tmp shape: {df_tmp_interp.shape}")
        #eprint(f"merged shape: {df_merged.shape}")
        if len(df_merged) == 0:
            df_merged = df_tmp_interp
        else:
            tmp = df_tmp_interp.drop('enh_loc-gene_id-sample_name', axis=1)
            if 'enh_loc' in tmp.columns:
                tmp.drop('enh_loc', axis=1, inplace=True)
            df_merged = df_merged.merge(tmp, on='gene_id-sample_name')
            del tmp
        df_merged.drop_duplicates(subset=None, keep='first', inplace=True)
        #df_merged.to_csv(f"{file}.{enh_loc}.test")
        del feat_cols
        del df_tmp_interp
    return df_merged
