
import numpy as np
import pandas as pd

def add_interp_header(file, header_list):
    with open(file, 'r') as FH:
        for line in FH:
            if line.startswith('#'):
                header_list.append(line)
    return header_list

def add_tss_interp(df_interp, file):
    df_merged = ''
    tss_interp = pd.read_csv(file, comment='#')
    if len(df_interp) == 0:
        df_merged = tss_interp
    else:
        #df_merged = pd.concat([df_interp, tss_interp], axis=1)
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
        
def add_enh_interp(df_interp, file):
    df_merged = ''
    enh_interp_all = pd.read_csv(file, comment='#')
    enh_interp_all['enh_loc'] = enh_interp_all['enh_loc-gene_id-sample_name'].apply(lambda x: x.split('-',1)[0])
    enh_interp_all['gene_id-sample_name'] = enh_interp_all['enh_loc-gene_id-sample_name'].apply(lambda x: x.split('-',1)[1])
    if len(df_interp) != 0:
        df_merged = df_interp
    for enh_loc in enh_interp_all['enh_loc'].unique():
        df_tmp_interp = enh_interp_all[enh_interp_all['enh_loc'] == enh_loc]
        df_tmp_interp = fix_column_headers(df_tmp_interp, enh_loc)
        if len(df_merged) == 0:
            df_merged = df_tmp_interp
        else:
            tmp = df_tmp_interp.drop('enh_loc-gene_id-sample_name', axis=1)
            tmp.drop('enh_loc', axis=1, inplace=True)
            df_merged = df_merged.merge(tmp, on='gene_id-sample_name')
            del tmp
        del df_tmp_interp
    return df_merged
