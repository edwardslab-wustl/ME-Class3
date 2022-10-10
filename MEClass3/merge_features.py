
import pandas as pd

def exec_merge_features(args):
    
    # Merge features.
    p_sidx = args.pst_inp # 1
    p_eidx = args.ped_inp # 2
    n_sidx = args.nst_inp # 1
    n_eidx = args.ned_inp # 2
    tss_interp_file = args.tif_inp
    
    # df interp for tss
    df_tss = pd.read_csv(tss_interp_file).set_index('gene_id-sample_name')

    # open df for proximity list of regulatory elements and merge
    df_merged = df_tss

    # Merge n-end
    if (n_sidx != -1):
        dict_dfn = {}
        for i in range(n_sidx, n_eidx+1):
            dict_dfn['n_end_'+str(i)] = pd.read_csv('re_n_'+str(i)+'_interp.csv').set_index('gene_id-sample_name')
            dict_dfn['n_end_'+str(i)] = dict_dfn['n_end_'+str(i)][ (dict_dfn['n_end_'+str(i)].index).isin(df_tss.index) ]
            #
            df_merged = pd.concat( [ df_merged, dict_dfn['n_end_'+str(i)] ], axis=1 ) #, join='inner')
            # clear memory
            dict_dfn['n_end_'+str(i)] = dict_dfn['n_end_'+str(i)].iloc[0:0]
            del dict_dfn['n_end_'+str(i)]


    # Merge p-end
    if (p_sidx != -1):
        dict_dfp = {}
        for i in range(p_sidx, p_eidx+1):
            dict_dfp['p_end_'+str(i)] = pd.read_csv('re_p_'+str(i)+'_interp.csv').set_index('gene_id-sample_name')
            dict_dfp['p_end_'+str(i)] = dict_dfp['p_end_'+str(i)][ (dict_dfp['p_end_'+str(i)].index).isin(df_tss.index) ]
            #
            df_merged = pd.concat( [ dict_dfp['p_end_'+str(i)], df_merged ], axis=1 ) #, join='inner')
            # clear memory
            dict_dfp['p_end_'+str(i)] = dict_dfp['p_end_'+str(i)].iloc[0:0]
            del dict_dfp['p_end_'+str(i)]

    #fill all NaN to 0.0 and assign name to index
    df_merged = df_merged.fillna(0)
    df_merged.index.name = 'gene_id-sample_name'

    #Write output to a file.
    df_merged.to_csv('interp_merged.csv', sep=',')
