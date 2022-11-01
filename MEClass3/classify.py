import argparse
import resource 

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold

from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import format_args_to_print
from MEClass3.io_functions import read_params_from_interp_2
from MEClass3.sample import read_sample_file
from MEClass3.classify_functions import read_interp_files
from MEClass3.classify_functions import normalize_labels
from MEClass3.classify_functions import plot_featureImportance


def exec_classify(args):
    pair_list = read_sample_file(args.input_list)
    df = read_interp_files(args.interp_files)
#    df = ( pd.read_csv(args.dfi_inp) ).set_index('gene_id-sample_name')
    df = df[~df.isin([np.nan, np.inf, -np.inf]).any(axis=1)]
    #set up columns to store classifier performance
    df['prob_dn'] = 0.0
    df['prob_up'] = 0.0
    df['expr_pred'] = np.nan
    # Add column to store a train and test flag and mark all null
    df['clf_flag'] = 'null'
    df.set_index('gene_id-sample_name')
    #set up dataframe for feature importances
    drop_col_list = ['gene_id-sample_name','sample_name','expr_value','expr_flag','prob_up','prob_dn','expr_pred','clf_flag','gene_id']
    feature_column_names = [ x for x in df.columns if x not in drop_col_list ]
    df_fi = pd.DataFrame(columns=feature_column_names)
    #===============================================
    with open(args.logfile, 'w') as log_FH:
        print_to_log(log_FH, format_args_to_print(args))
        shuffle=True
        if args.no_shuffle:
            shuffle=False
        if len(pair_list) == 1:
            if args.gnorm:
                df = normalize_labels(df)
            gene_id_list = np.asarray( list(set( df['gene_id'].tolist() ) ))
            kf = KFold(n_splits=args.folds, shuffle=shuffle ) 
            kf.get_n_splits(gene_id_list)
            for train_index, test_index in kf.split(gene_id_list): # K-Fold gene split
                I_train, I_test = gene_id_list[train_index], gene_id_list[test_index]
                df_kf = df.copy()
                idx = df_kf.index[ ~df_kf['gene_id'].isin(I_test) ]
                df_kf.loc[idx, 'clf_flag'] = 'train'
                del idx
                idx = df_kf.index[ df_kf['gene_id'].isin(I_test) ]
                df_kf.loc[idx, 'clf_flag'] = 'test'
                del idx
                # Select traing data set
                df_kf_train = df_kf[ df_kf['clf_flag']=='train' ].copy(deep=True)
                df_kf_test = df_kf[ df_kf['clf_flag']=='test' ].copy(deep=True)
                y_train = df_kf_train['expr_flag'].values
                #Setup classifier and run it
                df_kf_train.drop(drop_col_list, axis=1, inplace=True)
                df_kf_test.drop(drop_col_list, axis=1, inplace=True)
                clf = RandomForestClassifier(n_estimators=args.num_trees, n_jobs=args.threads) 
                clf.fit(df_kf_train, y_train)
                y_test_prob = clf.predict_proba(df_kf_test)
                y_test_pred = clf.predict(df_kf_test)
                # store performance values in df_out
                for i,idx in enumerate(df_kf_test.index):
                    df.loc[idx, 'expr_pred'] = y_test_pred[i]
                    df.loc[idx, 'prob_dn'] = y_test_prob[i][0]
                    df.loc[idx, 'prob_up'] = y_test_prob[i][1]
                # Feature importance
                if args.featureImportance:
                    df_fi = pd.concat( [df_fi, (pd.DataFrame( [(clf.feature_importances_)],  columns=feature_column_names))] )
                df_kf.drop(df_kf.index, inplace=True)
                df_kf_train.drop(df_kf_train.index, inplace=True)
                df_kf_test.drop(df_kf_test.index, inplace=True)
                del df_kf, df_kf_train, df_kf_test, clf
        #===============================================
        else: # loso loop
            for sample_pair in pair_list:
                #sample = sample_pair.name
                df_loso = df.copy()
                # test sample
                sample_name_test_list = []
                sample_name_test_list.append(sample_pair.name)
                print_to_log(log_FH, ", ".join(sample_name_test_list)+'\n')
                # train samples
                sample_name_train_list = []
                for sample_pair_2 in pair_list:
                    if sample_pair_2.name != sample_pair.name:
                        if ( sample_pair_2.tag1 == sample_pair.tag1 or \
                            sample_pair_2.tag2 == sample_pair.tag1 or \
                            sample_pair_2.tag1 == sample_pair.tag2 or \
                            sample_pair_2.tag2 == sample_pair.tag2  ):
                            continue
                        sample_name_train_list.append(sample_pair_2.name)
                #for item in sample_names:
                #    if any( i in item.split('_') for i in sample.split('_') ): # any(i in b for i in a)
                #        continue
                #    sample_name_train_list.append(item)
                print_to_log(log_FH, ", ".join(sample_name_train_list)+'\n')
                #----------------------------------------------------------------------
                if args.gnorm:
                    # Count number of up and down for test set
                    test_expr_flag_value_count = df_loso[ df_loso['sample_name'].isin(sample_name_test_list) ]['expr_flag'].value_counts()
                    num_up_test = pd.Series(test_expr_flag_value_count)[1]
                    num_dn_test = pd.Series(test_expr_flag_value_count)[-1]
                    print_to_log(log_FH,'# of up regulated in test set = '+str(num_up_test)+'\n')
                    print_to_log(log_FH,'# of down regulated in test set = '+str(num_dn_test)+'\n\n')
                    # Normalize up and down regulated for train samples
                    train_expr_flag_value_count = df_loso[ df_loso['sample_name'].isin(sample_name_train_list) ]['expr_flag'].value_counts()
                    try:
                        num_up = pd.Series(train_expr_flag_value_count)[1]
                    except KeyError:
                        num_up = 0
                    print_to_log(log_FH,'# of up regulated in train set = '+str(num_up)+'\n')
                    try:
                        num_dn = pd.Series(train_expr_flag_value_count)[-1]
                    except KeyError:
                        num_dn = 0
                    print_to_log(log_FH,'# of down regulated in train set= '+str(num_dn)+'\n\n')
                    if num_up > num_dn:
                        idx = df_loso.index[ (df_loso['expr_flag'] == 1) & df_loso['sample_name'].isin(sample_name_train_list) ]
                        tmp_idx = np.random.choice(idx, size=(num_up-num_dn), replace=False)
                        df_loso.loc[tmp_idx, 'expr_flag'] = 0
                        del idx
                    elif num_dn > num_up:
                        idx = df_loso.index[ (df_loso['expr_flag'] == -1) & df_loso['sample_name'].isin(sample_name_train_list) ]
                        tmp_idx = np.random.choice(idx, size=(num_dn-num_up), replace=False)
                        df_loso.loc[tmp_idx, 'expr_flag'] = 0
                        del idx            
                    # Count again after normalization
                    train_expr_flag_value_count_norm = df_loso[ df_loso['sample_name'].isin(sample_name_train_list) ]['expr_flag'].value_counts()
                    num_up_norm = pd.Series(train_expr_flag_value_count_norm)[1]
                    num_dn_norm = pd.Series(train_expr_flag_value_count_norm)[-1]
                    tot_items_norm = num_up_norm + num_dn_norm
                    print_to_log(log_FH,'# of up regulated after norm = '+str(num_up_norm)+'\n')
                    print_to_log(log_FH,'# of down regulated after norm = '+str(num_dn_norm)+'\n')
                    print_to_log(log_FH,'# of Total items in train set after norm = '+str(tot_items_norm)+'\n\n')
                #---------------------------------------------------------------------------
                ## Run K-Fold on gene level ##
                gene_id_list = np.asarray( list( set( df_loso[ df_loso['sample_name'].isin(sample_name_test_list) ]['gene_id'].tolist() ) ) ) # me-class
        #        gene_id_list = list(set( df_loso[ df_loso['sample_name'].isin(sample_name_test_list) ]['gene_id'].tolist() )) # me-class
        #        gene_id_list = np.asarray( set( df_loso[ ((df_kf.expr_flag==-1) | (df_kf.expr_flag==1)) \
        #                & df_loso['sample_name'].isin(sample_name_train_list) ]['gene_id'].tolist() ) ) # 
        #        gene_id_list = np.asarray( set( df_loso['gene_id'].tolist() ) )
                kf = KFold(n_splits=args.folds, shuffle=shuffle)
                kf.get_n_splits(gene_id_list)
                kf_run_idx = 0    
                for train_index, test_index in kf.split(gene_id_list): # K-Fold gene split
                    kf_run_idx += 1 # for output track
                    df_kf = df_loso.copy()
                    I_train, I_test = gene_id_list[train_index], gene_id_list[test_index]
                    # Mark train set 
                    idx = df_kf.index[ ((df_kf.expr_flag==-1) | (df_kf.expr_flag==1)) & \
                        df_kf['sample_name'].isin(sample_name_train_list) & ~df_kf['gene_id'].isin(I_test) ]
                    df_kf.loc[idx, 'clf_flag'] = 'train'
                    del idx
                    
                    # Mark test set
                    idx = df_kf.index[ ((df_kf.expr_flag==-1) | (df_kf.expr_flag==1)) & \
                        df_kf['sample_name'].isin(sample_name_test_list) & df_kf['gene_id'].isin(I_test) ]
                    df_kf.loc[idx, 'clf_flag'] = 'test'
                    del idx
                
                    # Select traing data set
                    df_kf_train = df_kf[ df_kf['clf_flag']=='train' ].copy(deep=True)
                    df_kf_test = df_kf[ df_kf['clf_flag']=='test' ].copy(deep=True)
                
                    y_train = df_kf_train['expr_flag'].values
                
                    #Setup classifier and run it
                    clf = RandomForestClassifier(n_estimators=args.num_trees, n_jobs=args.threads) # me-class 1001 default
                    #clf.fit(df_kf_train[features], y_train)
                    df_kf_train.drop(drop_col_list, axis=1, inplace=True)
                    df_kf_test.drop(drop_col_list, axis=1, inplace=True)
                    clf.fit(df_kf_train, y_train)
                
                    #Performance Evaluation
                    if len(df_kf_test) > 0:
                        y_test_prob = clf.predict_proba(df_kf_test)
                        y_test_pred = clf.predict(df_kf_test)
                        # store performance values in df_out
                        for i,idx in enumerate(df_kf_test.index):
                            df.loc[idx, 'expr_pred'] = y_test_pred[i]
                            df.loc[idx, 'prob_dn'] = y_test_prob[i][0]
                            df.loc[idx, 'prob_up'] = y_test_prob[i][1]
                
                    # Calculate feature importances
                    if args.featureImportance:
                        df_fi = pd.concat( [df_fi, (pd.DataFrame( [(clf.feature_importances_)],  columns=feature_column_names))] )
                
                    # Clear Memory
                    print_to_log(log_FH,'Memory use: '+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)+'MB\n')
                    df_kf.drop(df_kf.index, inplace=True)
                    df_kf_train.drop(df_kf_train.index, inplace=True)
                    df_kf_test.drop(df_kf_test.index, inplace=True)
                    del df_kf, df_kf_train, df_kf_test, clf
                
                df_loso.drop(df_loso.index, inplace=True)
                del df_loso
                
        # Print final dataframe without feature columns
        df_out = df.copy()
        out_cols=['gene_id-sample_name','gene_id','sample_name','expr_value','expr_flag','prob_dn','prob_up','expr_pred']
        df_out = df_out[ (df_out['expr_pred'].isnull()==False) ]
    #    df_out = df_out[ np.isnan(df_out['expr_pred']) == False ]    
        df_out[out_cols].to_csv(args.tag+'.RandomForestClassifier.csv', sep=',', index=False)
        
        # Feature importance output
        if args.featureImportance:
            df_fi.to_csv(args.tag+'.featureImportance.perFold.csv', sep=',', index=False)
            #dfi_items = [','.join(['feature','sum','mean'])]
            feat_imp_data = [ [x, df_fi[x].sum(), df_fi[x].mean()] for x in df_fi.columns ]
            feat_imp_df = pd.DataFrame(feat_imp_data, columns=['feature','sum','mean'])
            feat_imp_df.to_csv(args.tag + '.featureImportance.mean.csv', index=False)
            if not args.no_plot_featureImportance:
                param_data_dict = dict()
                param_dict = dict()
                for interp_file in args.interp_files:
                    param_data_dict, param_dict = read_params_from_interp_2(interp_file)
                    continue
                plot_featureImportance(feat_imp_df,param_dict,param_data_dict,args)
    return

def exec_classify_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('interp_files', metavar='interp_files', type=str, nargs='+',
        default=argparse.SUPPRESS,
        help='interpolation files for classification')
    parser_required.add_argument('-i', '--input_list', dest='input_list',
        default=argparse.SUPPRESS,
        required=True, help='Input list of sample names and file locations for pairings.')
    parser_classifier = parser.add_argument_group('classifier arguments')
    parser_classifier.add_argument('--num_trees',
        type=int, default=5001, help='Number of trees for Random Forest Classifier')
    parser_classifier.add_argument('-t', '--threads',
        type=int, default=8, help='Number of Processors for RF run')
    parser_classifier.add_argument('--folds',
        type=int, default=10, help='Number of folds for gene level cross-fold validation')
    parser_classifier.add_argument('--tag',
        default='classifier_results', help='Tag which will start all output filenames.')
    parser_classifier.add_argument('--no_shuffle',
        action='store_true', default=False, help='Do no shuffle data during kfold division.')
    parser_classifier.add_argument('--ngnorm',
        action='store_false', dest='gnorm', default=True, help='Normalize gene count or not') 
    parser_featureImportance = parser.add_argument_group('featureImportance arguments')
    parser_featureImportance.add_argument('--featureImportance',
        action='store_true', default=False, help='Compute feature importances') 
    parser_featureImportance.add_argument('--no_plot_featureImportance',
        action='store_true', default=False,
        help='Do not plot feature importances. Unused unless --featureImportance is set.')
    parser_featureImportance.add_argument('--featureImportance_max_y', default=0, type=float,
        help='max y-value for feature importance, set to 0 to auto-scale. Unused unless --featureImportance is set.')
    parser_general = parser.add_argument_group('general arguments')
    parser_general.add_argument('--logfile', action='store', dest='logfile',
        default='classify.log', help='log file')
    #parser._action_groups.reverse()
    return(parser)
