import sys
import resource 

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score


def exec_run_clf(args):
    # Load data csv file
    df = ( pd.read_csv(args.dfi_inp) ).set_index('gene_id-sample_name')
    df = df[~df.isin([np.nan, np.inf, -np.inf]).any(1)]
    # List of sample names
    sample_names =  list( set( df['sample_name'].tolist() ) )
    #Output info. This will store classifier performance
    df['prob_dw'] = 0.0
    df['prob_up'] = 0.0
    df['expr_pred'] = np.nan
    # Add train and test flag for classifier to use and mark all null
    df['clf_flag'] = 'null'
    # Feature importance
    if args.fsl_inp == 1: # TSS 
        feature_column_names = ( df.columns[pd.Series(df.columns).str.startswith('ftss')] ).tolist()
    if args.fsl_inp == 2: # TSS + RE
        feature_column_names = ( df.columns[pd.Series(df.columns).str.startswith('f')] ).tolist()
    if args.fsl_inp == 3: #RE Only
        feature_column_names = ( df.columns[pd.Series(df.columns).str.startswith('fre')] ).tolist()
    df_fi = pd.DataFrame(columns=feature_column_names)
    #===============================================
    if args.ss:
        if args.gnorm:
            # Normalize up and down
            expr_flag_value_count = df['expr_flag'].value_counts()
            num_up = pd.Series(expr_flag_value_count)[1]
            num_dn = pd.Series(expr_flag_value_count)[-1]
            if num_up > num_dn:
                idx = df.index[ (df['expr_flag'] == 1) ]
                df.at[np.random.choice(idx, size=(num_up-num_dn), replace=False), 'expr_flag'] = 0
                del idx
            elif num_dn > num_up:
                idx = df.index[ (df['expr_flag'] == -1) ]
                df.at[np.random.choice(idx, size=(num_dn-num_up), replace=False), 'expr_flag'] = 0
                del idx
  
        gene_id_list = np.asarray( list(set( df['gene_id'].tolist() ) ))
        kf = KFold(n_splits=10 ) # 10 fold #shuffle=True in me-class
        kf.get_n_splits(gene_id_list)
        for train_index, test_index in kf.split(gene_id_list): # K-Fold gene split
            I_train, I_test = gene_id_list[train_index], gene_id_list[test_index]
            df_kf = df.copy()
            idx = df_kf.index[ ~df_kf['gene_id'].isin(I_test) ]
            df_kf.at[idx, 'clf_flag'] = 'train'
            del idx
            idx = df_kf.index[ df_kf['gene_id'].isin(I_test) ]
            df_kf.at[idx, 'clf_flag'] = 'test'
            del idx
            # Select traing data set
            df_kf_train = df_kf[ df_kf['clf_flag']=='train' ]
            df_kf_test = df_kf[ df_kf['clf_flag']=='test' ]
            features = df_kf_train.columns[pd.Series(df_kf_train.columns).str.startswith('f')]
            y_train = df_kf_train['expr_flag'].values
            #Setup classifier and run it
            clf = RandomForestClassifier(n_estimators=args.ntr_inp, n_jobs=args.npr_inp) #
            clf.fit(df_kf_train[features], y_train)
            y_test_prob = clf.predict_proba(df_kf_test[features])
            y_test_pred = clf.predict(df_kf_test[features])
            # store performance values in df_out
            for i,idx in enumerate(df_kf_test[features].index):
                    df.at[idx, 'expr_pred'] = y_test_pred[i]
                    df.at[idx, 'prob_dw'] = y_test_prob[i][0]
                    df.at[idx, 'prob_up'] = y_test_prob[i][1]
            # Feature importance
            df_fi = pd.concat( [df_fi, (pd.DataFrame( [(clf.feature_importances_)],  columns=feature_column_names))] )
            df_kf.drop(df_kf.index, inplace=True)
            df_kf_train.drop(df_kf_train.index, inplace=True)
            df_kf_test.drop(df_kf_test.index, inplace=True)
            del df_kf, df_kf_train, df_kf_test, clf
    #===============================================
    if not args.ss: 
        # loso loop
        for sample in sample_names:
            df_loso = df.copy()
            # test sample
            sample_name_test_list = []
            sample_name_test_list.append(sample)
            sys.stdout.write(", ".join(sample_name_test_list)+'\n')
            # train samples
            sample_name_train_list = []
            for item in sample_names:
                if any( i in item.split('_') for i in sample.split('_') ): # any(i in b for i in a)
                    continue
                sample_name_train_list.append(item)
            sys.stdout.write(", ".join(sample_name_train_list)+'\n')
            #----------------------------------------------------------------------
            if args.gnorm:
                # Count number of up and down for test set
                test_expr_flag_value_count = df_loso[ df_loso['sample_name'].isin(sample_name_test_list) ]['expr_flag'].value_counts()
                num_up_test = pd.Series(test_expr_flag_value_count)[1]
                num_dn_test = pd.Series(test_expr_flag_value_count)[-1]
                sys.stdout.write('# of up regulated in test set = '+str(num_up_test)+'\n')
                sys.stdout.write('# of down regulated in test set = '+str(num_dn_test)+'\n\n')
                # Normalize up and down regulated for train samples
                train_expr_flag_value_count = df_loso[ df_loso['sample_name'].isin(sample_name_train_list) ]['expr_flag'].value_counts()
                try:
                    num_up = pd.Series(train_expr_flag_value_count)[1]
                except KeyError:
                    num_up = 0
                sys.stdout.write('# of up regulated in train set = '+str(num_up)+'\n')
                try:
                    num_dn = pd.Series(train_expr_flag_value_count)[-1]
                except KeyError:
                    num_dn = 0
                sys.stdout.write('# of down regulated in train set= '+str(num_dn)+'\n\n')
                if num_up > num_dn:
                    idx = df_loso.index[ (df_loso['expr_flag'] == 1) & df_loso['sample_name'].isin(sample_name_train_list) ]
                    df_loso.at[np.random.choice(idx, size=(num_up-num_dn), replace=False), 'expr_flag'] = 0
                #    df_loso.set_value(np.random.choice(idx, size=(num_up-num_dn), replace=False), 'expr_flag', 0)
                    del idx
                elif num_dn > num_up:
                    idx = df_loso.index[ (df_loso['expr_flag'] == -1) & df_loso['sample_name'].isin(sample_name_train_list) ]
                    df_loso.at[np.random.choice(idx, size=(num_dn-num_up), replace=False), 'expr_flag'] = 0
                #    df_loso.set_value(np.random.choice(idx, size=(num_dn-num_up), replace=False), 'expr_flag', 0)
                    del idx            
                # Count again after normalization
                train_expr_flag_value_count_norm = df_loso[ df_loso['sample_name'].isin(sample_name_train_list) ]['expr_flag'].value_counts()
                num_up_norm = pd.Series(train_expr_flag_value_count_norm)[1]
                num_dn_norm = pd.Series(train_expr_flag_value_count_norm)[-1]
                tot_items_norm = num_up_norm + num_dn_norm
                sys.stdout.write('# of up regulated after norm = '+str(num_up_norm)+'\n')
                sys.stdout.write('# of down regulated after norm = '+str(num_dn_norm)+'\n')
                sys.stdout.write('# of Total items in train set after norm = '+str(tot_items_norm)+'\n\n')
                #
            #---------------------------------------------------------------------------
            ## Run K-Fold on gene level ##
            gene_id_list = np.asarray( list( set( df_loso[ df_loso['sample_name'].isin(sample_name_test_list) ]['gene_id'].tolist() ) ) ) # me-class
    #        gene_id_list = list(set( df_loso[ df_loso['sample_name'].isin(sample_name_test_list) ]['gene_id'].tolist() )) # me-class
    #        gene_id_list = np.asarray( set( df_loso[ ((df_kf.expr_flag==-1) | (df_kf.expr_flag==1)) \
    #                & df_loso['sample_name'].isin(sample_name_train_list) ]['gene_id'].tolist() ) ) # 
    #        gene_id_list = np.asarray( set( df_loso['gene_id'].tolist() ) )
            kf = KFold(n_splits=10, shuffle=args.suf_inp) # 10 fold #shuffle=True in me-class
            kf.get_n_splits(gene_id_list)
            kf_run_idx = 0    
            for train_index, test_index in kf.split(gene_id_list): # K-Fold gene split
                kf_run_idx += 1 # for output track
                df_kf = df_loso.copy()
                I_train, I_test = gene_id_list[train_index], gene_id_list[test_index]
                # Mark train set 
                idx = df_kf.index[ ((df_kf.expr_flag==-1) | (df_kf.expr_flag==1)) & \
                        df_kf['sample_name'].isin(sample_name_train_list) & ~df_kf['gene_id'].isin(I_test) ]
                df_kf.at[idx, 'clf_flag'] = 'train'
                #df_kf.set_value(idx, 'clf_flag', 'train')
                del idx
                    
                # Mark test set
                idx = df_kf.index[ ((df_kf.expr_flag==-1) | (df_kf.expr_flag==1)) & \
                        df_kf['sample_name'].isin(sample_name_test_list) & df_kf['gene_id'].isin(I_test) ]
                df_kf.at[idx, 'clf_flag'] = 'test'
                #df_kf.set_value(idx, 'clf_flag', 'test')
                del idx
                
                # Select traing data set
                df_kf_train = df_kf[ df_kf['clf_flag']=='train' ]
                df_kf_test = df_kf[ df_kf['clf_flag']=='test' ]
                
                # Select features
                if args.fsl_inp == 1: # TSS
                    features = df_kf_train.columns[pd.Series(df_kf_train.columns).str.startswith('ftss')]
                if args.fsl_inp == 2: # TSS + RE
                    features = df_kf_train.columns[pd.Series(df_kf_train.columns).str.startswith('f')]
                if args.fsl_inp == 3: # RE
                    features = df_kf_train.columns[pd.Series(df_kf_train.columns).str.startswith('fre')]
                
                y_train = df_kf_train['expr_flag'].values
                
                #Setup classifier and run it
                clf = RandomForestClassifier(n_estimators=args.ntr_inp, n_jobs=args.npr_inp) # me-class 1001 default
                clf.fit(df_kf_train[features], y_train)
                
                #Performance Evaluation
                y_test_prob = clf.predict_proba(df_kf_test[features])
                y_test_pred = clf.predict(df_kf_test[features])
                # store performance values in df_out
                for i,idx in enumerate(df_kf_test[features].index):
                    df.at[idx, 'expr_pred'] = y_test_pred[i]
                    df.at[idx, 'prob_dw'] = y_test_prob[i][0]
                    df.at[idx, 'prob_up'] = y_test_prob[i][1]

                    # df.set_value(idx, 'expr_pred', y_test_pred[i]) 
                    # df.set_value(idx, 'prob_dw', y_test_prob[i][0]) 
                    # df.set_value(idx, 'prob_up', y_test_prob[i][1]) 
                    
                
                # Feature importance.
                #print(clf.feature_importances_)
                df_fi = pd.concat( [df_fi, (pd.DataFrame( [(clf.feature_importances_)],  columns=feature_column_names))] )
                
                sys.stdout.write('Memory use: '+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)+'MB\n')
                
                # Clear Memory
                
                df_kf.drop(df_kf.index, inplace=True)
                df_kf_train.drop(df_kf_train.index, inplace=True)
                df_kf_test.drop(df_kf_test.index, inplace=True)
                
                del df_kf, df_kf_train, df_kf_test, clf
                
            df_loso.drop(df_loso.index, inplace=True)
            del df_loso
    #============================
    # Final dataframe without feature column
    df_out = df.copy()
    df_out = df_out[ (df_out['expr_pred'].isnull()==False) ]
#    df_out = df_out[ np.isnan(df_out['expr_pred']) == False ]    
    df_out = df_out.drop(df_out.columns[pd.Series(df_out.columns).str.startswith('f')] , axis=1)
    df_out.to_csv(args.tag_inp+'.RandomForestClassifier.csv', sep=',')
    #
    # Feature importance output
    df_fi.to_csv(args.tag_inp+'.RF_fi.csv', sep=',')
    dfi_items = []
    for item in df_fi.columns:
        dfi_items.append( str(df_fi[item].sum())+'\t'+str(df_fi[item].mean())+'\n' )
    #
    open(args.tag_inp+'_fi_out.txt', 'w').write(''.join(dfi_items))
 
 
def exec_run_clf_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-dfi', action='store', dest='dfi_inp', required=True, help='Dataframe output from interpolation step')
    parser.add_argument('-ntr', action='store', dest='ntr_inp', type=int, default=5001, help='Number of trees for Random Forest Classifier')
    parser.add_argument('-npr', action='store', dest='npr_inp', type=int, default=8, help='Number of Processors for RF run')
    parser.add_argument('-tag', action='store', dest='tag_inp', default='test', help='Tag for Output Writing')
    parser.add_argument('-fsl', action='store', dest='fsl_inp', type=int, default=1, help='Feature Selection. 1: TSS; 2: TSS+RE')
    parser.add_argument('-suf', action='store', dest='suf_inp', type=bool, default=True, help='Shuffle true ot false')
    parser.add_argument('-ss', action='store_true', dest='ss', default=False, help='Single sample or not') 
    parser.add_argument('-ngnorm', action='store_false', dest='gnorm', default=True, help='Normalize gene count or not') 
    #parser._action_groups.reverse()
    return(parser)
