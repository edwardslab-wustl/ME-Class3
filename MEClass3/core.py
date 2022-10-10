#
# Authors: Manoj Kumar Singh (manoj@wustl.edu) and John R Edwards (jredwards@wustl.edu).
#
#---------------------------------------------------------------------------

import argparse

from MEClass3.process_sample import exec_proc_sample, exec_proc_sample_help
from MEClass3.filter import exec_filter, exec_filter_help
from MEClass3.setup_expr import exec_setup_expr, exec_setup_expr_help
from MEClass3.merge_features import exec_merge_features, exec_merge_features_help
from MEClass3.interpolation import exec_interp, exec_interp_help
from MEClass3.proximity_list import exec_proximity_list, exec_proximity_list_help
from MEClass3.process_regions import exec_proc_reg, exec_proc_reg_help
from MEClass3.classifier import exec_run_clf, exec_run_clf_help
from MEClass3.evaluation import exec_eval, exec_eval_help

def main():
    FUNCTION_MAP = {
            'proc_sample'       : exec_proc_sample,
            'filter'            : exec_filter,
            'proc_reg'          : exec_proc_reg,
            'proximity_list'    : exec_proximity_list,
            'interp'            : exec_interp,
            'merge_features'    : exec_merge_features,
            'setup_expr'        : exec_setup_expr,
            'run_clf'           : exec_run_clf,
            'eval'              : exec_eval,
            }

    parser = argparse.ArgumentParser()
    parse = setup_subparsers(parser)
    args = parser.parse_args()
    funct = FUNCTION_MAP[args.command]
    funct(args)

def setup_subparsers(parser):
    subparsers = parser.add_subparsers(dest='command')
    subparsers.required = True

    filter_desc=str("filter")
    sp_filter = subparsers.add_parser('filter', 
        description=filter_desc, help=filter_desc)
    sp_filter = exec_filter_help(sp_filter)
    
    proc_sample_desc=str("proc_sample")
    sp_proc_sample = subparsers.add_parser('proc_sample', 
        description=proc_sample_desc, help=proc_sample_desc)
    sp_proc_sample = exec_proc_sample_help(sp_proc_sample)
    
    proximity_list_desc=str("proximity_list")
    sp_proximity_list = subparsers.add_parser('proximity_list', 
        description=proximity_list_desc, help=proximity_list_desc)
    sp_proximity_list = exec_proximity_list_help(sp_proximity_list)
    
    merge_features_desc=str("merge_features")
    sp_merge_features = subparsers.add_parser('merge_features', 
        description=merge_features_desc, help=merge_features_desc)
    sp_merge_features = exec_merge_features_help(sp_merge_features)
    
    interp_desc=str("interp")
    sp_interp = subparsers.add_parser('interp', 
        description=interp_desc, help=interp_desc)
    sp_interp = exec_interp_help(sp_interp)
    
    setup_expr_desc=str("setup_expr")
    sp_setup_expr = subparsers.add_parser('setup_expr', 
        description=setup_expr_desc, help=setup_expr_desc)
    sp_setup_expr = exec_setup_expr_help(sp_setup_expr)
    
    run_clf_desc=str("run_clf")
    sp_run_clf = subparsers.add_parser('run_clf', 
        description=run_clf_desc, help=run_clf_desc)
    sp_run_clf = exec_run_clf_help(sp_run_clf)
    
    eval_desc=str("eval")
    sp_eval = subparsers.add_parser('eval', 
        description=eval_desc, help=eval_desc)
    sp_eval = exec_eval_help(sp_eval)
  
    proc_reg_desc=str("proc_reg")
    sp_proc_reg = subparsers.add_parser('proc_reg', 
        description=proc_reg_desc, help=proc_reg_desc)
    sp_proc_reg = exec_proc_reg_help(sp_proc_reg)
    
    return(parser)
 

if __name__ == '__main__':
    main()
