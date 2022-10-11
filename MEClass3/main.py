#
# Authors: Manoj Kumar Singh (manoj@wustl.edu) and John R Edwards (jredwards@wustl.edu).
#
#---------------------------------------------------------------------------

import argparse

from MEClass3.subcommand_parser_functions import setup_subparsers, Subcommand
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
    parser = argparse.ArgumentParser()
    subcommands = setup_subcommands()
    parser = setup_subparsers(parser, subcommands)
    args = parser.parse_args()
    args.func(args)
    
def setup_subcommands():       
    subcommands = []
    
    subcommands.append( Subcommand( name="filter", 
                                function=exec_filter,
                                help_function=exec_filter_help,
                                help = "Filters data",
                                desc = "Filter description"))

    subcommands.append( Subcommand( name="proc_sample", 
                                function=exec_proc_sample,
                                help_function=exec_proc_sample_help,
                                help = "process sample",
                                desc = "process sample description"))
    
    subcommands.append( Subcommand( name="proximity_list", 
                                function=exec_proximity_list,
                                help_function=exec_proximity_list_help,
                                help = "Proximity list",
                                desc = "proximity list description"))

    subcommands.append( Subcommand( name="merge_features", 
                                function=exec_merge_features,
                                help_function=exec_merge_features_help,
                                help = "Merge features",
                                desc = "Merge features description"))

    subcommands.append( Subcommand( name="setup_expr", 
                                function=exec_setup_expr,
                                help_function=exec_setup_expr_help,
                                help = "setup expr",
                                desc = "setup expr description"))

    subcommands.append( Subcommand( name="run_clf", 
                                function=exec_run_clf,
                                help_function=exec_run_clf_help,
                                help = "run classifier",
                                desc = "run classifier description"))

    subcommands.append( Subcommand( name="eval", 
                                function=exec_eval,
                                help_function=exec_eval_help,
                                help = "Evaluation",
                                desc = "Evaluation description"))

    subcommands.append( Subcommand( name="proc_reg", 
                                function=exec_proc_reg,
                                help_function=exec_proc_reg_help,
                                help = "Process region",
                                desc = "Process region description"))

    return subcommands

if __name__ == '__main__':
    main()
