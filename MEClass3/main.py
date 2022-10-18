#
# Authors: Manoj Kumar Singh (manoj@wustl.edu) and John R Edwards (jredwards@wustl.edu).
#
#---------------------------------------------------------------------------

import argparse

from MEClass3.subcommand_parser_functions import setup_subparsers, Subcommand
from MEClass3.preprocess import exec_preprocess, exec_preprocess_help
from MEClass3.interpolation import exec_interp, exec_interp_help
from MEClass3.overlap_genes_enhancers import exec_overlap_genes_enhancers, exec_overlap_genes_enhancers_help
from MEClass3.classifier import exec_classify, exec_classify_help
from MEClass3.evaluation import exec_eval, exec_eval_help
from MEClass3.merge_data import exec_merge_data, exec_merge_data_help
#from MEClass3.filter import exec_filter, exec_filter_help
#from MEClass3.setup_expr import exec_setup_expr, exec_setup_expr_help
#from MEClass3.merge_features import exec_merge_features, exec_merge_features_help
#from MEClass3.process_regions import exec_proc_reg, exec_proc_reg_help

def main():
    meclass3_description = '''Tool for integrative analysis of methylation and expression data'''
    parser = argparse.ArgumentParser(description=meclass3_description)
    subcommands = setup_subcommands()
    parser = setup_subparsers(parser, subcommands)
    args = parser.parse_args()
    args.func(args)
    
def setup_subcommands():       
    subcommands = []
    
    subcommands.append( Subcommand( name="preprocess", 
                                function=exec_preprocess,
                                help_function=exec_preprocess_help,
                                help = "Preprocess methylation data for sample pairs",
                                desc = "Preprocess methylation data for sample pairs."))
    
    subcommands.append( Subcommand( name="overlap_genes_enhancers", 
                                function=exec_overlap_genes_enhancers,
                                help_function=exec_overlap_genes_enhancers_help,
                                help = "Assigns enhancer annotations to genes",
                                desc = "Assigns enhancer annotations to genes"))
    
    subcommands.append( Subcommand( name="interp", 
                                function=exec_interp,
                                help_function=exec_interp_help,
                                help = "Interpolates methylation data to create methylation signatures",
                                desc = "Interpolates methylation data to create methylation signatures for each gene, enhancer, or region."))
    
    subcommands.append( Subcommand( name="merge_data", 
                                function=exec_merge_data,
                                help_function=exec_merge_data_help,
                                help = "Merge methylation signatures and expression data",
                                desc = "Merge interpolated methylation signatures and expression data to create input files for classify step."))

    subcommands.append( Subcommand( name="classify", 
                                function=exec_classify,
                                help_function=exec_classify_help,
                                help = "Train and run classifier",
                                desc = "Trains classifier to predict expression changes using methylation signatures for each gene/enhancer/region"))

    subcommands.append( Subcommand( name="eval", 
                                function=exec_eval,
                                help_function=exec_eval_help,
                                help = "Evaluate classifier performance",
                                desc = "Evaluate classifier performance."))

#    subcommands.append( Subcommand( name="merge_features", 
#                                function=exec_merge_features,
#                                help_function=exec_merge_features_help,
#                                help = "Merge features",
#                                desc = "Merge features description"))

#    subcommands.append( Subcommand( name="setup_expr", 
#                                function=exec_setup_expr,
#                                help_function=exec_setup_expr_help,
#                                help = "setup expr",
#                                desc = "setup expr description"))

#    subcommands.append( Subcommand( name="proc_reg", 
#                                function=exec_proc_reg,
#                                help_function=exec_proc_reg_help,
#                                help = "Process region",
#                                desc = "Process region description"))

#    subcommands.append( Subcommand( name="filter", 
#                                function=exec_filter,
#                                help_function=exec_filter_help,
#                                help = "Filters data",
#                                desc = "Filter description"))


    return subcommands

if __name__ == '__main__':
    main()
