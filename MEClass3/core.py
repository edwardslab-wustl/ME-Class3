#
# Authors: Manoj Kumar Singh (manoj@wustl.edu) and John R Edwards (jredwards@wustl.edu).
#
#---------------------------------------------------------------------------

import argparse

from MEClass3.process_sample import exec_proc_sample
from MEClass3.filter import exec_filter
from MEClass3.setup_expr import exec_setup_expr
from MEClass3.merge_features import exec_merge_features
from MEClass3.interpolation import exec_interp
from MEClass3.proximity_list import exec_proximity_list
from MEClass3.process_regions import exec_proc_reg
from MEClass3.classifier import exec_run_clf
from MEClass3.evaluation import exec_eval
from MEClass3.io_functions import print_tup

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
    subparsers = parser.add_subparsers(dest='command')

    sp_filter = subparsers.add_parser('filter', help='Filter argument')
    sp_filter.add_argument('-expr', action='store', dest='expr_inp', help='Name of first file')
    sp_filter.add_argument('-annot', action='store', dest='annot_inp', help='Name of second file')
    sp_filter.add_argument('-oex', action='store', dest='out_expr_inp', help='Name of output expr file')
    sp_filter.add_argument('-oan', action='store', dest='out_annot_inp', help='Name of output annot file')
    sp_filter.add_argument('-apt', action='store', dest='annot_type', default='RefSeq', help='Gencode/RefSeq')

    sp_proc_sample = subparsers.add_parser('proc_sample', help='Proc Sample Arguments')
    sp_proc_sample.add_argument('-ctp', action='store', dest='ctp_inp', help='Cell type Fofn')
    sp_proc_sample.add_argument('-pto', action='store', dest='pto_inp', default='.', help='Path to Output')
    sp_proc_sample.add_argument('-expr', action='store', dest='expr_inp', help='Name of expression file')

    sp_proximity_list = subparsers.add_parser('proximity_list', help='proximity list argument')
    sp_proximity_list.add_argument('-pst', action='store', dest='pst_inp', type=int, default=1, help='pend start index')
    sp_proximity_list.add_argument('-ped', action='store', dest='ped_inp', type=int, default=2, help='pend end index')
    sp_proximity_list.add_argument('-nst', action='store', dest='nst_inp', type=int, default=1, help='nend start index')
    sp_proximity_list.add_argument('-ned', action='store', dest='ned_inp', type=int, default=2, help='nend end index')
    sp_proximity_list.add_argument('-twn', action='store', dest='twin_inp', type=int, default=5000, help='nend end index')
    sp_proximity_list.add_argument('-glf', action='store', dest='glf_inp', help='gene list file')

    sp_merge_features = subparsers.add_parser('merge_features', help='merge features argument')
    sp_merge_features.add_argument('-pst', action='store', dest='pst_inp', type=int, default=1, help='pend start index')
    sp_merge_features.add_argument('-ped', action='store', dest='ped_inp', type=int, default=2, help='pend end index')
    sp_merge_features.add_argument('-nst', action='store', dest='nst_inp', type=int, default=1, help='nend start index')
    sp_merge_features.add_argument('-ned', action='store', dest='ned_inp', type=int, default=2, help='nend end index')
    sp_merge_features.add_argument('-tif', action='store', dest='tif_inp', help='tss interpolation file')

    sp_interp = subparsers.add_parser('interp', help='interp argument')
    sp_interp.add_argument('-sigma', action='store', dest='sigma_inp', type=int, default=50, help='Value of sigma')
    sp_interp.add_argument('-nip', action='store', dest='nip_inp', type=int, default=500, help='Number of interp points')
    sp_interp.add_argument('-ibin', action='store', dest='ibin_inp', type=int, default=5000, help='Size of bin around TSS/GB')
    sp_interp.add_argument('-ach', action='store', dest='anch_win', type=int, default=100000, help='Anchor window')
    sp_interp.add_argument('-rfl', action='store', dest='refl_inp', type=int, default=500, help='RE flnk length')
    sp_interp.add_argument('-rff', action='store', dest='reff_inp', type=int, default=25, help='RE flnk features')
    sp_interp.add_argument('-rfn', action='store', dest='rfn_inp', help='Region fofn')
    sp_interp.add_argument('-tag', action='store', dest='tag_inp', help='Tag for output')
    sp_interp.add_argument('-sfn', action='store', dest='sfn_inp', help='Sample file name')
    sp_interp.add_argument('-fln', action='store_true', dest='flankNorm', default=False, help='Flank Normalized interpolation')
    sp_interp.add_argument('-gsl', action='store_true', dest='geneSelect', default=False, help='Interpolation  for slected genes?')
    sp_interp.add_argument('-frf', action='store_true', dest='fixedReFlnk', default=False, help='Fixed features for RE flank')
    sp_interp.add_argument('-mmg', action='store', dest='min_gene_meth', type=int, default=40, help='Minimum CpGs assayed')
    sp_interp.add_argument('-mmr', action='store', dest='min_re_meth', type=int, default=2, help='Minimum CpGs assayed')

    sp_setup_expr = subparsers.add_parser('setup_expr', help='setup expr argument')
    sp_setup_expr.add_argument('-intrp', action='store', dest='intrp_inp', help='Interpolation CSV file')
    sp_setup_expr.add_argument('-expr', action='store', dest='expr_inp', help='Expression file')
    sp_setup_expr.add_argument('-fef', action='store', dest='floor_expr', type=bool, default=True, help='Floor expression value?')
    sp_setup_expr.add_argument('-efv', action='store', dest='efv_inp', type=float, default=5.0, help='Expression floor value')
    sp_setup_expr.add_argument('-def', action='store', dest='dexpr_flag', type=int, default=1, help='Method for differential expression')
    sp_setup_expr.add_argument('--expr_cutoff', action='store', dest='expr_cutoff', type=int, default=2, help='fold change in expression cutoff, must use -def 3')

    sp_run_clf = subparsers.add_parser('run_clf', help='Run classification argument')
    sp_run_clf.add_argument('-dfi', action='store', dest='dfi_inp', help='Dataframe output from interpolation step')
    sp_run_clf.add_argument('-ntr', action='store', dest='ntr_inp', type=int, default=5001, help='Number of trees for Random Forest Classifier')
    sp_run_clf.add_argument('-npr', action='store', dest='npr_inp', type=int, default=8, help='Number of Processors for RF run')
    sp_run_clf.add_argument('-tag', action='store', dest='tag_inp', default='test', help='Tag for Output Writing')
    sp_run_clf.add_argument('-fsl', action='store', dest='fsl_inp', type=int, default=1, help='Feature Selection. 1: TSS; 2: TSS+RE')
    sp_run_clf.add_argument('-suf', action='store', dest='suf_inp', type=bool, default=True, help='Shuffle true ot false')
    sp_run_clf.add_argument('-ss', action='store_true', dest='ss', default=False, help='Single sample or not') 
    sp_run_clf.add_argument('-ngnorm', action='store_false', dest='gnorm', default=True, help='Normalize gene count or not') 

    sp_eval = subparsers.add_parser('eval', help='Evaluation argument')
    sp_eval.add_argument('-dfi', action='store', dest='dfi_inp', help='Dataframe output from classification')
    sp_eval.add_argument('-nstp', action='store', dest='steps_inp', type=int, default=101, help='Number of steps')

    sp_proc_reg = subparsers.add_parser('proc_reg', help='Process RE regions')
    sp_proc_reg.add_argument('-gpj', action='store', dest='gapj_inp', type=int, default=500, help='RE closer than this will be joined together.')
    sp_proc_reg.add_argument('-msz', action='store', dest='max_inp', type=int, default=1000, help='Maximum size of RE to take into account.')
    sp_proc_reg.add_argument('-twn', action='store', dest='twin_inp', type=int, default=2000, help='RE falling +/- is this region around TSS will be ignored.')
    sp_proc_reg.add_argument('-glf', action='store', dest='glf_inp', help='Gene information file')
 
    args = parser.parse_args()
    funct = FUNCTION_MAP[args.command]
    funct(args)

if __name__ == '__main__':
    main()
