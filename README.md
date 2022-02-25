# ME-Class3

## Requirements
 * Python Packages
     *  numpy
     *  pysam
     *  biopython

## Installation
Clone and install with pip:

````
    git clone https://github.com/edwardslab-wustl/ME-Class3.git
    cd ME-Class3
    pip install .
````

## Usage

### Input files

### Example

## License information


## Parameters

```

1. Process samples: 
	This module processes fractional methylaton files and generates differentail methylation files. 

    -ctp	: Cell type Fofn (Please refer to example data). 
    -pto	: Path to Output.
    -expr	: Name of expression file. Samples without expression information will be omitted.

2. Generate proximity list:
	This module generates   

    sp_proximity_list = subparsers.add_parser('proximity_list', help='proximity list argument')
    sp_proximity_list.add_argument('-pst', action='store', dest='pst_inp', type=int, default=1, help='pend start index')
    sp_proximity_list.add_argument('-ped', action='store', dest='ped_inp', type=int, default=2, help='pend end index')
    sp_proximity_list.add_argument('-nst', action='store', dest='nst_inp', type=int, default=1, help='nend start index')
    sp_proximity_list.add_argument('-ned', action='store', dest='ned_inp', type=int, default=2, help='nend end index')
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

    sp_run_clf = subparsers.add_parser('run_clf', help='run clf argument')
    sp_run_clf.add_argument('-dfi', action='store', dest='dfi_inp', help='Dataframe output from interpolation step')
    sp_run_clf.add_argument('-ntr', action='store', dest='ntr_inp', type=int, default=5001, help='Number of trees for Random Forest Classifier')
    sp_run_clf.add_argument('-npr', action='store', dest='npr_inp', type=int, default=8, help='Number of Processors for RF run')
    sp_run_clf.add_argument('-tag', action='store', dest='tag_inp', default='test', help='Tag for Output Writing')
    sp_run_clf.add_argument('-fsl', action='store', dest='fsl_inp', type=int, default=1, help='Feature Selection. 1: TSS; 2: TSS+RE')
    sp_run_clf.add_argument('-suf', action='store', dest='suf_inp', type=bool, default=True, help='Shuffle true ot false')
