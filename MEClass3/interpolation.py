
import argparse
import resource
import time
from collections import defaultdict

from MEClass3.interpolation_functions import generate_out_header
from MEClass3.interpolation_functions import add_fail
from MEClass3.interpolation_functions import interp_list_sp
from MEClass3.interpolation_functions import interp_list_mp
#from MEClass3.interpolation_functions import generate_param_list
from MEClass3.interpolation_functions import generate_param_comment
from MEClass3.interpolation_functions import format_fail_dict
from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import read_anno_file 
from MEClass3.io_functions import eprint 
from MEClass3.io_functions import read_bed_file
from MEClass3.io_functions import format_args_to_print
from MEClass3.sample import read_sample_file


def exec_interp(args):
    anno_file = args.anno_file
    anno_type = args.anno_type
    out_header = generate_out_header(args.num_interp_points, anno_type, args.data_type)
    pair_list = read_sample_file(args.input_list)
    if args.sample:
        sample_set_flag = False
        for sample_pair in pair_list:
            if args.sample == sample_pair.name:
                pair_list = [ sample_pair ]
                sample_set_flag = True
        if not sample_set_flag:
            eprint("Can't find sample name " + args.sample + " in pair list. Check pair list and sample name and rerun.")
            exit()
    with open(args.logfile, 'w') as log_FH:
        print_to_log(log_FH, format_args_to_print(args))
        anno_list_prefilter = read_anno_file(anno_file, anno_type)
        anno_list_postfilter =[]
        fail_text = 'regions'
        print_to_log(log_FH, '\nSetting up Annotation Information\n')
        anno_fail_dict = defaultdict(list)
        if anno_type == 'tss':
            for gene in anno_list_prefilter:
                if gene.gene_length() < args.min_gene_length:
                    anno_fail_dict = add_fail( 'gene_length', gene.id, anno_fail_dict)
                elif args.filter_cdsStats and ( (gene.cdsStartStat != 'cmpl') or (gene.cdsEndStat != 'cmpl') ):
                    anno_fail_dict = add_fail( 'cdsStats', gene.id, anno_fail_dict)
                else:
                    anno_list_postfilter.append(gene)
            out_file_suffix = '_gene_interp'
            fail_text = 'genes'
            #param_data = generate_param_list(args.num_interp_points, args.ibin_inp, args.data_type, anno_type)
            param_data = generate_param_comment(args.num_interp_points, args.ibin_inp, args.data_type, anno_type)
        elif anno_type == 'enh':
            anno_list_postfilter = anno_list_prefilter
            out_file_suffix = '_enh_interp'
            fail_text = 'enhancers'
            #param_data = generate_param_list(args.num_interp_points, args.refl_inp, args.data_type, anno_type)
            param_data = generate_param_comment(args.num_interp_points, args.refl_inp, args.data_type, anno_type)
        else:
            eprint("Can't recognize anno_type: " + anno_type + "\nCheck --anno_type specification in help.")
            exit()
        print_to_log(log_FH, format_fail_dict(anno_fail_dict, fail_text) + '\n\n')
        for sample_pair in pair_list:
            sample_id = sample_pair.name 
            sample_file = args.output_path + '/' + sample_pair.name +'.bedgraph' 
            out_file = args.output_path + "/" + sample_pair.name + out_file_suffix + ".csv"
            #param_file = args.output_path + "/" + sample_pair.name + out_file_suffix + ".param"
            print_to_log(log_FH, "processing: " + sample_id + " -> " + sample_file + "\n")
            print_to_log(log_FH, "out_file: " + out_file + "\n")
            #print_to_log(log_FH, "param_file: " + param_file + "\n")
            dict_bed = {}
            bed_file_lines = read_bed_file(sample_file)
            chr_track = 'chr00'
            for line in bed_file_lines:
                if line.split()[0] != chr_track:
                    chr_track = line.split()[0]
                    dict_bed[chr_track] = {}
                    dict_bed[chr_track][int(line.split()[1])] = float(line.split()[3])
                    continue
                dict_bed[chr_track][int(line.split()[1])] = float(line.split()[3])
            del bed_file_lines[:]
            print_to_log(log_FH, 'Memory use: '+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)+'MB\n')
            gene_fail_dict = defaultdict(list)
            interp_start_time = time.perf_counter()
            if args.num_proc == 0:
                out_data, gene_fail_dict = interp_list_sp(anno_list_postfilter, dict_bed, sample_id, gene_fail_dict, args)
            else:
                out_data, gene_fail_dict = interp_list_mp(anno_list_postfilter, dict_bed, sample_id, gene_fail_dict, args)
            interp_end_time = time.perf_counter()
            print_to_log(log_FH, f"Interpolation time: {interp_end_time - interp_start_time:0.4f} seconds\n\n")
            print_to_log(log_FH, format_fail_dict(gene_fail_dict, fail_text))
            with open(out_file, 'w') as out_FH:
                out_FH.write(param_data + "\n")
                out_FH.write(out_header)
                out_FH.write(out_data)
            #with open(param_file, 'w') as out_FH:
            #    out_FH.write("\n".join(param_data))
            dict_bed.clear()
            del dict_bed


def exec_interp_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-a','--anno_file', default=argparse.SUPPRESS, required=True, help='region or gene annotation file')
    parser_required.add_argument('-t', dest='tag_inp', default=argparse.SUPPRESS, required=True, help='Tag for output')
    parser_required.add_argument('-i', '--input_list', default=argparse.SUPPRESS, required=True, help='List of sample pairs')
    parser_general = parser.add_argument_group('general arguments')
    parser_general.add_argument('--anno_type', default="tss", 
                                 choices=["tss", "enh"],
                                 help='region or gene annotation file')
    parser_general.add_argument('--data_type', default="mC", 
                                 choices=["mC", "hmC", "other"],
                                 help='type of data')
    parser_general.add_argument('--sample', default=None, help='Use to select specific sample from pair list and only run that one.')
    parser_general.add_argument('--num_proc', type=int, default=0, help='number of processes to run')
    parser_general.add_argument('-o', '--output_path', action='store', dest='output_path',
        default='intermediate_files', help='Path to Output')
    parser_general.add_argument('--logfile', dest='logfile', default='interpolation.log', help='log file')
    parser_filter = parser.add_argument_group('filtering arguments')
    parser_filter.add_argument('--min_gene_cpgs', type=int, default=40, help='Minimum CpGs assayed')
    parser_filter.add_argument('--min_gene_meth', type=float, default=0.2, help='Must have at least 1 CpG with at least this meth diff')
    parser_filter.add_argument('--min_gene_length', type=int, default=5000, help='Min. gene length, ideally shorter than window size')
    parser_filter.add_argument('--filter_cdsStats', action='store_true', default=True, help='Whether to filter on the cdsStart and cdsEnd as cmpl' )
    parser_filter.add_argument('--min_reg_cpgs', type=int, default=2, help='Minimum CpGs assayed in region')
    parser_interp = parser.add_argument_group('interpolation arguments')
    parser_interp.add_argument('--sigma', type=int, default=50, help='Value of sigma for Gaussian smoothing')
    parser_interp.add_argument('--num_interp_points', type=int, default=500, help='Number of interp points')
    parser_interp.add_argument('--flankNorm', action='store_true', default=False, help='Flank Normalized interpolation')
    parser_interp.add_argument('--ibin', dest='ibin_inp', type=int, default=5000, help='Size of bin around TSS/GB')
    parser_interp.add_argument('--anchor_window', dest='anch_win', type=int, default=100000, help='Anchor window')
    parser_interp.add_argument('--rfl', dest='refl_inp', type=int, default=500, help='RE flnk length')
    parser_interp.add_argument('--rff', dest='reff_inp', type=int, default=25, help='RE flnk features')
    parser_interp.add_argument('--frf', action='store_true', dest='fixedReFlnk', default=False, help='Fixed features for RE flank')
    #parser.add_argument('--geneSelect', action='store_true', default=False, help='Interpolation for slected genes?')
    #parser._action_groups.reverse()
    return(parser)



