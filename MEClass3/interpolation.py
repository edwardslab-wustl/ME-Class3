from collections import defaultdict
import sys
import re
import resource
import gzip
import time
import argparse

import multiprocessing as mp
from itertools import repeat

from MEClass3.Interpolation_class import Interpolation
from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import read_anno_file 
from MEClass3.io_functions import eprint 
from MEClass3.io_functions import read_bed_file
from MEClass3.io_functions import format_args_to_print
from MEClass3.sample import read_sample_file


def generate_out_header(num_pts, anno_type):
    header = ''
    if anno_type == 'gene_tss':
        header = 'gene_id-sample_name'
        for pos in range(num_pts):
            header = header + (',ftss_'+str(pos))
    elif anno_type == 'enh':
        header = 'enh_loc-gene_id-sample_name'
        for pos in range(num_pts):
            #header = header + (',f'+''.join(file_id.split('_'))+'_'+str(pos)) 
            header = header + (',fre_'+str(pos)) 
    else:
        eprint("Can't recognize anno_type. Check --anno_type specification in help.")
        exit()
    return header + "\n"

def add_fail (status, gene, dict):
    if status in dict:
        dict[status].append(gene)
    else:
        dict[status] = [ gene ]
    return dict

def exec_interp(args):
    anno_file = args.anno_file
    anno_type = args.anno_type
    out_header = generate_out_header(args.num_interp_points, anno_type)
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
    anchor_window = args.anch_win
    with open(args.logfile, 'w') as log_FH:
        print_to_log(log_FH, format_args_to_print(args))
        anno_list_prefilter = read_anno_file(anno_file, anno_type)
        anno_list_postfilter =[]
        dict_sample = {}
        for sample_pair in pair_list:
            sample_id = sample_pair.name 
            sample_file = args.output_path + '/' + sample_pair.name +'.bedgraph' 
            print_to_log(log_FH, "processing: " + sample_id + " -> " + sample_file + "\n")
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
            fail_text = 'regions'
            if anno_type == 'gene_tss':
                for gene in anno_list_prefilter:
                    if gene.gene_length() < args.min_gene_length:
                        gene_fail_dict = add_fail( 'gene_length', gene.id, gene_fail_dict)
                    elif args.filter_cdsStats and ( (gene.cdsStartStat != 'cmpl') or (gene.cdsEndStat != 'cmpl') ):
                        gene_fail_dict = add_fail( 'cdsStats', gene.id, gene_fail_dict)
                    else:
                        anno_list_postfilter.append(gene)
                out_file = args.output_path + "/" + sample_pair.name + '_gene_interp.csv'
                fail_text = 'genes'
            elif anno_type == 'enh':
                anno_list_postfilter = anno_list_prefilter
                out_file = args.output_path + "/" + sample_pair.name + '_enh_interp.csv'
                fail_text = 'enhancers'
            else:
                eprint("Can't recognize anno_type. Check --anno_type specification in help.")
                exit()
            interp_start_time = time.perf_counter()
            if args.num_proc == 0:
                out_data, gene_fail_dict = interp_list_sp(anno_list_postfilter, dict_bed, sample_id, gene_fail_dict, log_FH, args)
            else:
                out_data, gene_fail_dict = interp_list_mp(anno_list_postfilter, dict_bed, sample_id, gene_fail_dict, log_FH, args)
            interp_end_time = time.perf_counter()
            print_to_log(log_FH, f"Interpolation time: {interp_end_time - interp_start_time:0.4f} seconds\n")
            out_data = out_header + out_data
            for status in gene_fail_dict:
                num_fail = len(gene_fail_dict[status])
                print_to_log(log_FH, str(num_fail) + " " + fail_text + " failed to pass " + status + " filter.\n")
                print_to_log(log_FH, ",".join(gene_fail_dict[status]) + "\n")
            with open(out_file, 'w') as out_FH:
                out_FH.write(out_data)
            dict_bed.clear()
            del dict_bed
 
def interp_list_sp(item_list, dict_bed, sample_id, fail_dict, log_FH, args):
    final_results = []
    #fail_dict = defaultdict(list)
    for item in item_list:
        if args.anno_type == 'gene_tss':
            (status, result) = interp_gene(item, dict_bed, sample_id, args)
        elif args.anno_type == 'enh':
            (status, result) = interp_region(item, dict_bed, sample_id, args)
        else:
            eprint("Can't recognize anno_type. Check --anno_type specification in help.")
            exit()
        if status == 'pass':
            final_results.append(result)
        else:
            fail_dict = add_fail(status, result, fail_dict)
    return "\n".join(final_results), fail_dict
    
def interp_list_mp(item_list, dict_bed, sample_id, fail_dict, log_FH, args):
    arg_iterable = zip(item_list, repeat(dict_bed), repeat(sample_id), repeat(args))
    with mp.Pool(processes=args.num_proc) as pool:
        if args.anno_type == 'gene_tss':
            results = pool.starmap(interp_gene, arg_iterable)
        elif args.anno_type == 'enh':
            results = pool.starmap(interp_region, arg_iterable)
        else:
            eprint("Can't recognize anno_type. Check --anno_type specification in help.")
            exit()
        #res = pool.starmap_async(interp_gene, arg_iterable)
        #results = res.get()
    final_results = []
    for (status, result) in results:
        if status == 'pass':
            final_results.append(result)
        else:
            fail_dict = add_fail(status, result, fail_dict)
    return "\n".join(final_results), fail_dict
    
def interp_gene(gene, dict_bed, sample_id, args):
    result = 'na'
    interp_bin = args.ibin_inp
    reg_type='gene'
    cpos_dat, cpos_tmp, dmet_dat, dmet_tmp = [], [], [], []
    tss = gene.tss()
    tes = gene.tes()
    dict_cpg = dict_bed[gene.chr]
    # ----------------------------------------------------------------------------
    # Methylation based gene filters
    # 3. Genes with <40 CpGs assayed within +/-5kb of the TSS
    # 4. Genes with all CpGs within +/-5 kb of the TSS had <0.2 methylation change
    #-----------------------------------------------------------------------------    
    for cpos in range( (tss-interp_bin), (tss+interp_bin)+1 ):            
        if cpos in dict_cpg:
            dmet_tmp.append(dict_cpg[cpos])
            cpos_tmp.append(cpos)
    filter_flag = False
    if len(dmet_tmp) < args.min_gene_cpgs:
        result = ("cpg", gene.id)
        filter_flag = True
    elif max(dmet_tmp) < args.min_gene_meth:
        filter_flag = True
        result = ("meth_diff", gene.id)
    else:
        if not args.flankNorm: # Endpoint correction will be used in this case
            anchor_window = 0
    #    if args.flankNorm:
        for cpos in range( gene.txStart-(interp_bin+anchor_window), gene.txEnd+(interp_bin+anchor_window)+1 ):
            if cpos in dict_cpg:
                dmet_dat.append(dict_cpg[cpos])
                cpos_dat.append(cpos)
        # Missing anchor            
        if args.flankNorm and \
            (( cpos_dat[0] not in range( gene.txStart-(interp_bin+anchor_window), gene.txStart-interp_bin ) ) or \
            ( cpos_dat[-1] not in range( gene.txEnd+interp_bin+1, gene.txEnd+(interp_bin+anchor_window)+1 ) )):
            filter_flag = True
            result = ("flankNorm", gene.id)
    if not filter_flag:
        # Gather interpolated data
        #eprint(gene.id)
        interpolated_dmet_data = Interpolation(cpos_dat, dmet_dat, gene.txStart, gene.txEnd, gene.strand, reg_type, args).dat_proc()
        # cross-check number of interpolated features
        if len(interpolated_dmet_data) != args.num_interp_points:
            #eprint('Inconsistent number of interpolation features: ' + gene.id)
            result = ('interpolation inconsistency in number of features', gene.id)
            #exit() # Exit if number is inconsistent 
        # Write data
        else:
            result = ('pass', gene.id+'-'+sample_id+','+','.join(str(item) for item in interpolated_dmet_data) )
    del dmet_dat[:], cpos_dat[:], dmet_tmp[:], cpos_tmp[:]
    #eprint(result)
    #print_to_log(log_FH, result + "\n")
    return result

def interp_region(region, dict_bed, sample_id, args):
    result = 'na'
    reg_type = 're'
    interp_bin = 500 # for model gene plots
    dict_cpg = dict_bed[region.chr]
    # Regulatory Elements
    cpos_dat, cpos_tmp, dmet_dat, dmet_tmp = [], [], [], []
    # ----------------------------------------------------------------------------
    # Methylation based filters
    # 3. Genes with <2 CpGs assayed within enhancer
    #-----------------------------------------------------------------------------
    for cpos in range( (region.start-interp_bin), (region.end+interp_bin)+1 ):
        if cpos in dict_cpg:
            dmet_tmp.append(dict_cpg[cpos])
            cpos_tmp.append(cpos)
    filter_flag = False
    if len(dmet_tmp) < args.min_reg_cpgs:
        result = ("cpg", region.id)
        filter_flag = True
    else:
        if not args.flankNorm:
            anchor_window = 0
        for cpos in range( region.start-(interp_bin+anchor_window), region.end+(interp_bin+anchor_window)+1 ):
            if cpos in dict_cpg:
                dmet_dat.append(dict_cpg[cpos])
                cpos_dat.append(cpos)
                # Missing anchor
        if args.flankNorm and \
            (( cpos_dat[0] not in range( region.start-(interp_bin+anchor_window), region.start-interp_bin ) ) or \
            ( cpos_dat[-1] not in range( region.end+interp_bin+1, region.end+(interp_bin+anchor_window)+1 ) )):
            filter_flag = True
            result = ("flankNorm", region.id)
    if not filter_flag:
        # Gather interpolated data
        interpolated_dmet_data = Interpolation(cpos_dat, dmet_dat, region.start, region.end, region.strand, reg_type, args).dat_proc()
        # cross-check number of interpolated features
        if len(interpolated_dmet_data) != args.num_interp_points:
            result = ('interpolation inconsistency in number of features', region.id)
        else:
            result = ('pass', region.id+'-'+sample_id+','+','.join(str(item) for item in interpolated_dmet_data) )
    del dmet_dat[:], cpos_dat[:], dmet_tmp[:], cpos_tmp[:]
    return result 

def exec_interp_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-a','--anno_file', default=argparse.SUPPRESS, required=True, help='region or gene annotation file')
    parser_required.add_argument('-t', dest='tag_inp', default=argparse.SUPPRESS, required=True, help='Tag for output')
    parser_required.add_argument('-i', '--input_list', default=argparse.SUPPRESS, required=True, help='List of sample pairs')
    parser_general = parser.add_argument_group('general arguments')
    parser_general.add_argument('--anno_type', default="gene_tss", 
                                 choices=["gene_tss", "enh"],
                                 help='region or gene annotation file')
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



###############################
### OLD CODE
###############################
def interp_genes_OLD(gene_list, dict_bed, sample_id, log_FH, args):
    out_data = ''
    interp_bin = args.ibin_inp
    reg_type='gene'
    for gene in gene_list:
        cpos_dat, cpos_tmp, dmet_dat, dmet_tmp = [], [], [], []
        tss = gene.tss()
        tes = gene.tes()
        print_to_log(log_FH, gene.id + "\n")
        # ----------------------------------------------------------------------------
        # Methylation based gene filters
        # 3. Genes with <40 CpGs assayed within +/-5kb of the TSS
        # 4. Genes with all CpGs within +/-5 kb of the TSS had <0.2 methylation change
        for cpos in range( (tss-interp_bin), (tss+interp_bin)+1 ):            
            if gene.chr in dict_bed and cpos in dict_bed[gene.chr]:
                dmet_tmp.append(dict_bed[gene.chr][cpos])
                cpos_tmp.append(cpos)
        if len(dmet_tmp) < args.min_gene_cpgs:
            print_to_log(log_FH, gene.id + ' has < ' +  str(args.min_gene_cpgs)  + ' CpGs assayed in ' + sample_id + '\n')
            continue
        if max(dmet_tmp) < args.min_gene_meth:
            print_to_log(log_FH, gene.id + ' has < ' + str(args.min_gene_meth) +' maximum methylation change in ' + sample_id + '\n')
            continue
        #-----------------------------------------------------------------------------    
        if not args.flankNorm: # Endpoint correction will be used in this case
            anchor_window = 0
#        if args.flankNorm:
        for cpos in range( gene.txStart-(interp_bin+anchor_window), gene.txEnd+(interp_bin+anchor_window)+1 ):
            if gene.chr in dict_bed and cpos in dict_bed[gene.chr]:
                dmet_dat.append(dict_bed[gene.chr][cpos])
                cpos_dat.append(cpos)
        # Missing anchor            
        if args.flankNorm and \
            (( cpos_dat[0] not in range( gene.txStart-(interp_bin+anchor_window), gene.txStart-interp_bin ) ) or \
            ( cpos_dat[-1] not in range( gene.txEnd+interp_bin+1, gene.txEnd+(interp_bin+anchor_window)+1 ) )):
            print_to_log(log_FH, gene.id + ' has not CpGs in anchor windows in ' + sample_id + '\n')
            continue 
        # Gather interpolated data
        interpolated_dmet_data = Interpolation(cpos_dat, dmet_dat, gene.txStart, gene.txEnd, gene.strand, reg_type, args).dat_proc()
        # cross-check number of interpolated features
        if len(interpolated_dmet_data) != args.num_interp_points:
            eprint('Inconsistent number of interpolation features!')
            exit() # Exit if number is inconsistent 
        # Write data
        out_data = out_data + '\n'+gene.id+'-'+sample_id+','+','.join(str(item) for item in interpolated_dmet_data)
        del dmet_dat[:], cpos_dat[:], dmet_tmp[:], cpos_tmp[:]
    return out_data


def interp_regions_OLD(region_list, dict_bed, sample_id, log_FH, args):
    out_data=''
    reg_type = 're'
    interp_bin = 500 # for model gene plots
    # Regulatory Elements
    for region in region_list: # regulatory elements
        cpos_dat, cpos_tmp, dmet_dat, dmet_tmp = [], [], [], []
        print_to_log(log_FH, region.id+'\n')
        # ----------------------------------------------------------------------------
        # Methylation based filters
        # 3. Genes with <2 CpGs assayed within enhancer
        #
        for cpos in range( (region.start-interp_bin), (region.end+interp_bin)+1 ):
            if region.chr in dict_bed and cpos in dict_bed[region.chr]:
                dmet_tmp.append(dict_bed[region.chr][cpos])
                cpos_tmp.append(cpos)
        if len(dmet_tmp) < args.min_reg_cpgs:
            print_to_log(log_FH,region.id + ' has < ' +  str(args.min_reg_cpgs) \
                        + ' CpGs assayed in ' + sample_id + '\n')
            continue
        #-----------------------------------------------------------------------------
        if not args.flankNorm:
            anchor_window = 0
        for cpos in range( region.start-(interp_bin+anchor_window), region.end+(interp_bin+anchor_window)+1 ):
            if region.chr in dict_bed and cpos in dict_bed[region.chr]:
                dmet_dat.append(dict_bed[region.chr][cpos])
                cpos_dat.append(cpos)
                # Missing anchor
        if args.flankNorm and \
            (( cpos_dat[0] not in range( region.start-(interp_bin+anchor_window), region.start-interp_bin ) ) or \
            ( cpos_dat[-1] not in range( region.end+interp_bin+1, region.end+(interp_bin+anchor_window)+1 ) )):
            print_to_log(log_FH,region.id + ' has no CpGs in anchor windows in ' + sample_id + '\n')
            continue
        # Gather interpolated data
        interpolated_dmet_data = Interpolation(cpos_dat, dmet_dat, region.start, region.end, region.strand, reg_type, args).dat_proc()
        # cross-check number of interpolated features
        if len(interpolated_dmet_data) != args.num_interp_points:
            print_to_log(log_FH,str(len(interpolated_dmet_data)))
            print_to_log(log_FH,'Inconsistent number of interpolation features!')
            exit()
        # Write data 
        out_data = out_data + '\n'+region.id+'-'+sample_id+','+','.join(str(item) for item in interpolated_dmet_data)
        del dmet_dat[:], cpos_dat[:], dmet_tmp[:], cpos_tmp[:]
    return out_data 