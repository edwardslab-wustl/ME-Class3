
import multiprocessing as mp
from itertools import repeat

from MEClass3.Interpolation_class import Interpolation
from MEClass3.io_functions import eprint 

def generate_out_header(num_pts, anno_type, data_type):
    header = ''
    if anno_type == 'tss':
        header = 'gene_id-sample_name'
    elif anno_type == 'enh':
        header = 'enh_loc-gene_id-sample_name'
    else:
        eprint("Can't recognize anno_type: " + anno_type + "\nCheck --anno_type specification in help.")
        exit()
    for pos in range(num_pts):
        header +=  ',' + "-".join([data_type,anno_type,str(pos)])
    return header + "\n"

def generate_param_comment(num_pts, region_size, data_type, anno_type):
    anno_id = data_type + '-' + anno_type
    param_data = '# Info:' + ",".join([anno_id,str(region_size),str(num_pts)])
    param_data += ' Data:'
    step_size = int(2 * region_size / num_pts )
    if anno_type == 'tss':
        pos = -region_size + int(step_size/2)
    else:
        pos = -region_size + int(step_size/2) # THERE SEEMS TO BE AN OFFSET ISSUE HERE FOR ENHANCERS!
    for pnt in range(num_pts):
        tmp = "-".join([anno_id,str(pnt)])
        param_data += (tmp + ',' + str(pos) + ';')
        pos += step_size
    return param_data

def generate_param_list(num_pts, region_size, data_type, anno_type):
    param_data = []
    param_data.append('#' + ",".join([anno_type,str(region_size),str(num_pts)]))
    step_size = int(2 * region_size / num_pts )
    if anno_type == 'tss':
        pos = -region_size + int(step_size/2)
    else:
        pos = -region_size + int(step_size/2) # THERE SEEMS TO BE AN OFFSET ISSUE HERE FOR ENHANCERS!
    for pnt in range(num_pts):
        tmp = "-".join([data_type,anno_type,str(pnt)])
        param_data.append(tmp + ',' + str(pos))
        pos += step_size
    return param_data

def add_fail (status, gene, dict):
    if status in dict:
        dict[status].append(gene)
    else:
        dict[status] = [ gene ]
    return dict

def format_fail_dict (fail_dict, fail_text):
    result = ''
    for status in fail_dict:
       num_fail = len(fail_dict[status])
       result += str(num_fail) + " " + fail_text + " failed to pass " + status + " filter.\n"
       result += ",".join(fail_dict[status]) + "\n\n"
    return result

def interp_list_sp(item_list, dict_bed, sample_id, fail_dict, args):
    final_results = []
    for item in item_list:
        if args.anno_type == 'tss':
            (status, result) = interp_gene(item, dict_bed, sample_id, args)
        elif args.anno_type == 'enh':
            (status, result) = interp_region(item, dict_bed, sample_id, args)
        else:
            eprint("Can't recognize anno_type: " + args.anno_type + "\nCheck --anno_type specification in help.")
            exit()
        if status == 'pass':
            final_results.append(result)
        else:
            fail_dict = add_fail(status, result, fail_dict)
    return "\n".join(final_results), fail_dict
    
def interp_list_mp(item_list, dict_bed, sample_id, fail_dict, args):
    arg_iterable = zip(item_list, repeat(dict_bed), repeat(sample_id), repeat(args))
    with mp.Pool(processes=args.num_proc) as pool:
        if args.anno_type == 'tss':
            results = pool.starmap(interp_gene, arg_iterable)
        elif args.anno_type == 'enh':
            results = pool.starmap(interp_region, arg_iterable)
        else:
            eprint("Can't recognize anno_type: " + args.anno_type + "\nCheck --anno_type specification in help.")
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
        interpolated_dmet_data = Interpolation(cpos_dat, dmet_dat, gene.txStart, gene.txEnd, gene.strand, reg_type, args).dat_proc()
        # cross-check number of interpolated features
        if len(interpolated_dmet_data) != args.num_interp_points:
            result = ('interpolation inconsistency in number of features', gene.id)
        else:
            result = ('pass', gene.id+'-'+sample_id+','+','.join(str(item) for item in interpolated_dmet_data) )
    del dmet_dat[:], cpos_dat[:], dmet_tmp[:], cpos_tmp[:]
    return result

def interp_region(region, dict_bed, sample_id, args):
    result = 'na'
    reg_type = 're'
    interp_bin = 500 # for model gene plots
    dict_cpg = dict_bed[region.chr]
    cpos_dat, cpos_tmp, dmet_dat, dmet_tmp = [], [], [], []
    # ----------------------------------------------------------------------------
    # Regulatory Elements
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


###############################
### OLD CODE
###############################
# def interp_genes_OLD(gene_list, dict_bed, sample_id, log_FH, args):
#     out_data = ''
#     interp_bin = args.ibin_inp
#     reg_type='gene'
#    for gene in gene_list:
#        cpos_dat, cpos_tmp, dmet_dat, dmet_tmp = [], [], [], []
#        tss = gene.tss()
#        tes = gene.tes()
#        print_to_log(log_FH, gene.id + "\n")
#        # ----------------------------------------------------------------------------
#        # Methylation based gene filters
#        # 3. Genes with <40 CpGs assayed within +/-5kb of the TSS
#        # 4. Genes with all CpGs within +/-5 kb of the TSS had <0.2 methylation change
#        for cpos in range( (tss-interp_bin), (tss+interp_bin)+1 ):            
#            if gene.chr in dict_bed and cpos in dict_bed[gene.chr]:
#                dmet_tmp.append(dict_bed[gene.chr][cpos])
#                cpos_tmp.append(cpos)
#        if len(dmet_tmp) < args.min_gene_cpgs:
#            print_to_log(log_FH, gene.id + ' has < ' +  str(args.min_gene_cpgs)  + ' CpGs assayed in ' + sample_id + '\n')
#            continue
#        if max(dmet_tmp) < args.min_gene_meth:
#            print_to_log(log_FH, gene.id + ' has < ' + str(args.min_gene_meth) +' maximum methylation change in ' + sample_id + '\n')
#            continue
#        #-----------------------------------------------------------------------------    
#        if not args.flankNorm: # Endpoint correction will be used in this case
#            anchor_window = 0
##        if args.flankNorm:
#        for cpos in range( gene.txStart-(interp_bin+anchor_window), gene.txEnd+(interp_bin+anchor_window)+1 ):
#            if gene.chr in dict_bed and cpos in dict_bed[gene.chr]:
#                dmet_dat.append(dict_bed[gene.chr][cpos])
#                cpos_dat.append(cpos)
#        # Missing anchor            
#        if args.flankNorm and \
#            (( cpos_dat[0] not in range( gene.txStart-(interp_bin+anchor_window), gene.txStart-interp_bin ) ) or \
#            ( cpos_dat[-1] not in range( gene.txEnd+interp_bin+1, gene.txEnd+(interp_bin+anchor_window)+1 ) )):
#            print_to_log(log_FH, gene.id + ' has not CpGs in anchor windows in ' + sample_id + '\n')
#            continue 
#        # Gather interpolated data
#        interpolated_dmet_data = Interpolation(cpos_dat, dmet_dat, gene.txStart, gene.txEnd, gene.strand, reg_type, args).dat_proc()
#        # cross-check number of interpolated features
#        if len(interpolated_dmet_data) != args.num_interp_points:
#            eprint('Inconsistent number of interpolation features!')
#            exit() # Exit if number is inconsistent 
#        # Write data
#        out_data = out_data + '\n'+gene.id+'-'+sample_id+','+','.join(str(item) for item in interpolated_dmet_data)
#        del dmet_dat[:], cpos_dat[:], dmet_tmp[:], cpos_tmp[:]
#    return out_data
#
#
#def interp_regions_OLD(region_list, dict_bed, sample_id, log_FH, args):
#    out_data=''
#    reg_type = 're'
#    interp_bin = 500 # for model gene plots
#    # Regulatory Elements
#    for region in region_list: # regulatory elements
#        cpos_dat, cpos_tmp, dmet_dat, dmet_tmp = [], [], [], []
#        print_to_log(log_FH, region.id+'\n')
#        # ----------------------------------------------------------------------------
#        # Methylation based filters
#        # 3. Genes with <2 CpGs assayed within enhancer
#        #
#        for cpos in range( (region.start-interp_bin), (region.end+interp_bin)+1 ):
#            if region.chr in dict_bed and cpos in dict_bed[region.chr]:
#                dmet_tmp.append(dict_bed[region.chr][cpos])
#                cpos_tmp.append(cpos)
#        if len(dmet_tmp) < args.min_reg_cpgs:
#            print_to_log(log_FH,region.id + ' has < ' +  str(args.min_reg_cpgs) \
#                        + ' CpGs assayed in ' + sample_id + '\n')
#            continue
#        #-----------------------------------------------------------------------------
#        if not args.flankNorm:
#            anchor_window = 0
#        for cpos in range( region.start-(interp_bin+anchor_window), region.end+(interp_bin+anchor_window)+1 ):
#            if region.chr in dict_bed and cpos in dict_bed[region.chr]:
#                dmet_dat.append(dict_bed[region.chr][cpos])
#                cpos_dat.append(cpos)
#                # Missing anchor
#        if args.flankNorm and \
#            (( cpos_dat[0] not in range( region.start-(interp_bin+anchor_window), region.start-interp_bin ) ) or \
#            ( cpos_dat[-1] not in range( region.end+interp_bin+1, region.end+(interp_bin+anchor_window)+1 ) )):
#            print_to_log(log_FH,region.id + ' has no CpGs in anchor windows in ' + sample_id + '\n')
#            continue
#        # Gather interpolated data
#        interpolated_dmet_data = Interpolation(cpos_dat, dmet_dat, region.start, region.end, region.strand, reg_type, args).dat_proc()
#        # cross-check number of interpolated features
#        if len(interpolated_dmet_data) != args.num_interp_points:
#            print_to_log(log_FH,str(len(interpolated_dmet_data)))
#            print_to_log(log_FH,'Inconsistent number of interpolation features!')
#            exit()
#        # Write data 
#        out_data = out_data + '\n'+region.id+'-'+sample_id+','+','.join(str(item) for item in interpolated_dmet_data)
#        del dmet_dat[:], cpos_dat[:], dmet_tmp[:], cpos_tmp[:]
#    return out_data 