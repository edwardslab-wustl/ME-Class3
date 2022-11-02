
import argparse
from os import getcwd

from MEClass3.io_functions import read_enhancer_file
from MEClass3.io_functions import read_gene_file
from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import format_args_to_print

def exec_overlap_genes_enhancers(args):
    # This part of code generate proximity list of regulatory elements
    # add orientation relative to TSS, in addition to n and p.
    p_sidx = args.pst_inp # 1
    p_eidx = args.ped_inp # 2
    n_sidx = args.nst_inp # 1
    n_eidx = args.ned_inp # 2
    twin = args.twin_inp
    gene_file = args.anno_file
    enhancer_file = args.enh_file
    #path_to_input = getcwd()
    delimiter = 1000 #index size
 
    with open(args.logfile, 'w') as log_FH:   
        print_to_log(log_FH, format_args_to_print(args))
        dict_bed, dict_bed_idx =  read_enhancer_file(enhancer_file, delimiter=delimiter)
        gene_list = read_gene_file(gene_file)
        out_data = ''
        for gene in gene_list:
            print_to_log(log_FH, gene.id + "\n")
            cnt_ne = 0
            tss = gene.tss()
            txStart_idx = int(int(gene.txStart/delimiter))
            for info in dict_bed[gene.chr][ (dict_bed_idx[gene.chr].get( txStart_idx, 1))-1 : ]:
                if info[0] > int( tss+twin ):
                    cnt_ne += 1
                    if cnt_ne in range (n_sidx, n_eidx+1):
                        if gene.strand == '+':
                            enh_loc = 'n_end_'+str(cnt_ne)
                        elif gene.strand == '-':
                            enh_loc = 'p_end_'+str(cnt_ne)
                        new_data = "\t".join([enh_loc, gene.id, gene.chr, gene.strand, str(gene.txStart), str(gene.txEnd), str(info[0]),str(info[1])])
                        out_data = out_data + new_data + "\n"
                    if cnt_ne > n_eidx:
                        break
            # p-end
            # Get the index of first p_end reg.
            for info in dict_bed[gene.chr][ (dict_bed_idx[gene.chr].get( txStart_idx, 1))-1 : ]: 
                if info[1] < int( tss-twin ):
                    prev_index = dict_bed[gene.chr].index(info)
                    continue
                if info[1] >= int( tss-twin ):
                    p_search_index = prev_index
                    break
            cnt_pe = 0
            write_dict ={}
            for info in reversed( dict_bed[gene.chr][ : p_search_index+2 ] ):
                if info[1] < int( tss-twin ):
                    cnt_pe += 1
                    if cnt_pe in range(p_sidx, p_eidx+1):
                        if gene.strand == '+':
                            enh_loc = 'p_end_'+str(cnt_pe)
                        elif gene.strand == '-':
                            enh_loc = 'n_end_'+str(cnt_pe)
                    if (gene.id, enh_loc) not in write_dict: #p_end_5 genes were written twice, this is a quick hack to fix
                        new_data = "\t".join([enh_loc, gene.id, gene.chr, gene.strand, str(gene.txStart), str(gene.txEnd), str(info[0]),str(info[1])])
                        out_data = out_data + new_data + "\n"
                        write_dict[(gene.id,enh_loc)] = 1 
                    if cnt_pe > p_eidx:
                        break
        with open(args.outFile, 'w') as out_FH:
            out_FH.write(out_data)

def exec_overlap_genes_enhancers_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-a','--anno_file', required=True,
        default=argparse.SUPPRESS,
        help='gene annotation file')
    parser_required.add_argument('-e','--enh_file', required=True,
        default=argparse.SUPPRESS,
        help='enhancer annotation file')
    parser.add_argument('-o', '--outFile', default='enh_gene.overlap.txt', help='output file')
    parser.add_argument('--pst', action='store', dest='pst_inp', type=int, default=1, help='pend start index')
    parser.add_argument('--ped', action='store', dest='ped_inp', type=int, default=2, help='pend end index')
    parser.add_argument('--nst', action='store', dest='nst_inp', type=int, default=1, help='nend start index')
    parser.add_argument('--ned', action='store', dest='ned_inp', type=int, default=2, help='nend end index')
    parser.add_argument('--twn', action='store', dest='twin_inp', type=int, default=5000, help='nend end index')
    parser.add_argument('--logfile', dest='logfile', default='enh_gene_overlap.log', help='log file')
    parser._action_groups.reverse()
    return(parser)