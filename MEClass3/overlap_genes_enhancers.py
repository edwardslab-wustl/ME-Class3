
import argparse

from MEClass3.io_functions import read_enhancer_file
from MEClass3.io_functions import read_gene_file
from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import format_args_to_print

def exec_overlap_genes_enhancers(args):
    with open(args.logfile, 'w') as log_FH:   
        print_to_log(log_FH, format_args_to_print(args))
        #old: n_end is downstream, p_end is upstream of TSS
        up_start_idx = args.up_start_idx 
        up_end_idx = args.up_end_idx 
        dn_start_idx = args.dn_start_idx 
        dn_end_idx = args.dn_end_idx 
        promoter_region = args.promoter_region
        gene_file = args.gene_file
        enhancer_file = args.enh_file
        index_size = args.index_size
        up_label = 'up_'
        dn_label = 'dn_'
 
        dict_bed, dict_bed_idx =  read_enhancer_file(enhancer_file, delimiter=index_size)
        gene_list = read_gene_file(gene_file)
        out_data = ''
        for gene in gene_list:
            print_to_log(log_FH, gene.id + "\n")
            cnt_dn = 0
            tss = gene.tss()
            txStart_idx = int(int(gene.txStart/index_size))
            for info in dict_bed[gene.chr][ (dict_bed_idx[gene.chr].get( txStart_idx, 1))-1 : ]:
                if info[0] > int( tss + promoter_region ):
                    cnt_dn += 1
                    if cnt_dn in range (dn_start_idx, dn_end_idx+1):
                        if gene.strand == '+':
                            enh_loc = dn_label+str(cnt_dn)
                        elif gene.strand == '-':
                            enh_loc = up_label+str(cnt_dn)
                        new_data = "\t".join([enh_loc, gene.id, gene.chr, gene.strand, str(gene.txStart), str(gene.txEnd), str(info[0]),str(info[1])])
                        out_data = out_data + new_data + "\n"
                    if cnt_dn > dn_end_idx:
                        break
            # p-end
            # Get the index of first p_end reg.
            for info in dict_bed[gene.chr][ (dict_bed_idx[gene.chr].get( txStart_idx, 1))-1 : ]: 
                if info[1] < int( tss-promoter_region ):
                    prev_index = dict_bed[gene.chr].index(info)
                    continue
                if info[1] >= int( tss-promoter_region ):
                    up_search_index = prev_index
                    break
            cnt_up = 0
            write_dict ={}
            for info in reversed( dict_bed[gene.chr][ : up_search_index+2 ] ):
                if info[1] < int( tss-promoter_region ):
                    cnt_up += 1
                    if cnt_up in range(up_start_idx, up_end_idx+1):
                        if gene.strand == '+':
                            enh_loc = up_label+str(cnt_up)
                        elif gene.strand == '-':
                            enh_loc = dn_label+str(cnt_up)
                    if (gene.id, enh_loc) not in write_dict: #p_end_5 genes were written twice, this is a quick hack to fix
                        new_data = "\t".join([enh_loc, gene.id, gene.chr, gene.strand, str(gene.txStart), str(gene.txEnd), str(info[0]),str(info[1])])
                        out_data = out_data + new_data + "\n"
                        write_dict[(gene.id,enh_loc)] = 1 
                    if cnt_up > up_end_idx:
                        break
        with open(args.outFile, 'w') as out_FH:
            out_FH.write(out_data)

def exec_overlap_genes_enhancers_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-g','--gene_file', required=True,
        default=argparse.SUPPRESS,
        help='gene annotation file in refSeq format. Use refGene.txt file from the UCSC genome browser.')
    parser_required.add_argument('-e','--enh_file', required=True,
        default=argparse.SUPPRESS,
        help='enhancer annotation file')
    parser.add_argument('-o', '--outFile', default='enh_gene.overlap.txt', help='output file')
    parser.add_argument('--dn_start_idx', type=int, default=1, 
        help='downstream start index. 1 means to start with the closest enhancer downstream of the TSS.')
    parser.add_argument('--dn_end_idx', type=int, default=5, 
        help='downstream end index. 5 means stop at 5 enhancers downstream of the TSS.')
    parser.add_argument('--up_start_idx', type=int, default=1, 
        help='upstream start index. 1 means to start with the closest enhancer upstream of the TSS.')
    parser.add_argument('--up_end_idx', type=int, default=5, 
        help='upstream end index. 5 means to stop at 5 enhancers upstream of the TSS.')
    parser.add_argument('--promoter_region', type=int, default=5000,
        help='ignore all enhancers within +/- promoter_region bp from the TSS. By default excludes all enhancers +/- 5000bp from TSS.')
    parser.add_argument('--index_size', type=int, default=1000,
        help='Size of index used for overlapping enhancers and promoters. Set larger than enhancer/feature size.')
    parser.add_argument('--logfile',  default='enh_gene_overlap.log', help='log file')
    parser._action_groups.reverse()
    return(parser)