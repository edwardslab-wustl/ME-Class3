
import argparse

from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import format_args_to_print
from overlap_genes_enhancers_functions import overlap_enh_genes_by_distance
from overlap_genes_enhancers_functions import read_geneHancer
from overlap_genes_enhancers_functions import assign_enh_to_gene_by_score
from overlap_genes_enhancers_functions import format_enh_anno
from overlap_genes_enhancers_functions import add_gene_info

def exec_overlap_genes_enhancers(args):
    with open(args.logfile, 'w') as log_FH:   
        print_to_log(log_FH, format_args_to_print(args))
        results = []
        if args.type == 'distance':
            results = overlap_enh_genes_by_distance(args.gene_file,
                                                    args.enh_file,
                                                    args.up_start_idx,
                                                    args.up_end_idx,
                                                    args.dn_start_idx,
                                                    args.dn_end_idx,
                                                    args.promoter_region,
                                                    args.index_size)
        elif args.type == 'genehancer_score':
            enh_dict = read_geneHancer(args.enh_file)
            enh_dict_list = assign_enh_to_gene_by_score(enh_dict, args.score_num_enh)
            enh_dict_anno, total_num, convert_num = add_gene_info(enh_dict, args.gene_file)
            print_to_log(log_FH, f"converted ids for {convert_num} genes out of {total_num} total genes.")
            results = format_enh_anno(enh_dict_list, enh_dict_anno)
        out_data = [ "\t".join(result) for result in results ]
        with open(args.outFile, 'w') as out_FH:
            out_FH.write("\n".join(out_data) + "\n")

def exec_overlap_genes_enhancers_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-g','--gene_file', required=True,
        default=argparse.SUPPRESS,
        help='gene annotation file in refSeq format. Use refGene.txt file from the UCSC genome browser.')
    parser_required.add_argument('-e','--enh_file', required=True,
        default=argparse.SUPPRESS,
        help='enhancer annotation file')
    parser.add_argument('-t', '--type', default='distance',
        choices=['distance', 'genehancer_score'],
        help='type of enhancer file and matching metric. Distance expects bed \
            file of enhancer locations. genehancer_score expects genehancer \
            .csv file with genehancer gene-enhancer connectivity info.' )
    parser.add_argument('-o', '--outFile',
        default='enh_gene.overlap.txt',
        help='output file')
    parser.add_argument('--logfile',  default='enh_gene_overlap.log', help='log file')
    parser_distance = parser.add_argument_group('distance arguments')
    parser_distance.add_argument('--dn_start_idx', type=int, default=1, 
        help='downstream start index. 1 means to start with the closest enhancer downstream of the TSS.')
    parser_distance.add_argument('--dn_end_idx', type=int, default=5, 
        help='downstream end index. 5 means stop at 5 enhancers downstream of the TSS.')
    parser_distance.add_argument('--up_start_idx', type=int, default=1, 
        help='upstream start index. 1 means to start with the closest \
            enhancer upstream of the TSS.')
    parser_distance.add_argument('--up_end_idx', type=int, default=5, 
        help='upstream end index. 5 means to stop at 5 enhancers upstream of the TSS.')
    parser_distance.add_argument('--promoter_region', type=int, default=5000,
        help='ignore all enhancers within +/- promoter_region bp from the TSS. \
            By default excludes all enhancers +/- 5000bp from TSS.')
    parser_distance.add_argument('--index_size', type=int, default=1000,
        help='Size of index used for overlapping enhancers and promoters. Set larger than enhancer/feature size.')
    parser_score = parser.add_argument_group('score arguments')
    parser_score.add_argument('--score_num_enh', type=int, default=5,
        help='number of enhancers to use if using a scoring metric such as genehancer. Ignored if --type = distance.')
    parser._action_groups.reverse()
    return(parser)