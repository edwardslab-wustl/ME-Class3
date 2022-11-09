
from MEClass3.overlap_genes_enhancers_functions import read_geneHancer
from MEClass3.overlap_genes_enhancers_functions import assign_enh_to_gene_by_score
from MEClass3.overlap_genes_enhancers_functions import format_enh_anno


file = "meclass3_min_example_chr17/enh_db_test/GeneHancer.test.txt"

enh_dict = read_geneHancer(file)
enh_dict_list = assign_enh_to_gene_by_score(enh_dict, 1)
results = format_enh_anno(enh_dict_list, enh_dict)


print('done')
