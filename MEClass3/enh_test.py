
#from MEClass3.combine_enhancer_atlas import read_enh_atlas_file
from MEClass3.combine_enhancer_atlas import read_geneHancer_2
from MEClass3.io_functions import read_gene_file
from MEClass3.overlap_genes_enhancers_functions import index_gene_list

#from MEClass3.combine_enhancer_atlas import add_gene_atlas_data

#file_1="meclass3_min_example_chr17/enh_atlas/Liver_EP.chr17.txt"
#file_2="meclass3_min_example_chr17/enh_atlas/Skeletal_muscle_EP.chr17.txt"
gene_file="meclass3_min_example_chr17/refGene.gene.chr17"

#file_list = [file_1, file_2]
file_list = [ "meclass3_min_example_chr17/enh_db_test/GeneHancer.chr17" ]
results = dict()
for file in file_list:
    #results = read_enh_atlas_file(file, results)
    results = read_geneHancer_2(file, results)
    
gene_list = read_gene_file(gene_file)
gene_dict = index_gene_list(gene_list)

for symbol, result in results.items():
    result.add_gene_info(gene_dict)
    if result.num_enh() > 2:
        #print(f"gene: {result.gene_id}\nset: {result.enh_set}")
        count = 0
        for enh in result.top_n_enh(5):
            count += 1
            label = 'enh_' + str(count)
            print(f"{label}\t{enh.format_output()}\n")
#results = add_gene_atlas_data(results, data_1)
#results = add_gene_atlas_data(results, data_2)

print('done')


#from MEClass3.overlap_genes_enhancers_functions import read_geneHancer
#from MEClass3.overlap_genes_enhancers_functions import assign_enh_to_gene_by_score
#from MEClass3.overlap_genes_enhancers_functions import format_enh_anno
#
#
#file = "meclass3_min_example_chr17/enh_db_test/GeneHancer.test.txt"
#
#enh_dict = read_geneHancer(file)
#enh_dict_list = assign_enh_to_gene_by_score(enh_dict, 1)
#results = format_enh_anno(enh_dict_list, enh_dict)


#print('done')
