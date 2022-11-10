from collections import defaultdict

from dataclasses import dataclass

from MEClass3.io_functions import RegionAnno
from MEClass3.overlap_genes_enhancers_functions import add_dict_set


#### dataclasses ####
@dataclass()
class GeneEnhSet:
    gene_id: str
    enh_set: set
    
    def num_enh(self) -> int:
        return len(self.enh_set)
    
    def add_enh(self, enh):
        self.enh_set.add(enh)
        return
    
    def discard_enh(self, enh):
        self.enh_set.discard(enh)
        return
    
    def top_n_enh(self, n) -> list:
        return_list = sorted(self.enh_set, key=lambda x: x.connect_score)
        if self.num_enh() > n:
            return_list = return_list[0:n]
        return return_list
        
    def min_score_enh(self) -> RegionAnno:
        min_score = -1
        min_enh = 'na'
        for enh in self.enh_set:
            if min_enh == 'na':
                min_score = enh.connect_score
                min_enh = enh
            elif enh.connect_score < min_score:
                min_score = enh.connect_score
                min_enh = enh
        return min_enh        
    
    def add_gene_info(self, gene_dict):
        for enh in self.enh_set:
            if enh.gene in gene_dict:
                gene = gene_dict[enh.gene]
                new_enh = RegionAnno(enh.id, enh.region_id,
                                    gene.id,
                                    enh.chr, enh.start, enh.end,
                                    gene.strand,
                                    gene.txStart,
                                    gene.txEnd,
                                    enh.connect_score)
                self.discard_enh(enh)
                self.add_enh(new_enh)
        return
                
 
#def add_gene_atlas_data(data_dict, new_data_dict):
#    for symbol, enh_list in new_data_dict.items():
#        for enh in enh_list:
#            if symbol in data_dict:
#                data_dict[symbol].add_enh(enh)
#            else:
#                data_dict[symbol] = GeneEnhSet(symbol, {enh})
#    return data_dict
       
def read_geneHancer_2(file, result_dict, enh_score_threshold=0., connection_score_threshold=0.):
    #result_enh_list = []
    #result_enh_dict = dict()
    with open(file, 'r') as FH:
        line_count = 0
        for line in FH:
            line_count += 1
            if line_count == 1:
                continue  # skip header
            # chrom,source,feature name,start,end,score,strand,frame,attributes
            line_data = line.strip().split(',')
            enh_score = float(line_data[5])
            if enh_score >= enh_score_threshold:
                attribute_list = line_data[8].split(';')
                enh_id = "na"
                gene_id = "na"
                connection_score = "na"
                for attr in attribute_list:
                    (param, val) = attr.split('=')
                    if param == "genehancer_id":
                        enh_id = val
                    elif param == "connected_gene":
                        gene_id = val
                    elif param == "score":
                        connection_score = float(val)
                        if connection_score >= connection_score_threshold:
                            coord=f"{line_data[0]}:{line_data[3]}-{line_data[4]}"
                            enh = RegionAnno( enh_id,
                                              coord,
                                              gene_id,
                                              line_data[0],  # chrom
                                              line_data[3],  # start
                                              line_data[4],  # end
                                              'na',  # strand
                                              -1,  # txStart
                                              -1,  # txEnd 
                                              connect_score=connection_score
                                              )
                            #result_enh_dict[enh_id] = enhancer
                        if gene_id in result_dict:
                            result_dict[gene_id].add_enh(enh)
                        else:
                            result_dict[gene_id] = GeneEnhSet(gene_id, {enh})
    return result_dict

def read_enh_atlas_file(file, result_dict, label='enh_', score_threshold=0.):
    #chr1:1266920-1267010_ENSG00000242485$MRPL20$chr1$1342693$-      1.239509
    #result_dict = defaultdict(set)
    with open(file, 'r') as FH:
        count = 0
        for line in FH:
            count += 1
            line_data = line.strip().split()
            score = float(line_data[1])
            if score >= score_threshold:
                (coord, gene_info) = line_data[0].split('_')
                (ensg_id, symbol, gene_chrom, tss, strand) = gene_info.split('$')
                (enh_chr, position_range) = coord.split(":")
                (enh_start, enh_end) = position_range.split('-')
                enh = RegionAnno(label + str(count),
                                 coord,
                                 symbol,
                                 enh_chr,
                                 int(enh_start),
                                 int(enh_end),
                                 strand,
                                 tss,
                                 tss,
                                 score)
                #add_dict_set(result_dict, symbol, enh)
                if symbol in result_dict:
                    result_dict[symbol].add_enh(enh)
                else:
                    result_dict[symbol] = GeneEnhSet(symbol, {enh})
    return result_dict
                
                