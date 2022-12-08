
from collections import defaultdict
from dataclasses import dataclass

from MEClass3.io_functions import RegionAnno
from MEClass3.io_functions import read_enhancer_file
from MEClass3.io_functions import read_gene_file

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
        return_list = sorted(self.enh_set, reverse=True, key=lambda x: x.connect_score)
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
        new_enh_set = set()
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
                new_enh_set.add(new_enh)
            else:
                new_enh_set.add(enh)
        self.enh_set = new_enh_set       
        return
 

def overlap_enh_genes_by_score(gene_file, enhancer_file, enh_type, num_enh, base_label):
    results = []
    file_list = enhancer_file.split(',')
    tmp_results = dict()
    for file in file_list:
        if enh_type == "genehancer":
            tmp_results = read_geneHancer(file, tmp_results)
        elif enh_type == "enhancerAtlas":
            tmp_results = read_enh_atlas_file(file, tmp_results)
        gene_list = read_gene_file(gene_file)
        gene_dict = index_gene_list(gene_list)
    for symbol, result in tmp_results.items():
        result.add_gene_info(gene_dict)
        if result.num_enh() >= num_enh:
            count = 0
            for enh in result.top_n_enh(num_enh):
                count += 1
                label = base_label + '_' + str(count)
                out_list = [label]
                out_list.extend(enh.format_output())
                results.append(out_list)
                #results.append(f"{label}\t{enh.format_output()}\n")
    return results

def extract_enh_data_as_dict(enhancer_dict, index_range = 1000, unique_id = None):
    #cnt = 0
    enh_index = RegionIndex(defaultdict(dict), index_range)
    for symbol, enhancerSet in enhancer_dict.items():
        enh_index.add_set_to_index(enhancerSet.enh_set)
        
        #for enhancer in enhancerSet.enh_set:
        #    enh_index.add_to_index(enhancer)
        #    cnt += 1
    return enh_index

@dataclass
class RegionIndex:
    reg_index: defaultdict(dict)
    index_range: int = 1000
    
    def add_to_index(self, region):
        pos1 = int(region.start / self.index_range)
        pos2 = int(region.end / self.index_range) + 1
        for i in range(pos1, pos2):
            add_dict_set_2keys(self.reg_index, region.chr, i, region)
        return 
    
    def add_set_to_index(self, region_set):
        for region in region_set:
            self.add_to_index(region)
        return
    
    def merge_to_index_max(self, region_index):
        for chr in region_index.reg_index.keys():
            if chr in self.reg_index:
                for index, region_set in region_index.reg_index[chr].items():
                    if index in self.reg_index[chr]:
#                        new_set = overlap_region_sets_max(region_set, self.reg_index[chr][index])
#                        self.reg_index[chr][index] = new_set
                        self.add_set_to_index(region_set)
                    else:
                        self.add_set_to_index(region_set)
            else:
                for index, region_set in region_index.reg_index[chr].items():
                    self.add_set_to_index(region_set)
        return
    
def overlap_region_sets_max(region_set_new, region_set_ref):
    new_set = set()
    ref_set = set() 
    for region_new in region_set_new:
        match_flag = False
        for region_ref in region_set_ref:
            new_match_flag, new_regions = merge_coords_max(region_new, region_ref)
            ref_set.add(region_ref)
            new_set.update(new_regions)
            if new_match_flag:
                match_flag = True
        if not match_flag:
            new_set.add(region_new)
    new_set.update( {region for region in region_set_ref if region not in ref_set} )
    return new_set

def merge_coords_max(region_1, region_2):
    return_region_set = set()
    merge_flag = False
    if region_1.start <= region_2.start and region_2.end <= region_1.end:
        return_region_set.add(region_1)
        merge_flag = True
    elif region_2.start <= region_1.start and region_1.end <= region_2.end:
        return_region_set.add(region_2)
        merge_flag = True
    elif region_2.start <= region_1.start and region_1.start <= region_2.end:
        region_1.start = region_2.start
        return_region_set.add(region_1)
        merge_flag = True
    elif region_2.start <= region_1.end and region_1.end <= region_2.end:
        region_1.end = region_2.end
        return_region_set.add(region_1)
        merge_flag = True
    else:
        #return_region_set.add(region_2)
        pass
    return merge_flag, return_region_set
        
def add_dict_set_2keys(dictionary, key1, key2, value):
    if key1 in dictionary and key2 in dictionary[key1]:
        dictionary[key1][key2].add(value)
    else:
        dictionary[key1][key2] = {value}
    return
    
def add_dict_set(dictionary, key, value):
    if key in dictionary:
        dictionary[key].add(value)
    else:
        dictionary[key] = {value}
    return

def overlap_enh_genes_by_distance_2( gene_file, enhancer_file_list,
                                     up_start_index, up_end_index, dn_start_index, dn_end_index,
                                     promoter_region, index_size, max_distance):
    gene_list = read_gene_file(gene_file)
    enh_index = RegionIndex(defaultdict(dict))
    for enh_file in enhancer_file_list.split(','):
        enh_data_dict = dict()
        enh_data_dict = read_enh_atlas_file(enh_file, enh_data_dict)
        enh_obj_tmp = extract_enh_data_as_dict(enh_data_dict, index_range=index_size)
        enh_index.merge_to_index_max(region_index=enh_obj_tmp)
    overlap_results = find_enh_gene_list_overlap( gene_list, enh_index, up_end_index, dn_end_index, promoter_region, max_distance)
    results = format_distance_results(overlap_results, up_start_index, up_end_index, dn_start_index, dn_end_index)
    return results

def overlap_enh_genes_by_distance_new_bed( gene_file, enhancer_file_list,
                                     up_start_index, up_end_index, dn_start_index, dn_end_index,
                                     promoter_region, index_size, max_distance):
    gene_list = read_gene_file(gene_file)
    enh_index = RegionIndex(defaultdict(dict))
    for enh_file in enhancer_file_list.split(','):
        enh_data_dict = dict()
        enh_data_dict = read_enh_bed_file(enh_file, enh_data_dict)
        enh_obj_tmp = extract_enh_data_as_dict(enh_data_dict, index_range=index_size)
        enh_index.merge_to_index_max(region_index=enh_obj_tmp)
    overlap_results = find_enh_gene_list_overlap( gene_list, enh_index, up_end_index, dn_end_index, promoter_region, max_distance)
    results = format_distance_results(overlap_results, up_start_index, up_end_index, dn_start_index, dn_end_index)
    return results

def find_enh_gene_list_overlap(gene_list, enhancer_index, num_up, num_dn, promoter_region, max_distance = 1000000):
    results = defaultdict(list)
    for gene in gene_list:
        enh_list_up = gene.identify_nearest_up_objects(num_up, promoter_region, max_distance, enhancer_index)
        enh_list_dn = gene.identify_nearest_dn_objects(num_dn, promoter_region, max_distance, enhancer_index)
        if len(enh_list_up) >= num_up and len(enh_list_dn) >= num_dn:
            results[gene] = [enh_list_up, enh_list_dn]
    return results

def format_distance_results( distance_results, up_start_index, up_end_index, dn_start_index, dn_end_index ):
    results = []
    for gene, [enh_list_up, enh_list_dn] in distance_results.items():
        if len(enh_list_up) >= up_end_index and len(enh_list_dn) >= dn_end_index:
            if gene.strand == '-':
                enh_list_dn,enh_list_up = enh_list_up,enh_list_dn
            for i in range(up_start_index - 1, up_end_index):
                enh_id = f"enh_up_{i+1}"
                #results.append( "\t".join([ enh_id, gene.id, gene.chr, gene.strand,
                #                    str(gene.txStart), str(gene.txEnd),
                #                    str(enh_list_up[i].start), str(enh_list_up[i].end) ]) )
                results.append( [ enh_id, gene.id, gene.chr, gene.strand,
                                  str(gene.txStart), str(gene.txEnd),
                                  str(enh_list_up[i].start), str(enh_list_up[i].end) ] )
            for i in range(dn_start_index - 1, dn_end_index):
                enh_id = f"enh_dn_{i+1}"
                #results.append( "\t".join([ enh_id, gene.id, gene.chr, gene.strand,
                #                    str(gene.txStart), str(gene.txEnd),
                #                    str(enh_list_dn[i].start), str(enh_list_dn[i].end) ]) )
                results.append( [ enh_id, gene.id, gene.chr, gene.strand,
                                  str(gene.txStart), str(gene.txEnd),
                                  str(enh_list_dn[i].start), str(enh_list_dn[i].end) ] )
    return results
            
    
#def add_dict_list(dictionary, key, value):
#    if key in dictionary:
#        dictionary[key].append(value)
#    else:
#        dictionary[key] = [value]
#    return


def overlap_enh_genes_by_distance(gene_file, enhancer_file, up_start_idx, up_end_idx, dn_start_idx, dn_end_idx, promoter_region, index_size):
    results = []
    #old: n_end is downstream, p_end is upstream of TSS
    up_label = 'up_'
    dn_label = 'dn_'
    dict_bed, dict_bed_idx =  read_enhancer_file(enhancer_file, delimiter=index_size)
    gene_list = read_gene_file(gene_file)
    for gene in gene_list:
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
                    result = [enh_loc, gene.id, gene.chr, gene.strand,
                              str(gene.txStart), str(gene.txEnd), str(info[0]),str(info[1])]
                    results.append(result)
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
        write_dict = {}
        for info in reversed( dict_bed[gene.chr][ : up_search_index+2 ] ):
            if info[1] < int( tss-promoter_region ):
                cnt_up += 1
                if cnt_up in range(up_start_idx, up_end_idx+1):
                    if gene.strand == '+':
                        enh_loc = up_label+str(cnt_up)
                    elif gene.strand == '-':
                        enh_loc = dn_label+str(cnt_up)
                if (gene.id, enh_loc) not in write_dict:  # p_end_5 genes were written twice, this is a quick hack to fix
                    result = [enh_loc, gene.id, gene.chr, gene.strand,
                              str(gene.txStart), str(gene.txEnd), str(info[0]),str(info[1])]
                    results.append(result)
                    write_dict[(gene.id,enh_loc)] = 1 
                if cnt_up > up_end_idx:
                    break
    return results

def index_gene_list(gene_list):
    result_dict = dict()
    for gene in gene_list:
        result_dict[gene.symbol] = gene
    return result_dict


def read_geneHancer(file, result_dict, enh_score_threshold=0., connection_score_threshold=0.):
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
                        if gene_id in result_dict:
                            result_dict[gene_id].add_enh(enh)
                        else:
                            result_dict[gene_id] = GeneEnhSet(gene_id, {enh})
    return result_dict

def read_enh_atlas_file(file, result_dict, label='enh_', score_threshold=0.):
    #chr1:1266920-1267010_ENSG00000242485$MRPL20$chr1$1342693$-      1.239509
    with open(file, 'r') as FH:
        count = 0
        for line in FH:
            count += 1
            line_data = line.strip().split()
            score = float(line_data[1])
            if score >= score_threshold:
                (coord, gene_info) = line_data[0].split('_', 1)
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
                if symbol in result_dict:
                    result_dict[symbol].add_enh(enh)
                else:
                    result_dict[symbol] = GeneEnhSet(symbol, {enh})
    return result_dict

def read_enh_bed_file(file, result_dict, label='enh_'):
    #chr1:1266920-1267010_ENSG00000242485$MRPL20$chr1$1342693$-      1.239509
    with open(file, 'r') as FH:
        count = 0
        tss = -1
        symbol ='na'
        score = -1
        strand = "+"
        for line in FH:
            count += 1
            line_data = line.strip().split()
            enh_chr = line_data[0]
            enh_start = line_data[1]
            enh_end = line_data[2]
            region_id = f"{enh_chr}:{enh_start}-{enh_end}"
            symbol=enh_chr
            enh = RegionAnno(label + str(count),
                             region_id,
                             symbol,
                             enh_chr,
                             int(enh_start),
                             int(enh_end),
                             strand,
                             tss,
                             tss,
                             score)
            if symbol in result_dict:
                result_dict[symbol].add_enh(enh)
            else:
                result_dict[symbol] = GeneEnhSet(symbol, {enh})
    return result_dict

########
# OLD
########

#def add_gene_atlas_data(data_dict, new_data_dict):
#    for symbol, enh_list in new_data_dict.items():
#        for enh in enh_list:
#            if symbol in data_dict:
#                data_dict[symbol].add_enh(enh)
#            else:
#                data_dict[symbol] = GeneEnhSet(symbol, {enh})
#    return data_dict
 
#def format_enh_anno(enh_dict_list, enh_dict, label='enh'):
#    results = []
#    #enh_label = label
#    for gene_id in enh_dict_list:
#        count = 0
#        for enh_id in enh_dict_list[gene_id]:
#            count += 1
#            enh_label = label + "_" + str(count)
#            enh = enh_dict[enh_id]
#            result = [enh_label, enh.gene, enh.chr, enh.strand,
#                      str(enh.gene_txStart), str(enh.gene_txEnd),
#                      str(enh.start), str(enh.end)]
#            results.append(result)
#    return results      

#def assign_enh_to_gene_by_score(enh_dict, num_enh_per_gene=5, enh_score_thresh=0.):          
#    gene_enh_dict = defaultdict(set)
#    enh_scores = dict()
#    results_dict = defaultdict(list)
#    for enh in enh_dict.values():
#        if enh.connect_score >= enh_score_thresh:
#            uni_id = enh.gene + '-' + enh.id
#            enh_scores[uni_id] = enh.connect_score
#            if enh.gene in gene_enh_dict:
#                if len(gene_enh_dict[enh.gene]) >= num_enh_per_gene:
#                    min_enh_id, min_score = check_min_set_vals(gene_enh_dict[enh.gene], enh_scores, enh.gene)
#                    if enh.connect_score > min_score:
#                       gene_enh_dict[enh.gene].discard(min_enh_id)
#                       gene_enh_dict[enh.gene].add(enh.id)
#                else:
#                    add_dict_set(gene_enh_dict, enh.gene, enh.id)
#            else:
#                add_dict_set(gene_enh_dict, enh.gene, enh.id)
#    for gene_id, enh_set in gene_enh_dict.items():
#        if len(enh_set) >= num_enh_per_gene:
#            results_dict[gene_id] = sorted(enh_set, key=lambda x: enh_scores[gene_id + '-' + x])
#    return results_dict


#def check_min_set_vals(enh_set, enh_scores, gene_id):
#    min_enh_id = "na"
#    min_enh_score = -1
#    for enh_id in enh_set:
#        uni_id = gene_id + '-' + enh_id
#        if min_enh_score == -1:
#            min_enh_score = enh_scores[uni_id]
#            min_enh_id = enh_id
#        elif enh_scores[uni_id] < min_enh_score:
#            min_enh_score = enh_scores[uni_id]
#            min_enh_id = enh_id
#    return min_enh_id, min_enh_score
 
 
#def read_geneHancer_old(file, enh_score_threshold=0., connectionn_score_threshold=0.):
#    #result_enh_list = []
#    result_enh_dict = dict()
#    with open(file, 'r') as FH:
#        line_count = 0
#        for line in FH:
#            line_count += 1
#            if line_count == 1:
#                continue  # skip header
#            # chrom,source,feature name,start,end,score,strand,frame,attributes
#            line_data = line.strip().split(',')
#            enh_score = float(line_data[5])
#            if enh_score >= enh_score_threshold:
#                attribute_list = line_data[8].split(';')
#                enh_id = "na"
#                gene_id = "na"
#                connection_score = "na"
#                for attr in attribute_list:
#                    (param, val) = attr.split('=')
#                    if param == "genehancer_id":
#                        enh_id = val
#                    elif param == "connected_gene":
#                        gene_id = val
#                    elif param == "score":
#                        connection_score = float(val)
#                        if connection_score >= connectionn_score_threshold:
#                            enhancer = RegionAnno( enh_id,
#                                                   'na',
#                                                   gene_id,
#                                                   line_data[0],  # chrom
#                                                   line_data[3],  # start
#                                                   line_data[4],  # end
#                                                   'na',  # strand
#                                                   -1,  # txStart
#                                                   -1,  # txEnd 
#                                                   connect_score=connection_score
#                                                   )
#                            result_enh_dict[enh_id] = enhancer
#    return result_enh_dict

#def add_gene_info(enh_dict, gene_file):
#    gene_list = read_gene_file(gene_file)
#    gene_dict = index_gene_list(gene_list)
#    result_dict = dict()
#    total_num = 0
#    convert_num = 0
#    for enh_id, enh in enh_dict.items():
#        newEnh = ""
#        total_num += 1
#        if enh.gene in gene_dict:
#            convert_num += 1
#            new_gene = gene_dict[enh.gene]
#            newEnh = RegionAnno(enh.id,
#                                enh.region_id,
#                                new_gene.id,
#                                enh.chr,
#                                enh.start,
#                                enh.end,
#                                new_gene.strand,
#                                new_gene.txStart,
#                                new_gene.txEnd,
#                                enh.connect_score,
#                               )
#        else:
#            newEnh = enh
#        result_dict[enh_id] = newEnh
#    return result_dict, total_num, convert_num
 
 
