import sys
import os
import re
import gzip
from collections import defaultdict

from dataclasses import dataclass

#### dataclasses ####
@dataclass(frozen=True)
class Param:
    anno_type: str
    region_size: int
    num_pts: str
    lowerBound: int = 0
    upperBound: int = 0
    locationList: str = ''
    
    def size(self) -> int:
        size = 2 * self.region_size 
        if self.lowerBound !=0 or self.upperBound !=0:
            size = self.upperBound - self.lowerBound
        return size
    
    def left_label(self) -> str:
        if self.lowerBound !=0 or self.upperBound !=0:
            left_label = self.lowerBound
        else:
            left_label = -self.region_size
        return left_label
    
    def right_label(self) -> str:
        if self.lowerBound !=0 or self.upperBound !=0:
            right_label = self.upperBound
        else:
            right_label = self.region_size
        return right_label
        

@dataclass(frozen=True)
class GeneAnno:
    id: str
    txStart: int
    txEnd: int
    cdsStart: int
    cdsEnd: int
    cdsStartStat: str
    cdsEndStat: str
    exonCount: int
    strand: str
    chr: str
    
    def gene_length(self) -> int:
        return self.txEnd - self.txStart
    
    def tss(self) -> int:
        if self.strand == '+':
            tss = self.txStart
        else:
            tss = self.txEnd
        return tss
    
    def tes(self) -> int:
        if self.strand == '+':
            tes = self.txEnd
        else:
            tes = self.txStart
        return tes
                
@dataclass(frozen=True)
class RegionAnno:
    id: str
    region_id: str
    gene: str
    chr: str
    start: int
    end: int
    strand: str
    gene_txStart: int
    gene_txEnd: int
    
    def region_length(self) -> int:
        return self.end - self.start
    
    def mid_point(self) -> int:
        return round(float(self.end - self.start) / 2)
 
#### functions ####
def format_args_to_print( args ):
    result = '##### Input Params #####\n'
    for arg in vars(args):
        if arg != 'func':
            result = result + arg + ": " + str(getattr(args, arg)) + "\n"
    result = result + '########################\n\n'
    return result

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def print_tup(inp_tup, noftf, seperation, end_chr):
        noftf.write(seperation.join(str(dat) for dat in inp_tup) + end_chr)
        
def mk_output_dir(dir):
    if not os.path.exists(dir):
        os.mkdir(dir)

def print_to_log(FH, *args, **kwargs):
    FH.write(*args, **kwargs)
    FH.flush()
    
def read_anno_file(anno_file, anno_type):
    if anno_type == 'tss':
        feature_list = read_gene_file(anno_file)
    elif anno_type == 'enh':
        feature_list = read_region_file(anno_file)
    else:
        eprint("Can't recognize anno_type: " + anno_type + "\nCheck --anno_type specification in help.")
        exit()
    return feature_list

def read_param_file_list(file):
    result = list()
    with open(file, 'r') as FH:
        for line in FH:
            if line.startswith('#'):
                (anno_type,region_size,num_pts) = line.lstrip('#').strip().split(',')
                param = Param(anno_type, int(region_size), int(num_pts))
            else:
                (key, pos) = line.strip().split(',')
                result.append(int(pos))
    return result, param

def read_params_from_interp(file):
    param_result = dict()
    data_result = dict()
    with open(file, 'r') as FH:
        for line in FH:
            if line.startswith('#'):
                (comment, info, data) = line.strip().split()
                (anno_type,region_size,num_pts) = info.split(':')[1].split(',')
                param = Param(anno_type, int(region_size), int(num_pts))
                data_list = data.split(':')[1].split(';')
                data_pos_list = [ int(x.split(',')[1]) for x in data_list if len(x) > 0]
                if anno_type in param_result:
                    eprint(f"Found multiple entries for anno_type: {anno_type} in interp file: {file}\nOnly using first one.")
                else:
                    param_result[ anno_type ] = param
                    data_result [ anno_type ] = data_pos_list
    return data_result, param_result

def read_params_from_interp_2(file):
    param_result = dict()
    data_result = dict()
    with open(file, 'r') as FH:
        for line in FH:
            if line.startswith('#'):
                (comment, info, data) = line.strip().split()
                #(anno_type,region_size,num_pts) = info.split(':')[1].split(',')
                info_list = info.split(':')[1].split(',')
                anno_type = info_list[0]
                region_size = int(info_list[1])
                num_pts = int(info_list[2])
                param = ''
                if len(info_list) > 3:
                    region_start = int(info_list[3])
                    region_end = int(info_list[4])
                    enh_list = []
                    if len(info_list) > 5:
                        for i in range(5, len(info_list)):
                            enh_list.append(info_list[i])
                        param = Param(anno_type, region_size, num_pts, region_start, region_end, ','.join(enh_list))
                    else:
                        param = Param(anno_type, region_size, num_pts, region_start, region_end)
                else:
                    param = Param(anno_type, region_size, num_pts)
                data_list = data.split(':')[1].split(';')
                data_pos_list = [ [x.split(',')[0],int(x.split(',')[1])] for x in data_list if len(x) > 0]
                if anno_type in param_result:
                    eprint(f"Found multiple entries for anno_type: {anno_type} in interp file: {file}\nOnly using first one.")
                else:
                    param_result[ anno_type ] = param
                    data_result [ anno_type ] = data_pos_list
    return data_result, param_result

def index_raw_data (file, size):
    index = defaultdict(list)
    data = defaultdict(dict)
    idx_size = size * 4 
    with open(file, 'r') as FH:
        for line in FH:
            if not line.startswith('#'):
                (chrom, pos1, pos2, val) = line.strip().split()
                pos1 = int(pos1)
                idx = int(pos1 / idx_size)
                data[chrom][pos1] = float(val)
                if (chrom, idx) in index:
                    index[(chrom,idx)].append(pos1)
                else:
                    index[(chrom,idx)] = [pos1]
    return data, index, idx_size

def read_region_file(file):
    region_list = []
    with open(file, 'r') as FH:
        for line in FH:
            if not line.startswith('#'):
                items = line.strip().split()
                id = items[0] + '-' + items[1]
                region = RegionAnno( id, items[0], items[1], items[2],
                                    int(items[6]), int(items[7]), items[3],
                                    int(items[4]), int(items[5]) )
                region_list.append(region)
    return region_list

def read_enhancer_file(file, delimiter = 1000):
    dict_bed = defaultdict(list)
    dict_bed_idx = defaultdict(dict)
    with open(file, 'r') as FH:
        line_num=0
        prior_chrom = 'this is not a chromosome name'
        for line in FH:
            if not line.startswith('#'):
                data = line.strip().split()
                chrom = data[0]
                pos1 = int(data[1])
                pos2 = int(data[2])
                if chrom in dict_bed:
                    dict_bed[chrom].append( [pos1,pos2])
                else:
                    dict_bed[chrom] = [ [pos1,pos2] ]
                if chrom != prior_chrom:
                    line_num = 0
                    prior_chrom = chrom
                    pos_track = int( pos1 / delimiter )
                    dict_bed_idx[chrom][pos_track] = 1 #line_num to handle -1 for very first line
                curr_pos_track = int( pos1 / delimiter)
                if pos_track != curr_pos_track:
                    for i in range (pos_track, curr_pos_track):
                        if pos_track + 1 == curr_pos_track:
                            dict_bed_idx[chrom][curr_pos_track] = line_num
                            pos_track = curr_pos_track
                        if pos_track + 1 < curr_pos_track:
                            dict_bed_idx[chrom][pos_track+ 1] = dict_bed_idx[chrom][pos_track]
                            pos_track += 1
                line_num += 1
    return dict_bed, dict_bed_idx
                
def read_gene_file(file):
    gene_list = []
    with open(file, 'r') as FH:
        for line in FH:
            if not line.startswith('#'):
                items = line.strip().split()
                gene_list.append( GeneAnno( items[1], #geneID
                                            int(items[4]), #txStart
                                            int(items[5]), #txEnd
                                            int(items[6]), #cdsStart
                                            int(items[7]), #cdsEnd
                                            items[13], #cdsStartStat
                                            items[14], #cdsEndStat
                                            int(items[8]), #exonCount
                                            items[3], #strand
                                            items[2] )) #chr
    return gene_list

def read_bed_file(file):
    if re.search(".gz$", file):
        with gzip.open(file, 'rt') as FH:
            bed_file_lines = FH.readlines()
    else:
        with open(file, 'r') as FH:
            bed_file_lines = FH.readlines()
    return bed_file_lines