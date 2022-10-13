from dataclasses import dataclass

import pandas as pd

@dataclass(frozen=True)
class SamplePair:
    name: str
    tag1: str
    file1: str
    tag2: str
    file2: str

def read_sample_pair (file, num):
    num = str(num)
    df_bed = pd.read_table(file, index_col=False,
        na_values = 'NA', names = ['chrom'+num, 'start'+num, 'end'+num, 'value'+num])
    df_bed['idx_value'+num] = df_bed['chrom'+num]+'_'+df_bed['start'+num].astype(str)
    df_bed = df_bed.set_index('idx_value'+num)
    return(df_bed)

def read_sample_file (input_list_file):
    pair_list = []
    with open(input_list_file, 'r') as input_list_FH:
        for line in input_list_FH:
            if not line.startswith('#'):
                (tag1, file1, tag2, file2) = line.strip().split()
                name = "_".join((tag1, tag2))
                pair_list.append(SamplePair(name, tag1, file1, tag2, file2))
    return pair_list