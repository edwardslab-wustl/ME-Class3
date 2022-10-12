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
