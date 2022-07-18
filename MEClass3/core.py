#
# Authors: Manoj Kumar Singh (manoj@wustl.edu) and John R Edwards (jredwards@wustl.edu).
#
#---------------------------------------------------------------------------

import sys
import os
import argparse
import scipy
import scipy.interpolate
import scipy.ndimage.filters
import numpy as np
import pandas as pd
import resource
import gc
import math
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score
import gzip
import re

#---------------------------------------------------------------------------

class Interpolation:
    #
    def __init__(self, cpos_dat, dmet_dat, reg_start, reg_end, strand, reg_type, args):
        #
        self.cpos_raw = cpos_dat
        self.dmet_raw = dmet_dat
        self.sigma = args.sigma_inp # 50
        self.num_intrp_point = args.nip_inp #5
        self.anchor_window = args.anch_win
        self.re_flnk_len = args.refl_inp
        self.re_flnk_fetrs = args.reff_inp    
        #
        self.flankNorm = args.flankNorm
        self.fixedReFlnk = args.fixedReFlnk
        #
        self.reg = reg_type
        self.strand = strand
        self.reg_start = reg_start
        self.reg_end = reg_end
        #
        if self.reg == 'gene':
            #
            self.interp_bin = args.ibin_inp # to select region
            #
            self.tss = self.reg_start if self.strand == '+' else self.reg_end
            #
            self.interp_start = self.reg_start - self.interp_bin
            self.interp_end = self.reg_end + self.interp_bin
            #
            # if == 'TSS':
            self.output_reg_start = self.tss - self.interp_bin  
            self.output_reg_end = self.tss + self.interp_bin
            #    
        if self.reg == 're':
            #
            self.interp_bin = self.re_flnk_len # default for re
            #    
            self.interp_start = self.reg_start - self.interp_bin
            self.interp_end = self.reg_end + self.interp_bin
            #
            self.output_reg_start = self.reg_start - self.interp_bin
            self.output_reg_end = self.reg_end + self.interp_bin
            #
        #
        # this is used for flank_normalized only
        self.lower_lim = self.interp_start - self.anchor_window
        self.upper_lim = self.interp_end + self.anchor_window
        
    def flank_zscore_normalize_values(self):
        #
        self.flank_dmet = []
        for self.i, self.pos in enumerate(self.cpos_raw):
            if self.pos in range(self.lower_lim, self.interp_start+1) \
                    or self.pos in range(self.interp_end, self.upper_lim+1): 
                self.flank_dmet.append(self.dmet_raw[self.i])
            else:
                continue
        self.flank_dmet = np.array(self.flank_dmet)
        self.flank_dmet_mean = self.flank_dmet.mean()
        self.flank_dmet_std = self.flank_dmet.std()
        zscore_norm_dmet = ( self.dmet_raw - self.flank_dmet_mean ) / self.flank_dmet_std
        #
        return( zscore_norm_dmet )
        
    def endpoint_correction(self, p_end, n_end, cpos, dmet, sigma):
        #
        self.p_end = p_end
        self.n_end = n_end
        self.cpos_arr = cpos
        self.dmet_arr = dmet
        self.sigma_value = sigma
        #
        if self.cpos_arr[0] >= self.p_end + self.sigma_value:
            try:
                if self.dmet_arr[0] > 0 and self.dmet_arr[0] != self.dmet_arr[1]:
                        self.start_corr_dmet = max(0, ( self.dmet_arr[0] \
                        - abs(self.dmet_arr[0]-self.dmet_arr[1]) ))
                elif self.dmet_arr[0] < 0 and self.dmet_arr[0] != self.dmet_arr[1]:
                        self.start_corr_dmet = min(0, ( self.dmet_arr[0] \
                        + abs(self.dmet_arr[0]-self.dmet_arr[1]) ))
                elif self.dmet_arr[0] == 0 or self.dmet_arr[0] == self.dmet_arr[1]:
                    self.start_corr_dmet = self.dmet_arr[0]/2
            except IndexError:
                self.start_corr_dmet = self.dmet_arr[0]/2
        #
        self.value_correction = True
        self.counter = 0
        while self.value_correction == True:
            if self.cpos_arr[0] >= self.p_end+self.sigma_value and self.counter < 2:
                self.dmet_arr.insert(0, self.start_corr_dmet)
                self.cpos_arr.insert(0, self.cpos_arr[0]-self.sigma_value)
                self.start_corr_dmet = 0
                self.counter += 1
                continue
            if self.counter == 2:
                if self.cpos_arr[0] > self.p_end:
                    self.dmet_arr.insert(0, 0)
                    self.cpos_arr.insert(0, p_end)
                    self.value_correction = False
                    continue
            if self.cpos_arr[0] < self.p_end+self.sigma_value and self.cpos_arr[0] \
                                > self.p_end and self.counter == 1:
    #            dmet_arr.insert(0, 0)
    #            cpos_arr.insert(0, cpos_arr[0]-sigma_value)
                if self.cpos_arr[0] > self.p_end:
                    self.dmet_arr.insert(0, 0)
                    self.cpos_arr.insert(0, p_end)
                    self.value_correction = False
                    continue
            if self.cpos_arr[0] < self.p_end+self.sigma_value and self.cpos_arr[0] \
                                > self.p_end and self.counter == 0:
                if self.cpos_arr[0] > self.p_end:
                    self.dmet_arr.insert(0, self.dmet_arr[0])
                    self.cpos_arr.insert(0, self.p_end)
                    self.value_correction = False
                    continue
            if self.cpos_arr[0] == self.p_end:
                self.value_correction = False
                continue

        if self.cpos_arr[-1] <= self.n_end-self.sigma_value:
            try:
                if self.dmet_arr[-1] > 0:
                    self.end_corr_dmet = max(0, (self.dmet_arr[-1] \
                        - abs(self.dmet_arr[-1]-self.dmet_arr[-2])))
                elif self.dmet_arr[-1] < 0:
                    self.end_corr_dmet = min(0, (self.dmet_arr[-1] \
                        + abs(self.dmet_arr[-1]-self.dmet_arr[-2])))
                elif self.dmet_arr[-1] == 0 or self.dmet_arr[-1] == self.dmet_arr[-2]:
                    self.end_corr_dmet = self.dmet_arr[-1]/2
            except IndexError:
                self.end_corr_dmet = self.dmet_arr[-1]/2

        self.value_correction = True
        self.counter = 0
        while self.value_correction == True:
            if self.cpos_arr[-1] <= self.n_end-self.sigma_value and self.counter < 2:
                self.dmet_arr.append(self.end_corr_dmet)
                self.cpos_arr.append(self.cpos_arr[-1]+self.sigma_value)
                self.end_corr_dmet = 0
                self.counter += 1
                continue
            if self.counter == 2:
                if self.cpos_arr[-1] < self.n_end:
                    self.dmet_arr.append(0)
                    self.cpos_arr.append(n_end)
                    self.value_correction = False
                    continue
            if self.cpos_arr[-1] > self.n_end-self.sigma_value and self.cpos_arr[-1] \
                                < self.n_end and self.counter == 1:
    #            dmet_arr.append(0)
    #            cpos_arr.append(cpos_arr[-1]+sigma_value)
                if self.cpos_arr[-1] < self.n_end:
                    self.dmet_arr.append(0)
                    self.cpos_arr.append(n_end)
                    self.value_correction = False
                    continue
            if self.cpos_arr[-1] > self.n_end-self.sigma_value and self.cpos_arr[-1] \
                                < self.n_end and self.counter == 0:
                if self.cpos_arr[-1] < self.n_end:
                    self.dmet_arr.append(self.dmet_arr[-1])
                    self.cpos_arr.append(self.n_end)
                    self.value_correction = False
                    continue
            if self.cpos_arr[-1] == self.n_end:
                self.value_correction = False
                continue
    #
        self.cpos_dmet_arr = [self.cpos_arr, self.dmet_arr]
        return( self.cpos_dmet_arr )
        
    def interpolate_processed_data( self, cpos, dmet, x ):
        #
        self.cpos = cpos
        self.dmet = dmet
        self.x = x
        #
        if len( self.cpos ) >= 2:
            self.dmet_intrp_pregauss = scipy.interpolate.pchip_interpolate( self.cpos, self.dmet, self.x )
            self.dmet_intrp = scipy.ndimage.filters.gaussian_filter1d( self.dmet_intrp_pregauss, self.sigma )
        elif len( self.cpos ) == 1: # ?????
            self.dmet_intrp = [ float(self.dmet[0]) for item in x ]
        else:
            self.dmet_intrp = [ float(0) for item in x ]

        self.dict_cpos_dmet_intrp = dict( list( zip( self.x, self.dmet_intrp ) ) )
        #
        return( self.dict_cpos_dmet_intrp )

    def down_sample_data( self, dmet_intrp, num_point ):    

        self.pre_dwnsmpl_dmet_intrp = dmet_intrp
        self.num_intrp_point = num_point
        #
        # Downsampling. Each point represents average over a range.
                #self.pre_dwnsmpl_dmet_intrp = np.array( self.pre_dwnsmpl_dmet_intrp )
                #
        self.split_array = np.linspace( 0, len(self.pre_dwnsmpl_dmet_intrp), num=(self.num_intrp_point+1), dtype=int )
        #
        self.dwnsmpl_dmet_intrp_subarr = np.split( self.pre_dwnsmpl_dmet_intrp, self.split_array[1:] )
        #
        self.dwnsmpl_dmet_intrp = np.array( list( np.mean(item) for item in self.dwnsmpl_dmet_intrp_subarr[:-1] ) )
        #
        return( self.dwnsmpl_dmet_intrp )

    def dat_proc( self ):
        #
        if self.flankNorm:
            self.x = range( self.lower_lim, self.upper_lim+1, 1 )
            self.dmet_flnrm = self.flank_zscore_normalize_values()
            #print(self.dmet_flnrm)
            self.cpos_dmet_intrp = self.interpolate_processed_data( self.cpos_raw, self.dmet_flnrm, self.x )
        
        if not self.flankNorm:
            self.x = range( self.interp_start, self.interp_end+1, 1 )
            # Trim data beyond interp_start and interp_end 
            self.cpos_dmet_corr = self.endpoint_correction( self.interp_start, \
                        self.interp_end, self.cpos_raw, self.dmet_raw, self.sigma )
            self.cpos_dmet_intrp = self.interpolate_processed_data( \
                        self.cpos_dmet_corr[0], self.cpos_dmet_corr[1], self.x )            

        self.pre_dwnsmpl_dmet_intrp = []    
#        for pos in np.linspace(self.output_reg_start, self.output_reg_end, num=self.num_intrp_point, dtype=int):
        for pos in range(self.output_reg_start, self.output_reg_end+1): 
            self.pre_dwnsmpl_dmet_intrp.append(self.cpos_dmet_intrp[pos])
        #
        # Downsampling. Each point represents average over a range.
        if self.reg == 're' and self.fixedReFlnk:
            # This part fixes number of features from flanking regions of re.
            #
            self.pre_dwnsmpl_dmet_intrp = np.array( self.pre_dwnsmpl_dmet_intrp )
            #
            # number of features in re
            self.re_fetrs = self.num_intrp_point - ( 2*self.re_flnk_fetrs )
            #
            self.re_flnk_split_arr = np.split(self.pre_dwnsmpl_dmet_intrp, \
                     [self.re_flnk_len, (len(self.pre_dwnsmpl_dmet_intrp)-self.re_flnk_len)])
            #
            self.re_fln_p_pre_dwnsmpl_dmet_intrp = self.re_flnk_split_arr[0]
            self.re_fln_p_dwnsmpl_dmet_intrp = \
                    self.down_sample_data( self.re_fln_p_pre_dwnsmpl_dmet_intrp, self.re_flnk_fetrs )
            #
            self.re_pre_dwnsmpl_dmet_intrp = self.re_flnk_split_arr[1]
            self.re_dwnsmpl_dmet_intrp = \
                    self.down_sample_data( self.re_pre_dwnsmpl_dmet_intrp, self.re_fetrs)
            #
            self.re_fln_n_pre_dwnsmpl_dmet_intrp = self.re_flnk_split_arr[2]
            self.re_fln_n_dwnsmpl_dmet_intrp = \
                    self.down_sample_data( self.re_fln_n_pre_dwnsmpl_dmet_intrp, self.re_flnk_fetrs )
            # Merge all three nd-array
            self.dwnsmpl_dmet_intrp = np.concatenate([self.re_fln_p_dwnsmpl_dmet_intrp, \
                            self.re_dwnsmpl_dmet_intrp, self.re_fln_n_dwnsmpl_dmet_intrp])
            #
        else:
            #
            self.pre_dwnsmpl_dmet_intrp = np.array( self.pre_dwnsmpl_dmet_intrp )
            #
            self.dwnsmpl_dmet_intrp = self.down_sample_data( self.pre_dwnsmpl_dmet_intrp, self.num_intrp_point )
                
        #    
        if self.strand == '-': # reverse for reverse strand
            self.dwnsmpl_dmet_intrp = self.dwnsmpl_dmet_intrp[::-1]

        return( self.dwnsmpl_dmet_intrp ) #

#---------------------------------------------------------------------------
def print_tup(inp_tup, noftf, seperation, end_chr):
        noftf.write(seperation.join(str(dat) for dat in inp_tup) + end_chr)

#---------------------------------------------------------------------------
def exec_proc_sample(args):
    #
    ctp_fofn = args.ctp_inp
    path_to_output = args.pto_inp
    expr_file = args.expr_inp
    #    
    df_expr = pd.read_table(expr_file).set_index('gene_id')
    expr_cell_ids = list(df_expr.columns.values)
    #
    with open(ctp_fofn, 'r') as ctp_fofn_file:
        ctp_fofn_file_lines = ctp_fofn_file.readlines()
    ctp_fofn_file.close()
    #
    #ctp_fofn_file_lines_fltrd = list( item for item in ctp_fofn_file_lines if item.strip().split()[0] in expr_cell_ids )
    ctp_fofn_file_lines_fltrd = ctp_fofn_file_lines
    #
    # Wheel difference
    i = 0
    while i < len(ctp_fofn_file_lines_fltrd):
        first_info = ctp_fofn_file_lines_fltrd[i-1] 
        second_info = ctp_fofn_file_lines_fltrd[i]

        bg1 = first_info.strip().split()[1] # First .bedgraph
        bg2 = second_info.strip().split()[1] # Second .bedgraph
        tag = first_info.strip().split()[0]+'_'+second_info.strip().split()[0]

        df1_bed = pd.read_table(bg1, index_col=False, \
                na_values = 'NA', names = ['chrom1', 'start1', 'end1', 'value1'])
#        df1_bed.sort_values(by=['chrom1', 'start1']) # sort (doesn't change answer)    
        df1_bed['idx_value1'] = df1_bed['chrom1']+'_'+df1_bed['start1'].astype(str)
        df1_bed = df1_bed.set_index('idx_value1')

        df2_bed = pd.read_table(bg2, index_col=False, \
                na_values = 'NA', names = ['chrom2', 'start2', 'end2', 'value2'])
#        df2_bed.sort_values(by=['chrom2', 'start2']) # sort (doesn't change answer)
        df2_bed['idx_value2'] = df2_bed['chrom2']+'_'+df2_bed['start2'].astype(str)
        df2_bed = df2_bed.set_index('idx_value2')

        df_merged = df_merged = pd.concat( [ df1_bed, df2_bed ], axis=1, join='inner')
        df_merged['dcm'] = (df_merged['value2'] - df_merged['value1']).round(2)
#        df_merged['start1end'] = df_merged['start1'] + 1

        output_file = tag+'.bedgraph'
        cols_to_keep = ['chrom1', 'start1', 'end1', 'dcm']
        df_merged[cols_to_keep].to_csv(path_to_output+'/'+output_file, sep='\t', na_rep='NA', header=None, index=False)
        
        i += 1
        #
        df1_bed.drop(df1_bed.index, inplace=True)
        df2_bed.drop(df2_bed.index, inplace=True)
        df1_bed.drop(df1_bed.index, inplace=True)
        del df1_bed, df2_bed, df_merged

#---------------------------------------------------------------------------
def exec_filter(args):

    # Try keeping format fixed to ['gene_id', 'chrom', 'chromStart', 'chromEnd', 'strand', 'code_info', 'name', 'description']
    #
    expr_file = args.expr_inp
    annot_file = args.annot_inp
    output_annot_file_name = args.out_annot_inp 
    output_expr_file_name = args.out_expr_inp 
    annot_type = args.annot_type    
    #
    valid_chrom = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', \
             'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', \
             'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', \
             'chrX', 'chrY']
    #
    # Read expression data and set gene_id as index
    df_expr = ( pd.read_table(expr_file, index_col=False) ).set_index('gene_id')
    # Remove duplicate gene_id entries and keep first
    df_expr = df_expr.groupby(df_expr.index).first()
    # Read annotation data
    if annot_type == 'Gencode':
        df_annot = ( pd.read_table(annot_file, index_col=False, na_values = 'NA', names = \
                ['gene_id', 'chrom', 'chromStart', 'chromEnd', 'strand', 'code_info', 'name', 'description']) ).set_index('gene_id')
    if annot_type == 'RefSeq':
        df_annot = ( pd.read_table(annot_file, index_col=False, na_values = 'NA', names = \
                ['gene_id', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', \
                    'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']) ).set_index('gene_id')
    # Remove duplicate gene_id entries and keep firs
    df_annot = df_annot.groupby(df_annot.index).first()
    # Select only those for which expr daata is available
    df_filter = df_annot[(df_annot.index).isin(df_expr.index)]
    # add chr
    if annot_type == 'Gencode': 
        df_filter.loc[ :, 'chrom'] = 'chr' + df_filter['chrom'].astype(str)
    # chromosome filter
    df_filter = df_filter[(df_filter['chrom']).isin(valid_chrom)]
    # Write it to a file
    df_filter.to_csv(output_annot_file_name, sep='\t', na_rep='NA') #, header=None)
    #
    # filter and Write expression file after removing duplicates
    df_expr = df_expr[(df_expr.index).isin(df_filter.index)] 
    df_expr.to_csv(output_expr_file_name, sep='\t', na_rep='NA', float_format='%.3f') #, header=None)

#---------------------------------------------------------------------------
def exec_proximity_list(args):

    # This part of code generate proximity list of regulatory elements
    # add orientation relative to TSS, in addition to n and p.
    #
    p_sidx = args.pst_inp # 1
    p_eidx = args.ped_inp # 2
    n_sidx = args.nst_inp # 1
    n_eidx = args.ned_inp # 2
    twin = args.twin_inp
    gene_list_file = args.glf_inp
    path_to_input = os.getcwd()

    delimiter = 1000
    chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', \
            'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', \
            'chr22', 'chrX', 'chrY']
    dict_bed = {}
    for cnt in range( len(chr_list) ):
        with open(path_to_input+'/enhc_bed/'+chr_list[cnt]+'_all_cell_type_noTSS_cons.bed', 'r') as bed_file:
            bed_file_lines = bed_file.readlines()
        bed_file.close()

    #    print(chr_list[cnt])
        dict_bed[chr_list[cnt]] = []
        dict_bed[chr_list[cnt]+'_idx'] = {}
        for line_num, line in enumerate(bed_file_lines):
            dict_bed[chr_list[cnt]].append( [int(line.split()[1]), int(line.split()[2])] )
            if line_num == 0:
                pos_track = int( (int(line.split()[1]))/delimiter )
                dict_bed[chr_list[cnt]+'_idx'][int((int(line.split()[1]))/delimiter)] = 1 #line_num to handle -1 for very first line

            if pos_track != int( (int(line.split()[1]))/delimiter ):
                for i in range(pos_track, int( (int(line.split()[1]))/delimiter )):
                    if pos_track + 1 == int( (int(line.split()[1]))/delimiter ):
                        dict_bed[chr_list[cnt]+'_idx'][int((int(line.split()[1]))/delimiter)] = line_num
                        pos_track = int( (int(line.split()[1]))/delimiter )
                    if pos_track + 1 < int( (int(line.split()[1]))/delimiter ):
                        dict_bed[chr_list[cnt]+'_idx'][pos_track+1] = dict_bed[chr_list[cnt]+'_idx'][pos_track]
                        pos_track = pos_track + 1

        del bed_file_lines[:]

    # open files for proximity list of regulatory elements
    dict_flnm = {}
    for i in range(p_sidx, p_eidx+1):
        dict_flnm['p_end_'+str(i)] = open('re_p_'+str(i), 'w')
    for i in range(n_sidx, n_eidx+1):
        dict_flnm['n_end_'+str(i)] = open('re_n_'+str(i), 'w')

    # Readin gene info
    with open(gene_list_file, 'r') as gene_list_file_file:
        gene_list_file_lines = gene_list_file_file.readlines()
    gene_list_file_file.close()

    #gene_data = []
    for reg_info in gene_list_file_lines:
        reg_items = reg_info.strip().split()
        #
        if reg_items[1].startswith('N'): # Refseq
            gene_id, chrom, strand, tx_start, tx_end = \
                reg_items[1], reg_items[2], reg_items[3], int(reg_items[4]), int(reg_items[5])
        #
        if reg_items[0].startswith('ENSG'): # Ensembl
            gene_id, chrom, tx_start, tx_end = \
                reg_items[0], reg_items[1], int(reg_items[2]), int(reg_items[3])
            strand = '-' if reg_items[4] == '-1' else '+'
        #    
        if chrom not in chr_list:
            continue 
        #
        tss = tx_start if strand == '+' else tx_end
        #
        cnt_ne = 0
        for info in dict_bed[chrom][ (dict_bed[chrom+'_idx'].get( int(int(tx_start)/delimiter), 1))-1 : ]:
        #    if info[0] > int(tx_end):
            if info[0] > int( tss+twin ):
                cnt_ne += 1
                if cnt_ne in range (n_sidx, n_eidx+1):
                    if strand == '+':
                        dict_flnm['n_end_'+str(cnt_ne)].write( gene_id+'\t' \
                        +chrom+'\t'+strand+'\t'+str(tx_start)+'\t'+str(tx_end)+'\t' \
                        +str(info[0])+'\t'+str(info[1])+'\n' )
                    if strand == '-':
                        dict_flnm['p_end_'+str(cnt_ne)].write( gene_id+'\t' \
                        +chrom+'\t'+strand+'\t'+str(tx_start)+'\t'+str(tx_end)+'\t' \
                        +str(info[0])+'\t'+str(info[1])+'\n' )
                if cnt_ne > n_eidx:
                    break
        #
        # p-end
        # Get the index of first p_end reg.
        for info in dict_bed[chrom][ (dict_bed[chrom+'_idx'].get( int(int(tx_start)/delimiter), 1))-1 : ]: 
            # Notice -2 here: This may give key error -> fix it later # fixed
        #    if info[1] < int(tx_start):
            if info[1] < int( tss-twin ):
                prev_index = dict_bed[chrom].index(info)
                continue
        #    if info[1] >= int(tx_start):
            if info[1] >= int( tss-twin ):
                p_search_index = prev_index
                break
        cnt_pe = 0
        for info in reversed( dict_bed[chrom][ : p_search_index+2 ] ):
        #    if info[1] < int(tx_start):
            if info[1] < int( tss-twin ):
                cnt_pe += 1
                if cnt_pe in range(p_sidx, p_eidx+1):
                    if strand == '+':
                        dict_flnm['p_end_'+str(cnt_pe)].write( gene_id+'\t' \
                        +chrom+'\t'+strand+'\t'+str(tx_start)+'\t'+str(tx_end)+'\t' \
                        +str(info[0])+'\t'+str(info[1])+'\n' )
                    if strand == '-':
                        dict_flnm['n_end_'+str(cnt_pe)].write( gene_id+'\t' \
                        +chrom+'\t'+strand+'\t'+str(tx_start)+'\t'+str(tx_end)+'\t' \
                        +str(info[0])+'\t'+str(info[1])+'\n' )
                if cnt_pe > p_eidx:
                    break

    # close files for proximity list of regulatory elements
    for i in range(p_sidx, p_eidx+1):
        dict_flnm['p_end_'+str(i)].close()
    for i in range(n_sidx, n_eidx+1):
        dict_flnm['n_end_'+str(i)].close()
#---------------------------------------------------------------------------
def exec_interp(args):
    #
    num_intrp_point = args.nip_inp
    reg_fofn = args.rfn_inp
    sample_fofn = args.sfn_inp
    #
    interp_bin = args.ibin_inp
    #
    anchor_window = args.anch_win
    #
    # Readin region file info
    with open(reg_fofn,'r') as reg_fofn_file:
        reg_fofn_lines = reg_fofn_file.readlines()
    reg_fofn_file.close()
    #
    dict_gene, dict_gene_out = {}, {}
    dict_re, dict_re_out = {}, {}
    #
    for item in reg_fofn_lines:
        #
        if item.startswith('#'):
            continue
        #
        file_name = (item.strip().split()[1])
        file_id = (item.strip().split()[1]).strip('.gene')
        #
        with open(file_name,'r') as reg_file:
            file_lines_prefltrd = reg_file.readlines()
        reg_file.close()
        #
        # Filter Genes
        # 1. Ambiguous or incomplete TSS annotation
        # 2. Genes shorter than 5 kb.
        # 3. Genes with <40 CpGs assayed within +/-5kb of the TSS
        # 4. Genes with all CpGs within +/-5 kb of the TSS had <0.2 methylation change
        # 5. Alternative Promoters 
        # 6. cdsStartStat == cmpl; cdsEndStat == cmpl;
        # 7. differentially expressed genes >=2-Fold difference after floor of five.  
        # 
        if item.strip().split()[0] == 'gene':
            file_lines = []
            for line in file_lines_prefltrd:
                if line.startswith('#'): # header
                    continue
                #
                lines_items = line.strip().split()
                #
                txStart = int( lines_items[4] )
                txEnd = int( lines_items[5] )
                cdsStart = int( lines_items[6] )
                cdsEnd = int( lines_items[7] )
                exonCount = int( lines_items[8] )
                cdsStartStat = lines_items[13]
                cdsEndStat = lines_items[14]
                #    
                gene_length = txEnd - txStart 
                if gene_length < 5000: #Genes shorter than 5 kb.
                    continue
                #
                if (cdsStartStat != 'cmpl') or (cdsEndStat != 'cmpl'):
                    continue
                #
                file_lines.append(line)
            #
            open(file_id+'_'+args.tag_inp+'_fltrd.gene', 'w').write(''.join(file_lines))
        #
        #
        if item.strip().split()[0] == 'gene':
            dict_gene[file_id] = file_lines
#            dict_gene[file_id] = file_lines_prefltrd
            dict_gene_out[file_id] = []
            dict_gene_out[file_id].append('gene_id-sample_name') # write header
            for pos in range(num_intrp_point):
                dict_gene_out[file_id].append(',ftss_'+str(pos))
        #
        if item.strip().split()[0] == 're':
            dict_re[file_id] = file_lines_prefltrd
            dict_re_out[file_id] = []
            dict_re_out[file_id].append('gene_id-sample_name') # write header
            for pos in range(num_intrp_point):
                dict_re_out[file_id].append(',f'+''.join(file_id.split('_'))+'_'+str(pos))    
        #
#        del file_name, file_id, file_lines, file_lines_prefltrd
    #
    # Sample Information.
    with open(sample_fofn,'r') as sample_fofn_file:
        sample_fofn_lines = sample_fofn_file.readlines()
    sample_fofn_file.close()
    #
    dict_sample = {}
    if args.geneSelect: #gene_select: # for gene selection 
        sys.stdout.write('Gene Selection has been chosen. A third column in sample_fofn must have gene_ids')
        dict_sample_gene_list = {} 
    for item in sample_fofn_lines:
        if item.startswith('#'):
            continue
        dict_sample[item.strip().split()[0]] = item.strip().split()[1]
        if args.geneSelect: #gene_select:
            dict_sample_gene_list[item.strip().split()[0]] = item.strip().split()[2] #
    #
    #
    for sample_id, sample_file in dict_sample.items():
        #
        if args.geneSelect: #gene_select: #
            with open(dict_sample_gene_list[sample_id], 'r') as gene_list_file:
                gene_list_lines = gene_list_file.readlines()
            gene_list_file.close()
            #
            gene_id_list = []
            for item in gene_list_lines:
                gid = item.strip().split()[0]
                if gid.startswith('#'):
                    continue
                gene_id_list.append(gid)
            gene_id_list = set(gene_id_list)
        #
        dict_bed = {}
        sys.stdout.write(sample_id+'\n')
        #
        if re.search(".gz$", sample_file):
            with gzip.open(sample_file, 'rt') as bed_file:
                bed_file_lines = bed_file.readlines()
        else:
            with open(sample_file, 'r') as bed_file:
                bed_file_lines = bed_file.readlines()
        #
        chr_track = 'chr00'
        for line in bed_file_lines:
            if line.split()[0] != chr_track:
                chr_track = line.split()[0]
                dict_bed[chr_track] = {}
                dict_bed[chr_track][int(line.split()[1])] = float(line.split()[3])
                continue
            dict_bed[chr_track][int(line.split()[1])] = float(line.split()[3])
        del bed_file_lines[:]
        #
        # memory use
        sys.stdout.write('Memory use: '+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)+'MB\n')
        #
        for file_id, file_lines in dict_gene.items(): # gene
            #
            reg_type = 'gene'
            for reg_info in file_lines:
                reg_items = reg_info.strip().split()
                cpos_dat, cpos_tmp, dmet_dat, dmet_tmp = [], [], [], []
                #
                if reg_items[1].startswith('N'):
                    gene_id, chrom, strand, tx_start, tx_end = \
                        reg_items[1], reg_items[2], reg_items[3], int(reg_items[4]), int(reg_items[5])
                    tss = tx_start if strand == '+' else tx_end
                    tes = tx_end if strand == '+' else tx_start
                #
                if reg_items[0].startswith('ENSG'): # Ensembl
                    gene_id, chrom, tx_start, tx_end = \
                        reg_items[0], reg_items[1], int(reg_items[2]), int(reg_items[3])
                    strand = '-' if reg_items[4] == '-1' else '+'
                    tss = tx_start if strand == '+' else tx_end
                    tes = tx_end if strand == '+' else tx_start 
                #
                if args.geneSelect and ( gene_id not in gene_id_list ): # gene selection
                    continue
                #
                sys.stdout.write(gene_id+'\n')
                # ----------------------------------------------------------------------------
                # Methylation based gene filters
                # 3. Genes with <40 CpGs assayed within +/-5kb of the TSS
                # 4. Genes with all CpGs within +/-5 kb of the TSS had <0.2 methylation change
                #
                for cpos in range( (tss-interp_bin), (tss+interp_bin)+1 ):            
                    try:
                        dmet_tmp.append(dict_bed[chrom][cpos])
                        cpos_tmp.append(cpos)
                    except KeyError:
                        continue
                #                    
                if len(dmet_tmp) < args.min_gene_meth:
                    sys.stderr.write(gene_id + ' has < ' +  str(args.min_gene_meth)  + ' CpGs assayed in ' + sample_id + '\n')
                    continue
                #
                if max(dmet_tmp) < 0.2:
                    sys.stderr.write(gene_id + ' has < 0.2 maximum methylation change in ' + sample_id + '\n')
                    continue
                #    
                #-----------------------------------------------------------------------------    
                #
                if not args.flankNorm: # Endpoint correction will be used in this case
                    anchor_window = 0
                #
#                if args.flankNorm:
                for cpos in range( tx_start-(interp_bin+anchor_window), tx_end+(interp_bin+anchor_window)+1 ):
                    try:
                        dmet_dat.append(dict_bed[chrom][cpos])
                        cpos_dat.append(cpos)
                    except KeyError:
                        continue
                #
                # Missing anchor            
                if args.flankNorm and \
                    (( cpos_dat[0] not in range( tx_start-(interp_bin+anchor_window), tx_start-interp_bin ) ) or \
                    ( cpos_dat[-1] not in range( tx_end+interp_bin+1, tx_end+(interp_bin+anchor_window)+1 ) )):
                    sys.stderr.write(gene_id + ' has not CpGs in anchor windows in ' + sample_id + '\n')
                    continue 
                # 
                # Gather interpolated data
                interpolated_dmet_data = Interpolation(cpos_dat, dmet_dat, tx_start, tx_end, strand, reg_type, args).dat_proc()
                
                # cross-check number of interpolated features
                if len(interpolated_dmet_data) != num_intrp_point:
                    sys.stdout.write('Inconsistent number of interpolation features!')
                    sys.exit() # Exit if number is inconsistent 
                #
                # Write data
#                if len(interpolated_dmet_data) > 0:
                dict_gene_out[file_id].append('\n'+gene_id+'-'+sample_id+',')
                dict_gene_out[file_id].append(','.join(str(item) for item in interpolated_dmet_data))
                #
                del dmet_dat[:], cpos_dat[:], dmet_tmp[:], cpos_tmp[:]
        #        #
        # Regulatory Elements
        for file_id, file_lines in dict_re.items(): # regulatory elements
            #
            reg_type = 're'
            interp_bin = 500 # for model gene plots
            #
            for reg_info in file_lines:
                #
                reg_items = reg_info.strip().split()
                cpos_dat, cpos_tmp, dmet_dat, dmet_tmp = [], [], [], []
                #
                gene_id, chrom, strand, reg_start, reg_end = \
                        reg_items[0], reg_items[1], reg_items[2], int(reg_items[5]), int(reg_items[6])
                #
                sys.stdout.write(gene_id+'\n')
                #
                                # ----------------------------------------------------------------------------
                                # Methylation based gene filters
                                # 3. Genes with <2 CpGs assayed within +/-5kb of the TSS
                                #
                for cpos in range( (reg_start-interp_bin), (reg_end+interp_bin)+1 ):
                    try:
                        dmet_tmp.append(dict_bed[chrom][cpos])
                        cpos_tmp.append(cpos)
                    except KeyError:
                        continue
                                #
                if len(dmet_tmp) < args.min_re_meth:
                    sys.stderr.write(gene_id + ' has < ' +  str(args.min_re_meth) \
                                 + ' CpGs assayed in ' + sample_id + '\n')
                    continue
                                #
#                if max(dmet_tmp) < 0.2
#                    sys.stderr.write(gene_id + ' has < 0.2 maximum methylation change in ' + sample_id + '\n')
#                    continue
                                #
                                #-----------------------------------------------------------------------------
                #
                if not args.flankNorm:
                    anchor_window = 0
                #
#                if args.flankNorm:
                for cpos in range( reg_start-(interp_bin+anchor_window), reg_end+(interp_bin+anchor_window)+1 ):
                    try:
                        dmet_dat.append(dict_bed[chrom][cpos])
                        cpos_dat.append(cpos)
                    except KeyError:
                        continue

                # Missing anchor
                if args.flankNorm and \
                    (( cpos_dat[0] not in range( reg_start-(interp_bin+anchor_window), reg_start-interp_bin ) ) or \
                    ( cpos_dat[-1] not in range( reg_end+interp_bin+1, reg_end+(interp_bin+anchor_window)+1 ) )):
                    sys.stderr.write(gene_id + ' has not CpGs in anchor windows in ' + sample_id + '\n')
                    continue

                # Gather interpolated data
                interpolated_dmet_data = Interpolation(cpos_dat, dmet_dat, reg_start, reg_end, strand, reg_type, args).dat_proc()
                #
                # cross-check number of interpolated features
                if len(interpolated_dmet_data) != num_intrp_point:
                    sys.stdout.write(str(len(interpolated_dmet_data)))
                    sys.stdout.write('Inconsistent number of interpolation features!')
                    sys.exit()
                #
                # Write data
                #if len(interpolated_dmet_data) > 0:
                dict_re_out[file_id].append('\n'+gene_id+'-'+sample_id+',')
                dict_re_out[file_id].append(','.join(str(item) for item in interpolated_dmet_data))
                #
                del dmet_dat[:], cpos_dat[:], dmet_tmp[:], cpos_tmp[:]
                #
        dict_bed.clear()
        del dict_bed
    #
    # Write data to a file
    for file_id, out_lines in dict_gene_out.items(): 
         open(args.tag_inp+'_gene_interp.csv', 'w').write(''.join(out_lines))
    #    
    for file_id, out_lines in dict_re_out.items():
        open(file_id+'_'+args.tag_inp+'_re_interp.csv', 'w').write(''.join(out_lines))
#---------------------------------------------------------------------------
def exec_merge_features(args):
    
    # Merge features.
    p_sidx = args.pst_inp # 1
    p_eidx = args.ped_inp # 2
    n_sidx = args.nst_inp # 1
    n_eidx = args.ned_inp # 2
    tss_interp_file = args.tif_inp
    
    # df interp for tss
    df_tss = pd.read_csv(tss_interp_file).set_index('gene_id-sample_name')

    # open df for proximity list of regulatory elements and merge
    df_merged = df_tss

    # Merge n-end
    if (n_sidx != -1):
        dict_dfn = {}
        for i in range(n_sidx, n_eidx+1):
            dict_dfn['n_end_'+str(i)] = pd.read_csv('re_n_'+str(i)+'_interp.csv').set_index('gene_id-sample_name')
            dict_dfn['n_end_'+str(i)] = dict_dfn['n_end_'+str(i)][ (dict_dfn['n_end_'+str(i)].index).isin(df_tss.index) ]
            #
            df_merged = pd.concat( [ df_merged, dict_dfn['n_end_'+str(i)] ], axis=1 ) #, join='inner')
            # clear memory
            dict_dfn['n_end_'+str(i)] = dict_dfn['n_end_'+str(i)].iloc[0:0]
            del dict_dfn['n_end_'+str(i)]


    # Merge p-end
    if (p_sidx != -1):
        dict_dfp = {}
        for i in range(p_sidx, p_eidx+1):
            dict_dfp['p_end_'+str(i)] = pd.read_csv('re_p_'+str(i)+'_interp.csv').set_index('gene_id-sample_name')
            dict_dfp['p_end_'+str(i)] = dict_dfp['p_end_'+str(i)][ (dict_dfp['p_end_'+str(i)].index).isin(df_tss.index) ]
            #
            df_merged = pd.concat( [ dict_dfp['p_end_'+str(i)], df_merged ], axis=1 ) #, join='inner')
            # clear memory
            dict_dfp['p_end_'+str(i)] = dict_dfp['p_end_'+str(i)].iloc[0:0]
            del dict_dfp['p_end_'+str(i)]

    #fill all NaN to 0.0 and assign name to index
    df_merged = df_merged.fillna(0)
    df_merged.index.name = 'gene_id-sample_name'

    #Write output to a file.
    df_merged.to_csv('interp_merged.csv', sep=',')

#---------------------------------------------------------------------------
def exec_setup_expr(args):
    #
    interp_file = args.intrp_inp
    expr_file = args.expr_inp
    expr_floor_value = args.efv_inp
    #
    # Read in interpolation csv file.
    df_interp =  pd.read_csv(interp_file)
    # Add sample and gene_id column
    df_interp['gene_id'] = df_interp['gene_id-sample_name'].apply(lambda x: x.split('-')[0])
    df_interp['sample_name'] = df_interp['gene_id-sample_name'].apply(lambda x: x.split('-')[1])    
    # List of sample names
    sample_names =  list( set( df_interp['sample_name'].tolist() ) )
    cell_type_names = []
    for item in sample_names:
        cell_type_names.extend( item.strip().split('_') )
    #
    # Read expression data and set gene_id as index
    df_expr = ( pd.read_table(expr_file, index_col=False) ).set_index('gene_id')
    # Remove Duplicates
#    df_expr = df_expr.groupby(df_expr.index).first()
    # Floor expression values for selected cell types
    if args.floor_expr:    
        for cell_type in cell_type_names:
            df_expr.loc[ df_expr[cell_type] < expr_floor_value, cell_type ] = expr_floor_value
    #
    # Fold change in expression value [1]. [2].  [3]. log2(floored_Expr_B/floored_Expr_A)
    for sample in sample_names:
        if args.dexpr_flag == 0:
            df_expr[sample] = np.where( df_expr[ sample.split('_')[0] ] > df_expr[ sample.split('_')[1] ], \
                        -( df_expr[ sample.split('_')[0] ]/df_expr[ sample.split('_')[1] ] ), \
                        ( df_expr[ sample.split('_')[1] ]/df_expr[ sample.split('_')[0] ] ) )
            df_expr[sample] = np.where( df_expr[sample] == 1.0, 0.0, df_expr[sample] )
#            df_expr[sample] = np.where( np.logical_or((df_expr[sample] == 1.0), \
#                                (df_expr[sample] == -1.0)), 0.0, df_expr[sample] )
        #
        if args.dexpr_flag == 1: # 
            df_expr[sample] = np.where( df_expr[ sample.split('_')[0] ] > df_expr[ sample.split('_')[1] ], \
                        -( df_expr[ sample.split('_')[0] ]/df_expr[ sample.split('_')[1] ] ), \
                        ( df_expr[ sample.split('_')[1] ]/df_expr[ sample.split('_')[0] ] ) )
        #
        if args.dexpr_flag == 2:
            df_expr[sample] = np.log2( df_expr[ sample.split('_')[1] ]/df_expr[ sample.split('_')[0] ] )
            
        if args.dexpr_flag == 3: # 
            df_expr[sample] = np.where( df_expr[ sample.split('_')[0] ] > df_expr[ sample.split('_')[1] ], \
                        -( df_expr[ sample.split('_')[0] ]/df_expr[ sample.split('_')[1] ] ), \
                        ( df_expr[ sample.split('_')[1] ]/df_expr[ sample.split('_')[0] ] ) )
    #
    # processing of interf df
    #
    df_interp = df_interp.set_index('gene_id-sample_name')
    # Remove Duplicates
#    df_interp = df_interp.groupby(df_interp.index).first()
    #
    df_interp['expr_value'] = 0.0  # Float type
    df_interp['expr_flag'] = 0      # int type
    # 
    for gene_item in df_interp.index:
    #    print(gene_item)
        try:
            df_interp.loc[gene_item, 'expr_value'] = df_expr.loc[ gene_item.split('-')[0], gene_item.split('-')[1] ]
        except KeyError:
            sys.stdout.write(gene_item.split('-')[0]+' not found in expression file.\n')
            df_interp.loc[gene_item, 'expr_value'] = 0
        #
    if args.dexpr_flag == 0: # me-class demo
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] > 0, 1, df_interp['expr_flag'] )
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] < 0, -1, df_interp['expr_flag'] )
        df_interp['expr_flag'] =  df_interp['expr_flag'].astype(int)
    #
    if args.dexpr_flag == 1: # me-class paper
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] >= 2, 1, df_interp['expr_flag'] )
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] <= -2, -1, df_interp['expr_flag'] )
        df_interp['expr_flag'] =  df_interp['expr_flag'].astype(int)
    #
    if args.dexpr_flag == 2:
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] > 0.0, 1, df_interp['expr_flag'] )
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] < 0.0, -1, df_interp['expr_flag'] )
        df_interp['expr_flag'] = df_interp['expr_flag'].astype(int)
        
    if args.dexpr_flag == 3: # custom
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] >= args.expr_cutoff, 1, df_interp['expr_flag'] )
        df_interp['expr_flag'] = np.where( df_interp['expr_value'] <= -args.expr_cutoff, -1, df_interp['expr_flag'] )
        df_interp['expr_flag'] =  df_interp['expr_flag'].astype(int)
    #
    #
    # Finally write data for classifier
    df_interp.to_csv('interp_expr_data.csv', sep=',')
#---------------------------------------------------------------------------
def exec_run_clf(args):
    #
    # Load data csv file
    df = ( pd.read_csv(args.dfi_inp) ).set_index('gene_id-sample_name')
    #
    df = df[~df.isin([np.nan, np.inf, -np.inf]).any(1)]
    # List of sample names
    sample_names =  list( set( df['sample_name'].tolist() ) )
    #
    #Output info. This will store classifier performance
    df['prob_dw'] = 0.0
    df['prob_up'] = 0.0
    df['expr_pred'] = np.nan
    # Add train and test flag for classifier to use and mark all null
    df['clf_flag'] = 'null'
    #
    # Feature importance
    if args.fsl_inp == 1: # TSS 
        feature_column_names = ( df.columns[pd.Series(df.columns).str.startswith('ftss')] ).tolist()
        #
    if args.fsl_inp == 2: # TSS + RE
        feature_column_names = ( df.columns[pd.Series(df.columns).str.startswith('f')] ).tolist()
        #
    if args.fsl_inp == 3: #RE Only
        feature_column_names = ( df.columns[pd.Series(df.columns).str.startswith('fre')] ).tolist()
        #
    df_fi = pd.DataFrame(columns=feature_column_names)
    #===============================================
    if args.ss:
        #
        if args.gnorm:
            # Normalize up and down
            expr_flag_value_count = df['expr_flag'].value_counts()
            num_up = pd.Series(expr_flag_value_count)[1]
            num_dn = pd.Series(expr_flag_value_count)[-1]
            #
            if num_up > num_dn:
                idx = df.index[ (df['expr_flag'] == 1) ]
                df.at[np.random.choice(idx, size=(num_up-num_dn), replace=False), 'expr_flag'] = 0
                del idx
            #
            elif num_dn > num_up:
                idx = df.index[ (df['expr_flag'] == -1) ]
                df.at[np.random.choice(idx, size=(num_dn-num_up), replace=False), 'expr_flag'] = 0
                del idx
  
        gene_id_list = np.asarray( list(set( df['gene_id'].tolist() ) ))
        #
        kf = KFold(n_splits=10 ) # 10 fold #shuffle=True in me-class
        kf.get_n_splits(gene_id_list)
        for train_index, test_index in kf.split(gene_id_list): # K-Fold gene split
            #
            I_train, I_test = gene_id_list[train_index], gene_id_list[test_index]
            #
            df_kf = df.copy()
            #
            idx = df_kf.index[ ~df_kf['gene_id'].isin(I_test) ]
            df_kf.at[idx, 'clf_flag'] = 'train'
            del idx
            #
            idx = df_kf.index[ df_kf['gene_id'].isin(I_test) ]
            df_kf.at[idx, 'clf_flag'] = 'test'
            del idx
            # Select traing data set
            df_kf_train = df_kf[ df_kf['clf_flag']=='train' ]
            df_kf_test = df_kf[ df_kf['clf_flag']=='test' ]
            #
            features = df_kf_train.columns[pd.Series(df_kf_train.columns).str.startswith('f')]
            #
            y_train = df_kf_train['expr_flag'].values
            #
            #Setup classifier and run it
            clf = RandomForestClassifier(n_estimators=args.ntr_inp, n_jobs=args.npr_inp) #
            clf.fit(df_kf_train[features], y_train)
            #
            y_test_prob = clf.predict_proba(df_kf_test[features])
            y_test_pred = clf.predict(df_kf_test[features])
            # store performance values in df_out
            for i,idx in enumerate(df_kf_test[features].index):
                    df.at[idx, 'expr_pred'] = y_test_pred[i]
                    df.at[idx, 'prob_dw'] = y_test_prob[i][0]
                    df.at[idx, 'prob_up'] = y_test_prob[i][1]
            # Feature importance
            df_fi = pd.concat( [df_fi, (pd.DataFrame( [(clf.feature_importances_)],  columns=feature_column_names))] )
            #
            df_kf.drop(df_kf.index, inplace=True)
            df_kf_train.drop(df_kf_train.index, inplace=True)
            df_kf_test.drop(df_kf_test.index, inplace=True)
            #
            del df_kf, df_kf_train, df_kf_test, clf
            #
    #===============================================
    if not args.ss: 
        # loso loop
        for sample in sample_names:
            #
            df_loso = df.copy()
            # test sample
            sample_name_test_list = []
            sample_name_test_list.append(sample)
            sys.stdout.write(", ".join(sample_name_test_list)+'\n')
            #
            # train samples
            sample_name_train_list = []
            for item in sample_names:
                if any( i in item.split('_') for i in sample.split('_') ): # any(i in b for i in a)
                    continue
                sample_name_train_list.append(item)
            sys.stdout.write(", ".join(sample_name_train_list)+'\n')
            #
            #----------------------------------------------------------------------
            if args.gnorm:
                # Count number of up and down for test set
                test_expr_flag_value_count = df_loso[ df_loso['sample_name'].isin(sample_name_test_list) ]['expr_flag'].value_counts()
                num_up_test = pd.Series(test_expr_flag_value_count)[1]
                num_dn_test = pd.Series(test_expr_flag_value_count)[-1]
                sys.stdout.write('# of up regulated in test set = '+str(num_up_test)+'\n')
                sys.stdout.write('# of down regulated in test set = '+str(num_dn_test)+'\n\n')
                #
                # Normalize up and down regulated for train samples
                train_expr_flag_value_count = df_loso[ df_loso['sample_name'].isin(sample_name_train_list) ]['expr_flag'].value_counts()
                #    
                try:
                    num_up = pd.Series(train_expr_flag_value_count)[1]
                except KeyError:
                    num_up = 0
                sys.stdout.write('# of up regulated in train set = '+str(num_up)+'\n')
                #
                try:
                    num_dn = pd.Series(train_expr_flag_value_count)[-1]
                except KeyError:
                    num_dn = 0
                sys.stdout.write('# of down regulated in train set= '+str(num_dn)+'\n\n')
                #
                if num_up > num_dn:
                    idx = df_loso.index[ (df_loso['expr_flag'] == 1) & df_loso['sample_name'].isin(sample_name_train_list) ]
                    df_loso.at[np.random.choice(idx, size=(num_up-num_dn), replace=False), 'expr_flag'] = 0
                #    df_loso.set_value(np.random.choice(idx, size=(num_up-num_dn), replace=False), 'expr_flag', 0)
                    del idx
                #
                elif num_dn > num_up:
                    idx = df_loso.index[ (df_loso['expr_flag'] == -1) & df_loso['sample_name'].isin(sample_name_train_list) ]
                    df_loso.at[np.random.choice(idx, size=(num_dn-num_up), replace=False), 'expr_flag'] = 0
                #    df_loso.set_value(np.random.choice(idx, size=(num_dn-num_up), replace=False), 'expr_flag', 0)
                    del idx            
                #
                # Count again after normalization
                train_expr_flag_value_count_norm = df_loso[ df_loso['sample_name'].isin(sample_name_train_list) ]['expr_flag'].value_counts()
                num_up_norm = pd.Series(train_expr_flag_value_count_norm)[1]
                num_dn_norm = pd.Series(train_expr_flag_value_count_norm)[-1]
                tot_items_norm = num_up_norm + num_dn_norm
                sys.stdout.write('# of up regulated after norm = '+str(num_up_norm)+'\n')
                sys.stdout.write('# of down regulated after norm = '+str(num_dn_norm)+'\n')
                sys.stdout.write('# of Total items in train set after norm = '+str(tot_items_norm)+'\n\n')
                #
            #---------------------------------------------------------------------------
            ## Run K-Fold on gene level ##
            #
            gene_id_list = np.asarray( list( set( df_loso[ df_loso['sample_name'].isin(sample_name_test_list) ]['gene_id'].tolist() ) ) ) # me-class
    #        gene_id_list = list(set( df_loso[ df_loso['sample_name'].isin(sample_name_test_list) ]['gene_id'].tolist() )) # me-class
    #        gene_id_list = np.asarray( set( df_loso[ ((df_kf.expr_flag==-1) | (df_kf.expr_flag==1)) \
    #                & df_loso['sample_name'].isin(sample_name_train_list) ]['gene_id'].tolist() ) ) # 
    #        gene_id_list = np.asarray( set( df_loso['gene_id'].tolist() ) )
            kf = KFold(n_splits=10, shuffle=args.suf_inp) # 10 fold #shuffle=True in me-class
            kf.get_n_splits(gene_id_list)
            #
            kf_run_idx = 0    
            for train_index, test_index in kf.split(gene_id_list): # K-Fold gene split
                #
                kf_run_idx += 1 # for output track
                #            
                df_kf = df_loso.copy()
                #
                I_train, I_test = gene_id_list[train_index], gene_id_list[test_index]
                #
                # Mark train set 
                idx = df_kf.index[ ((df_kf.expr_flag==-1) | (df_kf.expr_flag==1)) & \
                        df_kf['sample_name'].isin(sample_name_train_list) & ~df_kf['gene_id'].isin(I_test) ]
                df_kf.at[idx, 'clf_flag'] = 'train'
                #df_kf.set_value(idx, 'clf_flag', 'train')
                del idx
                #    
                # Mark test set
                idx = df_kf.index[ ((df_kf.expr_flag==-1) | (df_kf.expr_flag==1)) & \
                        df_kf['sample_name'].isin(sample_name_test_list) & df_kf['gene_id'].isin(I_test) ]
                df_kf.at[idx, 'clf_flag'] = 'test'
                #df_kf.set_value(idx, 'clf_flag', 'test')
                del idx
                #
                # Select traing data set
                df_kf_train = df_kf[ df_kf['clf_flag']=='train' ]
                df_kf_test = df_kf[ df_kf['clf_flag']=='test' ]
                #
                # Select features
                if args.fsl_inp == 1: # TSS
                    features = df_kf_train.columns[pd.Series(df_kf_train.columns).str.startswith('ftss')]
                if args.fsl_inp == 2: # TSS + RE
                    features = df_kf_train.columns[pd.Series(df_kf_train.columns).str.startswith('f')]
                if args.fsl_inp == 3: # RE
                    features = df_kf_train.columns[pd.Series(df_kf_train.columns).str.startswith('fre')]
                #
                y_train = df_kf_train['expr_flag'].values
                #
                #Setup classifier and run it
                clf = RandomForestClassifier(n_estimators=args.ntr_inp, n_jobs=args.npr_inp) # me-class 1001 default
                clf.fit(df_kf_train[features], y_train)
                #
                #Performance Evaluation
                y_test_prob = clf.predict_proba(df_kf_test[features])
                y_test_pred = clf.predict(df_kf_test[features])
                # store performance values in df_out
                for i,idx in enumerate(df_kf_test[features].index):
                    df.at[idx, 'expr_pred'] = y_test_pred[i]
                    df.at[idx, 'prob_dw'] = y_test_prob[i][0]
                    df.at[idx, 'prob_up'] = y_test_prob[i][1]

                    # df.set_value(idx, 'expr_pred', y_test_pred[i]) 
                    # df.set_value(idx, 'prob_dw', y_test_prob[i][0]) 
                    # df.set_value(idx, 'prob_up', y_test_prob[i][1]) 
                    
                #
                # Feature importance.
                #print(clf.feature_importances_)
                df_fi = pd.concat( [df_fi, (pd.DataFrame( [(clf.feature_importances_)],  columns=feature_column_names))] )
                #
                sys.stdout.write('Memory use: '+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)+'MB\n')
                #
                # Clear Memory
                
                df_kf.drop(df_kf.index, inplace=True)
                df_kf_train.drop(df_kf_train.index, inplace=True)
                df_kf_test.drop(df_kf_test.index, inplace=True)
                #
                del df_kf, df_kf_train, df_kf_test, clf
                #
            df_loso.drop(df_loso.index, inplace=True)
            del df_loso
    #============================
    # Final dataframe without feature column
    df_out = df.copy()
    df_out = df_out[ (df_out['expr_pred'].isnull()==False) ]
#    df_out = df_out[ np.isnan(df_out['expr_pred']) == False ]    
    df_out = df_out.drop(df_out.columns[pd.Series(df_out.columns).str.startswith('f')] , axis=1)
    df_out.to_csv(args.tag_inp+'.RandomForestClassifier.csv', sep=',')
    #
    # Feature importance output
    df_fi.to_csv(args.tag_inp+'.RF_fi.csv', sep=',')
    dfi_items = []
    for item in df_fi.columns:
        dfi_items.append( str(df_fi[item].sum())+'\t'+str(df_fi[item].mean())+'\n' )
    #
    open(args.tag_inp+'_fi_out.txt', 'w').write(''.join(dfi_items))
    
#---------------------------------------------------------------------------
def exec_eval(args):

    df = pd.read_csv( args.dfi_inp, index_col=[0] ).reset_index()
    steps = args.steps_inp
    tag = ((args.dfi_inp).strip().split('.'))[0]
    #
    output_file_name = ((args.dfi_inp).strip('.csv'))+'_acc_rej.txt'
    #
    df['expr_pred'] = df['expr_pred'].astype(int)
    totalGenes = 0
    totalGenes_P = 0
    totalGenes_N = 0
    #
    totalGenes = df['gene_id-sample_name'].count()
    #
    expr_flag_value_count = df['expr_flag'].value_counts()
    try:
            totalGenes_P = pd.Series(expr_flag_value_count)[1]
    except KeyError:
            totalGenes_P = 0
    try:
            totalGenes_N = pd.Series(expr_flag_value_count)[-1]
    except KeyError:
            totalGenes_N = 0
    #
    output_file = open(output_file_name,'w')
    #
    sample_names = list( set( df['sample_name'].tolist() ) )
    init_gene_values = [0] * len( sample_names )
    sample_freq_90 = sample_freq_85 = sample_freq_80 = sample_freq_75  = pd.Series( init_gene_values, index=sample_names ) 
    gene_count_90_per_acc, gene_count_85_per_acc, gene_count_80_per_acc, gene_count_75_per_acc = 0, 0, 0, 0
    #
    for threshold in np.linspace(0,0.5,steps):
        #
        TP_P = 0
        TP_N = 0
        numGenes_P = 0
        numGenes_N = 0
        #
        df_tmp = df[ (df['prob_up'] >= 1-threshold) | (df['prob_dw'] >= 1-threshold) ]
        #
        #
        df_tmp_P = df_tmp[ (df_tmp['prob_up'] >= 1-threshold) ]
        numGenes_P = df_tmp_P['gene_id-sample_name'].count()
        df_tmp_P = df_tmp_P[ (df_tmp_P['expr_flag'] == 1) & (df_tmp_P['expr_pred'] == 1) ]
        TP_P = df_tmp_P['gene_id-sample_name'].count()
        #
        df_tmp_N = df_tmp[ (df_tmp['prob_dw'] >= 1-threshold) ]
        numGenes_N = df_tmp_N['gene_id-sample_name'].count()
        df_tmp_N = df_tmp_N[ ( df_tmp_N['expr_flag'] == -1 ) & (df_tmp_N['expr_pred'] == -1) ]
        TP_N = df_tmp_N['gene_id-sample_name'].count()
        #sys.stdout.write( str(threshold)+'\t'+str(numGenes_P)+'\t'+str(TP_P)+'\t'+str(numGenes_N)+'\t'+str(TP_N)+'\n')
        #
        accuracy = 0
        rejectRate = 0
        #
        if (numGenes_P + numGenes_N) > 0:
                accuracy = float(TP_P + TP_N) / float(numGenes_P + numGenes_N)
        if totalGenes > 0:
                rejectRate = float(numGenes_P + numGenes_N) / float(totalGenes)
        #
        if ( accuracy >= 0.90) :
                gene_count_90_per_acc = ( numGenes_P + numGenes_N )
                sample_freq_90 = df_tmp['sample_name'].value_counts()
        if ( accuracy >= 0.85) :
                gene_count_85_per_acc = ( numGenes_P + numGenes_N )
                sample_freq_85 = df_tmp['sample_name'].value_counts()
        if ( accuracy >= 0.80) :
                gene_count_80_per_acc = ( numGenes_P + numGenes_N )
                sample_freq_80 = df_tmp['sample_name'].value_counts()
        if ( accuracy >= 0.75) :
                gene_count_75_per_acc = ( numGenes_P + numGenes_N )
                sample_freq_75 = df_tmp['sample_name'].value_counts()
        #
        df_tmp = df_tmp.iloc[0:0]
        df_tmp_P = df_tmp_P.iloc[0:0]
        df_tmp_N = df_tmp_N.iloc[0:0]
        #
        output_file.write(str(threshold)+'\t'+str(rejectRate)+'\t'+str(accuracy)+'\t'+str( numGenes_P + numGenes_N )+'\n')

    output_file.close()
    #
    y_true = df['expr_flag'].tolist()
    y_pred = df['expr_pred'].tolist()
    #
    sys.stdout.write( 'F1 Score (macro):\t'+str( (f1_score(y_true, y_pred, average='macro') ))+'\n' )
    sys.stdout.write( 'F1 Score (micro):\t'+str( (f1_score(y_true, y_pred, average='micro') ))+'\n' )
    sys.stdout.write( 'Accuracy score:\t'+str( accuracy_score(y_true, y_pred) )+'\n' )
    #
    sys.stdout.write( '#Genes with >90% accuracy:\t'+str( gene_count_90_per_acc )+'\n' )
    with open( tag+'_genes_with_90_percent_accuracy.txt', 'w') as output_90per:
        for sample in sample_names:
            try:
                output_90per.write( str(pd.Series( sample_freq_90 )[sample]) +'\t'+ sample +'\n' )
            except KeyError:
                output_90per.write( str(0) +'\t'+ sample +'\n' )
    #
    sys.stdout.write( '#Genes with >85% accuracy:\t'+str( gene_count_85_per_acc )+'\n' )
    with open( tag+'_genes_with_85_percent_accuracy.txt', 'w') as output_85per:
        for sample in sample_names:
            try:
                output_85per.write( str(pd.Series( sample_freq_85 )[sample]) +'\t'+ sample +'\n' )
            except KeyError:
                output_85per.write( str(0) +'\t'+ sample +'\n' )
    #
    sys.stdout.write( '#Genes with >80% accuracy:\t'+str( gene_count_80_per_acc )+'\n' )
    with open( tag+'_genes_with_80_percent_accuracy.txt', 'w') as output_80per:
        for sample in sample_names:
            try:
                output_80per.write( str(pd.Series( sample_freq_80 )[sample]) +'\t'+ sample +'\n' )
            except KeyError:
                output_80per.write( str(0) +'\t'+ sample +'\n' )
    #
    sys.stdout.write( '#Genes with >75% accuracy:\t'+str( gene_count_75_per_acc )+'\n' )
    with open( tag+'_genes_with_75_percent_accuracy.txt', 'w') as output_75per:
        for sample in sample_names:
            try:
                output_75per.write( str(pd.Series( sample_freq_75 )[sample]) +'\t'+ sample +'\n' )
            except KeyError:
                output_75per.write( str(0) +'\t'+ sample +'\n' )

    sys.stdout.write( '#Total_Genes:\t'+str( totalGenes )+'\n' )

#---------------------------------------------------------------------------
def exec_proc_reg(args):

    gap_join = args.gapj_inp # 500
    max_size = args.max_inp # 10000
    tss_win = args.twin_inp # 2000
    input_gene_file = args.glf_inp
    #
    path_to_output = os.getcwd()
    if not os.path.isdir('enhc_bed'):
        os.makedirs('enhc_bed')
    #
    # Store TSS information.
    with open(input_gene_file, 'r') as refgene_file:
        refgene_file_lines = refgene_file.readlines()
    refgene_file.close()

    dict_refgene = {}
    for line in refgene_file_lines:
        refgene_items = line.strip().split()
        #
        if not refgene_items[1].startswith('NM_'):
            continue
        #
        chrom, strand, tx_start, tx_end = \
            refgene_items[2], refgene_items[3], int(refgene_items[4]), int(refgene_items[5])
        tss = tx_start if strand == '+' else tx_end
        #
        try:
            dict_refgene[chrom].append(tss-tss_win)
        except KeyError:
            dict_refgene[chrom] = []
            dict_refgene[chrom].append(tss-tss_win)
        #
    for chrom_name, chrom_items in dict_refgene.items():
        #
        dict_refgene[chrom_name] = sorted( list(set(dict_refgene[chrom_name])) , key=int )
    #
    #print(len(list(set(dict_refgene['chr1']))))
    #sys.exit()
    #
    enhc_fofn = 'enhc_list' #args.efn_inp
    with open(enhc_fofn, 'r') as enhc_fofn_file:
        enhc_fofn_file_lines = enhc_fofn_file.readlines()
    enhc_fofn_file.close()
    #
    # Assign file names to a dictionary
    df_type_name = {}
    for item in enhc_fofn_file_lines:
        if item.startswith('#'):
            continue
        df_type_name[item.strip().split()[0]] = pd.read_table(item.strip().split()[1], \
                 compression='gzip', names = ['chrom', 'chromStart', 'chromEnd', \
                'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', \
                'blockCount', 'blockSizes', 'blockStarts'])
        df_type_name[item.strip().split()[0]]['ctp'] = item.strip().split()[0]
    #
    df_chrom = {}
    df_ctp_chrom = {}
    df_chrom_cons = {}
    df_chrom_notss = {}
    #
    chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', \
                    '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
    for i in chr_list:
        chrom = 'chr'+i

        df_chrom[chrom] = pd.DataFrame(columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', \
            'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts', 'ctp'])

        df_chrom[chrom][['chromStart', 'chromEnd', 'thickStart', 'thickEnd']] = \
            df_chrom[chrom][['chromStart', 'chromEnd', 'thickStart', 'thickEnd']].astype(int) # change datatype to int

        # extract chromosome specific regions
        for item in enhc_fofn_file_lines:
            if item.startswith('#'):
                continue
            df_chrom[chrom] = df_chrom[chrom].append( \
                    df_type_name[item.strip().split()[0]].loc[ df_type_name[item.strip().split()[0]].chrom == chrom, : ] )

        # Short, reset index, reemove duplicates and print bed file
        df_chrom[chrom] = df_chrom[chrom].sort_values('thickStart')
        df_chrom[chrom] = df_chrom[chrom].reset_index(drop=True)
        df_ctp_chrom[chrom] = df_chrom[chrom].groupby(['chrom', 'thickStart', 'thickEnd'])['ctp'].apply(list).to_frame().reset_index()#.rename(columns={0:'count'})
        df_chrom[chrom] = df_chrom[chrom].groupby(['chrom', 'thickStart', 'thickEnd']).size().reset_index().rename(columns={0:'count'})
        df_chrom[chrom] = df_chrom[chrom].merge(df_ctp_chrom[chrom])
        df_chrom[chrom]['ctp'] = df_chrom[chrom]['ctp'].apply(lambda x: ','.join(x))
        df_chrom[chrom]['count'] = df_chrom[chrom]['count'].astype(int)
        df_chrom[chrom].to_csv(path_to_output+'/enhc_bed/'+chrom+'_all_cell_type.bed', sep='\t', index=False)
        #
        # Remove those falling in TSS region
        df_chrom_notss[chrom] = pd.DataFrame(columns = ['chrom', 'thickStart', 'thickEnd', 'count', 'ctp'])
        df_chrom_notss[chrom][['thickStart', 'thickEnd', 'count']] = \
                df_chrom_notss[chrom][['thickStart', 'thickEnd', 'count']].astype(int)
        #
        write_flag = 1
        for item in df_chrom[chrom].index:
            #
            if df_chrom[chrom].loc[item, 'thickStart'] > df_chrom[chrom].loc[item, 'thickEnd']:
                sys.stdout.write('thickStart > thickEnd \n')
                sys.exit()
            #
            for tss_item in dict_refgene[chrom]: #
                #
                if df_chrom[chrom].loc[item, 'thickStart'] > tss_item+(2*tss_win):
                    continue
                #
                if ( tss_item <= df_chrom[chrom].loc[item, 'thickStart'] <= tss_item+(2*tss_win) ) or \
                    ( tss_item <= df_chrom[chrom].loc[item, 'thickEnd'] <= tss_item+(2*tss_win) ) or \
                    ( df_chrom[chrom].loc[item, 'thickStart'] <= tss_item <=  df_chrom[chrom].loc[item, 'thickEnd']):
                    write_flag = 0
                    break
                #
                if tss_item > df_chrom[chrom].loc[item, 'thickEnd']:
                    break
            #
            if write_flag == 1:
                df_chrom_notss[chrom].loc[item] = df_chrom[chrom].loc[item, :]
            #
            write_flag = 1
        #
        df_chrom_notss[chrom] = df_chrom_notss[chrom].reset_index(drop=True)
        df_chrom_notss[chrom].to_csv(path_to_output+'/enhc_bed/'+chrom+'_all_cell_type_noTSS.bed', sep='\t', index=False)

        #Condolidate overlapping
        df_chrom_cons[chrom] = pd.DataFrame(columns = ['chrom', 'thickStart', 'thickEnd', 'count', 'ctp'])
        df_chrom_cons[chrom][['thickStart', 'thickEnd', 'count']] = \
                df_chrom_cons[chrom][['thickStart', 'thickEnd', 'count']].astype(int)
        prev_start = df_chrom_notss[chrom].loc[0, 'thickStart']
        prev_end = df_chrom_notss[chrom].loc[0, 'thickEnd']
        prev_count = 0 # df_chrom[chrom].loc[0, 'count']
        prev_ctp = '' #list(df_chrom[chrom].loc[0, 'ctp'])
        #
        for item in df_chrom_notss[chrom].index:
            if df_chrom_notss[chrom].loc[item, 'thickStart'] <= prev_end + gap_join:
                prev_count = prev_count + df_chrom_notss[chrom].loc[item, 'count']
                prev_ctp = prev_ctp +','+ df_chrom_notss[chrom].loc[item, 'ctp']
                if df_chrom_notss[chrom].loc[item, 'thickEnd'] > prev_end:
                    prev_end = df_chrom_notss[chrom].loc[item, 'thickEnd']
            if df_chrom_notss[chrom].loc[item, 'thickStart'] > prev_end + gap_join:
                if prev_end-prev_start <= max_size:
                    df_chrom_cons[chrom].loc[item] = [chrom, prev_start, prev_end, prev_count, prev_ctp]
                prev_start = df_chrom_notss[chrom].loc[item, 'thickStart']
                prev_end = df_chrom_notss[chrom].loc[item, 'thickEnd']
                prev_count = df_chrom_notss[chrom].loc[item, 'count']
                prev_ctp = df_chrom_notss[chrom].loc[item, 'ctp']

        df_chrom_cons[chrom].to_csv(path_to_output+'/enhc_bed/'+chrom+'_all_cell_type_noTSS_cons.bed', sep='\t', index=False, header=None)

#---------------------------------------------------------------------------

def main():
    FUNCTION_MAP = {
            'proc_sample'       : exec_proc_sample,
            'filter'            : exec_filter,
            'proc_reg'          : exec_proc_reg,
            'proximity_list'    : exec_proximity_list,
            'interp'            : exec_interp,
            'merge_features'    : exec_merge_features,
            'setup_expr'        : exec_setup_expr,
            'run_clf'           : exec_run_clf,
            'eval'              : exec_eval,
            }

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')

    sp_filter = subparsers.add_parser('filter', help='Filter argument')
    sp_filter.add_argument('-expr', action='store', dest='expr_inp', help='Name of first file')
    sp_filter.add_argument('-annot', action='store', dest='annot_inp', help='Name of second file')
    sp_filter.add_argument('-oex', action='store', dest='out_expr_inp', help='Name of output expr file')
    sp_filter.add_argument('-oan', action='store', dest='out_annot_inp', help='Name of output annot file')
    sp_filter.add_argument('-apt', action='store', dest='annot_type', default='RefSeq', help='Gencode/RefSeq')

    sp_proc_sample = subparsers.add_parser('proc_sample', help='Proc Sample Arguments')
    sp_proc_sample.add_argument('-ctp', action='store', dest='ctp_inp', help='Cell type Fofn')
    sp_proc_sample.add_argument('-pto', action='store', dest='pto_inp', default='.', help='Path to Output')
    sp_proc_sample.add_argument('-expr', action='store', dest='expr_inp', help='Name of expression file')

    sp_proximity_list = subparsers.add_parser('proximity_list', help='proximity list argument')
    sp_proximity_list.add_argument('-pst', action='store', dest='pst_inp', type=int, default=1, help='pend start index')
    sp_proximity_list.add_argument('-ped', action='store', dest='ped_inp', type=int, default=2, help='pend end index')
    sp_proximity_list.add_argument('-nst', action='store', dest='nst_inp', type=int, default=1, help='nend start index')
    sp_proximity_list.add_argument('-ned', action='store', dest='ned_inp', type=int, default=2, help='nend end index')
    sp_proximity_list.add_argument('-twn', action='store', dest='twin_inp', type=int, default=5000, help='nend end index')
    sp_proximity_list.add_argument('-glf', action='store', dest='glf_inp', help='gene list file')

    sp_merge_features = subparsers.add_parser('merge_features', help='merge features argument')
    sp_merge_features.add_argument('-pst', action='store', dest='pst_inp', type=int, default=1, help='pend start index')
    sp_merge_features.add_argument('-ped', action='store', dest='ped_inp', type=int, default=2, help='pend end index')
    sp_merge_features.add_argument('-nst', action='store', dest='nst_inp', type=int, default=1, help='nend start index')
    sp_merge_features.add_argument('-ned', action='store', dest='ned_inp', type=int, default=2, help='nend end index')
    sp_merge_features.add_argument('-tif', action='store', dest='tif_inp', help='tss interpolation file')

    sp_interp = subparsers.add_parser('interp', help='interp argument')
    sp_interp.add_argument('-sigma', action='store', dest='sigma_inp', type=int, default=50, help='Value of sigma')
    sp_interp.add_argument('-nip', action='store', dest='nip_inp', type=int, default=500, help='Number of interp points')
    sp_interp.add_argument('-ibin', action='store', dest='ibin_inp', type=int, default=5000, help='Size of bin around TSS/GB')
    sp_interp.add_argument('-ach', action='store', dest='anch_win', type=int, default=100000, help='Anchor window')
    sp_interp.add_argument('-rfl', action='store', dest='refl_inp', type=int, default=500, help='RE flnk length')
    sp_interp.add_argument('-rff', action='store', dest='reff_inp', type=int, default=25, help='RE flnk features')
    sp_interp.add_argument('-rfn', action='store', dest='rfn_inp', help='Region fofn')
    sp_interp.add_argument('-tag', action='store', dest='tag_inp', help='Tag for output')
    sp_interp.add_argument('-sfn', action='store', dest='sfn_inp', help='Sample file name')
    sp_interp.add_argument('-fln', action='store_true', dest='flankNorm', default=False, help='Flank Normalized interpolation')
    sp_interp.add_argument('-gsl', action='store_true', dest='geneSelect', default=False, help='Interpolation  for slected genes?')
    sp_interp.add_argument('-frf', action='store_true', dest='fixedReFlnk', default=False, help='Fixed features for RE flank')
    sp_interp.add_argument('-mmg', action='store', dest='min_gene_meth', type=int, default=40, help='Minimum CpGs assayed')
    sp_interp.add_argument('-mmr', action='store', dest='min_re_meth', type=int, default=2, help='Minimum CpGs assayed')

    sp_setup_expr = subparsers.add_parser('setup_expr', help='setup expr argument')
    sp_setup_expr.add_argument('-intrp', action='store', dest='intrp_inp', help='Interpolation CSV file')
    sp_setup_expr.add_argument('-expr', action='store', dest='expr_inp', help='Expression file')
    sp_setup_expr.add_argument('-fef', action='store', dest='floor_expr', type=bool, default=True, help='Floor expression value?')
    sp_setup_expr.add_argument('-efv', action='store', dest='efv_inp', type=float, default=5.0, help='Expression floor value')
    sp_setup_expr.add_argument('-def', action='store', dest='dexpr_flag', type=int, default=1, help='Method for differential expression')
    sp_setup_expr.add_argument('--expr_cutoff', action='store', dest='expr_cutoff', type=int, default=2, help='fold change in expression cutoff, must use -def 3')

    sp_run_clf = subparsers.add_parser('run_clf', help='Run classification argument')
    sp_run_clf.add_argument('-dfi', action='store', dest='dfi_inp', help='Dataframe output from interpolation step')
    sp_run_clf.add_argument('-ntr', action='store', dest='ntr_inp', type=int, default=5001, help='Number of trees for Random Forest Classifier')
    sp_run_clf.add_argument('-npr', action='store', dest='npr_inp', type=int, default=8, help='Number of Processors for RF run')
    sp_run_clf.add_argument('-tag', action='store', dest='tag_inp', default='test', help='Tag for Output Writing')
    sp_run_clf.add_argument('-fsl', action='store', dest='fsl_inp', type=int, default=1, help='Feature Selection. 1: TSS; 2: TSS+RE')
    sp_run_clf.add_argument('-suf', action='store', dest='suf_inp', type=bool, default=True, help='Shuffle true ot false')
    sp_run_clf.add_argument('-ss', action='store_true', dest='ss', default=False, help='Single sample or not') 
    sp_run_clf.add_argument('-ngnorm', action='store_false', dest='gnorm', default=True, help='Normalize gene count or not') 

   
    sp_eval = subparsers.add_parser('eval', help='Evaluation argument')
    sp_eval.add_argument('-dfi', action='store', dest='dfi_inp', help='Dataframe output from classification')
    sp_eval.add_argument('-nstp', action='store', dest='steps_inp', type=int, default=101, help='Number of steps')

    sp_proc_reg = subparsers.add_parser('proc_reg', help='Process RE regions')
    sp_proc_reg.add_argument('-gpj', action='store', dest='gapj_inp', type=int, default=500, help='RE closer than this will be joined together.')
    sp_proc_reg.add_argument('-msz', action='store', dest='max_inp', type=int, default=1000, help='Maximum size of RE to take into account.')
    sp_proc_reg.add_argument('-twn', action='store', dest='twin_inp', type=int, default=2000, help='RE falling +/- is this region around TSS will be ignored.')
    sp_proc_reg.add_argument('-glf', action='store', dest='glf_inp', help='Gene information file')
 
    args = parser.parse_args()
    funct = FUNCTION_MAP[args.command]
    funct(args)

if __name__ == '__main__':
    main()
