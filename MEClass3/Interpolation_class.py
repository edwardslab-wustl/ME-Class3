import scipy
import numpy as np

class Interpolation:
    #
    def __init__(self, cpos_dat, dmet_dat, reg_start, reg_end, strand, reg_type, args):
        #
        self.cpos_raw = cpos_dat
        self.dmet_raw = dmet_dat
        self.sigma = args.sigma # 50
        self.num_intrp_point = args.num_interp_points #5
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
