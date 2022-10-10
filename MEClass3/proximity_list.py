from os import getcwd

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
    #path_to_input = os.getcwd()
    path_to_input = getcwd()

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