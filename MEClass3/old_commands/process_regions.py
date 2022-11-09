import os
import sys

import pandas as pd

   
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

def exec_proc_reg_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('-glf', action='store', dest='glf_inp', required=True, help='Gene information file')
    parser.add_argument('-gpj', action='store', dest='gapj_inp', type=int, default=500, help='RE closer than this will be joined together.')
    parser.add_argument('-msz', action='store', dest='max_inp', type=int, default=1000, help='Maximum size of RE to take into account.')
    parser.add_argument('-twn', action='store', dest='twin_inp', type=int, default=2000, help='RE falling +/- is this region around TSS will be ignored.')
    parser._action_groups.reverse()
    return(parser)