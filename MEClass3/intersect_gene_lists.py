import os

from MEClass3.io_functions import print_to_log
from MEClass3.io_functions import format_args_to_print
from MEClass3.io_functions import mk_output_dir

def exec_intersect_gene_lists(args):
    with open(args.logfile, 'w') as log_FH:   
        print_to_log(log_FH, format_args_to_print(args))
        mk_output_dir(args.outDir)
        header_dict = dict()
        gene_dict = dict()
        file_count = 0
        for file in args.input_files:
            file_count += 1
            header_dict, gene_dict = add_interp_to_dict(file, header_dict, gene_dict)
        gene_set = prune_num_dict(gene_dict, file_count)
        for file in args.input_files:
            result = get_set_data_from_interp(file, gene_set)
            out_file_name = args.outDir + "/" + os.path.basename(file)
            with open(out_file_name, "w") as out_FH:
                out_FH.write(result)
    return
 
 
def get_set_data_from_interp(file, gene_set):       
    result = ""
    with open(file, 'r') as FH:
        for line in FH:
            if line.startswith("#") or line.startswith("gene_id-sample_name"):
                result += line
            else:
                id, data = line.split(",", 1)
                if id in gene_set:
                    result += line
    return result
     
     
def prune_num_dict(num_dict, num):
    return_set = set()
    for key, val in num_dict.items():
        if val >= num:
            return_set.add(key)
    return return_set
        
           
def add_interp_to_dict(file, header_dict, gene_dict):
    with open(file, 'r') as FH:
        for line in FH:
            if line.startswith("#") or line.startswith("gene_id-sample_name"):
                if file in header_dict:
                    header_dict[file] += line
                else:
                    header_dict[file] = line
            else:
                id, data = line.split(",", 1)
                gene_dict = increment_dict(gene_dict, id)
    return header_dict, gene_dict

def increment_dict(num_dict, key):
    if key in num_dict:
        num_dict[key] += 1
    else:
        num_dict[key] = 1
    return num_dict
                    
        
def exec_intersect_gene_lists_help(parser):
    parser_required = parser.add_argument_group('required arguments')
    parser_required.add_argument('input_files',
        nargs='+',help="output files from merge_data in csv format")
    parser.add_argument('-o', '--outDir', default="filtered_files",
        help="files will have same name, but be stored in the specified directory.")
    parser.add_argument('--logfile',
        default='filtered_files.log', help='log file')
    return parser