
import sys
import os
import argparse
from datetime import datetime

def exec_evaluate_performance(args):
        #
	inp_ref_file = args.ref_file_inp
	inp_pred_file = args.pred_file_inp
	inp_isz = args.isz_inp
	output_file_name = args.out_file_inp
	#
	with open(inp_ref_file, 'r') as ref_file:
		ref_file_lines = ref_file.readlines()
	ref_file.close()
	#
	with open(inp_pred_file, 'r') as pred_file:
		pred_file_lines = pred_file.readlines()
	pred_file.close()
	#
	dict_ref = {}
	dict_ref_keys = []
	for line in ref_file_lines:
		items = line.strip().split()
		te_type = items[0]
		chrom = items[1]
		ip = items[2]
		temp_ref_ip_id = te_type+'-'+chrom+'-'+ip
		dict_ref[ temp_ref_ip_id ] = []
		dict_ref_keys.append(temp_ref_ip_id ) # definite order compared to doing {for key in dict_ref}:
	#
	track_pred_ids = []
	#
	for key in dict_ref_keys:
		#
		pred_te_type = key.strip().split('-')[0]
		pred_chrom = key.strip().split('-')[1]
		pred_ip = int(key.strip().split('-')[2])
		#
		for line in pred_file_lines:
			items = line.strip().split()
			te_type = items[0]
			chrom = items[1]
			ip = int(items[2])
			if ( pred_te_type == te_type ) and ( pred_chrom == chrom ) and ( pred_ip in range(ip-inp_isz, ip+inp_isz+1) ) and ( line not in track_pred_ids):
				dict_ref[ key ].append( (te_type, chrom, ip ) )
				track_pred_ids.append( line )

	# Count true positive and false negative
	true_positive = 0
	true_positive_list = []
	false_negative = 0
	false_negative_list = []
	#
	for key in dict_ref_keys:
		if len(dict_ref[key]) == 0:
			false_negative += 1
			false_negative_list.append( ( key.strip().split('-')[0], key.strip().split('-')[1], key.strip().split('-')[2] )  ) 
		else:
			true_positive += 1
			true_positive_list.append( ( key.strip().split('-')[0], key.strip().split('-')[1], key.strip().split('-')[2] )  ) 

	# Count False positive
	false_positive_list = []
	false_positive = len( pred_file_lines ) - len( list(set( track_pred_ids ) ))
	for line in pred_file_lines:
		if line in track_pred_ids:
			continue
		false_positive_list.append( ( line.strip().split()[0], line.strip().split()[1], line.strip().split()[2] )   )
	#
	output_file = open( output_file_name, 'w')
	output_file.write('# '+ str( datetime.now() ) + '\n')
	output_file.write('# Ref file: '+ args.ref_file_inp + '; Pred file: ' + args.pred_file_inp + '; Insert size: ' + str(args.isz_inp)  + '\n')
	# Write output
	#sys.stdout.write('True Positive = '+ str( true_positive ) +'\n')
	output_file.write('True Positive = '+ str( true_positive ) +'\n')
	for item in true_positive_list:
		output_file.write('TP\t'+item[1]+'\t'+str(item[2])+'\n')
	#
	#sys.stdout.write('False Positive = '+ str( false_positive ) +'\n')
	output_file.write('False Positive = '+ str( false_positive ) +'\n')
	for item in false_positive_list:
		output_file.write('FP\t'+item[1]+'\t'+str(item[2])+'\n')
	#
	#sys.stdout.write('False Negative = '+ str(false_negative ) +'\n')
	output_file.write('False Negative = '+ str(false_negative ) +'\n')
	for item in false_negative_list:
		output_file.write('FN\t'+item[1]+'\t'+str(item[2])+'\n')

	output_file.close()
	#

def main(): 
#	FUNCTION_MAP = {
#			'evaluate_performance' : exec_evaluate_performance,
#			}
#
#	parser = argparse.ArgumentParser()
#	parser.add_argument('module', choices=FUNCTION_MAP.keys())
#	subparsers = parser.add_subparsers()
#
#	sp_evaluate_performance = subparsers.add_parser('evaluate_performance_opt', help="extract reads argument")
#	sp_evaluate_performance.add_argument('-ref', action='store', dest='ref_file_inp', required=True, help='Bam(.bam) file with full path')
#	sp_evaluate_performance.add_argument('-pred', action='store', dest='pred_file_inp', required=True, help='region file with full path')
#	sp_evaluate_performance.add_argument('-out', action='store', dest='out_file_inp', type=str, default='performance_evaluation_output.txt', help='region file with full path')
#	sp_evaluate_performance.add_argument('-isz', action='store', dest='isz_inp', type=int, default=360, help='Read range')

	parser = argparse.ArgumentParser()
	parser.add_argument('-ref', action='store', dest='ref_file_inp', required=True, help='Bam(.bam) file with full path')
	parser.add_argument('-pred', action='store', dest='pred_file_inp', required=True, help='region file with full path')
	parser.add_argument('-out', action='store', dest='out_file_inp', type=str, default='performance_evaluation_output.txt', help='region file with full path')
	parser.add_argument('-isz', action='store', dest='isz_inp', type=int, default=360, help='Read range')

	args = parser.parse_args()
#	funct = FUNCTION_MAP[args.module]
#	funct(args)
	exec_evaluate_performance(args)

if __name__ == '__main__':
	import sys
	import os
	import argparse
	from datetime import datetime
	main()

