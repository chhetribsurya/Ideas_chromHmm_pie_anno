import os
import scipy.stats
import json
import time
import pandas as pd
import pybedtools


start_time = time.time()


directory = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/ENCODE_poster/chromHMM_overlap_encode")
if not os.path.exists(directory):
    os.makedirs(directory)
    #os.rmdir(prior_dir)
    #shutil.rmtree()
    

####Input files for script to run are attached in the dir input_file_for script:
chromHMM_hepg2_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/chrom_impute_chromHMM_data/corrected/E118_25_imputed12marks_dense.bed")  
#chromHMM_hepg2_file = os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/chrom_impute_chromHMM_data/corrected/head_200_E118_25_imputed12marks_dense.bed")  

peaks_bedfile1= os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/encode_org_chip_data/Pol2/myers/ENCFF001UKU.bed")
peaks_bedfile2= os.path.expanduser("~/Dropbox/local_miscellaneous_data/test_data/encode_org_chip_data/Pol2/myers/ENCFF002CMI.bed")

####Output files from the script run (provide the path location to save the file):
#Provide the file path location to save the output of ENCODE_TF_cluster file as JSON format (could be useful for future usage):
JSON_dict_file = os.path.expanduser(os.path.join(directory,"ENCODE_chromHMM_JSON_1.txt"))


master_dict_return = {}

def parse_chromHMM(chromHMM_file, Json_output_file):
	with open(chromHMM_file,"r") as file:
		with open(JSON_dict_file, 'w') as outfile:
			for line in file.readlines():
				#print line
				splitted = line.strip("\n").split("\t",4)
				chrom, start, end, ch_state = splitted[0], splitted[1], splitted[2], splitted[3]
				cell_line = splitted[4].split()[1]
				line_info = [chrom, str(start), str(end), str(ch_state)]

				if ch_state in master_dict_return.keys():
					if chrom in master_dict_return[ch_state].keys():
						master_dict_return[ch_state][chrom].append((chrom, int(start), int(end), "\t".join(line_info)))
					else:
						master_dict_return[ch_state][chrom] = [(chrom, int(start), int(end), "\t".join(line_info))]
				else:
					master_dict_return[ch_state] = {chrom:[(chrom, int(start), int(end), "\t".join(line_info))]}
					master_dict_return[ch_state].update({"A_bed_count_hits":0})
					master_dict_return[ch_state].update({"B_bed_count_hits":0})
					master_dict_return[ch_state].update({"custom_overlap_list":[]})
					master_dict_return[ch_state].update({"As_overlap_list":[]})
					master_dict_return[ch_state].update({"Bs_overlap_list":[]})

			json.dump(master_dict_return, outfile)
				#print master_dict_return
			return(master_dict_return)

#master_chromHMM_dict = parse_chromHMM(chromHMM_hepg2_file, JSON_dict_file)
#print "\n\nTime to parse chromHMM_segmentation based on region and chromosome = ", time.time()-start_time


def parse_chromHMM_based_on_region(chromHMM_file, Json_output_file):
	with open(chromHMM_file,"r") as file:
		with open(JSON_dict_file, 'w') as outfile:
			for line in file.readlines():
				#print line
				splitted = line.strip("\n").split("\t",4)
				chrom, start, end, ch_state = splitted[0], splitted[1], splitted[2], splitted[3]
				cell_line = splitted[4].split()[1]
				line_info = [chrom, str(start), str(end), str(ch_state)]

				if ch_state in master_dict_return.keys():
					master_dict_return[ch_state]["loci"].append([chrom, int(start), int(end)])
					#master_dict_return[ch_state]["loci"].append([chrom, int(start), int(end), "\t".join(line_info)])
					
				else:
					master_dict_return[ch_state] = {"loci" : [[chrom, int(start), int(end)]]}
					#master_dict_return[ch_state] = {"loci" : [[chrom, int(start), int(end), "\t".join(line_info)]]}
					master_dict_return[ch_state].update({"A_bed_count_hits":0})
					master_dict_return[ch_state].update({"B_bed_count_hits":0})
					master_dict_return[ch_state].update({"custom_overlap_list":[]})
					master_dict_return[ch_state].update({"As_overlap_list":[]})
					master_dict_return[ch_state].update({"Bs_overlap_list":[]})

			json.dump(master_dict_return, outfile)
				#print master_dict_return
			return(master_dict_return)

master_chromHMM_dict = parse_chromHMM_based_on_region(chromHMM_hepg2_file, JSON_dict_file)
print "\n\nTime to parse chromHMM_segmentation based on region = ", time.time()-start_time




#################
#################




start_time = time.time()
print "\n\n Overlap analysis starting..."


#Check if file exists to prevent the append mode for each keys/ch_state value being written to the file:
#prior_file = directory + "/" +file_name
file_name_list = [key+".txt" for key in master_chromHMM_dict.iterkeys()]
for file_name in file_name_list:
	prior_file = os.path.join(directory, file_name)
	if os.path.exists(prior_file):
		os.remove(prior_file) 


def load_stringf_for_peak_files(input_file):

	df = pd.read_csv(input_file, sep= "\t", header= None)
	df = df.iloc[:,0:6]
	df.columns = ["chrom", "start", "end", "name", "score", "strand"]
	TF_peak_file = df.loc[:,["chrom","start", "end"]]
	#df = df.ix[:,0:6]
	#df.loc[:,"end"] = df["end"] + 1
	#df.loc[:,"strand"] = "."
	#Cpg_bed_filter =	TF_peak_file.loc[(df["chrom"] != "*")]
	#Cpg_list = Cpg_bed_filter.values.tolist()
	TF_peak_list =	TF_peak_file.values.tolist()
	TF_peak_string_list = [ "\t".join(map(str, line)) for line in TF_peak_list]
	TF_peak_string_bed_format = "\n".join(TF_peak_string_list)

	TF_bed_file = pybedtools.BedTool(TF_peak_string_bed_format,from_string=True)

	return(TF_bed_file)

peak_file = load_stringf_for_peak_files(peaks_bedfile1)



for key,value in master_chromHMM_dict.iteritems():
	# name the output file with its keyword:
	file_name = key.replace("/","-") + ".txt"
	outfile_path = os.path.join(directory, file_name)
	
	each_tf_dict = master_chromHMM_dict[key]
	print "processing....",  key, "final......."
	TF_loci_list = each_tf_dict.get("loci")
	TF_loci_string_list = [ "\t".join(map(str, line)) for line in TF_loci_list]
	TF_loci_string_bed_format = "\n".join(TF_loci_string_list)
	chromHMM_loci_bed_file = pybedtools.BedTool(TF_loci_string_bed_format,from_string=True)

	# BEDTool intersect operation:
	bed_intersect = peak_file.intersect(chromHMM_loci_bed_file, wa = True, wb = True)
	intersect_count = bed_intersect.count()

	if intersect_count >= 1:
		# Bedtool object converted to pandas dataframe:
		df_bed = pd.read_table(bed_intersect.fn, sep="\t", header=None)
		As_overlap_list = df_bed.iloc[:,0:3]
		Bs_overlap_list = df_bed.iloc[:,3:6]
		custom_overlap_list = df_bed

		print key, "appending As overlap_list"
		master_chromHMM_dict[key]["As_overlap_list"].append(As_overlap_list)
		master_chromHMM_dict[key]["As_overlap_list"][0] = master_chromHMM_dict[key]["As_overlap_list"][0].drop_duplicates()
		master_chromHMM_dict[key]["A_bed_count_hits"] = len(master_chromHMM_dict[key]["As_overlap_list"][0].index)


		master_chromHMM_dict[key]["Bs_overlap_list"].append(Bs_overlap_list)
		master_chromHMM_dict[key]["Bs_overlap_list"][0] = master_chromHMM_dict[key]["Bs_overlap_list"][0].drop_duplicates()
		master_chromHMM_dict[key]["B_bed_count_hits"] = len(master_chromHMM_dict[key]["Bs_overlap_list"][0].index)

		master_chromHMM_dict[key]["custom_overlap_list"].append(df_bed)

		df_bed.to_csv(outfile_path, sep ="\t", header = True, index = False)
		outfile_path_prefix = os.path.splitext(outfile_path)
		master_chromHMM_dict[key]["As_overlap_list"][0].to_csv(outfile_path_prefix[0] + "_with_A_overlap.bed")




print "\n overlap b/w bed files completed!!!!!"
print "Time for overlap operation = ", time.time()-start_time

out_file = directory + "/" + "intersectbed_counts_with_chromHMM.txt"
with open(out_file, "w") as outfile:
	for key,value in master_chromHMM_dict.items():
		result = "%s\t%s" %(key, master_chromHMM_dict[key]["A_bed_count_hits"])
		outfile.write(result + "\n")
		print result


print "\nEnd of job!!\n"

