#!/bin/bash
import os
import re
import pandas as pd
import numpy as np
import glob
import pybedtools
import json
import time
import subprocess
import warnings

xls_file = os.path.expanduser("~/for_chris/HepG2_SLs.xls.xlsx")
dir_path = os.path.expanduser("~/for_chris/batch_I/idr_passed_peaks/SL*") 
chip_analysis_path = os.path.expanduser("~/for_chris/batch_I/chip_analysis/IDR*")
file_suffix = "intersectbed_counts_with_ideas.txt"


###################################################
###################################################
###################################################

read_xls =  pd.read_excel(xls_file, sheetname = 0)
df_xls =  read_xls.iloc[:,[0,1,2,3,4,6]]

Target_col_list = [u'Exp ID rep1', u'Control ID rep1', u'Target', u'Exp ID rep2', u'Control ID rep2']
df_xls_select = df_xls.loc[:, df_xls.columns.isin(Target_col_list)]
df_xls_select.columns = ["rep1", "control1", "tf_name", "rep2", "control2",]
df_xls_ordered = df_xls_select.loc[:,["rep1","rep2","control1","control2","tf_name"]]

df_xls["combined_col"] = df_xls_ordered.dropna().apply(lambda x: '_'.join(x.astype(str).values), axis=1)
SL_all_split_list = df_xls[u'combined_col'].dropna().str.split("_", 4).tolist()

SL_cross_check_list = [[],[]]
for each in SL_all_split_list:
    SL_cross_check_list[0].append("_".join(each[0:2]))
    SL_cross_check_list[1].append(each[4])


idr_file_list = glob.glob(dir_path)

count = 0
tf_filepath_metadata = [[],[],[]]
for each_idr in idr_file_list:
    each_idr_basename= os.path.basename(each_idr)
    test_str = each_idr_basename   #test_str = 'SL101267.filt.nodup.srt.SE_VS_SL115566.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak.gz'
    #idr_string = filter(None, re.split("[.,_]", test_str))   #idr_string = re.split("[.,_]", test_str)
    idr_string = re.split("[.,_]", test_str)  #idr_string = re.split("[.,_]", test_str)
    idr_SL_list = [ idr_string[i] for i in [0,6]]
    idr_check_string= "_".join(idr_SL_list)

    # compare the SL string from the idr_files to the SL no. from the excel file:
    for i,item in enumerate(SL_cross_check_list[0]):
        if idr_check_string == item :
            count += 1
            new_file_name =  each_idr + "_" + SL_cross_check_list[1][i]

            tf_filepath_metadata[0].append(SL_cross_check_list[1][i])
            tf_filepath_metadata[1].append(each_idr)
            tf_filepath_metadata[2].append(idr_check_string)
            #print count, SL_cross_check_list[1][i]

#####################################################
#####################################################
#####################################################

my_analysis_list = df_xls["Prelim-Analysis on"].dropna().tolist()
start_time = time.time()


output_dir = os.path.expanduser("~/for_chris/batch_I/chromHMM_overlap_encode")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

####Input files for script to run are attached in the dir input_file_for script:
chromHMM_hepg2_file = os.path.expanduser("~/for_chris/batch_I/E118_25_imputed12marks_dense.bed")  

#Provide the file path location to save the output of ENCODE_TF_cluster file as JSON format (could be useful for future usage):
JSON_dict_file = os.path.expanduser(os.path.join(output_dir,"ENCODE_chromHMM_JSON_1.txt"))



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

def load_stringf_for_centerOfpeak_files(input_file):

	df = pd.read_csv(input_file, sep= "\t", header= None)
	df = df.iloc[:,0:6]
	df.columns = ["chrom", "start1", "end1", "name", "score", "strand"]
	df["mid_point"] = (df["start1"] + df["end1"])/2
	df["start"] = df["mid_point"].astype(int)
	df["end"] = df["mid_point"].astype(int) + 1
	
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


for each in my_analysis_list:
    index = tf_filepath_metadata[0].index(each)
    tf_name = tf_filepath_metadata[0][index]
    file_path = tf_filepath_metadata[1][index]

    #final_output_file = output_dir + "/" + tf_name + "_" + "intersectbed_counts_with_ideas.txt"
    final_output_file = os.path.join(output_dir, (tf_name + "_" + file_suffix))
	if not os.path.exists(final_output_file):  
	    peaks_bedfile1 = file_path
	    master_dict_return = {}

	    master_chromHMM_dict = parse_chromHMM_based_on_region(chromHMM_hepg2_file, JSON_dict_file)
	    print "\n\nTime to parse chromHMM_segmentation based on region = ", time.time()-start_time

	    start_time = time.time()
	    print "\n\n Overlap analysis starting..."

	    #Check if file exists to prevent the append mode for each keys/ch_state value being written to the file:
	    #prior_file = output_dir + "/" +file_name
	    file_name_list = [key+".txt" for key in master_chromHMM_dict.iterkeys()]
	    for file_name in file_name_list:
		prior_file = os.path.join(output_dir, file_name)
		if os.path.exists(prior_file):
		    os.remove(prior_file)

	    peak_file = load_stringf_for_centerOfpeak_files(peaks_bedfile1)


	    for key,value in master_chromHMM_dict.iteritems():
		# name the output file with its keyword:
		file_name = key.replace("/","-") + ".txt"
		outfile_path = os.path.join(output_dir, file_name)

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

	    out_file = output_dir + "/" + tf_name + "_" + "intersectbed_counts_with_chromHMM.txt"
	    with open(out_file, "w") as outfile:
		for key,value in master_chromHMM_dict.items():
		    result = "%s\t%s" %(key, master_chromHMM_dict[key]["A_bed_count_hits"])
		    outfile.write(result + "\n")
		    print result

	    print "\nIntersection job for %s completed....!!\n" %(tf_name)

    #chromHMM_file = output_dir + "/" + tf_name+ "_" + "intersectbed_counts_with_chromHMM.txt"
    chromHMM_file = os.path.join(output_dir, (tf_name + "_" + file_suffix))
    os.environ["chromHMM_file"] =  chromHMM_file
    os.environ["tf_name"] = tf_name
    os.environ["out_dir_name"] = output_dir
    os.system('Rscript ./ChromHMM_piechart_morgan.R $chromHMM_file $tf_name $out_dir_name')
    #subprocess.call("Rscript " + "ChromHMM_piechart_morgan.R " +  chromHMM_file + tf_name + output_dir], shell=True)
    #subprocess.call("Rscript ChromHMM_piechart_morgan.R --args chromHMM_file tf_name output_dir", shell=True)
    print "\nRunning of Rscript for plotting completed!!!....\n"
    print "\nCheck your plots in %s dir\n" %(output_dir)
    
