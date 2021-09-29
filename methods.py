import os.path
import subprocess

""" 
HAS TO BE RUN IN PYTHON 3!!!
"""

def filter_strain(strain_name,n_bp,fam_name,type_name):
    # setup for "input_for_bedtools.bed"
    # setup for "input_for_bedtools.bed" file location
    save_path = "/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/inputs/input_for_bedtools"
    name_of_file = "input_for_bedtools_"+strain_name
    completeName = os.path.join(save_path,name_of_file+".bed")

    # open "strain_name.fa.out"
    f = open("/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/files/"+strain_name+".fa.out","r")
    lines = f.readlines()[3:]

    for line in lines:
        x = line.split()
        d = int(x[6]) - int(x[5])
        if (d>=n_bp) and (fam_name == x[9]) and (type_name == x[10]):
            with open(completeName,"a") as bedfile:
                bedfile.write(x[4] + "\t")
                bedfile.write(x[5] + "\t")
                bedfile.write(x[6] + "\t")
                bedfile.write(x[9] + "\t")
                bedfile.write(x[10] + "\t")
                bedfile.write("\n")
            bedfile.close()
    f.close()

def bedtools_runner(strain_name,fam_name,type_name):    
    # will output a file named: Unknown-fam-33.fa
    # runs bedtools via bashcript, outputs to /Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/outputs/output_from_bedtools

    # has to cd into /Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/inputs/input_for_bedtools
    subprocess.run(["cd","/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/inputs/input_for_bedtools"])
    
    # ex: fam_name = "rnd-1_family-33", then stripped_fam_name = "family-33"
    stripped_fam_name = fam_name.partition("_")[-1]

    # gets put into subprocess.run()
    subprocess_argument_frag = ["bedtools", "getfasta","-fi","/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/files/"+strain_name+".fa","-bed","/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/inputs/input_for_bedtools/input_for_bedtools_"+strain_name+".bed"]

    # will return the output
    output = subprocess.run(subprocess_argument_frag,stdout = subprocess.PIPE)
    
    # append to other file location
    save_path = "/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/outputs/output_from_bedtools"
    name_of_file = strain_name+"-"+type_name+"-"+stripped_fam_name
    completeName = os.path.join(save_path,name_of_file+".fa")
    with open(completeName,"a") as fa_file:
        fa_file.write(output.stdout.decode("utf-8"))
    fa_file.close()

def sixpack_primer(strain_name,fam_name,type_name):
    # will prime bedtools output so that it can be sent through sixpack
    stripped_fam_name = fam_name.partition("_")[-1]
    f = open("/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/outputs/output_from_bedtools/"+strain_name+"-"+type_name+"-"+stripped_fam_name+".fa","r")
    file_list = []

    # counter for sequence files (ex: "sequence_<a>")
    a = 1

    for line in f:
        file_list.append(line)

    f.close() 
    
    # setup for file save location
    save_path = "/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/inputs/input_for_six_pack/"

    # counter for iterator
    i = 0
    while i<len(file_list):
        with open(save_path+"input_for_sixpack_"+strain_name+"_sequence_"+str(a),"a") as primer:
            primer.write(file_list[i])
            primer.write(file_list[i+1])
        primer.close()
        a+=1
        i+=2

def sixpack(strain_name,fam_name,type_name):
    # will run sixpack through subprocess and output .seq and .translated files for the corresponding "input_for_six_pack" file in "output_for_six_pack" folder
    # how many sequence files are there
    stripped_fam_name = fam_name.partition("_")[-1]
    f = open("/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/outputs/output_from_bedtools/"+strain_name+"-"+type_name+"-"+stripped_fam_name+".fa","r")
    text = f.read()
    num_caret = text.count(">")
    f.close()

    # paths into output and input folder for sixpack
    sixpack_output_path = "/Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/outputs/output_from_six_pack/"
    sixpack_input_path = "/Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/inputs/input_for_six_pack/"
    
    # utilizing for loop in order to run sixpack on each "input_for_six_pack" file
    for n in range(1,num_caret+1):
        input_file = sixpack_input_path+"input_for_sixpack_"+strain_name+"_sequence_"+str(n)
        output_translated = sixpack_output_path+"output_for_sixpack_"+strain_name+str(n)+".translated"
        output_seq = sixpack_output_path+"output_for_sixpack_"+strain_name+str(n)+".seq"
        subprocess_argument_frag = ["sixpack","-mstart","-table","0","-sequence",input_file,"-outfile",output_translated,"-outseq",output_seq]
        subprocess.run(subprocess_argument_frag)

"""
def longest_orf(strain_name,fam_name,type_name):
    # will find the longest_orf for each .seq file in "output_for_six_pack"
    
    # how many sequence files are there
    stripped_fam_name = fam_name.partition("_")[-1]
    f = open("/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/outputs/output_from_bedtools/"+strain_name+"-"+type_name+"-"+stripped_fam_name+".fa","r")
    text = f.read()
    num_caret = text.count(">")
    f.close()
    
    # counter represents 1 - 56 for the .seq files
    a = 1

    #going from 0 ~ num_caret for .seq files
    with open("/Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/outputs/longest_orf.fa","a") as l:
        while a < num_caret+1:
            current_file_name = "/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/outputs/output_from_six_pack/output_for_sixpack_"+strain_name+str(a)+".seq"

            with open(current_file_name) as f:
                first_line = f.readline().rstrip()
            
            l.write(first_line)
            l.write(longest_orf(current_file_name)) 
            l.write("\n")

            #increment a by 1
            a+=1
"""