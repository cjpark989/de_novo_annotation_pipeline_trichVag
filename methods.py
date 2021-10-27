import os.path
import subprocess

""" 
HAS TO BE RUN IN PYTHON 3!!!
"""

def categorize(strain_name,type_name):
    # will output a list with all the fam_name listed 
    # these fam_name would be associated with the type_name
    save_path = "/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/files"
    name_of_file = "category_"+strain_name+"_"+type_name
    #completeName = os.path.join(save_path,name_of_file)

    # open strain_name.fa.out and read from line 3
    f = open("/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/files/"+strain_name+".fa.out","r")
    lines = f.readlines()[3:]

    # create list to put in all the fam_name
    l = []

    for line in lines:
        x = line.split()
        if (x[10]==type_name) and (x[9] not in l):
            l.append(x[9])
    f.close()
    return l

###########################################################################################################################################################################################################################

def filter_strain(strain_name,n_bp,fam_name,type_name):
    # setup for "input_for_bedtools.bed"
    # setup for "input_for_bedtools.bed" file location
    save_path = "/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/inputs/input_for_bedtools"
    name_of_file = "input_for_bedtools_"+strain_name+"_"+fam_name
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
    file_name = "/Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/inputs/input_for_bedtools/" + "input_for_bedtools_"+strain_name+"_"+fam_name+".bed"
    # has to cd into /Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/inputs/input_for_bedtools
    subprocess.run(["cd","/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/inputs/input_for_bedtools"])
        
        
    # gets put into subprocess.run()
    subprocess_argument_frag = ["bedtools", "getfasta","-fi","/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/files/"+strain_name+".fa","-bed","/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/inputs/input_for_bedtools/input_for_bedtools_"+strain_name+"_"+fam_name+".bed"]

    # will return the output
    output = subprocess.run(subprocess_argument_frag,stdout = subprocess.PIPE)
    
    # append to other file location
    save_path = "/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/outputs/output_from_bedtools"
    name_of_file = strain_name+"-"+type_name+"-"+fam_name
    completeName = os.path.join(save_path,name_of_file+".fa")
    with open(completeName,"a") as fa_file:
        fa_file.write(output.stdout.decode("utf-8"))
    fa_file.close()
  
def sixpack_primer(strain_name,fam_name,type_name):
    # will prime bedtools output so that it can be sent through sixpack
    f = open("/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/outputs/output_from_bedtools/"+strain_name+"-"+type_name+"-"+fam_name+".fa","r")
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

    # could switch the locations of while loop and with file open statement
    while i<len(file_list):
        with open(save_path+"input_for_sixpack_"+strain_name+"_sequence_"+fam_name+"_"+str(a),"a") as primer:
            primer.write(file_list[i])
            primer.write(file_list[i+1])
        primer.close()
        a+=1
        i+=2

def sixpack_runner(strain_name,fam_name,type_name):
    # will run sixpack through subprocess and output .seq and .translated files for the corresponding "input_for_six_pack" file in "output_for_six_pack" folder
    # how many sequence files are there
    stripped_fam_name = fam_name.partition("_")[-1]
    f = open("/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/outputs/output_from_bedtools/"+strain_name+"-"+type_name+"-"+fam_name+".fa","r")
    text = f.read()
    num_caret = text.count(">")
    f.close()

    # paths into output and input folder for sixpack
    sixpack_output_path = "/Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/outputs/output_from_six_pack/"
    sixpack_input_path = "/Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/inputs/input_for_six_pack/"
    
    # utilizing for loop in order to run sixpack on each "input_for_six_pack" file
    for n in range(1,num_caret+1):
        # strain_name = "T_Stableri_BTPI"
        # fam_name = "rnd-1_family-33"
        # n is just a counter
        input_file = sixpack_input_path+"input_for_sixpack_"+strain_name+"_sequence_"+fam_name+"_"+str(n)
        output_translated = sixpack_output_path+"output_from_sixpack_"+strain_name+"_"+fam_name+"_"+str(n)+".translated"
        output_seq = sixpack_output_path+"output_from_sixpack_"+strain_name+"_"+fam_name+"_"+str(n)+".seq"
        subprocess_argument_frag = ["sixpack","-mstart","-table","0","-sequence",input_file,"-outfile",output_translated,"-outseq",output_seq]
        subprocess.run(subprocess_argument_frag)

# example input = longest_orf("T_stableri_BTPI","rnd-1_family-33","Unknown")
# only refer to the "ouput_from_six_pack" file:
    # file name will be: "output_for_sixpack_"+"T_Stableri_BTPI" (strain_name)+"_"+"rnd-1_family-33" (fam_name)+"_"+(1~57)+".seq"
# look for those files.
# open those files, and then, go through the .seq file and find the longest orf. 
# make an longest_orf_T_stableri_btpi_rnd-1_family-33 file in "files". There should be "n" amount of them per family in the unknown type of strain_name
# will get a shit ton of files.

def longest_orf(strain_name,fam_name,type_name):
    # one fam_name will go through X amounts of iterations
    # DEFINE X
    f = open("/Users/chungjunepark/Documents/"+strain_name+"/files_and_outputs/outputs/output_from_bedtools/"+strain_name+"-"+type_name+"-"+fam_name+".fa","r")
    text = f.read()
    x = text.count(">")
    f.close()

    for n in range(1,x+1):
        # list that will contain all of fam_name_n.seq's outputs
        l = []
        # list that will contain the indexes of all the lines that start with ">"
        caret = []
        # file path + name for input file
        file_name = "/Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/outputs/output_from_six_pack/output_from_sixpack_"+strain_name+"_"+fam_name+"_"+str(n)+".seq"
        # read file_name and get the sequences
        # f_2 is a list of all the separate lines in the file as its elements
        f_2 = open(file_name).readlines()

        # iterate through the lines list and find the indexes of the ">" occurring ones
        for element in f_2:
            if element.startswith(">"):
                caret.append(f_2.index(element))
        # iterate through the index list and concatinate lines that are inbetween those indexes
        for index,element in enumerate(caret):
            # element is what we're on right now
            # caret[index+1] is the next element
            if index+1 < len(caret):
                l.append("".join(f_2[element+1:caret[index+1]]))
            else:
                break

        # longest_orf for specific n in specific fam_name
        longest_orf = max(l, key = len)
        while not longest_orf.startswith("M"):
            l.remove(longest_orf)
            longest_orf = max(l, key=len)
        str1 = ""
        longest_orf = str1.join(longest_orf)

        # create a file in longest_orf to save .txt file of longest_orf
        file_path = "/Users/chungjunepark/Documents/T_stableri_BTPI/files_and_outputs/longest_orf/longest_orf_"+strain_name+"_"+fam_name+"_"+str(n)+".txt"

        with open(file_path,"a") as f_3:
            f_3.write(">"+fam_name+"_"+str(n)+"\n")
            f_3.write(longest_orf)
        
        f_3.close()

#filter_strain("T_Stableri_BTPI",1000,"rnd-1_family-33","Unknown")
#bedtools_runner("T_Stableri_BTPI","rnd-1_family-33","Unknown")
#sixpack_primer("T_Stableri_BTPI","rnd-1_family-33","Unknown")
#sixpack_runner("T_Stableri_BTPI","rnd-1_family-33","Unknown")
#longest_orf("T_Stableri_BTPI","rnd-1_family-33","Unknown")