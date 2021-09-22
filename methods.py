import itertools as it


def LTR(genome_name,n,family_name,type_name):
    #takes genome name (ex: "T_tenax.fa.out"), length specification, and keyword (ex: "Line") as an input
    #what do you want to filter from the .fa.out file
    #exports .bed file named "1000.bed"

    f = open(genome_name+ ".fa.out", "r")

    for lines in f:
        x = lines.split()
        d = int(x[6]) - int(x[5])
        
        if (d > n) and (family_name in x[9]) and (type_name in x[10]):
            with open("1000.bed","a") as bedfile:
                bedfile.write(x[4] + "\t")
                bedfile.write(x[5] + "\t")
                bedfile.write(x[6] + "\t")
                bedfile.write(x[9] + "\t")
                bedfile.write(x[10] + "\t")
                bedfile.write("\n")
    bedfile.close()

def genome(genome_name):
    #takes genome name (ex: "T_tenax.fa") as an input
    #exports .bed file named "1000_2.bed"

    f = open(genome_name + ".fa", "r")
    n = 1

    for lines in f:
        with open("1000_2.bed","a") as bedfile:
            if "LINE" in lines:
                bedfile.write(lines)
                n*=-1
                continue
            if n < 0:
                n*=-1
                bedfile.write(lines)
                bedfile.write("\n")
        
    bedfile.close()

def categorize():
    #opens "1000.bed" file, and produces a "category_1000.bed" file 
    #"category_1000.bed" file will include repeat/class family types, and how many of those were each found

    f = open("1000.bed", "r")
    l = {}

    for lines in f:
        x = lines.split()
        
        if x[3] not in l.keys():
            l[x[3]] = 1
        else:
            l[x[3]] += 1
        
    with open("category_1000.bed","a") as bedfile:
        for keys in l.keys():
            bedfile.write(keys + "\t" + str(l[keys]))
            bedfile.write("\n")
        
    bedfile.close()

def sixpack_primer(file_name):
    #takes one .fa file as an input, and outputs multiple .fa files (each containing one fasta line)
    #sixpack_primer would prime the .fa file for sixpack -sequence -output
    
    #first open .fa file in "r" mode
    f = open(file_name + ".fa", "r")
    file_list = []
    #unique .fa file name attribute
    a = 0

    for lines in f:
        file_list.append(lines)

    i = 0
    while i < len(file_list):
        with open("sequence_" + str(a),"a") as primer:
            primer.write(file_list[i])
            primer.write(file_list[i+1])
        primer.close()
        
        a+=1
        i+=2          

def longest_orf(file_name):    
    #WARNING: uses "aa" in .seq file as separator
    with  open(file_name + ".seq") as fp:
        contents = fp.read()
    fp.close()
    
    l = []
    
    for entry in contents.split('>'):
        #do something with entry  
        #entry = 8345757-8347785_6_ORF70  Translation of 8345757-8347785 in frame 6, ORF 70, threshold 1, 31aa ///// LFEVHLNVHRLHLYQVHNRKAFWFALLCYSX
        l.append(entry)
    
    l.pop(0)

    l_dict = {}
    
    for line in l:
        x = line.split("aa")
        l_dict[x[0] + "aa"] = x[1]

    longest_orf = max(l_dict.values(), key=len)
    
    return longest_orf

