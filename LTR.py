f = open("T_tenax.fa.out", "r")

for lines in f:
    x = lines.split()
    d = int(x[6]) - int(x[5])
    
    if (d > 1000) and ("Unknown" not in x[10]):
        with open("1000.bed","a") as bedfile:
            bedfile.write(x[4] + "\t")
            bedfile.write(x[5] + "\t")
            bedfile.write(x[6] + "\t")
            bedfile.write(x[10] + "\t")
            bedfile.write("\n")
bedfile.close()

