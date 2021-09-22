f = open("1000.fa", "r")
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