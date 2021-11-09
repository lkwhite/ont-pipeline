import pandas as pd
from statistics import mean

'''
get_all_quals processes the output of samtools mpileup by splitting each line into a list with elements:
0. chromosome
1. position (a single nucleotide)
2. reference base at position
3. # of reads at position
4. "read bases": a string of characters describing the alignment at the position. Different notation
for matches, mismatches, indels, strand, mapping quality, and starts and ends of reads, not necessarily 1 char/read.
5. Base qualities, encoded as ASCII characters (1/ position).
6. Alignment mapping qualities, encoded as ASCII characters.
'''

def get_all_quals(file1, file2):
    with open(file1, 'r') as in_pileup, open(file2, 'w') as out_pileup:

        for line in in_pileup:
            line = line.rstrip("\n") # pull that special character off the end
            pileup_data=line.split("\t") # split the string up by tabs, obviously
            base_qual_holder = "" # we're going to hold all the numeric quality values in this string
            align_qual_holder = "" # we're going to hold all the numeric quality values in this string

            for value in pileup_data[5]:
                base_qual_holder += str(ord(value)-33) + ","
            
            for value in pileup_data[6]:
                align_qual_holder += str(ord(value)-33) + ","

            base_qual_holder = base_qual_holder.rstrip(",") # remove the final comma from the string of quality values
            align_qual_holder = align_qual_holder.rstrip(",") # remove the final comma from the string of quality values
            out_pileup.write("\t".join(pileup_data[0:4]) + "\t" + base_qual_holder + "\t" + align_qual_holder + "\n")
            
'''
get_qual_means takes a file (such as the processed pileup output of get_all_quals)
opens the file and adds the appropriate header
splits the base and alignment quality scores into lists using commas as the split value
converts each object in the lists into an integer
and adds new columns containing single values, the result of doing some math on the integers in each list
    for the moment, that math is just a mean
    but taking the arithmetic mean of qscores is generally bad practice because they're logrithmic
and then writes out the data frame as a new tab separated file for downstream processing in R
''' 
def get_qual_means(file1, file2):
    with open(file1, 'r') as in_file, open(file2, 'w') as out_file:
        header_list = ["chrom", "pos", "ref_base", "read_num", "base_quals", "align_quals"]
        df = pd.read_csv(in_file, sep='\t', names = header_list)
        df["base_quals"] = df["base_quals"].str.split(",").apply(lambda x: [int(i) for i in x])
        df["base_qual_mean"] = df.apply(lambda x: mean(x.base_quals), axis = 1)
        df["align_quals"] = df["align_quals"].str.split(",").apply(lambda x: [int(i) for i in x])
        df["align_qual_mean"] = df.apply(lambda x: mean(x.align_quals), axis = 1)
        df.to_csv(out_file, sep='\t')
