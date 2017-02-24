#!/usr/bin/python

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
# WHAT: Parse fasta nucleotide files for microbial genomes
# WHY: Programming Midterm Assignment
# HOW:
# WHO: agm/2016 provided assignment and base code, jhr/2016 modified code for asignment

import re, sys
from collections import defaultdict

# - - - - - U S E R    V A R I A B L E S - - - - - - - -
# NOTE: 
# The location of the input file is defined by the inFolder variable.
# In the class examples, there is a folder 01-Genomes located in the same
#     folder as the script (i.e. the current working directory).
# So the path for the interpreter to find the file is:
#             /01-Genomes/AcinetoBaum-uid59271.ffn
inGenome = "AcinetoBaum-uid59271.ffn"
inFolder = "01-Genomes"

# - - - - - G L O B A L  D E C L A R A T I O N S  - - - - - -
inFILE = open(inFolder + "/" + inGenome, 'r')
GeneNukes = defaultdict(lambda: defaultdict(lambda: 'unk' ))
fourNucleotides = ['A', 'G', 'T', 'C']

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\nStarting to process genome file: ", inGenome

# - - --  - -----   - ---- --- -- -- ------- -  - -- - -  - ---  - --- -- ---- 
# TASK 1 . . . . . . . 
# Input fasta seqs from file . . . . . . . . . . . 
#print "     1",
seq = ''
for line in inFILE:
	if re.match(r'^>', line):
		if len(seq) > 1:
		# Save current working sequence . . . . 
			GeneNukes[geneid]['seq'] = seq
		# Begin processing new sequences . . . . .
		# HEADER:  >gi|215481761|ref|NC_011595.1|:95-1492 Acinetobacter baumannii AB307-0294, complete genome
		# TARGET:                               |:95-1492 Acinetobacter baumannii AB307-0294,
		head = re.search(r'\|\:([\w\d\-]+) ([\w \d\-]+),', line )
		seq = ''
		if head is not None:
			geneid = head.group(2)+"|pos:"+head.group(1)
		else:
			print line
			print "There was an error parsing header info . . . . \n"
			sys.exit()
	else:
		seq += line.rstrip()	
GeneNukes[geneid]['seq'] = seq

#print "There are %d of genes in the input file." % len(GeneNukes.keys())

# - - --  - -----   - ---- --- -- -- ------- -  - -- - -  - ---  - --- -- ---- 
# TASK 2 . . . . . . . 
# Calculate the Fraction %A+T composition of the genome as a mean across genes, and
#    keep track of the minimum and maximum values.
# NOTE: This task is complete. The code is provided as an example for the last task section.
# NOTE: The 'A' and the 'T' counts are stored in the GeneNukes dictionary.
# NOTE: The ATfrac value for each gene is stored in the GeneNukes dictionary.
#print "     2"
fracAVG  = 0                                    # variable to store the average
fracMIN  = 999                                  # variable to store the minimum
fracMAX  = -1                                   # variable to store the maximum
for gid in GeneNukes.iterkeys():                # iterate the unique Gene IDs using a variable named 'gid'
	gene = GeneNukes[gid]['seq']                # copy the NT sequence to a variable named 'gene' 
	N = len(gene)
	GeneNukes[gid]['length'] = N
	Acount = gene.count('A')                    # count A and save in Acount
	Tcount = gene.count('T')                    # count T and save in Tcount
	GeneNukes[gid]['A'] = Acount                # save Acount to GeneNukes dictionary
	GeneNukes[gid]['T'] = Tcount                # save Tcount to GeneNukes dictionary
	ATfrac = 100 * (Acount + Tcount)/float(N)   # calculate the A+T fraction as a percentage
	GeneNukes[gid]['A+T'] = ATfrac              # save each gene ATfrac value back into the GeneNukes dictionary
	fracAVG += ATfrac 							# calculate the gene mean A+T percentage
	if fracMIN > ATfrac: fracMIN = ATfrac       # find the ATfrac gene minimum
	if fracMAX < ATfrac: fracMAX = ATfrac       # find the ATfrac gene maximum

fracAT = fracAVG/len(GeneNukes.keys())			# I added this because fracAVG was a huge number that was total 
												#of all the ATfrac's
# - - --  - -----   - ---- --- -- -- ------- -  - -- - -  - ---  - --- -- ---- 
# TASK 3 . . . . . . . 
# Calculate the OBSERVED and EXPECTED frequencies of "AT" motifs across genes, and
#    keep track of the minimum and maximum values.
# a. For the OBSERVED calculation, use the count('AT') function to find the number of
#    'AT' motifs and then divide by the length of the gene sequence. Leave the result as
#    a fractional frequency (do not multiply by 100).
# b. For the EXPECTED calculation, use the A and T counts already stored in the GeneNukes dictionary
#    and divide by gene length to get fractional composition of each. Multiplying the frequency of
#    A times the frequency of T yields a product that is the random expectation for the frequency
#    of A and T appearing side by side:  (A/N)*(T/N). Leave this number as a frequency value (do not
#    multiply by 100) so that it is comparable to the OBSERVED value calculated above.
# c. Store the OBS and EXP values back into the GeneNukes dictionary.
# d. Find the MEAN, MIN and MAX values for both OBS and EXP calculations across all genes. 
#print "     3"
ATOBStotal = 0 
ATEXPtotal = 0

ATOBSMIN = 9999
ATEXPMIN = 9999
	
ATOBSMAX = -1
ATEXPMAX = -1

# count AT motifs . . . . . 
for gid in GeneNukes.iterkeys():  	
	
	#count # of AT motifs
	ATCount = GeneNukes[gid]['seq'].count('AT')
	
	#put AT count into dict
	GeneNukes[gid]['ATCount'] = ATCount
	
	#get gene length 
	N = len(gene) 
	
	#get AT frequency - #ATmotifs / gene length
	ATOBS = ATCount /float(N)	
	
	#put observed AT motif frequency into dict
	GeneNukes[gid]['ATOBS'] = ATOBS	
	
	#calculate expected AT motif frequency
	# formula is (A/N)*(T/N), broke this into three lines of code to make it simpler
	AEXP = (GeneNukes[gid]['A'])/float(GeneNukes[gid]['length'])
	TEXP = (GeneNukes[gid]['T'])/float(GeneNukes[gid]['length'])
	ATEXP = AEXP * float(TEXP)
	
	#put expected AT frequency into dict
	GeneNukes[gid]['ATEXP'] = ATEXP
		
	#AT mean calculations, need to get total of all AT frequencies
	#and then divide by number of total genes. 
	#Division done outside of the for loop
	ATOBStotal += ATOBS
	ATEXPtotal += ATEXP
	
	#AT min calculations
	if ATOBSMIN > ATOBS: ATOBSMIN = ATOBS
	if ATEXPMIN > ATEXP: ATEXPMIN = ATEXP
		
	# AT MAX calculations
	if ATOBSMAX < ATOBS: ATOBSMAX = ATOBS
	if ATEXPMAX < ATEXP: ATEXPMAX = ATEXP
	
	
#This is last step of mean calculations, have total and divide by number of genes	
ATOBSMEAN = ATOBStotal/len(GeneNukes.keys())
ATEXPMEAN = ATEXPtotal/len(GeneNukes.keys())	

# - - --  - -----   - ---- --- -- -- ------- -  - -- - -  - ---  - --- -- ---- 
# TASK 4 . . . . . . .
# Print to screen in the following format:
# !!--------------------			       < 2 exclaims, 20 hyphens
# Calculated Results:                             < title
#         mean OBS %A+T:    52.31 %                   < 8 spaces, var_string, colon, 4 spaces, %0.2f %%
# 		  min  OBS %A+T:    22.22 %
#         max  OBS %A+T:    66.66 %
#         mean OBS  ApT:    12.54 %                     < 8 spaces, var_string, colon, 4 spaces, %0.2f 
#         min  OBS  ApT:    4 %                         < 8 spaces, var_string, colon, 4 spaces, %d
#         max  OBS  ApT:    22 %
#         mean EXP  ApT:    12.54 %                     < 8 spaces, var_string, colon, 4 spaces, %0.2f 
#         min  EXP  ApT:    4 %                         < 8 spaces, var_string, colon, 4 spaces, %d
#         max  EXP  ApT:    22 %
# !!--------------------			       < 2 exclaims, 20 hyphens

print "!!--------------------"
print "Calculated Results:"  
print "        mean OBS %%A+T:    %0.2f %%" % fracAT
print "        min  OBS %%A+T:    %0.2f %%" % fracMIN
print "        max  OBS %%A+T:    %0.2f %%" % fracMAX
print "        mean OBS  ApT:    %0.4f" % ATOBSMEAN #used %0.4f instead of %0.2f because the values are small
print "        min  OBS  ApT:    %0.4f" % ATOBSMIN  #when I use %d it prints 0 instead of decimal 
print "        max  OBS  ApT:    %0.4f" % ATOBSMAX
print "        mean EXP  ApT:    %0.4f" % ATEXPMEAN
print "        min  EXP  ApT:    %0.4f" % ATEXPMIN
print "        max  EXP  ApT:    %0.4f" % ATEXPMAX
print "!!--------------------"


#print "     4"

# - - --  - -----   - ---- --- -- -- ------- -  - -- - -  - ---  - --- -- ---- 
# TASK 5 . . . . . . .
# Write to file a data tab-delimited data table that has three columns of data
#    for each gene in the genome:
#           geneID      ATobs    ATexp

OUT = open("77-Rhoades.Midterm.OUTPUT.txt", 'w')
for gid in GeneNukes.iterkeys():
	OUT.write("%s\t%0.4f\t%0.4f\n" % (gid, GeneNukes[gid]['ATOBS'], GeneNukes[gid]['ATEXP']  ) )

#print "     5"


print "\n\n* * * * *   D O N E   * * * * * \n\n"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - E n d   O f   F i l e - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


