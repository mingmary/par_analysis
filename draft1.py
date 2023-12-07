"""
usage: python3 draft1.py
---in progress---
"""

import re

gap1file = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_filtered.gap1.sorted.sam"
outputfile = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep2_ACATCG_gaptypes_total.sam"

"""
def gap1torri - transforms gap1_filtered sorted SAM file to rri format: each read is repeated on two lines where each line represents one of the two interacting fragments.
logic:
temporary output SAM file is created
gap1_filtered file is filtered line by line so that:
- header is written to the output file 
- reads are filtered for mapq 255
- reads with more than 2 matches in cigar pairs are filtered (<1%)
- reads that are not recognized by STAR as chimeric per read are processed as following:
	- One of max. two possible alternatives is selected. Alternatives are mapped to the same genomic region but have minor differences in cigar strings
	- Cigar string is processed so that it is divided by two strings with one match only. Example: '6S20M22N25M' -> ['6S20M47N', 48N35M']
	- A position for the second record is created as follows: POS + X, where X - number of bases in cigar strings before the second match 
	- Two modified lines derived from 1 original line per read are written in the temporary output file
- reads that are recognized by STAR as chimeric are processed as following:
	- The original line is written to the output file without chimeric tags
	- The coordinates in the original line are changed to the coordinates with the coordinates from SA:Z: tag
	- Modified line is written in the output file without 
"""

def gap1torri (gap1file):
	cigar1 = "dummy" #auxillary variable
	cigar2 = "dummy" #auxillary variable
	pos2 = "dummy" #auxillary variable
	curname = "dummy" #auxillary variable
	gap1reformed = open("tmp.gap1.sam", "w") #temporary gap1 file
	with open(gap1file) as gap1:
		for line in gap1:
			if line.startswith("@"): #keeping the header
				gap1reformed.write(line)
			else:
				curline = line.split("\t")
				if curline[4] == '255': #mapq filter
					if curline[5].count('M') == 2: #filtering artefacts from gapm (<1%)
						if len(curline) == 19: #19 for standart gap1 or 21 for chimeric
							if curline[0] != curname: #choosing the first occurence out of 2 alternatives
								curname == curline[0]
								cigar = curline[5]
								cigarlist = re.findall(r'\d+[A-Z]', cigar) #string to list
								pos1 = curline[3]
								matcheslist = [index for index, element in enumerate(cigarlist) if element[-1] == 'M']
								if (matcheslist)[0] == 0: #cigar1 creation
									cigar1 = cigarlist[0]
									nonmatches = sum([int(element) for element in ([element[:-1] for element in cigarlist[1:]])])
									cigar1 = cigar1 + str(nonmatches) + "N" #is N a correct letter to put???
								else:
									befmatch = cigarlist[:(matcheslist[0] + 1)]
									aftmatch = cigarlist[(matcheslist[0] + 1):]
									aftsum = sum([int(element) for element in ([element[:-1] for element in aftmatch])])
									cigar1 = ''.join(map(str, befmatch)) + str(aftsum) + "N"
								befomatch = cigarlist[:(matcheslist[1])] #cigar2 creation
								befosum = sum([int(element) for element in ([element[:-1] for element in befomatch])])
								pos2 = int(pos1) + befosum
								if matcheslist[1] != (len(cigarlist) - 1):
									cigar2 = str(befosum) + "N" + ''.join(map(str,  cigarlist[(matcheslist[1]):]))
								else:	
									cigar2 = str(befosum) + "N" + cigarlist[-1]
								curline[5] = cigar1	
								gap1reformed.write('\t'.join(map(str, curline))) #writing two lines with non-overlapping cigar strings for further annotation
								curline[5] = cigar2
								curline[3] = str(pos2)
								gap1reformed.write('\t'.join(map(str, curline)))
						elif len(curline) == 21: #19 for standart gap1 or 21 for chimeric
							gap1reformed.write('\t'.join(map(str, curline[:-2])) + "\n")
							chimera = (curline[-1]).split(",") #SA:Z:chr1,179085836,+,7S29M48H,255,0;
							curline[2] = ((chimera[0]).split(":"))[-1] #chromosome
							curline[3] = chimera[1] #position
							curline[5] = chimera[3] #cigar
							curline[15] = "NM:i:" + (chimera[-1])[:-2]#NM
							curline[14] = "nM:i:" + (chimera[-1])[:-2]#nm
							gap1reformed.write('\t'.join(map(str, curline[:-2])) +"\n")
						else: 
							print(curline)	

gap1torri(gap1file)							
