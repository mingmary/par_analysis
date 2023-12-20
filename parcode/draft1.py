"""
usage: python3 draft1.py
---in progress---
"""

import re
import subprocess

gap1file = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_filtered.gap1.name.sorted.sam"
homofile = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_homo.name.sorted.sam"
rrifile = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_rri.name.sorted.sam"
outputfilesam = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_total.sam"
outputfilebed = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_total.bed"

"""
def gap1torri - transforms gap1_filtered sorted SAM file to rri format: each read is repeated on two lines where each line represents one of the two interacting fragments.
logic:
temporary output SAM file is created
gap1_filtered file is filtered line by line so that:
- header is ommited
- reads are filtered for mapq 255
- reads that are not recognized by STAR as chimeric per read are processed as following:
	- One of max. two possible alternatives is selected. Alternatives are mapped to the same genomic region but have minor differences in cigar strings
	- Cigar string is processed so that it is divided by two strings with one match only. Example: '6S20M22N25M' -> ['20M', 35M']. D/I are also taken in consideration
	- A position for the second record is created as follows: POS + X, where X - number of bases in cigar strings before the second match 
	- Two modified lines derived from 1 original line per read are written in the temporary output file
- reads that are recognized by STAR as chimeric are processed as following:
	- The original line is written to the output file with match only in cihar and second full cigar pair
	- The coordinates in the original line are changed to the coordinates from SA:Z: tag, match is taken from chimeric cigar. original pos and cigar are written into SA:Z tag
- gap1 tag is added

def homotorri - transforms homo sorted SAM file to rri format
logic:
temporary SAM file is created
homo file is filtered line by line so that:
- header is ommited
- reads (recognized by STAR as chimeric only) are processed as following:
	- The original line is written to the output file 
- homo tag is added

def rriaddtag - adds "rri" tag to all the lines of the original sorted SAM file
logic:
temporary SAM file is created
rri file is being processed line by line so that:
- header is ommited
- reads are rewritten as in the original file but the "rri" tag is added to the end of each line 

"""

def createsam(outputfilesam):
	outputfilesam = open(outputfilesam, "w")

def gap1torrinew (gap1file):
	mapqdict = {"2048":"0", "0":"2048", "16":"2064", "2064":"16"}
	pos2 = "dummy" #auxillary variable
	curname = "dummy" #auxillary variable
	cigarpair1 = ""
	cigarpair2 = ""
	gap1reformed = open("tmp.gap1.sam", "w") #temporary gap1 file
	outpfile = open(outputfilesam, "a")
	with open(gap1file) as gap1:
		for line in gap1:
			if line.startswith("@"): #keeping the header
				pass
			else:
				curline = line.split("\t")
				if curline[4] == '255': #mapq filter
					if len(curline) == 19: #19 for standart gap1 or 21 for chimeric
						if curline[0] != curname: #choosing the first occurence out of 2 alternatives
							curname = curline[0]
							cigar = curline[5]
							cigarlist = re.findall(r'\d+[A-Z]', cigar) #string to list
							pos1 = curline[3]
							cigarletterslist = [element[-1] for element in cigarlist]
							Mindex = cigarletterslist.index('M')
							Nindex = cigarletterslist.index('N')
							cigarpair1 = ""	
							while Mindex < Nindex:
								cigarpair1 = cigarpair1 + cigarlist[Mindex]
								Mindex += 1
							lastNindex = len(cigarletterslist) - 1 - cigarletterslist[::-1].index('N')
							if cigarletterslist[-1] != "M":
								newcigarlist = cigarlist[(lastNindex+1):-1]
							else:
								newcigarlist = cigarlist[(lastNindex+1):]	
							cigarpair2 = ''.join(map(str, newcigarlist))
							if cigarletterslist[0] != "M":
								cigarforpos2 = cigarlist[1:(lastNindex+1)]
							else:
								cigarforpos2 = cigarlist[:(lastNindex+1)]	
							if curline[1] == "0":
								sign  = "+"
							elif curline[1] == "16":
								sign = "-"
							pos2 = int(pos1)
							
							for el in cigarforpos2:
								if el[:-1] != "I":
									pos2 = pos2 + int(el[:-1])
							curline[5] = cigarpair1
							gap1reformed.write(('\t'.join(map(str, curline)))[:-1] + "\tch:A:1\tSA:Z:" + curline[2] + "," + str(pos2) + "," + sign + "," + cigar + "," + (curline[15])[-1] + ";\tCO:Z:gap1\n")
							outpfile.write(('\t'.join(map(str, curline)))[:-1] + "\tch:A:1\tSA:Z:" + curline[2] + "," + str(pos2) + "," + sign + "," + cigar + "," + (curline[15])[-1] + ";\tCO:Z:gap1\n")
							curline[5] = cigarpair2
							curline[3] = pos2
							curline[1] = mapqdict[(curline[1])]
							gap1reformed.write(('\t'.join(map(str, curline)))[:-1] + "\tch:A:1\tSA:Z:" + curline[2] + "," + str(pos1) + "," + sign + "," + cigar + "," + (curline[15])[-1] + ";\tCO:Z:gap1\n")
							outpfile.write(('\t'.join(map(str, curline)))[:-1] + "\tch:A:1\tSA:Z:" + curline[2] + "," + str(pos1) + "," + sign + "," + cigar + "," + (curline[15])[-1] + ";\tCO:Z:gap1\n")
					elif len(curline) == 21: #19 for standart gap1 or 21 for chimeric
						cigar = curline[5]
						cigarlist = re.findall(r'\d+[A-Z]', cigar) #string to list
						cigarletterslist = [element[-1] for element in cigarlist]
						Mindex = cigarletterslist.index('M')
						Nindex = cigarletterslist.index('N')
						cigarpair1 = ""	
						while Mindex < Nindex:
							cigarpair1 = cigarpair1 + cigarlist[Mindex]
							Mindex += 1
						if curline[1] == "0" or curline[1] == "2048":
							sign = "+"
						elif curline[1] == "16" or curline[1] == "2064":
							sign = "-"
						chimera = curline[-1].split(",")
						pos2 = chimera[1]
						cigar2 = chimera[3]
						nm = ((chimera[5]).split(";"))[0]
						curline[5] = cigarpair1
						gap1reformed.write(('\t'.join(map(str, curline)))[:-1] + "\tCO:Z:gap1\n")
						outpfile.write(('\t'.join(map(str, curline)))[:-1] + "\tCO:Z:gap1\n")
						satag = "SA:Z" + curline[2] + "," + curline[3] + "," + sign + "," + cigar + "," + curline[4] + "," + ((curline[15]).split(":"))[-1] + ";"
						cigarlistchim2 = re.findall(r'\d+[A-Z]', cigar2)
						matcheslist = [index for index, element in enumerate(cigarlistchim2) if element[-1] == 'M']
						cigarpairchim = ""
						cont = 0
						if len(matcheslist) == 1:
							cigarpairchim = cigarlistchim2[matcheslist[0]]
						else:
							mind = matcheslist[0]
							while mind < matcheslist[-1] + 1:	#this allows to take in consideration I/D in cigar
								if (cigarlistchim2[mind])[-1] != "I":
									cont = cont + int((cigarlistchim2[mind])[:-1])
								elif (cigarlistchim2[mind])[-1] == "I":
									cont = cont - int((cigarlistchim2[mind])[:-1])	
								if (cigarlistchim2[mind])[-1] == "N":
									cigarpairchim = ""
									pos2 = int(pos2) + int((cigarlistchim2[mind])[:-1]) + cont
								else:
									cigarpairchim = cigarpairchim + cigarlistchim2[mind]
								mind += 1		
						curline[3] = str(pos2)
						curline[5] = cigarpairchim
						curline[1] = mapqdict[(curline[1])]
						curline[14] = curline[14] + nm
						curline[15] = curline[15] + nm
						curline[-1] = satag
						gap1reformed.write(('\t'.join(map(str, curline))) + "\tgap1\n")
						outpfile.write(('\t'.join(map(str, curline))) + "\tgap1\n")
						

def homotorri(homofile):
	homoreformed = open("tmp.homo.sam", "w") #temporary homo file
	outpfile = open(outputfilesam, "a")
	with open(homofile) as homo:
		for line in homo:
			if line.startswith("@"): #keeping the header
				pass
			else:
				homoreformed.write(line[:-1] + "\tCO:Z:homo\n")
				outpfile.write(line[:-1] + "\tCO:Z:homo\n")


def rriaddtag(rrifile):		
	outpfile = open(outputfilesam, "a")
	rrireformed = open("tmp.rri.sam", "w") #temporary homo file
	with open(rrifile) as rri:
		for line in rri:
			if line.startswith("@"): #keeping the header
				pass
			else:
				rrireformed.write(line[:-1] + "\tCO:Z:rri\n")
				outpfile.write(line[:-1] + "\tCO:Z:rri\n")

def createbed(inputsamfile, outputbedfile):
	outpfile = open(outputbedfile, "w")
	outpfile.write("chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tcigar\n")
	with open(inputsamfile, "r") as samfile:
		for line in samfile:
			alignment = line.split("\t")
			chrom = alignment[2]
			chromStart = alignment[3]
			cigarlist = re.findall(r'\d+[A-Z]', alignment[5])
			cont = sum([int(element) for element in ([element[:-1] for element in cigarlist])])
			#chromEnd = 0
			#strand = "dummy"
			if alignment[1] == "0" or alignment[1] == "2048":
				strand = "+"
			elif alignment[1] == "16" or alignment[1] == "2064":
				strand = "-"
			chromEnd = int(chromStart) + cont	
			name = alignment[0]	
			#print (alignment[-2])
			cigar = (alignment[-2].split(","))[3]			
			outpfile.write(f"{chrom}\t{chromStart}\t{chromEnd}\t{name}\t255\t{strand}\t{cigar}\n")

createsam(outputfilesam)
gap1torrinew(gap1file)	
homotorri(homofile)
rriaddtag(rrifile)
createbed(outputfilesam, outputfilebed)
