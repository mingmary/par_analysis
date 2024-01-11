"""
usage: python3 draft1.py
---in progress---
"""

import re
import sys
import subprocess

gap1file = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_filtered.gap1.name.sorted.sam"
homofile = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_homo.name.sorted.sam"
rrifile = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_rri.name.sorted.sam"
outputfilesam = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_total.sam"
outputfilebed = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_total.bed"
annogtffile = "/data_01/mmintsev/ref/table5pENSTpCAT.gtf"
annobedfile = "/data_01/mmintsev/ref/table5pENSTpCAT.bed"
intersectedbed =  "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_wo.annotation.bed"
selectedintersected = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_gaptypes_selected.annotation.bed"
stats = "/data_01/mmintsev/paris/iPSC/Set18-8_iPSC_rep1_CGTGAT_interactions.stats.txt"
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
						gap1reformed.write(('\t'.join(map(str, curline))) + "\tCO:Z:gap1\n")
						outpfile.write(('\t'.join(map(str, curline))) + "\tCO:Z:gap1\n")
						

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
	outpfile.write("chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tcigar\torigin\n")
	with open(inputsamfile, "r") as samfile:
		for line in samfile:
			alignment = line.split("\t")
			chrom = alignment[2]
			chromStart = alignment[3]
			cigarlist = re.findall(r'\d+[A-Z]', alignment[5])
			cont = sum([int(element) for element in ([element[:-1] for element in cigarlist])])
			#chromEnd = 0
			#strand = ("dummy"
			if alignment[1] == "0" or alignment[1] == "2048":
				strand = "+"
			elif alignment[1] == "16" or alignment[1] == "2064":
				strand = "-"
			chromEnd = int(chromStart) + cont	
			name = alignment[0]	
			#print (alignment[-2])
			cigar = (alignment[-2].split(","))[3]
			origin = (alignment[-1].split(":"))[2]
			print(origin)
			outpfile.write(f"{chrom}\t{chromStart}\t{chromEnd}\t{name}\t255\t{strand}\t{cigar}\t{origin}")

def gtftobed(inputgtffile, outputbedfile):
	featuresvocabulary = {}
	transdictionary = {}
	genedictionary = {}
	outpfile = open(outputbedfile, "w")
	outpfile.write("chrom\tchromStart\tchromEnd\tname\tscore\tstrand\ttranscript_id\ttranscript_type\tgene_name\tgene_id\texon_number\ttranscript_novelty\tgene_novelty\tsource\n")
	with open(inputgtffile, "r") as gtffile:
		for line in gtffile:
			feature = line.split("\t")
			ent = feature[2]
			if ent == "gene":
				gene_id = (((((feature[-1]).split(";"))[0]).split(" "))[-1]).replace('"', '')
				novelty = (((((feature[-1])[:-2]).split(";"))[-1]).replace('"', '')).split(" ")
				if novelty[-2] == "gene_novelty":
					genedictionary[gene_id] = novelty[-1]
				else:	
					genedictionary[gene_id] = "novel"
	with open(inputgtffile, "r") as gtffile:
		for line in gtffile:
			feature = line.split("\t")
			ent = feature[2]
			if ent == "transcript":
				transcript_id = (((((feature[-1]).split(";"))[1]).split(" "))[-1]).replace('"', '')
				novelty = (((((feature[-1])[:-2]).split(";"))[-1]).replace('"', '')).split(" ")
				if novelty[-2] == "transcript_novelty":
					transdictionary[transcript_id] = novelty[-1]
				else:	
					transdictionary[transcript_id] = "novel"
			if ent == "exon":
				chrom = feature[0]
				chromStart = feature[3]
				chromEnd = feature[4]
				strand = feature[6]
				annofeatures = ((feature[-1])[:-2]).split(";")
				for feat in annofeatures:
					feat = feat.replace('"', '')
					featlist = feat.split(" ")
					featuresvocabulary[(featlist[-2])] = featlist[-1]
				transcript_id = featuresvocabulary['transcript_id']
				source = feature[1]
				try:
					transcript_type = featuresvocabulary['transcript_type']
				except KeyError:
					transcript_type = "NA"
				try:
					gene_name = featuresvocabulary['gene_name']
				except KeyError:
					gene_name = "NA"
				gene_id = featuresvocabulary['gene_id']
				try:
					exon_number = featuresvocabulary['exon_number']
				except KeyError:
					exon_number = "NA"	
				transcript_novelty = transdictionary[transcript_id]
				gene_novelty = genedictionary[gene_id]	
				outpfile.write(f"{chrom}\t{chromStart}\t{chromEnd}\t.\t.\t{strand}\t{transcript_id}\t{transcript_type}\t{gene_name}\t{gene_id}\t{exon_number}\t{transcript_novelty}\t{gene_novelty}\t{source}\n")
				featuresvocabulary = {}

def prioritisation(dummylist, count, typeslist, trnoveltylist):
	while count > 1:	
		if (dummylist[-1])[15] != (dummylist[-2])[15]:
			if typeslist[((dummylist[-1])[15])] > typeslist[((dummylist[-2])[15])]:
				dummylist.pop(-1)
			else:
				dummylist.pop(-2)	
		elif (dummylist[-1])[19] != (dummylist[-2])[19]:
			if trnoveltylist[((dummylist[-1])[19])] > trnoveltylist[((dummylist[-2])[19])]:
				dummylist.pop(-1)
			else:
				dummylist.pop(-2)
		elif (dummylist[-1])[-1] != (dummylist[-2])[-1]:
			if (dummylist[-1])[-1] > (dummylist[-2])[-1]:
				dummylist.pop(-2)
			else:
				dummylist.pop(-1)
		elif (dummylist[-1])[-9] != (dummylist[-2])[-9]:
			if (dummylist[-1])[-9] > (dummylist[-2])[-9]:
				dummylist.pop(-2)
			else:
				dummylist.pop(-1)
		else:
			dummylist.pop(-1)	
		count = count - 1
	return dummylist[0]	

def selecttranscripts(inputbed, outputbed, statsfile):
	output = open(outputbed, "w")
	stats = open(statsfile, "w")
	#input = open(inputbed, "r")
	#lines = input.readlines()
	count = 0
	read = ""
	chr = ""
	start = ""
	strand = ""
	dummy = False
	readcount = False
	dummylist = []
	typeslist = {
        "IG_C_gene": 3,
        "IG_C_pseudogene" : 3,
        "IG_D_gene" : 3,
        "IG_J_gene" : 3,
        "IG_J_pseudogene" : 3,
        "IG_pseudogene" : 3,
        "IG_V_gene" : 3,
        "IG_V_pseudogene" : 3,
        "lncRNA" : 3,
        "miRNA" : 3,
        "misc_RNA" : 3,
        "Mt_rRNA" : 3,
        "Mt_tRNA" : 1,
        "NA" : 4,
        "ncRNA" : 3,
        "nonsense_mediated_decay" : 3,
        "non_stop_decay" : 3,
        "others" : 4,
        "polymorphic_pseudogene" : 3,
        "processed_pseudogene" : 3,
        "processed_transcript" : 3,
        "protein_coding" : 2,
        "protein_coding_CDS_not_defined" : 3,
        "protein_coding_LoF" : 3,
        "pseudogene" : 3,
        "retained_intron" : 3,
        "ribozyme" : 3,
        "rRNA" : 3,
        "rRNA_pseudogene" : 3,
        "scaRNA" : 3,
        "scRNA" : 3,
        "snoRNA" : 1,
        "snRNA" : 3,
        "sRNA" : 3,
        "TEC" : 3,
        "transcribed_processed_pseudogene" : 3,
        "transcribed_unitary_pseudogene" : 3,
        "transcribed_unprocessed_pseudogene" : 3,
        "translated_processed_pseudogene" : 3,
        "translated_unprocessed_pseudogene" : 3,
        "TR_C_gene" : 3,
        "TR_D_gene" : 3,
        "TR_J_gene" : 3,
        "TR_V_gene" : 3,
        "TR_J_pseudogene" : 3,
        "TR_V_pseudogene" : 3,
        "tRNA" : 1,
        "TR_V_gene" : 3,
        "TR_V_pseudogene" : 3,
        "unitary_pseudogene" : 3,
        "unprocessed_pseudogene" : 3,
        "vault_RNA" : 3,
    }

	trnoveltylist = {
		"novel" : 3,
		"GENCODEv39" : 1,
		"GENCODE_updated" : 2,
	}
	dummylist1 = []
	with open(inputbed, "r") as file0:
		for line in file0:
			ents = line.split("\t")
			if ents[7] == "gap1":
				if ents[3] == read and ents[1] == start:
					dummylist.append(ents)
				else:					
					if read == ents[3]:
						readcount = True
						dummylist1 = dummylist
					else:
						read = ents[3]
						if readcount == True:
							ent1 = prioritisation(dummylist1, len(dummylist1), typeslist, trnoveltylist)
							ent2 = prioritisation(dummylist, len(dummylist), typeslist, trnoveltylist)
							output.write('\t'.join(ent1))
							output.write('\t'.join(ent2))
							stats.write(f"{ent1[14]}\t{ent2[14]}\t{max((ent1[14]), (ent2[14]))}+{min((ent1[14]), (ent2[14]))}\t{ent1[15]}\t{ent2[15]}\t{max((ent1[15]), (ent2[15]))}+{min((ent1[15]), (ent2[15]))}\tgap1\n")
						readcount = False
					dummylist = []
					start = ents[1]
					dummylist.append(ents)				
			if ents[7] == "homo":
				if ents[3] == read:
					dummylist.append(ents)
				else:
					read = ents[3]
					ent1 = prioritisation(dummylist, len(dummylist), typeslist, trnoveltylist)
					output.write('\t'.join(ent1))
					stats.write(f"{ent1[14]}\t{ent1[14]}\t{ent1[14]}+{ent1[14]}\t{ent1[15]}\t{ent1[15]}\t{ent1[15]}+{ent1[15]}\thomo\n")
					dummylist = []
					dummylist.append(ents)
			if ents[7] == "rri":	
				if ents[3] == read and ents[0] == chr and ents[5] == strand:
					dummylist.append(ents)
				else:
					if read == ents[3]:
						readcount = True
						dummylist1 = dummylist
					else:
						read = ents[3]
						if readcount == True:
							ent1 = prioritisation(dummylist1, len(dummylist1), typeslist, trnoveltylist)
							ent2 = prioritisation(dummylist, len(dummylist), typeslist, trnoveltylist)
							output.write('\t'.join(ent1))
							output.write('\t'.join(ent2))
							stats.write(f"{ent1[14]}\t{ent2[14]}\t{max((ent1[14]), (ent2[14]))}+{min((ent1[14]), (ent2[14]))}\t{ent1[15]}\t{ent2[15]}\t{max((ent1[15]), (ent2[15]))}+{min((ent1[15]), (ent2[15]))}\trri\n")
						readcount = False
					dummylist = []
					chr = ents[0]
					strand = ents[5]
					dummylist.append(ents)	
					





#createsam(outputfilesam)
#gap1torrinew(gap1file)	
#homotorri(homofile)
#rriaddtag(rrifile)
#createbed(outputfilesam, outputfilebed)
#gtftobed(annogtffile, annobedfile)
selecttranscripts(intersectedbed, selectedintersected, stats)