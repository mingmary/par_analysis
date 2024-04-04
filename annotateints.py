import re
import subprocess

scratchdir = "/scratch/mmintsev/paris2/"
#scratchdir = "/net/valis-01/scratch/mmintsev/paris2/"

name = "iPSC/Set18-8_iPSC_rep1_CGTGAT_"
#name = "iPSC/Set18-8_iPSC_rep2_ACATCG_"
#name = "iPSC/Set18-8_iPSC-Cyto_rep0_TCAAGT_"
#name = "iPSC/Set18-8_iPSC-Nuc_rep0_GATCTG_"
#name = "NSC/Set18-8_NSC_rep1_GCCTAA_"
#name = "NSC/Set18-8_NSC_rep2_TGGTCA_"
#name = "Neurons/Set18-8_Neuron_rep1_CACTGT_"
#name = "Neurons/Set18-8_Neuron_rep2_ATTGGC_"

#tpm enabled
tpmstatus = True
#Check tpmfile
tpmfile = "/data_01/mmintsev/raw/rnaseq/ipsc.quant.sf"
gap1file = "/data_01/mmintsev/paris/" + name + "gaptypes_filtered.gap1.name.sorted.sam"
gapmfile = "/data_01/mmintsev/paris/" + name + "gaptypes_filtered.gapm.name.sorted.sam"
homofile = "/data_01/mmintsev/paris/" + name + "gaptypes_homo.name.sorted.sam"
rrifile = "/data_01/mmintsev/paris/" + name + "gaptypes_rri.name.sorted.sam"
outputfilesam = "/data_01/mmintsev/paris/" + name + "gaptypes_total.sam"
outputfilebed = "/data_01/mmintsev/paris/" + name + "gaptypes_total.bed"
intersectedbed =  scratchdir + name + "gaptypes_wo.annotation.bed"
selectedintersected = scratchdir + name + "gaptypes_selected.annotation.bed"
stats = "/data_01/mmintsev/paris/" + name + "interactions.stats.txt"
filewithseq = "/data_01/mmintsev/paris/" + name + "gaptypes_selected.annotation.seq.bed"
uniqfile = scratchdir + name + "uniq.interactions.stats.txt"
countsfile = scratchdir + name + "interactions.types.stats.txt"
counts2 = scratchdir + name + "interactions.counts.txt"
counts3 = scratchdir + name + "interactions.intermed.counts.txt"
stats2 = scratchdir + name + "interactions.intermed.stats.txt"
fastqfile = "/data_01/mmintsev/raw/paris_" + name + "L1+L2comb_trim3.fastq"

annogtffile = "/data_01/mmintsev/ref/table5pENSTpCAT.gtf"
annobedfile = "/data_01/mmintsev/ref/table5pENSTpCAT.bed"

print(name)

#The function creates outputfile for all the gaptypes
def createsam(outputfilesam):
	outputfilesam = open(outputfilesam, "w")
	print("sam created")

#The function processes gapm file to convert it into rri format. The alignments with multiple gaps longer than 5 bp are classified as gapm, with one gap longer than 5 bp are classified as gap1. Other alignments are discarded.
#The first itiration through the file aims to find the maximum score for each read (important for multimappers)
#As the result of the function each match part of the alignment is written in the separate line. The original cigar string is written to SA:Z tag. 
#The output is written into outputfilesam
def gapmtorri(gapmfile):
	print("gapmstarted...")
	Nlist = []
	fakegap1 = 0
	truegapm = 0
	failed = 0
	curread = ""
	linecount = 0
	linestowrite = []
	newcigarlist =[]
	multimapcount = 0
	scorelist = []
	scoredict = {}
	name_id = ""
	output = open("tmpgapm.txt", "w")
	outpfile = open(outputfilesam, "a")
	with open(gapmfile) as gapm:
		for liney in gapm:
			curline = liney.split("\t")
			if liney.startswith("@"):
				continue
			if curline[0] == name_id:
				scorelist.append(int(((curline[13]).split(":"))[-1]))
				scoredict[name_id] = max(scorelist)
			else:
				if name_id != "":
					scoredict[name_id] = max(scorelist)
				name_id = curline[0]
				scorelist = []
				scorelist.append(int(((curline[13]).split(":"))[-1]))
				scoredict[name_id] = max(scorelist)	
		print("max assigned")		
	with open(gapmfile) as gapm:					
		for line in gapm:
			failedif = False
			Nlist = []
			if line.startswith("@"):
				pass
			else:
				curline = line.split("\t")
				if int(((curline[13]).split(":"))[-1]) != scoredict[(curline[0])]:
					continue #maxscore filtering
				if len(curline) == 21 and ((((curline[-1]).split(":"))[2]).split(","))[0] != curline[2]:
					continue #in case of gapm alignments marked with chimeric tag SA:Z are very few and they are incorrect (cigar pairs show they have alignment overlap)
				if curline[0] != curread:						
					curread = curline[0]
					for line in linestowrite:
						weight = 1 / (multimapcount)
						output.write(line + "\t" + str(weight) + "\n")
						outpfile.write(line + "\t" + str(weight) + "\n")
					linestowrite = []	
					multimapcount = 0
				multimapcount +=1						
				cigar = curline[5]
				cigarlist = re.findall(r'\d+[A-Z]', cigar)
				for el in cigarlist:
					if el[-1] == "N":
						Nlist.append(int(el[:-1]))	
				Mindexes = [index for index, element in enumerate(cigarlist) if element.endswith("M")]
				count = sum(match > 5 for match in Nlist)	
				if count == 0:
					failed += 1
					failedif = True
				elif count == 1:
					fakegap1 +=1
				elif count > 1:
					truegapm +=1

				if failedif == False:
					curpos = int(curline[3])
					skip = False
					addtomatches = 0 
					poslist = []
					mind = Mindexes[0]
					nextM = 0
					while mind <= Mindexes[-1]:
						if (cigarlist[mind])[-1] == "M":
							if skip == False:
								poslist.append(curpos)
								curpos = curpos + int((cigarlist[mind])[:-1])
								newcigarlist.append(cigarlist[mind])
								nextM += 1 
							if skip == True:
								newcigarlist[-1] = str(int((newcigarlist[-1])[:-1]) + addtomatches + int((cigarlist[mind])[:-1])) + "M"
								curpos = curpos + int((cigarlist[mind])[:-1])
								skip = False
								addtomatches += 0
						elif (cigarlist[mind])[-1] == "N":
							curpos = curpos + int((cigarlist[mind])[:-1])
							if int((cigarlist[mind])[:-1]) <= 5:
								addtomatches += int((cigarlist[mind])[:-1])
								skip = True
						elif (cigarlist[mind])[-1] == "I":								
							curpos = curpos - int((cigarlist[mind])[:-1])
							skip = True
						elif (cigarlist[mind])[-1] == "D":			
							curpos = curpos + int((cigarlist[mind])[:-1])
							addtomatches += int((cigarlist[mind])[:-1])
							skip = True							
						mind +=1
					if curline[-2].startswith("jI:B:i"):
						tag = curline[-2]
					else:
						tag = "jI:B:i,NA"	

					if curline[1] in ["0", "256", "2048"]:
						strand = "+"
					else:
						strand = "-"		
					linecount = 0
					while linecount < len(newcigarlist):
						linestowrite.append(f"{curline[0]}_{str(multimapcount)}\t{curline[1]}\t{curline[2]}\t{poslist[linecount]}\t{curline[4]}\t{newcigarlist[linecount]}\t" + ('\t').join(curline[5:18]) + f"\t{tag}\tSA:Z:{curline[2]},{curline[3]},{strand},{curline[5]},{curline[4]},{(curline[15])[-1]};\tCO:Z:gapm\tCO:Y:{str(len(newcigarlist))}")
						linecount+=1
					cigarlist = []
					newcigarlist = []
					poslist = []	

	print(F"failed {str(failed)} true gapm {str(truegapm)} gap1 {str(fakegap1)}")	


#The function processes gap1 file to convert it into rri format.
#As the result of the function each match part of the alignment is written in the separate line. The original cigar string is written to SA:Z tag in case of multimappers (alignment with indexes _X).
#In case of chimeric alignments (non-multimappers) cigar string is processed as in rri format
#The output is written into outputfilesam
#Weight is equal to multimappers number
def gap1torrinew(gap1file):
	mapqdict = {"2048":"0", "0":"2048", "16":"2064", "2064":"16"}
	mapqdictmapqstr = {"256":"0", "0":"0", "272":"16", "16":"16"}	
	listoflists = []
	pos2 = "dummy" #auxillary variable
	curname = "dummy" #auxillary variable
	scores = []
	readcur = ""
	listsmall = []
	cigarpair1 = ""
	cigarpair2 = ""
	readindex = 0
	outpfile = open(outputfilesam, "a")
	linecount = 0
	with open(gap1file) as gap1:
		for line in gap1:
			if line.startswith("@"): #keeping the header
				pass
			else:
				linecount = linecount + 1
				curline = line.split("\t")		
				if curline[4] != '255': #mapq filter
						if curline[0] == readcur:
							listoflists.append(curline)
							scores.append(int(((curline[13]).split(":"))[-1]))
						elif curline[0] != readcur:
							if readcur != "":		
								maxscore = max(scores)	
								maxscores = ([index for index, value in enumerate(scores) if value == maxscore])
								for el in maxscores:
									listsmall.append(listoflists[el])
								weight = str(1 / len(listsmall))	
								for listik in listsmall:
									curname = listik[0]
									cigar = listik[5]
									cigarlist = re.findall(r'\d+[A-Z]', cigar) #string to list
									pos1 = listik[3]
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
									if listik[1] == "0":
										sign  = "+"
									elif listik[1] == "16":
										sign = "-"
									pos2 = int(pos1)
									
									for el in cigarforpos2:
										if el[:-1] != "I":
											pos2 = pos2 + int(el[:-1])
									listik[5] = cigarpair1
									listik[0] = listik[0] + "_" + str(readindex)
									readindex = readindex + 1
									listik[1] = mapqdictmapqstr[(listik[1])]
									outpfile.write(('\t'.join(map(str, listik)))[:-1] + "\tch:A:1\tSA:Z:" + listik[2] + "," + str(pos2) + "," + sign + "," + cigar + "," + (listik[15])[-1] + ";\tCO:Z:gap1\tCO:Y:2\t" + weight + "\n")
									listik[5] = cigarpair2
									listik[3] = pos2
									listik[1] = mapqdict[(listik[1])]
									outpfile.write(('\t'.join(map(str, listik)))[:-1] + "\tch:A:1\tSA:Z:" + listik[2] + "," + str(pos1) + "," + sign + "," + cigar + "," + (listik[15])[-1] + ";\tCO:Z:gap1\tCO:Y:2\t" + weight + "\n")
							listoflists = []
							listsmall = []
							listoflists.append(curline)
							readcur = curline[0]
							readindex = 0
							scores = []
							scores.append(int(((curline[13]).split(":"))[-1]))
				if curline[4] == '255': #mapq filter
					weight = str(1)
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
							outpfile.write(('\t'.join(map(str, curline)))[:-1] + "\tch:A:1\tSA:Z:" + curline[2] + "," + str(pos2) + "," + sign + "," + cigar + "," + (curline[15])[-1] + ";\tCO:Z:gap1\tCO:Y:2\t" + weight + "\n")
							curline[3] = pos2
							curline[1] = mapqdict[(curline[1])]
							outpfile.write(('\t'.join(map(str, curline)))[:-1] + "\tch:A:1\tSA:Z:" + curline[2] + "," + str(pos1) + "," + sign + "," + cigar + "," + (curline[15])[-1] + ";\tCO:Z:gap1\tCO:Y:2\t" + weight + "\n")
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
						weight = str(1)
						outpfile.write(('\t'.join(map(str, curline)))[:-1] + f"\tCO:Z:gap1\tCO:Y:2\t{weight}\n")
						satag = "SA:Z" + curline[2] + "," + curline[3] + "," + sign + "," + cigar + "," + curline[4] + "," + ((curline[15]).split(":"))[-1] + ";"
						cigarlistchim2 = re.findall(r'\d+[A-Z]', cigar2)
						matcheslist = [index for index, element in enumerate(cigarlistchim2) if element[-1] == 'M']
						cigarpairchim = ""
						cont = 0
						if len(matcheslist) == 1:
							cigarpairchim = cigarlistchim2[matcheslist[0]]
						else:
							mind = matcheslist[0]
							while mind < matcheslist[-1] + 1: #this allows to take in consideration I/D in cigar
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
						outpfile.write(('\t'.join(map(str, curline))) + "\tCO:Z:gap1\tCO:Y:2\t" + weight + "\n")
	print("gap1 done")						

#The function processes homo file to convert it into rri format. 
#The line remains intact but the homo tag and weight are added
#The output is written into outputfilesam
def homotorri(homofile):
	weight = str(1)
	outpfile = open(outputfilesam, "a")
	with open(homofile) as homo:
		for line in homo:
			if line.startswith("@"): #keeping the header
				pass
			else:
				outpfile.write(line[:-1] + "\tCO:Z:homo\tCO:Y:2\t" + weight + "\n")
	print("homo done")				

#The function processes rri file. 
#The line remains intact but the rri tag and weight are added
#The output is written into outputfilesam
def rriaddtag(rrifile):	
	weight = str(1)	
	outpfile = open(outputfilesam, "a")
	with open(rrifile) as rri:
		for line in rri:
			if line.startswith("@"): #keeping the header
				pass
			else:
				outpfile.write(line[:-1] + "\tCO:Z:rri\tCO:Y:2\t" + weight + "\n")
	print("rri done")				

#The functions processes sam file with all the alignments in rri format and converts it into bed file keeping all the necessary information
def createbed(inputsamfile, outputbedfile):
	outpfile = open(outputbedfile, "w")
	outpfile.write("chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tcigar\torigin\tpairsnumber\tweight\n")
	with open(inputsamfile, "r") as samfile:
		for line in samfile:
			alignment = line.split("\t")
			chrom = alignment[2]
			chromStart = alignment[3]
			cigarlist = re.findall(r'\d+[A-Z]', alignment[5])
			cont = sum([int(element) for element in ([element[:-1] for element in cigarlist])])

			if alignment[1] == "0" or alignment[1] == "2048" or alignment[1] == "256":
				strand = "+"
			elif alignment[1] == "16" or alignment[1] == "2064" or alignment[1] == "272":
				strand = "-"
			chromEnd = int(chromStart) + cont	
			name = alignment[0]	
			cigar = (alignment[-4].split(","))[3]
			origin = (alignment[-3].split(":"))[2]
			pairsnumber = ((alignment[-2]).split(":"))[-1] #Is needed for prioritisation process. Number of matches per read (different from 2 just in case of the gapm)
			weight = alignment[-1]
			mapq =  alignment[4]
			outpfile.write(f"{chrom}\t{chromStart}\t{chromEnd}\t{name}\t{mapq}\t{strand}\t{cigar}\t{origin}\t{pairsnumber}\t{weight}")

#The function processed annotation gtf file, extracts exonic features, keeping some gene and transcript information and converts it into bed format. It keeps necessary information for annotation.
#Fantom novel genes are assigned with genetype according to ensembl genetypes (see typedict variable)
def gtftobed(inputgtffile, outputbedfile):
	typedict = {
		"coding_mRNA" : "protein_coding",
		"lncRNA_antisense": "lncRNA",
		"lncRNA_divergent": "lncRNA",
		"lncRNA_intergenic": "lncRNA",
		"lncRNA_sense_intronic": "lncRNA",
		"sense_overlap_RNA": "lncRNA",
		"short_ncRNA": "ncRNA",
		"uncertain_coding": "others",
	}
	fantomgenes = {}
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
				if novelty[-2] == "gene_class":
					fantomgenes[gene_id] = typedict[(novelty[-1])]
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
				gene_id = featuresvocabulary['gene_id']
				source = feature[1]
				try:
					transcript_type = fantomgenes[gene_id]
				except KeyError:
					try:
						transcript_type = featuresvocabulary['transcript_type']
					except KeyError:
						transcript_type = "NA"
				try:
					gene_name = featuresvocabulary['gene_name']
				except KeyError:
					gene_name = "NA"
				try:
					exon_number = featuresvocabulary['exon_number']
				except KeyError:
					exon_number = "NA"	
				try:
					transcript_novelty = transdictionary[transcript_id]
				except KeyError:
					transcript_novelty = "GENCODEv39" #only assigned manually to tRNA 
				try:
					gene_novelty = genedictionary[gene_id]
				except KeyError:
					gene_novelty = "GENCODEv39" #only assigned manually to tRNA
				outpfile.write(f"{chrom}\t{chromStart}\t{chromEnd}\t.\t.\t{strand}\t{transcript_id}\t{transcript_type}\t{gene_name}\t{gene_id}\t{exon_number}\t{transcript_novelty}\t{gene_novelty}\t{source}\n")
				featuresvocabulary = {}

#The function intersects gaptypes total bed file with the annotation bed file
def bedtoolsintersect(annofile, intfile, outfile):
	command = "bedtools intersect -a " + intfile + " -b " + annofile + " -wo -s > " + outfile
	result = subprocess.run(command, shell=True)			

#The function is a prioritisation scheme used in a transcript selection (further function).
#The prioritisation scheme is needed to ovrcome transcript redundancy problem. From the list of the transcripts assigned to a read it selects the most "suitable one", comparing pairs of transcripts till only one possible transcript remains.
def prioritisation(dummylist, count, typeslist, trnoveltylist, tpm, tpmlist):
	while count > 1:	
		if (dummylist[-1])[-9] == (dummylist[-2])[-9]: #if assigned to the same transcript
			dummylist.pop(-1)
		elif (dummylist[-1])[17] != (dummylist[-2])[17]: #biotype selection
			if typeslist[((dummylist[-1])[17])] > typeslist[((dummylist[-2])[17])]: 
				dummylist.pop(-1)
			else:
				dummylist.pop(-2)		
		elif tpm == True and tpmlist[((dummylist[-1])[-9])] != tpmlist[((dummylist[-2])[-9])]: #tpm selection
			if tpmlist[((dummylist[-1])[-9])] < tpmlist[((dummylist[-2])[-9])]:
				dummylist.pop(-1)
			else:
				dummylist.pop(-2)			
		elif (dummylist[-1])[21] != (dummylist[-2])[21]: #novelty selection (less new is selected)
			if trnoveltylist[((dummylist[-1])[21])] > trnoveltylist[((dummylist[-2])[21])]:
				dummylist.pop(-1)
			else:
				dummylist.pop(-2)
		elif (dummylist[-1])[-1] != (dummylist[-2])[-1]: #overlap selection
			if (dummylist[-1])[-1] > (dummylist[-2])[-1]:
				dummylist.pop(-2)
			else:
				dummylist.pop(-1)
		elif (dummylist[-1])[-9] != (dummylist[-2])[-9]: #alphabetic order
			if (dummylist[-1])[-9] > (dummylist[-2])[-9]:
				dummylist.pop(-2)
			else:
				dummylist.pop(-1)
		else: #random selection
			dummylist.pop(-1)	
		count = count - 1
	return dummylist[0]	

#The function collects all the possible annotation features for each read and through the proiritisation scheme selects one transcript for one match in cigar (in the interaction).
def selecttranscripts(inputbed, outputbed, statsfile, tpm, countsint, statsintd):
	output = open(outputbed, "w")
	stats = open(statsfile, "w")
	countsinter = open(countsint, "w")
	statsinter = open(statsintd, "w")
	read = ""
	chr = ""
	start = ""
	strand = ""
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
	gapmlists = []
	number = 0
	gapments = []
	tpmlist = {}
	if tpm == True:
		with open(tpmfile) as tpms:
			tpms.readline()
			for line in tpms:
				linetpm = line.split("\t")
				tpmlist[(linetpm[0])] = float(linetpm[1])
		print("tpm created")
	with open(inputbed, "r") as file0:
		for line in file0:
			ents = line.split("\t")
			if ents[7] == "gapm" or ents[7] == "gap1":			
				if ents[8] == "2":
					ents[7] = "gap1"
				if ents[3] == read	and ents[1] == start:
					dummylist.append(ents)
				elif ents[3] == read:
					start = ents[1]
					gapmlists.append(dummylist)
					dummylist = []
					dummylist.append(ents)
				else:
					gapmlists.append(dummylist)
					read = ents[3]
					if len(gapmlists) == int(number) and number != 0:
						for listm in gapmlists:
							gapments.append(prioritisation(listm, len(listm), typeslist, trnoveltylist, tpm, tpmlist))
						if tpm == True:
							for ent in gapments:								
								output.write(('\t'.join(ent))[:-1] + "\t" + str(tpmlist[(ent[16])]) + "\n")
						else:
							for ent in gapments:
								output.write(('\t'.join(ent))[:-1] + "\n")
						part1 = ('\t').join([ent[16] for ent in gapments]) + "\t"
						part2 = ('+').join(sorted([ent[16] for ent in gapments])) + "\t"
						part3 = ('\t').join([ent[17] for ent in gapments]) + "\t"
						part4 = ('+').join(sorted([ent[17] for ent in gapments])) + "\t"
						countscounts = len(gapments)
						while countscounts < 4:
							part1 += ".\t"
							part3 += ".\t"
							countscounts += 1
						stats.write(f"{part1}{part2}{part3}{part4}{(gapments[0])[9]}\t{ents[7]}\n")	
						countsinter.write((gapments[0])[9] + "\t" + part2 + "\n")
						statsinter.write((gapments[0])[9] + "\t" + part4 + "\n")
						gapments = []
					start = ents[1]		
					number = ents[8]
					dummylist = []
					dummylist.append(ents)
					gapmlists = []
			if ents[7] == "homo":
				if ents[3] == read:
					dummylist.append(ents)
				else:
					read = ents[3]
					ent1 = prioritisation(dummylist, len(dummylist), typeslist, trnoveltylist, tpm, tpmlist)
					if tpm == True:
						output.write(('\t'.join(ent1))[:-1] + "\t" + str(tpmlist[(ent1[16])]) + "\n")
					else:
						output.write(('\t'.join(ent1)))
					stats.write(f"{ent1[16]}\t{ent1[16]}\t.\t.\t{ent1[16]}+{ent1[16]}\t{ent1[17]}\t{ent1[17]}\t.\t.\t{ent1[17]}+{ent1[17]}\t1\thomo\n")
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
							ent1 = prioritisation(dummylist1, len(dummylist1), typeslist, trnoveltylist, tpm, tpmlist)
							ent2 = prioritisation(dummylist, len(dummylist), typeslist, trnoveltylist, tpm, tpmlist)
							if tpm == True:
								output.write(('\t'.join(ent1))[:-1] + "\t" + str(tpmlist[(ent1[16])]) + "\n")
								output.write(('\t'.join(ent2))[:-1] + "\t" + str(tpmlist[(ent1[16])]) + "\n")
							else:
								output.write(('\t'.join(ent1)))
								output.write(('\t'.join(ent2)))
							stats.write(f"{ent1[16]}\t{ent2[16]}\t.\t.\t{max((ent1[16]), (ent2[16]))}+{min((ent1[16]), (ent2[16]))}\t{ent1[17]}\t{ent2[17]}\t.\t.\t{max((ent1[17]), (ent2[17]))}+{min((ent1[17]), (ent2[17]))}\t1\trri\n")
						readcount = False
					dummylist = []
					chr = ents[0]
					strand = ents[5]
					dummylist.append(ents)	

gtftobed(annogtffile, annobedfile)

createsam(outputfilesam)
gapmtorri(gapmfile)
gap1torrinew(gap1file)	
homotorri(homofile)
rriaddtag(rrifile)
createbed(outputfilesam, outputfilebed)
bedtoolsintersect(annobedfile, outputfilebed, intersectedbed)
selecttranscripts(intersectedbed, selectedintersected, stats, tpmstatus, counts3, stats2)
