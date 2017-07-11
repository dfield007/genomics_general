#!/usr/bin/env python
# Essential modules
import re, sys, getopt, itertools
print "                                  "
print "----------------------------------"
print "-  Sequence extraction for SNPs  -"
print "-     For 8 pools of WGS         -"
print "-      David Field 23/10/2013    -"
print "----------------------------------"
print "                                  "
print " * to run script, you need (i) script name, (ii) a list of snps to extract, (iii) a sync file, (iv) an output file name,"
print "     (v) the min read depth for sequence, (vi) min copies for an allele call (i.e. > 1, no singletons)"
print " * e.g, SNPextract8poolNew.py SNPlistPilot2test.txt sync.txt snpoutputtest.txt 15 2"
print " * e.g, SNPextract8poolFast.py cpSNPlistNewtest.txt yp4_yelPool20_yp2_yp1_mp2_mp4_mgtPool20_mp11.cp.stampy.markdup.sync snpoutput_cp_test.txt 15 100"
print " * Note: only works if SNP list file has line breaks saved as Unix (LF). "
print " "
print "Importing data... "
print 'Snp list (Input file 1):', sys.argv[1]
print 'Sync file (Input file 2):', sys.argv[2]
print 'Output file:', sys.argv[3]
print 'Min depth:', sys.argv[4]
print 'Min allele call:', sys.argv[5]

# read in the whole SNP list
InFile = open(sys.argv[1], 'r')
SNPlist  = [i.split()[0:] for i in InFile.readlines()]
InFile.close()
SNPlist = SNPlist[1:len(SNPlist)]
numLines = sum(1 for line in SNPlist)
numLines = float(int(numLines))
# read in just the scaffold names
InFile = open(sys.argv[1], 'r')
scaffs_looking_for  = [i.split()[0:1] for i in InFile.readlines()]
scaffs_looking_for = scaffs_looking_for[1:len(scaffs_looking_for)]
InFile.close() 
scaffs_looking_for = list(itertools.chain.from_iterable(scaffs_looking_for))
scaffs_looking_for = set(scaffs_looking_for)

#print "Matrix of SNP list:", SNPlist
#print "Scaffold Names:", scaffs_looking_for
# Open and import the SNP list file
#InFileTest = open(sys.argv[1], 'r')
#print InFileTest
#numLines = sum(1 for line in InFileTest)
#InFileTest  = [i.split()[0:] for i in InFile.readlines()]
# lengthInFile=len(InFileTest)
# print lengthInFile
# numLines = (sum(1 for line in InFile))-1
#InFile.close()
print 'Number of SNPs to process:', (numLines)
print " "
print "Processing... "
print " "

OutFile = open(sys.argv[3], 'w')
# Initialise 
LineNumber = 0
thisScaff = []
thisPos = []
theDataSummary = []
theDataFull = []
minDepthSeqInfo=int(sys.argv[4])
minAlleleCount=int(sys.argv[5])
# Synch file setup
startFlag = 0
stopcounter = 0
switchOn = 0
numNonBiallele=0
SNPProcessCounter=0
theDataTemp=[]
theDataFull=[]
# open big sync file
with open(sys.argv[2]) as InFile2:
	for line in InFile2:
		if LineNumber == 0:
			OutputString = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % \
				("scaffold", "position","locus","refallele","alleles","seqIUPAC","Depth","leftSeq","rightSeq","minSeqEdge","seqAvail","lenIUPAC","numNonBiallele")
			# Outputting to file
			print OutputString
			OutFile.write(OutputString+"\n")
			theDataFull.append(OutputString+"\n")
			LineNumber=LineNumber+1
		if line.split(None, 1)[0] in scaffs_looking_for:
			thisLine = line.split(None, 2)[0:2]
			#print thisLine
			thisScaff = thisLine[0]
			currentPos = int(thisLine[1])
			#print "Scaff & Pos:", thisScaff, currentPos
			res=[n for n, (s, i, p) in enumerate(SNPlist) if s==thisScaff and currentPos>=(int(i)-50) and currentPos<=(int(i)+50)]
			if res:
				#print currentPos, res
				switchOn = 1
				LineList2 = line.split()
				currentScaff = LineList2[0]
				currentPos = LineList2[1]
				currentScaff = str(currentScaff)
				currentPos = int(currentPos)
				res = res[0]
				thisSnpPos = SNPlist[res][1]
				MarkerType = SNPlist[res][2]
				print "This SNP position:", thisSnpPos
				print currentScaff, currentPos, res
				if (startFlag==0):
					startFlag = 1
				stopcounter = stopcounter + 1
				refAllele = str(LineList2[2])
				print "refAllele:", refAllele
				seq = "NA"
				minDepth = "NA"
				focal = "N"
				if (int(currentPos)==int(thisSnpPos)):
					focal = "Y"
					refAlleleFocalSNP=refAllele
				# pool1
				alleles1=str(LineList2[3])
				alleles1=alleles1.split(":")
				alleles1=alleles1[0:4]
				# pool2
				alleles2=str(LineList2[4])
				alleles2=alleles2.split(":")
				alleles2=alleles2[0:4]
				# pool3
				alleles3=str(LineList2[5])
				alleles3=alleles3.split(":")
				alleles3=alleles3[0:4]
				# pool4
				alleles4=str(LineList2[6])
				alleles4=alleles4.split(":")
				alleles4=alleles4[0:4]
				# pool5
				alleles5=str(LineList2[7])
				alleles5=alleles5.split(":")
				alleles5=alleles5[0:4]
				# pool6
				alleles6=str(LineList2[8])
				alleles6=alleles6.split(":")
				alleles6=alleles6[0:4]
				# pool7
				alleles7=str(LineList2[9])
				alleles7=alleles7.split(":")
				alleles7=alleles7[0:4]
				# pool8
				alleles8=str(LineList2[10])
				alleles8=alleles8.split(":")
				alleles8=alleles8[0:4]
				# seq code transform
				seqCode=["A","T","C","G"]
				print "Alleles in each pool:", alleles1, alleles2, alleles3, alleles4, alleles5, alleles6, alleles7, alleles8
				# pool 1
				counter=0
				for val in alleles1:
					alleles1[counter]=int(val)
					counter=counter+1
				counter=0
				alleles1present = [1 if x>=minAlleleCount else 0 for x in alleles1]
				for allele in alleles1present:
					alleles1present[counter]=int(allele)*seqCode[counter]
					counter=counter+1
				counter=0
				# pool 2
				for val in alleles2:
					alleles2[counter]=int(val)
					counter=counter+1
				alleles2present = [1 if x>=minAlleleCount else 0 for x in alleles2]
				counter=0
				for allele in alleles2present:
					alleles2present[counter]=int(allele)*seqCode[counter]
					counter=counter+1
				counter=0	
				# pool 3
				for val in alleles3:
					alleles3[counter]=int(val)
					counter=counter+1
				alleles3present = [1 if x>=minAlleleCount else 0 for x in alleles3]
				counter=0
				for allele in alleles3present:
					alleles3present[counter]=int(allele)*seqCode[counter]
					counter=counter+1
				counter=0	
				# pool 4
				for val in alleles4:
					alleles4[counter]=int(val)
					counter=counter+1
				alleles4present = [1 if x>=minAlleleCount else 0 for x in alleles4]
				counter=0
				for allele in alleles4present:
					alleles4present[counter]=int(allele)*seqCode[counter]
					counter=counter+1
				counter=0				
				# pool 5
				for val in alleles5:
					alleles5[counter]=int(val)
					counter=counter+1
				alleles5present = [1 if x>=minAlleleCount else 0 for x in alleles5]
				counter=0
				for allele in alleles5present:
					alleles5present[counter]=int(allele)*seqCode[counter]
					counter=counter+1
				counter=0					
				# pool 6
				for val in alleles6:
					alleles6[counter]=int(val)
					counter=counter+1
				alleles6present = [1 if x>=minAlleleCount else 0 for x in alleles6]
				counter=0
				for allele in alleles6present:
					alleles6present[counter]=int(allele)*seqCode[counter]
					counter=counter+1
				counter=0	
				# pool 7
				for val in alleles7:
					alleles7[counter]=int(val)
					counter=counter+1
				alleles7present = [1 if x>=minAlleleCount else 0 for x in alleles7]
				counter=0
				for allele in alleles7present:
					alleles7present[counter]=int(allele)*seqCode[counter]
					counter=counter+1
				counter=0	
				# pool 8
				for val in alleles8:
					alleles8[counter]=int(val)
					counter=counter+1
				alleles8present = [1 if x>=minAlleleCount else 0 for x in alleles8]
				counter=0
				for allele in alleles8present:
					alleles8present[counter]=int(allele)*seqCode[counter]
					counter=counter+1
				counter=0	
				print "Alleles <min removed: ", alleles1present, alleles2present, alleles3present, alleles4present, alleles5present, alleles6present, alleles7present, alleles8present
				# depths
				depthPop1=sum(alleles1)
				depthPop2=sum(alleles2)
				depthPop3=sum(alleles3)
				depthPop4=sum(alleles4)
				depthPop5=sum(alleles5)
				depthPop6=sum(alleles6)
				depthPop7=sum(alleles7)
				depthPop8=sum(alleles8)
				if (MarkerType=="D"):
					minDepth=min([depthPop1,depthPop8], key=int)
				if (MarkerType=="P"):
					minDepth=min([depthPop3,depthPop4,depthPop5,depthPop6], key=int)
				# IUPAC conversion
				seqIUPAC="NA"
				allelesBothPops = alleles1present + alleles2present + alleles3present + alleles4present + alleles5present + alleles6present + alleles7present + alleles8present
				allelesBothPopsList = set(allelesBothPops)
				allelesBothPopsList = list(allelesBothPopsList)
				allelesBothPopsList.sort()
				allelesBothPopsList = ''.join(allelesBothPopsList) 
				#print "allelesBothPopsList:",allelesBothPopsList
				OutputList = (thisScaff,currentPos,refAllele,allelesBothPopsList,minDepth)
				print "OutputList:",OutputList
				seqIUPAC = ''.join(allelesBothPopsList) 
				allelesBothPopsListJoin=''.join(allelesBothPopsList) 
				if (len(seqIUPAC)>2):
					numNonBiallele=numNonBiallele+1
				if (int(currentPos)==int(thisSnpPos)):
					OutputListFocal = OutputList
				if (seqIUPAC=='AG'):
					seqIUPAC="R"
				if (seqIUPAC=='CT'):
					seqIUPAC="Y"
				if (seqIUPAC=='CG'):
					seqIUPAC="S"
				if (seqIUPAC=='AT'):
					seqIUPAC="W"
				if (seqIUPAC=='GT'):
					seqIUPAC="K"
				if (seqIUPAC=='AC'):
					seqIUPAC="M"
				print "Current & thiPos:", currentPos, thisSnpPos
				print "(currentPos==thisSnpPos)?:", (currentPos==thisSnpPos)
				if (int(currentPos)==int(thisSnpPos)):
					seqIUPAC = "[",seqIUPAC,"]"
					seqIUPAC=''.join(seqIUPAC)
				# empty readings use ref allele instead
				print "seqIUPAC: ", seqIUPAC
				if (allelesBothPopsListJoin.__len__() is 0):
					 allelesBothPopsListJoin=refAllele
					 seqIUPAC=refAllele
				# Save to current dataset
				theDataTemp.append([thisScaff,currentPos,refAllele,allelesBothPopsListJoin,seqIUPAC,minDepth,focal])
				theDataFull.append([thisScaff,currentPos,refAllele,allelesBothPopsListJoin,seqIUPAC,minDepth,focal])
				# 100 bp sequence details
				if (stopcounter == 101):
					# turn off switch
					switchOn = 0
					seqDepthCompress = []
					seqIUPACcompress = []
					rangeVals=range(0,101)
					for val in rangeVals:
						thisOne=theDataTemp[val][4]
						thisOne2=str(theDataTemp[val][5])
						seqIUPACcompress.append(thisOne)
						seqDepthCompress.append(thisOne2)
					# How close to depth dropoff (i.e. end of short read)
					leftSeq = 0
					rightSeq = 0
					counter = 0
					# Extent of good Seq depth around SNP
					for val in seqDepthCompress:
						if (int(val) >= minDepthSeqInfo):
							if (counter < 50):
								leftSeq = leftSeq+1
							if (counter > 50):
								rightSeq = rightSeq+1
						counter = counter+1
					minSeqAvail=min([leftSeq,rightSeq], key=int)
					seqIUPACcompressNew=''.join(seqIUPACcompress)
					seqDepthCompressNew=','.join(seqDepthCompress)
					edge=len(thisScaff)
					numCodeThisScaff=thisScaff[8:edge]
					locus=["s",str(numCodeThisScaff),"_",str(thisSnpPos)]
					locus=''.join(locus)
					lengthIUPAC = (len(seqIUPACcompressNew)-2)
					theDataSummary.append([thisScaff,thisSnpPos,locus,OutputListFocal[2],OutputListFocal[3],\
					seqIUPACcompressNew,seqDepthCompressNew,minDepth,leftSeq,rightSeq,minSeqAvail,100])
					OutputString3 = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (thisScaff,thisSnpPos,locus,OutputListFocal[2],\
						OutputListFocal[3],seqIUPACcompressNew,seqDepthCompressNew,leftSeq,rightSeq,minSeqAvail,100,lengthIUPAC,numNonBiallele)
					print OutputString3 
					OutFile.write(OutputString3+"\n")
					SNPProcessCounter=SNPProcessCounter+1
					LineNumber=LineNumber+1
					# clear out
					dataOutputFile=1
					theDataTemp=[]
					stopcounter = 0
					startFlag = 0
				# If goes outside range but switch still on, need to collect data for shorter region than 100bp
				if (currentPos>(int(i)+50) and switchOn == 1):
					switchOn = 0
					seqDepthCompress = []
					seqIUPACcompress = []
					maxReads = len(theDataTemp)
					rangeVals=range(0,(maxReads-1))
					for val in rangeVals:
						thisOne=theDataTemp[int(val)][4]
						thisOne2=str(theDataTemp[int(val)][5])
						seqIUPACcompress.append(thisOne)
						seqDepthCompress.append(thisOne2)
					# How close to depth edge
					leftSeq = 0
					rightSeq = 0
					counter = 0
					# Extent of good Seq depth around SNP
					for val in seqDepthCompress:
						if (int(val) >= minDepthSeqInfo):
							if (counter < 50):
								leftSeq = leftSeq+1
							if (counter > 50):
								rightSeq = rightSeq+1
						counter = counter+1
					minSeqAvail=min([leftSeq,rightSeq], key=int)
					seqIUPACcompressNew=''.join(seqIUPACcompress)
					seqDepthCompressNew=','.join(seqDepthCompress)
					edge=len(thisScaff)
					numCodeThisScaff=thisScaff[8:edge]
					locus=["s",str(numCodeThisScaff),"_",str(thisSnpPos)]
					locus=''.join(locus)
					lengthIUPAC = (len(seqIUPACcompressNew)-2)
					theDataSummary.append([thisScaff,thisSnpPos,locus,OutputListFocal[2],OutputListFocal[3],\
						seqIUPACcompressNew,seqDepthCompressNew,leftSeq,rightSeq,minSeqAvail,maxReads])
					OutputString3 = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (thisScaff,thisSnpPos,locus,OutputListFocal[2],\
						OutputListFocal[3],seqIUPACcompressNew,seqDepthCompressNew,leftSeq,rightSeq,minSeqAvail,maxReads,lengthIUPAC,numNonBiallele)
					print OutputString3
					SNPProcessCounter=SNPProcessCounter+1
					LineNumber=LineNumber+1
					print " "
					sys.stdout.write(OutputString3)
					sys.stdout.flush()
					#print OutputString3
					OutFile.write(OutputString3+"\n")
					# clear out
					theDataTemp=[]
					stopcounter = 0
					dataOutputFile=1
				startFlag = 0
		PercentComplete = round((SNPProcessCounter/numLines)*100,2)
		PercentComplete = str(float(PercentComplete))
		sys.stdout.write("\rPercent Complete: %s" % PercentComplete)
		sys.stdout.flush()
#OutFile.close()
print " "
print " "
print "Processing complete"