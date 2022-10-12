import sys

from optparse import OptionParser

usage = "usage: %prog [options] <SAM formatted data from stdin>"
parser = OptionParser(usage=usage)
parser.add_option("--name", action="store", type="string", dest="name",help="Name of the fasta sequence (individual)",default="X")
parser.add_option("--minsupport", action="store", type="float", dest="minsupport",help="Minimum consensus support",default=0.8)
parser.add_option("--mindepth", action="store", type="int", dest="mindepth",help="Minimum depth",default=1)
parser.add_option("--noheader", action="store_true", dest="noheader",help="Do not print header describing the columns in the output",default=False)
(options, args) = parser.parse_args()




threshold=options.minsupport
header=options.name

nucleotides=['A','C','T','G']
def basefun(ref,inp):
	string =''
	for thebase in inp:
		if thebase.isalpha() and thebase.upper() in nucleotides:
			string += thebase.upper()
		else:
			if thebase in ['.',',']:
				string += ref
	return string

newseq=''
for line in sys.stdin:
	col=line.split()
	refbase=col[2]
	depth=int(col[3])
	if depth < options.mindepth:
		newallele='N'
		newseq += newallele
		continue

	pileupcol=col[4]
	pileup=basefun(refbase,pileupcol)	
		
	alleles=list(set(pileup))
	allalleles=list(set(pileup+refbase))
	#find major base
	alleledict={}
	refcount=0
	for allele in alleles:
		alleledict[allele]=pileup.count(allele)

	counts=alleledict.values()
	if len(counts) <1:
		newallele='N'
		newseq += newallele
		continue
	maxcount=max(counts)

	chosenallele=[]	
	for allele in alleledict.keys():
		if alleledict[allele] == maxcount:
			chosenallele.append(allele)
	if len(chosenallele) > 1:
		newallele='N'
		newseq += newallele
		continue
	majorallele=chosenallele[0]
	majF=1.0*maxcount/depth
	
	if majF >= threshold:
		newallele=majorallele
		#newseq +=majorallele
	else:
		#newseq += 'N'
		newallele='N'
	newseq += newallele

	#print newallele,line, 
		
print '>'+header		
print newseq