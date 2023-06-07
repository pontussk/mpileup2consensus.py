import sys

from optparse import OptionParser

usage = "usage: %prog [options] <SAM formatted data from stdin>"
parser = OptionParser(usage=usage)
parser.add_option("--name", action="store", type="string", dest="name",help="Name of the fasta sequence (individual)",default="X")
parser.add_option("--minsupport", action="store", type="float", dest="minsupport",help="Minimum consensus support",default=0.8)
parser.add_option("--mindepth", action="store", type="int", dest="mindepth",help="Minimum depth",default=1)
parser.add_option("--maxdepth", action="store", type="int", dest="maxdepth",help="Maximum depth",default=10000000000)
parser.add_option("--noheader", action="store_true", dest="noheader",help="Do not print header describing the columns in the output",default=False)
parser.add_option("--ssDNAlib", action="store_true", dest="ssDNAlib",help="genotype ssDNA libraries by dropping all Ts that aligned to the forward strand and all As that aligned to the reverse strand",default=False)
parser.add_option("--ssDNAlib_refcall", action="store_true", dest="ssDNAlib_refcall",help="genotype ssDNA libraries by dropping all Ts (when the ref is C) that aligned to the forward strand and all As (when the ref is G) that aligned to the reverse strand",default=False)

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

def ssDNAbasefun(ref,inp):
	string =''
	for thebase in inp:
		if   thebase=='T':continue
		elif thebase=='a':continue
		elif thebase=='.' and ref=='T': continue
		elif thebase==',' and ref=='A': continue
		elif thebase.isalpha() and thebase.upper() in nucleotides:
			string += thebase.upper()
		else:
			if thebase in ['.',',']:
				string += ref
	return string
	
def ssDNArefbasefun(ref,inp):
	string =''
	for thebase in inp:
		if   thebase=='T' and ref.upper()=='C':continue
		elif thebase=='a' and ref.upper()=='G':continue
		elif thebase.isalpha() and thebase.upper() in nucleotides:
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

	if depth > options.maxdepth:
		newallele='N'
		newseq += newallele
		continue

	pileupcol=col[4]
	
	if options.ssDNAlib==True:
		pileup=ssDNAbasefun(refbase,pileupcol)
		if len(pileup) <1:
			newallele='N'
			newseq += newallele
			continue
	elif options.ssDNAlib_refcall==True:
		pileup=ssDNArefbasefun(refbase,pileupcol)
		if len(pileup) <1:
			newallele='N'
			newseq += newallele
			continue
	else:
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
