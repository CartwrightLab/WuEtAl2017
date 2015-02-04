#! /usr/local/bin/python
import sys
import gzip
import os
import errno
import time

class info:
    def __init__(self, callby, ref, alt, variable, altalt, exome):
        self.callby = callby        #called by (1)single caller, (3) trio caller, (2) both
        self.ref = ref              #reference allele
        self.alt = alt              #alternate allele (same as ref for homo ref)
        self.variable = variable    #(1) found as variable in the 1000 genomes data - what we're using as our standard for true snps
        self.altalt = altalt        #(1) a different alterate allele was found in the 1000 genomes data
        self.exome = exome          #(1) in WEx data - not sure if this is the whole exome, but within these data we can compare N and S snp freq
        
#################################
chr=sys.argv[1]
person=sys.argv[2]
variable_site_file=sys.argv[3] # should be able to take both gz and unzip
working_dir = sys.argv[4]
pileupFile = sys.argv[5]
exome=str(sys.argv[6])

allcallsdict=dict()     #dict contains dicts of (1)hets, (2)Hr, (3)Ha
allcallsdict['hets']=dict()         #each dict contains calls as class info
allcallsdict['homos_ref']=dict()
allcallsdict['homos_alt']=dict()

#make dictionaries for hr, ha, homo for trio sites    
print "Parse calls from trio"
#filein = open('./calls_trio_vcfs/vcf_file_list.txt', 'r')
filein = open(working_dir+"vcf_file_trio_list.txt", 'r')
filenames=filein.readlines()
for filename in filenames:
    vcf = open(filename.rstrip(), 'r')
    for line in vcf:
        if line.startswith(chr):
            splitline=line.split()
            if((len(splitline[3])>1) or (len(splitline[4])>1)):     #skip indels
                continue
            elif (splitline[4] is '.'): #homo ref
                allcallsdict['homos_ref'][splitline[1]]=info(3, splitline[3], splitline[3], 0, 0, 0)        #record
            else:
                geno = (splitline[9]).split(':') #if kid has alt allele
                if (geno[0][2]) is '1':
                    if (geno[0][0]) is '0':     #het
                        allcallsdict['hets'][splitline[1]]=info(3, splitline[3], splitline[4], 0, 0, 0)        #record
                    if (geno[0][0]) is '1':     #homo alt
                        #allcallsdict['hets'][splitline[1]]=info(3, splitline[3], splitline[4], 0, 0, 0)        #record                  
                        allcallsdict['homos_alt'][splitline[1]]=info(3, splitline[3], splitline[4], 0, 0, 0)        #record                  
    vcf.close()

#get single called sites - add to prev dictionary
print "Parse calls from single"
#filein = open('./calls_single_vcfs/vcf_file_list.txt', 'r')
filein = open(working_dir+"vcf_file_single_list.txt", 'r')
filenames=filein.readlines()
for filename in filenames:
    vcf = open(filename.rstrip(), 'r')
    for line in vcf:
        if line.startswith(chr):
            splitline=line.split()
            if((len(splitline[3])>1) or (len(splitline[4])>1)):     #skip indels
                continue
            elif (splitline[4] is '.'): #homo ref
                if splitline[1] in allcallsdict['homos_ref']:     #already found in trio caller
                    allcallsdict['homos_ref'][splitline[1]].callby = 2
                else:
                    allcallsdict['homos_ref'][splitline[1]]=info(1, splitline[3], splitline[3], 0, 0, 0)
            else:
                geno = (splitline[9]).split(':') #if kid has alt allele
                if (geno[0][2]) is '1':
                    if (geno[0][0]) is '0':     #het
                        if splitline[1] in allcallsdict['hets']:     #already found in trio caller
                            allcallsdict['hets'][splitline[1]].callby = 2
                        else:
                            allcallsdict['hets'][splitline[1]]=info(1, splitline[3], splitline[4], 0, 0, 0)
                    if (geno[0][0]) is '1':     #homo alt
                        if splitline[1] in allcallsdict['homos_alt']:
                            allcallsdict['homos_alt'][splitline[1]].callby = 2
                        else:                            
                            #allcallsdict['hets'][splitline[1]]=info(1, splitline[3], splitline[4], 0, 0, 0)
                            allcallsdict['homos_alt'][splitline[1]]=info(1, splitline[3], splitline[4], 0, 0, 0)

#go through variable sites file and note variable in each sub-dictionary 
## TODO sholud implement a test on file type here, do this later
print "Get variable sites"
ref_data = open(variable_site_file, 'r')
#ref_data = gzip.open(variable_site_file, 'rb')
for line in ref_data:      #read line in file
    if line.startswith(chr):            #skip all initial info lines
        splitline=line.split()
        if(len(splitline[3])>1 or len(splitline[4])>1):     #skip indels
            continue
        else:
            if splitline[1] in allcallsdict['hets']:  #found by us as het             #dictionaries are non-unique because a site found by the trio as one type could be found by single as another
                allcallsdict['hets'][splitline[1]].variable = 1        #note prev found
                if allcallsdict['hets'][splitline[1]].alt is not splitline[4]:        #check alt allele
                    allcallsdict['hets'][splitline[1]].altalt = 1 
            if splitline[1] in allcallsdict['homos_alt']:  #found by us as homoalt
                allcallsdict['homos_alt'][splitline[1]].variable = 1        #note prev found
                if allcallsdict['homos_alt'][splitline[1]].alt is not splitline[4]:        #check alt allele
                    allcallsdict['homos_alt'][splitline[1]].altalt = 1     
ref_data.close()

#go through exome data and note if found in exome
if exome is not '0':
    print "Pares exome file %s" % exome
    #filein = open('chr'+chr+'Ex_'+person+'.pileups', 'r')        #check if in exome    
    filein = open(exome, 'r')        #check if in exome    
    for site in filein:
        splitline = site.split()
        if len(splitline)>4:  
            chr,pos,ref,count,bases,qual = site.split()
            if pos in allcallsdict['hets']:
                if count is not '0':
                    allcallsdict['hets'][pos].exome = 1
            if pos in allcallsdict['homos_alt']:
                if count is not '0':
                    allcallsdict['homos_alt'][pos].exome = 1
            if pos in allcallsdict['homos_ref']:
                if count is not '0':
                    allcallsdict['homos_ref'][pos].exome = 1
    filein.close()

outfilenames=('homos_ref','homos_alt','hets')

result_dir = working_dir+"base_count/"
try:
  os.mkdir(result_dir)
except OSError as e:
    if not (e.errno == errno.EEXIST and os.path.isdir(result_dir)):
        raise

for filename in outfilenames:
    print "Process %s" % filename
    fileout1 = open(result_dir+'base_counts_'+filename+'_'+person+'_flag_filtered.txt','w')
    fileout1.write('pos'+"\t"+'ref'+"\t"+'alt'+"\t"+'As'+"\t"+'Cs'+"\t"+'Gs'+"\t"+'Ts'+"\t"+'callby'+"\t"+'snp'+"\t"+'snpdif'+"\t"+'Ex'+"\n")
    fileout2 = open(result_dir+'base_counts_'+filename+'_'+person+'_byref_flag_filtered.txt','w')
    fileout2.write('pos'+"\t"+'ref'+"\t"+'alt'+"\t"+'refs'+"\t"+'alts'+"\t"+'e1s'+"\t"+'e2s'+"\t"+'callby'+"\t"+'snp'+"\t"+'snpdif'+"\t"+'Ex'+"\n")
    
    #pileups=open('chr'+chr+'_'+person+'.pileups','r')
    pileups=open(pileupFile,'r')
    for line in pileups:
        splitline = line.split()
        if len(splitline)>4:            
            chr,pos,ref,count,bases,qual = line.split()
            if pos in allcallsdict[filename]:
                if ref is not 'N':                    
                    basecounts=dict()
                    basecounts['ref']=bases.count('.')+bases.count(',')
                    basecounts['A']=bases.count('A')+bases.count('a')
                    basecounts['C']=bases.count('C')+bases.count('c')
                    basecounts['G']=bases.count('G')+bases.count('g')
                    basecounts['T']=bases.count('T')+bases.count('t')
                    basecounts[ref]=basecounts['ref']
                
                    fileout1.write(pos+"\t"+ref+"\t"+allcallsdict[filename][pos].alt+"\t"+str(basecounts['A'])+"\t"+str(basecounts['C'])+"\t"+str(basecounts['G'])+"\t"+str(basecounts['T'])+"\t"+str(allcallsdict[filename][pos].callby)+"\t"+str(allcallsdict[filename][pos].variable)+"\t"+str(allcallsdict[filename][pos].altalt)+"\t"+str(allcallsdict[filename][pos].exome)+"\n")
                
                    baselist=list()
                    baselist.append(basecounts[ref])
                    del basecounts[ref]
                    del basecounts['ref']
                    if allcallsdict[filename][pos].alt is not ref:
                        baselist.append(basecounts[allcallsdict[filename][pos].alt])
                        del basecounts[allcallsdict[filename][pos].alt]
                    for k,v in basecounts.iteritems():
                        baselist.append(v)        
                    fileout2.write(pos+"\t"+ref+"\t"+allcallsdict[filename][pos].alt+"\t"+str(baselist[0])+"\t"+str(baselist[1])+"\t"+str(baselist[2])+"\t"+str(baselist[3])+"\t"+str(allcallsdict[filename][pos].callby)+"\t"+str(allcallsdict[filename][pos].variable)+"\t"+str(allcallsdict[filename][pos].altalt)+"\t"+str(allcallsdict[filename][pos].exome)+"\n")
    pileups.close()   
    fileout1.close()
    fileout2.close()
