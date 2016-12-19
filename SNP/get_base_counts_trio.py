#!/usr/bin/env python2
import sys
import gzip
import os
import errno
import time
# from collections import namedtuple
# info = namedtuple('info', ['callby','ref', 'alt', 'variable', 'altalt', 'exome'])
# def summary_info(info):
#     s = '%i\t%i\t%i\t%i' % (
#         info.callby, info.variable, info.altalt, info.exome)

# class info:
#     __slots__ = ('callby','ref', 'alt', 'variable', 'altalt', 'exome')
#     def __init__(self, callby, ref, alt, variable=0, altalt=0, exome=0):
#         self.callby = callby        #called by (1)single caller, (3) trio caller, (2) both
#         self.ref = ref              #reference allele
#         self.alt = alt              #alternate allele (same as ref for homo ref)
#         self.variable = variable    #(1) found as variable in the 1000 genomes data - what we're using as our standard for true snps
#         self.altalt = altalt        #(1) a different alterate allele was found in the 1000 genomes data
#         self.exome = exome          #(1) with exome data


def build_info (callby, ref, alt, variable=0, altalt=0, exome=0):
    return {'callby':callby, 'ref':ref, 'alt':alt, 'variable':variable, 'altalt':altalt, 'exome':exome}

def summary_info(info):
        s = '%i\t%i\t%i\t%i' % (
            info['callby'], info['variable'], info['altalt'], info['exome'])
        return s


#################################
# Python 2.6 doesn't have .major .minor
if sys.version_info[0] == 3:
    raise SystemError("Python 3.x is not supported")

if sys.version_info[1] < 7:
    raise SystemError("Must be using Python 2.7 or greater")


PREFIX_HET="0/1:"
PREFIX_HOMO_ALT="1/1:"

print('===== Parse output files =====')
chr = sys.argv[1]
person = sys.argv[2]
variable_site_file = sys.argv[3]  # should be able to take both gz and unzip
working_dir = sys.argv[4]
pileup_file = sys.argv[5]
exome = str(sys.argv[6])

all_calls = {'hets': dict(), 'homos_ref':dict(), 'homos_alt':dict()}
#outfile_names = ['hets'] #Use this if you only want hets
outfile_names=['homos_ref','homos_alt','hets']
result_dir = working_dir + 'base_count/'

try:
    os.mkdir(result_dir)
except OSError as e:
    if not (e.errno == errno.EEXIST and os.path.isdir(result_dir)):
        raise

#make dictionaries for hr, ha, homo for trio sites
print ('Parse calls from trio')
sys.stdout.flush()
with open(working_dir+'vcf_file_trio_list.txt', 'r') as filein:
    vcf_files = filein.readlines()

for filename in vcf_files:
    with open(filename.rstrip(), 'r') as vcf:
        for line in vcf:
            if line.startswith(chr):
                splitline=line.split()
                if((len(splitline[3])>1) or (len(splitline[4])>1)):     #skip indels
                    continue
                elif (splitline[4] is '.'): #homo ref
                    all_calls['homos_ref'][splitline[1]] = build_info(3, splitline[3], splitline[3])
                else:
                    geno = splitline[9][0:4] #if kid has alt allele
                    if PREFIX_HET == geno: #het
                        all_calls['hets'][splitline[1]]=build_info(3, splitline[3], splitline[4])
                    elif PREFIX_HOMO_ALT == geno : #homo alt
                        all_calls['homos_alt'][splitline[1]]=build_info(3, splitline[3], splitline[4])

#get single called sites - add to prev dictionary
print('Parse calls from single')
sys.stdout.flush()
with open(working_dir+'vcf_file_single_list.txt', 'r') as filein:
    vcf_files = filein.readlines()

for filename in vcf_files:
    with open(filename.rstrip(), 'r') as vcf:
        for line in vcf:
            if line.startswith(chr):
                splitline=line.split()
                if((len(splitline[3])>1) or (len(splitline[4])>1)):     #skip indels
                    continue
                elif (splitline[4] is '.'): #homo ref
                    try:
                        all_calls['homos_ref'][splitline[1]]['callby'] = 2
                    except KeyError:
                        all_calls['homos_ref'][splitline[1]]=build_info(1, splitline[3], splitline[3])

                else:
                    geno = splitline[9][0:4] #if kid has alt allele
                    if PREFIX_HET == geno: #het
                        try:
                            all_calls['hets'][splitline[1]]['callby'] = 2
                        except KeyError:
                            all_calls['hets'][splitline[1]]=build_info(1, splitline[3], splitline[4])

                    elif PREFIX_HOMO_ALT == geno : #homo alt
                        try:
                            all_calls['homos_alt'][splitline[1]]['callby'] = 2
                        except KeyError:
                            all_calls['homos_alt'][splitline[1]]=build_info(1, splitline[3], splitline[4])


#go through variable sites file and note variable in each sub-dictionary
#TODO(SW): Implement a test for file type, gzip or not
print('Get variable sites')
sys.stdout.flush()
# with open(variable_site_file, 'r') as ref_data:
with gzip.open(variable_site_file, 'rb') as ref_data:
    for line in ref_data:      #read line in file
        if line.startswith(chr):            #skip all initial info lines
            splitline=line.split()
            if(len(splitline[3])>1 or len(splitline[4])>1):     #skip indels
                continue
            else:
                if splitline[1] in all_calls['hets']:  #found by us as het
                    #dictionaries are non-unique because a site found by the trio as one type could be found by single as another
                    all_calls['hets'][splitline[1]]['variable'] = 1        #note prev found
                    if all_calls['hets'][splitline[1]]['alt'] is not splitline[4]:        #check alt allele
                        all_calls['hets'][splitline[1]]['altalt'] = 1
                if splitline[1] in all_calls['homos_alt']:  #found by us as homoalt
                    all_calls['homos_alt'][splitline[1]]['variable'] = 1        #note prev found
                    if all_calls['homos_alt'][splitline[1]]['alt'] is not splitline[4]:        #check alt allele
                        all_calls['homos_alt'][splitline[1]]['altalt'] = 1

#go through exome data and note if found in exome
if exome is not '0':
    print ('Pares exome file %s' % exome)
    sys.stdout.flush()
    #filein = open('chr'+chr+'Ex_'+person+'.pileups', 'r')        #check if in exome
    with open(exome, 'r') as filein:        #check if in exome
        for site in filein:
            splitline = site.split()
            if len(splitline)==6:
                chr,pos,ref,count,bases,qual = site.split()
                if (pos in all_calls['hets'] and count is not '0'):
                    all_calls['hets'][pos]['exome'] = 1
                if (pos in all_calls['homos_alt'] and count is not '0'):
                    all_calls['homos_alt'][pos]['exome'] = 1
                if (pos in all_calls['homos_ref'] and count is not '0'):
                    all_calls['homos_ref'][pos]['exome'] = 1


label=('pos', 'ref', 'alt', 'As', 'Cs', 'Gs', 'Ts', 'callby', 'snp', 'snpdif', 'Ex')
title1= '\t'.join(label)+'\n'
label= ('pos', 'ref', 'alt', 'refs', 'alts', 'e1s', 'e2s', 'callby', 'snp', 'snpdif', 'Ex')
title2= '\t'.join(label)+'\n'

for filename in outfile_names:
    print ('Process %s' % filename)
    sys.stdout.flush()
    prefix = result_dir+'counts_'+filename+'_'+person

    fileout1 = open(prefix+'_flag.txt','w')
    fileout1.write(title1)

    fileout2 = open(prefix+'_byref.txt','w')
    fileout2.write(title2)

    with gzip.open(pileup_file, 'rb') as pileups:
        for line in pileups:
            splitline = line.split()
            if len(splitline)==6:
                chr,pos,ref,count,bases,qual = line.split()
                if (pos in all_calls[filename] and ref is not 'N'):
                    basecounts=dict()
                    basecounts['ref']=bases.count('.')+bases.count(',')
                    basecounts['A']=bases.count('A')+bases.count('a')
                    basecounts['C']=bases.count('C')+bases.count('c')
                    basecounts['G']=bases.count('G')+bases.count('g')
                    basecounts['T']=bases.count('T')+bases.count('t')
                    basecounts[ref]=basecounts['ref']
                    basecounts_summary = '%d\t%d\t%d\t%d' % (basecounts['A'], basecounts['C'], basecounts['G'], basecounts['T'])

                    fileout1.write(pos+'\t'+ref+'\t'+all_calls[filename][pos]['alt']+'\t'+
                                   basecounts_summary +'\t'+
                                   summary_info(all_calls[filename][pos])+'\n')

                    baselist=list()
                    baselist.append(basecounts[ref])
                    del basecounts[ref]
                    del basecounts['ref']
                    if all_calls[filename][pos]['alt'] is not ref:
                        baselist.append(basecounts[all_calls[filename][pos]['alt']])
                        del basecounts[all_calls[filename][pos]['alt']]
                    for k,v in basecounts.iteritems():
                        baselist.append(v)
                    baselist_summary = '%d\t%d\t%d\t%d' % (baselist[0], baselist[1], baselist[2], baselist[3])
                    fileout2.write(pos+'\t'+ref+'\t'+all_calls[filename][pos]['alt']+'\t'+
                                   baselist_summary+'\t'+
                                   summary_info(all_calls[filename][pos])+'\n')

    fileout1.close()
    fileout2.close()
