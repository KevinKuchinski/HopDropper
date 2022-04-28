import argparse
from collections import Counter

parser=argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True, help="List of file names of interleaved FASTQ data to analyze")
parser.add_argument("-l", "--len", type=int, required=True, help="Length of UMIs")
parser.add_argument("-q", "--qual", type=int, required=True, help="Minimum PHRED base quality within UMI")
parser.add_argument("-m", "--min", type=int, required=True, help="Minimum number of reads pairs containing UMI pair")
parser.add_argument("-o", "--output", type=str, required=True, help="Job name to append to hop_dropper output files")
args=parser.parse_args()

'''
hop_dropper analyzes interleaved FASTQ files and removes index-hopped read pairs. 
It identifies all reads originating from the same template molecule by creating an intrinsic unique molecular identifier (iUMI) 
from the first n bases of the read sequence (assuming all n bases have a minimnum PHRED score). hop_dropper filters each provided 
FASTQ file to remove read pairs where iUMIs from either read are clearly more abundant in another FASTQ file. To be considered clearly 
most abundant in a particular library, the instances of an iUMI in that library must be at least twice the number of instances in any 
other library (ie one PCR cycle). hop_dropper also removes read pairs where the total instances of either read's iUMI do not exceed a 
certain number, eg to remove singletons. 
'''

UMILength=args.len ## The number of bases at the start of each read that is taken as the UMI for that read
minUMIQual=args.qual ## The minimum phred score all bases in the UMI must have for the read not to be rejected

## Create a list containing the file names of the reads to be analyzed
with open(args.input) as inputFile:
    fileList=inputFile.readlines()
fileList=[fileName.rstrip() for fileName in fileList if fileName.rstrip()!='']

R1UMIs=[]
R2UMIs=[]
print()
print('Reading iUMIs in input files...')
for fileName in fileList:
    print('  ',fileName,'...')
    with open(fileName) as inputFASTQ:
        inputFASTQLines=inputFASTQ.readlines()
    read1Seqs=inputFASTQLines[1::8]
    read1Quals=inputFASTQLines[3::8]
    read2Seqs=inputFASTQLines[5::8]
    read2Quals=inputFASTQLines[7::8]
    for read1Seq,read2Seq,read1Qual,read2Qual in zip(read1Seqs,read2Seqs,read1Quals,read2Quals):
        UMI1Seq=read1Seq.rstrip()[0:UMILength]
        UMI2Seq=read2Seq.rstrip()[0:UMILength]
        UMI1Qual=read1Qual.rstrip()[0:UMILength]
        UMI2Qual=read2Qual.rstrip()[0:UMILength]
        if all([ord(character)-33>=minUMIQual for character in UMI1Qual])==True:
            if all([ord(character)-33>=minUMIQual for character in UMI2Qual])==True:
                R1UMIs.append(UMI1Seq)
                R2UMIs.append(UMI2Seq)

R1UMIs=set(R1UMIs) ## Dereplicate UMI lists
R2UMIs=set(R2UMIs) 

## Initialize empty lists
R1UMICount={R1UMI:0 for R1UMI in R1UMIs}
R2UMICount={R2UMI:0 for R2UMI in R2UMIs}
R1UMILibCount={R1UMI:[] for R1UMI in R1UMIs}
R2UMILibCount={R2UMI:[] for R2UMI in R2UMIs}

print()
print('Counting occurences of iUMIs in each library ...')
for fileName in fileList:
    print('  ',fileName,'...')
    with open(fileName) as inputFASTQ:
        inputFASTQLines=inputFASTQ.readlines()
    read1Seqs=inputFASTQLines[1::8]
    read1Quals=inputFASTQLines[3::8]
    read2Seqs=inputFASTQLines[5::8]
    read2Quals=inputFASTQLines[7::8]
    for read1Seq,read2Seq,read1Qual,read2Qual in zip(read1Seqs,read2Seqs,read1Quals,read2Quals):
        UMI1Seq=read1Seq.rstrip()[0:UMILength]
        UMI2Seq=read2Seq.rstrip()[0:UMILength]
        UMI1Qual=read1Qual.rstrip()[0:UMILength]
        UMI2Qual=read2Qual.rstrip()[0:UMILength]
        if all([ord(character)-33>=minUMIQual for character in UMI1Qual])==True:
            if all([ord(character)-33>=minUMIQual for character in UMI2Qual])==True:
                R1UMICount[UMI1Seq]+=1
                R2UMICount[UMI2Seq]+=1
                R1UMILibCount[UMI1Seq].append(fileName)
                R2UMILibCount[UMI2Seq].append(fileName)

print()
print('Determining top library for each iUMI...')
print('  ', end='')
UMIReport=open(args.output+'_UMI_report.tsv','w')
UMIReport.write('UMI\tUMI_count\tTop_lib\tTop_lib_count\tSecond_lib\tSecond_lib_count\tBest_lib\n')

R1UMIBestLib={R1UMI:[] for R1UMI in R1UMIs}
R2UMIBestLib={R2UMI:[] for R2UMI in R2UMIs}

i=0
for UMI in R1UMIs:
    i+=1
    if i%100==0:
        print('.', end='')
    ## Determine the libraries in which this UMI makes its most and second-most appearances
    libCount=Counter(R1UMILibCount[UMI])
    topLib=libCount.most_common(1)[0][0]
    topLibCount=libCount.most_common(1)[0][1]
    secondLib=libCount.most_common(2)[-1][0]
    secondLib=secondLib if secondLib!=topLib else None
    ## Determine if the library in which this library makes its most appeareances is the best choice (by having at least twice the appearances as the second-most common library)
    if secondLib!=None:
        secondLibCount=libCount.most_common(2)[-1][1]
    else:
        secondLibCount=0
    if topLibCount>secondLibCount*2:
        R1UMIBestLib[UMI]=topLib
    else:
        R1UMIBestLib[UMI]=None
    UMIReport.write(UMI+'\tR1\t'+str(R1UMICount[UMI])+'\t'+topLib+'\t'+str(topLibCount)+'\t'+str(secondLib)+'\t'+str(secondLibCount)+'\t'+str(R1UMIBestLib[UMI])+'\n')

for UMI in R2UMIs:
    i+=1
    if i%100==0:
        print('.', end='')
    ## Determine the libraries in which this UMI makes its most and second-most appearances                                                                                                                
    libCount=Counter(R2UMILibCount[UMI])
    topLib=libCount.most_common(1)[0][0]
    topLibCount=libCount.most_common(1)[0][1]
    secondLib=libCount.most_common(2)[-1][0]
    secondLib=secondLib if secondLib!=topLib else None
    ## Determine if the library in which this library makes its most appeareances is the best choice (by having at least twice the appearances as the second-most common library)                           
    if secondLib!=None:
        secondLibCount=libCount.most_common(2)[-1][1]
    else:
        secondLibCount=0
    if topLibCount>secondLibCount*2:
        R2UMIBestLib[UMI]=topLib
    else:
        R2UMIBestLib[UMI]=None
    UMIReport.write(UMI+'\tR2\t'+str(R2UMICount[UMI])+'\t'+topLib+'\t'+str(topLibCount)+'\t'+str(secondLib)+'\t'+str(secondLibCount)+'\t'+str(R2UMIBestLib[UMI])+'\n')

UMIReport.close()
print()

print()
print('Filtering read pairs in input files and writing to separate FASTQs for valid, chimeric, and index hopped read pairs...')
readsReport=open(args.output+'_reads_report.tsv','w') ## Create a TSV file describing the number of reads, read pairs, and UMI pairs in each file
readsReport.write('Lib\tTotal_R1\tTotal_R2\tHi_qual_R1\tHi_qual_R2\tTotal_pairs\tHi_qual_pairs\tValid_Pairs\tHop_pairs\tUniq_valid_pairs\tUniq_hop_pairs\n')

for fileName in fileList:
    print('  ',fileName,'...')
    with open(fileName) as inputFASTQ:
        inputFASTQLines=inputFASTQ.readlines()
    read1Names=inputFASTQLines[0::8]
    read1Seqs=inputFASTQLines[1::8]
    read1Quals=inputFASTQLines[3::8]
    read2Names=inputFASTQLines[4::8]
    read2Seqs=inputFASTQLines[5::8]
    read2Quals=inputFASTQLines[7::8]
    validFile=open(''.join(fileName.split('.')[:-1])+'_valid_'+args.output+'.fq','w')
    hopFile=open(''.join(fileName.split('.')[:-1])+'_hop_'+args.output+'.fq','w')
    totalR1s=0
    totalR2s=0
    totalReadPairs=0
    hiQualR1s=0
    hiQualR2s=0
    hiQualReadPairs=0
    validReadPairs=0
    validUMIPairs=[]
    hopReadPairs=0
    hopUMIPairs=[]
    for read1Name,read2Name,read1Seq,read2Seq,read1Qual,read2Qual in zip(read1Names,read2Names,read1Seqs,read2Seqs,read1Quals,read2Quals):
        totalR1s+=1
        totalR2s+=1
        totalReadPairs+=1
        UMI1Seq=read1Seq.rstrip()[0:UMILength]
        UMI2Seq=read2Seq.rstrip()[0:UMILength]
        UMI1Qual=read1Qual.rstrip()[0:UMILength]
        UMI2Qual=read2Qual.rstrip()[0:UMILength]
        R1HiQual=True if all([ord(character)-33>=minUMIQual for character in UMI1Qual]) else False
        if R1HiQual==True:
            hiQualR1s+=1
        R2HiQual=True if all([ord(character)-33>=minUMIQual for character in UMI2Qual]) else False
        if R2HiQual==True:
            hiQualR2s+=1
        if R1HiQual==True and R2HiQual==True:
            hiQualReadPairs+=1
            if all([R1UMICount[UMI1Seq]>args.min,R2UMICount[UMI2Seq]>args.min]):
                if all([R1UMIBestLib[UMI1Seq]==fileName,R2UMIBestLib[UMI2Seq]==fileName]):
                    validReadPairs+=1
                    validUMIPairs.append(tuple(sorted([UMI1Seq,UMI2Seq])))
                    validFile.write(read1Name)
                    validFile.write(read1Seq)
                    validFile.write('+\n')
                    validFile.write(read1Qual)
                    validFile.write(read2Name)
                    validFile.write(read2Seq)
                    validFile.write('+\n')
                    validFile.write(read2Qual)
                else:
                    hopReadPairs+=1
                    hopUMIPairs.append(tuple(sorted([UMI1Seq,UMI2Seq])))
                    hopFile.write(read1Name)
                    hopFile.write(read1Seq)
                    hopFile.write('+\n')
                    hopFile.write(read1Qual)
                    hopFile.write(read2Name)
                    hopFile.write(read2Seq)
                    hopFile.write('+\n')
                    hopFile.write(read2Qual)
    validFile.close()
    hopFile.close()
    validUMIPairs=len(set(validUMIPairs))
    hopUMIPairs=len(set(hopUMIPairs))
    readsReport.write(fileName+'\t')
    readsReport.write(str(totalR1s)+'\t'+str(totalR2s)+'\t'+str(hiQualR1s)+'\t'+str(hiQualR2s)+'\t'+str(totalReadPairs)+'\t')
    readsReport.write(str(hiQualReadPairs)+'\t'+str(validReadPairs)+'\t'+str(hopReadPairs)+'\t'+str(validUMIPairs)+'\t'+str(hopUMIPairs)+'\n')

readsReport.close()

print()
print('Done.')
