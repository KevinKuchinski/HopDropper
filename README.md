# HopDropper
HopDropper is a bioinformatic tool that analyzes paired-end genomic sequencing data to: 1) identify distinct fragments of genomic material, and 2) generate pairs of consensus sequences that describe both ends of each fragment. This single-molecule level of resolution allows genomic fragments to be counted and individually characterized. As part of this analysis, HopDropper also identifies chimeric sequence artefacts generated during library construction and targeted enrichment. Removal of these artefacts improves confidence in the genomic fragments identified by HopDropper and reduces false detections.

HopDropper does the above using unique molecular indices (UMIs). UMIs are distinct nucleotide sequences that become associated with each end of each genomic fragment during library construction. These UMIs are introduced before amplification and target enrichment, so all subsequent molecular copies of the same template fragment will share the same pair of UMIs. When using paired-end sequencing chemistry, this allows all read pairs describing the same template fragment to be grouped together. These groups of read pairs can then be collapsed into a single pair of consensus sequences describing both ends of a single template fragment.

HopDropper UMIs are composed of an extrinsic portion (eUMI) and an intrinsic portion (iUMI). eUMIs are incorporated into specialized sequencing adapter oligomers that are ligated onto input genomic material during library construction. Some adapters contain random eUMI sequences while others contain one eUMI from a fixed repertoire. iUMIs are naturally created by the random fragmentation of genomic material which creates distinct sequence motifs at the ends of the fragments.

HopDropper also uses UMIs to discard chimeric molecules created during library construction and target enrichment. Each UMI should only co-occur with one other UMI (its so-called best mate). Infrequent co-occurrences between UMIs that are not each other’s best mates are removed. Similarly, HopDropper uses UMIs to discard index hops. Index hops are a special type of chimera wherein genomic material originating from one specimen acquires library barcodes from a different specimen. Each UMI should only appear in one library (its so-called best library). Infrequent appearances of a UMI in libraries other than its best library are removed by HopDropper.

## HopDropper input
HopDropper analysis is performed on a group of FASTQ files. This group contains the first and second read files (generated by paired-end sequencing) for one or multiple libraries. If multiple libraries are analyzed, they should form a group of libraries among which there was a risk of chimera formation and index hopping, e.g. libraries pooled for PCR or hybridization probe capture-based enrichment. 

The FASTQ file pairs included in a HopDropper analysis are indicated in a TSV input file. This file associates each FASTQ file pair with a library name. The input TSV file should not have a header line, and each line must be ordered as follows:

```library_name R1_file_name    R2_file_name```

FASTQ file pairs must be unzipped and located in the same directory as the TSV input file.

## HopDropper Algorithm
HopDropper performs the following sequence of operations on the provided FASTQ files:
1.	HopDropper extracts a UMI from each read. The user specifies the length of eUMI present in the adapter oligo and the length of the desired iUMI. HopDropper analyzes the first <i>n</i> base quality scores of each read, where <i>n</i> is the combined eUMI and iUMI length. If all of those base quality scores exceed a minimum value (default is PHRED score of 30), the first <i>n</i> bases of the read sequence are extracted as its UMI.

2.	HopDropper assigns a best library to each UMI. A UMI’s best library is the library in which it occurs most frequently. If a UMI occurs in more than one library, the number of occurrences in its most frequent library must be significantly greater than the number of occurrences in its second most frequent library. Significance is determined by converting frequency counts into proportions then comparing the proportions of the top two libraries with the two proportion Z-test (default alpha level is 5%). If a UMI occurs in only one library, the two proportion Z-test compares the proportions 1.0 and 0.0, and the significance relies solely on the sample size, i.e. the number of times the UMI was observed. If a UMI does not have a clear best library after this analysis, it is assigned no best library and read pairs containing this UMI will be discarded when valid reads are outputted. 

3.	HopDropper assigns a best UMI mate to each UMI. A UMI’s best mate is the other UMI with which it co-occurs in read pairs the most frequently. If a UMI co-occurs with more than one other UMI, the number of co-occurrences with its most frequent mate must be significantly greater than the number of co-occurrences with its second most frequent mate. Significance is determined by converting frequency counts into proportions then comparing the proportions of the top two UMI mates with the two proportion Z-test (default alpha level is 5%). If a UMI co-occurs with only one other UMI, the two proportion Z-test compares the proportions 1.0 and 0.0, and the significance relies solely on the sample size, i.e. the number of times the UMI was observed. If a UMI does not have a clear best mate after this analysis, it is assigned no best mate and read pairs containing this UMI will be discarded when valid reads are outputted.

4.	HopDropper writes valid read pairs to new FASTQ files. For each read pair, HopDropper considers the UMIs for both reads. The pair is only valid if both UMIs have the same best library and that library is the one to which the read pair was assigned. The read pair is also valid only if both UMIs are each other’s best UMI mate. Finally, the read pair is only valid if its UMI pair has been observed a minimum number of times (default is 30 occurrences).

5.	HopDropper creates pairs of consensus sequences representing both ends of each distinct input genomic molecule. All read pairs in a library are grouped based on their UMI pair. Each group of read pairs is randomly sub-sampled for faster computation (default is up to 200 read pairs). Low quality bases in these reads are masked with ‘N’ ambiguity characters (default threshold is a PHRED score of 30), then a multiple sequence alignment is conducted for each read in the pair. These multiple sequence alignments are used to generate consensus sequences for both reads in the pair, i.e. both ends of the input genomic molecule described by that group of read pairs. If a position in the multiple sequence alignment is not described by a sufficient number of non-ambiguous bases (default is 20 non-ambiguous bases), then it is masked with an ‘N’ ambiguity character. If the non-ambiguous bases describing a position in the multiple sequence alignment do not converge on a single nucleotide (default consensus threshold is 0.75), then it is masked with an ‘N’ ambiguity character. Consensus sequences are truncated after the first appearance of a run of ‘N’ ambiguity character (default run length is 5).

6.	HopDropper writes each pair of consensus sequences to a FASTA file. Headers are formatted as follows: ```>experiment_name|library_name|fragment_number|fragment_end|fragment_copies|fragment_UMI_pair``` 

Experiment_name is a name given to a HopDropper analysis at runtime. Library_name is the name of the library in which the fragment was detected.      Fragment_number is a sequential identifier for all fragments from the same library. The combination of library_name and fragment_number is unique for each fragment within experiment_name. Fragment_copies is the number of read pairs describing the fragment, i.e. the number of PCR copies of the input genomic molecule that were sequenced. Fragment_UMI_pair is the UMI pair describing the fragment. This is useful for identifying the same fragment across different experiments (when the same fragment may be assigned a different fragment_number within its library).

## Important Notes
A.	Do not trim bases from the beginning of reads! HopDropper must consider those bases as-sequenced for UMI assignment.

B.	UMIs are not always unique. Sequencing adapter oligomers with fixed eUMIs are obviously limited in number. Sequencing adapter oligomers with random eUMIs have a limited number of permutations. iUMIs are limited by the number of positions in the genome where fragmentation can occur, and repetitive genome structures or simple chance may result in redundant sequence motifs at certain fragmentation locations. Taken together, this means that some input fragments may share UMIs with other input fragments on one or both of their ends by chance. This coincidental sharing of UMIs is called a UMI collision.
The rate of UMI collisions depends on the diversity of eUMIs available on adapter molecules. It also depends on the diversity of iUMIs residing within the genomes of the organisms in the specimen which, in turn, depends on total genome length and the presence and frequency of repetitive structural elements. The UMI collision rate also depends on the number of input molecules; fewer input molecules require fewer UMIs, so there is a lower probability of collision. UMI collision rates are easily underestimated, even when the diversity of eUMIs and iUMIs seems extensive (the same cognitive bias behind the classic “birthday problem” from statistics and probability theory applies here).
Unfortunately, UMI collisions are alike in appearance to chimeras and index hops. HopDropper does not attempt to distinguish them. If a UMI collision occurs and one of the input fragments containing that UMI is significantly more abundant, the other input fragment will be discarded. If a UMI collision occurs and none of the input fragments containing that UMI are significantly more abundant than the rest, all of these input fragments will be discarded. This means that conditions favouring UMI collisions will increase data attrition by HopDropper, potentially resulting in false negatives. Users should understand that HopDropper’s priority is to remove false positives.

C.	HopDropper expects the eUMIs to be in-line with the input molecule, i.e. directly adjacent to the end of the input molecule. If the eUMI’s position within the sequencing adapter oligomer does not place it directly adjacent to the input molecule upon ligation, eUMI sequences and their corresponding base quality scores must be prepended to read entries in the FASTQ files.

D.	Users may provide a list of fixed eUMIs if there is a specific repertoire of eUMI sequences in the sequencing adapter oligomers they use. If this list is provided, HopDropper will reject any read whose eUMI does not appear in the list. 

E.	Users may choose not to use an eUMI, relying solely on iUMIs. This requires a library construction method where input genomic material is randomly fragmented before any adapters are added and before any amplification occurs. Thus, amplicon enrichment is not conducive for analysis by HopDropper using iUMIs alone.

F.	Users may choose not to use an iUMI, relying solely on eUMIs in library adapter oligos. If there are a limited number of eUMIs provided, this may increase the rate of UMI collisions and, consequently, data attrition.

## HopDropper usage
```
$ python hop_dropper.py -i <input file> -o <output dir> -I <iUMI length> -E <eUMI length> [-q <min UMI qual> -f <fixed UMIs file> -m <min UMI count> -Q <min read qual> -s <subsample size> -d <min read depth> -c <consensus threshold> -N <max consecutive Ns>]
```
<b>Arguments:</b>

    -i : path to TSV input file
    -o : output directory
    -I : iUMI length
    -E : eUMI length
    -q : minimum PHRED quality for bases in UMI (default = 30)
    -f : path to TXT file listing fixed eUMIs (one eUMI per line)
    -m : minimum occurences of UMI pair for retention (UMI pairs that do not appear at least this many times are not outputted, default = 30)
    -Q : minimum PHRED quality of bases in reads (bases below this quality are masked with Ns while generating fragment end consensus sequences, default = 30)
    -s : sub-sample size (the maximum number of read pairs used for generating fragment end consensus sequences, default = 200)
    -d : minimum read depth (the minimum number of reads with non-ambiguous nucleotides required for generating fragment end consensus sequences, default = 20)
    -c : consensus threshold (the proportion of reads that must contain the same nucleotide when generating fragment end consensus sequences, default = 0.75)
    -N : maximum consecutive Ns (reads are trimmed after this number of consecutive Ns, default = 5)
