# Third year project
# Antimicrobial resistance gene composition of UK environmental isolates of E. coli: The genetic context of AMR genes.
![alt text](https://github.com/BadgerRob/Ecoli_assembly_Jack/blob/master/Vic_plasmid_mapping.png)
## NB The python script called ‘wrench.py’ is an automated assembly pipeline.

## Introduction

Escherichia coli is a common bacterial species that has large diversity at the strain level. The genetic and phenotypic variation among these strains lead to a large diversity of functional niches with different strains having vastly different pathogenicity, virulence and environmental requirements. Some strains of E coli can be highly pathogenic to humans (e.g. Enterohemorrhagic Escherichia coli O157:H7) and a major food borne pathogen, while others fulfil a commensal role within the human gut. Furthermore, certain linages (e.g. ST131, clade V etc) have members that are able to persist, replicate and survive in the environment and are often referred to as Extra Intestinal Pathogenic E. coli (ExPEC). ExPEC strains are a substantial threat to human health due to their ability to survive in the environment and their potential pathenogisity to humans. ExPEC E. coli are responsible for the majority of UTI infections in Low to middle income countries, where sanitation and waste management can be poor or non-existent. 

The presence of antimicrobial resistance in pathogenic ecoil is a confounding factor when considering the threat of these strains with respect to human heath and morbidity. Antimicrobial (Antibiotic) resistance originates in the environmental settings where competition between species of microorganisms is in flux. The constant selection for antibiotic compounds and their corresponding resistance mechanisms can be described as an 'evolutionary arms race' where competition and selection drives the diversity of both antibiotic variation and resistance genes. As selection is an optimising process, resistance genes under selection are expected to move to fixation within a population and thus become synonymous with specific strains and their chromosomal context.

Many E. coli strain in a community setting contain plasmids that consist of extra chromosomal DNA that can frequently move between bacterial hosts. Such plasmids can be thought of as independent units of selection that contain functional genes that can move horizontally through a community. While these plasmids can contain numerous genes with function, the presence of anti-microbial resistance genes is of particular concern to human health due to their rapid dissemination through a natural community under selection, and furthermore, the chance of re-entering the human population. 

The interplay between the environment, selection and the genetic composition of resistance genes within a bacterial community is complex and remains far from resolved. The additional anthropogenic selection pressures on the community acquired resistance profile, often referred to as ‘the resistome’, therefore remain elusive. It is hypothesised that the increased use of clinically prescribed antibiotics has resulted an increase of antimicrobial compounds being released into the environment. One such source of these compounds is human waste and waste water effluents entering the UK water system. While it is likely that an increase in antimicrobial compounds in the environment is having a substantial effect on the selective pressure of a the microbial communities exposed to such compounds, the mechanistic response is difficult and tenuous to ascertain due to the inability to describe the direction of travel of AMR genes between the clinic and the environment. 

One such difficulty in quantifying the risk of the environmental resistome, and the consequences of anthropogenic selection, is the restriction in accurately partitioning the mobile elements of a resistome; termed the mobilome. Defining the genomic context of individual genes in a sequenced isolate is posable using short read technology, however defining the wider genetic context of such genes, i.e. the entire plasmid / chromosomal structure that the gene of interest resides on, is limited by the length of individual sequencing reads.


Imagen you have an E. coli strain which contains five large environmental plasmids of similar structure but differing evolutionary histories. The approximate size of a environmental E. coli strain is 5 mega bases and the approximate size of these plasmids is 150-200 kb. When sequencing this isolate you would extract the entire DNA from multiple bacterial cells in a culture, fragment the DNA to ~150 base pairs, sequence everything and try and assemble an unknown number of plasmids and a genome using many overlapping pairs of reads of 150 base pairs in length. This is normally done using a "k-mer clustering" approach.

As an analogy, it is like trying to put six different jigsaws back together, only each complete jigsaw comes from a million different copies of that jigsaw, cut up in different ways, all mixed up together with the other five jigsaws. To make things more interesting, the picture on the jigsaws are very similar. What you end up with is one or two big jigsaws which are representative of all six jigsaws together. This makes figuring out which jigsaw parts belong to which individual jigsaws quite a challenge. One way to make things easier is to use larger jigsaw pieces... really large jigsaw pieces... DNA reads… OK I got carried away with the jigsaw analogy.        


Long read sequencing provides a means to directly sequence long bits of extracted DNA from a sample. Short read sequencing use a sequencing by synthesis approach which sequences a complementary strand as it is synthesised from the target. The length of such sequences are limited by the chemistry used in this method. Therefore the ability to directly sequence long bits of DNA removes the length limitations of sequencing by synthesis. There are many advantages and disadvantages to both technologies and I encourage you to explore this further in your discussion. 

When assembling a genome without a reference genome we undertake what is known as 'de-novo assembly.' There are a number of programs that are available to assemble long read genomes (e.g. miniasm, canu, flye etc.) and we will be using Flye for this pipeline. 

Onne important note is the error structure associated with the long read platform we are using. Thus while we are able to get very long sequences which produces a much more contiguous assembly using long reads, the accuracy of a final assembly is limited by the raw read accuracy. It is important to note that the accuracy of your final assembly will impact the validity of your results depending on the nature of your research question. For example, if you need to know the order of genes and their classes relative to a genetic position (e.g. a promotor type), then highly contiguous assemblies may be more important than individual base pair accuracy. However, if you are looking at SNP or variant calling then single base pair accuracy may be more important that sequence contiguity. In this study we require both highly contiguous assemblies to identify multiple plasmids and the genetic context of genes of interest, as well as individual base pair accuracy due to the need to resolve multiple isolates.

One way to improve this draft assembly is to 'polish' out the errors using the original raw long reads. This is reliant on having a random error structure in the raw reads and completing an alignment pile up against the draft assembly to calculate the most likely correct base for a given position. There are a number of tools we can use to do this (e.g. nanopolish, racon, medaka etc.) however we will be using racon and medaka in this pipeline. 

Another approach to improving a draft assembly is to use highly accurate short read illumina data to complement the assembly. This can be done in tandem with the long read assembly or used as a further polishing step to the long read draft. Unicycler and Pilon both have the capability of undertaking a 'hybrid assembly' and we will use Pilon for this pipeline.

The aim of this study is to resolve, assemble and annotate up to 20 E. coli isolates and their associated plasmids from a UK water system using a hybrid assembly approach. This will allow us to explore the mobility of antibiotic resistance genes in the environment which will elucidate the selective mechanisms for persistance and highlight potentual human health risks.

## Workflow
![alt text](https://github.com/BadgerRob/Ecoli_assembly_Jack/blob/master/Assembly.png)

A flowchart outlining the assembly and polishing process.

## Basecalling Fast5 files

Fast5 files store the raw data from the ONT sequencing platform. This data, turmed 'squiggle data' needs to be transformed into the fastq seqence data we know and love. This is undertaken using trained nural networks which predict the most likely base in a given context of the squiggle pattern (read up on this if interested, it is quite facinating, it's like having a 4-6 bp open reading frame). To do this we use the ONT program guppy set to use the high accracy configuration mode. I have completed this bit for you as it is computationaly intensive and takes some time. The commands were as follows.

```
guppy_basecaller -r --input_path path/to/fast5/dir --save_path path/to /out/dir --config dna_r9.4.1_450bps_hac.cfg  --qscore_filtering --min_qscore 7 --cpu_threads_per_caller 4 --num_callers 2 --flowcell FLO-MIN106 --kit SQK-LSK109
```
|Flag                      | Description               | 
| -------------------------|:-------------------------:| 
| `guppy_basecaller`       |calls Guppy                | 
| `-r`                     |recursive mode             | 
| `--input-path`           |path to fast5 dir/         |
| `--save-path`            |path to output fastq files |
| `--qscore_filtering`     |filter for quality score   |
| `--min_qscore`           |minimum qscore for pass    |
| `cpu_threads_per_caller` |number of threads          |
| `--num_callers`          |number of basecallers      |
| `--config`               |configuration file         |



Reads are output as `.fastq` files containing 4000 reads and quality data per file. Sequences contain a `@` followed by a header then sequence. This is separated from quality data by `+`. 

(_Optional_) If you to watch in real time how many sequences are being written you can change to the directory where your fastq files are being written (/pass) and enter the bash one-liner:

```

watch -n 10 'find . -name "*.fastq" -exec grep 'read=' -c {} \; | paste -sd+ | bc'

```

## Sanitation, Demultiplexing and trimming

Once reads have been basecalled a quick check and fix for broken reads is usually a good idea before further manipulations.

```
seqkit sana path/to/reads.fastq | filtlong –min_length 400 –keep_percent 95 –min_mean_q 8 - | qcat –f - --detect_middle --min_read-length 400 –trim –b path/to/barcode_out_dir/

```

A further round of trimming can be used to remove missed adapters and discard middle adaptors resulting from chimeric reads. Now the files have been split into barcode files, a loop can be used to speed things up.

```
for file in barcode*.fa
do

	stub=${file%.fa}
	stub=${stub#barcode}
	echo ${stub}
	porechop –i $file –o path/to/${stub}_porechop_out.fastq –-discard_middle –discard_unnassigned
done
```

A final step to sort by length can be used on each barcode dir for completion. This can also be made into a loop as before to process multiple files.

```
Seqkit sort –lr path/to/reads.fastq | gzip > reads.sorted.fastq.gz

```
## Draft assembly
Now the first draft assembly can be undertaken using Flye.

```
flye --nano-raw path/to/reads.sorted.gz -–plasmids --out-dir path/to /flye_assembly/ --genome-size 5mb  --threads 8

```

|Flag                      | Description               | 
| -------------------------|:-------------------------:| 
| `flye`	            |calls flye                 | 
| `-r`                     |recursive mode             | 
| --nano-raw  `           |path to fastq/             |
| -–plasmids             |retain plasmid drafts      |
| `--genome_size           |estimate of genome size    |
| `--threads`              |number of threads to run   |
