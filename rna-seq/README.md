## Get data

Paper: http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0093338

GEO: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54413

Put these files into a text file called `sra-files.txt`

```
# UVB
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450265/SRR1145042/SRR1145042.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450266/SRR1145043/SRR1145043.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450267/SRR1145044/SRR1145044.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450268/SRR1145045/SRR1145045.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450269/SRR1145046/SRR1145046.sra

# Control
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450270/SRR1145047/SRR1145047.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450271/SRR1145048/SRR1145048.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450272/SRR1145049/SRR1145049.sra
ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX450%2FSRX450273/SRR1145050/SRR1145050.sra
```

Then use wget to get all the files

```bash
wget sra-files.txt
```

Now use the sra-toolkit to unpack all the files.


```bash
find *sra | parallel fastq-dump {}
```

Use a perl script to convert all to Illumina 1.9

```bash
find *fastq | parallel "illumina15-to-sanger.pl {} > {}.fq"
```

FASTQC on all to make sure the quality encoding worked well.

```bash
find *.f*q | parallel fastqc {} --outdir .
```

```bash
mv SRR1145042.fastq.fq uvb1.fq
mv SRR1145043.fastq.fq uvb2.fq
mv SRR1145044.fastq.fq uvb3.fq
mv SRR1145045.fastq.fq uvb4.fq
mv SRR1145046.fastq.fq uvb5.fq

mv SRR1145047.fastq.fq ctl1.fq
mv SRR1145048.fastq.fq ctl2.fq
mv SRR1145049.fastq.fq ctl3.fq
mv SRR1145050.fastq.fq ctl4.fq


## Alignment
NCPUS=8
MEM="40GB"
GENOMEDIR="/home/sdt5z/genomes/star/hg19/"
QSUB="qsub -l select=1:ncpus=$NCPUS:mem=$MEM,walltime=24:00:00 -q uvabx -W group_list=uvabx -V -j oe -m bae -M vustephen+fir@gmail.com"
find `pwd` -name "*.fq" | sed 's/\.fq$//g' | sort | xargs -i echo $QSUB -- `which time` `which STAR` --genomeDir $GENOMEDIR --runThreadN $NCPUS --outFileNamePrefix {}.star. --readFilesIn {}.fq > runstar.sh

## Count
featureCounts -a ~/genomes/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf -o counts.txt -T 12 -t exon -g gene_id *.sam

## Extract chr4:60000000-160000000
NCPUS=8
MEM="40GB"
QSUB="qsub -l select=1:ncpus=$NCPUS:mem=$MEM,walltime=24:00:00 -q uvabx -W group_list=uvabx -V -j oe -m bae -M vustephen+fir@gmail.com"
find *.sam | sed 's/.star.Aligned.out.sam//g' | sort | parallel --dry-run 'samtools view -Sb {}.star.Aligned.out.sam | samtools sort -o - | bedtools intersect -abam - -b roi.bed | bedtools bamtofastq -i - -fq {}.fastq'

```
