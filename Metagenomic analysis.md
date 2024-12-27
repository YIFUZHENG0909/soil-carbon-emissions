1.Quality filtering

trim Galore (v0.6.10)

```
#!/bin/bash
for i in $(cat sampleID); do
/public/home/wuq8022600160/anaconda3/envs/genomic/bin/trim_galore -o 01.cleandata --gzip --paired 00.Rawdata/${i}_val_1.fq.gz  00.Rawdata/${i}_val_2.fq.gz
done
#quality check
for i in $(cat sampleID); do
fastqc -t 32 -o QC -f fastq 01.cleandata/${i}_val_1.fq.gz
fastqc -t 32 -o QC -f fastq 01.cleandata/${i}_val_2.fq.gz
done
```

2.Assembly

MEGAHIT (v1.2.9)

```
#!/bin/bash
mkdir 03.Assembly_megahit
for i in $(cat sampleID); do
megahit -1 01.cleandata/${i}_val_1.fq.gz -2 01.cleandata/${i}_val_2.fq.gz -o 03.Assembly_megahit/${i} -t 64 --tmp-dir tmp_megahit --presets meta-large
done
##quality check
mkdir quast_megahit
for i in $(cat sampleID); do
/public/home/wuq8022600160/anaconda3/envs/quast/bin/quast.py -1 01.cleandata/${i}_val_1.fq.gz -2 01.cleandata/${i}_val_2.fq.gz -o quast_megahit/${i} -t 64 03.Assembly_megahit/${i}/final.contigs.fa
done
##remove short contigs (<500 bp)
mkdir 03.Assembly_megahit/contig_500
for i in $(cat sampleID); do
seqtk seq -L 500 03.Assembly_megahit/${i}/final.contigs.fa > 03.Assembly_megahit/contig_500/${i}.fa
done
```

3.predict the open reading frame (ORF)

Prodigal (v2.6.3) 

```
#!/bin/bash
mkdir 04.prodigal
for i in $(cat sampleID); do
/public/home/wuq8022600160/anaconda3/envs/checkm/bin/prodigal -i 03.Assembly_megahit/contig_500/${i}.fa -d 04.prodigal/${i}.fasta
done
```

4.obtain nonredundant gene sets

CD-HIT (v4.8.1) 

```
#!/bin/bash
cat 04.prodigal/*.fasta > 04.prodigal/all.fa
/public/home/wuq8022600160/cd-hit/cd-hit-est -c 0.95 -d 400 -T 64 -M 0 -n 8 -aS 0.9 -G 0 -g 1 -i 04.prodigal/all.fa -o 04.prodigal/uni_gene.fa
```

5.taxonomy annotation

Kraken2 (v 2.1.2) 

EasyMicrobiome (https://github.com/YongxinLiu/EasyMicrobiome)

```
# Generate report in default taxid output
mkdir 05.annotation/kraken2
mkdir 05.annotation/kraken2/temp
kraken2 --db /mnt/hpc/home/wuq8022600160/db/kraken2/pluspf \
04.prodigal/uni_gene.fa \
  --threads 32 \
  --report 05.annotation/kraken2/temp/NRgene.report \
  --output 05.annotation/kraken2/temp/NRgene.output
# Genes & taxid list
grep '^C' 05.annotation/kraken2/temp/NRgene.output|cut -f 2,3|sed '1 i Name\ttaxid' \
  > 05.annotation/kraken2/temp/NRgene.taxid
# Add taxonomy
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$0} NR>FNR{print $1,a[$2]}' \
  /mnt/hpc/home/wuq8022600160/zyf/tools/EasyMicrobiome/kraken2/taxonomy.txt \
  05.annotation/kraken2/temp/NRgene.taxid \
  > 05.annotation/kraken2/nucleotide.tax
```

6.function annotation

emapper-2.1.12 / eggNOG DB version: 5.0.2 /diamond version 2.1.6 / MMseqs2 version 14.7e284/ Compatible novel families DB version: 1.0.1

```
#!/bin/bash
mkdir 05.annotation/eggnog
python3 /public/home/wuq8022600160/anaconda3/envs/eggnog/bin/emapper.py -i 04.prodigal/uni_gene.fasta --itype CDS --translate --cpu 64 --data_dir /public/home/wuq8022600160/xulei/metagenome/eggnog/ --dmnd_db /public/home/wuq8022600160/xulei/metagenome/eggnog/eggnog_proteins.dmnd -o 05.annotation/eggnog/FUNCTION
done
```

7.gene quantitfy

```
mkdir -p temp/salmon
salmon -v # 1.10.1
# Build index, -t sequence, -i index
mkdir 06.quantify/salmon
mkdir 06.quantify/salmon/index 
salmon index -t 04.prodigal/uni_gene.fasta \
-i 06.quantify/salmon/index
for j in $(cat sampleID); do
salmon quant -i 06.quantify/salmon/index -l A --meta \
    -1 01.cleandata/${j}_val_1.fq.gz -2 01.cleandata/${j}_val_2.fq.gz \
    -o 06.quantify/salmon/${j}.quant
done
# Quantitative, l library type automatic selection, p thread, --meta metagenomic model
# Merge
mkdir 06.quantify/salmon/result
salmon quantmerge --quants 06.quantify/salmon/*.quant \
    -o 06.quantify/salmon/result/gene.TPM
salmon quantmerge --quants 06.quantify/salmon/*.quant \
    --column NumReads -o 06.quantify/salmon/result/gene.count
sed -i '1 s/.quant//g' 06.quantify/salmon/result/gene.*
```



