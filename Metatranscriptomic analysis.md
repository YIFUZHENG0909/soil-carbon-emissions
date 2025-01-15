1.Quality filtering

trim Galore (v0.6.10)

```
#!/bin/bash
for i in $(cat sampleID); do
/public/home/wuq8022600160/anaconda3/envs/genomic/bin/trim_galore -o 01.cleandata --gzip --paired 00.Rawdata/${i}.R1.raw_val_1.fq.gz   00.Rawdata/${i}.R2.raw_val_2.fq.gz 
done
#quality check
for i in $(cat sampleID); do
fastqc -t 32 -o QC -f fastq 01.cleandata/${i}.R1.raw_val_1.fq.gz 
fastqc -t 32 -o QC -f fastq 01.cleandata/${i}.R2.raw_val_2.fq.gz 
done
```

2. remove rRNA

```
#INDEX
##workdir public/home/wuq8022600160/database/sortmerna
sortmerna --index 1 --threads 32 \
--ref ~/sortmerna4.3/rRNA_databases/silva-bac-16s-id90.fasta \
--ref ~/sortmerna4.3/rRNA_databases/silva-bac-23s-id98.fasta \
--ref ~/sortmerna4.3/rRNA_databases/silva-arc-16s-id95.fasta \
--ref ~/sortmerna4.3/rRNA_databases/silva-arc-23s-id98.fasta \
--ref ~/sortmerna4.3/rRNA_databases/silva-euk-18s-id95.fasta \
--ref ~/sortmerna4.3/rRNA_databases/silva-euk-28s-id98.fasta \
--ref ~/sortmerna4.3/rRNA_databases/rfam-5s-database-id98.fasta \
--ref ~/sortmerna4.3/rRNA_databases/rfam-5.8s-database-id98.fasta \
--workdir ./
#remove rRNA
mkdir 2.assembly-rna
for i in $(cat sampleID); do
mkdir 2.assembly-rna/${i}
mkdir 2.assembly-rna/${i}/idx
cp /public/home/wuq8022600160/database/sortmerna/idx/* 2.assembly-rna/${i}/idx
sortmerna --index 0 --threads 32 \
--ref /public/home/wuq8022600160/database/sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta \
--ref /public/home/wuq8022600160/database/sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta \
--ref /public/home/wuq8022600160/database/sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta \
--ref /public/home/wuq8022600160/database/sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta \
--ref /public/home/wuq8022600160/database/sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta \
--ref /public/home/wuq8022600160/database/sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta \
--ref /public/home/wuq8022600160/database/sortmerna/data/rRNA_databases/rfam-5s-database-id98.fasta \
--ref /public/home/wuq8022600160/database/sortmerna/data/rRNA_databases/rfam-5.8s-database-id98.fasta \
--workdir 2.assembly-rna/${i} \
--reads 01.cleandata/${i}.R1.raw_val_1.fq.gz \
--reads 01.cleandata/${i}.R2.raw_val_2.fq.gz \
--aligned "2.assembly-rna/${i}N_rRNA" --other "2.assembly-rna/${i}N_clean" --paired_in --fastx --out2
done
```

3.Assembly

MEGAHIT (v1.2.9)

```
mkdir 03.Assembly_megahit
for i in $(cat sampleID); do
megahit -1  2.assembly-rna/${i}N_clean_fwd.fq.gz -2 2.assembly-rna/${i}N_clean_rev.fq.gz -o 03.Assembly_megahit/${i} -t 64 --tmp-dir tmp_megahit23 --presets meta-large
done
```

4.predict the open reading frame (ORF)

Prodigal (v2.6.3) 

```
#!/bin/bash
mkdir 04.prodigal
for i in $(cat sampleID); do
/public/home/wuq8022600160/anaconda3/envs/checkm/bin/prodigal -i 03.Assembly_megahit/${i}/final.contigs.fa -d 04.prodigal/${i}.fasta
done
```

5.obtain nonredundant gene sets

CD-HIT (v4.8.1) 

```
#!/bin/bash
cat 04.prodigal/*.fasta > 04.prodigal/all.fa
/public/home/wuq8022600160/cd-hit/cd-hit-est -c 0.95 -d 400 -T 64 -M 0 -n 8 -aS 0.9 -G 0 -g 1 -i 04.prodigal/all.fa -o 04.prodigal/uni_gene.fa
```

6.taxonomy annotation

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

7.function annotation

emapper-2.1.12 / eggNOG DB version: 5.0.2 /diamond version 2.1.6 / MMseqs2 version 14.7e284/ Compatible novel families DB version: 1.0.1

```
#!/bin/bash
mkdir 05.annotation/eggnog
python3 /public/home/wuq8022600160/anaconda3/envs/eggnog/bin/emapper.py -i 04.prodigal/uni_gene.fasta --itype CDS --translate --cpu 64 --data_dir /public/home/wuq8022600160/xulei/metagenome/eggnog/ --dmnd_db /public/home/wuq8022600160/xulei/metagenome/eggnog/eggnog_proteins.dmnd -o 05.annotation/eggnog/FUNCTION
done
```

8.gene quantitfy

salmon -v # 1.10.1

```
mkdir -p temp/salmon
# Build index, -t sequence, -i index
mkdir 06.quantify/salmon
mkdir 06.quantify/salmon/index 
salmon index -t 04.prodigal/uni_gene.fasta \
-i 06.quantify/salmon/index
for j in $(cat sampleID); do
salmon quant -i 06.quantify/salmon/index -l A --meta \
    -1 01.cleandata/${j}R1.raw_val_1.fq.gz  -2 01.cleandata/${j}R2.raw_val_2.fq.gz \
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



