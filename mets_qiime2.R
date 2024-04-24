if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R") # current version is 0.99.20

install.packages("ape")
install.packages("Biostrings")
install.packages("biomformat")
install.packages("phyloseq")
install.packages("Hmisc")
install.packages("yaml")

library(ape)
library(Biostrings)
library(biomformat)
library(phyloseq)
library(Hmisc)
library(yaml)
library(tidyr)
library(dplyr)
library(stats)
library(utils)

library(qiime2R)

file <- "/home/wheelerlab3/Data/METS/METS_2022-10-14_gut_microbiome/mets-rep-seqs.qza"

data <- read_qza(file)

data$data[1:5,1:5]


system('mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh')

system('wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-linux-conda.yml')


system('export PATH="/home/jacob/miniconda3/bin:$PATH"')

system('/home/jacob/miniconda3/bin/conda list')


#use wget to download greengenes reference reads and reference taxonomy for mets microbiome classification

wget \
-O 'gg-13-8-99-nb-classifier.qza' \
'https://docs.qiime2.org/jupyterbooks/cancer-microbiome-intervention-tutorial/data/030-tutorial-downstream/020-taxonomy/gg-13-8-99-nb-classifier.qza'


system('/home/jacob/miniconda3/bin/conda env create -n qiime2 --file /home/jacob/qiime2-2022.2-py38-linux-conda.yml')
system('/home/jacob/miniconda3/bin/conda activate qiime2')

/home/wheelerlab3/anaconda2/bin/conda env create -n qiime2 --file /home/jacob/qiime2-2022.2-py38-linux-conda.yml

conda activate qiime2

qiime feature-classifier classify-sklearn \
--i-classifier /home/jacob/gg-13-8-99-nb-classifier.qza \
--i-reads /home/wheelerlab3/Data/METS/METS_2022-10-14_gut_microbiome/mets-rep-seqs.qza \
--o-classification taxonomy.qza


system("qiime feature-classifier classify-consensus-blast --i-query /home/wheelerlab3/Data/METS/METS_2022-10-14_gut_microbiome/mets-rep-seqs.qza --i-reference-reads /var/DB/greengenes/gg_99_reference_seqs.qza --i-reference-taxonomy /var/DB/greengenes/gg_99_reference_taxonomy.qza ")
#--p-perc-identity 0.97
#--o-classification {BLAST-TAXONOMY}.qza \
#--verbose

qiime metadata tabulate --m-input-file /home/jacob/taxonomy.qza --o-visualization /home/jacob/taxonomy.qzv

qiime rescript merge-taxa --i-data /home/camilla/taxonomy.qza /home/camilla/taxonomy_silva.qza --p-mode score --o-merged-data taxonomy_merged.qza

qiime metadata tabulate --m-input-file /home/camilla/taxonomy_silva.qza --o-visualization /home/jacob/taxonomy_silva.qzv

iconv -t UTF-8 -f ISO-8859-1 taxonomy_merged.qza > taxonomy_merged.txt
file -I taxonomy_merged.txt
taxonomy_merged.txt: text/plain; charset=utf-8

qiime metadata tabulate --m-input-file out.txt --o-visualization mapping_skinABiL.qzv

qiime metadata tabulate --m-input-file taxonomy_merged.txt --o-visualization taxonomy_merged.qzv

qiime taxa barplot --i-table /home/wheelerlab3/Data/METS/METS_2022-10-14_gut_microbiome/mets-table.qza --i-taxonomy /home/jacob/taxonomy.qza --m-metadata-file /home/wheelerlab3/Data/METS/METS_2022-10-14_gut_microbiome/mets-metadata.tsv --o-visualization taxa-bar-plots.qzv

tax_out <- read_qza("/home/jacob/taxonomy.qza")

#create a feature table that has taxonomy instead of feature IDs

qiime taxa collapse \
--i-table /home/wheelerlab3/Data/METS/METS_2022-10-14_gut_microbiome/mets-table.qza \
--i-taxonomy taxonomy.qza \
--p-level 7 \
--o-collapsed-table phyla-table.qza 

#convert new frequency to relative frequency

qiime feature-table relative-frequency \
--i-table phyla-table.qza \
--o-relative-frequency-table rel-phyla-table.qza

qiime tools export --input-path rel-phyla-table.qza \
--output-path rel-table

biom convert -i /home/wheelerlab3/Data/METS/METS_2022-10-14_gut_microbiome/extracted-mets-table/41a5cfa7-22c4-4cdf-ac2d-374b073c494d/data/feature-table.biom -o extracted-mets-table.tsv --to-tsv

biom convert -i rel-table/feature-table.biom -o rel-phyla-table.tsv --to-tsv

/home/jacob/mets_qiime.R

biom add-metadata \
--input-fp /home/wheelerlab3/Data/METS/METS_2022-10-14_gut_microbiome/extracted-mets-table/41a5cfa7-22c4-4cdf-ac2d-374b073c494d/data/feature-table.biom \
--observation-metadata-fp /home/camilla/metadata.tsv \
--output-fp biom-with-taxonomy.biom

biom convert \
--input-fp biom-with-taxonomy.biom \
--output-fp biom-with-taxonomy.tsv \
--to-tsv \


--output-observation-metadata \
--header-key taxonomy

colnames(taxonomy.qzv)[which(names(dataframe) == "columnName")] <- "newColumnName"

merge_taxa <- read.table(file = "/home/camilla/METS/merged_biom_taxonomy_taxalastcolumn.tsv", sep = '\t', header = TRUE)


library(dplyr)

merge_taxa %>% 
  rowwise() %>% 
  mutate(
    Zero = sum(c_across(starts_with("X"))==0),
    NonZero = sum(c_across(starts_with("X"))!=0)
  )

merge_taxa$non_zero <- rowSums(merge_taxa != 0)

merge_taxa <- filter(merge_taxa, non_zero > 213)

merge_taxa <- merge_taxa %>% arrange(desc(non_zero))
