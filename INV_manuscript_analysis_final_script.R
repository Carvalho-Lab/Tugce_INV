library(dplyr)
library(data.table)
library(ggplot2)
library(rJava)
library(UpSetR)
library(ggrepel)

# Datasets ----------------------------------------------------------------
#Date:15June2024 - 9August2024
#######GNOMAD data######
gnomad_hg38_invs <- read.table(file = "Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/gnomad_v4/merged_all_invs_gnomad_v4_sorted_nochr.bed", header = F, sep = "\t")
final_gnomad_hg38_invs <- gnomad_hg38_invs[, c(1,2,3,4,5)]

colnames(final_gnomad_hg38_invs) <- c("inv_chr", "inv_start", "inv_end", "inv_id", "inv_length")

nrow(final_gnomad_hg38_invs)
final_gnomad_hg38_invs$resource <- "gnomad_v4"

write.table(final_gnomad_hg38_invs, file = "gnomad_hg38_data.bed", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
length(unique(final_gnomad_hg38_invs$inv_id))
summary(final_gnomad_hg38_invs$inv_length)

######DGV data######
hg38_dgv_data <- read.table(file = "Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/dgv/dgv_hg38_variants_2020-02-25.txt", header = T, sep = "\t")
hg38_dgv_invs <- subset(hg38_dgv_data, hg38_dgv_data$variantsubtype == "inversion")
nrow(hg38_dgv_invs)
table(hg38_dgv_invs$reference)

final_dgv_hg38_invs <- hg38_dgv_invs[!hg38_dgv_invs$chr %in% c(""),]
nrow(final_dgv_hg38_invs)

final_dgv_hg38_invs$inv_id <- paste(final_dgv_hg38_invs$chr, final_dgv_hg38_invs$start, final_dgv_hg38_invs$end, sep = "_")
final_dgv_hg38_invs$inv_length <- final_dgv_hg38_invs$end - final_dgv_hg38_invs$start + 1

final_dgv_hg38_invs <- final_dgv_hg38_invs[, c(2,3,4,21,22,7)]

colnames(final_dgv_hg38_invs) <- c("inv_chr", "inv_start", "inv_end", "inv_id", "inv_length", "resource")

nrow(final_dgv_hg38_invs)
length(unique(final_dgv_hg38_invs$inv_id))
table(final_dgv_hg38_invs$inv_chr)
final_dgv_hg38_invs <- final_dgv_hg38_invs[!final_dgv_hg38_invs$inv_chr %in% c("4_GL000008v2_random"),]
nrow(final_dgv_hg38_invs)
length(unique(final_dgv_hg38_invs$inv_id))

final_dgv_hg38_invs <- final_dgv_hg38_invs[!duplicated(final_dgv_hg38_invs$inv_id),]
nrow(final_dgv_hg38_invs)
table(final_dgv_hg38_invs$resource)
write.table(final_dgv_hg38_invs, file = "dgv_hg38_data.bed", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

summary(final_dgv_hg38_invs$inv_length)


######1000G data######
hg38_thog_data <- read.table(file = "Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/1000g/hg38_1KGP_inv.bed", header = F, sep = "\t")
colnames(hg38_thog_data) <- c("inv_chr", "inv_start", "inv_end", "inv_id", "inv_length", "NHet", "NHom")
final_thog_hg38_invs <- select(hg38_thog_data, inv_chr, inv_start, inv_end, inv_id, inv_length)

final_thog_hg38_invs$inv_chr <- gsub("chr","",as.character(final_thog_hg38_invs$inv_chr))
final_thog_hg38_invs$resource <- "Byrska_Bishop_et_al"

write.table(final_thog_hg38_invs, file = "thog_hg38_data.bed", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
nrow(final_thog_hg38_invs)
length(final_thog_hg38_invs$inv_id)

summary(final_thog_hg38_invs$inv_length)

######Ebert et al data######
hg38_ebert_data <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/ebert/ebert_hg38.tsv")
table(hg38_ebert_data$type)
hg38_ebert_data$`#CHROM` <- gsub("chr","",as.character(hg38_ebert_data$`#CHROM`))

final_ebert_hg38_invs <- hg38_ebert_data[, c(2,3,4,1,6,36)]
colnames(final_ebert_hg38_invs) <- c("inv_chr", "inv_start", "inv_end", "inv_id", "inv_length", "resource")

final_ebert_hg38_invs <- final_ebert_hg38_invs[!final_ebert_hg38_invs$inv_chr %in% c("14_GL000225v1_random", "Un_KI270743v1"),]

write.table(final_ebert_hg38_invs, file = "ebert_hg38_data.bed", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

nrow(final_ebert_hg38_invs)
length(unique(final_ebert_hg38_invs$inv_id))

summary(final_ebert_hg38_invs$inv_length)

######porubsky et al data######
hg38_porubsky_data <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/porubsky/porubsky_invs.csv")

final_porubsky_hg38_invs <- hg38_porubsky_data[, c(1,2,3,5,4,9)]
colnames(final_porubsky_hg38_invs ) <- c("inv_chr", "inv_start", "inv_end","inv_id", "inv_length","resource")
final_porubsky_hg38_invs$inv_chr <- gsub("chr","",as.character(final_porubsky_hg38_invs$inv_chr))
write.table(final_porubsky_hg38_invs, file = "porubsky_hg38_data.bed", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
nrow(final_porubsky_hg38_invs)

summary(final_porubsky_hg38_invs$inv_length)

#UCSC and Ensemble biomart
gencode_v46_knowncanonical <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/gencode_knownCanonical")
nrow(gencode_v46_knowncanonical)

ensembl_biomart <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/mart_export.txt")
nrow(ensembl_biomart)

#merge gencode and ensemble biomart
gencode_v46_ensembl_biomart_merged <- merge(gencode_v46_knowncanonical, ensembl_biomart, by.x = "transcript", by.y = "Transcript stable ID version", all.x = TRUE)
nrow(gencode_v46_ensembl_biomart_merged)
length(unique(gencode_v46_ensembl_biomart_merged$`Gene stable ID`))

######OMIM data######
library(readxl)
omim <- readLines("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/omim/genemap2.txt")
omim_hg38_data <- read.delim("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/omim/genemap2.txt", skip = 3, header = TRUE, nrows=length(omim) - 80)
nrow(omim_hg38_data)
length(omim_hg38_data$MIM.Number) 

final_omim_hg38_data <- select(omim_hg38_data, X..Chromosome, Genomic.Position.Start, Genomic.Position.End, MIM.Number, Approved.Gene.Symbol, Gene.Locus.And.Other.Related.Symbols, Gene.Name, Entrez.Gene.ID, Ensembl.Gene.ID, Phenotypes, Comments)
colnames(final_omim_hg38_data) <- c("chr", "start", "end", "Gene_stable_ID", "approved_gene_symbol", "other_gene_symbol", "gene_name", "entrez_gene_id", "Ensembl_gene_ID", "Phenotypes", "Comments")


final_omim_hg38_data$Inheritance_pattern <- ""

final_omim_hg38_data$Inheritance_pattern[grep("Autosomal dominant", final_omim_hg38_data$Phenotypes)] <- "Autosomal dominant"
final_omim_hg38_data$Inheritance_pattern[grep("Autosomal recessive", final_omim_hg38_data$Phenotypes)] <- paste0(final_omim_hg38_data$Inheritance_pattern[grep("Autosomal recessive", final_omim_hg38_data$Phenotypes)], ", ", "Autosomal recessive")
final_omim_hg38_data$Inheritance_pattern[grep("Digenic recessive", final_omim_hg38_data$Phenotypes)] <- "Digenic recessive"
final_omim_hg38_data$Inheritance_pattern[grep("X-linked recessive", final_omim_hg38_data$Phenotypes)] <- "X-linked recessive"
final_omim_hg38_data$Inheritance_pattern[grep("X-linked dominant", final_omim_hg38_data$Phenotypes)] <- "X-linked dominant"
final_omim_hg38_data$Inheritance_pattern[grep("Y-linked", final_omim_hg38_data$Phenotypes)] <- "Y-linked"
final_omim_hg38_data$Inheritance_pattern[which(final_omim_hg38_data$Inheritance_pattern  == ", Autosomal recessive")] = "Autosomal recessive"
final_omim_hg38_data$Inheritance_pattern[which(final_omim_hg38_data$Inheritance_pattern  == "")] = "Other"
table(final_omim_hg38_data$Inheritance_pattern)

final_omim_hg38_data$chr <- gsub("chr","",as.character(final_omim_hg38_data$chr))

hg38_omim_wo_phenotype <- subset(final_omim_hg38_data, final_omim_hg38_data$Phenotypes == "")
nrow(hg38_omim_wo_phenotype)
length(unique(hg38_omim_wo_phenotype$MIM_number))
 
hg38_omim_with_phenotype <- final_omim_hg38_data[-which(final_omim_hg38_data$Phenotypes == ""), ]
nrow(hg38_omim_with_phenotype)
length(unique(hg38_omim_with_phenotype$Ensembl_gene_ID))

#### merging gencode data with omim
merged_final_gencode_v46_ensembl_biomart_omim <- merge(gencode_v46_ensembl_biomart_merged, final_omim_hg38_data, by.x = "Gene stable ID", by.y = "Ensembl_gene_ID", all.x = TRUE)

sorted_merged_final_gencode_v46_ensembl_biomart_omim <- merged_final_gencode_v46_ensembl_biomart_omim[, c(3,4,5,1,10,7,2,8,9,12,13,14,15,17,18,19,20,21,24,26,28)]

nrow(sorted_merged_final_gencode_v46_ensembl_biomart_omim)

unique_merged_final_gencode_v46_ensembl_biomart_omim <- unique(sorted_merged_final_gencode_v46_ensembl_biomart_omim)
nrow(unique_merged_final_gencode_v46_ensembl_biomart_omim)

length(unique(unique_merged_final_gencode_v46_ensembl_biomart_omim$`Gene stable ID`))

# orphadata ---------------------------------------------------------------
all_orpha_genes <- read_excel("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/orphadata/orphanet_disease_genes.xlsx")
length(unique(all_orpha_genes$Gene_Symbol))

orpha_genes <- all_orpha_genes[,1]
orpha_genes <- unique(orpha_genes)
nrow(orpha_genes)
length(orpha_genes$Gene_Symbol)
orpha_genes$gene_class <- "Orphanet_gene"

#add orphanet annotation
merged_orphanet_omim_gencode_ensembl_hg38_data <- merge(unique_merged_final_gencode_v46_ensembl_biomart_omim, orpha_genes, by.x = "Gene name", by.y = "Gene_Symbol", all.x = TRUE)

merged_orphanet_omim_gencode_ensembl_hg38_data$gene_class[is.na(merged_orphanet_omim_gencode_ensembl_hg38_data$gene_class)] <- "Not_Orphanet_gene"

merged_orphanet_omim_gencode_ensembl_hg38_data$Phenotypes[merged_orphanet_omim_gencode_ensembl_hg38_data$Phenotypes == ""] <- "No_OMIM_phenotype"

merged_orphanet_omim_gencode_ensembl_hg38_data$chromStart <- merged_orphanet_omim_gencode_ensembl_hg38_data$chromStart + 1

merged_orphanet_omim_gencode_ensembl_hg38_data$chrom <- gsub("chr","",as.character(merged_orphanet_omim_gencode_ensembl_hg38_data$chrom))

#keep the regular chromosomes
table(merged_orphanet_omim_gencode_ensembl_hg38_data$chrom)
chromosomes <- c("1","2", "3","4","5","6","7","8","9","10","11","12","13","14","15","16",
                 "17","18","19","20","21","22","X","Y")
chr_final_gencode_orpha_omim_ensembl_biomart_merged <- merged_orphanet_omim_gencode_ensembl_hg38_data[which(merged_orphanet_omim_gencode_ensembl_hg38_data$chrom %in% chromosomes),]

#keep the protein-coding ones
ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged <- subset(chr_final_gencode_orpha_omim_ensembl_biomart_merged, chr_final_gencode_orpha_omim_ensembl_biomart_merged$`Transcript type` == "protein_coding")
nrow(ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged)


sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged <- ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged[, c(2,3,4,5,6,7,8,9,10,11,12,13,14,20,21,22)]
nrow(sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged)
length(unique(sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged$`Gene stable ID`))

duplicated(sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged$`Gene stable ID`)

unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged <- unique(sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged)

nrow(unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged)

duplicated(unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged$`Gene stable ID`)

#remove the duplicated rows
unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged %>% group_by(`Gene stable ID`) %>% filter(n() > 1)
filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged <- unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged[-c(4003,15173), ]

nrow(filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged)

write.table(filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged, file = "merged_orphanet_omim_gencodeV46_ensembl_hg38_data.bed", row.names = FALSE, quote = FALSE, col.names = F, sep = "\t")

####

hg38_genes_wo_omim_phenotype <- subset(filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged, filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged$Phenotypes == "No_OMIM_phenotype")
nrow(hg38_genes_wo_omim_phenotype)
length(unique(hg38_genes_wo_omim_phenotype$`Gene stable ID`))

hg38_genes_with_NA_phenotype <- filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged[which(is.na(filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged$Phenotypes)), ]
nrow(hg38_genes_with_NA_phenotype)
length(unique(hg38_genes_with_NA_phenotype$`Gene stable ID`))

hg38_genes_with_omim_phenotype <- filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged[-which(filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged$Phenotypes == "No_OMIM_phenotype" | is.na(filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged$Phenotypes)), ]
nrow(hg38_genes_with_omim_phenotype)
length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))

hg38_orphanet_genes <- subset(filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged, filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged$gene_class == "Orphanet_gene")
length(unique(hg38_orphanet_genes$`Gene stable ID`))

hg38_not_orphanet_genes <- subset(filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged, filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged$gene_class == "Not_Orphanet_gene")
length(unique(hg38_not_orphanet_genes$`Gene stable ID`))

#check the difference between omim data
setdiff(hg38_omim_with_phenotype$Ensembl_gene_ID, hg38_genes_with_omim_phenotype$`Gene stable ID`)
setdiff(hg38_omim_wo_phenotype$Ensembl_gene_ID, hg38_genes_wo_omim_phenotype$`Gene stable ID`)

######GENE INTERSECTIONS######
###bedtools commands
##sort the files
#sort -k1,1V -k2,2n ~/gnomad_hg38_data.bed > sorted_gnomad_hg38_data.bed
#sort -k1,1V -k2,2n ~/dgv_hg38_data.bed > sorted_dgv_hg38_data.bed
#sort -k1,1V -k2,2n ~/thog_hg38_data.bed > sorted_thog_hg38_data.bed
#sort -k1,1V -k2,2n ~/ebert_hg38_data.bed > sorted_ebert_hg38_data.bed
#sort -k1,1V -k2,2n ~/porubsky_hg38_data.bed > sorted_porubsky_hg38_data.bed
#sort -k1,1V -k2,2n ~/merged_orphanet_omim_gencodev46_ensembl_hg38_data.bed > sorted_merged_orphanet_omim_gencodev46_ensembl_hg38_data.bed
##intersect inversions with protein-coding genes
#bedtools intersect -a ~/gnomad_hg38_data.bed -b ~/sorted_merged_orphanet_omim_gencodev46_ensembl_hg38_data.bed -e -sorted -wo -g ~/hg38.chrom.sizes > gnomad_genes_hg38.bed
#bedtools intersect -a ~/sorted_dgv_hg38_data.bed -b ~/sorted_merged_orphanet_omim_gencodev46_ensembl_hg38_data.bed -e -sorted -wo -g ~/hg38.chrom.sizes > dgv_genes_hg38.bed
#bedtools intersect -a ~/sorted_thog_hg38_data.bed -b ~/sorted_merged_orphanet_omim_gencodev46_ensembl_hg38_data.bed -e -sorted -wo -g ~/hg38.chrom.sizes > thog_genes_hg38.bed
#bedtools intersect -a ~/sorted_ebert_hg38_data.bed -b ~/sorted_merged_orphanet_omim_gencodev46_ensembl_hg38_data.bed -e -sorted -wo -g ~/hg38.chrom.sizes > ebert_genes_hg38.bed
#bedtools intersect -a ~/sorted_porubsky_hg38_data.bed -b ~/sorted_merged_orphanet_omim_gencodev46_ensembl_hg38_data.bed -e -sorted -wo -g ~/hg38.chrom.sizes > porubsky_genes_hg38.bed

#gnomad
gnomad_genes <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/invs_genes_intersections/gnomad_genes_hg38.bed")
colnames(gnomad_genes) <- c("inv_chr", "inv_start", "inv_end", "inv_id","inv_length", "resource", "gene_chr", "gene_start", "gene_end", "Gene_stable_ID", "Gene_stable_ID_version", "Transcript_ID", "Transcript_stable_ID", "Protein_ID", "Protein_stable_ID", "Ensenbl_canonical", "Transcript_type", "Refseq_match_transcript", "approved_gene_symbol", "Phenotypes", "Inheritance_pattern", "gene_class", "length_of_intersected_region")
nrow(gnomad_genes)

length(unique(gnomad_genes$inv_id))

gnomad_genes$dataset <- "gnomad"
length(unique(gnomad_genes$inv_id))

gnomad_genes$gene_length <- gnomad_genes$gene_end - gnomad_genes$gene_start + 1


length(unique(gnomad_genes$inv_id))/nrow(final_gnomad_hg38_invs)*100
length(unique(gnomad_genes$Gene_stable_ID))/nrow(filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged)*100

inh_pat_gnomad<- select(gnomad_genes, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_gnomad <- inh_pat_gnomad %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_gnomad$Inheritance_pattern)

gnomad_genes$intersection_category <- ""

for(i in 1:nrow(gnomad_genes)) {       
  if(gnomad_genes$gene_start[i] > gnomad_genes$inv_start[i] & gnomad_genes$gene_end[i] < gnomad_genes$inv_end[i]) {
    gnomad_genes$intersection_category[i] <- 1
    
  }else if(gnomad_genes$gene_start[i] > gnomad_genes$inv_start[i] & gnomad_genes$gene_end[i] > gnomad_genes$inv_end[i]){
    gnomad_genes$intersection_category[i] <- 2 
    
  }else if(gnomad_genes$gene_start[i] < gnomad_genes$inv_start[i] & gnomad_genes$gene_end[i] < gnomad_genes$inv_end[i]){
    gnomad_genes$intersection_category[i] <- 2 
    
  }else if(gnomad_genes$gene_start[i] < gnomad_genes$inv_start[i] & gnomad_genes$gene_end[i] > gnomad_genes$inv_end[i]){
    gnomad_genes$intersection_category[i] <- 3
  }
  
}


gnomad_1 <- subset(gnomad_genes, gnomad_genes$intersection_category == "1")
length(unique(gnomad_1$inv_id))

invs_in_gnomad_1 <- c(unique(gnomad_1$inv_id))

nrow(gnomad_1)/nrow(gnomad_genes)*100
inh_pat_gnomad_1 <- select(gnomad_1, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_gnomad_1 <- inh_pat_gnomad_1 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_gnomad_1$inheritance_pattern)/nrow(gnomad_1)
write.table(gnomad_1, file = "gnomad_1_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

gnomad_2 <- subset(gnomad_genes, gnomad_genes$intersection_category == "2")
nrow(gnomad_2)/nrow(gnomad_genes)*100

invs_in_gnomad_2 <- c(unique(gnomad_2$inv_id))
length(invs_in_gnomad_2)

genes_in_gnomad_2 <- c(unique(gnomad_2$Gene_stable_ID))
length(genes_in_gnomad_2)

inh_pat_gnomad_2 <- select(gnomad_2, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_gnomad_2 <- inh_pat_gnomad_2 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_gnomad_2$Inheritance_pattern)/nrow(gnomad_2)
write.table(gnomad_2, file = "gnomad_2_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

gnomad_3 <- subset(gnomad_genes, gnomad_genes$intersection_category == "3")
nrow(gnomad_3)/nrow(gnomad_genes)*100

invs_in_gnomad_3 <- c(unique(gnomad_3$inv_id))
length(invs_in_gnomad_3)

inh_pat_gnomad_3 <- select(gnomad_3, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_gnomad_3 <- inh_pat_gnomad_3 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_gnomad_3$Inheritance_pattern)/nrow(gnomad_3)
write.table(gnomad_3, file = "gnomad_3_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

#dgv
dgv_genes <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/invs_genes_intersections/dgv_genes_hg38.bed")
colnames(dgv_genes) <- c("inv_chr", "inv_start", "inv_end", "inv_id","inv_length", "resource", "gene_chr", "gene_start", "gene_end", "Gene_stable_ID", "Gene_stable_ID_version", "Transcript_ID", "Transcript_stable_ID", "Protein_ID", "Protein_stable_ID", "Ensenbl_canonical", "Transcript_type", "Refseq_match_transcript", "approved_gene_symbol", "Phenotypes", "Inheritance_pattern", "gene_class", "length_of_intersected_region")
nrow(dgv_genes)
dgv_genes$dataset <- "dgv"
length(unique(dgv_genes$inv_id))

dgv_genes$gene_length <- dgv_genes$gene_end - dgv_genes$gene_start + 1

length(unique(dgv_genes$inv_id))/nrow(final_dgv_hg38_invs)*100

length(unique(dgv_genes$Gene_stable_ID))/nrow(filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged)*100

inh_pat_dgv <- select(dgv_genes, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_dgv <- inh_pat_dgv %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_dgv$Inheritance_pattern)

dgv_genes$intersection_category <- ""

for(i in 1:nrow(dgv_genes)) {       
  if(dgv_genes$gene_start[i] > dgv_genes$inv_start[i] & dgv_genes$gene_end[i] < dgv_genes$inv_end[i]) {
    dgv_genes$intersection_category[i] <- 1
    
  }else if(dgv_genes$gene_start[i] > dgv_genes$inv_start[i] & dgv_genes$gene_end[i] > dgv_genes$inv_end[i]){
    dgv_genes$intersection_category[i] <- 2 
    
  }else if(dgv_genes$gene_start[i] < dgv_genes$inv_start[i] & dgv_genes$gene_end[i] < dgv_genes$inv_end[i]){
    dgv_genes$intersection_category[i] <- 2 
    
  }else if(dgv_genes$gene_start[i] < dgv_genes$inv_start[i] & dgv_genes$gene_end[i] > dgv_genes$inv_end[i]){
    dgv_genes$intersection_category[i] <- 3
  }
  
}

dgv_1 <- subset(dgv_genes, dgv_genes$intersection_category == "1")
nrow(dgv_1)/nrow(dgv_genes)*100
inh_pat_dgv_1 <- select(dgv_1, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_dgv_1 <- inh_pat_dgv_1 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_dgv_1$Inheritance_pattern)/nrow(dgv_1)
write.table(dgv_1, file = "dgv_1_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

dgv_2 <- subset(dgv_genes, dgv_genes$intersection_category == "2")
nrow(dgv_2)/nrow(dgv_genes)*100
inh_pat_dgv_2 <- select(dgv_2, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_dgv_2 <- inh_pat_dgv_2 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_dgv_2$Inheritance_pattern)/nrow(dgv_2)
write.table(dgv_2, file = "dgv_2_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

dgv_3 <- subset(dgv_genes, dgv_genes$intersection_category == "3")
nrow(dgv_3)/nrow(dgv_genes)*100
inh_pat_dgv_3 <- select(dgv_3, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_dgv_3 <- inh_pat_dgv_3 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_dgv_3$Inheritance_pattern)/nrow(dgv_3)
write.table(dgv_3, file = "dgv_3_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

#thog
thog_genes <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/invs_genes_intersections/thog_genes_hg38.bed")
colnames(thog_genes) <- c("inv_chr", "inv_start", "inv_end", "inv_id","inv_length", "resource", "gene_chr", "gene_start", "gene_end", "Gene_stable_ID", "Gene_stable_ID_version", "Transcript_ID", "Transcript_stable_ID", "Protein_ID", "Protein_stable_ID", "Ensenbl_canonical", "Transcript_type", "Refseq_match_transcript", "approved_gene_symbol", "Phenotypes", "Inheritance_pattern", "gene_class", "length_of_intersected_region")
nrow(thog_genes)
thog_genes$dataset <- "thog"
length(unique(thog_genes$inv_id))

thog_genes$gene_length <- thog_genes$gene_end - thog_genes$gene_start + 1

length(unique(thog_genes$inv_id))/nrow(final_thog_hg38_invs)*100

length(unique(thog_genes$Gene_stable_ID))/nrow(filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged)*100

inh_pat_thog<- select(thog_genes, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_thog <- inh_pat_thog %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_thog$Inheritance_pattern)

thog_genes$intersection_category <- ""

for(i in 1:nrow(thog_genes)) {       
  if(thog_genes$gene_start[i] > thog_genes$inv_start[i] & thog_genes$gene_end[i] < thog_genes$inv_end[i]) {
    thog_genes$intersection_category[i] <- 1
    
  }else if(thog_genes$gene_start[i] > thog_genes$inv_start[i] & thog_genes$gene_end[i] > thog_genes$inv_end[i]){
    thog_genes$intersection_category[i] <- 2 
    
  }else if(thog_genes$gene_start[i] < thog_genes$inv_start[i] & thog_genes$gene_end[i] < thog_genes$inv_end[i]){
    thog_genes$intersection_category[i] <- 2 
    
  }else if(thog_genes$gene_start[i] < thog_genes$inv_start[i] & thog_genes$gene_end[i] > thog_genes$inv_end[i]){
    thog_genes$intersection_category[i] <- 3
  }
  
}
thog_1 <- subset(thog_genes, thog_genes$intersection_category == "1")
nrow(thog_1)/nrow(thog_genes)*100
inh_pat_thog_1 <- select(thog_1, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_thog_1 <- inh_pat_thog_1 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_thog_1$Inheritance_pattern)/nrow(thog_1)
write.table(thog_1, file = "thog_1_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

thog_2 <- subset(thog_genes, thog_genes$intersection_category == "2")
nrow(thog_2)/nrow(thog_genes)*100
inh_pat_thog_2 <- select(thog_2, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_thog_2 <- inh_pat_thog_2 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_thog_2$Inheritance_pattern)/nrow(thog_2)
write.table(thog_2, file = "thog_2_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

thog_3 <- subset(thog_genes, thog_genes$intersection_category == "3")
nrow(thog_3)/nrow(thog_genes)*100
inh_pat_thog_3 <- select(thog_3, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_thog_3 <- inh_pat_thog_3 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_thog_3$Inheritance_pattern)/nrow(thog_3)
write.table(thog_3, file = "thog_3_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

#ebert
ebert_genes <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/invs_genes_intersections/ebert_genes_hg38.bed")
colnames(ebert_genes) <- c("inv_chr", "inv_start", "inv_end", "inv_id","inv_length", "resource", "gene_chr", "gene_start", "gene_end", "Gene_stable_ID", "Gene_stable_ID_version", "Transcript_ID", "Transcript_stable_ID", "Protein_ID", "Protein_stable_ID", "Ensenbl_canonical", "Transcript_type", "Refseq_match_transcript", "approved_gene_symbol", "Phenotypes", "Inheritance_pattern", "gene_class", "length_of_intersected_region")
nrow(ebert_genes)
ebert_genes$dataset <- "Ebert"
length(unique(ebert_genes$inv_id))

ebert_genes$gene_length <- ebert_genes$gene_end - ebert_genes$gene_start + 1

length(unique(ebert_genes$inv_id))/nrow(final_ebert_hg38_invs)*100

length(unique(ebert_genes$Gene_stable_ID))

inh_pat_ebert<- select(ebert_genes, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_ebert <- inh_pat_ebert %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_ebert$Inheritance_pattern)

ebert_genes$intersection_category <- ""

for(i in 1:nrow(ebert_genes)) {       
  if(ebert_genes$gene_start[i] > ebert_genes$inv_start[i] & ebert_genes$gene_end[i] < ebert_genes$inv_end[i]) {
    ebert_genes$intersection_category[i] <- 1
    
  }else if(ebert_genes$gene_start[i] > ebert_genes$inv_start[i] & ebert_genes$gene_end[i] > ebert_genes$inv_end[i]){
    ebert_genes$intersection_category[i] <- 2 
    
  }else if(ebert_genes$gene_start[i] < ebert_genes$inv_start[i] & ebert_genes$gene_end[i] < ebert_genes$inv_end[i]){
    ebert_genes$intersection_category[i] <- 2 
    
  }else if(ebert_genes$gene_start[i] < ebert_genes$inv_start[i] & ebert_genes$gene_end[i] > ebert_genes$inv_end[i]){
    ebert_genes$intersection_category[i] <- 3
  }
  
}
ebert_1 <- subset(ebert_genes, ebert_genes$intersection_category == "1")
nrow(ebert_1)/nrow(ebert_genes)*100
inh_pat_ebert_1 <- select(ebert_1, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_ebert_1 <- inh_pat_ebert_1 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_ebert_1$Inheritance_pattern)/nrow(ebert_1)
write.table(ebert_1, file = "ebert_1_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

ebert_2 <- subset(ebert_genes, ebert_genes$intersection_category == "2")
nrow(ebert_2)/nrow(ebert_genes)*100
inh_pat_ebert_2 <- select(ebert_2, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_ebert_2 <- inh_pat_ebert_2 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_ebert_2$Inheritance_pattern)/nrow(ebert_2)
write.table(ebert_2, file = "ebert_2_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

ebert_3 <- subset(ebert_genes, ebert_genes$intersection_category == "3")
nrow(ebert_3)/nrow(ebert_genes)*100
inh_pat_ebert_3 <- select(ebert_3, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_ebert_3 <- inh_pat_ebert_3 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_ebert_3$Inheritance_pattern)/nrow(ebert_3)
write.table(ebert_3, file = "ebert_3_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

#porubsky
porubsky_genes <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/invs_genes_intersections/porubsky_genes_hg38.bed")
colnames(porubsky_genes) <- c("inv_chr", "inv_start", "inv_end", "inv_id","inv_length", "resource", "gene_chr", "gene_start", "gene_end", "Gene_stable_ID", "Gene_stable_ID_version", "Transcript_ID", "Transcript_stable_ID", "Protein_ID", "Protein_stable_ID", "Ensenbl_canonical", "Transcript_type", "Refseq_match_transcript", "approved_gene_symbol", "Phenotypes", "Inheritance_pattern", "gene_class", "length_of_intersected_region")
nrow(porubsky_genes)
porubsky_genes$dataset <- "Porubsky"
length(unique(porubsky_genes$inv_id))

porubsky_genes$gene_length <- porubsky_genes$gene_end - porubsky_genes$gene_start + 1

length(unique(porubsky_genes$inv_id))/nrow(final_porubsky_hg38_invs)*100

length(unique(porubsky_genes$Gene_stable_ID))
length(unique(porubsky_genes$Gene_stable_ID))/length(unique(filtered_unique_sorted_ProteinCoding_orphanet_omim_gencode_v46_ensembl_biomart_merged$`Gene stable ID`))*100

inh_pat_porubsky<- select(porubsky_genes, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_porubsky <- inh_pat_porubsky %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_porubsky$Inheritance_pattern)


porubsky_genes$intersection_category <- ""

for(i in 1:nrow(porubsky_genes)) {       
  if(porubsky_genes$gene_start[i] > porubsky_genes$inv_start[i] & porubsky_genes$gene_end[i] < porubsky_genes$inv_end[i]) {
    porubsky_genes$intersection_category[i] <- 1
    
  }else if(porubsky_genes$gene_start[i] > porubsky_genes$inv_start[i] & porubsky_genes$gene_end[i] > porubsky_genes$inv_end[i]){
    porubsky_genes$intersection_category[i] <- 2 
    
  }else if(porubsky_genes$gene_start[i] < porubsky_genes$inv_start[i] & porubsky_genes$gene_end[i] < porubsky_genes$inv_end[i]){
    porubsky_genes$intersection_category[i] <- 2 
    
  }else if(porubsky_genes$gene_start[i] < porubsky_genes$inv_start[i] & porubsky_genes$gene_end[i] > porubsky_genes$inv_end[i]){
    porubsky_genes$intersection_category[i] <- 3
  }
  
}

porubsky_1 <- subset(porubsky_genes, porubsky_genes$intersection_category == "1")
nrow(porubsky_1)/nrow(porubsky_genes)*100
inh_pat_porubsky_1 <- select(porubsky_1, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_porubsky_1 <- inh_pat_porubsky_1 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_porubsky_1$inheritance_pattern)/nrow(porubsky_1)
write.table(porubsky_1, file = "porubsky_1_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

porubsky_2 <- subset(porubsky_genes, porubsky_genes$intersection_category == "2")
nrow(porubsky_2)/nrow(porubsky_genes)*100
inh_pat_porubsky_2 <- select(porubsky_2, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_porubsky_2 <- inh_pat_porubsky_2 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_porubsky_2$inheritance_pattern)/nrow(porubsky_2)
write.table(porubsky_2, file = "porubsky_2_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

porubsky_3 <- subset(porubsky_genes, porubsky_genes$intersection_category == "3")
nrow(porubsky_3)/nrow(porubsky_genes)*100
inh_pat_porubsky_3 <- select(porubsky_3, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_porubsky_3 <- inh_pat_porubsky_3 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_porubsky_3$Inheritance_pattern)/nrow(porubsky_3)
write.table(porubsky_3, file = "porubsky_3_genes.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

# Main Figures -----------------------------------------------------------------

#Figure 3
intersection_category_1_df <- data_frame(dataset = c("gnomAD", "DGV", "1KGP", "Ebert et al.", "Porubsky et al."),
                                         percentage_of_intersection = c(round(nrow(gnomad_1)/nrow(gnomad_genes)*100,1), round(nrow(dgv_1)/nrow(dgv_genes)*100,1), round(nrow(thog_1)/nrow(thog_genes)*100,1), round(nrow(ebert_1)/nrow(ebert_genes)*100,1), round(nrow(porubsky_1)/nrow(porubsky_genes)*100,1)))

intersection_category_2_df <- data_frame(dataset = c("gnomAD", "DGV", "1KGP", "Ebert et al.", "Porubsky et al."),
                                         percentage_of_intersection = c(round(nrow(gnomad_2)/nrow(gnomad_genes)*100,1), round(nrow(dgv_2)/nrow(dgv_genes)*100,1), round(nrow(thog_2)/nrow(thog_genes)*100,1), round(nrow(ebert_2)/nrow(ebert_genes)*100,1), round(nrow(porubsky_2)/nrow(porubsky_genes)*100,1)))

intersection_category_3_df <- data_frame(dataset = c("gnomAD", "DGV", "1KGP", "Ebert et al.", "Porubsky et al."),
                                         percentage_of_intersection = c(round(nrow(gnomad_3)/nrow(gnomad_genes)*100,1), round(nrow(dgv_3)/nrow(dgv_genes)*100,1), round(nrow(thog_3)/nrow(thog_genes)*100,1), round(nrow(ebert_3)/nrow(ebert_genes)*100,1), round(nrow(porubsky_3)/nrow(porubsky_genes)*100,1)))

ggplot(data=intersection_category_1_df, aes(x=dataset, y=percentage_of_intersection, fill = dataset)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=percentage_of_intersection), vjust=0.5, color="black", size=16) +
  coord_cartesian(ylim = c(0, 100)) + 
  labs(x ="", y = "% of Overlaps") +
  theme_update(text = element_text(size=32)) +
  theme_classic() +
  theme(legend.title=element_text(size=32),
        axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 34),
        axis.title = element_text(size = 34),
        legend.position = "none") +
  theme(axis.text.x = element_text(colour = c("black")),
        axis.text.y = element_text(colour = c("black")))


ggplot(data=intersection_category_2_df, aes(x=dataset, y=percentage_of_intersection, fill = dataset)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=percentage_of_intersection), vjust=1, color="black", size=16) +
  coord_cartesian(ylim = c(0, 100)) + 
  labs(x ="", y = "% of Overlaps") +
  theme_update(text = element_text(size=32)) +
  theme_classic() +
  theme(legend.title=element_text(size=32),
        axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 34),
        axis.title = element_text(size = 34),
        legend.position = "none") +
  theme(axis.text.x = element_text(colour = c("black")),
        axis.text.y = element_text(colour = c("black")))

ggplot(data=intersection_category_3_df, aes(x=dataset, y=percentage_of_intersection, fill = dataset)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=percentage_of_intersection), vjust=1, color="black", size=16) +
  coord_cartesian(ylim = c(0, 100)) + 
  labs(x ="", y = "% of Overlaps") +
  theme_update(text = element_text(size=32)) +
  theme_classic() +
  theme(legend.title=element_text(size=32),
        axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32),
        axis.title = element_text(size = 32),
        legend.position = "none") +
  theme(axis.text.x = element_text(colour = c("black")),
        axis.text.y = element_text(colour = c("black")))


#Figure 4
intersection_category_1_all_datasets <- rbind(gnomad_1, dgv_1, thog_1, ebert_1, porubsky_1)
length(unique(intersection_category_1_all_datasets$inv_id))
intersection_category_1_all_datasets$inv_name <- paste(intersection_category_1_all_datasets$inv_chr, intersection_category_1_all_datasets$inv_start, intersection_category_1_all_datasets$inv_end, sep = "_")
length(unique(intersection_category_1_all_datasets$inv_name))

intersection_category_2_all_datasets <- rbind(gnomad_2, dgv_2, thog_2, ebert_2, porubsky_2)
intersection_category_2_all_datasets$inv_name <- paste(intersection_category_2_all_datasets$inv_chr, intersection_category_2_all_datasets$inv_start, intersection_category_2_all_datasets$inv_end, sep = "_")
length(unique(intersection_category_2_all_datasets$inv_name))

length(unique(intersection_category_2_all_datasets$Gene_stable_ID))
write.table(intersection_category_2_all_datasets, file = "intersection_category_2_all_datasets.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")


intersection_category_3_all_datasets <- rbind(gnomad_3, dgv_3, thog_3, ebert_3, porubsky_3)
intersection_category_3_all_datasets$inv_name <- paste(intersection_category_3_all_datasets$inv_chr, intersection_category_3_all_datasets$inv_start, intersection_category_3_all_datasets$inv_end, sep = "_")
length(unique(intersection_category_3_all_datasets$inv_name))

length(unique(intersection_category_3_all_datasets$Gene_stable_ID))
write.table(intersection_category_3_all_datasets, file = "intersection_category_3_all_datasets.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")


#category2
intersection_category_2_all_datasets$OMIM_Phenotype <- c()
intersection_category_2_all_datasets$OMIM_Phenotype[intersection_category_2_all_datasets$Phenotype == "No_OMIM_phenotype"] <- "No_OMIM_phenotype"
intersection_category_2_all_datasets$OMIM_Phenotype[is.na(intersection_category_2_all_datasets$OMIM_Phenotype)] <- "NA"
intersection_category_2_all_datasets$OMIM_Phenotype[-which(intersection_category_2_all_datasets$Phenotypes == "No_OMIM_phenotype" | is.na(intersection_category_2_all_datasets$Phenotypes))] <- "OMIM_phenotype"

unique_invs_intersection_category_2_all_datasets <- intersection_category_2_all_datasets %>% distinct(inv_chr, inv_start, inv_end)
nrow(unique_invs_intersection_category_2_all_datasets)

unique_genes_intersection_category_2_all_datasets <- intersection_category_2_all_datasets %>% distinct(Gene_stable_ID, OMIM_Phenotype, gene_class)
nrow(unique_genes_intersection_category_2_all_datasets)

table(unique_genes_intersection_category_2_all_datasets$OMIM_Phenotype)
table(unique_genes_intersection_category_2_all_datasets$OMIM_Phenotype)/nrow(unique_genes_intersection_category_2_all_datasets)*100

178/length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))*100
490/length(unique(hg38_genes_wo_omim_phenotype$`Gene stable ID`))*100
162/length(unique(hg38_genes_with_NA_phenotype$`Gene stable ID`))*100

#category3
intersection_category_3_all_datasets$OMIM_Phenotype <- c()
intersection_category_3_all_datasets$OMIM_Phenotype[intersection_category_3_all_datasets$Phenotype == "No_OMIM_phenotype"] <- "No_OMIM_phenotype"
intersection_category_3_all_datasets$OMIM_Phenotype[is.na(intersection_category_3_all_datasets$OMIM_Phenotype)] <- "NA"
intersection_category_3_all_datasets$OMIM_Phenotype[-which(intersection_category_3_all_datasets$Phenotypes == "No_OMIM_phenotype" | is.na(intersection_category_3_all_datasets$Phenotypes))] <- "OMIM_phenotype"

unique_invs_intersection_category_3_all_datasets <- intersection_category_3_all_datasets %>% distinct(inv_chr, inv_start, inv_end)
nrow(unique_invs_intersection_category_3_all_datasets)

unique_genes_intersection_category_3_all_datasets <- intersection_category_3_all_datasets %>% distinct(Gene_stable_ID, OMIM_Phenotype, gene_class)
nrow(unique_genes_intersection_category_3_all_datasets)

table(intersection_category_3_all_datasets$OMIM_Phenotype)/nrow(intersection_category_3_all_datasets)*100

table(unique_genes_intersection_category_3_all_datasets$OMIM_Phenotype)
332/length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))*100
599/length(unique(hg38_genes_wo_omim_phenotype$`Gene stable ID`))*100
99/length(unique(hg38_genes_with_NA_phenotype$`Gene stable ID`))*100

invs_omim_plot_df <- data.frame(Category = c("category 2", "category 2", "category 2", "category 3", "category 3", "category 3"),
                                Group = c("Genes in OMIM without any phenotype", "Genes in OMIM with at least one phenotype", "Genes not presented in OMIM", "Genes in OMIM without any phenotype", "Genes in OMIM with at least one phenotype", "Genes not presented in OMIM"),
                                Percentage = c(4.3, 3.6, 4.6, 5.3, 6.8,2.9)
)


ggplot(data=invs_omim_plot_df, aes(x=Category, y=Percentage, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette = "Dark2")+
  labs(x ="Category", y = "% of genes") +
  theme_update(text = element_text(size=32)) +
  theme_classic() +
  theme(legend.title=element_text(size=32),
        legend.text=element_text(size=32),
        axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32),
        axis.title = element_text(size = 32),
        ) +
  theme(axis.text.x = element_text(colour = c("black")),
        axis.text.y = element_text(colour = c("black")))+
  theme(legend.position="bottom", legend.direction="vertical") 


# pie charts of inheritance patterns
#cat2
intersection_category_2_all_datasets_omim_pheno <- intersection_category_2_all_datasets[-which(intersection_category_2_all_datasets$Phenotypes == "No_OMIM_phenotype" | is.na(intersection_category_2_all_datasets$Phenotypes)), ]
length(unique(intersection_category_2_all_datasets_omim_pheno$Gene_stable_ID))
table(intersection_category_2_all_datasets_omim_pheno$Inheritance_pattern) / nrow(intersection_category_2_all_datasets_omim_pheno) * 100

all_invs_category_2_omim_inh_pattern <- data.frame(
  Inheritance_pattern = c("AD","AR", "AR/AD", "Digenic recessive", "Other", "X-linked dominant", "X-linked recessive", "Y-linked"),
  value = c(29.4, 41.5, 19.8, 0.8, 4.4, 0.8, 3.2, 0)
)
df_cat2 <- all_invs_category_2_omim_inh_pattern %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(all_invs_category_2_omim_inh_pattern, aes(x = "" , y = value, fill = Inheritance_pattern)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set2") +
  geom_label_repel(data = df_cat2,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 6, nudge_x = 1, show.legend = F) +
  guides(fill = guide_legend(title = "Inheritance pattern")) +
  theme_void()

#cat3
intersection_category_3_all_datasets_omim_pheno <- intersection_category_3_all_datasets[-which(intersection_category_3_all_datasets$Phenotypes == "No_OMIM_phenotype" | is.na(intersection_category_3_all_datasets$Phenotypes)), ]
length(unique(intersection_category_3_all_datasets_omim_pheno$Gene_stable_ID))
round(table(intersection_category_3_all_datasets_omim_pheno$Inheritance_pattern) / nrow(intersection_category_3_all_datasets_omim_pheno) * 100, 1)

all_invs_category_3 <- data.frame(
  Inheritance_pattern = c("AD","AR", "AR/AD", "Digenic recessive", "Other", "X-linked dominant", "X-linked recessive", "Y-linked"),
  value = c(28.9, 47.9, 10.5, 0.7, 7.7, 1.4, 2.7, 0.2)
)
df_cat3 <- all_invs_category_3 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))
ggplot(all_invs_category_3, aes(x = "" , y = value, fill = Inheritance_pattern)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set2") +
  geom_label_repel(data = df_cat3,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 6, nudge_x = 1, show.legend = F) +
  guides(fill = guide_legend(title = "Inheritance pattern")) +
  theme_void()

#Supplementary Figures ---------------------------------------------------

#Supplementary Figure 1
#barplot of transcript types
table(merged_orphanet_omim_gencode_ensembl_hg38_data$`Transcript type`)
transcript_type_barplot <- ggplot(merged_orphanet_omim_gencode_ensembl_hg38_data, aes(x=`Transcript type`, fill = `Transcript type`))+
  geom_bar(stat="count", width=0.7) +
  theme_classic() 
transcript_type_barplot + 
  theme(axis.text.x = element_text(angle = 90, size = 22),
        axis.title = element_text(size = 22),
        axis.text.y = element_text(size = 22)) +
  geom_text(stat='count', aes(label=..count..), vjust=-1, size = 7) +
  theme(axis.text.x = element_text(colour = c("black")),
        axis.text.y = element_text(colour = c("black"))) +
  theme(legend.position = "none") 


#Supplementary Figure 4
# gnomaAD frequency analyses ---------------------------------------------------------------
gnomad_hg38_invs$SVLEN <- as.character(sapply(gnomad_hg38_invs$V6, function(x) strsplit(strsplit(x, "SVLEN=")[[1]][2], ";")[[1]][1]))

gnomad_hg38_invs$AF <- as.character(sapply(gnomad_hg38_invs$V6, function(x) strsplit(strsplit(x, "AF=")[[1]][2], ";")[[1]][1]))
gnomad_hg38_invs$AF <- as.numeric(gnomad_hg38_invs$AF)

gnomad_hg38_invs$N_HET <- as.character(sapply(gnomad_hg38_invs$V6, function(x) strsplit(strsplit(x, "N_HET=")[[1]][2], ";")[[1]][1]))
gnomad_hg38_invs$N_HOMALT <- as.character(sapply(gnomad_hg38_invs$V6, function(x) strsplit(strsplit(x, "N_HOMALT=")[[1]][2], ";")[[1]][1]))
gnomad_hg38_invs$FREQ_HET <- as.character(sapply(gnomad_hg38_invs$V6, function(x) strsplit(strsplit(x, "FREQ_HET=")[[1]][2], ";")[[1]][1]))
gnomad_hg38_invs$FREQ_HOMALT <- as.character(sapply(gnomad_hg38_invs$V6, function(x) strsplit(strsplit(x, "FREQ_HOMALT=")[[1]][2], ";")[[1]][1]))
gnomad_hg38_invs$FREQ_HOMALT <- as.numeric(gnomad_hg38_invs$FREQ_HOMALT)
gnomad_hg38_invs$FREQ_HET <- as.numeric(gnomad_hg38_invs$FREQ_HET)


gnomad_hg38_invs$var_type <- ""

gnomad_hg38_invs$var_type[which(gnomad_hg38_invs$AF < 0.05)] <- "rare(AF<0.05)"
gnomad_hg38_invs$var_type[which(gnomad_hg38_invs$AF >= 0.05)] <- "common"

table(gnomad_hg38_invs$var_type)

final_gnomad_hg38_detailed_invs <- gnomad_hg38_invs[, c(1,2,3,4,7,8,9,10,11,12,13)]

write.table(final_gnomad_hg38_detailed_invs, file = "gnomad_hg38_invs_detailed.bed", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

rare_gnomad_hg38_detailed_invs <- subset(final_gnomad_hg38_detailed_invs, final_gnomad_hg38_detailed_invs$var_type == "rare(AF<0.05)")
write.table(rare_gnomad_hg38_detailed_invs, file = "rare_gnomad_hg38_detailed_invs.bed", row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")

#frequency pie chart
table(final_gnomad_hg38_detailed_invs$var_type)/nrow(final_gnomad_hg38_detailed_invs)*100

freq_table <- data.frame(
  group = c("common", "rare(<0.05)"),
  value = c(1.1, 98.9)
)

bp <- ggplot(freq_table, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")
bp
pie <- bp + coord_polar("y", start=0)
library(scales)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_blank(),
  )

pie + scale_fill_brewer(palette="Set2") + blank_theme +
  theme(axis.text.x=element_blank())+
  geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), 
                label = percent(value/100)), size=9) +
  theme(legend.title=element_text(size=22),
        legend.text=element_text(size=22))


######INTERSECTIONS######
###bedtools commands
##sort
#sort -k1,1V -k2,2n ~/gnomad_hg38_invs_detailed.bed > sorted_gnomad_hg38_invs_detailed.bed
#bedtools intersect -a ~/sorted_gnomad_hg38_invs_detailed.bed -b ~/sorted_merged_orphanet_omim_gencodev46_ensembl_hg38_data.bed -e -sorted -wo -g ~/hg38.chrom.sizes > gnomad_detailed_genes_hg38.bed

#gnomad 
gnomad_detailed_genes <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/gnomad_v4_freq/gnomad_detailed_genes_hg38.bed")
colnames(gnomad_detailed_genes) <- c("inv_chr", "inv_start", "inv_end", "inv_id","inv_length", "AF", "N_HET", "N_HOMALT", "FREQ_HET", "FREQ_HOMALT", "var_type", "gene_chr", "gene_start", "gene_end", "Gene_stable_ID", "Gene_stable_ID_version", "Transcript_ID", "Transcript_stable_ID", "Protein_ID", "Protein_stable_ID", "Ensenbl_canonical", "Transcript_type", "Refseq_match_transcript", "approved_gene_symbol", "Phenotypes", "Inheritance_pattern", "gene_class", "length_of_intersected_region")    
nrow(gnomad_detailed_genes)

gnomad_detailed_genes$gene_length <- gnomad_detailed_genes$gene_end - gnomad_detailed_genes$gene_start + 1

gnomad_detailed_genes$intersection_category <- ""

for(i in 1:nrow(gnomad_detailed_genes)) {       
  if(gnomad_detailed_genes$gene_start[i] > gnomad_detailed_genes$inv_start[i] & gnomad_detailed_genes$gene_end[i] < gnomad_detailed_genes$inv_end[i]) {
    gnomad_detailed_genes$intersection_category[i] <- 1
    
  }else if(gnomad_detailed_genes$gene_start[i] > gnomad_detailed_genes$inv_start[i] & gnomad_detailed_genes$gene_end[i] > gnomad_detailed_genes$inv_end[i]){
    gnomad_detailed_genes$intersection_category[i] <- 2 
    
  }else if(gnomad_detailed_genes$gene_start[i] < gnomad_detailed_genes$inv_start[i] & gnomad_detailed_genes$gene_end[i] < gnomad_detailed_genes$inv_end[i]){
    gnomad_detailed_genes$intersection_category[i] <- 2 
    
  }else if(gnomad_detailed_genes$gene_start[i] < gnomad_detailed_genes$inv_start[i] & gnomad_detailed_genes$gene_end[i] > gnomad_detailed_genes$inv_end[i]){
    gnomad_detailed_genes$intersection_category[i] <- 3
  }
  
}

gnomad_detailed_genes$OMIM_Phenotype <- c()
gnomad_detailed_genes$OMIM_Phenotype[gnomad_detailed_genes$Phenotype == "No_OMIM_phenotype"] <- "No_OMIM_phenotype"
gnomad_detailed_genes$OMIM_Phenotype[is.na(gnomad_detailed_genes$OMIM_Phenotype)] <- "NA"
gnomad_detailed_genes$OMIM_Phenotype[-which(gnomad_detailed_genes$Phenotypes == "No_OMIM_phenotype" | is.na(gnomad_detailed_genes$Phenotypes))] <- "OMIM_phenotype"

#RARE inversions
rare_gnomad_detailed_genes <- subset(gnomad_detailed_genes, gnomad_detailed_genes$var_type == "rare(AF<0.05)")
nrow(rare_gnomad_detailed_genes)
table(rare_gnomad_detailed_genes$intersection_category)

#filter out category 1
rare_gnomad_detailed_genes_category2_3 <- rare_gnomad_detailed_genes[-which(rare_gnomad_detailed_genes$intersection_category == "1")]
nrow(rare_gnomad_detailed_genes_category2_3)

length(rare_gnomad_detailed_genes_category2_3$inv_id)
unique_genes_in_rare_gnomad_detailed_genes <- rare_gnomad_detailed_genes_category2_3 %>% distinct(Gene_stable_ID, OMIM_Phenotype, gene_class)
nrow(unique_genes_in_rare_gnomad_detailed_genes)

table(unique_genes_in_rare_gnomad_detailed_genes$OMIM_Phenotype)
table(unique_genes_in_rare_gnomad_detailed_genes$OMIM_Phenotype, unique_genes_in_rare_gnomad_detailed_genes$gene_class)

write.table(rare_gnomad_detailed_genes_category2_3, file = "rare_gnomad_detailed_genes_category2_3.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

rare_gnomad_plot_df <- data.frame(
                                Group = c("Genes in OMIM without any phenotype", "Genes in OMIM with at least one phenotype", "Genes not presented in OMIM"),
                                Count = c(521, 247, 105)
)

library(RColorBrewer)
ggplot(data=rare_gnomad_plot_df, aes(x=Group, y=Count, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#FFAB91","#FF8A65","#FF6F00"))+
  theme_classic() +
  theme(legend.title=element_text(size=24),
        legend.text=element_text(size=24),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = c("black")),
        axis.title = element_text(size = 24),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.position="bottom", legend.direction="vertical") 
  
rare_gnomad_detailed_genes_category2_3_omim_pheno <- subset(rare_gnomad_detailed_genes_category2_3,rare_gnomad_detailed_genes_category2_3$OMIM_Phenotype == "OMIM_phenotype")
length(unique(rare_gnomad_detailed_genes_category2_3_omim_pheno$Gene_stable_ID))
length(unique(rare_gnomad_detailed_genes_category2_3_omim_pheno$inv_id))
a <- rare_gnomad_detailed_genes_category2_3_omim_pheno %>% distinct(Gene_stable_ID, Inheritance_pattern)

table(a$Inheritance_pattern)/nrow(a) * 100

rare_gnomad_detailed_wo_hom_freq_genes_category2_3_omim_pheno <- subset(rare_gnomad_detailed_genes_category2_3_omim_pheno,rare_gnomad_detailed_genes_category2_3_omim_pheno$N_HOMALT == "0")
length(unique(rare_gnomad_detailed_wo_hom_freq_genes_category2_3_omim_pheno$inv_id))

rare_gnomad_detailed_wo_hom_freq_genes_category2_3_omim_AR_pheno <- subset(rare_gnomad_detailed_wo_hom_freq_genes_category2_3_omim_pheno,rare_gnomad_detailed_wo_hom_freq_genes_category2_3_omim_pheno$Inheritance_pattern == "Autosomal recessive")
length(unique(rare_gnomad_detailed_wo_hom_freq_genes_category2_3_omim_AR_pheno$Gene_stable_ID))

write.table(rare_gnomad_detailed_wo_hom_freq_genes_category2_3_omim_AR_pheno, file = "rare_gnomad_detailed_wo_hom_freq_genes_category2_3_omim_AR_pheno.txt", row.names = FALSE, quote = FALSE, col.names = T, sep = "\t")

#COMMON inversions
common_gnomad_detailed_genes <- subset(gnomad_detailed_genes, gnomad_detailed_genes$var_type == "common")

nrow(common_gnomad_detailed_genes)
table(common_gnomad_detailed_genes$intersection_category)

unique_genes_in_common_gnomad_detailed_genes <- common_gnomad_detailed_genes %>% distinct(Gene_stable_ID, OMIM_Phenotype, gene_class)
nrow(unique_genes_in_common_gnomad_detailed_genes)

table(unique_genes_in_common_gnomad_detailed_genes$OMIM_Phenotype)
table(unique_genes_in_common_gnomad_detailed_genes$OMIM_Phenotype, unique_genes_in_common_gnomad_detailed_genes$gene_class)

common_gnomad_plot_df <- data.frame(
  Group = c("Genes in OMIM without any phenotype", "Genes in OMIM with at least one phenotype", "Genes not presented in OMIM"),
  Count = c(6, 7, 5)
)
ggplot(data=common_gnomad_plot_df, aes(x=Group, y=Count, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c("#B2DFDB","#4DB6AC", "#00838F"))+
  scale_y_continuous(limits = c(0,500)) +
  theme_classic() +
  theme(legend.title=element_text(size=24),
        legend.text = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24, colour = c("black")),
        axis.title = element_text(size = 24),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(legend.position="bottom", legend.direction="vertical") 


#Supplementary Figure 5
lplot_gnomad <- final_gnomad_hg38_invs
lplot_gnomad$dataset <- "gnomAD"

lplot_dgv <- final_dgv_hg38_invs
lplot_dgv$dataset <- "DGV"

lplot_thog <-  final_thog_hg38_invs
lplot_thog$dataset <- "1KGP"

lplot_ebert <- final_ebert_hg38_invs
lplot_ebert$dataset <- "Ebert et al."

lplot_porubsky <- final_porubsky_hg38_invs
lplot_porubsky$dataset <- "Porubsky et al."

all_invs <- rbind(lplot_gnomad, lplot_dgv, lplot_thog, lplot_ebert, lplot_porubsky)

length_plot <- ggplot(all_invs, aes(x= dataset, y=log2(inv_length), fill= dataset)) +
  geom_violin(trim=FALSE, show.legend = FALSE)

length_plot + labs(x ="", y = "Inversion length (bp,log2 transformed)", size = 18) +
  theme_update(text = element_text(size=24)) +
  theme_classic() +
  theme(legend.title=element_text(size=24),
        axis.text.x = element_text(size = 24, colour = c("black")),
        axis.text.y = element_text(size = 24, colour = c("black")),
        axis.title = element_text(size = 24),
        ) 

#Supplementary Figure 6
venn_gnomad <- final_gnomad_hg38_invs
venn_gnomad$inv_name <- paste(venn_gnomad$inv_chr, venn_gnomad$inv_start, venn_gnomad$inv_end, sep = "_")
length(unique(venn_gnomad$inv_name))

venn_dgv <- final_dgv_hg38_invs
venn_dgv$inv_name <- paste(venn_dgv$inv_chr, venn_dgv$inv_start, venn_dgv$inv_end, sep = "_")
length(unique(venn_dgv$inv_name))

venn_thog <- final_thog_hg38_invs
venn_thog$inv_name <- paste(venn_thog$inv_chr, venn_thog$inv_start, venn_thog$inv_end, sep = "_")
length(unique(venn_thog$inv_name))

venn_ebert <- final_ebert_hg38_invs
venn_ebert$inv_name <- paste(venn_ebert$inv_chr, venn_ebert$inv_start, venn_ebert$inv_end, sep = "_")
length(unique(venn_ebert$inv_name))

venn_porubsky <- final_porubsky_hg38_invs
venn_porubsky$inv_name <- paste(venn_porubsky$inv_chr, venn_porubsky$inv_start, venn_porubsky$inv_end, sep = "_")
length(unique(venn_porubsky$inv_name))


listInput <- list(gnomAD = venn_gnomad$inv_name, DGV = venn_dgv$inv_name, '1KGP' = venn_thog$inv_name, 'Ebert et al.' = venn_ebert$inv_name, 'Porubsky et al.'= venn_porubsky$inv_name)

upset(fromList(listInput), 
      order.by = "freq",
      point.size=5,
      nsets = 5, 
      mainbar.y.label = "Number of inversions", sets.x.label = "", 
      text.scale = c(2, 2, 2, 2, 2, 1.5),
      sets.bar.color=c("#A3A500","#00B0F6","#F8766D", "#00BF7D", "#E76BF3"),
)

#Supplementary Figure 7
######Inversion intersections
#bedtools commands
#bedtools intersect -a ~/sorted_gnomad_hg38_data.bed -b ~/sorted_dgv_hg38_data.bed ~/sorted_thog_hg38_data.bed ~/sorted_ebert_hg38_data.bed ~/sorted_porubsky_hg38_data.bed -f 0.5 -F 0.5 -e -sorted -names dgv thog ebert porubsky -wo > gnomad_others.bed
#bedtools intersect -a ~/sorted_dgv_hg38_data.bed -b ~/sorted_gnomad_hg38_data.bed ~/sorted_thog_hg38_data.bed ~/sorted_ebert_hg38_data.bed ~/sorted_porubsky_hg38_data.bed -f 0.5 -F 0.5 -e -sorted -names gnomad thog ebert porubsky  -wo > dgv_others.bed
#bedtools intersect -a ~/sorted_thog_hg38_data.bed -b ~/sorted_gnomad_hg38_data.bed ~/sorted_dgv_hg38_data.bed ~/sorted_ebert_hg38_data.bed ~/sorted_porubsky_hg38_data.bed -f 0.5 -F 0.5 -e -sorted -names gnomad dgv ebert porubsky -wo > thog_others.bed
#bedtools intersect -a ~/sorted_ebert_hg38_data.bed -b ~/sorted_gnomad_hg38_data.bed ~/sorted_dgv_hg38_data.bed ~/sorted_thog_hg38_data.bed ~/sorted_porubsky_hg38_data.bed -f 0.5 -F 0.5 -e -sorted -names gnomad dgv thog porubsky -wo > ebert_others.bed
#bedtools intersect -a ~/sorted_porubsky_hg38_data.bed -b ~/sorted_gnomad_hg38_data.bed ~/sorted_dgv_hg38_data.bed ~/sorted_thog_hg38_data.bed ~/sorted_ebert_hg38_data.bed -f 0.5 -F 0.5 -e -sorted -names gnomad dgv thog ebert -wo > porubsky_others.bed
#
#gnomad
gnomad_others <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/inversion_intersections/gnomad_others.bed")
colnames(gnomad_others) <- c("inv_chr", "inv_start", "inv_end", "inv_id","inv_length", "resource", "database", "others_chr", "others_start", "others_end", "others_inv_id","others_inv_length", "resource","length_of_intersected_region")


gnomad_dgv <- subset(gnomad_others, gnomad_others$database == "dgv")
length(unique(gnomad_dgv$inv_id))/nrow(final_gnomad_hg38_invs)*100

gnomad_thog <- subset(gnomad_others, gnomad_others$database == "thog")
length(unique(gnomad_thog$inv_id))/nrow(final_gnomad_hg38_invs)*100

gnomad_ebert <- subset(gnomad_others, gnomad_others$database == "ebert")
length(unique(gnomad_ebert$inv_id))/nrow(final_gnomad_hg38_invs)*100

gnomad_porubsky <- subset(gnomad_others, gnomad_others$database == "porubsky")
length(unique(gnomad_porubsky$inv_id))/nrow(final_gnomad_hg38_invs)*100

gnomad_plot_df <- data.frame(Database = c("DGV", "1KGP", "Ebert et al.", "Porubsky et al."), Percentage = c(49.4, 27.3, 11.2, 11.1))

gnomad_overlapped_invs_plot <- ggplot(data=gnomad_plot_df, aes(x = Database, y=Percentage, fill = Database)) +
  geom_bar(stat = "identity", width=0.5) +
  xlab("") + ylab("% Inversions") +
  theme_classic()+
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  geom_text(aes(label=Percentage), vjust=1.6, color="white", size=12)+
  theme(legend.title=element_text(size=24),
        axis.text.x = element_text(size = 24, colour = c("black")),
        axis.text.y = element_text(size = 24, colour = c("black")),
        axis.title = element_text(size = 24),
        legend.position = "none") +
  scale_fill_manual(values=c("#00A5FF","#0D47A1","#90CAF9","#00B8E5"))
gnomad_overlapped_invs_plot + coord_cartesian(ylim=c(0,80))

#dgv
dgv_others <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/inversion_intersections/dgv_others.bed")
colnames(dgv_others) <- c("inv_chr", "inv_start", "inv_end", "inv_id","inv_length", "resource", "database", "others_chr", "others_start", "others_end", "others_inv_id","others_inv_length", "resource","length_of_intersected_region")

dgv_gnomad <- subset(dgv_others, dgv_others$database == "gnomad")
length(unique(dgv_gnomad$inv_id))/nrow(final_dgv_hg38_invs)*100

dgv_thog <- subset(dgv_others, dgv_others$database == "thog")
length(unique(dgv_thog$inv_id))/nrow(final_dgv_hg38_invs)*100

dgv_ebert <- subset(dgv_others, dgv_others$database == "ebert")
length(unique(dgv_ebert$inv_id))/nrow(final_dgv_hg38_invs)*100

dgv_porubsky <- subset(dgv_others, dgv_others$database == "porubsky")
length(unique(dgv_porubsky$inv_id))/nrow(final_dgv_hg38_invs)*100

dgv_plot_df <- data.frame(Database = c("gnomAD", "1KGP", "Ebert et al.", "Porubsky et al."), Percentage = c(77.2, 20.7, 24.1, 24.3))

dgv_overlapped_invs_plot <- ggplot(data=dgv_plot_df, aes(x = Database, y=Percentage, fill = Database)) +
  geom_bar(stat = "identity", width=0.5) +
  xlab("") + ylab("% Inversions") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  geom_text(aes(label=Percentage), vjust=1.6, color="white", size=12)+
  theme(legend.title=element_text(size=24),
        axis.text.x = element_text(size = 24, colour = c("black")),
        axis.text.y = element_text(size = 24, colour = c("black")),
        axis.title = element_text(size = 24),
        legend.position = "none") +
  scale_fill_manual(values=c("#9E9D24","#827717","#FFE082","#BB9D00"))
dgv_overlapped_invs_plot + coord_cartesian(ylim=c(0,80))

#thog
thog_others <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/inversion_intersections/thog_others.bed")
colnames(thog_others) <- c("inv_chr", "inv_start", "inv_end", "inv_id","inv_length", "resource", "database", "others_chr", "others_start", "others_end", "others_inv_id","others_inv_length", "resource","length_of_intersected_region")

thog_gnomad <- subset(thog_others, thog_others$database == "gnomad")
length(unique(thog_gnomad$inv_id))/nrow(final_thog_hg38_invs)*100

thog_dgv <- subset(thog_others, thog_others$database == "dgv")
length(unique(thog_dgv$inv_id))/nrow(final_thog_hg38_invs)*100

thog_ebert <- subset(thog_others, thog_others$database == "ebert")
length(unique(thog_ebert$inv_id))/nrow(final_thog_hg38_invs)*100

thog_porubsky <- subset(thog_others, thog_others$database == "porubsky")
length(unique(thog_porubsky$inv_id))/nrow(final_thog_hg38_invs)*100


thog_plot_df <- data.frame(Database = c("gnomAD", "DGV", "Ebert et al.", "Porubsky et al."), Percentage = c(78.3, 43.5, 12.5, 15.2))

thog_overlapped_invs_plot <- ggplot(data=thog_plot_df, aes(x = Database, y=Percentage, fill = Database)) +
  geom_bar(stat = "identity", width=0.5) +
  xlab("") + ylab("% Inversions") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  geom_text(aes(label=Percentage), vjust=1.6, color="white", size=12)+
  theme(legend.title=element_text(size=24),
        axis.text.x = element_text(size = 24, colour = c("black")),
        axis.text.y = element_text(size = 24, colour = c("black")),
        axis.title = element_text(size = 24),
        legend.position = "none") +
  scale_fill_manual(values=c("#FF7043","#FF5742","#EF9A9A","#F8766D"))
thog_overlapped_invs_plot + coord_cartesian(ylim=c(0,80))

#ebert
ebert_others <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/inversion_intersections/ebert_others.bed")
colnames(ebert_others) <- c("inv_chr", "inv_start", "inv_end", "inv_id","inv_length", "resource", "database", "others_chr", "others_start", "others_end", "others_inv_id","others_inv_length", "resource","length_of_intersected_region")

ebert_gnomad <- subset(ebert_others, ebert_others$database == "gnomad")
length(unique(ebert_gnomad$inv_id))/nrow(final_ebert_hg38_invs)*100

ebert_dgv<- subset(ebert_others, ebert_others$database == "dgv")
length(unique(ebert_dgv$inv_id))/nrow(final_ebert_hg38_invs)*100

ebert_thog<- subset(ebert_others, ebert_others$database == "thog")
length(unique(ebert_thog$inv_id))/nrow(final_ebert_hg38_invs)*100

ebert_porubsky<- subset(ebert_others, ebert_others$database == "porubsky")
length(unique(ebert_porubsky$inv_id))/nrow(final_ebert_hg38_invs)*100

ebert_plot_df <- data.frame(Database = c("gnomAD", "DGV", "1KGP", "Porubsky et al."), Percentage = c(64.0, 62.3, 20.5, 69.1))

ebert_overlapped_invs_plot <- ggplot(data=ebert_plot_df, aes(x = Database, y=Percentage, fill = Database)) +
  geom_bar(stat = "identity", width=0.5) +
  xlab("") + ylab("% Inversions") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  geom_text(aes(label=Percentage), vjust=1.6, color="white", size=12)+
  theme(legend.title=element_text(size=24),
        axis.text.x = element_text(size = 24, colour = c("black")),
        axis.text.y = element_text(size = 24, colour = c("black")),
        axis.title = element_text(size = 24),
        legend.position = "none") +
  scale_fill_manual(values=c("#00BC59","#00695C","#A5D6A7","#00BF7D"))
ebert_overlapped_invs_plot + coord_cartesian(ylim=c(0,80))


#porubsky
porubsky_others <- fread("Z:/Members/clun/manuscripts/tugce_inv/inversion_manuscript_datasets/inversion_intersections/porubsky_others.bed")
colnames(porubsky_others) <- c("inv_chr", "inv_start", "inv_end", "inv_id","inv_length", "resource", "database", "others_chr", "others_start", "others_end", "others_inv_id","others_inv_length", "resource","length_of_intersected_region")

porubsky_gnomad <- subset(porubsky_others, porubsky_others$database == "gnomad")
length(unique(porubsky_gnomad$inv_id))/nrow(final_porubsky_hg38_invs)*100

porubsky_dgv <- subset(porubsky_others, porubsky_others$database == "dgv")
length(unique(porubsky_dgv$inv_id))/nrow(final_porubsky_hg38_invs)*100

porubsky_thog <- subset(porubsky_others, porubsky_others$database == "thog")
length(unique(porubsky_thog$inv_id))/nrow(final_porubsky_hg38_invs)*100

porubsky_ebert <- subset(porubsky_others, porubsky_others$database == "ebert")
length(unique(porubsky_ebert$inv_id))/nrow(final_porubsky_hg38_invs)*100

porubsky_plot_df <- data.frame(Database = c("gnomAD", "DGV", "1KGP", "Ebert et al."), Percentage = c(62.9, 70.4, 21.6, 76.4))

porubsky_overlapped_invs_plot <- ggplot(data=porubsky_plot_df, aes(x = Database, y=Percentage, fill = Database)) +
  geom_bar(stat = "identity", width=0.5) +
  xlab("") + ylab("% Inversions") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  geom_text(aes(label=Percentage), vjust=1.6, color="white", size=12)+
  theme(legend.title=element_text(size=24),
        axis.text.x = element_text(size = 24, colour = c("black")),
        axis.text.y = element_text(size = 24, colour = c("black")),
        axis.title = element_text(size = 24),
        legend.position = "none") +
  scale_fill_manual(values=c("#AB47BC","#8E24AA","#E1BEE7","#DC71FA"))
porubsky_overlapped_invs_plot + coord_cartesian(ylim=c(0,80))

#Supplementary Figure 8
######Percentage of OMIM genes overlapping with inversions
#category 2
#gnomad
gnomad_2_with_phenotype <- gnomad_2[-which(gnomad_2$Phenotypes == "No_OMIM_phenotype" | is.na(gnomad_2$Phenotypes)), ]
length(unique(gnomad_2_with_phenotype$Gene_stable_ID))/length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))*100

#dgv
dgv_2_with_phenotype <- dgv_2[-which(dgv_2$Phenotypes == "No_OMIM_phenotype" | is.na(dgv_2$Phenotypes)), ]
length(unique(dgv_2_with_phenotype$Gene_stable_ID))/length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))*100

#1kgp
thog_2_with_phenotype <- thog_2[-which(thog_2$Phenotypes == "No_OMIM_phenotype" | is.na(thog_2$Phenotypes)), ]
length(unique(thog_2_with_phenotype$Gene_stable_ID))/length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))*100

#ebert
ebert_2_with_phenotype <- ebert_2[-which(ebert_2$Phenotypes == "No_OMIM_phenotype" | is.na(ebert_2$Phenotypes)), ]
length(unique(ebert_2_with_phenotype$Gene_stable_ID))/length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))*100

#porubsky
porubsky_2_with_phenotype <- porubsky_2[-which(porubsky_2$Phenotypes == "No_OMIM_phenotype" | is.na(porubsky_2$Phenotypes)), ]
length(unique(porubsky_2_with_phenotype$Gene_stable_ID))/length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))*100

invs_omim_plot_df_cat2 <- data.frame(Database = c("gnomAD", "DGV", "1KGP", "Ebert et al.", "Porubsky et al."), 
                                     Percentage = c(2.1, 1.6, 0.4, 0.2, 0.1))

ggplot(invs_omim_plot_df_cat2, aes(x = Database, y =Percentage, fill = Database)) +
  geom_bar(stat = "identity", width=0.5) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 5, 1)) +
  geom_text(aes(label=Percentage), vjust=0.5, color="black", size=14) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x ="", y = "% of OMIM phenotype-related genes") +
  theme(legend.title=element_text(size=30),
        axis.text.x = element_text(size = 30, colour = c("black")),
        axis.text.y = element_text(size = 30, colour = c("black")),
        axis.title = element_text(size = 30),
        legend.position = "none") +
  coord_cartesian(ylim=c(0,5)) 

#category 3
#gnomad
gnomad_3_with_phenotype <- gnomad_3[-which(gnomad_3$Phenotypes == "No_OMIM_phenotype" | is.na(gnomad_3$Phenotypes)), ]
length(unique(gnomad_3_with_phenotype$Gene_stable_ID))/length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))*100

#dgv
dgv_3_with_phenotype <- dgv_3[-which(dgv_3$Phenotypes == "No_OMIM_phenotype" | is.na(dgv_3$Phenotypes)), ]
length(unique(dgv_3_with_phenotype$Gene_stable_ID))/length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))*100

#1kgp
thog_3_with_phenotype <- thog_3[-which(thog_3$Phenotypes == "No_OMIM_phenotype" | is.na(thog_3$Phenotypes)), ]
length(unique(thog_3_with_phenotype$Gene_stable_ID))/length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))*100

#ebert
ebert_3_with_phenotype <- ebert_3[-which(ebert_3$Phenotypes == "No_OMIM_phenotype" | is.na(ebert_3$Phenotypes)), ]
length(unique(ebert_3_with_phenotype$Gene_stable_ID))/length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))*100

#porubsky
porubsky_3_with_phenotype <- porubsky_3[-which(porubsky_3$Phenotypes == "No_OMIM_phenotype" | is.na(porubsky_3$Phenotypes)), ]
length(unique(porubsky_3_with_phenotype$Gene_stable_ID))/length(unique(hg38_genes_with_omim_phenotype$`Gene stable ID`))*100

invs_omim_plot_df_cat3 <- data.frame(Database = c("gnomAD", "DGV", "1KGP", "Ebert et al.", "Porubsky et al."), 
                                     Percentage = c(3.2, 3.7, 1.7, 0.5, 0.2))

ggplot(invs_omim_plot_df_cat3, aes(x = Database, y =Percentage, fill = Database)) +
  geom_bar(stat = "identity", width=0.5) +
  theme_classic() +
  scale_y_continuous(breaks = seq(0, 5, 1)) +
  geom_text(aes(label=Percentage), vjust=1.6, color="black", size=14) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x ="", y = "% of OMIM phenotype-related genes") +
  theme(legend.title=element_text(size=30),
        axis.text.x = element_text(size = 30, colour = c("black")),
        axis.text.y = element_text(size = 30, colour = c("black")),
        axis.title = element_text(size = 30),
        legend.position = "none") +
  coord_cartesian(ylim=c(0,5)) 


#Pie charts
#gnomad
inh_pat_gnomad_2 <- select(gnomad_2_with_phenotype, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_gnomad_2 <- inh_pat_gnomad_2 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_gnomad_2$Inheritance_pattern)/nrow(gnomad_2)
round(table(inh_pat_gnomad_2$Inheritance_pattern)/length(unique(gnomad_2_with_phenotype$Gene_stable_ID))*100,1)



df_gnomad_2 <- data.frame(
  Inheritance_pattern = c("AD","AR", "AR/AD", "Digenic recessive", "Other", "X-linked dominant", "X-linked recessive", "Y-linked"),
  value = c(24, 46.2, 26, 1, 2.9, 0, 0, 0)
)

df2_gnomad <- df_gnomad_2 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df_gnomad_2, aes(x = "" , y = value, fill = Inheritance_pattern)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Blues") +
  geom_label_repel(data = df2_gnomad,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 8, nudge_x = 1, show.legend = F) +
  guides(fill = guide_legend(title = "Inheritance pattern")) +
  theme_void() +
  theme(legend.text=element_text(size=14),
        legend.title = element_text(size = 14))


#dgv
inh_pat_dgv_2 <- select(dgv_2_with_phenotype, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_dgv_2 <- inh_pat_dgv_2 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_dgv_2$Inheritance_pattern)/nrow(dgv_2)
round(table(inh_pat_dgv_2$Inheritance_pattern)/length(unique(dgv_2_with_phenotype$Gene_stable_ID))*100,1)

df_dgv_2 <- data.frame(
  Inheritance_pattern = c("AD","AR", "AR/AD", "Digenic recessive", "Other", "X-linked dominant", "X-linked recessive", "Y-linked"),
  value = c(33.8, 41.6, 11.7, 1.3, 5.2, 2.6, 3.9,0)
)

df2_dgv <- df_dgv_2 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df_dgv_2, aes(x = "" , y = value, fill = Inheritance_pattern)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#F0F4C3","#DCE775", "#CDDC39","#AFB42B","#9E9D24","#BB9D00","#827717","#8D6E63")) +
  geom_label_repel(data = df2_dgv,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 8, nudge_x = 1, show.legend = F) +
  guides(fill = guide_legend(title = "Inheritance pattern")) +
  theme_void() +
  theme(legend.text=element_text(size=14),
        legend.title = element_text(size = 14))


#thog
inh_pat_thog_2 <- select(thog_2_with_phenotype, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_thog_2 <- inh_pat_thog_2 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_thog_2$Inheritance_pattern)/nrow(thog_2)
round(table(inh_pat_thog_2$Inheritance_pattern)/length(unique(thog_2_with_phenotype$Gene_stable_ID))*100,1)

df_thog_2 <- data.frame(
  Inheritance_pattern = c("AD","AR", "AR/AD", "Digenic recessive", "Other", "X-linked dominant", "X-linked recessive","Y-linked"),
  value = c(27.8, 44.4, 11.1, 0, 5.6, 0, 11.1,0)
)

df2_thog <- df_thog_2 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df_thog_2, aes(x = "" , y = value, fill = Inheritance_pattern)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Reds") +
  geom_label_repel(data = df2_thog,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 8, nudge_x = 1, show.legend = F) +
  guides(fill = guide_legend(title = "Inheritance pattern")) +
  theme_void() +
  theme(legend.text=element_text(size=14),
        legend.title = element_text(size = 14))

#ebert
inh_pat_ebert_2 <- select(ebert_2_with_phenotype, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_ebert_2 <- inh_pat_ebert_2 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_ebert_2$Inheritance_pattern)/nrow(ebert_2)
round(table(inh_pat_ebert_2$Inheritance_pattern)/length(unique(ebert_2_with_phenotype$Gene_stable_ID))*100,1)

df_ebert_2 <- data.frame(
  Inheritance_pattern = c("AD","AR", "AR/AD", "Digenic recessive", "Other", "X-linked dominant", "X-linked recessive", "Y-linked"),
  value = c(16.7, 41.7, 16.7, 0, 16.7, 0, 8.3, 0)
)

df2_ebert <- df_ebert_2 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df_ebert_2, aes(x = "" , y = value, fill = Inheritance_pattern)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#C8E6C9","#A5D6A7", "#81C784","#00BC59","#00BF7D", "#26A69A","#00897B","#00695C"))  +
  geom_label_repel(data = df2_ebert,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 8, nudge_x = 1, show.legend = F) +
  guides(fill = guide_legend(title = "Inheritance pattern")) +
  theme_void() +
  theme(legend.text=element_text(size=14),
        legend.title = element_text(size = 14))


#porubsky
inh_pat_porubsky_2 <- select(porubsky_2_with_phenotype, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_porubsky_2 <- inh_pat_porubsky_2 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
table(inh_pat_porubsky_2$Inheritance_pattern)/nrow(porubsky_2)
round(table(inh_pat_porubsky_2$Inheritance_pattern)/length(unique(porubsky_2_with_phenotype$Gene_stable_ID))*100,1)

df_porubsky_2 <- data.frame(
  Inheritance_pattern = c("AD","AR", "AR/AD", "Digenic recessive", "Other", "X-linked dominant", "X-linked recessive", "Y-linked"),
  value = c(16.7, 50, 0, 0, 16.7, 0, 16.7, 0)
)

df2_porubsky <- df_porubsky_2 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df_porubsky_2, aes(x = "" , y = value, fill = Inheritance_pattern)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "RdPu") +
  geom_label_repel(data = df2_porubsky,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 8, nudge_x = 1, show.legend = F) +
  guides(fill = guide_legend(title = "Inheritance pattern")) +
  theme_void() +
  theme(legend.text=element_text(size=14),
        legend.title = element_text(size = 14))


#gnomad3
inh_pat_gnomad_3 <- select(gnomad_3_with_phenotype, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_gnomad_3 <- inh_pat_gnomad_3 %>% 
  distinct(Gene_stable_ID, .keep_all = T)

round(table(inh_pat_gnomad_3$Inheritance_pattern)/length(unique(gnomad_3_with_phenotype$Gene_stable_ID))*100,1)

df_gnomad_3 <- data.frame(
  Inheritance_pattern = c("AD","AR", "AR/AD", "Digenic recessive", "Other", "X-linked dominant", "X-linked recessive", "Y-linked"),
  value = c(31.4, 45.9, 11.3, 0, 8.2, 1.3, 1.9,0)
)

df3_gnomad <- df_gnomad_3 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df_gnomad_3, aes(x = "" , y = value, fill = Inheritance_pattern)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Blues") +
  geom_label_repel(data = df3_gnomad,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 8, nudge_x = 1, show.legend = F) +
  guides(fill = guide_legend(title = "Inheritance pattern")) +
  theme_void() +
  theme(legend.text=element_text(size=14),
        legend.title = element_text(size = 14))


#dgv3
inh_pat_dgv_3 <- select(dgv_3_with_phenotype, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_dgv_3 <- inh_pat_dgv_3 %>% 
  distinct(Gene_stable_ID, .keep_all = T)

round(table(inh_pat_dgv_3$Inheritance_pattern)/length(unique(dgv_3_with_phenotype$Gene_stable_ID))*100,1)

df_dgv_3 <- data.frame(
  Inheritance_pattern = c("AD","AR", "AR/AD", "Digenic recessive", "Other", "X-linked dominant", "X-linked recessive", "Y-linked"),
  value = c(28.3, 49.5, 12.5, 0, 6, 1.1, 2.2,0)
)

df3_dgv <- df_dgv_3 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df_dgv_3, aes(x = "" , y = value, fill = Inheritance_pattern)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#F0F4C3","#DCE775", "#CDDC39","#AFB42B","#9E9D24","#BB9D00","#827717","#8D6E63")) +
  geom_label_repel(data = df3_dgv,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 8, nudge_x = 1, show.legend = F) +
  guides(fill = guide_legend(title = "Inheritance pattern")) +
  theme_void() +
  theme(legend.text=element_text(size=14),
        legend.title = element_text(size = 14))

#thog3
inh_pat_thog_3 <- select(thog_3_with_phenotype, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_thog_3 <- inh_pat_thog_3 %>% 
  distinct(Gene_stable_ID, .keep_all = T)

round(table(inh_pat_thog_3$Inheritance_pattern)/length(unique(thog_3_with_phenotype$Gene_stable_ID))*100,1)

df_thog_3 <- data.frame(
  Inheritance_pattern = c("AD","AR", "AR/AD", "Digenic recessive", "Other", "X-linked dominant", "X-linked recessive", "Y-linked"),
  value = c(32.1, 41.7, 8.3, 1.2, 8.3, 1.2, 6, 1.2)
)

df3_thog <- df_thog_3 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df_thog_3, aes(x = "" , y = value, fill = Inheritance_pattern)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Reds") +
  geom_label_repel(data = df3_thog,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 8, nudge_x = 1, show.legend = F) +
  guides(fill = guide_legend(title = "Inheritance pattern")) +
  theme_void() +
  theme(legend.text=element_text(size=14),
        legend.title = element_text(size = 14))


#ebert
inh_pat_ebert_3 <- select(ebert_3_with_phenotype, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_ebert_3 <- inh_pat_ebert_3 %>% 
  distinct(Gene_stable_ID, .keep_all = T)
round(table(inh_pat_ebert_3$Inheritance_pattern)/length(unique(ebert_3_with_phenotype$Gene_stable_ID))*100,1)

df_ebert_3 <- data.frame(
  Inheritance_pattern = c("AD","AR", "AR/AD", "Digenic recessive", "Other", "X-linked dominant", "X-linked recessive", "Y-linked"),
  value = c(18.2, 40.9, 18.2, 4.5, 13.6, 4.5, 0, 0)
)

df3_ebert <- df_ebert_3 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df_ebert_3, aes(x = "" , y = value, fill = Inheritance_pattern)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values=c("#C8E6C9","#A5D6A7", "#81C784","#00BC59","#00BF7D", "#26A69A","#00897B","#00695C"))  +
  geom_label_repel(data = df3_ebert,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 8, nudge_x = 1, show.legend = F) +
  guides(fill = guide_legend(title = "Inheritance pattern")) +
  theme_void() +
  theme(legend.text=element_text(size=14),
        legend.title = element_text(size = 14))


#porubsky
inh_pat_porubsky_3 <- select(porubsky_3_with_phenotype, c("Gene_stable_ID", "Inheritance_pattern"))
inh_pat_porubsky_3 <- inh_pat_porubsky_3 %>% 
  distinct(Gene_stable_ID, .keep_all = T)

round(table(inh_pat_porubsky_3$Inheritance_pattern)/length(unique(porubsky_3_with_phenotype$Gene_stable_ID))*100,1)

df_porubsky_3 <- data.frame(
  Inheritance_pattern = c("AD","AR", "AR/AD", "Digenic recessive", "Other", "X-linked dominant", "X-linked recessive", "Y-linked"),
  value = c(16.7, 50, 8.3, 0, 16.7, 8.3, 0, 0)
)

df3_porubsky <- df_porubsky_3 %>% 
  mutate(csum = rev(cumsum(rev(value))), 
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

ggplot(df_porubsky_3, aes(x = "" , y = value, fill = Inheritance_pattern)) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "RdPu") +
  geom_label_repel(data = df3_porubsky,
                   aes(y = pos, label = paste0(value, "%")),
                   size = 8, nudge_x = 1, show.legend = F) +
  guides(fill = guide_legend(title = "Inheritance pattern")) +
  theme_void() +
  theme(legend.text=element_text(size=14),
        legend.title = element_text(size = 14))

