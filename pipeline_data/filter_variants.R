library("VariantAnnotation")
library("stringr")
library("dplyr")

##### Filter variants ######

  args = commandArgs(trailingOnly = TRUE)
  vcf_file = args[1]
  ref = "hg38"

  #Load VCF file and convert to dataframe
  vcf <- readVcf(vcf_file, ref)
  vcf_info <- data.frame(info(vcf))

  #Load epilepsy gene list info
  gene_file <- args[2] 
  gene_info <- read.table(file = gene_file, sep = '\t', header = TRUE)
  
  
  vcf_info$CHR <- as.character(seqnames(vcf))
  vcf_info$POS <- as.character(start(vcf))
  vcf_info$REF <- as.character(ref(vcf))
  vcf_info$ALT <- data.frame(alt(vcf))
  vcf_info$ALT <- vcf_info$ALT$value
  vcf_info$Variant <- result <- gsub("\\s+", "", (paste(vcf_info$CHR,":",vcf_info$POS, ":", vcf_info$REF, ">", vcf_info$ALT)))
  

  genotypes <- data.frame(geno(vcf)$GT)
  genotype_summary <- c()
  
  # Iterate variants and add sample names of variants with alt alleles to a new column
  for (i in 1:nrow(genotypes)) {
    sample_genotypes <- c()
    for (sample in colnames(genotypes)) {
      if (!(genotypes[i, sample] %in% c("0/0", "0|0", "./.", ".|."))) {
        sample_genotypes <- c(sample_genotypes, paste0(sample, " (", genotypes[i, sample], ")"))
      }
    }
    # Concatenate sample names and genotypes into string and store in list
    sample_string <- paste(sample_genotypes, collapse = "; ")
    genotype_summary[i] <- sample_string
  }
  
  vcf_info$Genotype_Summary <- genotype_summary
  
  #Replicate variant rows have multiple VEP entries (i.e. variants that impact multiple genes)
  vcf_info <- vcf_info %>% 
    unnest(CSQ) %>% 
    distinct(CSQ, .keep_all = TRUE)
  
  ##VEP CSQ##
  #VEP consequences are split into individual fields
  if ("CSQ" %in% colnames(vcf_info)){
    print('Processing Ensembl VEP CSQ')
    csq = info(header(vcf))['CSQ',]$Description
    csq_col = unlist(strsplit(csq, ':'))[2]
    csq_fields = unlist(strsplit(csq_col, '\\|'))
    csq_data <- data.frame(matrix(nrow = nrow(vcf_info), ncol = length(csq_fields)))
    colnames(csq_data) <- csq_fields
    for (i in 1:nrow(vcf_info)){
      transcript <- unlist(vcf_info[i,]$CSQ)
      fields <- unlist(strsplit(transcript, split = '\\|'))
      for (col in 1:length(csq_fields)){
        csq_data[i,col] <- fields[col]}
    }
    vcf_info <- cbind(csq_data, vcf_info)}
  
  ## gnomAD ##
  if ("AC_gnomad3.0" %in% colnames(vcf_info)){
    print('Processing AC_gnomad3.0')
    for (i in 1:length(vcf_info$AC_gnomad3.0)){
      entry = vcf_info$AC_gnomad3.0[i]
      if (length(unlist(entry)) > 1) {
        max = 0
        for (j in 1:length(unlist(entry))){
          if (!is.na(unlist(entry)[j])){
            count = as.numeric(unlist(entry)[j])
            if (count > max){
              max = count}
          }
          entry = max}
      } else if (is.na(entry)){
        entry = 0
      }
      vcf_info$AC_gnomad3.0[i] = as.numeric(entry)
    }
    vcf_info$AC_gnomad3.0 <- as.numeric(vcf_info$AC_gnomad3.0)
  }
  
  if ("Homozygous_gnomad3.0" %in% colnames(vcf_info)){
    print('Processing Homozygous_gnomad3.0')
    for (i in 1:length(vcf_info$Homozygous_gnomad3.0)){
      entry = vcf_info$Homozygous_gnomad3.0[i]
      if (length(unlist(entry)) > 1) {
        max = 0
        for (j in 1:length(unlist(entry))){
          if (!is.na(unlist(entry)[j])){
            count = as.numeric(unlist(entry)[j])
            if (count > max){
              max = count}
          }
          entry = max}
      } else if (is.na(entry)){
        entry = 0
      }
      vcf_info$Homozygous_gnomad3.0[i] = as.numeric(entry)
    }
    vcf_info$Homozygous_gnomad3.0 <- as.numeric(vcf_info$Homozygous_gnomad3.0)
  }
  
  #Match VCF genes with Epilepsy Gene List genes
  matches <- match(vcf_info$Gene, gene_info$Ensemble_ID)
  merged_matches <- ifelse(is.na(matches), match(vcf_info$SYMBOL, gene_info$Gene), matches)
  
  gene_summary <- gene_info[merged_matches,]
  
  # create a new column with the matching Phenotypes
  vcf_info_plus <- cbind(vcf_info, gene_summary)
  
  #Add intron distance column
  vcf_info_plus$intron_dist <- as.numeric(str_extract(vcf_info_plus$HGVSc, "(?<=\\d)[+-](\\d+)(?=A|T|G|C)"))
  
  #Remove variants with no annotated gene
  filter1 <- vcf_info_plus[complete.cases(vcf_info_plus$Ensemble_ID), ]
  
  #Protein coding transcripts only
  filter2 <- filter1[(filter1$BIOTYPE == "protein_coding"),]
  
  #Intronic only
  filter3 <- filter2[grepl("splice_region_variant|intron_variant|splice_acceptor_variant|splice_donor_variant",
                           filter2$Consequence),]
  
  #Remove repeat entries (i.e. same variant across multiple transcripts). 
  # Find all duplicated variants and their corresponding rows
  rownames(filter3) <- seq(nrow(filter3))
  duplicate_variant_index_1 <- which(duplicated(filter3$Variant) | duplicated(filter3$Variant, fromLast = TRUE))
  
  #First, remove any duplicates that have a blank SYMBOL column (i.e. no gene name)
  filter4 <- subset(filter3, !((rownames(filter3) %in% duplicate_variant_index_1) & (SYMBOL == "")))
  
  vcf_info_final <- filter4
  
####################

vcf_info_final$Pangolin_Max <- ifelse(is.na(vcf_info_final$Pangolin_Max), 0, as.numeric(vcf_info_final$Pangolin_Max))
vcf_info_final$SpliceAI_Max <- ifelse(is.na(vcf_info_final$SpliceAI_Max), 0, as.numeric(vcf_info_final$SpliceAI_Max))
vcf_info_final$CADD_Score <- ifelse(is.na(vcf_info_final$CADD_Score), 0, as.numeric(vcf_info_final$CADD_Score))
vcf_info_final$SQUIRLS_Max <- ifelse(is.na(vcf_info_final$SQUIRLS_Max), 0, as.numeric(vcf_info_final$SQUIRLS_Max))

#Add column for basic variant type (SNV or Other) - required for model predictions
vcf_info_final$Variant_Type <- as.character(vcf_info_final$VARIANT_CLASS)
vcf_info_final$Variant_Type <- ifelse(vcf_info_final$Variant_Type == "SNV", "SNV", "Other")
vcf_info_final$Variant_Type <- factor(vcf_info_final$Variant_Type)

#Pathogenic threshold for each tool
A <- quote(SpliceAI_Max >= 0.4475)
B <- quote(Pangolin_Max >= 0.2175)
C <- quote(CADD_Score >= 21.175)
D <- quote(SQUIRLS_Max >= 0.023)

#Load models
reg_model <- readRDS("INSERT_LINK") 
tree <- readRDS("INSERT_LINK")
rf_model <- readRDS("INSERT_LINK")

#Multiple Tools
{vcf_info_final$Four_Tools <- with(vcf_info_final, ifelse(((eval(A) & eval(B) & eval(C) & eval(D))), 1, 0))
  
  vcf_info_final$Three_Tools <- with(vcf_info_final, ifelse(((eval(A) & eval(B) & eval(C)) | (eval(A) & eval(B) & eval(D)) | (eval(A) & eval(C) & eval(D)) | 
                                                                   (eval(B) & eval(C) & eval(D))), 1, 0))
  
  vcf_info_final$Two_Tools <- with(vcf_info_final, ifelse(((eval(A) & eval(B)) | (eval(A) & eval(C)) | (eval(A) & eval(D)) | (eval(B) & eval(C)) | 
                                                                 (eval(B) & eval(D)) | (eval(C) & eval(D))), 1, 0))
  
  vcf_info_final$One_Tool <- with(vcf_info_final, ifelse(eval(A) | eval(B) | eval(C) | eval(D) , 1, 0))
  
}

#Regression
{reg_pred <- predict(reg_model, vcf_info_final, type = 'response')
  vcf_info_final$Reg_Model_Pred <- ifelse(reg_pred >= 0.5, 'Pathogenic', 'Benign')
}

#Decision Tree
{vcf_info_final$Decision_Tree_Pred <- predict(tree, vcf_info_final, type = 'class')
}

#Random Forest
rf_vcf_info_final <- vcf_info_final[, c('Pangolin_Max', 'SpliceAI_Max', 'SQUIRLS_Max','CADD_Score', 'Variant_Type')] 
vcf_info_final$Random_Forest_Pred <- predict(rf_model, rf_vcf_info_final, type = 'class')
