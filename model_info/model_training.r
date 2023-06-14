library("cluster")
library("dplyr")
library("UpSetR")
library("stringr")
library("VariantAnnotation")
library("ggplot2")
library("tidyr")
library("caret")
library("randomForest")
library("rpart")

ref = "hg38"

### SNP ###
{ benign_snp_file = "INSERT FILE"
  benign_snp_vcf <- readVcf(benign_snp_file, ref)
  benign_snp_info <- data.frame(info(benign_snp_vcf))
  benign_snp_info$Status <- "Benign"
  benign_snp_info$CHR <- as.character(seqnames(benign_snp_vcf))
  benign_snp_info$POS <- as.character(start(benign_snp_vcf))
  benign_snp_info$REF <- as.character(ref(benign_snp_vcf))
  benign_snp_info$ALT <- data.frame(alt(benign_snp_vcf))
  benign_snp_info$ALT <- benign_snp_info$ALT$value
  benign_snp_info$Variant <- result <- gsub("\\s+", "", (paste(benign_snp_info$CHR,":",benign_snp_info$POS, ":", benign_snp_info$REF, ">", benign_snp_info$ALT)))
  
  #benign_snp_info$SQUIRLS_FILTER <- grepl("SQUIRLS", rowRanges(benign_snp_vcf)$FILTER)
  
  #Replicate variant rows that have multiple VEP entries (i.e. impact multiple genes)
  benign_snp_info <- benign_snp_info %>% 
    unnest(CSQ) %>% 
    distinct(CSQ, .keep_all = TRUE)
  
  ##VEP CSQ##
  if ("CSQ" %in% colnames(benign_snp_info)){
    print('Processing Ensembl VEP CSQ')
    csq = info(header(benign_snp_vcf))['CSQ',]$Description
    csq_col = unlist(strsplit(csq, ':'))[2]
    csq_fields = unlist(strsplit(csq_col, '\\|'))
    csq_data <- data.frame(matrix(nrow = nrow(benign_snp_info), ncol = length(csq_fields)))
    colnames(csq_data) <- csq_fields
            
    #Update CSQ fields
    for (i in 1:nrow(benign_snp_info)){
      transcript <- unlist(benign_snp_info[i,]$CSQ)
      fields <- unlist(strsplit(transcript, split = '\\|'))
      for (col in 1:length(csq_fields)){
        csq_data[i,col] <- fields[col]}
    }
    benign_snp_info <- cbind(csq_data, benign_snp_info)}  
  
  path_snp_file = ""
  path_snp_vcf <- readVcf(path_snp_file, ref)
  path_snp_info <- data.frame(info(path_snp_vcf))
  path_snp_info$Status <- "Pathogenic"
  path_snp_info$CHR <- as.character(seqnames(path_snp_vcf))
  path_snp_info$POS <- as.character(start(path_snp_vcf))
  path_snp_info$REF <- as.character(ref(path_snp_vcf))
  path_snp_info$ALT <- data.frame(alt(path_snp_vcf))
  path_snp_info$ALT <- path_snp_info$ALT$value
  path_snp_info$Variant <- result <- gsub("\\s+", "", (paste(path_snp_info$CHR,":",path_snp_info$POS, ":", path_snp_info$REF, ">", path_snp_info$ALT)))
  
  #path_snp_info$SQUIRLS_FILTER <- grepl("SQUIRLS", rowRanges(path_snp_vcf)$FILTER)
  
  #Replicate variant rows that have multiple VEP entries (i.e. impact multiple genes)
  path_snp_info <- path_snp_info %>% 
    unnest(CSQ) %>% 
    distinct(CSQ, .keep_all = TRUE)
  
  ##VEP CSQ##
  if ("CSQ" %in% colnames(path_snp_info)){
    print('Processing Ensembl VEP CSQ')
    csq = info(header(path_snp_vcf))['CSQ',]$Description
    csq_col = unlist(strsplit(csq, ':'))[2]
    csq_fields = unlist(strsplit(csq_col, '\\|'))
    csq_data <- data.frame(matrix(nrow = nrow(path_snp_info), ncol = length(csq_fields)))
    colnames(csq_data) <- csq_fields
    
    #Update CSQ fields
    for (i in 1:nrow(path_snp_info)){
      transcript <- unlist(path_snp_info[i,]$CSQ)
      fields <- unlist(strsplit(transcript, split = '\\|'))
      for (col in 1:length(csq_fields)){
        csq_data[i,col] <- fields[col]}
    }
    path_snp_info <- cbind(csq_data, path_snp_info)} 
    
  #Combine dataframes
  snp_info <- rbind(benign_snp_info, path_snp_info)

  #Add intron dist column
  snp_info$intron_dist <- as.numeric(str_extract(snp_info$HGVSc, "(?<=\\d)[+-](\\d+)(?=A|T|G|C)"))
  
  #Protein coding transcripts only
  snp_info_1 <- snp_info[(snp_info$BIOTYPE == "protein_coding"),]
  
  #Intronic only
  snp_info_2 <- snp_info_1[!(is.na(snp_info_1$intron_dist)),]
  
  # <5000 bp into intron
  snp_info_3 <- snp_info_2[abs(snp_info_2$intron_dist)<5000,]
  
  #Remove repeat entries (i.e. same variant across multiple transcripts). 
  # Find all duplicated variants and their corresponding rows
  rownames(snp_info_3) <- seq(nrow(snp_info_3))
  duplicate_variant_index_1 <- which(duplicated(snp_info_3$Variant) | duplicated(snp_info_3$Variant, fromLast = TRUE))
  
  #First, remove any duplicates that have a blank SYMBOL column (i.e. no gene name)
  snp_info_4 <- subset(snp_info_3, !((rownames(snp_info_3) %in% duplicate_variant_index_1) & (SYMBOL == "")))
  
  #Next, keep the row that has the smallest abs(intron_dist) value
  rownames(snp_info_4) <- seq(nrow(snp_info_4))
  duplicate_variant_index_2 <- duplicated(snp_info_4$Variant) | duplicated(snp_info_4$Variant, fromLast = TRUE)
  
  # Identify rows that are duplicated and do not have the minimum abs(intron_dist) value
  rows_to_remove <- which(duplicate_variant_index_2 & snp_info_4$intron_dist != tapply(abs(snp_info_4$intron_dist), snp_info_4$Variant, min)[snp_info_4$Variant])
  snp_info_5 <- snp_info_4[-rows_to_remove, ]
  
  #Keep first occurence of any remaining duplicates
  snp_info_6 <- snp_info_5[!duplicated(snp_info_5$Variant), ]
  
  #Separate intron regions
  snp_info_splice_site <- snp_info_6[(abs(snp_info_6$intron_dist) <= 8),]
  table(snp_info_splice_site$Status)
  
  snp_info_deep_intronic <- snp_info_6[(abs(snp_info_6$intron_dist) > 8),]
  table(snp_info_deep_intronic$Status)
}

### INDEL ###

{ref = "hg38"
  
  benign_indels_file = "INSERT FILE"
  benign_indels_vcf <- readVcf(benign_indels_file, ref)
  benign_indels_info <- data.frame(info(benign_indels_vcf))
  benign_indels_info$Status <- "Benign"
  benign_indels_info$CHR <- as.character(seqnames(benign_indels_vcf))
  benign_indels_info$POS <- as.character(start(benign_indels_vcf))
  benign_indels_info$REF <- as.character(ref(benign_indels_vcf))
  benign_indels_info$ALT <- data.frame(alt(benign_indels_vcf))
  benign_indels_info$ALT <- benign_indels_info$ALT$value
  benign_indels_info$Variant <- result <- gsub("\\s+", "", (paste(benign_indels_info$CHR,":",benign_indels_info$POS, ":", benign_indels_info$REF, ">", benign_indels_info$ALT)))
  
  #benign_indels_info$SQUIRLS_FILTER <- grepl("SQUIRLS", rowRanges(benign_indels_vcf)$FILTER)
  
 #Replicate variant rows that have multiple VEP entries (i.e. impact multiple genes)
  benign_indels_info <- benign_indels_info %>% 
    unnest(CSQ) %>% 
    distinct(CSQ, .keep_all = TRUE)
  
  
  ##VEP CSQ##
  if ("CSQ" %in% colnames(benign_indels_info)){
    print('Processing Ensembl VEP CSQ')
    csq = info(header(benign_indels_vcf))['CSQ',]$Description
    csq_col = unlist(strsplit(csq, ':'))[2]
    csq_fields = unlist(strsplit(csq_col, '\\|'))
    csq_data <- data.frame(matrix(nrow = nrow(benign_indels_info), ncol = length(csq_fields)))
    colnames(csq_data) <- csq_fields
    
    
    
    #Update CSQ fields
    for (i in 1:nrow(benign_indels_info)){
      transcript <- unlist(benign_indels_info[i,]$CSQ)
      fields <- unlist(strsplit(transcript, split = '\\|'))
      for (col in 1:length(csq_fields)){
        csq_data[i,col] <- fields[col]}
    }
    benign_indels_info <- cbind(csq_data, benign_indels_info)} 
  
  path_indels_file = "INSERT FILE"
  path_indels_vcf <- readVcf(path_indels_file, ref)
  path_indels_info <- data.frame(info(path_indels_vcf))
  path_indels_info$Status <- "Pathogenic"
  path_indels_info$CHR <- as.character(seqnames(path_indels_vcf))
  path_indels_info$POS <- as.character(start(path_indels_vcf))
  path_indels_info$REF <- as.character(ref(path_indels_vcf))
  path_indels_info$ALT <- data.frame(alt(path_indels_vcf))
  path_indels_info$ALT <- path_indels_info$ALT$value
  path_indels_info$Variant <- result <- gsub("\\s+", "", (paste(path_indels_info$CHR,":",path_indels_info$POS, ":", path_indels_info$REF, ">", path_indels_info$ALT)))
  
  #path_indels_info$SQUIRLS_FILTER <- grepl("SQUIRLS", rowRanges(path_indels_vcf)$FILTER)
  
  #Replicate variant rows that have multiple VEP entries (i.e. impact multiple genes)
  path_indels_info <- path_indels_info %>% 
    unnest(CSQ) %>% 
    distinct(CSQ, .keep_all = TRUE)
  
  
  ##VEP CSQ##
  if ("CSQ" %in% colnames(path_indels_info)){
    print('Processing Ensembl VEP CSQ')
    csq = info(header(path_indels_vcf))['CSQ',]$Description
    csq_col = unlist(strsplit(csq, ':'))[2]
    csq_fields = unlist(strsplit(csq_col, '\\|'))
    csq_data <- data.frame(matrix(nrow = nrow(path_indels_info), ncol = length(csq_fields)))
    colnames(csq_data) <- csq_fields
    
     
    #Update CSQ fields
    for (i in 1:nrow(path_indels_info)){
      transcript <- unlist(path_indels_info[i,]$CSQ)
      fields <- unlist(strsplit(transcript, split = '\\|'))
      for (col in 1:length(csq_fields)){
        csq_data[i,col] <- fields[col]}
    }
    path_indels_info <- cbind(csq_data, path_indels_info)} 
  
  
  #Combine dataframes
  indels_info <- rbind(benign_indels_info, path_indels_info)
  
  #Protein coding transcripts only
  indels_info_1 <- indels_info[(indels_info$BIOTYPE == "protein_coding"),]
  
  #Intronic only
  indels_info_2 <- indels_info[indels_info$Consequence %in% c("intron_variant", "splice_acceptor_variant", "spice_donor_variant", 
                                                              "splice_acceptor_variant&intron_variant",  "splice_donor_variant&intron_variant",
                                                              "splice_donor_variant&splice_acceptor_variant&frameshift_variant&stop_lost&intron_variant",
                                                              "splice_region_variant&intron_variant"),]
  
  #Add intron dist column
  indels_info_2$intron_dist <- as.numeric(str_extract(indels_info_2$HGVSc, "(?<=-|\\+)\\d+"))
  
  # <5000 bp into intron
  indels_info_3 <- indels_info_2[abs(indels_info_2$intron_dist)<5000,]
  
  #Remove repeat entries (i.e. same variant across multiple transcripts). 
  # Find all duplicated variants and their corresponding rows
  rownames(indels_info_3) <- seq(nrow(indels_info_3))
  duplicate_variant_index_1 <- which(duplicated(indels_info_3$Variant) | duplicated(indels_info_3$Variant, fromLast = TRUE))
  
  #First, remove any duplicates that have a blank SYMBOL column (i.e. no gene name)
  indels_info_4 <- subset(indels_info_3, !((rownames(indels_info_3) %in% duplicate_variant_index_1) & (SYMBOL == "")))
  
  #Next, keep the row that has the smallest abs(intron_dist) value
  rownames(indels_info_4) <- seq(nrow(indels_info_4))
  duplicate_variant_index_2 <- duplicated(indels_info_4$Variant) | duplicated(indels_info_4$Variant, fromLast = TRUE)
  
  # Identify rows that are duplicated and do not have the minimum abs(intron_dist) value
  rows_to_remove <- which(duplicate_variant_index_2 & indels_info_4$intron_dist != tapply(abs(indels_info_4$intron_dist), indels_info_4$Variant, min)[indels_info_4$Variant])
  indels_info_5 <- indels_info_4[-rows_to_remove, ]
  
  #Keep first occurence of any remaining duplicates
  indels_info_6 <- indels_info_5[!duplicated(indels_info_5$Variant), ]
  indels_info_6 <- indels_info_6[!(is.na(indels_info_6$Variant)), ]
  
  #Separate intron regions
  indels_info_splice_site <- indels_info_6[(abs(indels_info_6$intron_dist) <= 8),]
  indels_info_deep_intronic <- indels_info_6[(abs(indels_info_6$intron_dist) > 8),]
  
}

#### Combine variant sets ######
{all_variants <- rbind(indels_info_6, snp_info_6)

all_variants$Pangolin_Max <- as.numeric(all_variants$Pangolin_Max)
all_variants$SpliceAI_Max <- as.numeric(all_variants$SpliceAI_Max)
all_variants$CADD_Score <- as.numeric(all_variants$CADD_Score)
all_variants$SQUIRLS_Max <- as.numeric(all_variants$SQUIRLS_Max)
all_variants$VARIANT_CLASS <- as.factor(all_variants$VARIANT_CLASS)
all_variants$Status <- as.factor(all_variants$Status)

#Add binary variant type column (i.e. SNV or Other) - required for models that will encounter unseen variant classes in future
all_variants$Variant_Type <- as.character(all_variants$VARIANT_CLASS)
all_variants$Variant_Type <- ifelse(all_variants$Variant_Type == "SNV", "SNV", "Other")
all_variants$Variant_Type <- factor(all_variants$Variant_Type)
table(all_variants$Variant_Type)
}

### Split based on intron distance #####
{all_splice_site <- all_variants[(abs(all_variants$intron_dist) <= 8),]
all_deep_intronic <- all_variants[(abs(all_variants$intron_dist) > 8),]
}

### Create graphs of intron dist ###

{### Graph 1 = Deep Intronic ###
  bin_size <- 50
  all_deep_intronic_binned <- all_deep_intronic %>%
    mutate(bin = cut(intron_dist, seq(-5000, 5000, bin_size), include.lowest = TRUE, labels = FALSE)) %>%
    mutate(bin = ifelse(bin == 0, -5000, bin)) %>%
    mutate(bin = ifelse(bin > 0, (-5000 + (bin - 1) * bin_size), (-5000 - bin * bin_size))) %>%
    group_by(Status, bin) %>%
    summarise(Count = n(), .drop = FALSE) %>%
    ungroup() %>%
    complete(Status, bin, fill = list(Count = 0))
  
  #Set fill color for zero values to transparent
  fill_colors <- c("deepskyblue", "red")
  zero_color <- "transparent"
  names(fill_colors) <- levels(all_deep_intronic_binned$Status)
  fill_colors[is.na(fill_colors)] <- zero_color
  
  #Plot binned data on log10 scale with truncated counts
  deep_p1 <- ggplot(all_deep_intronic_binned, aes(x = bin, y = Count, fill = Status)) +
    geom_col(position = "dodge") +
    scale_x_continuous(limits = c(-5005, 0)) +
    scale_y_continuous(trans = "log10", limits = c(-1, 10000), expand = c(0, 0),
                       breaks = c(1, 10, 100, 1000, 10000), labels = c("1", "10", "100", "1,000", "10,000")) +
    labs(x = "Distance from intron (nucleotides)", y = "Variant count") +
    scale_fill_manual(values = fill_colors) +
    theme_classic()
  
  ggsave("deep_p1.png", plot = deep_p1, width = 8, height = 6, dpi = 300)
  
  
  deep_p2 <- ggplot(all_deep_intronic_binned, aes(x = bin, y = Count, fill = Status)) +
    geom_col(position = "dodge") +
    scale_x_continuous(limits = c(0, 5005)) +
    scale_y_continuous(trans = "log10", limits = c(-1, 10000), expand = c(0, 0),
                       breaks = c(1, 10, 100, 1000, 10000), labels = c("1", "10", "100", "1,000", "10,000")) +
    labs(x = "Distance from intron (nucleotides)", y = "Variant count") +
    scale_fill_manual(values = fill_colors) +
    theme_classic()
  
  ggsave("deep_p2.png", plot = deep_p2, width = 8, height = 6, dpi = 300)
  
  
  ### Graph 2 = Splice Site ###
  #Create a vector of all possible distances
  all_distances <- seq(-8, 8)
  
  #Summarize the counts, include zero counts
  all_splice_site_binned <- all_splice_site %>%
    group_by(Status, intron_dist) %>%
    summarise(Count = n(), .drop = FALSE) %>%
    ungroup() %>%
    complete(Status, intron_dist, fill = list(Count = 0))
  
  #Set fill color for zero values to transparent
  fill_colors <- c("deepskyblue", "red")
  zero_color <- "transparent"
  names(fill_colors) <- levels(all_splice_site_binned$Status)
  fill_colors[is.na(fill_colors)] <- zero_color
  
  #Create the plot
  splice_p1 <- ggplot(all_splice_site_binned, aes(x = intron_dist, y = Count, fill = Status)) +
    geom_col(position = "dodge") +
    scale_x_continuous(limits = c(-8, 0), breaks = all_distances) +
    scale_y_continuous(trans = "log10", limits = c(-1, 1000), expand = c(0, 0)) +
    scale_fill_manual(values = fill_colors) +
    labs(x = "Distance from intron/exon junction (nucleotides)", y = "Variant count") +
    theme_classic()
  
  ggsave("splice_p1.png", plot = splice_p1, width = 8, height = 6, dpi = 300)
  
   
  splice_p2 <- ggplot(all_splice_site_binned, aes(x = intron_dist, y = Count, fill = Status)) +
    geom_col(position = "dodge") +
    scale_x_continuous(limits = c(0, 8), breaks = all_distances) +
    scale_y_continuous(trans = "log10", limits = c(-1, 1000), expand = c(0, 0)) +
    scale_fill_manual(values = fill_colors) +
    labs(x = "Distance from intron/exon junction (nucleotides)", y = "Variant count") +
    theme_classic()
  
  ggsave("splice_p2.png", plot = splice_p2, width = 8, height = 6, dpi = 300)
  
}

table(all_splice_site$Status, all_splice_site$VARIANT_CLASS)
table(all_deep_intronic$Status, all_deep_intronic$VARIANT_CLASS)

##### Split Data (80% train, 20% validation) #####
{set.seed(123)
  benign_model_data <- all_deep_intronic[all_deep_intronic$Status == 'Benign', c('Pangolin_Max', 'SpliceAI_Max', 'SQUIRLS_Max','CADD_Score', 'Variant_Type', 'Status')]
  b_idx <- sample(seq(1, 2), size = nrow(benign_model_data), replace = TRUE, prob = c(.8, .2))
  benign_train <- benign_model_data[b_idx == 1, ]
  benign_val <- benign_model_data[b_idx == 2, ]
  
  path_model_data <- all_deep_intronic[all_deep_intronic$Status == 'Pathogenic', c('Pangolin_Max', 'SpliceAI_Max', 'SQUIRLS_Max','CADD_Score', 'Variant_Type', 'Status')]
  p_idx <- sample(seq(1, 2), size = nrow(path_model_data), replace = TRUE, prob = c(.8, .2))
  path_train <- path_model_data[p_idx == 1, ]
  path_val <- path_model_data[p_idx == 2, ]
  
  train_data <- rbind(benign_train, path_train)
  train_data[is.na(train_data)] <- 0
  
  val_data <- rbind(benign_val, path_val)
  val_data[is.na(val_data)] <- 0
}

table(train_data$Status)
#423/(423+12033)

table(val_data$Status)
#94/(94+2986)

#####################################

### Optimize Individual Tools via Cross-Validation ###
{ set.seed(123)
  
  file_path <- "results_summary.txt"
  file_conn <- file(file_path, open = "w")

  # Function to calculate precision, recall and F1 scores
  {F1_Score <- function(true_labels, predicted_labels, positive = "Pathogenic") {
    tp <- sum(true_labels == positive & predicted_labels == positive)
    fp <- sum(true_labels != positive & predicted_labels == positive)
    fn <- sum(true_labels == positive & predicted_labels != positive)
    
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    f1_score <- 2 * precision * recall / (precision + recall)
    
    return(c(precision, recall, f1_score))
  }
    
    # Create list of parameter grids to search over for each tool
    params <- list(
      SpliceAI_Max = c(data.frame(thresh = seq(0, 1, by = 0.0025)), default = 0.5),
      Pangolin_Max = c(data.frame(thresh = seq(0, 1, by = 0.0025)), default = 0.2),
      CADD_Score = c(data.frame(thresh = seq(0, 40, by = 0.025)), default = 20),
      SQUIRLS_Max = c(data.frame(thresh = seq(0, 0.1, by = 0.0001)), default = 0)
    )
    
    # Create folds
    set.seed(123) # Set seed for reproducibility
    folds <- createFolds(train_data$Status, k = 5, list = TRUE, returnTrain = FALSE)
    
    # Initialize results dataframe
    results <- data.frame(Tool = character(),
                          Fold = integer(),
                          Threshold = double(),
                          Precision = double(),
                          Recall = double(),
                          F1_Score = double(),
                          stringsAsFactors = FALSE)
    
    # Loop over folds and parameter grids
    for (fold in seq_along(folds)) {
      for (tool in names(params)) {
        # Subset data to current fold and tool
        fold_data <- train_data[folds[[fold]], ]
        fold_data <- fold_data[, c(tool, "Status")]
        
        # Loop over thresholds
        for (thresh in params[[tool]]$thresh) {
          # Make predictions based on current threshold
          pred_labels <- ifelse(fold_data[[tool]] > thresh, "Pathogenic", "Benign")
          
          # Calculate precision, recall and F1 scores
          scores <- F1_Score(fold_data$Status, pred_labels)
          
          # Add results to dataframe
          results <- rbind(results, data.frame(Tool = tool,
                                               Fold = fold,
                                               Threshold = thresh,
                                               Precision = scores[1],
                                               Recall = scores[2],
                                               F1_Score = scores[3],
                                               stringsAsFactors = FALSE))
        }
      }
    }
    
    # Group results by Tool, Threshold, and Fold and calculate mean values
    mean_results <- results %>%
      group_by(Tool, Threshold) %>%
      summarize(mean_precision = mean(Precision),
                mean_recall = mean(Recall),
                mean_f1_score = mean(F1_Score))
    
    # Iterate over each tool and create a plot for each one
    plots <- lapply(unique(mean_results$Tool), function(tool_name) {
      
      # Extract relevant information from mean_results
      tool_df <- mean_results %>% filter(Tool == tool_name)
      threshold <- tool_df$Threshold
      mean_precision <- tool_df$mean_precision
      mean_recall <- tool_df$mean_recall
      mean_f1_score <- tool_df$mean_f1_score
      
      max_f1 <- max(tool_df$mean_f1_score, na.rm = TRUE)
      opti_thresh <- max(tool_df$Threshold[tool_df$mean_f1_score == max_f1], na.rm = TRUE)
      
      # Write the results to summary text file
      writeLines(paste("Results for", tool_name), file_conn)
      writeLines(paste("Max F1 score =", round(max_f1, 2)), file_conn)
      writeLines(paste("Optimum Threshold =", round(opti_thresh, 5)), file_conn)
      writeLines("", file_conn)
      
      # Create the plot
      p <- ggplot(data = tool_df, aes(x = threshold)) +
        geom_smooth(aes(y = mean_precision, color = "Precision"), size = 1.5, se = FALSE, method = "loess", span = 0.3) +
        geom_smooth(aes(y = mean_recall, color = "Recall"), size = 1.5, se = FALSE, method = "loess", span = 0.3) +
        geom_smooth(aes(y = mean_f1_score, color = "F1 Score"), size = 1.5, se = FALSE, method = "loess", span = 0.3) +
        labs(title = paste(tool_name, "Metrics vs. Threshold"), x = "Threshold", y = "Score") +
        scale_color_manual(name = "Metric", values = c("Precision" = "#619CFF", "Recall" = "#F8766D", "F1 Score" = "#00BA38")) +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
        geom_point(aes(x = opti_thresh, y = max_f1), shape = "*", size = 5) +
        geom_vline(xintercept = params[[tool_name]]$default, linetype = "dotted", color = "black", alpha = 0.75) +
        theme_minimal() +
        theme(panel.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      
      filename <- paste0(tool_name, ".png")
      
      ggsave(filename, plot = p, width = 10, height = 6, dpi = 300)
      
    })
    
  }
  
  close(file_conn)
  
}

{ A <- quote(SpliceAI_Max >= 0.4475)
  B <- quote(Pangolin_Max >= 0.2175)
  C <- quote(CADD_Score >= 21.175)
  D <- quote(SQUIRLS_Max >= 0.023)
  
  train_data$SpliceAI_Class <- with(train_data, ifelse(eval(A), 1, 0))
  train_data$Pangolin_Class <- with(train_data, ifelse(eval(B), 1, 0))
  train_data$CADD_Class <- with(train_data, ifelse(eval(C), 1, 0))
  train_data$SQUIRLS_Class <- with(train_data, ifelse(eval(D), 1, 0))

upset(train_data[train_data$Status == "Pathogenic",c("SpliceAI_Class", "Pangolin_Class", "CADD_Class", "SQUIRLS_Class")], 
      order.by = "freq")
}


### Check performance on validation data ###
# Sample dataframe
performance_summary <- data.frame(
  Tool = c("SpliceAI_Class", "Pangolin_Class", "CADD_Class", "SQUIRLS_Class", "One_Tool", "Two_Tools", 
           "Three_Tools", "Four_Tools", "Reg_Model_Pred", "Decision_Tree_Pred", "Random_Forest_Pred"),
  Precision = numeric(length = 11),
  Recall = numeric(length = 11),
  F1 = numeric(length = 11)
)

{data <- val_data
  A <- quote(SpliceAI_Max >= 0.4475)
  B <- quote(Pangolin_Max >= 0.2175)
  C <- quote(CADD_Score >= 21.175)
  D <- quote(SQUIRLS_Max >= 0.023)
  
  data$SpliceAI_Class <- with(data, ifelse(eval(A), 1, 0))
  data$Pangolin_Class <- with(data, ifelse(eval(B), 1, 0))
  data$CADD_Class <- with(data, ifelse(eval(C), 1, 0))
  data$SQUIRLS_Class <- with(data, ifelse(eval(D), 1, 0))
  
  p <- data[data$Status == 'Pathogenic',]
  b <- data[data$Status == 'Benign',]
  
  for (col in c('SpliceAI_Class', 'Pangolin_Class', 'CADD_Class', 'SQUIRLS_Class')){
    print(col)
    TP <- (sum(p[,col] == 1))
    FN <- sum(p[,col] == 0)
    TN <- sum(b[,col] == 0)
    FP <- sum(b[,col] == 1)
    Precision <- TP / (TP + FP)
    Recall <- TP / (TP + FN)
    F1 <- ( 2 * ((Precision * Recall) / (Precision + Recall)))
    performance_summary[performance_summary$Tool == col, 'Precision'] = Precision
    performance_summary[performance_summary$Tool == col, 'Recall'] = Recall
    performance_summary[performance_summary$Tool == col, 'F1'] = F1
    print(TP)
    print(FN)
    print(TN)
    print(FP)
  }
  
  upset(data[data$Status == "Pathogenic",c("SpliceAI_Class", "Pangolin_Class", "CADD_Class", "SQUIRLS_Class")], 
        order.by = "freq")
}

##### Multiple Tools ######
{data$Four_Tools <- with(data, ifelse(((eval(A) & eval(B) & eval(C) & eval(D))), 1, 0))

data$Three_Tools <- with(data, ifelse(((eval(A) & eval(B) & eval(C)) | (eval(A) & eval(B) & eval(D)) | (eval(A) & eval(C) & eval(D)) | 
                                         (eval(B) & eval(C) & eval(D))), 1, 0))

data$Two_Tools <- with(data, ifelse(((eval(A) & eval(B)) | (eval(A) & eval(C)) | (eval(A) & eval(D)) | (eval(B) & eval(C)) | 
                                       (eval(B) & eval(D)) | (eval(C) & eval(D))), 1, 0))

data$One_Tool <- with(data, ifelse(eval(A) | eval(B) | eval(C) | eval(D) , 1, 0))

p <- data[data$Status == 'Pathogenic',]
b <- data[data$Status == 'Benign',]


for (col in c('Four_Tools', 'Three_Tools', 'Two_Tools', 'One_Tool')){
  print(col)
  TP <- (sum(p[,col] == 1))
  FN <- sum(p[,col] == 0)
  TN <- sum(b[,col] == 0)
  FP <- sum(b[,col] == 1)
  Precision <- TP / (TP + FP)
  Recall <- TP / (TP + FN)
  F1 <- ( 2 * ((Precision * Recall) / (Precision + Recall)))
  performance_summary[performance_summary$Tool == col, 'Precision'] = Precision
  performance_summary[performance_summary$Tool == col, 'Recall'] = Recall
  performance_summary[performance_summary$Tool == col, 'F1'] = F1
}

}


#### Decision Tree #####
{set.seed(123)
  tree <- rpart(Status ~ Pangolin_Max + SpliceAI_Max + SQUIRLS_Max + CADD_Score + Variant_Type, 
                data = train_data, method = "class", 
                parms = list(split = 'information'),
                minsplit = 2)
  
  fancyRpartPlot(tree, main = "Decision Tree for Deep-Intronic Variants",
                 sub = "", yesno = 2, type= 5)
  
  data$Decision_Tree_Pred <- predict(tree, val_data, type = 'class')
  
  p <- data[data$Status == 'Pathogenic',]
  b <- data[data$Status == 'Benign',]
  
  for (col in c('Decision_Tree_Pred')){
    print(col)
    TP <- (sum(p[,col] == 'Pathogenic'))
    FN <- sum(p[,col] == 'Benign')
    TN <- sum(b[,col] == 'Benign')
    FP <- sum(b[,col] == 'Pathogenic')
    Precision <- TP / (TP + FP)
    Recall <- TP / (TP + FN)
    F1 <- ( 2 * ((Precision * Recall) / (Precision + Recall)))
    print(F1)
    performance_summary[performance_summary$Tool == col, 'Precision'] = Precision
    performance_summary[performance_summary$Tool == col, 'Recall'] = Recall
    performance_summary[performance_summary$Tool == col, 'F1'] = F1
  }
}


#### Random Forest ####
{rf_model <- randomForest(Status ~ Pangolin_Max + SpliceAI_Max + CADD_Score + SQUIRLS_Max + Variant_Type, 
                          data=train_data, importance = TRUE)

data$Random_Forest_Pred <- predict(rf_model, val_data, type = 'class')

p <- data[data$Status == 'Pathogenic',]
b <- data[data$Status == 'Benign',]

for (col in c('Random_Forest_Pred')){
  print(col)
  TP <- (sum(p[,col] == 'Pathogenic'))
  FN <- sum(p[,col] == 'Benign')
  TN <- sum(b[,col] == 'Benign')
  FP <- sum(b[,col] == 'Pathogenic')
  Precision <- TP / (TP + FP)
  Recall <- TP / (TP + FN)
  F1 <- ( 2 * ((Precision * Recall) / (Precision + Recall)))
  performance_summary[performance_summary$Tool == col, 'Precision'] = Precision
  performance_summary[performance_summary$Tool == col, 'Recall'] = Recall
  performance_summary[performance_summary$Tool == col, 'F1'] = F1
}

rf_model$importance
print(F1)
}

###### Regression Model ######
# Fit a logistic regression model using glm()
{set.seed(123)
  train_data$Status_Num <- as.numeric(train_data$Status) - 1
  train_data$Variant_Num <- as.numeric(train_data$Variant_Type) - 1
  
  val_data$Variant_Num <- as.numeric(val_data$Variant_Type) - 1
  
  reg_data <- train_data[,c("Pangolin_Max", "SpliceAI_Max", "SQUIRLS_Max", "CADD_Score", "Status_Num", "Variant_Num")]
  
  reg_model <- glm(Status_Num ~ ., data = reg_data, family = binomial)
  reg_pred <- predict(reg_model, val_data, type = 'response')
  
  data$Reg_Model_Pred <- ifelse(reg_pred >= 0.5, 'Pathogenic', 'Benign')
  
  table(data$Status, data$Reg_Model_Pred)
  
  p <- data[data$Status == 'Pathogenic',]
  b <- data[data$Status == 'Benign',]
  
  for (col in c('Reg_Model_Pred')){
    print(col)
    TP <- (sum(p[,col] == 'Pathogenic'))
    FN <- sum(p[,col] == 'Benign')
    TN <- sum(b[,col] == 'Benign')
    FP <- sum(b[,col] == 'Pathogenic')
    Precision <- TP / (TP + FP)
    Recall <- TP / (TP + FN)
    F1 <- ( 2 * ((Precision * Recall) / (Precision + Recall)))
    performance_summary[performance_summary$Tool == col, 'Precision'] = Precision
    performance_summary[performance_summary$Tool == col, 'Recall'] = Recall
    performance_summary[performance_summary$Tool == col, 'F1'] = F1
    print(Precision)
    print(Recall)
    print(F1)
  }
}

#Create performance_summary graph
{# Reshape the dataframe into long format
df_long <- pivot_longer(performance_summary, cols = c(Precision, Recall, F1), names_to = "Metric", values_to = "Value")

#Define custom colors for each metric
metric_colors <- c("Precision" = "#619CFF", "Recall" = "#F8766D", "F1" = "#00BA38")

#Specify the desired order of the metrics
metric_order <- c("Precision", "Recall", "F1")
df_long$Metric <- factor(df_long$Metric, levels = metric_order)

tool_order <- c("SQUIRLS_Class", "SpliceAI_Class", "Pangolin_Class", "CADD_Class", 
                "One_Tool", "Two_Tools", "Three_Tools", "Four_Tools",   
                "Reg_Model_Pred", "Decision_Tree_Pred", "Random_Forest_Pred")

#Convert the Tool variable to a factor with the desired order
df_long$Tool <- factor(df_long$Tool, levels = tool_order)

#Plot the bar chart with custom colors, specified order, and y-axis range
ggplot(df_long, aes(x = Tool, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
  scale_fill_manual(values = metric_colors) +
  labs(x = "Tool", y = "Value") +
  ggtitle("Metrics Comparison by Tool") +
  scale_x_discrete(labels = c("SQUIRLS", "SpliceAI", "Pangolin", "CADD-Splice", 
                              "One Tool", "Two Tools", "Three Tools", "Four Tools",   
                              "Regression Model", "Decision Tree", "Random Forest")) +
  ylim(0, 1) + 
  theme_minimal() +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line())
}

### Splice Site Data Performance ####
{splice_site_data <- all_splice_site[, c('Pangolin_Max', 'SpliceAI_Max', 'SQUIRLS_Max','CADD_Score', 'Variant_Type', 'Status')]
splice_site_data$Variant_Num <- as.numeric(splice_site_data$Variant_Type) - 1
str(splice_site_data)
splice_site_data[is.na(splice_site_data)] <- 0
data <- splice_site_data
}
#Individual Tools
{ A <- quote(SpliceAI_Max >= 0.4475)
  B <- quote(Pangolin_Max >= 0.2175)
  C <- quote(CADD_Score >= 21.175)
  D <- quote(SQUIRLS_Max >= 0.023)
  
  data$SpliceAI_Class <- with(data, ifelse(eval(A), 1, 0))
  data$Pangolin_Class <- with(data, ifelse(eval(B), 1, 0))
  data$CADD_Class <- with(data, ifelse(eval(C), 1, 0))
  data$SQUIRLS_Class <- with(data, ifelse(eval(D), 1, 0))
  
  p <- data[data$Status == 'Pathogenic',]
  b <- data[data$Status == 'Benign',]
  
  for (col in c('SpliceAI_Class', 'Pangolin_Class', 'CADD_Class', 'SQUIRLS_Class')){
    print(col)
    TP <- (sum(p[,col] == 1))
    FN <- sum(p[,col] == 0)
    TN <- sum(b[,col] == 0)
    FP <- sum(b[,col] == 1)
    Precision <- TP / (TP + FP)
    Recall <- TP / (TP + FN)
    F1 <- ( 2 * ((Precision * Recall) / (Precision + Recall)))
  }
  
  upset(data[data$Status == "Pathogenic",c("SpliceAI_Class", "Pangolin_Class", "CADD_Class", "SQUIRLS_Class")], 
        order.by = "freq")
}

#Multiple Tools
{data$Four_Tools <- with(data, ifelse(((eval(A) & eval(B) & eval(C) & eval(D))), 1, 0))
  
  data$Three_Tools <- with(data, ifelse(((eval(A) & eval(B) & eval(C)) | (eval(A) & eval(B) & eval(D)) | (eval(A) & eval(C) & eval(D)) | 
                                           (eval(B) & eval(C) & eval(D))), 1, 0))
  
  data$Two_Tools <- with(data, ifelse(((eval(A) & eval(B)) | (eval(A) & eval(C)) | (eval(A) & eval(D)) | (eval(B) & eval(C)) | 
                                         (eval(B) & eval(D)) | (eval(C) & eval(D))), 1, 0))
  
  data$One_Tool <- with(data, ifelse(eval(A) | eval(B) | eval(C) | eval(D) , 1, 0))
  
  p <- data[data$Status == 'Pathogenic',]
  b <- data[data$Status == 'Benign',]
  
  
  for (col in c('Four_Tools', 'Three_Tools', 'Two_Tools', 'One_Tool')){
    print(col)
    TP <- (sum(p[,col] == 1))
    FN <- sum(p[,col] == 0)
    TN <- sum(b[,col] == 0)
    FP <- sum(b[,col] == 1)
    Precision <- TP / (TP + FP)
    Recall <- TP / (TP + FN)
    F1 <- ( 2 * ((Precision * Recall) / (Precision + Recall)))
    print(TP)
    print(FN)
    print(TN)
    print(FP)
    print(F1)
  }
}

#Regression
{reg_pred <- predict(reg_model, splice_site_data, type = 'response')
  
  data$Reg_Model_Pred <- ifelse(reg_pred >= 0.5, 'Pathogenic', 'Benign')
  
  p <- data[data$Status == 'Pathogenic',]
  b <- data[data$Status == 'Benign',]
  
  for (col in c('Reg_Model_Pred')){
    print(col)
    TP <- (sum(p[,col] == 'Pathogenic'))
    FN <- sum(p[,col] == 'Benign')
    TN <- sum(b[,col] == 'Benign')
    FP <- sum(b[,col] == 'Pathogenic')
    Precision <- TP / (TP + FP)
    Recall <- TP / (TP + FN)
    F1 <- ( 2 * ((Precision * Recall) / (Precision + Recall)))
    print(TP)
    print(FN)
    print(TN)
    print(FP)
    print(F1)
  }
}

#Decision Tree
{data$Decision_Tree_Pred <- predict(tree, splice_site_data, type = 'class')
  
  p <- data[data$Status == 'Pathogenic',]
  b <- data[data$Status == 'Benign',]
  
  for (col in c('Decision_Tree_Pred')){
    print(col)
    TP <- (sum(p[,col] == 'Pathogenic'))
    FN <- sum(p[,col] == 'Benign')
    TN <- sum(b[,col] == 'Benign')
    FP <- sum(b[,col] == 'Pathogenic')
    Precision <- TP / (TP + FP)
    Recall <- TP / (TP + FN)
    F1 <- ( 2 * ((Precision * Recall) / (Precision + Recall)))
    print(TP)
    print(FN)
    print(TN)
    print(FP)
    print(F1)
  }
}

#Random Forest
{data$Random_Forest_Pred <- predict(rf_model, splice_site_data, type = 'class')
  
  p <- data[data$Status == 'Pathogenic',]
  b <- data[data$Status == 'Benign',]
  
  for (col in c('Random_Forest_Pred')){
    print(col)
    TP <- (sum(p[,col] == 'Pathogenic'))
    FN <- sum(p[,col] == 'Benign')
    TN <- sum(b[,col] == 'Benign')
    FP <- sum(b[,col] == 'Pathogenic')
    Precision <- TP / (TP + FP)
    Recall <- TP / (TP + FN)
    F1 <- ( 2 * ((Precision * Recall) / (Precision + Recall)))
    print(TP)
    print(FN)
    print(TN)
    print(FP)
    print(F1)
  }
}
