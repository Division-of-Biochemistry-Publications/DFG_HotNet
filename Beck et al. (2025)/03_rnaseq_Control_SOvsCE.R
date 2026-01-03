## SET DIRECTORIES ######################################################################################################
########################################################################################################################
inputPath <- "../RNA_Seq"
outputPath <- file.path(inputPath, "03_Results_CONTROL_comparison")

dir.create(outputPath, recursive = TRUE)



## LIBRARIES ###########################################################################################################
########################################################################################################################
library(DESeq2)
library(tidyverse)
library(openxlsx)



## IMPORT DATA #########################################################################################################
########################################################################################################################
featurecounts <- read.delim(file.path(inputPath, "01_Preprocessed_data", "07_FeatureCounts_gene_UNIQUE.txt"), row.names = 1, comment.char = "#") # & 07_FeatureCounts_gene_UNIQUE.txt was previously generated with “01_raw_data_preprocessing.sh” or see supplementary data
meta_data <- read.csv(file.path(inputPath, "01_Preprocessed_data", "meta_data.csv"), sep = ";")# see supplementary data

anno <- read.delim(file.path(inputPath, "00_Meta_Data_etc", "DM_1-3_516_R44_potato.v6.1.working_models.func_anno.txt"), header=FALSE) # downloaded from SPUD DB
primary <- read.delim(file.path(inputPath, "00_Meta_Data_etc", "DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3"), header=FALSE) # downloaded from SPUD DB



## ADJUST DATAFRAMES ###################################################################################################
########################################################################################################################
meta_data <- meta_data %>% filter(treatment == "Control")

# only keep samples, which are also present in meta_data
featurecounts <- featurecounts[,meta_data$Proben_ID] 
featurecounts <- featurecounts[,meta_data$Proben_ID] # fiter for all columns in "featurecounts" which are present in column "Proben_ID" in "meta_data"
featurecounts <- rownames_to_column(featurecounts, var = "locusName")
featurecounts <- featurecounts %>% mutate(locusName = sub("\\.v6\\.1$", "", locusName)) # removes ".v6.1"
featurecounts <- column_to_rownames(featurecounts, var = "locusName")

primary <- primary %>%  # to filter for primary transcript
                  dplyr::filter(
                    V3 == "mRNA",
                    V1 %in% paste0("chr", sprintf("%02d", 1:12)) ) %>%  # generates chr01 … chr12
                  dplyr::mutate(locusName = stringr::str_extract(V9, "(?<=ID=)[^;]+")) %>% 
                  dplyr::select(locusName) #  locusName is primary transcript

anno <- anno %>% 
              dplyr::rename("locusName"= "V1",
                            "v6.1_description" = "V2") %>% 
              inner_join(primary, by = "locusName") %>% 
              mutate(locusName = substr(locusName, 1, nchar(locusName) - 2))



## DIFFERENTIAL EXPRESSION ANALYSES ####################################################################################
########################################################################################################################

# create dds matrix which connects meta_data & featurecounts #
dds = DESeqDataSetFromMatrix(
              countData = featurecounts,
              colData = meta_data,
              design = ~ genotype
            )


# DEGS CALCULATION #
deseq_CE_SO <- DESeq(dds, test = c("Wald")) # compares every "genotype_dps_treatment" with every other using Wald test (see above)

deseq_results = as.data.frame(results(
  object = deseq_CE_SO,  contrast=c("genotype", "Solara", "Cecile"), pAdjustMethod = "BH"
))

deseq_results <- deseq_results %>% filter(abs(log2FoldChange) > 1, padj <= 0.05)
deseq_results <- rownames_to_column(deseq_results, var = "locusName")

SO_higher <- deseq_results %>% filter(log2FoldChange >= 1.0) # pos log2FC means higher expressed in SO
SO_lower <- deseq_results %>% filter(log2FoldChange <= -1.0) # neg log2FC means lower expressed in SO


# COUNT NORMALIZATION
dds_factor <- estimateSizeFactors(dds) # calculates scaling factor (normalization) with which samples are calculated
norm_counts <- counts(dds_factor, normalized = TRUE) # calculates original data with scaling factor (= normalizes data)
norm_counts <- as.data.frame(norm_counts)



################# REMOVE LOW EXPRESSED GENES ############################################################################
########################################################################################################################
norm_filtered <- norm_counts %>%
              mutate(row_name = rownames(norm_counts)) %>% # add row name with seperate column
              rowwise() %>%
              filter(mean(c_across(-row_name)) >= 10) %>%  # filer based on rowmean
              ungroup() %>% 
              column_to_rownames(var = "row_name") # set column "row_name" again as row name

# scale normalized & filtered data
norm_filtered <- as.data.frame(t(norm_filtered)) # transponse is necessary to scale via columns (scale over all samples for each gene)
norm_filtered_Z <- as.data.frame(scale(norm_filtered))
print(colMeans(norm_filtered_Z)) # scales via columns & check if colMean ~0 & SD: 1 
print(apply(norm_filtered_Z, 2, sd))
norm_filtered_scaled <- as.data.frame(t(norm_filtered_Z)) # transpose back
norm_filtered_scaled <- rownames_to_column(norm_filtered_scaled, var = "locusName")

# remove DEGs which are low expressed
SO_higher_DEGs_filtered <- SO_higher %>% 
                                  semi_join(norm_filtered_scaled, by = "locusName") %>% 
                                  left_join(anno, by = "locusName")

SO_lower_DEGs_filtered <- SO_lower %>% 
                                  semi_join(norm_filtered_scaled, by = "locusName") %>% 
                                  left_join(anno, by = "locusName")


# save all DEG lists in one Excel file 
wb <- createWorkbook()

# Liste aller DataFrames und ihrer gewünschten Sheet-Namen
dfs <- list("SO_higher_than_CE_DEGs" = SO_higher_DEGs_filtered, 
            "SO_lower_than_CE_DEGs" = SO_lower_DEGs_filtered)

# add each dfs as a seperate sheet
for (name in names(dfs)) {
  addWorksheet(wb, name)           # add sheet
  writeData(wb, name, dfs[[name]]) 
}

# save excel file in the desired path
saveWorkbook(wb, file = file.path(outputPath, "DEGs_SOvsCE_CONTROL.xlsx"), overwrite = TRUE)
