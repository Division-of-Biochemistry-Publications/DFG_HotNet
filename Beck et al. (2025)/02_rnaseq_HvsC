## SET DIRECTORIES ######################################################################################################
########################################################################################################################
inputPath <- "../RNA_Seq"
outputPath <- file.path(inputPath, "02_Results")
figPath <- file.path(inputPath, "02_Results/00_Figures")

dir.create(outputPath, recursive = TRUE)
dir.create(figPath, recursive = TRUE)



## LIBRARIES ###########################################################################################################
########################################################################################################################
library(DESeq2)
library(tidyverse)
library(readxl)
library(openxlsx)
library(svglite)
library(UpSetR)
library(clusterProfiler)



## IMPORT DATA #########################################################################################################
########################################################################################################################
#Note folder structure: “01_Preprocessed_data” must be a subfolder of the ‘inputPath’ variable 
featurecounts <- read.delim(file.path(inputPath, "01_Preprocessed_data", "07_FeatureCounts_gene_UNIQUE.txt"), row.names = 1, comment.char = "#") # & 07_FeatureCounts_gene_UNIQUE.txt was previously generated with “01_raw_data_preprocessing.sh” or see supplementary data
meta_data <- read.csv(file.path(inputPath, "01_Preprocessed_data", "meta_data.csv"), sep = ";")# see supplementary data

#Note folder structure: documents need to be part of “00_Meta_Data_etc” 
anno <- read.delim(file.path(inputPath, "00_Meta_Data_etc", "DM_1-3_516_R44_potato.v6.1.working_models.func_anno.txt"), header=FALSE) # downloaded from SPUD DB
primary <- read.delim(file.path(inputPath, "00_Meta_Data_etc", "DM_1-3_516_R44_potato.v6.1.repr_hc_gene_models.gff3"), header=FALSE) # downloaded from SPUD DB
potato.go <- read.delim(file.path(inputPath, "00_Meta_Data_etc", "DM_1-3_516_R44_potato.v6.1.working_models.go.txt"), header=FALSE) # downloaded from SPUD DB
AT_GOSLIM <- read.delim(file.path(inputPath, "00_Meta_Data_etc", "ATH_GO_GOSLIM.txt"), header=FALSE) # downloaded from TAIR

# import GO term assignment to higher-ordered groups 
Pathway_groups <- read_excel(file.path(inputPath, "00_Meta_Data_etc", "GO_pathways_grouping.xlsx"), sheet = "groups") # see supplementary data



## ADJUST DATAFRAMES ###################################################################################################
########################################################################################################################
featurecounts <- featurecounts %>%
  dplyr::select(any_of(meta_data$Proben_ID)) %>%
  tibble::rownames_to_column("locusName") %>%
  dplyr::mutate(locusName = sub("\\.v6\\.1$", "", locusName)) %>%
  tibble::column_to_rownames("locusName")

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



## VST data, norm counts, etc. #########################################################################################
########################################################################################################################

# create dds matrix which connects meta_data & featurecounts #
dds = DESeqDataSetFromMatrix(
  countData = featurecounts,
  colData = meta_data,
  design = ~ genotype_dps_treatment
)


# VST TRANSFORMATION + SCALING #
vst_data <- as.data.frame(assay(vst(dds)))
# Z-Score scaling
vst_data_Z <- vst_data %>%
  as.matrix() %>%
  { t(scale(t(.), center = TRUE, scale = TRUE)) } %>%
  as.data.frame() %>%
  tibble::rownames_to_column("locusName")
# check if mean ~0 and SD ~1 (for each gene)
checkup <- as.matrix(vst_data_Z[,-1])
summary(rowMeans(checkup))
summary(matrixStats::rowSds(checkup))


# COUNT NORMALIZATION #
dds_factor <- estimateSizeFactors(dds) # calculates scaling factor (normalization) with which samples are calculated
norm_counts <- counts(dds_factor, normalized = TRUE) # calculates original data with scaling factor (= normalized data)


# REMOVE LOW EXPRESSED GENES (since they are probably not relevant for further analysis) #
norm_filtered10 <- as.data.frame(norm_counts[rowMeans(norm_counts) >= 10,]) # used as input for PCA + DEG calculation
norm_filtered50 <- as.data.frame(norm_counts[rowMeans(norm_counts) >= 50,]) # stricter filter for analysis of individual pathways


# SCALE NORMALIZED + EXPRESSED DATA (for PCA) #
norm_filtered10_t <- as.data.frame(t(norm_filtered10)) # transpose is necessary to scale via columns 
norm_filtered10_Z <- as.data.frame(scale(norm_filtered10_t)) # scale over all samples (in column) for each gene)
print(colMeans(norm_filtered10_Z)) # scales via columns & check if mean ~0 & SD= 1 
print(apply(norm_filtered10_Z, 2, sd))


# DF SAVING #
write.csv(vst_data_Z,file = file.path(outputPath, "vst_Z.csv"), row.names = FALSE)
write.csv(norm_filtered10,file = file.path(outputPath, "norm_filtered10.csv"), row.names = FALSE)
write.csv(norm_filtered50,file = file.path(outputPath, "norm_filtered50.csv"), row.names = FALSE)




################# PRINCIPAL COMPONENT ANALYSIS #########################################################################
########################################################################################################################
## PC CALCULATION
pca_result <- prcomp(norm_filtered10_Z) # calculate PCs, x should contain a 47 x 47 matrix with PCs 
PCs <- as.data.frame(pca_result$x) 
PCs <- rownames_to_column(PCs, var = "Proben_ID")
PCs <- left_join(PCs, meta_data, by = "Proben_ID")

## Scree plot to visualize PC variance contribution (in order to manually ad values for PC1 & PC2 in PCA plot below) #
# generate new df with all 47 PC + variance of the PCs + calculate cum. Variance #
pca_var <- data.frame(
  PC = paste0("PC",1:length(pca_result$sdev)),
  expl_variance = pca_result$sdev**2)
pca_var$pct_expl_var = (pca_var$expl_variance/sum(pca_var$expl_variance))*100
pca_var$cum_expl_var = cumsum(pca_var$pct_expl_var)
pca_var <- pca_var %>% mutate_if(is.numeric, round, digits=1)

# Scree plot
ggplot(pca_var)  +
  geom_col(aes(c(1:nrow(pca_var)), pct_expl_var)) + 
  geom_line(aes(c(1:nrow(pca_var)), cum_expl_var), color = "red") +
  geom_line(aes(c(1:nrow(pca_var)), pct_expl_var), color = "blue") +
  geom_text(aes(label=pct_expl_var, x=c(1:nrow(pca_var)), y=round(pct_expl_var, digits=1)), size = 3.5, colour="black") +
  labs(x = "47 PC" , y = "Explained variance [%]", title = "Scree plot of leaf RNA Seq data (47 samples)") 

# Plot PC1 + PC2
PCs$dps <- factor(PCs$dps, levels = c(3, 7, 13))

ggplot(PCs, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = genotype_treatment, shape = dps), size= 3, alpha = 0.7) +
  scale_color_manual(values = c("#005F73", "#84190E", "#6DC1A5", "#F4C161")) +
  scale_shape_manual(values=c(16, 17, 15)) + 
  ylab("PC2 (13% of explained variance)") + # based on scree plot results 
  xlab("PC1 (34% of explained variance)") + # based on scree plot results 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=8),
        axis.title = element_text(size = 8, face = "bold"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.1,"cm"),
        legend.position = "bottom") + 
  guides(color = guide_legend(title = NULL))

ggsave("PCA.svg", path = figPath, device = "svg",
       width = 7.5, height = 6.5, units = "cm")


### DEGS CALCULATION  ##################################################################################################
########################################################################################################################
deseq <- DESeq(dds, test = c("Wald")) 

log2FC_H_C <- list()
# loop creates a new element in the log2FC_H_C list for each unique genotype_dps pair (i, e.g. “SO_3”, “CE_7”)
# each of these is a df containing the locusName, the genotype-dps pair and the results of differential gene expression between the conditions Heat (_H) and Control (_C)
for(i in unique(meta_data$genotype_dps)) {
  print(i)
  log2FC_H_C[[i]] = data.frame(
    locusName = rownames(featurecounts),
    genotype_dps = i,
    as.data.frame(
      results(
        deseq,
        contrast = c("genotype_dps_treatment",  paste0(i,"_H"),  paste0(i,"_C")),
        pAdjustMethod = "BH"
      )
    )
  )
}
log2FC_H_C <- bind_rows(log2FC_H_C)

DEGs_incl_low_expressed <- log2FC_H_C %>% filter(abs(log2FoldChange) > 1, padj <= 0.05) # filter for statistical significant genes between Heat & Control for both genotypes at each dps (log2FC>1: higher in H; log2FC<1: lower in H)
DEGs_incl_low_expressed <- remove_rownames(DEGs_incl_low_expressed)

# DF SAVING #
write.csv(log2FC_H_C,file = file.path(outputPath, "log2FC_HvsC.csv"), row.names = FALSE)



################# REMOVE LOW EXPRESSED DEGs ############################################################################
########################################################################################################################
# remove DEGs which are low expressed
norm_filtered10 <- rownames_to_column(norm_filtered10, var = "locusName")
DEGs_filtered_10 <- DEGs_incl_low_expressed %>% 
  semi_join(norm_filtered10, by = "locusName") %>% 
  left_join(anno, by = "locusName")

genotype_dps_unique <- unique(DEGs_filtered_10$genotype_dps) # define unique genotype_dps

for (genotype_dps_unique in genotype_dps_unique) { # extracts DEG lists for each unique genotype_dps (independent of up/down with heat)
  assign(genotype_dps_unique, DEGs_filtered_10 %>% 
           filter(genotype_dps == genotype_dps_unique))}



################# UPSET PLOT ###########################################################################################
########################################################################################################################
DEGs_H_C <- list(
  SO_3dps = as.vector(SO_3$locusName),
  SO_7dps = as.vector(SO_7$locusName),
  SO_13dps = as.vector(SO_13$locusName),
  CE_3dps = as.vector(CE_3$locusName), 
  CE_7dps = as.vector(CE_7$locusName), 
  CE_13dps = as.vector(CE_13$locusName))

svglite(file.path(figPath, "Upset.svg"),
        width = 8, height = 6.5)

# plots interesting sections
upset(fromList(DEGs_H_C), 
      intersections = list(
        "CE_13dps", "CE_7dps", "CE_3dps",
        "SO_13dps", "SO_7dps", "SO_3dps",
        c("SO_3dps", "SO_7dps", "SO_13dps", "CE_3dps", "CE_7dps", "CE_13dps")),
      sets = c("SO_13dps", "SO_7dps", "SO_3dps",
               "CE_13dps", "CE_7dps", "CE_3dps"), 
      keep.order = TRUE, 
      set.metadata = NULL, 
      matrix.color = "#2D0D8F", main.bar.color = "#46bac2",
      mainbar.y.label = "number of genes within intersections", mainbar.y.max = NULL,
      sets.bar.color = "#B7E4B2", sets.x.label = "number of DEGs (H/C)",
      point.size = 2, line.size = 1, mb.ratio = c(0.7, 0.3),
      expression = NULL, att.pos = NULL, 
      order.by = "freq",  
      decreasing = TRUE,
      show.numbers = "yes", number.angles = 0, 
      text.scale = c(1.3, 1.2, 1.3, 1.1, 1.2, 1.1), # here the different text sizes can be adapted
      cutoff = NULL, queries = NULL, query.legend = "none",
      shade.color = "gray88", shade.alpha = 0.25, matrix.dot.alpha = 0.5,
      empty.intersections = NULL, color.pal = 1, boxplot.summary = NULL,
      attribute.plots = NULL, scale.intersections = "identity",
      scale.sets = "identity", 
      set_size.angles = 0,
      set_size.show = TRUE, set_size.numbers_size = NULL,
      set_size.scale_max = NULL)

dev.off()



################# Extract DEG lists & different intersections from UPSET PLOT ###########################################
########################################################################################################################
# put all DEG into one list and assign names
all_dfs <- list(
  SO_3  = SO_3,
  SO_7  = SO_7,
  SO_13 = SO_13,
  CE_3  = CE_3,
  CE_7  = CE_7,
  CE_13 = CE_13
)

# inner join across all data frames (to extract common DEGs) 
joined <- Reduce(function(x, y) inner_join(x, y, by = "locusName"), all_dfs)

# only keep locusName + log2FoldChange-column 
lfc_cols <- grep("log2FoldChange", names(joined), value = TRUE)
joined <- joined %>%
  select(locusName, all_of(lfc_cols)) %>%
  setNames(c("locusName", paste0(rep(c("SO","SO","SO","CE","CE","CE"), # adjust name of columns
                                     each = 1),
                                 "_",
                                 c("3","7","13","3","7","13"),
                                 "dps"))) %>%
  left_join(anno, by = "locusName") # add annotation 

# common up- / down-regulated
all_common_up <- joined %>%
  filter(if_all(ends_with("dps"), ~ .x >= 1))

all_common_down <- joined %>%
  filter(if_all(ends_with("dps"), ~ .x <= -1))

# save all DEG lists in one Excel file 
wb <- createWorkbook()

# Liste aller DataFrames und ihrer gewünschten Sheet-Namen
dfs <- list("SO_3" = SO_3, "SO_7" = SO_7, "SO_13" = SO_13, 
            "CE_3" = CE_3, "CE_7" = CE_7, "CE_13" = CE_13,
            "common_up"   = all_common_up,
            "common_down" = all_common_down)

# add each dfs as a seperate sheet
for (name in names(dfs)) {
  addWorksheet(wb, name)           # add sheet
  writeData(wb, name, dfs[[name]]) 
}

# save excel file in the desired path
saveWorkbook(wb, file = file.path(outputPath, "DEGs_H_C.xlsx"), overwrite = TRUE)


################# Preparations for GO Enrichments #######################################################################
########################################################################################################################
DEGs_filtered <- bind_rows(SO_3, SO_7, SO_13, CE_3, CE_7, CE_13)

# Loop over each value in “genotype_dps”: extract up/ down regulated genes within DEG lists
genotype_dps_unique <- unique(DEGs_filtered$genotype_dps)

for (genotype in genotype_dps_unique) {
  # df for log2FoldChange >= 1 (up-regulated)
  assign(paste0(genotype, "_up"), 
         DEGs_filtered %>% 
           filter(genotype_dps == genotype & log2FoldChange >= 1))
  # df for log2FoldChange <= -1 (down-regulated)
  assign(paste0(genotype, "_down"), 
         DEGs_filtered %>% 
           filter(genotype_dps == genotype & log2FoldChange <= -1))
}

# GO_ID - TERM - ONTOLOGY
BP_goterms <- AT_GOSLIM %>%
  filter(row_number() > 4) %>% 
  dplyr::rename(
    TERM = V5,
    GOID = V6, 
    ONTOLOGY = V8) %>%
  filter(ONTOLOGY == "P") %>%  # filter for ONTOLOGY "BP" = biological process
  dplyr::select(-ONTOLOGY, -V1, -V2, -V3, -V4, -V7, -V9, -V10, -V11, -V12, -V13, -V14, -V15)  

BP_goterms <- unique(BP_goterms)

# only keep important columns
potato.go <- potato.go %>% 
  dplyr::rename(locusName = V2,
                GOID = V5) %>%
  dplyr::select("locusName", "GOID")

# filter for genes within potato.go that are primary & delete other isoforms .1 etc. & keep unique ones 
potato.go <- inner_join(potato.go, primary, by = "locusName")
potato.go <- potato.go %>% mutate(locusName = substr(locusName, 1, nchar(locusName) - 2))
potato.go <- unique(potato.go)

# join Soltu identifier with BP_goterms file (contains: TERM & GOID; filtering for BP was done before)
gene_pathway <- inner_join(potato.go, BP_goterms, by = "GOID")
gene_pathway <- unique(gene_pathway)

## TERM2GENE
gene_pathID <- gene_pathway[,c("GOID", "locusName")] # id2gene 
gene_pathID <- na.omit(gene_pathID) # removes rows with NAs 
gene_pathID <- unique(gene_pathID) # unique: remove doubles

## TERM2NAME
gene_path <- gene_pathway[,c("GOID", "TERM")] # id2name
gene_path <- na.omit(gene_path) # removes rows with NAs 
gene_path <- unique(gene_path) # unique: remove doubles



## PERFORM GO ENRICHMENTS ###############################################################################################
########################################################################################################################
# define function to perform enrichments
perform_enrichment <- function(df, gene_pathID, gene_path) {
  # enrichment analysis
  enr_df <- as.data.frame(enricher(df$locusName,
                                   pvalueCutoff = 0.05,
                                   pAdjustMethod = "bonferroni",
                                   TERM2GENE = gene_pathID,
                                   TERM2NAME = gene_path))
  
  # splitting the geneID column
  enr_geneID <- str_split_fixed(enr_df$geneID, "/", n=Inf)
  enr_df <- cbind(enr_df, enr_geneID)
  enr_df$ID <- NULL  # remove ID column
  
  # adding new columns for plotting
  enr_df$description_factor <- factor(enr_df$Description)
  enr_df$calculated_values <- NA  # add new column for values calculated in the following
  
  # calculate GeneRatio and save them in the new column
  enr_df$calculated_values <- apply(enr_df, 1, function(row) {
    ratio_parts <- strsplit(row["GeneRatio"], "/")[[1]]
    numerator <- as.numeric(ratio_parts[1])
    denominator <- as.numeric(ratio_parts[2])
    return(numerator / denominator)
  })
  
  # sort by the calculated values (descending)
  enr_df <- enr_df[order(-enr_df$calculated_values), ]
  
  # Add genotype_dps column (reformat from df's genotype_dps)
  enr_df$genotype_dps <- gsub("_", " ", unique(df$genotype_dps))  # Replace underscore with space
  enr_df$genotype_dps <- paste0(enr_df$genotype_dps, "dps")       # Add "dps" at the end
  
  # Add regulation column (based on the DataFrame's regulation)
  if (grepl("_up", deparse(substitute(df)))) {
    enr_df$regulation <- "up"
  } else if (grepl("_down", deparse(substitute(df)))) {
    enr_df$regulation <- "down"
  } else {
    enr_df$regulation <- NA  # just in case, although this should never happen if "_up" or "_down" are always present
  }
  
  return(enr_df)
}

# Calling the function for all dfs
SO_3_up_enr <- perform_enrichment(SO_3_up, gene_pathID, gene_path)
SO_3_down_enr <- perform_enrichment(SO_3_down, gene_pathID, gene_path)
SO_7_up_enr <- perform_enrichment(SO_7_up, gene_pathID, gene_path)
SO_7_down_enr <- perform_enrichment(SO_7_down, gene_pathID, gene_path)
SO_13_up_enr <- perform_enrichment(SO_13_up, gene_pathID, gene_path)
SO_13_down_enr <- perform_enrichment(SO_13_down, gene_pathID, gene_path)
CE_3_up_enr <- perform_enrichment(CE_3_up, gene_pathID, gene_path)
CE_3_down_enr <- perform_enrichment(CE_3_down, gene_pathID, gene_path)
CE_7_up_enr <- perform_enrichment(CE_7_up, gene_pathID, gene_path)
CE_7_down_enr <- perform_enrichment(CE_7_down, gene_pathID, gene_path)
CE_13_up_enr <- perform_enrichment(CE_13_up, gene_pathID, gene_path)
CE_13_down_enr <- perform_enrichment(CE_13_down, gene_pathID, gene_path)



## PLOTTING PREPARATIONS ###############################################################################################
########################################################################################################################
# join list of all enriched pathways together as a list
df_list <- list(SO_3_up_enr, SO_3_down_enr, SO_7_up_enr, SO_7_down_enr, 
                SO_13_up_enr, SO_13_down_enr, CE_3_up_enr, CE_3_down_enr, 
                CE_7_up_enr, CE_7_down_enr, CE_13_up_enr, CE_13_down_enr)

cols_to_keep <- c("description_factor", "Count", "p.adjust", "calculated_values", "genotype_dps", "regulation") # columns which should be kept

filtered_dfs <- lapply(df_list, function(df) { # loop to keep relevant columns of each df 
  df[, cols_to_keep, drop = FALSE]})

summary_df <- do.call(rbind, filtered_dfs) # connect all dfs

summary_df_groups <- summary_df %>% 
  filter(Count >= 15) %>% # for plotting: only keep enriched pathways which contain 15 or more DEGs 
  mutate(description_factor = as.character(description_factor)) %>% 
  add_count(description_factor) %>%  # add column "n" with frequency of pathways
  filter(n >= 2) %>%  # only keep rows with pathways which occur twice
  select(-n) %>%   
  left_join(Pathway_groups, by = "description_factor") # add higher-ordered groups for a better overview

colors_groups <- c(
  "abiotic stress" = "#08306B",
  "biotic stress" = "#3D7C99",
  "hormone" =  "#A4D7E1",
  "light perception" = "#F9EC9B",
  "photosynthesis" = "#67876E",
  "development" = "#A3C973",
  "protein" = "#46BAC2",
  "various" =  "#D9D9D9",
  "signalling" = "#B4E4B2",
  "cell wall" = "#B3C6E7")

## SUMMARY ENRICHMENT PLOT ##############################################################################################
########################################################################################################################
# define plotting order
summary_df_groups$genotype_dps <- factor(summary_df_groups$genotype_dps,
                                         levels = c("SO 3dps", "SO 7dps", "SO 13dps", 
                                                    "CE 3dps", "CE 7dps", "CE 13dps"))
summary_df_groups$regulation <- factor(summary_df_groups$regulation, levels = c("up", "down"))

ggplot(data = summary_df_groups, 
       aes(x = calculated_values, 
           y = reorder(description_factor, calculated_values),  # reorder based on calculated_values
           fill = Group)) + 
  geom_col(width = 0.4) +
  geom_text(aes(label = Count), hjust = 0, vjust = 0.5, size = 1.7, color = "#4D4D4D") +  
  scale_fill_manual(values = colors_groups) +
  facet_grid(~regulation ~Group ~ genotype_dps , scales = "free", space = "free") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 0.065)) + 
  labs(
    x = "GeneRatio",
    y = "Enriched Pathway (GO: BP)",
    fill = NULL
  ) + 
  theme_bw() +
  theme(
    panel.spacing = unit(0.1, "lines"),  
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill = NA, colour = NA),
    strip.text.x = element_text(size = 8, face = "bold"),
    strip.text.y = element_text(size = 8, angle = 0, hjust = 0),  
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    axis.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.position = "bottom") +
  guides(color = guide_legend(title = NULL))  

ggsave(file.path(figPath, "GO_enrichmnt_HvsC.svg"), 
       device = "svg", width = 22, height = 13, units = "cm")



## Extract genes from enriched pathways ################################################################################
########################################################################################################################
# define function to extract genes from enriched pathways + join with annotationsfile 
process_enrichment <- function(enr_df,
                               anno,
                               gene_pattern = "Soltu.DM\\.",   # keep only these loci
                               desc_col = "Description",
                               anno_desc_col = "v6.1_description") {
  # add GOID from rownames
  enr_df <- tibble::rownames_to_column(enr_df, var = "GOID")
  
  # extract genes per pathway (Description)
  meta <- enr_df %>% dplyr::select(GOID, !!sym(desc_col))
  
  res_list <- lapply(unique(meta[[desc_col]]), function(path){
    ID <- meta$GOID[meta[[desc_col]] == path][1]  # if multiple entries, take the first one
    
    current <- enr_df %>%
      dplyr::filter(GOID == ID, .data[[desc_col]] == path) %>%
      # remove metric columns -> only gene columns remain, which we then transpose
      dplyr::select(-GOID, -GeneRatio, -BgRatio, -RichFactor, -FoldEnrichment,
                    -zScore, -pvalue, -p.adjust, -qvalue, -geneID, -Count)
    
    current <- as.data.frame(t(current), stringsAsFactors = FALSE)
    colnames(current) <- "locusName"
    
    current %>%
      dplyr::filter(stringr::str_detect(locusName, gene_pattern)) %>%
      dplyr::mutate(Pathway = path, GOID = ID)
  })
  
  genes_df <- dplyr::bind_rows(res_list)
  
  # join annotation and reorder columns
  out <- genes_df %>%
    dplyr::left_join(anno, by = "locusName") %>%
    dplyr::relocate(locusName, all_of(anno_desc_col), Pathway, GOID)
  
  return(out)
}

# call function
result_SO_3_up_enr_SUM  <- process_enrichment(SO_3_up_enr,  anno)
result_SO_7_up_enr_SUM  <- process_enrichment(SO_7_up_enr,  anno)
result_SO_13_up_enr_SUM <- process_enrichment(SO_13_up_enr, anno)
result_SO_3_down_enr_SUM  <- process_enrichment(SO_3_down_enr,  anno)
result_SO_7_down_enr_SUM  <- process_enrichment(SO_7_down_enr,  anno)
result_SO_13_down_enr_SUM <- process_enrichment(SO_13_down_enr, anno)

result_CE_3_up_enr_SUM  <- process_enrichment(CE_3_up_enr,  anno)
result_CE_7_up_enr_SUM  <- process_enrichment(CE_7_up_enr,  anno)
result_CE_13_up_enr_SUM <- process_enrichment(CE_13_up_enr, anno)
result_CE_3_down_enr_SUM  <- process_enrichment(CE_3_down_enr,  anno)
result_CE_7_down_enr_SUM  <- process_enrichment(CE_7_down_enr,  anno)
result_CE_13_down_enr_SUM <- process_enrichment(CE_13_down_enr, anno)

# save all gene lists in one Excel file 
wb <- createWorkbook()

dfs <- list("SO_3_up" = result_SO_3_up_enr_SUM, 
            "SO_7_up" = result_SO_7_up_enr_SUM, 
            "SO_13_up" = result_SO_13_up_enr_SUM, 
            "SO_3_down" = result_SO_3_down_enr_SUM,
            "SO_7_down" = result_SO_7_down_enr_SUM, 
            "SO_13_down" = result_SO_13_down_enr_SUM, 
            
            "CE_3_up" = result_CE_3_up_enr_SUM, 
            "CE_7_up" = result_CE_7_up_enr_SUM, 
            "CE_13_up" = result_CE_13_up_enr_SUM, 
            "CE_3_down" = result_CE_3_down_enr_SUM,
            "CE_7_down" = result_CE_7_down_enr_SUM, 
            "CE_13_down" = result_CE_13_down_enr_SUM
)


for (name in names(dfs)) {
  addWorksheet(wb, name)           
  writeData(wb, name, dfs[[name]]) 
}

# save excel file in desired path
saveWorkbook(wb, file = file.path(outputPath, "Genes_from_overrepresented_pathways.xlsx"), overwrite = TRUE)
