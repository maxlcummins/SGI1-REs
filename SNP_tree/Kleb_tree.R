library(ggtree)
library(pheatmap)
library(magrittr)
library(dplyr)
library(readr)
library(reshape2)

#read in card data
df <- read_delim("~/Dropbox/Doctorate/Manuscripts/SGI1_BIGSI/final_assemblies/Kleb_Abricate_card/Kleb_card_abricate.txt", 
                 "\t", escape_double = FALSE, trim_ws = TRUE)

#change colnames to be more informative
colnames(df)[c(1,9:10)] <- c("name","perc_coverage", "perc_identity")

#store colnames for reuse
df_colnames <- colnames(df)

#read in card data
df <- read_delim("~/Dropbox/Doctorate/Manuscripts/SGI1_BIGSI/final_assemblies/Kleb_Abricate_card/Kleb_card_abricate.txt", 
                 "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE, skip = 1)

#filter all duplicate header rows by looking for keyword
df <- df %>% filter(X2 != "SEQUENCE")

#apply column names to the dataframe
colnames(df) <- df_colnames

#change data type to numeric
df$perc_coverage <- as.numeric(df$perc_coverage)

#change data type to numeric
df$perc_identity <- as.numeric(df$perc_identity)

#filter for gene hits with >90 cov and identity
df <- df %>% filter(perc_coverage > 90) %>% filter( perc_identity > 90)

#filter out non-card data
df <- df %>% filter(DATABASE %in% c("card"))

#select columns containg sample name and genes identified
df2 <- df %>%select(contains("name"), contains("GENE"))

#recast dataframe to generate a presence absence table
df3 <- dcast(data = df2, name ~ GENE, length, drop = FALSE)

#replace NAs with 0
df3[is.na(df3)] <- 0

#trim everything after "." in the sample names (to remove ".fasta")
rownames(df3) <- gsub("\\..*","",df3$name)

#add genus name to start of sample name and set rownames to this value
rownames(df3) <- gsub("E","Klebsiella_E",rownames(df3))

#remove column with row names
df3 <- df3[,2:ncol(df3)]

#colwise sum
colSums(df3) -> sums

#add the number of samples in dataset positive for a gene, colwise, as well as this number as a percentage
colnames(df3) <- paste0(colnames(df3), " ", sums, "/", nrow(df3), " (", round(sums/nrow(df3)*100), "%)")

#change gene names from APH to other more commonly used names
colnames(df3) <- gsub("APH.*Ib", "strA", colnames(df3))
colnames(df3) <- gsub("APH.*Ia", "aphA1", colnames(df3))
colnames(df3) <- gsub("APH.*Id", "strB", colnames(df3))

#change column order to be alphabetical
df3 <- df3[,order(names(df3))]

#read in the SNP tree
tree <- read.tree(file = "~/Dropbox/Doctorate/Manuscripts/SGI1_BIGSI/November/Kleb_tree.tree")

#create a preliminary tree
p2 <- ggtree(tree) %<+% 
  df3 +
  geom_tiplab(size = 3, align = FALSE) 

#change tip labels to trim off ".fasta"
tree$tip.label <- gsub("\\..*","",tree$tip.label)

#create a data frame from the tree tip labels for downstream analysis
tips <- as.data.frame(tree$tip.label)

#change col names to allow table joining
colnames(tips) <- "name"

#read in spreadsheet with genus designations for the strains
genus_table <- read_table2("~/Dropbox/Doctorate/Manuscripts/SGI1_BIGSI/November/Accession_bioproject_genus")

#join the genus table and the tip labels
new_names <- left_join(tips, genus_table)

#create a new column in the format "Genus_AccesionNo."
new_names$new <- paste0(new_names$Genus, "_", new_names$name)

#Fix the name for reference so it doesnt start with "NA"
new_names$new <- gsub("NA_", "", new_names$new)

#apply new tip labels to the tree
tree$tip.label <- new_names$new

#generate final tree
p <- ggtree(tree) %<+% 
  df3 +
  geom_tiplab(size = 3.5, align = TRUE) +
  geom_treescale()

#plot tree
p

#place next to our previously generated tree a genotypic table
a <- gheatmap(p = p, data = df3,
              colnames_offset_y = 0,
              low = "white",
              high = "blue",
              font.size = 3,
              hjust = 0,
              colnames_position = "top",
              colnames_angle = 90,
              offset = .6,
              width = .8) + ggplot2::ylim(NA, 20) +
  theme(legend.position = "none")

#plot the... plot
a




