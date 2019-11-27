#Processing of blast output using BLASTlord (script available at https://github.com/maxlcummins/BLASTlord/tree/master/R)
BLASTlord::blastlord("/Users/maxcummins/Dropbox/Doctorate/Manuscripts/SGI1_BIGSI/References_CDHIT/CD-hit_SGI1-REs.txt", "SGI1_CDhit", 90, 90)

#Trims column names to remove excess characters from gene names
colnames(SGI1_CDhit_full_summary_N90L90) <- gsub(pattern = ".*gene=", "", colnames(SGI1_CDhit_full_summary_N90L90))
colnames(SGI1_CDhit_full_summary_N90L90) <- gsub(pattern = "\\]\\[protein_id.*", "", colnames(SGI1_CDhit_full_summary_N90L90))
colnames(SGI1_CDhit_full_summary_N90L90) <- gsub(pattern = "\\] \\[protein_id.*", "", colnames(SGI1_CDhit_full_summary_N90L90))
colnames(SGI1_CDhit_full_summary_N90L90) <- gsub(pattern = "\\]\\[.*", "", colnames(SGI1_CDhit_full_summary_N90L90))
colnames(SGI1_CDhit_full_summary_N90L90) <- gsub(pattern = "\\] \\[.*", "", colnames(SGI1_CDhit_full_summary_N90L90))
colnames(SGI1_CDhit_full_summary_N90L90) <- gsub(pattern = ".*=", "", colnames(SGI1_CDhit_full_summary_N90L90))

#clone dataframe to one with a smaller name
SGI1_CDhit_full_summary_N90L90 -> df

#change column names so that they are unique (e.g. gene, gene -> gene_1, gene_2)
colnames(df) <- make.unique(colnames(df))

#remove strains that were identified as Salmonella, Proteus or Acinetobacter
df <- df %>% filter(!grepl('Salmonella|Proteus|Acinetobacter', sample_name))

#remove strains that were identified as contaminated
df <- df %>% filter(!grepl('ERR433421|]|ERR473388', sample_name))

#set sample names to rownames so they appear as the tip labels on the dendrogram
df$sample_name -> rownames(df)

#remove the column with sample names
df <- df[,2:ncol(df)]

#if a sample has more than one copy of a given gene change the value to 1
#this is because samples can have multiple integrons which would change their
#genotypic profile and make them seem distinct when they aren't necessarily
df[df > 1] <- 1

#load pheatmap
library(pheatmap)

#used to generate a preliminary dendrogram - not used in final publicaton
pheatmap(df, fontsize_col = 3, fontsize_row = 5)

#create a distance matrix for the starins based on their genotypic profiles
dd <- dist(scale(df), method = "euclidean")

#cluster the strains based on the distance matrix
hc <- hclust(dd, method = "ward.D2")

#create a dendrogram
hcd <- as.dendrogram(hc)

par(mar = c(10, 4.1, 4.1, 2.1))

#change the font size and other cosmetic aspects of the dendrogram
nodePar <- list(lab.cex = .7, pch = c(NA, 19), 
                cex = 0.3, col = "blue")

#plot the dendrogram
plot(hcd, nodePar = nodePar)

#generate a spreadsheet that assists in determining which samples had SGI1-REs
#that assembled to a single scaffold (used for exploratory analysis,
#not published in final manuscript)
`SGI1_CDhit_co-occurence_N90L90` -> df_cooccurence

df_cooccurence$same_scaff <- gsub(pattern = ".*gene=", "", df_cooccurence$same_scaff)
df_cooccurence$same_scaff <- gsub(pattern = "\\]\\[protein_id.*", "", df_cooccurence$same_scaff)
df_cooccurence$same_scaff <- gsub(pattern = "\\] \\[protein_id.*", "", df_cooccurence$same_scaff)
df_cooccurence$same_scaff <- gsub(pattern = "\\]\\[.*", "", df_cooccurence$same_scaff)
df_cooccurence$same_scaff <- gsub(pattern = "\\] \\[.*", "", df_cooccurence$same_scaff)
df_cooccurence$same_scaff <- gsub(pattern = ".*=", "", df_cooccurence$same_scaff)

