#Load Packages
library(tidyverse)
library(plyr)
library(dplyr)
library(gplots)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(stats)
library(GOSemSim)

################################################################
# GO analysis of DEPs with FC 1.2 and p<0.05 in the YNHH Cohort
################################################################

# SET THE DESIRED ORGANISM HERE
library(org.Hs.eg.db)
hs = org.Hs.eg.db
universe = read.csv("Universe.csv", header=TRUE)
universe.list= select(hs,  keys= universe$Protein_ID, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
universe.list
View(universe.list)
write.csv(universe.list, "Universe.list.csv")

#Upregulated genes
up.geneset = read.csv("Up.geneset.csv", header=TRUE)
up.geneset.list= select(hs,  keys= up.geneset$Protein_ID, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
up.geneset.list
View(up.geneset.list)
write.csv(up.geneset.list, "Up.geneset.list.csv")

up.geneset.list= read.csv("Up.geneset.list.csv", header= TRUE)
universe.list= read.csv("Universe.list.csv", header=TRUE)

#GO Enrichment for upregulated proteins
up.go.enrich= enrichGO(up.geneset.list$ENTREZID, universe= universe.list$ENTREZID, OrgDb= "org.Hs.eg.db", keyType = "ENTREZID", 
                       ont = "all", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE, pool = TRUE)

up.go.enrich.df= as.data.frame(up.go.enrich)
write.csv(up.go.enrich.df, "Up_GO_Enrich_LVO_YNHH.csv")
  
#Upset plot
X11()
upsetplot(up.go.enrich)

#Barplot
X11()
barplot(up.go.enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Pathways in Up Proteins- LVO vs. no LVO (YNHH Cohort)",
        font.size = 8) + facet_grid(ONTOLOGY ~ ., scales="free")

#Dotplot
X11()
dotplot(up.go.enrich, title = "GO Pathways of Up Proteins- LVO vs. no LVO (YNHH Cohort)", font.size = 8)

##############################################################
#Dotplot of Upregulated Concordant Pathways (Filtered GO IDs)
##############################################################

#Define the GO IDs we want
go_list <- c("GO:0062023", "GO:0097120", "GO:0050867", "GO:0099072", "GO:0002696", "GO:0050817", "GO:0018209", "GO:0021782", "GO:0010001",
             "GO:0002699", "GO:0002274", "GO:0018108")

#Filter the enrichment results
filtered_res <- up.go.enrich
filtered_res@result <- filtered_res@result %>%
  dplyr::filter(ID %in% go_list)

# Dot plot of filtered results
showCategory = length(go_list)

p <- dotplot(filtered_res, showCategory = length(go_list), title = "GO Pathways for Upregulated Proteins- LVO vs. no LVO (YNHH Cohort)",
  font.size = 8)

# Save as PNG
ggsave(
  filename = "Concordant GO_Dotplot_Up Proteins_YNHH.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

# Save as PDF
ggsave(
  filename = "Concordant GO_Dotplot_Up Proteins_YNHH.pdf",
  plot = p,
  width = 8,
  height = 6
)

################################################################################################################################################

#Downregulated genes
down.geneset = read.csv("Down.geneset.csv", header=TRUE)
down.geneset.list= select(hs,  keys= down.geneset$Protein_ID, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
down.geneset.list
View(down.geneset.list)
write.csv(down.geneset.list, "Down.geneset.list.csv")

down.geneset.list= read.csv("Down.geneset.list.csv", header= TRUE)
universe.list= read.csv("Universe.list.csv", header=TRUE)

#GO Enrichment for downregulated proteins
down.go.enrich= enrichGO(down.geneset.list$ENTREZID, universe= universe.list$ENTREZID, OrgDb= "org.Hs.eg.db", keyType = "ENTREZID", 
                         ont = "all", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE, pool = TRUE)

down.go.enrich.df= as.data.frame(down.go.enrich)
write.csv(down.go.enrich.df, "Down_GO_Enrich_LVO_YNHH.csv")

#Upset plot
X11()
upsetplot(down.go.enrich)

#Barplot
X11()
barplot(down.go.enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Pathways in Down Proteins- LVO vs. no LVO (YNHH Cohort)",
        font.size = 8) + facet_grid(ONTOLOGY ~ ., scales="free")

#Dotplot
X11()
dotplot(down.go.enrich, title = "GO Pathways in Down Proteins- LVO vs. no LVO (YNHH Cohort)", font.size = 8)

################################################################
#Dotplot of Downregulated Concordant Pathways (Filtered GO IDs)
################################################################

#Define the GO IDs we want
go_list <- c("GO:0005179", "GO:0031640", "GO:0141061", "GO:0141060", "GO:0043368")

#Filter the enrichment results
filtered_res <- down.go.enrich
filtered_res@result <- filtered_res@result %>%
  dplyr::filter(ID %in% go_list)

# Dot plot of filtered results
showCategory = length(go_list)

p <- dotplot(filtered_res, showCategory = length(go_list), title = "GO Pathways for Downregulated Proteins- LVO vs. no LVO (YNHH Cohort)",
             font.size = 8)

# Save as PNG
ggsave(
  filename = "Concordant GO_Dotplot_Down Proteins_YNHH.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

# Save as PDF
ggsave(
  filename = "Concordant GO_Dotplot_Down Proteins_YNHH.pdf",
  plot = p,
  width = 8,
  height = 6
)