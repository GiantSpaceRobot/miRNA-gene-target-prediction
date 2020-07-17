## ---------------------------
##
## Script name:  Common-Gene-Targets.R
##
## Purpose of script: Given miRNA gene target predictions for three miRNAs, compare
##
## Author: Dr. Paul Donovan
##
## Date Created: 17-07-20
##
## Email: pauldonovan@rcsi.com
##
## ---------------------------
##
## Notes: 
##   
##
## ---------------------------

library("ggplot2")
library("VennDiagram")
library("gplots")
library("dplyr")
library("stringr")
library("reshape2")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") # Suppress VennDiagram messages/log file generation

### Use command line arguments
args = commandArgs(trailingOnly=TRUE)

miR16.genes.df <- read.csv(args[1],
                       sep = "\t", 
                       header = F)
miR126.genes.df <- read.csv(args[2], 
                  sep = "\t", 
                  header = F)
miR335.genes.df <- read.csv(args[3], 
                    sep = ",", 
                    header = F)

# miR16.genes.df <- read.csv("/home/paul/Documents/Acute_Ischaemic_Stroke/miR-16-5p/miR-16-5p_all-predicted-gene-targets.tsv",
#                        sep = "\t",
#                        header = F)
# miR126.genes.df <- read.csv("/home/paul/Documents/Acute_Ischaemic_Stroke/miR-126-5p/miR-126-5p_all-predicted-gene-targets.tsv",
#                   sep = "\t",
#                   header = F)
# miR335.genes.df <- read.csv("/home/paul/Documents/Acute_Ischaemic_Stroke/miR-335-5p/miR-335-5p_all-predicted-gene-targets.tsv",
#                     sep = ",",
#                     header = F)

my.output <- args[4]

# Venn Diagram
geneLists <- c("miR-16-5p" = miR16.genes.df,
               "miR-126-5p" = miR126.genes.df,
               "miR-335-5p" = miR335.genes.df)

all.genes <- unique(c(as.vector(miR16.genes.df$V1), as.vector(miR126.genes.df$V1), as.vector(miR335.genes.df$V1)))

my.venn <- venn.diagram(geneLists, 
                        filename = NULL,
                        #fill=c("darkmagenta", "darkblue", "red"), 
                        fill=c("#D90092", "#0758C8", "#FFA500"),
                        alpha=c(0.4,0.4,0.4), 
                        cex = 1.4, 
                        #cat.fontface=1,
                        imagetype = "png",
                        cat.cex = 1.4,
                        cat.dist = c(0.08,0.08,0.03),
                        category.names=c(paste0("miR-16-5p\n(",nrow(miR16.genes.df), ")"), 
                                         paste0("miR-126-5p\n(",nrow(miR126.genes.df), ")"), 
                                         paste0("miR-335-5p (",nrow(miR335.genes.df), ")")),
                        main.cex = 1.4,
                        main=paste0("Total no. of genes predicted = ", 
                                    length(all.genes)))

pdf(file = paste0(args[4], "_VennDiagram.pdf"))
grid.newpage() # Create new grid
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc"))) # Resize this grid (80%)
grid.draw(my.venn) # Draw venn diagram in new smaller plot
dev.off()

### Get Venn intersections
a <- venn(geneLists, show.plot=FALSE)

# By inspecting the structure of the a object created, 
# you notice two attributes: 1) dimnames 2) intersections
# We can store the intersections in a new object named inters
inters <- attr(a,"intersections")

#my.gene.list <- c(inters$`TargetScan.unique.as.character.targetscan.Target.gene..:miRDB.unique.as.character.miRDB.Gene.Symbol..`, 
#                    inters$`TargetScan.unique.as.character.targetscan.Target.gene..:miRWalk.unique.as.character.miRWalk.genesymbol..`,
#                    inters$`TargetScan.unique.as.character.targetscan.Target.gene..:miRDB.unique.as.character.miRDB.Gene.Symbol..:miRWalk.unique.as.character.miRWalk.genesymbol..`,
#                    inters$`miRDB.unique.as.character.miRDB.Gene.Symbol..:miRWalk.unique.as.character.miRWalk.genesymbol..`)


### Create binary dataframe based on yes/no gene target prediction for each miRNA
miR16.genes.df$miR.16.5p <- 1
miR126.genes.df$miR.126.5p <- 1
merged.df <- merge(miR16.genes.df, miR126.genes.df, by = "V1", all = T)
miR335.genes.df$miR.335.5p <- 1
merged.df <- merge(merged.df, miR335.genes.df, by = "V1", all = T)
merged.df[is.na(merged.df)] <- 0
#rownames(merged.df) <- merged.df$V1
#merged.df <- merged.df %>% select(-c("V1"))

# Cluster
merged.df <- merged.df[rowSums(merged.df == 0) <= 1, ]   # Remove rows with more than two zeros
dat <- merged.df[,2:4]
rownames(dat) <- merged.df[,1]
row.order <- hclust(dist(dat))$order # Clustering
col.order <- hclust(dist(t(dat)))$order
dat_new <- dat[row.order, col.order] #Re-order maxtrix according to clustering
df_molten_dat <- melt(as.matrix(dat_new))
names(df_molten_dat)[c(1:2)] <- c("Genes", "miRNA")

# Plot
pdf(width = 4, height = 15, paste0(args[4], "_Heatmap.pdf"))
ggplot(df_molten_dat, aes(miRNA, Genes))+
  geom_tile(aes(fill = value), colour = "white", size=0.05) +
  scale_fill_gradient(low = "white", high = "darkorchid1") +
  #theme_minimal() +
  labs(x = "", y = "") +
  theme(legend.position = "none", axis.text.y = element_text(size=2)) +
  #theme(legend.position = "none", axis.text.y = element_blank()) +
  ggtitle("Clustered miRNA Gene Target Prediction Heatmap")
dev.off()

