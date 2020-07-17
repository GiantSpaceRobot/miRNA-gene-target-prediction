## ---------------------------
##
## Script name:  miRNA-gene-target-prediction-analysis.
##
## Purpose of script: Given raw TargetScan, miRDB, and miRWalk output, process data and analyse using EnrichR
##
## Author: Dr. Paul Donovan
##
## Date Created: 15-07-20
##
## Email: pauldonovan@rcsi.com
##
## ---------------------------
##
## Notes: 
##   
##
## ---------------------------

library("enrichR")
library("ggplot2")
library("VennDiagram")
library("gplots")
library("stringr")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger") # Suppress VennDiagram messages/log file generation

### Use command line arguments
args = commandArgs(trailingOnly=TRUE)

targetscan <- read.csv(args[1],
                       sep = "\t", 
                       header = T)
miRDB <- read.csv(args[2], 
                  sep = "\t", 
                  header = T)
miRWalk <- read.csv(args[3], 
                    sep = ",", 
                    header = T)


# targetscan <- read.csv("/home/paul/Documents/Acute_Ischaemic_Stroke/miR-126-5p/TargetScan7.2__miR-126-5p.predicted_targets.txt",
#                        sep = "\t", 
#                        header = T)
# miRDB <- read.csv("/home/paul/Documents/Acute_Ischaemic_Stroke/miR-126-5p/miRDB.tsv", 
#                   sep = "\t", 
#                   header = T)
# miRWalk <- read.csv("/home/paul/Documents/Acute_Ischaemic_Stroke/miR-126-5p/miRWalk_miRNA_Targets.csv", 
#                     sep = ",", 
#                     header = T)

my.output <- args[4]

# Get unique genes predicted by each method
my.targetscan.genes.df <- data.frame(unique(as.character(targetscan$Target.gene)))
my.miRDB.genes.df <- data.frame(unique(as.character(miRDB$Gene.Symbol)))
my.miRWalk.genes.df <- data.frame(unique(as.character(miRWalk$genesymbol)))

my.targetscan.genes <- unique(as.character(targetscan$Target.gene))
my.miRDB.genes <- unique(as.character(miRDB$Gene.Symbol))
my.miRWalk.genes <- unique(as.character(miRWalk$genesymbol))

# Venn Diagram
geneLists <- c("TargetScan" = my.targetscan.genes.df,
               "miRDB" = my.miRDB.genes.df,
               "miRWalk" = my.miRWalk.genes.df)

total.no.of.features <- unique(c(my.targetscan.genes, my.miRDB.genes, my.miRWalk.genes))
write.table(x = total.no.of.features, file = paste0(my.output, "_all-predicted-gene-targets.tsv"), 
            sep = "\t", 
            quote = F, 
            row.names = F, 
            col.names = F)

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
                        category.names=c(paste0("TargetScan\n(",length(my.targetscan.genes), ")"), 
                                         paste0("miRDB\n(",length(my.miRDB.genes), ")"), 
                                         paste0("miRWalk (",length(my.miRWalk.genes), ")")),
                        main.cex = 1.4,
                        main=paste0("Total no. of genes predicted = ", 
                                    length(total.no.of.features)))

png(file = paste0(args[4], "_VennDiagram.png"))
#png(file = "/home/paul/My_VennDiagram.png")
grid.newpage() # Create new grid
pushViewport(viewport(width=unit(0.8, "npc"), height = unit(0.8, "npc"))) # Resize this grid (80%)
grid.draw(my.venn) # Draw venn diagram in new smaller plot
dev.off()

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

### Get genes for GO analysis (use all or use ones prediced by two methods)
if (identical(args[5],"yes")) {
  ### Use all genes
  my.gene.list <- total.no.of.features
} else if (identical(args[5],"no")) {
  my.gene.list <- c(inters$`TargetScan.unique.as.character.targetscan.Target.gene..:miRDB.unique.as.character.miRDB.Gene.Symbol..`, 
                    inters$`TargetScan.unique.as.character.targetscan.Target.gene..:miRWalk.unique.as.character.miRWalk.genesymbol..`,
                    inters$`TargetScan.unique.as.character.targetscan.Target.gene..:miRDB.unique.as.character.miRDB.Gene.Symbol..:miRWalk.unique.as.character.miRWalk.genesymbol..`,
                    inters$`miRDB.unique.as.character.miRDB.Gene.Symbol..:miRWalk.unique.as.character.miRWalk.genesymbol..`)
}

### Get all genes predicted by two or more methods

### GO Analysis
dbs <- c("Enrichr_Libraries_Most_Popular_Genes",
         "GO_Molecular_Function_2018", 
         "GO_Cellular_Component_2018", 
         "GO_Biological_Process_2018",
         "WikiPathways_2019_Human",
         "KEGG_2019_Human")
my.enrichR <- enrichr(my.gene.list, databases = dbs)

my.summary <- printEnrich(data = my.enrichR, 
            file = paste0(my.output, "_Summary.tsv"), 
            sep = "\t", 
            columns = c(1:9))

write.table(x = my.enrichR[["Enrichr_Libraries_Most_Popular_Genes"]], file = paste0(my.output, "_Enrichr_Libraries_Most_Popular_Genes.tsv"))
write.table(x = my.enrichR[["GO_Molecular_Function_2018"]], file = paste0(my.output, "_GO_Molecular_Function_2018.tsv"))
write.table(x = my.enrichR[["GO_Cellular_Component_2018"]], file = paste0(my.output, "_GO_Cellular_Component_2018.tsv"))
write.table(x = my.enrichR[["GO_Biological_Process_2018"]], file = paste0(my.output, "_GO_Biological_Process_2018.tsv"))
#write.table(x = my.enrichR[["Chromosome_Location_hg19"]], file = paste0(my.output, "_Chromosome_Location_hg19.tsv"))
#write.table(x = my.enrichR[["OMIM_Disease"]], file = paste0(my.output, "_OMIM_Disease.tsv"))
#write.table(x = my.enrichR[["OMIM_Expanded"]], file = paste0(my.output, "_OMIM_Expanded.tsv"))
#write.table(x = my.enrichR[["ENCODE_Histone_Modifications_2015"]], file = paste0(my.output, "_ENCODE_Histone_Modifications_2015.tsv"))
#write.table(x = my.enrichR[["CORUM"]], file = paste0(my.output, "_CORUM.tsv"))
#write.table(x = my.enrichR[["InterPro_Domains_2019"]], file = paste0(my.output, "_InterPro_Domains_2019.tsv"))
#write.table(x = my.enrichR[["KEA_2015"]], file = paste0(my.output, "_KEA_2015.tsv"))
#write.table(x = my.enrichR[["miRTarBase_2017"]], file = paste0(my.output, "_miRTarBase_2017.tsv"))
#write.table(x = my.enrichR[["NCI-Nature_2016"]], file = paste0(my.output, "_NCI-Nature_2016.tsv"))
#write.table(x = my.enrichR[["NURSA_Human_Endogenous_Complexome"]], file = paste0(my.output, "_NURSA_Human_Endogenous_Complexome.tsv"))
#write.table(x = my.enrichR[["Panther_2016"]], file = paste0(my.output, "_Panther_2016.tsv"))
#write.table(x = my.enrichR[["Pfam_Domains_2019"]], file = paste0(my.output, "_Pfam_Domains_2019.tsv"))
#write.table(x = my.enrichR[["Pfam_InterPro_Domains"]], file = paste0(my.output, "_Pfam_InterPro_Domains.tsv"))
#write.table(x = my.enrichR[["Phosphatase_Substrates_from_DEPOD"]], file = paste0(my.output, "_Phosphatase_Substrates_from_DEPOD.tsv"))
#write.table(x = my.enrichR[["PPI_Hub_Proteins"]], file = paste0(my.output, "_PPI_Hub_Proteins.tsv"))
#write.table(x = my.enrichR[["Rare_Diseases_GeneRIF_Gene_Lists"]], file = paste0(my.output, "_Rare_Diseases_GeneRIF_Gene_Lists.tsv"))
#write.table(x = my.enrichR[["Reactome_2016"]], file = paste0(my.output, "_Reactome_2016.tsv"))
#write.table(x = my.enrichR[["TargetScan_microRNA_2017"]], file = paste0(my.output, "_TargetScan_microRNA_2017.tsv"))
write.table(x = my.enrichR[["WikiPathways_2019_Human"]], file = paste0(my.output, "_WikiPathways_2019_Human.tsv"))
write.table(x = my.enrichR[["KEGG_2019_Human"]], file = paste0(my.output, "_KEGG_2019_Human.tsv"))

df.mpg <- my.enrichR[["Enrichr_Libraries_Most_Popular_Genes"]]
df.mf <- my.enrichR[["GO_Molecular_Function_2018"]]
df.cc <- my.enrichR[["GO_Cellular_Component_2018"]]
df.bp <- my.enrichR[["GO_Biological_Process_2018"]]
#df.chr <- my.enrichR[["Chromosome_Location_hg19"]]
#df.disease <- my.enrichR[["OMIM_Disease"]]
#df7.disease.exp <- my.enrichR[["OMIM_Expanded"]]
#df8.mods <- my.enrichR[["ENCODE_Histone_Modifications_2015"]]
#df9.CORUM <- my.enrichR[["CORUM"]] 
#df10.InterPro <- my.enrichR[["InterPro_Domains_2019"]]
#df11.kea <- my.enrichR[["KEA_2015"]] 
#df12.miRTarBase <- my.enrichR[["miRTarBase_2017"]]
#df13.NCI <- my.enrichR[["NCI-Nature_2016"]]
#df14.NURSA <- my.enrichR[["NURSA_Human_Endogenous_Complexome"]]
#df15.Panther <- my.enrichR[["Panther_2016"]]
#df16.Pfam <- my.enrichR[["Pfam_Domains_2019"]]
#df17.Pfam.InterPro <- my.enrichR[["Pfam_InterPro_Domains"]]
#df18.Phos <- my.enrichR[["Phosphatase_Substrates_from_DEPOD"]]
#df19.PPI <- my.enrichR[["PPI_Hub_Proteins"]]
#df20.RareDisease <- my.enrichR[["Rare_Diseases_GeneRIF_Gene_Lists"]]
#df21.reactome <- my.enrichR[["Reactome_2016"]]
#df22.TargetScan <- my.enrichR[["TargetScan_microRNA_2017"]]
df23.WikiPathways <- my.enrichR[["WikiPathways_2019_Human"]]
df24.KEGG <- my.enrichR[["KEGG_2019_Human"]]

df.to.plot <- function(my.dataframe, analysis = "") {
  #print(paste0("Rows in ", deparse(substitute(my.dataframe)), ": ", nrow(my.dataframe)))
  sub.df <- head(my.dataframe, n = 10)
  sub.df <- sub.df[order(-sub.df$Adjusted.P.value),]
  sub.df$log10.of.padj <- abs(log10(sub.df$Adjusted.P.value))
  sub.df <- subset(sub.df, log10.of.padj>1.3)
  number.of.terms <- nrow(sub.df)
  if(number.of.terms<1){
    return(plot.new())
  } else {
  #text.size <- -log10(number.of.terms) ####
  #text.size <- ((2/number.of.terms)+1.18)^3 ####
  #print(text.size)
  level_order <- as.list(sub.df$Term) # Create factor of terms to reorder ggplot2
  my.plot <- ggplot(data=sub.df, aes(x=factor(Term, level = level_order), y=log10.of.padj)) +
    geom_bar(stat="identity", fill="darkred") +  #steelblue
    #geom_text(
    #  aes(label=Term), 
    #  hjust = 1,
    #  color="black", 
    #  size=text.size) +
    #theme_minimal() +
    theme(axis.text.y=element_text(size=10))+#, 
    #      axis.ticks.y=element_blank(), 
    #      panel.grid.major = element_blank(), 
    #      panel.grid.minor = element_blank()) +
    scale_x_discrete(name = "", 
                     labels=function(x) sub("(", "\n(", x, fixed=TRUE)) + # Add newline to labels with GO numbers
    #scale_x_discrete(name = "", 
    #                 labels=function(x) 
    #                   if(str_count(x, " ")>4){
    #                     space.count <- (str_count(x, " "))/2;
    #                     round(space.count);
    #                     str_replace(x, paste0("(.*?"," ",".*?) "), paste0("\\", space.count))}) + 
    scale_y_continuous(name = "-Log10 of Adjusted P value") +
    labs(title=analysis) +
    coord_flip()
  # my.plot <- 
  #   ggplot(data=sub.df, aes(x=factor(Term, level = level_order), y=log10.of.padj)) +
  #     geom_bar(stat="identity", fill="darkred") +  #steelblue
  #     geom_text(
  #       aes(label=Term), 
  #       hjust = 1,
  #       color="black", 
  #       size=text.size) +
  #     theme_minimal() +
  #     theme(axis.text.y=element_blank(), 
  #           axis.ticks.y=element_blank(), 
  #           panel.grid.major = element_blank(), 
  #           panel.grid.minor = element_blank()) +
  #     scale_x_discrete(name = "") +
  #     scale_y_continuous(name = "Log10 of Adjusted P value") +
  #     labs(title=analysis) +
  #     coord_flip()
  return(my.plot)
  }
}

pdf(height = 5, file = paste0(my.output, "_EnrichR_Barplots.pdf"))
df.to.plot(df.mf, "GO Molecular Function 2018")
df.to.plot(df.cc, "GO Cellular Component 2018")
df.to.plot(df.bp, "GO Biological Process 2018")
df.to.plot(df.mpg, "Enrichr Libraries Most Popular Genes")
#df.to.plot(df.chr, "Chromosome Location hg19")
#df.to.plot(df.disease, "OMIM Disease")
#df.to.plot(df7.disease.exp, "OMIM Expanded")
#df.to.plot(df8.mods, "ENCODE Histone Modifications 2015")
#df.to.plot(df9.CORUM, "CORUM") 
#df.to.plot(df10.InterPro, "InterPro Domains 2019")
#df.to.plot(df11.kea, "KEA 2015") 
#df.to.plot(df12.miRTarBase, "miRTarBase 2017")
#df.to.plot(df13.NCI, "NCI-Nature 2016")
#df.to.plot(df14.NURSA, "NURSA Human Endogenous Complexome")
#df.to.plot(df15.Panther, "Panther 2016")
#df.to.plot(df16.Pfam, "Pfam Domains 2019")
#df.to.plot(df17.Pfam.InterPro, "Pfam InterPro Domains")
#df.to.plot(df18.Phos, "Phosphatase Substrates from DEPOD")
#df.to.plot(df19.PPI, "PPI Hub Proteins")
#df.to.plot(df20.RareDisease, "Rare Diseases GeneRIF Gene Lists")
#df.to.plot(df21.reactome, "Reactome 2016")
#df.to.plot(df22.TargetScan, "TargetScan microRNA 2017")
df.to.plot(df23.WikiPathways, "WikiPathways 2019 Human")
df.to.plot(df24.KEGG, "KEGG 2019 Human")
dev.off()

