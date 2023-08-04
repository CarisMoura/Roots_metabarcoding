#-----------------------------------------------------------------------------------------------------------------------------------------#
#                                                                                                                                         #
#   Codes for: Moura et al. "Insights on root-fungal association in forest conversion systems through  network and indicator analysis"    # 
#                                                                                                                                         #
#-----------------------------------------------------------------------------------------------------------------------------------------#
#Contents: 
#1. Figure S1 - Accumulation curve using phyloseq object and gamma diversity barplot
#2. Figure S2 - Radar plot Chao1 and Shannon alpha diversity
#3. Table S1  - PERMANOVA
#4. Figure S3 - Weighted UPGMA tree based on Unifrac distances and barplot of taxa composition 
#5. Figure 1  - Bipartite network based on indicator species
#6. Table 1   - Network metrics of roots and fungi co-occurrence network
#7. Table S2  - Network metric results


#------------------------------------------------------------------------------#
#                                 Packages                                     #
#------------------------------------------------------------------------------#
library(readxl)
library(phyloseq)
library(limma)
library(dplyr)
library(ape)
library(vegan)
library(ggplot2)
library(BiodiversityR) # also loads vegan
library(ggsci)
library(hrbrthemes)
library(fmsb)
library(data.table)
library(ggradar)
#library(palmerpenguins)
library(tidyverse)
library(scales)
library(showtext)
library(ggfortify)
library(dendextend)
library(ggpubr)
library(qiime2R)

#------------------------------------------------------------------------------#
#                         Data Import and Preparation                          #
#------------------------------------------------------------------------------#

#-------------------------#
#     Plant data set      #
#-------------------------#
##### Load OTU table for root sequences #####
roots.df <- "C:\\Users\\carin\\Documents\\Network_R\\Table_roots_revised_CM.xlsx"
dim(roots.df)
roots.otux <- as.data.frame(readxl::read_excel(roots.df, sheet=1))
roots.tax <- as.data.frame(readxl::read_excel(roots.df, sheet=2))
roots.env <- as.data.frame(readxl::read_excel(roots.df, sheet=3))

dim(roots.otu1)


##### Load OTU table for root sequences #####
roots.otu2 <- roots.otux%>%
  tibble::column_to_rownames("OTU") 
roots.tax <- roots.tax %>% 
  tibble::column_to_rownames("OTU")
roots.env <- roots.env %>% 
  tibble::column_to_rownames("Sample") 


##### Transform into matrix #####
roots.otux <- as.matrix(roots.otu2)
roots.tax <- as.matrix(roots.tax)


##### Transform to phyloseq objects for rbcL #####
OTU_roots_prefilterx = otu_table(roots.otux, taxa_are_rows = TRUE)
dim(roots.env)
TAX_roots = tax_table(roots.tax)
samples_roots = sample_data(roots.env)

roots_ob_prefilter <- phyloseq(OTU_roots_prefilterx, TAX_roots, samples_roots)


##### Filter less 1% of reads per sample #####
roots.otu.filtered <-filter_taxa(roots_ob_prefilter, function(x) sum(x > 3) > (0.01*length(x)), TRUE)   


##### Add tree to phyloseq object #####
random_tree = rtree(ntaxa(roots.otu.filtered), rooted=TRUE, tip.label=taxa_names(roots.otu.filtered))
plot(random_tree)
otu_roots_filt <- otu_table(roots.otu.filtered)
physeq2 = phyloseq(otu_roots_filt, TAX_roots, samples_roots, random_tree)
plot_tree(physeq2, color="type", label.tips="Sample", ladderize="left", plot.margin=0.3)


##### Rarefy OTU table #####
roots.rarefied = rarefy_even_depth(physeq2, rngseed=1, sample.size = min(sample_sums(physeq2)),
                                   replace = F)
otu.rare <- otu_table(roots.rarefied)
otu_table(physeq2)


##### Estimate richness #####
estimate_richness(roots.rarefied)


##### Load plant dataset for accumulation curve figure #####
##This dataframe is the OTU table, with plots as lines (but plot ID is omitted) and species as columns, NA=0
rbclacc.1 <- as.data.frame(otu.rare)
dim(rbclacc.1)
rbclacc.1 <- t(rbclacc.1)


##This dataframe is the envionmental data, with plots as lines and colummn as variables
rbclacc.env <- as.data.frame(roots.env)
rbclacc.env 
rbclacc.env$site.total <- apply(rbclacc.1, 1, sum)


##### Calculate gamma diversity #####
gamma <- with(roots.env, specnumber(rbclacc.1, type))


##### Test if richness is significantly different #####
richness_roots <- read_excel("C:\\Users\\carin\\Documents\\Network_R\\Biodiv_roots_v1.xlsx")
richness_roots <- as.data.frame(richness_roots)

res.aov_rt <- aov(Shannon ~ Landuse, data= richness_roots)
summary(res.aov_rt)
TukeyHSD(res.aov_rt)


#-------------------------#
#     Fungi data set      #
#-------------------------#
##Load Fungi dataset 
##### Phyloseq object for fungi #####
fungi.otu<- read_excel("C:\\Users\\carin\\Documents\\Network_R\\Table_fungi_revised_CM.xlsx", 
                       sheet = 1)
fungi.tax<- read_excel("C:\\Users\\carin\\Documents\\Network_R\\Table_fungi_revised_CM.xlsx", 
                       sheet = 2)
fungi.env<- read_excel("C:\\Users\\carin\\Documents\\Network_R\\Table_fungi_revised_CM.xlsx", 
                       sheet = 3)


##### Define row names from the otu column #####
fungi.otu<- fungi.otu %>%
  tibble::column_to_rownames("OTU") 
fungi.tax <- fungi.tax %>% 
  tibble::column_to_rownames("OTU")
fungi.env <- fungi.env %>% 
  tibble::column_to_rownames("Sample") 

dim(fungi.otu)


##### Transform into matrixes otu and tax tables (sample table can be left as data frame) #####
fungi.otu1 <- as.matrix(fungi.otu)
fungi.tax <- as.matrix(fungi.tax)

dim(fungi.otu)


##### Transform to phyloseq objects #####
OTU_fungi = otu_table(fungi.otu1, taxa_are_rows = TRUE)
TAX_fungi = tax_table(fungi.tax)
SAMP_fungi = sample_data(fungi.env)

FungiOB_prefiltered <- phyloseq(OTU_fungi, TAX_fungi, SAMP_fungi)
FungiOB_prefiltered


##### Filter less 1% of reads per sample #####
fungi.otu.filtered <-filter_taxa(FungiOB_prefiltered, function(x) sum(x > 3) > (0.01*length(x)), TRUE)   


##### Add tree to phyloseq object #####
random_tree_fg = rtree(ntaxa(fungi.otu.filtered), rooted=TRUE, tip.label=taxa_names(fungi.otu.filtered))
plot(random_tree_fg)
otu_fg_filt <- otu_table(fungi.otu.filtered)
physeq_fg = phyloseq(otu_fg_filt, TAX_fungi, SAMP_fungi, random_tree_fg)
plot_tree(physeq_fg, color="type", label.tips="Sample", ladderize="left", plot.margin=0.3)


##### Rarefy OTU table #####
fg.rarefied = rarefy_even_depth(physeq_fg, rngseed=1, sample.size = min(sample_sums(physeq_fg)),
                                replace = F)


##### Estimate richness #####
estimate_richness(fg.rarefied)


##### Test if richness is significantly different #####
richness_fungi <- read_excel("C:\\Users\\carin\\Documents\\Network_R\\Biodiv_fungi.xlsx")
richness_fungi <- as.data.frame(richness_fungi)

res.aov <- aov(Shannon ~ Landuse, data= richness_fungi)
summary(res.aov)
TukeyHSD(res.aov)


##### Save otu table within an object #####
otu_fg <- otu_table(fg.rarefied)
fungi.otu.norm <- as.matrix(otu_fg)


##This dataframe is the OTU table, with plots as lines (but plot ID is omitted) and species as columns, NA=0
fungiacc.1 <- as.data.frame(fungi.otu.norm)
fungiacc.1

head(fungiacc.1)


##This dataframe is the envionmental data, with plots as lines and colummn as variables
fungiacc.1 <- t(fungiacc.1)
fungiacc.env <- as.data.frame(fungi.env)
dim(fungiacc.1)
fungiacc.env$site.total <- apply(fungiacc.1, 1, sum)


##### Calculate gamma diversity #####
gamma_fg <- with(fungi.env, specnumber(fungiacc.1, type))


#----------------------------------------------------#
# 1. Figure S1 Accumulation curve and OTU richness   #
#----------------------------------------------------#

## obs. Sometimes functions accumcomp or accumcomp.long fails to load
## functions can be extracted from here https://rdrr.io/cran/BiodiversityR/src/R/sites.long.R
##Run accumulation curve of rbcL with number of reads as x axis

##### Plot accumulation curves per landuse type (factor= type) and Add grid to the plot #####
Accum.2 <- accumcomp(rbclacc.1, y=rbclacc.env, factor='type', method='exact', legend=FALSE, conditioned = TRUE, scale= 'site.total')
Accum.2
grid(nx = NULL, ny= NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)

Accum.3 <- accumcomp(rbclacc.1, y=rbclacc.env, factor='type', method='collector', legend=FALSE, conditioned = TRUE, scale= 'site.total')
Accum.3
grid(nx = NULL, ny= NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)


##### Plot accumulation curves per plot (scale= site.total) and Add grid to the plot #####
Accum.rbcldf <- accumresult(rbclacc.1, y=rbclacc.env, scale='site.total', method='exact', conditioned = TRUE, ylim= c(0,150))
Accum.rbcldf
accumplot(Accum.rbcldf)

## accumcomp calculates the species accumulation curve for all levels of a selected environmental variable separatedly.
Accum.rbcldf1 <- accumcomp(rbclacc.1, y=rbclacc.env, factor='type', level= 'type', method='collector', legend=FALSE, conditioned = TRUE, ylim= c(0,200))
Accum.rbcldf1
grid(nx =NULL, ny= NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) 


Accum.rbcldf2 <- accumcomp(rbclacc.1, y=rbclacc.env, factor='type', method='collector', 
                           legend=FALSE, conditioned = TRUE, scale= 'site.total')


##### Run accumulation curve of fungi with number of reads as x axis #####
fungiacc.1 <- t(fungiacc.1)
fungiacc.env$site.total <- apply(fungiacc.1, 1, sum)


##### Plot accumulation curves per landuse type (factor= type) and Add grid to the plot #####
fungiAccum.2 <- accumcomp(fungiacc.1, y=fungiacc.env, factor='type', method='exact', legend=FALSE, conditioned = TRUE, scale= 'site.total')
fungiAccum.2
grid(nx = NULL, ny= NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)


##### FigureS1A - Plot accumulation curves using ggplot for roots #####
accum.long1 <- accumcomp.long(Accum.2, ci=NA, label.freq=1)
head(accum.long1)


plot_accum_roots <- ggplot(data=accum.long1, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=1) +
  geom_point(data=subset(accum.long1, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=3) + 
  theme_bw() +
  scale_color_manual(values = c("#006a00", "#7e9828","#ffdb6d", "#f7a436")) +
  labs(x = expression("Sequence sample size"), y = "Plant OTU richness", colour = "Grouping", shape = "Grouping")

plot_accum_roots <- plot_accum_roots + theme_light()
FigureS1A <- plot_accum_roots

##### FigureS1B - Accumulation curve for plant OTUs detected in roots using plots as x-axis #####
accum.long_plot <- accumcomp.long(Accum.rbcldf1, ci=NA, label.freq=1)

plot_accum_roots_plot <- ggplot(data=accum.long_plot, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=1) +
  geom_point(data=subset(accum.long_plot, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=3) + 
  theme_bw() +
  scale_color_manual(values = c("#006a00", "#7e9828","#ffdb6d", "#f7a436")) +
  labs(x = expression("Sample size"), y = "Plant OTU richness", colour = "Grouping", shape = "Grouping")

plot_accum_roots_plot <- plot_accum_roots_plot + theme_light()
FigureS1B <- plot_accum_roots_plot +  scale_x_continuous(breaks=seq(0, 10, 2))


##### FigureS1D - Plot accumulation curves using ggplot for fungi #####
accum.long_fg <- accumcomp.long(fungiAccum.2, ci=NA, label.freq=1)   ##label.freq indicates the frequence of plotting sample points to the accumulation curve
head(accum.long_fg)

plot_accum_fungi <- ggplot(data=accum.long_fg, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=1) +
  geom_point(data=subset(accum.long_fg, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=3) + 
  scale_color_manual(values = c("#006a00", "#7e9828","#ffdb6d", "#f7a436")) +
  labs(x = expression("Sequence sample size"), y = "Fungi OTU richness", colour = "Grouping", shape = "Grouping")

plot_accum_fungi1 <- plot_accum_fungi+ theme_light()
FigureS1D <- plot_accum_fungi1


##### FigureS1E - Plot accumulation curves using plots as x-axis #####
Accum.fg1 <- accumcomp(fungiacc.1, y=fungiacc.env, factor='type', level= 'type', method='collector', legend=FALSE, conditioned = TRUE, ylim= c(0,1200))
Accum.fg1

accum.long_plot_fg <- accumcomp.long(Accum.fg1, ci=NA, label.freq=1)

plot_accum_fg_plot <- ggplot(data=accum.long_plot_fg, aes(x = Sites, y = Richness)) + 
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  geom_line(aes(colour=Grouping), size=1) +
  geom_point(data=subset(accum.long_plot_fg, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=3) + 
  theme_bw() +
  scale_color_manual(values = c("#006a00", "#7e9828","#ffdb6d", "#f7a436")) +
  labs(x = expression("Sample size"), y = "Fungi OTU richness", colour = "Grouping", shape = "Grouping")

plot_accum_fg_plot <- plot_accum_fg_plot + theme_light()
plot_accum_fg_plot +  scale_x_continuous(breaks=seq(0, 10, 2))

FigureS1E <- plot_accum_fg_plot +  scale_x_continuous(breaks=seq(0, 10, 2))


##### Add rarified table to an object for plant sequences #####
x1 <- t(otu_table(roots.rarefied))
Accum.s1 <- accumcomp(x1, y=rbclacc.env, factor='type', method='collector', legend=FALSE, conditioned = TRUE, scale= 'site.total')
Accum.s1
grid(nx = NULL, ny= NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)

accum.longx <- accumcomp.long(Accum.s1, ci=NA, label.freq=4)

##### FigureS1C - Plot bar and standard deviation #####
my_sum <- accum.long1 %>%
  group_by(Grouping) %>%
  summarise( 
    n=n(),
    mean=mean(Richness),
    sd=sd(Richness)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

level_order <- c('Forest', 'Jungle', 'Rubber', 'OilPalm')
cores <- c("#006a00", "#7e9828", "#ffdb6d", "#f7a436")

p2 <- ggplot(my_sum) +
  geom_bar( aes(x= factor(Grouping, level= level_order), y=mean), stat="identity", fill=cores, alpha=0.7, width=0.6) + 
  geom_errorbar( aes(x=Grouping, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.9, size=1) +
  scale_color_manual("legend", values = c("#006a00", "#7e9828","#ffdb6d", "#f7a436"))+
  xlab("") +
  ylab("Plant OTU richness")

rt_richness_plot <- p2 +   theme_light()
FigureS1C <- rt_richness_plot


##### Add rarified table to an object for fungi sequences ######
f1 <- t(otu_table(fg.rarefied))
Accum.f1 <- accumcomp(f1, y=fungiacc.env, factor='type', method='collector', legend=FALSE, conditioned = TRUE, scale= 'site.total')
Accum.f1
grid(nx = NULL, ny= NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)

accum.longfgx <- accumcomp.long(Accum.f1, ci=NA, label.freq=4)

##### FigureS1F - Plot bar and standard deviation #####
my_sum_fg <- accum.longfgx %>%
  group_by(Grouping) %>%
  summarise( 
    n=n(),
    mean=mean(Richness),
    sd=sd(Richness)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

level_order <- c('Forest', 'Jungle', 'Rubber', 'OilPalm')
cores <- c("#006a00", "#7e9828", "#ffdb6d", "#f7a436")

f2 <- ggplot(my_sum_fg) +
  geom_bar( aes(x= factor(Grouping, level= level_order), y=mean), stat="identity", fill=cores, alpha=0.7, width=0.6) + 
  geom_errorbar( aes(x=Grouping, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="black", alpha=0.9, size=1) +
  scale_color_manual("legend", values = c("#006a00", "#7e9828","#ffdb6d", "#f7a436"))+
  xlab("") +
  ylab("Fungi OTU richness")


fg_richness_plot <- f2 +   theme_light()
FigureS1F <- fg_richness_plot


##### Place plots together #####
library(cowplot)
plot_grid(FigureS1A, FigureS1B, FigureS1C, FigureS1D, FigureS1E, FigureS1F, labels = "AUTO", align = "h")


#----------------------------------------------------------#
# 2. Figure S2 Diversity of functional groups and guilds   #
#----------------------------------------------------------#
##### Radar plot based on Shannon diversity of functional groups and guilds #####

##### Subset taxa #####
AMF <- subset_taxa(fg.rarefied, groups_new== "AM")
ECM <- subset_taxa(fg.rarefied, groups== "ECM")
Symbiotroph <- subset_taxa(fg.rarefied, Guild...9== "Symbiotroph")
Pathotroph <- subset_taxa(fg.rarefied, Guild...9== "Pathotroph")
Saprotroph <- subset_taxa(fg.rarefied, Guild...9== "Saprotroph")
Trees_p <- subset_taxa(roots.rarefied, Category== "Tree")
Shrub_p <- subset_taxa(roots.rarefied, Category== "Shrub")
Herb_p <- subset_taxa(roots.rarefied, Category== "Herb")
Alien_p <- subset_taxa(roots.rarefied, Native_Alien== "Non-native")
Native_p <- subset_taxa(roots.rarefied, Native_Alien== "Native")


##### Estimate richness for subset of functional groups and guilds #####
AMF_rich <- estimate_richness(AMF)
AMF_rich$Landuse <- roots.env$type
AMF_rich_gp <- aggregate(AMF_rich$Shannon, list(AMF_rich$Landuse), FUN=mean)

ECM_rich <- estimate_richness(ECM, measure= "Shannon")
ECM_rich$Landuse <- roots.env$type
ECM_rich_gp <- aggregate(ECM_rich$Shannon, list(ECM_rich$Landuse), FUN=mean)

Symbiotroph_rich <- estimate_richness(Symbiotroph, measure= "Shannon")
Symbiotroph_rich$Landuse <- roots.env$type
Symb_rich_gp <- aggregate(Symbiotroph_rich$Shannon, list(Symbiotroph_rich$Landuse), FUN=mean)

Pathotroph_rich <- estimate_richness(Pathotroph, measure= "Shannon")
Pathotroph_rich$Landuse <- roots.env$type
Patho_rich_gp <- aggregate(Pathotroph_rich$Shannon, list(Pathotroph_rich$Landuse), FUN=mean)

Saprotroph_rich <- estimate_richness(Saprotroph, measure= "Shannon")
Saprotroph_rich$Landuse <- roots.env$type
Sapro_rich_gp <- aggregate(Saprotroph_rich$Shannon, list(Saprotroph_rich$Landuse), FUN=mean)

Trees_rich <- estimate_richness(Trees_p, measure= "Shannon")
Trees_rich$Landuse <- roots.env$type
Trees_rich_gp <- aggregate(Trees_rich$Shannon, list(Trees_rich$Landuse), FUN=mean)

Shrub_rich <- estimate_richness(Shrub_p, measure= "Shannon")
Shrub_rich$Landuse <- roots.env$type
Shrub_rich_gp <- aggregate(Shrub_rich$Shannon, list(Shrub_rich$Landuse), FUN=mean)

Herb_rich <- estimate_richness(Herb_p, measure= "Shannon")
Herb_rich$Landuse <- roots.env$type
Herb_rich_gp <- aggregate(Herb_rich$Shannon, list(Herb_rich$Landuse), FUN=mean)

Alien_rich <- estimate_richness(Alien_p, measure= "Shannon")
Alien_rich$Landuse <- roots.env$type
Alien_rich_gp <- aggregate(Alien_rich$Shannon, list(Alien_rich$Landuse), FUN=mean)

Nat_rich <- estimate_richness(Native_p, measure= "Shannon")
Nat_rich$Landuse <- roots.env$type
Nati_rich_gp <- aggregate(Nat_rich$Shannon, list(Nat_rich$Landuse), FUN=mean)


##### Group diversity in one data frame #####
Shannon_scores <- data.frame(
  row.names= c("Rainforest", "Jungle", "Rubber", "Oilpalm"),
  AMF_rich_gp, ECM_rich_gp, Symb_rich_gp, Patho_rich_gp, Sapro_rich_gp,Trees_rich_gp,
  Shrub_rich_gp, Herb_rich_gp, Alien_rich_gp, Nati_rich_gp)


##### Remove unnecessary columns #####
df_radar = Shannon_scores[ -c(1,3, 5, 7, 9, 11, 13, 15, 17, 19) ]


##### Set min and max #####
max_min1 <- data.frame(
  x = c(3, 0), x.1 = c(3, 0), x.2 = c(3, 0), x.3 = c(3, 0), x.4 = c(3, 0), x.5 = c(3, 0), x.6 = c(3,0), x.7 = c(3,0),
  x.8 = c(3,0), x.9 = c(3,0))
rownames(max_min1) <- c("Max", "Min")


df_rad <- rbind(max_min1, df_radar)


##### Rename columns #####
setnames(df_rad, old = c('x','x.1','x.2', 'x.3', 'x.4', 'x.5', 'x.6', 'x.7', 'x.8', 'x.9'), 
         new = c('AMF','ECM','Symbiotroph','Pathotroph', ' Saprotroph', 'Tree', 'Shrub', 'Herb', 'Nonnative', 'Native'))


##### Subset only rainforest data #####
data <- df_rad[c("Max", "Min", "Rainforest"), ]


##### Plot radar plot #####
create_beautiful_radarchart <- function(data, color = "#006a00", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

op <- par(mar = c(1, 2, 2, 1))
create_beautiful_radarchart(data, caxislabels = c(0, 1, 2, 3, 4))
par(op)


##### Reduce plot margin using par() ######
op <- par(mar = c(1, 2, 2, 2))
dev.off()


##### Create the radar charts ######
create_beautiful_radarchart(
  data = df_rad, caxislabels = c(0, 1, 2 ,3),
  color = c("#006a00", "#7e9828", "#f7a436", "#ffdb6d"),
  cex = 1
) 


##### Add an horizontal legend #####
legend(
  x = "right", legend = rownames(df[-c(1,2),]), horiz = FALSE,
  bty = "n", pch = 20 , col = c("#006a00", "#7e9828", "#f7a436", "#ffdb6d"),
  text.col = "black", cex = 1, pt.cex = 1.5
)
par(op)


##### Define colors and titles #####
colors <- c("#006a00", "#7e9828", "#f7a436", "#ffdb6d")
titles <- c("Rainforest", "Jungle", "Rubber", "Oilpalm")


##### Reduce plot margin using par() #####
##### Split the screen in 4 parts #####
op <- par(mar = c(1, 1, 1, 1))
par(mfrow = c(1,4))


##### Create the radar chart #####
for(i in 1:4){
  create_beautiful_radarchart(
    data = df_rad[c(1, 2, i+2), ], caxislabels = c(0, 1, 2, 3, 4),
    color = colors[i], title = titles[i]
  )
}
par(op)


##### Estimate richness for subset of functional groups and guilds - focus on Chao1 #####
AMF_rich <- estimate_richness(AMF)
AMF_rich$Landuse <- roots.env$type
AMF_chao_gp <- aggregate(AMF_rich$Chao1, list(AMF_rich$Landuse), FUN=mean)

ECM_chao <- estimate_richness(ECM, measure= "Chao1")
ECM_chao$Landuse <- roots.env$type
ECM_chao_gp <- aggregate(ECM_chao$Chao1, list(ECM_chao$Landuse), FUN=mean)

Symbiotroph_chao <- estimate_richness(Symbiotroph, measure= "Chao1")
Symbiotroph_chao$Landuse <- roots.env$type
Symb_chao_gp <- aggregate(Symbiotroph_chao$Chao1, list(Symbiotroph_chao$Landuse), FUN=mean)

Pathotroph_chao <- estimate_richness(Pathotroph, measure= "Chao1")
Pathotroph_chao$Landuse <- roots.env$type
Patho_chao_gp <- aggregate(Pathotroph_chao$Chao1, list(Pathotroph_chao$Landuse), FUN=mean)

Saprotroph_chao <- estimate_richness(Saprotroph, measure= "Chao1")
Saprotroph_chao$Landuse <- roots.env$type
Sapro_chao_gp <- aggregate(Saprotroph_chao$Chao1, list(Saprotroph_chao$Landuse), FUN=mean)

Trees_chao <- estimate_richness(Trees_p, measure= "Chao1")
Trees_chao$Landuse <- roots.env$type
Trees_chao_gp <- aggregate(Trees_chao$Chao1, list(Trees_chao$Landuse), FUN=mean)

Shrub_chao <- estimate_richness(Shrub_p, measure= "Chao1")
Shrub_chao$Landuse <- roots.env$type
Shrub_chao_gp <- aggregate(Shrub_chao$Chao1, list(Shrub_chao$Landuse), FUN=mean)

Herb_chao <- estimate_richness(Herb_p, measure= "Chao1")
Herb_chao$Landuse <- roots.env$type
Herb_chao_gp <- aggregate(Herb_chao$Chao1, list(Herb_chao$Landuse), FUN=mean)

Alien_chao <- estimate_richness(Alien_p, measure= "Chao1")
Alien_chao$Landuse <- roots.env$type
Alien_chao_gp <- aggregate(Alien_chao$Chao1, list(Alien_chao$Landuse), FUN=mean)

Nat_chao <- estimate_richness(Native_p, measure= "Chao1")
Nat_chao$Landuse <- roots.env$type
Nati_chao_gp <- aggregate(Nat_chao$Chao1, list(Nat_chao$Landuse), FUN=mean)


##### Group diversity in one data frame #####
Shannon_scores_chao <- data.frame(
  row.names= c("Rainforest", "Jungle", "Rubber", "Oilpalm"),
  AMF_chao_gp, ECM_chao_gp, Symb_chao_gp, Patho_chao_gp, Sapro_chao_gp,Trees_chao_gp,
  Shrub_chao_gp, Herb_chao_gp, Alien_chao_gp, Nati_chao_gp)


##### Remove unnecessary columns #####
df_radar_chao = Shannon_scores_chao[ -c(1,3, 5, 7, 9, 11, 13, 15, 17, 19) ]


##### Set min and max #####
max_minx <- data.frame(
  x = c(100, 0), x.1 = c(100, 0), x.2 = c(100, 0), x.3 = c(100, 0), x.4 = c(100, 0), x.5 = c(100, 0), x.6 = c(100,0), x.7 = c(100,0),
  x.8 = c(100,0), x.9 = c(100,0))
rownames(max_minx) <- c("Max", "Min")


df_rad_chao <- rbind(max_minx, df_radar_chao)


##### Rename columns #####
setnames(df_rad_chao, old = c('x','x.1','x.2', 'x.3', 'x.4', 'x.5', 'x.6', 'x.7', 'x.8', 'x.9'), 
         new = c('AMF','ECM','Symbiotroph','Pathotroph', ' Saprotroph', 'Tree', 'Shrub', 'Herb', 'Nonnative', 'Native'))


##### Subset only rainforest data for quick visualization #####
data_chao <- df_rad_chao[c("Max", "Min", "Rainforest"), ]


##### Plot radar plot #####
create_beautiful_radarchart1 <- function(data_chao, color = "#006a00", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

op <- par(mar = c(1, 2, 2, 1))
create_beautiful_radarchart(data_chao, caxislabels = c(0, 25,50, 75, 100))
par(op)


##### Reduce plot margin using par() #####
op <- par(mar = c(1, 2, 2, 2))

dev.off()


##### Create the radar charts #####
radar_chao <- create_beautiful_radarchart(
  data = df_rad_chao, caxislabels = c(0, 25,50,75, 100),
  color = c("#006a00", "#7e9828", "#f7a436", "#ffdb6d"),
  cex = 1
) 


##### Add legend #####
legend(
  x = "right", legend = rownames(df[-c(1,2),]), horiz = FALSE,
  bty = "n", pch = 20 , col = c("#006a00", "#7e9828", "#f7a436", "#ffdb6d"),
  text.col = "black", cex = 1, pt.cex = 1.5
)
par(op)


#----------------------------------------------------#
# 3. Table S1 PERMANOVA                              #
#----------------------------------------------------#
##### Estimate bray curtis distance for plant roots data #####
bray_d <- distance(roots.rarefied, "bray") 


##### Estimate PERMANOVA for plant roots dataset #####
wunifrac_dist = phyloseq::distance(roots.rarefied, method="unifrac", weighted=F)
adonis2(wunifrac_dist ~ sample_data(roots.rarefied)$type)
rtbio.disper <- betadisper(bray_d, sample_data(roots.rarefied)$type)
permutest(rtbio.disper, pairwise = TRUE)


##### Estimate bray curtis distance for fungi data #####
bray_d_fg <- distance(fg.rarefied, "bray") 


##### Estimate PERMANOVA for fungi dataset #####
wunifrac_dist_fg = phyloseq::distance(fg.rarefied, method="unifrac", weighted=F)
adonis2(wunifrac_dist ~ sample_data(fg.rarefied)$type)
fgbio.disper <- betadisper(bray_d_fg, sample_data(fg.rarefied)$type)
permutest(fgbio.disper, pairwise = TRUE)


#----------------------------------------------------#
# 4. Figure S3 - Dendrogram and barplot              #
#----------------------------------------------------#

##### Estimate weightened unifrac distance for plants detected in roots #####
rtUnifrac <- UniFrac(roots.rarefied,  weighted=TRUE)
biod.hclust <- hclust(rtUnifrac, method="ward")
plot(biod.hclust)
plot(hclust(rtUnifrac, method="ward")) 


##### Convert clust to dendrogram #####
rt_dend <- as.dendrogram(biod.hclust)


##### Add symbols to tips #####
nodePar <- list(lab.cex = 0.8, pch = c(NA, 19), 
                cex = 1.2, col = "blue")


##### Plot #####
plot(rt_dend, ylab = "Height", nodePar = nodePar, leaflab = "perpendicular")


##### Specify different point types and colors for each leave #####
cores <- c("#7e9828", "#f7a436", "#7e9828", "#7e9828", "#006a00",
           "#f7a436", "#f7a436", "#f7a436", "#f7a436", "#f7a436", "#f7a436",
           "#ffdb6d", "#7e9828", "#ffdb6d", "#ffdb6d", "#ffdb6d", "#ffdb6d", "#ffdb6d", "#ffdb6d",
           "#7e9828", "#006a00", "#7e9828", "#7e9828", "#006a00", "#006a00", "#006a00",
           "#006a00", "#006a00",  "#7e9828", "#006a00", "#006a00")


cores_ward <- c("#ffdb6d", "#ffdb6d", "#ffdb6d", "#ffdb6d", "#ffdb6d","#ffdb6d", "#ffdb6d", 
                "#f7a436", "#f7a436", "#f7a436", "#f7a436", "#f7a436", "#f7a436", "#7e9828",
                "#7e9828", "#f7a436", "#006a00", "#006a00", "#7e9828", "#006a00", "#006a00",
                "#006a00", "#006a00", "#7e9828", "#7e9828", "#7e9828", "#7e9828", "#7e9828",
                "#006a00", "#006a00", "#006a00")


rt_dend %>% set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 1.5) %>%  # node point size
  set("leaves_col", cores_ward) %>%
  set("labels_cex", 1)%>% #node point color
  ladderize %>% 
  sort(type= "nodes", isReverse = TRUE)%>%
  ladderize %>% 
  sort() %>% 
   plot(main = "Plants - cluster dendrogram")

p_dend <- recordPlot()


##### Plot 10 most abundant orders in the dataset #####
##### Set the order of the dendrogram and barplot #####
top10 <- names(sort(taxa_sums(roots.rarefied), decreasing=TRUE)[1:20])

dat.aglo1 = tax_glom(roots.rarefied, taxrank = "Order")
dat.trans2 = transform_sample_counts(dat.aglo1, function(x) x/sum(x))
prune.dat.two1 = prune_taxa(top10, dat.trans2)
dat.dataframe1 = psmelt(prune.dat.two1)
dat.agr1 = aggregate(Abundance~Sample+type+Order, data=dat.dataframe1, FUN=mean)



rt_bar_order <- ggplot(dat.agr1, aes(x=Sample, y=Abundance, fill=Order)) + geom_bar(stat="identity", position="fill") + scale_fill_manual(values = c("#E9E681", "#F4D03F", "#D0AD32", "#967E2C", "#76D7C4", 
                                                                                                                                           "#45B39D", "#52BE80", "#229954", "#43856C", "#186A3B", 
                                                                                                                                           "#F4C9E0", "#E6B0AA", "#EC7063", "#AA6D66", "#8F4239", 
                                                                                                                                           "#C6AED7", "#8E74A1")) + theme(legend.position="bottom", legend.key.size = unit(0.4, 'cm'))

##### Plot 10 most abundant Genus in the dataset #####
##### Set the order of the dendrogram and barplot #####
top10 <- names(sort(taxa_sums(roots.rarefied), decreasing=TRUE)[1:10])

dat.aglo_genus = tax_glom(roots.rarefied, taxrank = "Genus")
dat.trans_genus = transform_sample_counts(dat.aglo_genus, function(x) x/sum(x))
prune.dat.two_genus = prune_taxa(top10, dat.trans_genus)
dat.dataframe_genus = psmelt(prune.dat.two_genus)
dat.agr_genus = aggregate(Abundance~Sample+type+Genus, data=dat.dataframe_genus, FUN=mean)



rt_bar_genus <- ggplot(dat.agr_genus, aes(x=Sample, y=Abundance, fill=Genus)) + geom_bar(stat="identity", position="fill") + scale_fill_manual(values = c("#E9E681", "#F4D03F", "#D0AD32", "#967E2C", "#76D7C4", 
                                                                                                                                                     "#45B39D", "#52BE80", "#229954", "#43856C", "#186A3B", 
                                                                                                                                                     "#F4C9E0", "#E6B0AA", "#EC7063", "#AA6D66", "#8F4239", 
                                                                                                                                                     "#C6AED7", "#8E74A1")) + theme(legend.position="bottom", legend.key.size = unit(0.4, 'cm'))

sample_order_Rt <- c( "BF2b", "BF4b", "BF3c","BJ4b","BJ2b", "HJ4b", "HJ3b", "HJ2b", "BF1b", "HF1b", "HF3c", "BJ5b",
                      "HF4b", "HF2b", "HF3b", "BJ3c", "HR1b", "BR1b", "BR4b", "HR3b", "HR4b", "BR2b", "HJ1b", "HR2b", 
                      "BO3b", "HO2b", "BO5b", "BO2b", "HO3c", "HO1c", "HO4b")
    
sample_order_Rt <- sample_names(roots.rarefied)[order.dendrogram(rt_dend)]


##### Plot barplot of plant orders composition in the dataset #####                                          
data_rt <- psmelt(roots.rarefied)
p <- ggplot(data=data_rt, aes(x=Sample, y=Abundance, fill=Order))
barplot_roots <- p + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#E9E681", "#F4D03F", "#D0AD32", "#967E2C", "#6b451fff", 
                                 "#82edd8ff", "#45B39D", "#52BE80", "#428345ff", "#13522eff", 
                                 "#F4C9E0", "#E6B0AA", "#7e4c4699", "#EC7063", "#a16159ff", 
                                 "#d6e6f5ff", "#bad2ebff", "#8ebedaff", "#357fb9ff", "#0a3d8aff", 
                                 "#ede5e5ff", "#d3d3d3ff", "#8c7c7c9a", "#6e6f6fff", "#4e4e4eff")) +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.key.size = unit(0.4, 'cm')) + guides(fill=guide_legend(nrow=3)) +
  scale_x_discrete(limits=sample_order_Rt) 
  

barplot_roots <- barplot_roots + theme_bw()
barplot_roots <- barplot_roots + theme(legend.position="bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(fill=guide_legend(nrow=4)) +
  scale_x_discrete(limits=sample_order_Rt) 

bottom_row <- plot_grid(p_dend, barplot_roots, ncol = 1, labels = "hv", align= "l")

plot_grid(p_dend, barplot_roots,
          nrow = 2,
          labels = "AUTO",
          label_size = 12,
          align = "v", axis = "rl"
)


##### Estimate weightened unifrac distance for fungi detected in roots #####
fgUnifrac <- UniFrac(fg.rarefied,  weighted=TRUE)
biod.hclust_fg <- hclust(fgUnifrac, method="ward")
plot(biod.hclust_fg)


##### Convert clust to dendrogram #####
fg_dend <- as.dendrogram(biod.hclust_fg)


##### Add symbols to tips #####
nodePar <- list(lab.cex = 0.8, pch = c(NA, 19), 
                cex = 1.2, col = "blue")


##### Plot #####
plot(fg_dend, ylab = "Height", nodePar = nodePar, leaflab = "perpendicular")


##### Specify different point types and colors for each leave #####
cores_fg <- c("#7e9828","#f7a436","#7e9828","#f7a436","#006a00", "#7e9828","#006a00","#7e9828","#006a00", "#006a00", "#006a00",
              "#006a00", "#f7a436", "#006a00", "#f7a436","#7e9828", "#f7a436","#ffdb6d", "#ffdb6d","#006a00", "#f7a436","#7e9828",
              "#f7a436", "#ffdb6d", "#ffdb6d", "#006a00", "#7e9828","#7e9828","#ffdb6d", "#ffdb6d", "#ffdb6d" ) 
  
              
fg_dend %>% set("leaves_pch", 19) %>%  # node point type
  set("leaves_cex", 1.5) %>%  # node point size
  set("leaves_col", cores_fg) %>%
  set("labels_cex",  1)%>% #node point color
  plot(main = "Fungi - cluster dendrogram")

fp_dend <- recordPlot()

                                          
##### Plot orders in the dataset #####
##### Set the order of the dendrogram and barplot #####
sample_order_fg <- sample_names(fg.rarefied)[order.dendrogram(fg_dend)]

fungi.tax
data_fg <- psmelt(fg.rarefied)
fp <- ggplot(data=data_fg, aes(x=Sample, y=Abundance, fill=phylum))
barplot_fg <- fp + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#E9E681", "#F4D03F", "#D0AD32", "#967E2C", "#6b451fff", 
                               "#82edd8ff", "#45B39D", "#52BE80", "#428345ff", "#13522eff", 
                               "#F4C9E0", "#E6B0AA", "#7e4c4699", "#EC7063", "#a16159ff", 
                               "#d6e6f5ff", "#bad2ebff", "#8ebedaff", "#357fb9ff", "#0a3d8aff", 
                               "#ede5e5ff", "#d3d3d3ff", "#8c7c7c9a", "#6e6f6fff", "#4e4e4eff")) +
  theme(legend.position="bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
                                          guides(fill=guide_legend(nrow=5)) + scale_x_discrete(limits=sample_order_fg) 


barplot_fg <- barplot_fg + theme_bw()
barplot_fg <- barplot_fg + theme(legend.position="bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(fill=guide_legend(nrow=5)) +
scale_x_discrete(limits=sample_order_fg) 


##### Plot most abundant orders for fungi data #####
top25_fg <- names(sort(taxa_sums(fg.rarefied), decreasing=TRUE)[1:106])

dat.aglo_ord = tax_glom(fg.rarefied, taxrank = "order")
dat.trans_ord = transform_sample_counts(dat.aglo_ord, function(x) x/sum(x))
prune.dat.two_ord = prune_taxa(top25_fg, dat.trans_ord)
dat.dataframe_ord = psmelt(prune.dat.two_ord)
dat.agr_ord = aggregate(Abundance~Sample+type+order, data=dat.dataframe_ord, FUN=mean)



fg_bar_order <- ggplot(dat.agr_ord, aes(x=Sample, y=Abundance, fill=order)) + geom_bar(stat="identity", position="fill")  + 
                                        theme(legend.position="bottom", legend.key.size = unit(0.4, 'cm')) + 
                                        scale_fill_discrete(breaks=c("Agaricales", "Archaeorhizomycetales", "Chaetosphaeriales", "Chaetothyriales", "Diaporthales",
                                                                     "Geastrales", "Glomerales", "Glomerellales", "Helotiales", "Hymenochaetales", "Hypocreales", 
                                                                     "Mortierellales", "Mucorales", "Phallales", "Pleosporales", "Polyporales", "Russulales", "Sarrameanales", 
                                                                     "Sordariales", "Thelephorales", "Trechisporales", "Trichosporonales", "Xylariales", "unidentified"))
                        
fg_bar_order
barplot_fg_order <- fg_bar_order + theme_bw() + scale_fill_manual(values =c("#E9E681", "#F4D03F", "#D0AD32", "#967E2C", "#6b451fff", 
                                                                            "#82edd8ff", "#45B39D", "#52BE80", "#428345ff", "#13522eff", 
                                                                            "#F4C9E0", "#E6B0AA", "#7e4c4699", "#EC7063", "#a16159ff", 
                                                                            "#d6e6f5ff", "#bad2ebff", "#8ebedaff", "#357fb9ff", "#0a3d8aff", 
                                                                            "#ede5e5ff", "#d3d3d3ff", "#8c7c7c9a", "#6e6f6fff", "#4e4e4eff")) 

barplot_fg_order <- barplot_fg_order + theme(legend.position="bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + guides(fill=guide_legend(nrow=5)) +
  scale_x_discrete(limits=sample_order_fg) 

barplot_fg_order <- barplot_fg_order +  scale_fill_discrete(breaks=c("Agaricales", "Archaeorhizomycetales", "Chaetosphaeriales", "Chaetothyriales",
                                                                     "Diaporthales", "Geastrales", "Glomerales", "Glomerellales", "Helotiales", "Hymenochaetales",
                                                                     "Hypocreales", "Mortierellales", "Mucorales", "Phallales", "Pleosporales", "Polyporales",
                                                                     "Russulales", "Sarrameanales", "Sordariales", "Thelephorales", "Trechisporales", "Trichosporonales",
                                                                     "Xylariales", "unidentified"))
                                        
barplot_fg_order <- barplot_fg_order +  scale_fill_discrete(breaks=c("Agaricales", "Archaeorhizomycetales", "Chaetosphaeriales", "Chaetothyriales",
                                                                       "Diaporthales", "Geastrales", "Glomerales", "Glomerellales", "Helotiales", "Hymenochaetales",
                                                                       "Hypocreales", "Mortierellales", "Mucorales", "Phallales", "Pleosporales", "Polyporales",
                                                                       "Russulales", "Sarrameanales", "Sordariales", "Thelephorales", "Trechisporales", "Trichosporonales",
                                                                       "Xylariales", "unidentified")) + 
                                        scale_fill_manual(values =c("#E9E681", "#F4D03F", "#D0AD32", "#967E2C", "#6b451fff", 
                                                                   "#82edd8ff", "#45B39D", "#52BE80", "#428345ff", "#13522eff", 
                                                                   "#F4C9E0", "#E6B0AA", "#7e4c4699", "#EC7063", "#a16159ff", 
                                                                   "#d6e6f5ff", "#bad2ebff", "#8ebedaff", "#357fb9ff", "#0a3d8aff", 
                                                                   "#ede5e5ff", "#d3d3d3ff", "#8c7c7c9a", "#6e6f6fff", "#4e4e4eff"))  
plot_grid(fp_dend, barplot_fg_order,
          nrow = 2,
          labels = "AUTO",
          label_size = 12,
          align = "v", axis = "rl"
)



#----------------------------------------------------------#
# 5. Figure 1 Bipartite network based on indicator species #
#----------------------------------------------------------#
                                        ##### Estimating land use system responsive OTUs with indicator species analysis for plants #####
indic_roots <- as.data.frame(t(otu.rare))
indic_roots_groups <- roots.env$type
length(unique(indic_roots_groups))


## set.seed(5971)
indicatorsp_roots <- multipatt(indic_roots,indic_roots_groups,func = "IndVal.g",control=how(nperm=9999))

summary(indicatorsp_roots,alpha=1,indvalcomp=T)
indic_roots_df <- indicatorsp_roots$sign
write.table(indic_roots_df,paste0("indic_roots_df.txt"),sep="\t",quote=F)


##### Select OTUs classified as indicator species with significant p-value #####
forest_roots <- as.matrix(indic_roots_df[which(indic_roots_df$s.Forest == 1 & indic_roots_df$p.value < 0.05),])
jungle_roots <- as.matrix(indic_roots_df[which(indic_roots_df$s.Jungle == 1 & indic_roots_df$p.value < 0.05),])
rubber_roots <- as.matrix(indic_roots_df[which(indic_roots_df$s.Rubber == 1 & indic_roots_df$p.value < 0.05),])
oil_roots <- as.matrix(indic_roots_df[which(indic_roots_df$s.OilPalm == 1 & indic_roots_df$p.value < 0.05),])

roots_r_values <- rbind(forest_roots,jungle_roots,rubber_roots,oil_roots)
colnames(roots_r_values)[1:4] <-c("Forest","Jungle","OilPalm","Rubber")

dim(roots_r_values)
##### Check range of correlation coefficients #####
range(roots_r_values[,"stat"])


##### Check total number of indicator OTUS #####
length(unique(rownames(roots_r_values)))


##### Proportion of roots OTUs correlated to each land use system #####
length(unique(rownames(roots_r_values)))/nrow(indic_roots)
roots_ra <- t(t(otu.rare)/colSums(otu.rare)) * 100
sum(colSums(roots_ra[unique(rownames(roots_r_values)),]))/sum(colSums(roots_ra))


##### Construct node table for plant communities from indicator species data #####
plants_bipartite <- data.frame(from= c(rep("Forest",length(which(roots_r_values[,"Forest"]==1))),
                                      rep("Jungle",length(which(roots_r_values[,"Jungle"]==1))),
                                      rep("OilPalm",length(which(roots_r_values[,"OilPalm"]==1))),
                                      rep("Rubber",length(which(roots_r_values[,"Rubber"]==1)))),
                              to= c(rownames(roots_r_values)[which(roots_r_values[,"Forest"]==1)],
                                    rownames(roots_r_values)[which(roots_r_values[,"Jungle"]==1)],
                                    rownames(roots_r_values)[which(roots_r_values[,"OilPalm"]==1)],
                                    rownames(roots_r_values)[which(roots_r_values[,"Rubber"]==1)]),
                              r= c(roots_r_values[which(roots_r_values[,"Forest"]==1),"stat"],
                                   roots_r_values[which(roots_r_values[,"Jungle"]==1),"stat"],
                                   roots_r_values[which(roots_r_values[,"OilPalm"]==1),"stat"],
                                   roots_r_values[which(roots_r_values[,"Rubber"]==1),"stat"]))


##### Plot bipartite network for plant data #####
r_bipartite <- graph.data.frame(plants_bipartite, directed=FALSE)
bipartite.mapping(r_bipartite)

V(r_bipartite)$type <- bipartite_mapping(r_bipartite)$type
plot(r_bipartite)

plot(r_bipartite, vertex.label.cex = 0.8, vertex.label.color = "black")
plot(r_bipartite, layout=layout.bipartite, vertex.size=7, vertex.label.cex=0.6)


##### Define attributes of network plots #####
## Color and shape
V(r_bipartite)$color <- ifelse(V(r_bipartite)$type, "darkgreen", "yellow")
V(r_bipartite)$shape <- ifelse(V(r_bipartite)$type, "circle", "square")
E(r_bipartite)$color <- "lightgray"

plot(r_bipartite, vertex.label.cex = 0.8, vertex.label.color = "black")
plot(r_bipartite, node.size= .4, vertex.label.cex = 0.4, vertex.label.color = "black")
plot(r_bipartite, node.size= .4, node.label= V(r_bipartite)$family, vertex.label.cex = 0.4, vertex.label.color = "black")


##### Some network metrics #####
types <- V(r_bipartite)$type                 
deg <- degree(r_bipartite)
bet <- betweenness(r_bipartite)
clos <- closeness(r_bipartite)
eig <- eigen_centrality(r_bipartite)$vector

cent_df <- data.frame(types, deg, bet, clos, eig)
cent_df[order(cent_df$type, decreasing = TRUE),] 


##### Sizing Vertices by Centrality #####
V(r_bipartite)$size <- degree(r_bipartite)
V(r_bipartite)$label.cex <- degree(r_bipartite) * 0.03
plot(r_bipartite, layout = layout_with_fr)


##### Prepare node attribute table for each OTU #####
root_bipartite_attrib <- data.frame(node=unique(rownames(roots_r_values)),indicgroup=0)

for (i in as.character(root_bipartite_attrib$node))
{
  root_bipartite_attrib[root_bipartite_attrib$node==i,"indicgroup"] <- paste(colnames(roots_r_values)[which(roots_r_values[i,1:4]==1)],collapse = "_")
}

root_bipartite_attrib <- cbind(root_bipartite_attrib,roots.tax[as.character(root_bipartite_attrib$node),])


##### Set labels for nodes #####
V(r_bipartite)$label <- V(r_bipartite)$name


##### Set node sizes #####
V(r_bipartite)$size <- 6


##### Set node shapes #####
V(r_bipartite)$shape <- c(rep("circle",4),rep("circle",158))


##### Link taxonomy able to otu table #####
tax_filter_roots <- roots.tax[rownames(otu.rare),]


##### Define node colors based on order assignment #####
V(r_bipartite)$color <- V(r_bipartite)$name
V(r_bipartite)$color[1:4] <- root_bipartite_attrib[V(r_bipartite)$name, ]$name
V(r_bipartite)$color <- root_bipartite_attrib[ V(r_bipartite)$color, ]$Order
V(r_bipartite)$color[1:4] <- "white"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Arecales"] <- "#522d05ff"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Ericales"] <- "#8c510aff"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Fabales"] <- "#c87137ff"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Gentianales"] <- "#edd1bfff"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Lamiales"] <- "#a1cca5"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Laurales"] <- "#90a955"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Malpighiales"] <- "#01665eff"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Myrtales"] <- "#8ebedaff"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Poales"] <- "#0c47a0ff"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Rosales"] <- "#e8ddb5"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Santalales"] <- "#d5c8e0"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Sapindales"] <- "#cccccc"
V(r_bipartite)$color[V(r_bipartite)$color %in% "Zingiberales"] <- "#595959"


##### Figure 1 A: Plot Plant Bipartite Network #####
plot(r_bipartite, node.size= 1,edge.curved=0.2,layout=layout_with_fr(r_bipartite),
     edge.arrow.size=3, vertex.label= NA, vertex.frame.color="#d2d2d2ff",vertex.frame.width= 0.5, vertex.label.cex = NA,
     vertex.size=11,vertex.color=V(r_bipartite)$color, vertex.shape="circle") 

Figure1A <-recordPlot()


##### Figure 1 B: Relative abundance of indicator plant OTUs per order#####
ind_plants <- V(r_bipartite)$label
dat.trans_ord = transform_sample_counts(dat.aglo_ord, function(x) x/sum(x))
prune.dat.two_ord = prune_taxa(ind_plants, roots.rarefied)
dat.dataframe_ind = psmelt(prune.dat.two_ord)
dat.agr_ind = aggregate(Abundance~Sample+type+Order, data=dat.dataframe_ind, FUN=mean)

p <- ggplot(data=dat.agr_ind, aes(x=type, y=Abundance, fill=Order))
Figure1B <- p + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values =  c("#522d05ff",  "#8c510aff",  "#c87137ff", "#edd1bfff", "#a1cca5",  "#90a955",
                                "#01665eff", "#8ebedaff",  "#0c47a0ff", "#e8ddb5",  "#b4a3cbff",  "#cccccc", "#595959")) +
  theme(legend.position="bottom", axis.text.x = element_text(vjust = 0.5, hjust=1), legend.key.size = unit(0.2, 'cm')) + guides(fill=guide_legend(nrow=9))


Figure1B <- Figure1B + theme_bw() + scale_x_discrete(limits = c("Forest","Jungle","Rubber", "OilPalm"))
Figure1B <- Figure1B + theme(legend.position="right", axis.text.x = element_text(vjust = 0.5, hjust=1))  


##### Figure 1 C: Relative abundance of indicator plant OTUs functional groups #####
dat.agr_ind_cat = aggregate(Abundance~Sample+type+Category, data=dat.dataframe_ind, FUN=mean)

p1 <- ggplot(data=dat.agr_ind_cat, aes(x=type, y=Abundance, fill=Category))
Figure1C <- p1 + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#018571ff", "#80cdc1ff", "#dfc27dff")) +
  theme(legend.position="bottom", axis.text.x = element_text(vjust = 0.5, hjust=1), legend.key.size = unit(0.2, 'cm')) + guides(fill=guide_legend(nrow=3))


Figure1C <- Figure1C + theme_bw()
Figure1C <- Figure1C + theme(legend.position="bottom", axis.text.x = element_text(vjust = 0.5, hjust=1))  


##### Relative abundance of indicator plant OTUs native/non-native #####
dat.agr_ind_nat = aggregate(Abundance~Sample+type+Native_Alien, data=dat.dataframe_ind, FUN=mean)

nat_plot <- ggplot(data=dat.agr_ind_nat, aes(x=type, y=Abundance, fill=Native_Alien))
nat_plot <- nat_plot + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#018571ff", "#80cdc1ff")) +
  theme(legend.position="bottom", axis.text.x = element_text(vjust = 0.5, hjust=1), legend.key.size = unit(0.2, 'cm')) + guides(fill=guide_legend(nrow=3))


nat_plot <- nat_plot + theme_bw()


##### Estimating land use system responsive OTUs with indicator species analysis for fungi #####
indic_fungi <- as.data.frame(t(otu_fg))
indic_fg_groups <- fungi.env$type
length(unique(indic_fg_groups))


## set.seed(5971)
indicatorsp_fg <- multipatt(indic_fungi,indic_fg_groups,func = "IndVal.g",control=how(nperm=9999))

summary(indicatorsp_fg,alpha=1,indvalcomp=T)
indic_fg_df <- indicatorsp_fg$sign
write.table(indic_fg_df,paste0("indic_fg_df.txt"),sep="\t",quote=F)


##### Select OTUs classified as indicator species with significant p-value #####
forest_fg <- as.matrix(indic_fg_df[which(indic_fg_df$s.Forest == 1 & indic_fg_df$p.value < 0.05),])
jungle_fg <- as.matrix(indic_fg_df[which(indic_fg_df$s.Jungle == 1 & indic_fg_df$p.value < 0.05),])
rubber_fg <- as.matrix(indic_fg_df[which(indic_fg_df$s.Rubber == 1 & indic_fg_df$p.value < 0.05),])
oil_fg <- as.matrix(indic_fg_df[which(indic_fg_df$s.OilPalm == 1 & indic_fg_df$p.value < 0.05),])

fg_r_values <- rbind(forest_fg,jungle_fg,rubber_fg,oil_fg)
colnames(fg_r_values)[1:4] <-c("Forest","Jungle","OilPalm","Rubber")


##### Check range of correlation coefficients #####
range(fg_r_values[,"stat"])


##### Check total number of indicator OTUS #####
length(unique(rownames(fg_r_values)))


##### Proportion of roots OTUs correlated to each land use system #####
length(unique(rownames(fg_r_values)))/nrow(indic_fungi)
fungi_ra <- t(t(otu_fg)/colSums(otu_fg)) * 100
sum(colSums(fungi_ra[unique(rownames(fg_r_values)),]))/sum(colSums(fungi_ra))


##### Construct node table for plant communities from indicator species data #####
fungi_bipartite <- data.frame(from= c(rep("Forest",length(which(fg_r_values[,"Forest"]==1))),
                                       rep("Jungle",length(which(fg_r_values[,"Jungle"]==1))),
                                       rep("OilPalm",length(which(fg_r_values[,"OilPalm"]==1))),
                                       rep("Rubber",length(which(fg_r_values[,"Rubber"]==1)))),
                               to= c(rownames(fg_r_values)[which(fg_r_values[,"Forest"]==1)],
                                     rownames(fg_r_values)[which(fg_r_values[,"Jungle"]==1)],
                                     rownames(fg_r_values)[which(fg_r_values[,"OilPalm"]==1)],
                                     rownames(fg_r_values)[which(fg_r_values[,"Rubber"]==1)]),
                               r= c(fg_r_values[which(fg_r_values[,"Forest"]==1),"stat"],
                                    fg_r_values[which(fg_r_values[,"Jungle"]==1),"stat"],
                                    fg_r_values[which(fg_r_values[,"OilPalm"]==1),"stat"],
                                    fg_r_values[which(fg_r_values[,"Rubber"]==1),"stat"]))


##### Plot bipartite network for plant data #####
fg_bipartite <- graph.data.frame(fungi_bipartite, directed=FALSE)
bipartite.mapping(fg_bipartite)

V(fg_bipartite)$type <- bipartite_mapping(fg_bipartite)$type
plot(fg_bipartite)

plot(fg_bipartite, vertex.label.cex = 0.8, vertex.label.color = "black")
plot(fg_bipartite, layout=layout.bipartite, vertex.size=7, vertex.label.cex=0.6)


##### Define attributes of network plots #####
## Color and shape
V(fg_bipartite)$color <- ifelse(V(fg_bipartite)$type, "darkgreen", "yellow")
V(fg_bipartite)$shape <- ifelse(V(fg_bipartite)$type, "circle", "square")
E(fg_bipartite)$color <- "lightgray"

plot(fg_bipartite, vertex.label.cex = 0.8, vertex.label.color = "black")
plot(fg_bipartite, node.size= .4, vertex.label.cex = 0.4, vertex.label.color = "black")
plot(fg_bipartite, node.size= .4, node.label= V(fg_bipartite)$family, vertex.label.cex = 0.4, vertex.label.color = "black")


##### Some network metrics #####
types_fg <- V(fg_bipartite)$type                 
deg_fg <- degree(fg_bipartite)
bet_fg <- betweenness(fg_bipartite)
clos_fg <- closeness(fg_bipartite)
eig_fg <- eigen_centrality(fg_bipartite)$vector

cent_df_fg <- data.frame(types_fg, deg_fg, bet_fg, clos_fg, eig_fg)
cent_df_fg[order(cent_df_fg$types_fg, decreasing = TRUE),] 


##### Sizing Vertices by Centrality #####
V(fg_bipartite)$size <- degree(fg_bipartite)
V(fg_bipartite)$label.cex <- degree(fg_bipartite) * 0.03
plot(fg_bipartite, layout = layout_with_fr)


##### Prepare node attribute table for each OTU #####
fg_bipartite_attrib <- data.frame(node=unique(rownames(fg_r_values)),indicgroup=0)

for (i in as.character(fg_bipartite_attrib$node))
{
  fg_bipartite_attrib[fg_bipartite_attrib$node==i,"indicgroup"] <- paste(colnames(fg_r_values)[which(fg_r_values[i,1:4]==1)],collapse = "_")
}

fg_bipartite_attrib <- cbind(fg_bipartite_attrib,fungi.tax[as.character(fg_bipartite_attrib$node),])
write_xlsx(fg_bipartite_attrib,"fg_bipartite_attrib.xlsx")


##### Set labels for nodes #####
V(fg_bipartite)$label <- V(fg_bipartite)$name


##### Set node sizes #####
V(fg_bipartite)$size <- 6


##### Set node shapes #####
V(fg_bipartite)$shape <- c(rep("circle",4),rep("circle",158))


##### Link taxonomy able to otu table #####
tax_filter_fg <- fungi.tax[rownames(otu_fg),]


##### Define node colors based on phylum assignment #####
V(fg_bipartite)$color <- V(fg_bipartite)$label
V(fg_bipartite)$color[1:4] <- fg_bipartite_attrib[V(fg_bipartite)$name, ]$name
V(fg_bipartite)$color <- fg_bipartite_attrib[ V(fg_bipartite)$color, ]$phylum
V(fg_bipartite)$color[1:4] <- "white"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "Ascomycota"] <- "#522d05ff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "Basidiomycota"] <- "#8f4f06ff"
#V(fg_bipartite)$color[V(fg_bipartite)$color %in% "Chytridiomycota"] <- "#d5bd8cff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "Calcarisporiellomycota"] <- "#8a802cff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "Glomeromycota"] <-"#f6e7c2ff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "Kickxellomycota"] <- "#f4f4f4ff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "Mortierellomycota"] <- "#cbe3dfff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "Mucoromycota"] <- "#86c7bcff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "Rozellomycota"] <- "#559a93ff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "unidentified"] <- "#01665eff"
V(fg_bipartite)$name[5:268] <- NA


##### Figure 1 D: Plot fungi Bipartite Network #####
plot(fg_bipartite, node.size= 1,edge.curved=0.2,layout=layout_with_fr(fg_bipartite),
     edge.arrow.size=2, vertex.label= V(fg_bipartite)$name, vertex.frame.color="#d2d2d2ff",vertex.frame.width= 0.5, vertex.label.cex = NA,
     vertex.size=11,vertex.color=V(fg_bipartite)$color, vertex.shape=V(fg_bipartite)$shape) 

Figure1D <-recordPlot()


##### Define node colors based on functional guilds assignment #####
V(fg_bipartite)$color <- V(fg_bipartite)$label
V(fg_bipartite)$color[1:4] <- fg_bipartite_attrib[V(fg_bipartite)$name, ]$name
V(fg_bipartite)$color <- fg_bipartite_attrib[ V(fg_bipartite)$color, ]$groups
V(fg_bipartite)$color[1:4] <- "white"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "AM"] <- "#522d05ff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "ECM"] <- "#a6611aff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "saprotroph"] <- "#80cdc1ff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "plant_pathogen"] <-"#f5f5f5ff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "unknown"] <- "#01665eff"
V(fg_bipartite)$color[V(fg_bipartite)$color %in% "other"] <- "#01665eff"

V(fg_bipartite)$name[5:268] <- NA

plot(fg_bipartite, node.size= 1,edge.curved=0.2,layout=layout_with_fr(fg_bipartite),
     edge.arrow.size=2, vertex.label= V(fg_bipartite)$name, vertex.frame.color="#d2d2d2ff",vertex.frame.width= 0.5, vertex.label.cex = NA,
     vertex.size=11,vertex.color=V(fg_bipartite)$color, vertex.shape=V(fg_bipartite)$shape) 

##### Figure 1 E: Relative abundance of indicator plant OTUs per phylum for fungi#####
ind_fungi <- V(fg_bipartite)$label

dat.trans_phy = transform_sample_counts(fg.rarefied, function(x) x/sum(x))
prune.dat.two_phy = prune_taxa(ind_fungi, dat.trans_phy)
dat.dataframe_phy = psmelt(prune.dat.two_phy)
dat.agr_phy = aggregate(Abundance~Sample+type+phylum, data=dat.dataframe_phy, FUN=mean)

fg_phy <- ggplot(data=dat.agr_phy, aes(x=type, y=Abundance, fill=phylum))
Figure1E <- fg_phy + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#522d05ff", "#8f4f06ff", "#8a802cff", "#f6e7c2ff", 
                               "#f4f4f4ff", "#cbe3dfff", "#86c7bcff", 
                               "#559a93ff", "#076260ff")) +
  theme(legend.position="right", axis.text.x = element_text(vjust = 0.5, hjust=1), legend.key.size = unit(0.2, 'cm')) + guides(fill=guide_legend(nrow=15))


Figure1E <- Figure1E + theme_bw() + scale_x_discrete(limits = c("Forest","Jungle","Rubber", "OilPalm"))


##### Figure 1 F: Relative abundance of indicator plant OTUs functional groups #####
dat.agr_guild = aggregate(Abundance~Sample+type+groups, data=dat.dataframe_phy, FUN=mean)
dat.agr_guild$groups[dat.agr_guild$groups %in% "other"] <- "unknown"
fg_phy_guild <- ggplot(data=dat.agr_guild, aes(x=type, y=Abundance, fill=groups))



p_guild <- ggplot(data=dat.agr_guild, aes(x=type, y=Abundance, fill=groups))
Figure_1F <- p_guild + geom_bar(aes(), stat="identity", position="fill") +
  scale_fill_manual(values = c("#522d05ff", "#a6611aff", "#f5f5f5ff", "#80cdc1ff", 
                               "#076260ff")) +
  theme(legend.position="bottom", axis.text.x = element_text(vjust = 0.5, hjust=1), legend.key.size = unit(0.2, 'cm')) + guides(fill=guide_legend(nrow=6))


Figure_1F <- Figure_1F + theme_bw() +  scale_x_discrete(limits = c("Forest","Jungle","Rubber", "OilPalm"))
Figure_1F

                                        
#---------------------------------------------------------------------#
# 6.Table 1 Network metrics of roots and fungi co-occurrence network  #
#---------------------------------------------------------------------#
##### Combine otu tables from plants and fungi #####
otu_norm_combine <- rbind(otu.rare, otu_fg)


##### Perform Spearman correlation of all OTU pairs #####
all_cor <- rcorr(t(otu_norm_combine), type=c("spearman"))


##### Extract the correlation coefficients and p values #####
CorrDF <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    from = rownames(cormat)[col(cormat)[ut]],
    to = rownames(cormat)[row(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


all_cor_df <- CorrDF(all_cor$r, all_cor$P)


##### Adjust P-values for Multiple Comparisons #####
all_cor_df$padj <- p.adjust(all_cor_df$p, method="none")


##### Extract only correlations above 0.7 and p value 0.001 #####
all_cor_df_padj <- all_cor_df[which(all_cor_df$cor > 0.7),]
all_cor_df_padj <-all_cor_df_padj[which(all_cor_df_padj$padj < 0.001),]
write_xlsx(all_cor_df_padj, "all_cor_df_padj.xlsx")


##### Modularity and topological roles #####
##### Prepare input file for netcarto (modularity) analysis #####
## It takes time, I run it the cluster - check netcarto.sh file 
nd1 = all_cor_df_padj$from
nd2 = all_cor_df_padj$to
weight = all_cor_df_padj$cor
web = list(nd1,nd2,weight)

library(rnetcarto)
library(readr)
library(readxl)

netcarto_out098_1707 <- netcarto(
  web,
  seed = 2987,
  iterfac = 1,
  coolingfac = 0.98,
  bipartite = TRUE
)

write.table(netcarto_out098_1707,paste0("netcarto_out098_1707.txt"),sep="\t",quote=F)


##### Estimate OTUs responsive to landuse systemm using edge R #####
## Fungi OTU responsive to land uses with likelihood ratio testing in edgeR
model_fungi <- model.matrix(~landuse, data=fungi.env)
fungi_tax_rare <- tax_table(fg.rarefied)
fungi_tax_rare <- as.data.frame(fungi_tax_rare)
otu_fg_df <- as.data.frame(otu_fg)

edgeR_fungi_system <- DGEList(counts=otu_fg_df, group=fungi.env$landuse, genes=fungi_tax_rare)
edgeR_fungi_system_factors <- calcNormFactors(edgeR_fungi_system)

dge_fungifsyst <- estimateGLMRobustDisp(edgeR_fungi_system_factors, design=model_fungi)

fit_fungi <- glmFit(dge_fungifsyst, design = model_fungi)
lrt_fungi <- glmLRT(fit_fungi, coef=2:4)
fsyst_fungi <- topTags(lrt_fungi, n=Inf, p.value=0.05)
fsyst_fungi <- fsyst_fungi$table


## Plant OTU responsive to land uses with likelihood ratio testing in edgeR
model_plants <- model.matrix(~type, data=roots.env)
plants_tax_rare <- tax_table(roots.rarefied)
plants_tax_rare <- as.data.frame(plants_tax_rare)
otu_pl_df <- as.data.frame(otu.rare)

edgeR_plants_fsyst <- DGEList(counts=otu_pl_df, group=roots.env$type, genes=plants_tax_rare)
edgeR_plants_fsyst <- calcNormFactors(edgeR_plants_fsyst)


dge_plants_fsyst <- estimateGLMRobustDisp(edgeR_plants_fsyst, design=model_plants)

fit_plants <- glmFit(dge_plants_fsyst, design = model_plants)
lrt_plants <- glmLRT(fit_plants, coef=2:4)
fsyst_plants <- topTags(lrt_plants, n=Inf, p.value=0.05)
fsyst_plants <- fsyst_plants$table


##### Identify land use sensitive taxa #####
indic_edge_plants <- intersect(rownames(roots_r_values), rownames(fsyst_plants))
indic_edge_fungi <- intersect(rownames(fg_r_values), rownames(fsyst_fungi))

indic_edge_combine <- c(indic_edge_plants, indic_edge_fungi)


##### Create attribute table for node values #####
nodeattrib_combine <- data.frame(node=union(all_cor_df_padj$from,all_cor_df_padj$to))
r_values_combine <- rbind(roots_r_values,fg_r_values)


nodeattrib_combine$indicgroup <- 0

for (i in as.character(nodeattrib_combine$node))
{
  if (i %in% indic_edge_combine == TRUE)
  {nodeattrib_combine[nodeattrib_combine$node==i,"indicgroup"] <- paste(colnames(r_values_combine)[which(r_values_combine[i,1:4]==1)],collapse = "_")}
  else
  {nodeattrib_combine[nodeattrib_combine$node==i,"indicgroup"]<- "NA"}
}

rownames(nodeattrib_combine) <- as.character(nodeattrib_combine$node)

dim(indic_edge_combine)


#####  Create network #####
all_net <- graph_from_data_frame(all_cor_df_padj,direct=F,vertices=nodeattrib_combine)


##### Count and Characterize nodes #####    

length(V(all_net))
length(grep("OTU_97_*",names(V(all_net))))
rtall <- sub("OTU_97_*", "f_97_", names(V(all_net)))
rtall1 <- rename("OTU_97_*", "f_97_", all_cor_df_padj$from)
length(grep("OTU*", rtall))


df_renamed <- all_cor_df_padj %>%
  mutate(across(starts_with("from"), ~str_replace(., "^OTU_97_", "fungi_")))
df_renamed <- df_renamed %>%
  mutate(across(starts_with("to"), ~str_replace(., "^OTU_97_", "fungi_")))


##### Connections #####
fg_occur <- droplevels(df_renamed[with(df_renamed, grepl("fungi_*",from) & grepl("fungi_*",to)),])
nrow(fg_occur)

rt_occur <- droplevels(df_renamed[with(df_renamed, grepl("OTU*",from) & grepl("OTU*",to)),])
nrow(rt_occur)

rtfg_occur <- droplevels(df_renamed[with(df_renamed, grepl("fungi_*",from) & grepl("OTU*",to)),])
nrow(rtfg_occur)

write_xlsx(rt_occur, "rt_occur.xlsx")
dim(rt_occur)


##### Estimate node degree values #####
deg <- degree(all_net, mode="all")
all_deg <- sort(degree(all_net,mode="all"),decr=T)
mean(all_deg)


##### Community structure and definition of modules that contain indicspecies in the landuses #####
Forest_nodes <- rownames(nodeattrib_combine[nodeattrib_combine$indicgroup=="Forest",])
Jungle_nodes <- rownames(nodeattrib_combine[nodeattrib_combine$indicgroup=="Jungle",])
Rubber_nodes <- rownames(nodeattrib_combine[nodeattrib_combine$indicgroup=="Rubber",])
OilPalm_nodes <- rownames(nodeattrib_combine[nodeattrib_combine$indicgroup=="OilPalm",])


cs_nodes_all <- c(Forest_nodes,Jungle_nodes,Rubber_nodes,OilPalm_nodes)
fg_nodes_all <- cs_nodes_all[grep("OTU_97_*",cs_nodes_all)]

cs_nodes_all1 <- sub("OTU_97_*", "f_97_",cs_nodes_all)
rt_nodes_all <- cs_nodes_all[grep("OTU*",cs_nodes_all1)]

                                        
#---------------------------------------------------------------------#
# 7.Table S2 Network metric results                                   #
#---------------------------------------------------------------------#
##### Merge indicator table results with taxonomic table and rnetcarto results #####
##### Upload netcarto results table #####
TableS2_Network <- read_excel("~/Network_R/TableS2_Network.xlsx")

colnames(roots.tax)[1] <- "OTU"
TableS2 <- merge(TableS2_Network, roots.tax, by="OTU")
write_xlsx(TableS2, "TableS2_plant.xlsx")

colnames(roots.tax)[1] <- "OTU"
TableS2_fg <- merge(TableS2_Network, fungi.tax, by="OTU")
write_xlsx(TableS2_fg, "TableS2_fg.xlsx")

colnames(fg_r_values)[1] <- "Forest"
head(fg_r_values)

colnames(nodeattrib_combine)[1] <- "OTU"
TableS2_ind <- merge(TableS2_Network, nodeattrib_combine, by="OTU")
write_xlsx(TableS2_ind, "TableS2_ind.xlsx")
          
