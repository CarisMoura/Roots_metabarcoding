#-----------------------------------------------------------------------------------------------------------------------------------------#
#                                                                                                                                         #
#   Codes for: Moura et al. "Insights on root-fungal association in forest conversion systems through  network and indicator analysis"    # 
#                                                                                                                                         #
#-----------------------------------------------------------------------------------------------------------------------------------------#

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
BiocManager::install("limma")

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

