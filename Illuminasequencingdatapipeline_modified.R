---
title: "Iluminasequencingdatapipeline"
author: "Sydney Schmitter and Tanvi Dutta Gupta"
date: "10/5/2022"
output: html_document
---
#Import Libraries
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(censusapi)
#library(sf)
library(tigris)
library(mapview)
library(readxl)
library(leaflet)
library(ggplot2)
library(shiny)
library(shinydashboard)
library(vegan)
library(sjPlot)
library(ggeffects)
library(scales)
```

#Anacapa Pipeline
```{r}
# __Anacapa__ is an eDNA toolkit that allows users to build comprehensive reference databases and assign taxonomy to raw multilocus metabarcode sequence data. It address longstanding needs of the eDNA for modular informatics tools, comprehensive and customizable reference databases, flexibility across high-throughput sequencing platforms, fast multilocus metabarcode processing, and accurate taxonomic assignment. __Anacapa__ toolkit processes eDNA reads and assigns taxonomy using existing software or modifications to existing software. This modular toolkit is designed to analyze multiple samples and metabarcodes simultaneously from any Ilumina sequencing platform. A significant advantage of the __Anacapa__ toolkit is that it does not require that paired reads overlap, or that both reads in a pair pass QC. Taxonomy results are generated for all read types and the user can decide which read types they wish to retain for downstream analysis.

# __Anacapa__ includes four modules:
# 1) building reference libraries using [__CRUX__](https://github.com/limey-bean/CRUX_Creating-Reference-libraries-Using-eXisting-tools)
# 2) running quality control (QC) and assigning Amplicon Sequence Variants (ASV) using Dada2 (Sequence QC and ASV Parsing),   
# 3) assigning taxonomy using Bowtie 2 and a Bowtie 2 specific Bayesian Least Common Ancestor (BLCA) (Assignment) and   
# 4) Running exploratory data analysis and generating ecological diversity summary statistics for a set of samples using [ranacapa](https://f1000research.com/articles/7-1734/v1).  
```

#Step 1: CRUX
```{r}
# This first part of the Anacapa toolkit uses __CRUX__ to generate reference libraries needed for taxonomic assignment.  Briefly,the output of __CRUX__ consists of two reference libraries, either unfiltered or filtered. Unfiltered libraries contain every dereplicated read found during BLAST searches. The filtered library contains only reads with robust taxonomic assignments. Specifically we refer to robust taxonomic assignments as any reads that do not have the following in their taxonomic path: 'uncultured', 'environmental', 'sample', or 'NA;NA;NA;NA'. Prebuilt __CRUX__ reference libraries (12S - MiFish, 16S - EMP, 18S V4, 18S V8-9, 18S - EMP, PITS - Plant ITS2, CO1 and FITS - Fungal ITS) [see Table 1] are available on [Google Drive](https://drive.google.com/drive/folders/0BycoA83WF7aNOEFFV2Z6bC1GM1E?usp=sharing)), and will be uploaded to Dryad upon acceptence of the Anacapa manuscript. Each library contains unique metabarcode specific reads that correspond to NCBI accession version numbers. Libraries consist of fasta files, taxonomy files, and a Bowtie 2 index library.  
```

# Step 2: Sequence QC and ASV Parsing using dada2
```{r}
# This next step of the toolkit aims to conduct standard sequence QC and generate amplicon sequence variants (ASV) from Illumina data using **dada2** (Callahan et al. 2016). ASVs are a novel solution to identifying biologically informative unique sequences in metabarcoding samples that replaces the operational taxonomic unit (OTU) framework. Unlike OTUs, which cluster sequences using an arbitrary sequence similarity (ex 97%), ASVs are unique sequence reads determined using Bayesian probabilities of known sequencing error. These unique sequences can be as little as 2 bp different, providing improved taxonomic resolution and an increase in observed diversity. Please see (Callahan et al. 2016, Amir et al. 2017) for further discussion.
```

# Step 3: Taxonomic Assignment using Bowtie 2 and BLCA
```{r}
# This next module of the Anacapa toolkit assigns taxonomy to ASVs using **Bowtie 2** and a Bowtie 2 specific **Bayesian Least Common Ancestor** (**BLCA**) algorithm.
```

# Step 4: ranacapa: Data exploration
```{r}
# This portion generates ecological diversity summary statistics for a set of samples (Kandlikar et al. 2018).
# 
# The last step of the **Anacapa** Pipeline conducts exploratory data analysis to provide a first pass look at sequencing depth, taxonomic assignments, and generated data tables. This analysis is not meant for publication, but solely as a first stab at visualization of your data. This is helpful in identifying potential glaring errors or contamination, and identifying patterns worth investigating further through more robust analysis. We highly encourage data exploration before further analysis as different parameters within the **Anacapa** pipeline may produce differences in downstream results and these parameters will vary by project, stringency of taxonomic assignment, and users opinions.
```


#Crop eDNA Metadata
```{r}
#Read metadata
all.metadata <- read.csv("metadata.csv")
time.div.hours <- as.integer(all.metadata$Start.Time/100)
time.div <- as.integer(time.div.hours/6)
all.metadata$TimeGroup <- time.div
# 11:38 removes the first 3 nearshore sites
all.metadata2 <- all.metadata[11:38,]

# Order is a count of the number of reads per order per sample
order.taxa <- read.csv("corder_taxa.csv")
# Convert all N/A values of .csv file to 0
order.taxa[is.na(order.taxa)] <- 0  
# 11:38 removes the first 3 nearshore sites
order.taxa <- order.taxa[11:38,]

col.orders <- colnames(order.taxa)
order.zero.removed <- data.frame(order.taxa$Acanthuriformes)

for (i in 1:length(order.taxa)) {
  order.name <- col.orders[i]
  if (sum(order.taxa[i]) == 0) {
    
  } else {
    order.zero.removed[order.name] <- order.taxa[i]
  }
}

colSums(order.zero.removed)

```

```{r}

# PCA analysis with depth as the factor
PCA = prcomp(order.zero.removed,scale. = T) #you don't need to center and scale before this cause this function does this for you, center is T by default and you set scale. = T


PCs = PCA$x #extracting PC scores
to_plot = data.frame(PC1 = PCs[,1],
                     PC2 = PCs[,2],
                     type = as.factor(all.metadata2$DEPTH)) #put your groups here as type
ggplot(to_plot, aes(PC1, PC2, color = type)) +
  geom_point() + ggtitle("PCA Analysis of Similarity of Sample Site Order Composition, factor deptg")

to_plot2 = data.frame(PC1 = PCs[,1],
                     PC2 = PCs[,2],
                     type = as.factor(all.metadata2$TimeGroup)) #put your groups here as type
ggplot(to_plot2, aes(PC1, PC2, color = type)) +
  geom_point() + ggtitle("PCA Analysis of Similarity of Sample Site Order Composition, factor time")

library(vegan)
# Use a species accumulation curve to find the estimated total diversity
order.specaccum <- specaccum(order.taxa)

# I presume we had some logic for trimming it in this way
order.taxa.cut <- order.taxa[11:38,]
rich.adcp <- specaccum(order.taxa.cut)
rich.adcp <- rich.adcp$richness
rich.adcp <- rich.adcp[11:38]

```

```{r}
# Figure 2
# Significant relationship bewteen zooplankton diversity and eDNA diversity
order.v.class <- lm(SWDivClassUncut ~ SWDivOrderUncut, all.metadata2)
summary(order.v.class)
e <- plot_model(order.v.class, type = "pred", term = "SWDivOrderUncut") +
  labs(x = "Order diversity of eDNA samples (Shannon-Weiner diversity index)", 
       y = "Class diversity of eDNA samples (Shannon-Weiner diversity index)")
plot(e)

#Figure 4A
# Significant relationship bewteen zooplankton diversity and eDNA diversity
meter.div.sw.div <- lm(Meter.Net.Diversity ~ SWDivOrderUncut, all.metadata2)
plot(all.metadata2$SWDivOrderUncut, all.metadata2$Meter.Net.Diversity)
summary(meter.div.sw.div)
e <- plot_model(meter.div.sw.div, type = "pred", term = "SWDivOrderUncut") + 
  ggtitle("Sample Biodiversity and Meter Net Zooplankton Diversity") + 
  labs(x = "Order diversity of eDNA samples (Shannon-Weiner diversity index)", 
       y = "Meter net ooplankton diversity (Shannon-Weiner diversity index)")
plot(e)

#Figure 4B
# Significant relationship - potentially demonstrates that as biovolume increases,
# species composition becomes more uniform
order.zoop.uncut <- lm(MeterNetZooplanktonVol ~ SWDivOrderUncut, all.metadata2)
summary(order.zoop.uncut)
e <- plot_model(order.zoop.uncut, type = "pred", term = "SWDivOrderUncut") + 
  ggtitle("Sample Biodiversity and Meter Net Zooplankton Volume") + 
  labs(x = "Order diversity of eDNA samples (Shannon-Weiner diversity index)", y = "Meter net zooplankton volume (mL)")
plot(e)

#Figure 4C, Echo Amplitude, include to compare significance against Ilumina, is insignificant with removing nearshore stations
# This code lets us analyze the relationship between echo amplitude and richness
adcp.df <- data.frame("Echo Amplitude" = all.metadata2$Echo.Amplitude, "Richness" = order.specaccum$richness)
adcp.richness <- lm(Richness ~ Echo.Amplitude, adcp.df)
plot(all.metadata2$Echo.Amplitude, order.specaccum$richness)
summary(adcp.richness)
adcp.richness
g <- plot_model(adcp.richness, type = "pred", term = "Echo.Amplitude") + 
  ggtitle("Echo Amplitude and Species Richness") + 
  labs(x = "Order richness of eDNA samples", y = "Echo amplitude")
plot(g)

#Figure 5A
#Salinity and Order Richness
df.specrich.sal.order <- data.frame(order.specaccum$richness, all.metadata2$Sal..psu.)
specaccum.lm.order.sal <- lm(order.specaccum.richness ~ all.metadata2.Sal..psu., df.specrich.sal.order)
summary(specaccum.lm.order.sal)
g <- plot_model(specaccum.lm.order.sal, type = "pred", term = "all.metadata2.Sal..psu.") + ggtitle("Salinity and Order Richness") + labs(x = "Salinity (psu)", y = "Order richness")
plot(g)

#Figure 5B
#Oxygen Concentration and Order Richness
df.specrich.o2.order <- data.frame(order.specaccum$richness, all.metadata2$O2.Seapoint..mL.L.)
specaccum.lm.order.o2 <- lm(order.specaccum.richness ~ all.metadata2.O2.Seapoint..mL.L., df.specrich.o2.order)
summary(specaccum.lm.order.o2)
g <- plot_model(specaccum.lm.order.o2, type = "pred", term = "all.metadata2.O2.Seapoint..mL.L.") + ggtitle("Oxygen Concentration and Order Richness") + labs(x = "Oxygen concentration (mL/L)", y = "Order richness")
plot(g)
```

```{r}


#bring in taxonomy and meta data from working directory. 
#CO1_raw_taxonomy <- read.csv(file = "CO1_ASV_raw_taxonomy_100.csv", header = TRUE, row.names = 1)
#or
#CO1_filtered_taxonomy <- read.csv(file = "CO1_ASV_sum_by_taxonomy_100_filtered.csv", header = TRUE, row.names = 1)
library(vegan)
library(ggplot2)

correlation_metadata <- read.csv("metadata.csv")
correlation_metadata <- correlation_metadata[11:38,]
time.div.hours <- as.integer(correlation_metadata$Start.Time/100)
time.div <- as.integer(time.div.hours/6)
correlation_metadata$TimeGroup <- time.div
order.taxa <- read.csv("corder_taxa.csv")
order.taxa[is.na(order.taxa)] <- 0  # Replace NaN & Inf with NA
order.taxa.trimmed <- order.taxa[11:38,]

pca_taxonomy <- order.taxa.trimmed
pca_metadata <- correlation_metadata
order.taxa.trimmed$Depths <- correlation_metadata$DEPTH
order.taxa.trimmed$Station <- correlation_metadata$STATION.NUMBER
#coi_metadata <- read.csv("coi-metadata.csv", row.names = 1)

```

#Correspondance Analysis
```{r}
# Calculating NMDS correlation

#Converting Absolute Abundance to Relative Abundance, using decostand funciton in vegan
# Calculating relative abundance and creating new dataframe with relative abundance data
pca_taxonomy_rel <- decostand(pca_taxonomy, method = "total")

# Calculate bray distance matrix
pca_taxonomy_dist <- vegdist(pca_taxonomy_rel, method = "bray")

# Running NMDS in vegan (metaMDS), 
pca_taxonomy_NMS <-
  metaMDS(pca_taxonomy_dist,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500)


pca_taxonomy_NMS$sample.num <- pca_metadata$SAMPLE..
pca_taxonomy_NMS$depth <- pca_metadata$DEPTH
pca_taxonomy_NMS$station <- pca_metadata$STATION.NUMBER
pca_taxonomy_NMS$time <- pca_metadata$TimeGroup


plot(pca_taxonomy_NMS)
points(pca_taxonomy_NMS)
ordi.plot <- ordiplot(pca_taxonomy_NMS, ordiellipse(pca_taxonomy_NMS$points, pca_taxonomy_NMS$sample.num))

#ordicluster or ordiplot or hclust for plotting
cl <- hclust(pca_taxonomy_dist)
pca_taxonomy_NMS_df <- data.frame(pca_taxonomy_NMS$points)

NMS.plot <- plot(pca_taxonomy_NMS
     , type="n")
points(pca_taxonomy_NMS, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
text(pca_taxonomy_NMS)
summary(pca_taxonomy_NMS)

ordination <- ordicluster(pca_taxonomy_NMS, cl)
ordiplot(ordination)

factor(pca_taxonomy_NMS$points, levels = all.metadata2$SAMPLE..)

order.taxa$Depths <- all.metadata2$DEPTH
order.taxa$Station <- all.metadata2$STATION.NUMBER



#bring in taxonomy and meta data from working directory. 
#CO1_raw_taxonomy <- read.csv(file = "CO1_ASV_raw_taxonomy_100.csv", header = TRUE, row.names = 1)
#or
#CO1_filtered_taxonomy <- read.csv(file = "CO1_ASV_sum_by_taxonomy_100_filtered.csv", header = TRUE, row.names = 1)

```

```{r}
#pca_taxonomy <- order.taxa[11:38,]
#pca_metadata <- all.metadata2[11:38,]
#coi_metadata <- read.csv("coi-metadata.csv", row.names = 1)

plot(all.metadata$DEPTH, all.metadata$Start.Time)

avg.depth <- aggregate(x = all.metadata$DEPTH,
                          by = list(all.metadata$STATION.NUMBER),
                          FUN = mean)
avg.time <- aggregate(x = all.metadata$Start.Time,
                          by = list(all.metadata$STATION.NUMBER),
                          FUN = mean)

avg.depth$Time <- avg.time$x

plot(avg.depth$x, avg.depth$Time)


#Converting Absolute Abundance to Relative Abundance, using decostand funciton in vegan
# Calculating relative abundance and creating new dataframe with relative abundance data
pca_taxonomy_rel <- decostand(pca_taxonomy, method = "total")

# Calculate bray distance matrix
pca_taxonomy_dist <- vegdist(pca_taxonomy_rel, method = "bray")

# Running NMDS in vegan (metaMDS), 
pca_taxonomy_NMS <-
  metaMDS(pca_taxonomy_dist,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)


#ordicluster or ordiplot or hclust for plotting
cl <- hclust(pca_taxonomy_dist)

plot(pca_taxonomy_NMS
     , type="n")
points(pca_taxonomy_NMS, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
text(pca_taxonomy_NMS)

pca_taxonomy_NMS$points

ordiellipse(pca_taxonomy_NMS,iidx, col=1:4, kind = "ehull", lwd=3)

ordination <- ordicluster(pca_taxonomy_NMS, cl)
ordiplot(ordination)

factor(pca_taxonomy_NMS$points, levels = all.metadata2$SAMPLE..)

# PCA analysis with depth as the factor
PCA = prcomp(order.taxa,scale. = T) #you don't need to center and scale before this cause this function does this for you, center is T by default and you set scale. = T
PCs = PCA$x #extracting PC scores
to_plot = data.frame(PC1 = PCs[,1],
                     PC2 = PCs[,2],
                     type = as.factor(all.metadata2$DEPTH)) #put your groups here as type
ggplot(to_plot, aes(PC1, PC2, color = type)) +
  geom_point() + ggtitle("PCA Analysis of Similarity of Sample Site Order Composition")

# PCA analysis with time as the factor
to_plot2 = data.frame(PC1 = PCs[,1],
                     PC2 = PCs[,2],
                     type = as.factor(all.metadata2$TimeGroup)) #put your groups here as type
ggplot(to_plot2, aes(PC1, PC2, color = type)) +
  geom_point() + ggtitle("PCA Analysis of Similarity of Sample Site Order Composition by Time")

df.nmds <- data.frame(pca_taxonomy_NMS$points, 
                      type = as.factor(all.metadata2$STATION.NUMBER))
ggplot(df.nmds, aes(MDS1, MDS2, color = type)) + geom_point() + 
  ggtitle("NMDS Analysis of Similarity of Sample Site Order Composition")

```

```{r}
# PCA analysis
PCA = prcomp(order.zero.removed,scale. = T) #you don't need to center and scale before this cause this function does this for you, center is T by default and you set scale. = T
PCs = PCA$x #extracting PC scores
to_plot = data.frame(PC1 = PCs[,1],
                     PC2 = PCs[,2],
                     type = as.factor(all.metadata2$STATION.NUMBER)) #put your groups here as type
ggplot(to_plot, aes(PC1, PC2, color = type)) +
  geom_point() + ggtitle("PCA Analysis of Similarity of Sample Site Order Composition")

df.nmds <- data.frame(pca_taxonomy_NMS$points, 
                      type = as.factor(all.metadata2$STATION.NUMBER))
ggplot(df.nmds, aes(MDS1, MDS2, color = type)) + geom_point() + 
  ggtitle("NMDS Analysis of Similarity of Sample Site Order Composition")


# CCA analysis
mod.time <- cca(order.taxa.trimmed ~ TimeGroup, pca_metadata)
plot(mod.time, type="n", scaling = "symmetric") 
with(pca_metadata, ordihull(mod.time, TimeGroup, scaling = "symmetric", label = TRUE))
title("Constrained Correspondence Analysis by Time")

mod.depth <- cca(order.taxa.trimmed ~ DEPTH, pca_metadata)
plot(mod.depth, type="n", scaling = "symmetric") 
with(pca_metadata, ordihull(mod.depth, DEPTH, scaling = "symmetric", label = TRUE))
title("Constrained Correspondence Analysis by Depth")

mod.station <- cca(order.taxa.trimmed ~ STATION.NUMBER, pca_metadata)
plot(mod.station, type="n", scaling = "symmetric") 
with(pca_metadata, ordihull(mod.station, STATION.NUMBER, scaling = "symmetric", label = TRUE))
title("Constrained Correspondence Analysis by Station Number")


## use more colours and add ellipsoid hulls
plot(mod.time, type = "n")
with(pca_metadata, ordicloud(mod.time, TimeGroup, scaling = "symmetric",
                           kind = "ehull", col = 1:4, lwd=1, label = TRUE))

# Trying anova dissimilarity index
station.ano.dissim <- anosim(order.taxa.trimmed, grouping = all.metadata2$STATION.NUMBER) # Insignificant difference
depth.ano.dissim <- anosim(order.taxa.trimmed, grouping = all.metadata2$DEPTH) # Insignificant difference
time.ano.dissim <- anosim(order.taxa.trimmed, grouping = all.metadata2$TimeGroup) # Significant difference - but not super explanatory
ocean.ano.dissim <- anosim(order.taxa.trimmed, grouping = all.metadata2$OceanType) # Significant difference - not super explanatory
lat.ano.dissim <- anosim(order.taxa.trimmed, grouping = all.metadata2$LatDeg) # Insignificant difference

# Trying mrpp dissimilarity index
station.mrpp.dissim <- mrpp(order.taxa.trimmed, grouping = all.metadata2$STATION.NUMBER, distance = "bray") # Significant
station.mrpp.dissim
depth.mrpp.dissim <- mrpp(order.taxa.trimmed, grouping = all.metadata2$DEPTH, distance = "bray") # Significant
depth.mrpp.dissim
time.mrpp.dissim <- mrpp(order.taxa.trimmed, grouping = all.metadata2$TimeGroup, distance = "bray") # Significant 
time.mrpp.dissim
ocean.mrpp.dissim <- mrpp(order.taxa.trimmed, grouping = all.metadata2$OceanType, distance = "bray") # Insignificant
ocean.mrpp.dissim
lat.mrpp.dissim <- mrpp(order.taxa.trimmed, grouping = all.metadata2$LatDeg, distance = "bray") # Significant
lat.mrpp.dissim

# Trying mrpp dissimilarity index with Jaccard 
# (see https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/edn3.239 for suggestion)
station.mrpp.dissim <- mrpp(order.taxa.trimmed, grouping = all.metadata2$STATION.NUMBER, distance = "jaccard") # Significant
station.mrpp.dissim
depth.mrpp.dissim <- mrpp(order.taxa.trimmed, grouping = all.metadata2$DEPTH, distance = "jaccard") # Significant
depth.mrpp.dissim
time.mrpp.dissim <- mrpp(order.taxa.trimmed, grouping = all.metadata2$TimeGroup, distance = "jaccard") # Insignificant 
time.mrpp.dissim
ocean.mrpp.dissim <- mrpp(order.taxa.trimmed, grouping = all.metadata2$OceanType, distance = "jaccard") # Insignificant
ocean.mrpp.dissim
lat.mrpp.dissim <- mrpp(order.taxa.trimmed, grouping = all.metadata2$LatDeg, distance = "jaccard") # Significant
lat.mrpp.dissim

# Tried for DNA again
dna.mrpp.dissim <- mrpp(order.taxa.trimmed, grouping = all.metadata2$DNA.extracted., distance = "bray") # Insignificant
dna.mrpp.dissim

# Trying adonis dissimilarity index
station.adonis.dissim <- adonis2(order.taxa.trimmed ~ STATION.NUMBER, all.metadata2) # Residual explains more than station num
depth.adonis.dissim <- adonis2(order.taxa.trimmed ~ DEPTH, all.metadata2) # Depth has strongly explanatory power than residual
time.adonis.dissim <- adonis2(order.taxa.trimmed ~ TimeGroup, all.metadata2) # Time group has poor (but mildly significant) explanatory power
ocean.adonis.dissim <- adonis2(order.taxa.trimmed ~ OceanType, all.metadata2) # Insignificant explanatory poer
lat.adonis.dissim <- adonis2(order.taxa.trimmed ~ LatDeg, all.metadata2) # Insignificant explanatory poer
lat.adonis.dissim

# Checking our data to see if we did it right lmao
lib.ano.dissim <- anosim(order.taxa.trimmed, grouping = all.metadata2$Library.No) # Insignificant difference
pcr.ano.dissim <- anosim(order.taxa.trimmed, grouping = all.metadata2$PCR.No.) # Insignificant difference
dna.ano.dissim <- anosim(order.taxa.trimmed, grouping = all.metadata2$DNA.extracted.) # Significant difference
niskin.ano.dissim <- anosim(order.taxa.trimmed, grouping = all.metadata2$Niskin.Bottle) # Significant difference

ocean.type <- c("Cyclonic Eddy",         "Cyclonic Eddy",         "Cyclonic Eddy",         "Cyclonic Eddy",        
"Eddy convergence",      "Eddy convergence",      "Anticyclonic Eddy",    "Anticyclonic Eddy",   
"Anticyclonic eddy",    "Anticyclonic eddy",    "Anticyclonic eddy",     "Anticyclonic eddy",    
"Anticyclonic eddy",     "Anticyclonic eddy",     "Edge of Cyclonic Eddy", "Edge of Cyclonic Eddy",
"Edge of Cyclonic Eddy", "Edge of Cyclonic Eddy", "Cyclonic Eddy",         "Cyclonic Eddy",        
"Subtropical Gyre",      "Subtropical Gyre", "Open Ocean", "Open Ocean",                                                 
"OMZ Compression",       "OMZ Compression",       "Divergence",            "Divergence") 
all.metadata2$OceanType <- ocean.type

mod.ocean <- cca(order.taxa.trimmed ~ ocean.type, pca_metadata)
plot(mod.ocean, type="n", scaling = "symmetric") 
with(pca_metadata, ordihull(mod.ocean, ocean.type, scaling = "symmetric", label = TRUE))
title("Constrained Correspondence Analysis by Ocean Type")

```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
