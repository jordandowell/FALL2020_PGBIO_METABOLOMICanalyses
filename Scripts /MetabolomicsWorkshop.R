
#metabolomics workshop 
#install and library the below packages


if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot")
library(EnhancedVolcano)
library(ggbiplot)








#check proper working directory

#import data 
metabolome<- read.csv("Data/S.enterica_metabolome_rida.csv",
                      header = T,
                      row.names = 1)
#data engineering
View(metabolome)
META<-(metabolome[-1])
View(META)

##Construct Volcano plots based on t-test & logfold change
#transpose data
VOLCANO_META<-t(META)
View(VOLCANO_META)
# lets run t.test for xylulose the first compound 
t.test(VOLCANO_META[1, 1:6], VOLCANO_META[1, 7:12])
#lets run a t.test for all compounds between wild type and single KO
ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)
  results = t.test(x, y)
  results$p.value
}
#get the pvalue for each compound 
rawpvalue <- apply(VOLCANO_META,
                  1,
                  ttestRat,
                  grp1 = c(1:6),
                  grp2 = c(7:12))


# a multiple testing correction would be the alpha(0.05)/number of tests(e.g. number of compounds 170)
#how many of our compounds have a pvalue less than our alpha 
sum(rawpvalue <= (0.05/170))



##transform our data into log10 base.
VOLCANO_META_log = log(VOLCANO_META)

#calculate the mean of each compound per control group
WILD <- apply(VOLCANO_META_log[, 1:6], 1, mean)

#calulate the standard deviation of each compound per control group
WILD_SD<-apply(VOLCANO_META_log[, 1:6], 1, sd)
WILD_SE<-WILD_SD/sqrt(6)

#set the fold change cutoff as 2* the average standard error 
#this is also the average 95% confidence interval 
hist(WILD_SE)
FoldchangeCutoff<-2*mean(WILD_SE)

#calcuate the mean of each compound per test group
SINGLE_KO <-apply(VOLCANO_META_log[, 7:12], 1, mean)




WILD_SINGLE_FOLDCHANGE <- WILD - SINGLE_KO
#View a histogram of fold changes
hist(WILD_SINGLE_FOLDCHANGE, xlab = "log2 Fold Change (WILD vs SINGLE_KO)")


#create a results table

results <- cbind(WILD_SINGLE_FOLDCHANGE, rawpvalue)
results <- as.data.frame(results)
results$probename <- rownames(results)

View(results)
#visualize with a volcano plot
dev.off()
DifferentialMETS <-
  EnhancedVolcano(
    toptable = results,
    x = "WILD_SINGLE_FOLDCHANGE",
    y = "rawpvalue",
    lab = rownames(results),
    xlab = bquote( ~ Log[10] ~ "fold change"),
    pCutoff = 0.05 / length(results[,1]),
    FCcutoff = FoldchangeCutoff,
    xlim = c(min(results[, 1], na.rm = TRUE),
             max(results[, 1], na.rm = TRUE) + 1),
    ylim = c(0, max(-log10(results[, 2]), na.rm = TRUE) + 1),
    title = 'Volcano plot',
    subtitle = 'Differentially Produced Metabolites in single RidA KO',
    caption = paste0('Total = ', nrow(results) + 1, ' metabolites'),
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    # boxedlabels = TRUE,
    legendLabels = c("NS", "Log10 FC", "P", "P & Log10 FC"),
    legendPosition = 'right',
    legendLabSize = 12,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.2,
    colConnectors = 'grey30'
  )


#View plot
DifferentialMETS




#PCA to visualize clusters

View(metabolome)

#PCA run
RIDA.pca <- prcomp(metabolome[, 2:ncol(metabolome)],
                   center = TRUE,
                   scale. = TRUE)

#assess the scree plot  look for inflection point.
#Visually how many components? As a rule of thumb any axis with an eigenvalue > 1 should be explored 
#eigenvalues are equal to the standard deviation of an axis squared
plot(RIDA.pca$sdev^2, ylab = "Eigenvalues")
#add a line to make the traditional cutoff of 1 
abline(h=1)

#assess components quatitatively
summary(RIDA.pca)

#PC 1&2
ggbiplot(
  RIDA.pca,
  choices = c(1, 2),
  ellipse.prob = 0.95,
  obs.scale = 1,
  var.scale = 1,
  groups = metabolome$treatment,
  labels = rownames(metabolome),
  var.axes = FALSE,
  ellipse = TRUE,
  circle = FALSE
) +
  scale_color_discrete(name = '') +
  theme_bw()+theme(legend.direction = 'horizontal', legend.position = 'top')

#lets investigate what cases the differences between groups using PC 2
#rotation gives us the orientation of individual compounds
View(RIDA.pca$rotation[,2])

PC2.Compounds<- RIDA.pca$rotation[,2][order(RIDA.pca$rotation[,2],decreasing = F)]

#lets look at the compounds relating to this separation 
View(PC2.Compounds[1:6])
View(tail(PC2.Compounds))



#lets assess another group of PCs to see if we can find a dosage effect 
#PC 3&4
ggbiplot(
  RIDA.pca,
  choices = c(3, 4),
  obs.scale = 1,
  var.scale = 1,
  ellipse.prob = 0.95,
  groups = metabolome$treatment,
  labels = rownames(metabolome),
  var.axes = FALSE,
  ellipse = TRUE,
  circle = FALSE
) +
  scale_color_discrete(name = '') +
  theme_bw()+theme(legend.direction = 'horizontal', legend.position = 'top')

# shows that there is a dosage effect, however not a major contributor of variation 

#lets investigate what cases the differences between groups using pcs 3 &4
#rotation gives us the orientation of individual compounds
View(RIDA.pca$rotation[,3:4])

#lets reduce this list to the compounds with equal directionality( e.g. both positve or both negative )

Directional.compounds<- RIDA.pca$rotation[,3:4]
Directional.compounds<-subset(Directional.compounds,  rowSums(sign(Directional.compounds))!=0)
#calculate variable importance for the joint vector 
VIP<- as.data.frame(Directional.compounds[,1]+Directional.compounds[,2])
Directional.compounds<-cbind(Directional.compounds,VIP)
Directional.compounds<- Directional.compounds[order(Directional.compounds[,3]),]


#View compounds ordered by variable importance
View(Directional.compounds)

#get the top 5 contributing compounds
View(head(Directional.compounds))
View(tail(Directional.compounds))


#lets investigate what cases the differences between groups using pcs 3 &4
#lets reduce this list to the compounds with equal directionality( e.g. both positve or both negative )

Directional.compounds<- RIDA.pca$rotation[,3:4]
Directional.compounds<-subset(Directional.compounds,  rowSums(sign(Directional.compounds))!=0)

VIP<- as.data.frame(Directional.compounds[,1]+Directional.compounds[,2])
Directional.compounds<-cbind(Directional.compounds,VIP)
Directional.compounds<- Directional.compounds[order(Directional.compounds[,3]),]


#View compounds ordered by variable importance
View(Directional.compounds)

#get the top 5 contributing compounds
View(head(Directional.compounds))
View(tail(Directional.compounds))





