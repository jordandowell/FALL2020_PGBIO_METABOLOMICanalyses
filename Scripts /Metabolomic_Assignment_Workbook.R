#assignment 
#metabolomics workshop 
#install and library the below packages



#dataset citation:
#Zytynska, Sharon E. et al. (2019), 
#Data from: Effect of plant chemical variation and mutualistic ants on the local population genetic structure of an aphid herbivore, Dryad, Dataset, 
#https://doi.org/10.5061/dryad.mm7bj56

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
metabolome<- read.csv("Data/Assignment_Dataset.csv",
                      header = T,
                      row.names = 1)
#data engineering
View(metabolome)
#change Naphids to a factor
metabolome$Naphids<-as.factor(metabolome$Naphids)
 


#remove the factor so we have a data set of just the compounds
META<-(metabolome[-1])
View(META)

#table will give you a list of how many individuals in each factor
#in this case we know that the first 145 lines are 1 aphid
#next 103 lines are 2 aphids etc. 
table(metabolome[1])
##Construct Volcano plots based on t-test & logfold change
#transpose data
VOLCANO_META<-t(META)
View(VOLCANO_META)
# run t.test for the first compound  [row,sequence of columns]
t.test(VOLCANO_META[1, 1:145], VOLCANO_META[1, 146:209])
#Run a t.test for all compounds between 2 groups 
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
                   grp1 = c(1:145),
                   grp2 = c(146:209))


# a multiple testing correction would be the alpha(0.05)/number of tests(e.g. number of compounds 170)
#how many of our compounds have a pvalue less than our alpha (alpha/number of compounds)
sum(rawpvalue <= (0.05/22))


######   Question 1: How many compounds have a pvalue less than the bonferroni correction and less than the alpha?






##This data has already been log transformed so we do not need to transform it.
VOLCANO_META_log <- VOLCANO_META

#calculate the mean of each compound per control group (the group you are comparing a "treatment to" )
WILD <- apply(VOLCANO_META_log[, 1:145], 1, mean)

#calulate the standard deviation of each compound per control group
WILD_SD<-apply(VOLCANO_META_log[, 1:145], 1, sd)
#calculate the standared error Standard deviaiton/ sqrt(N samples)
WILD_SE<-WILD_SD/sqrt(145)

#set the fold change cutoff as 2* the average standard error 
#this is also the average 95% confidence interval 
hist(WILD_SE)
FoldchangeCutoff<-2*mean(WILD_SE)

#calcuate the mean of each compound per test group
SINGLE_KO <-apply(VOLCANO_META_log[, 146:209], 1, mean)



#calculate the fold change 
WILD_SINGLE_FOLDCHANGE <- WILD - SINGLE_KO
#View a histogram of fold changes
# you can change the title of the x lable in the quotations
hist(WILD_SINGLE_FOLDCHANGE, xlab = "log2 Fold Change (WILD vs SINGLE_KO)")


#create a results table

results <- cbind(WILD_SINGLE_FOLDCHANGE, rawpvalue)
results <- as.data.frame(results)
results$probename <- rownames(results)

View(results)
#visualize with a volcano plot 
#change titles to represent you analysis
dev.off()
DifferentialMETS <-
  EnhancedVolcano(
    toptable = results,
    x = "WILD_SINGLE_FOLDCHANGE",
    y = "rawpvalue",
    lab = rownames(results),
    xlab = bquote( ~ Log[10] ~ "fold change"),
    #pvalue can be altered theres a few number of compounds so multiple testing correction is not truly needed
    pCutoff = 0.05 / length(results[,1]),
    #you can change the fold change cut off to what you think is importnt 
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
#this plot will need to be turned in 
#write a summary including the identity of any compounds above the foldchange threshold & below the p-value cutoff 

#######Question 2: Which compounds are above the fold change  & below the p-value thresholds? (plot must be included)





#PCA to visualize clusters



#PCA run
#for ease of examination reduce the dataset to the two groups you wish to compare. remember [rows to include, columns to include ] 
#include all the compound columns so you should only need to change the rows to include
View(metabolome)
table(metabolome$Naphids)


RIDA.pca <- prcomp(metabolome[1:nrow(metabolome) , 2:ncol(metabolome)],
                   center = TRUE,
                   scale. = TRUE)

#assess the scree plot 
#Visually how many components?
#save the following plot and describe how many axes of variaiton should be explored based on the eigenvalues
plot(RIDA.pca$sdev^2, ylab = "Eigenvalues")
#add a line to make the traditional cutoff of 1 
abline(h=1)


#assess components quatitatively
summary(RIDA.pca)

#how much cumulative variaiton is explained by the number of axes which should be explored


####### Question 3: How many components should be explored? how much cumulative variance is explained by the components that should be explored? (must include plot)



 ####### Question 4: Using the following bits of code identify which axes best separate out your group/s of interest and identify which compounds are separating your group 

#explore Pertinent PCs that can separate your groups of interest by changing choices 

#PC  change choices 1&2 to the axes of interest
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

#investigate what causes the differences between groups along your PC or PCs of interest 
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
#change != to == if you want to subset the compounds going in opposite directions 
Directional.compounds<-subset(Directional.compounds,  rowSums(sign(Directional.compounds))!=0)

VIP<- as.data.frame(Directional.compounds[,1]+Directional.compounds[,2])
Directional.compounds<-cbind(Directional.compounds,VIP)
Directional.compounds<- Directional.compounds[order(Directional.compounds[,3]),]


#View compounds ordered by variable importance
View(Directional.compounds)

#get the top 5 contributing compounds
View(head(Directional.compounds))
View(tail(Directional.compounds))



# summary should be submitted in a word-like document with figures and descriptions



