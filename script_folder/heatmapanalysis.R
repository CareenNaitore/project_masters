##################SIGNIFICANTLY EXPRESSED GENES HEATMAP VISUAL PRESENTATION####################
setwd("~/DE_analysis/smple/")

library(gplots)
library(RColorBrewer)
library(ggplot2)
library(reshape2)


data<-read.csv("normalized_readcounts.csv",header = T, sep=",")
head(data)
desc<-(data[,1])
#cbind the row mean of the average of the replicates in one data3  
data3<- cbind(desc,female=c(rowMeans(norms[1:3])),male=c(rowMeans(norms[4:6])),Larvae=c(rowMeans(norms[7:9])),Gravidfemale=c(rowMeans(norms[10:12])),Pupae=c(rowMeans(norms[13:15])),Teneral_females=c(rowMeans(norms[16:18])),Teneral_males=c(rowMeans(norms[19:21])), Wildtype_males=c(rowMeans(norms[22:24])))
View(data3)
as.data.frame(data3)
as.data.frame(desc)
head(desc)
data<-data[,2:length(data)]
head(data)
data_matrix<-as.matrix(data)

norms=data
head(norms)
hdata<-log2(norms)
head(hdata)
hdata<-cbind(desc,hdata)

head(hdata)
write.csv(hdata, file = "lifestage4heatmap.csv")

data <- read.csv("lifestage4heatmap.csv")
head(data)
rnames <- data[,2]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
head(mat_data)
rownames(mat_data) <- rnames  # assign row names
head(mat_data)
mat_data[mat_data==-Inf]=-2
head(mat_data)


sapply(mat_data, class)
dtx<-is.na(mat_data)
# creates a own color palette from red to green
my_palette <- colorRampPalette(c("blue1", "floralwhite", "red"))(n = 299)
ColAnn <- data.frame(colnames(heat))
colnames(ColAnn) <- c("Sample")
ColAnn <- HeatmapAnnotation(df=ColAnn, which="col")

RowAnn <- data.frame(df$Family)
colnames(RowAnn) <- c("Gene family")
colours <- list("Gene family"=c("ncRNA"="royalblue","pseudogene"="red3"))
RowAnn <- HeatmapAnnotation(df=RowAnn, col=colours, which="row")

# (optional) defines the color breaks manually for a "skewed" color transition

col_breaks = c(seq(-0.3,3,length=100),  # for purple
               seq(4,9,length=100),    # for white
               seq(10,18,length=100))    # for yellow/ non overlapping breaks!!


# creates a 5 x 5 inch image
png("heatmap3.png",    # create PNG for the heat map        
    width = 15*300,        # 15 x 300 pixels
    height = 12*300,
    res = 300,            # 300 pixels per inch
    pointsize =15)        # smaller font size


heatmap.2(mat_data,cexRow=0.3,
          main = "miRNA in developmental stages ",# heat map title
          # change font color of cell labels to black
          
          key = T,
          keysize = 1.0,
          #xlab = "sample names", 
          density.info="histogram",  # turns off density plot inside color legend(none or histogram)
          trace="none",         # turns off trace lines inside the heat map (none, column, row or both)
          margins =c(12,20),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram, put both to have row and column
          Colv="Rowv")      # turn off column clustering, put RowV to have both


dev.off()               # close the PNG device

