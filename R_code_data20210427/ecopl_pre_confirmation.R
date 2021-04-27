#Statistical recipe for quantifying microbial functional diversity from EcoPlate metabolic profiling
#Takeshi Miki, Taichi Yokokawa, Po-Ju Ke, I Fang Hsieh, Chih-hao Hsieh, Tomonori Kume, Kinuyo Yoney, Kazuaki Matsui
#kmatsui@civileng.kindai.ac.jp


sessionInfo()
install.packages("vegan")
library(vegan)

sample_EC <-read.csv("sample.csv", header=FALSE, skip=0) #sample ecoplate color paatten (Functional matrix E_C)
format_d<-read.csv("format_xitou_pattern.csv", header=T) #sample metadata

#Prepare the list of consistance name of chemicals included in ecoplate
chem_name<-c("Pyruvic-Acid-Methyl-Ester", "Tween-40", "Tween-80","alpha-Cyclodextrin", "Glycogen", "D-Cellobiose","alpha-D-Lactose", "beta-Methyl-D-Glucoside", "D-Xylose", "i-Erythritol", "D-Mannitol","N-Acetyl-D-Glucosamine", "D-Glucosaminic-Acid", "Glucose-1-Phosphate","alpha-Glycerol-Phosphate","D-Galactonic-Acid-gamma-Lactone", "D-Galacturonic-Acid", "2-Hydroxy-Benzoic-Acid", "4-Hydroxy-Benzoic-Acid", "gamma-Hydroxybutyric-Acid", "Itaconic-Acid", "alpha-Ketobutyric-Acid", "D-Malic-Acid", "L-Arginine", "L-Asparagine", "L-Phenylalanine", "L-Serine", "L-Threonine", "Glycyl-L-Glutamic-Acid", "Phenylethyl-amine", "Putrescine")

######The step to prepare Functional matrix in Figure 3
sample_EC<-sample_EC[,c(-1)] #exclusing water data (control, zero)
row.names(sample_EC)<-format_d$sample  #Add name (sample ID)
colnames(sample_EC)<-chem_name #Add name (substrate name)
View(sample_EC)  #-->shoud check the worksheet "sample_EC" of "sample_pre_data.xlsx".


######The step to obdain 31 binary values in Figure 3
#Convert data into binary values by quantile (quantile-based multifunctionality)
ebc<-list()  #prepare the output list
thres = 0.7 #assing the value to the threshold T
for(i in 1:31) ebc[[i]] <-(sample_EC[,i] > quantile(sample_EC[,i], thres))  # return true if the color value qunatile is greater than threshold
E_BC<-as.data.frame(ebc)  #converted into dataframe
row.names(E_BC)<-format_d$sample
colnames(E_BC)<-chem_name
E_BC[E_BC==TRUE]<-1  #converting T,F to 1 and 0; ts is now dataframe to represent the binary matrix (E_BC in Figure 3)
View(E_BC) #-->should check the worksheet "sample_EB_C" of "sample_pre_data.xlsx"

######The step to calculate functional diversity in Figure 3
mf<-list()
for(i in 1:24) mf[[i]]<-sum(E_BC[i,])
MF<-t(as.data.frame(mf))
row.names(MF)<-format_d$sample
View(MF) #-->should check the worksheet "MF" of "sample_pre_data.xlsx"

######The step to calculate the functional dissimilarity matrix DF from EC
#Calculating the default continuous dissimilarity (without chemical-simirality information) [baseline of the analysis]
sample_EC[sample_EC <0]<-0 #convert negative values to zero, which is necessary to calculate the bray-curtis distance
DF_UW<-vegdist(sample_EC, method="bray")
View(as.matrix(DF_UW)) #-->shoud check the whorksheet "DF_UW" of "sample_pre_data.xlsx"

######The step to calculate the functional dissimilarity matrix DF from EBC
#Calculating the default continuous dissimilarity (without chemical-simirality information) [baseline of the analysis]
DF_BUW<-vegdist(E_BC, method="jaccard")
View(as.matrix(DF_BUW)) #-->shoud check the whorksheet "DF_BUW" of "sample_pre_data.xlsx"


#write.csv(as.matrix(DF_BUW), "sample4.csv", quote=FALSE, row.names=TRUE)

