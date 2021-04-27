#Statistical recipe for quantifying microbial functional diversity from EcoPlate metabolic profiling
#Takeshi Miki, Taichi Yokokawa, Po-Ju Ke, I Fang Hsieh, Chih-hao Hsieh, Tomonori Kume, Kinuyo Yoney, Kazuaki Matsui
#kmatsui@civileng.kindai.ac.jp

###########################PART 0 Setting Environment################################################
sessionInfo()
version()


###Support > MacOSX 10.7.3, > Windows 7, Ubuntu14 
###Results from Java packages should be saved for the case when Java packages do not work

install.packages("rJava")
install.packages("rcdk")

install.packages("rJava", type='source')
install.packages("rcdk", type='source')

#install.packages("./library/rJava_0.9-8.tar.gz", repos = NULL, type = "source")
#install.packages("./library/rcdk_3.3.8.tar.gz", repos = NULL, type = "source")


install.packages("vegan")
ninstall.packages("picante")
install.packages("devetools")
devtools::install_github("hoxo-m/pforeach")  #This is necessary one for installing pforeach (library for parallel computing)
install.packages("GUniFrac")
install.packages("ggplot2")
install.packages("labdsv")
install.packages("SYNCSA")
##########################End of PART 0###########################################################

################################PART 1 Analysis of Chemical Similarity###############################################################################
{
library("rJava")    #load this library for checking the compatibility between R and Java
library("rcdk")     #library for java chemical analysis

##NOTE: there emerges the incompatibility between rcdk and Java installed in your OS. The potential solution is:
#i) check the version of java 
#>java -version
#If the version shown is <= 1.7, then update your java (via command line in ubuntu or visit Oracle website for Win/OSX)
#Even if the version show is > 1.8, the R tells that the version was older. In this case, please check /usr/lib/jvm/ via ls -la. Then the symbolic link (e.g., default-java) may link to the older version of Java. In this case you need to update the symbolic link
#ii) After updating the java, first make R recognize the updated vrsion of Java via >sudo R CMD javareconf, and then you need to install library rJava from R toolbar.
#iii) install Chemistry Development Kit (CDK) Java libraries via Ubuntu Software Center or >sudo apt-get install libcdk-java
#iv) Try reinstall rcdk.
#v) It would be better to reboot your PC/MAC/WS or Rstudio-server

library("picante")  #library for generating phylogenetic trees of which functions are used for generating chemical similarity trees.


#Load the list of 2D chemical structure of 31 carbon substrates in ecoplate
mol_list<-load.molecules("./sdf/list_ecoplate2.sdf")
#show the graphics of 2D structure of these chemicals 
view.molecule.2d(mol_list)

#Preparation for chemical similarity tree-X (X=a,b,or c)
#For tree-a & tree-b
fp.list_a <- lapply(mol_list, get.fingerprint)
fp.list_b <- lapply(mol_list, get.fingerprint, type="extended") #see 2771.pdf for the comparison of different fingerprinting


fp.dist_a<-fp.sim.matrix(fp.list_a, method='tanimoto')
fp.dist_b<-fp.sim.matrix(fp.list_b, method='tanimoto')

fp.dist_a<-as.data.frame(1-fp.dist_a)  #from Similarity to dissimilarity
fp.dist_b<-as.data.frame(1-fp.dist_b)  #from Similarity to dissimilarity


#write.csv(fp.dist_a, "./fp.dist_a.csv", row.names=F, quote=F)
#write.csv(fp.dist_b, "./fp.dist_b.csv", row.names=F, quote=F)

#When rJava and/or rcdk installation were failed
fp.dist_a<-read.csv("./fp.dist_a.csv")
fp.dist_b<-read.csv("./fp.dist_b.csv")

#Prepare the list of consistance name of chemicals included in ecoplate
chem_name<-c("Pyruvic-Acid-Methyl-Ester", "Tween-40", "Tween-80","alpha-Cyclodextrin", "Glycogen", "D-Cellobiose","alpha-D-Lactose", "beta-Methyl-D-Glucoside", "D-Xylose", "i-Erythritol", "D-Mannitol","N-Acetyl-D-Glucosamine", "D-Glucosaminic-Acid", "Glucose-1-Phosphate","alpha-Glycerol-Phosphate","D-Galactonic-Acid-gamma-Lactone", "D-Galacturonic-Acid", "2-Hydroxy-Benzoic-Acid", "4-Hydroxy-Benzoic-Acid", "gamma-Hydroxybutyric-Acid", "Itaconic-Acid", "alpha-Ketobutyric-Acid", "D-Malic-Acid", "L-Arginine", "L-Asparagine", "L-Phenylalanine", "L-Serine", "L-Threonine", "Glycyl-L-Glutamic-Acid", "Phenylethyl-amine", "Putrescine")

colnames(fp.dist_a)<-chem_name  #Add chemical labels first
colnames(fp.dist_b)<-chem_name  #Add chemical labels first 
fp.dist_a<-as.dist(fp.dist_a)   #converting to the format distance
fp.dist_b<-as.dist(fp.dist_b)   #converting to the format distance

#For tree-c From Online tool  http://chemmine.ucr.edu/
online_mat<-read.table("./online.dissimilarity")
fp.dist_c<-as.dist(online_mat)  #converting to the format distance

#Processes to generate the clustring tree
#We have concluded the method ="average" gives the best correlation between distance matrix and tree-shape (excute the next block)
#Check the difference of the clustring method : https://www1.doshisha.ac.jp/~mjin/R/28/28.html

clus.hier_a<-hclust(fp.dist_a, method="average")
plot(clus.hier_a, hang=-1, cex=1.0)
tree_a.ph<-as.phylo(clus.hier_a)  #convert to phytogenetic tree
plot(tree_a.ph)   ###Used for Fig.S1a

clus.hier_b<-hclust(fp.dist_b, method="average")
plot(clus.hier_b, hang=-1, cex=1.0)
tree_b.ph<-as.phylo(clus.hier_b)  #convert to phytogenetic tree
plot(tree_b.ph)   ###Used for Fig.S1b

clus.hier_c<-hclust(fp.dist_c, method="average")
plot(clus.hier_c, cex=1.0)
tree_c.ph<-as.phylo(clus.hier_c)
plot(tree_c.ph) ###Used for Fig.S1c

#choosing the best clustering by cophenetic correlation (see Petchey OL and Gaston KJ (2006) Functional diversity:back to basics and looking forward. Ecology Letters 9: 741-758)
{
#With different methods from the identical input (distance matrix), different tree shapes are obtained. Of-course, the some fraction of information from the distance matrix was lost or distorted during the generation of the clustering tree. Therefore, it is important to choose the method with which the distortion is at minimum

#For tree-a, generate the tree by different methods
clus_1<-hclust(fp.dist_a, method="ward.D")
clus_2<-hclust(fp.dist_a, method="ward.D2")
clus_3<-hclust(fp.dist_a, method="single")
clus_4<-hclust(fp.dist_a, method="complete")
clus_5<-hclust(fp.dist_a, method="average") #this returns the highest mantel r
clus_6<-hclust(fp.dist_a, method="mcquitty")
clus_7<-hclust(fp.dist_a, method="median")
clus_8<-hclust(fp.dist_a, method="centroid")

#Calculate the distance between nodes of the tree
coph_dist_1a<-cophenetic(clus_1)
coph_dist_2a<-cophenetic(clus_2)
coph_dist_3a<-cophenetic(clus_3)
coph_dist_4a<-cophenetic(clus_4)
coph_dist_5a<-cophenetic(clus_5)
coph_dist_6a<-cophenetic(clus_6)
coph_dist_7a<-cophenetic(clus_7)
coph_dist_8a<-cophenetic(clus_8)

#Mantel correlation test: original distance matrix (fp.dist_a) vs distance matrix generated from the tree
mantel(fp.dist_a, coph_dist_1a)
mantel(fp.dist_a, coph_dist_2a)
mantel(fp.dist_a, coph_dist_3a)
mantel(fp.dist_a, coph_dist_4a)
mantel(fp.dist_a, coph_dist_5a)  #this returns the highest mantel r
mantel(fp.dist_a, coph_dist_6a)
mantel(fp.dist_a, coph_dist_7a)
mantel(fp.dist_a, coph_dist_8a)

#For tree-b
clus_1<-hclust(fp.dist_b, method="ward.D")
clus_2<-hclust(fp.dist_b, method="ward.D2")
clus_3<-hclust(fp.dist_b, method="single")
clus_4<-hclust(fp.dist_b, method="complete")
clus_5<-hclust(fp.dist_b, method="average") #this returns the highest mantel r
clus_6<-hclust(fp.dist_b, method="mcquitty")
clus_7<-hclust(fp.dist_b, method="median")
clus_8<-hclust(fp.dist_b, method="centroid")

coph_dist_1b<-cophenetic(clus_1)
coph_dist_2b<-cophenetic(clus_2)
coph_dist_3b<-cophenetic(clus_3)
coph_dist_4b<-cophenetic(clus_4)
coph_dist_5b<-cophenetic(clus_5)
coph_dist_6b<-cophenetic(clus_6)
coph_dist_7b<-cophenetic(clus_7)
coph_dist_8b<-cophenetic(clus_8)

mantel(fp.dist_b, coph_dist_1b)
mantel(fp.dist_b, coph_dist_2b)
mantel(fp.dist_b, coph_dist_3b)
mantel(fp.dist_b, coph_dist_4b)
mantel(fp.dist_b, coph_dist_5b)  #this returns the highest r
mantel(fp.dist_b, coph_dist_6b)
mantel(fp.dist_b, coph_dist_7b)
mantel(fp.dist_b, coph_dist_8b)

#For tree-c
clus_1<-hclust(fp.dist_c, method="ward.D")
clus_2<-hclust(fp.dist_c, method="ward.D2")
clus_3<-hclust(fp.dist_c, method="single")
clus_4<-hclust(fp.dist_c, method="complete")
clus_5<-hclust(fp.dist_c, method="average") #this returns the highest r = UPGMA
clus_6<-hclust(fp.dist_c, method="mcquitty")
clus_7<-hclust(fp.dist_c, method="median")
clus_8<-hclust(fp.dist_c, method="centroid")

coph_dist_1c<-cophenetic(clus_1)
coph_dist_2c<-cophenetic(clus_2)
coph_dist_3c<-cophenetic(clus_3)
coph_dist_4c<-cophenetic(clus_4)
coph_dist_5c<-cophenetic(clus_5)
coph_dist_6c<-cophenetic(clus_6)
coph_dist_7c<-cophenetic(clus_7)
coph_dist_8c<-cophenetic(clus_8)

mantel(fp.dist_c, coph_dist_1c)
mantel(fp.dist_c, coph_dist_2c)
mantel(fp.dist_c, coph_dist_3c)
mantel(fp.dist_c, coph_dist_4c)
mantel(fp.dist_c, coph_dist_5c) #this returns the highest r
mantel(fp.dist_c, coph_dist_6c)
mantel(fp.dist_c, coph_dist_7c)
mantel(fp.dist_c, coph_dist_8c)

}

############Trial to generate random-tree######
#This will be used for statistical significance test by permutating trees
leng<-length(fp.dist_c)
z<-sample(1:leng, leng, replace=FALSE)
temp_d<-fp.dist_c
for(i in 1:leng) temp_d[i]<-fp.dist_c[z[i]]
#temp_d
temp_chem<-hclust(temp_d, method="average")
temp.ph<-as.phylo(temp_chem)
plot(temp.ph)
###############################################


###This is important part to compare different trees!! (For Supporting Information)####################
#Correlation between three trees
mantel(coph_dist_5a,coph_dist_5b) #Tree-a vs Tree-b
mantel(coph_dist_5a,coph_dist_5c) #Tree-a vs Tree-c
mantel(coph_dist_5b,coph_dist_5c) #Tree-b vs Tree-c

}#end of PART1
#######################End of PART1###########################################################


##############PART2 Ecopalte data loading ecoplate analysis######################
library("vegan")


####################Definitions of functions#########################
#[1] Definition of functions to load data from Xitou Ecoplate and those from the microcosm experiments
#This part should be modified depending of the format of data files in your computer

#This is an example of the function to load data:
#Setting: let's say folder A includes this script file & subfolder text_file, and the subforder text_file includes ecoplate data from multiple sampling months but the data from different sampling months are saved in different folders. 
#Here, this function specifies the relative path of the folder posision by the first parameter: relative_path1, and 2) the subfolder within it by the 2nd paramter: relative_path2. 
#Other parameters are set, depending on the file name of ecoplate data. Here, the file name is e.g. "20150228Eco_05_T1.TXT" where "20150228" is the sampling date (=innoculation date), "Eco" is the type of biolog (Ecoplate), "05" means the measurement day after innoculation date, and "N1" is the site name (=T) with replication ID (=1). The parameters: ini_date represents the sampling date, plate_id represent the site name with replication ID, ini_no1 is the first day of the incubation, and end_no2 is the final day of the incubation. The reason why we also have other two parameters (end_no1, ini_no2) is we have different format for the file name from 10th day. 
#The final parameter (no_skip) is necessary because the raw text data from microplate reader includes the header information (e.g. date information), which should be skipped. In our format, the first 5 lines including the empty lines should be skipped for data loading.

load_ecoplate_data_demo01 <-  function(relative_path1="./text_file/", relative_path2="2015_Feb", ini_date=20150228, ini_no1=1, end_no1=9, ini_no2=10, end_no2=20, plate_id="T1", no_skip=5)    #The parameters that already have values in the functional definition are recognized as default setting.
{
  
  data_list <-list() #generate empty list, ecoplate data are saved as a list
  
  for(i in ini_no1:end_no1) {
    file_name <- paste(relative_path1,relative_path2,"/", as.character(ini_date), "Eco_0", as.character(as.integer(i)), "_", plate_id, ".TXT",sep="")   #the function past is used to combine multiple character sequences. 
    e <-try(read.table(file_name, skip=no_skip), silent=FALSE)   #error management
    if(class(e) == "try-error") next  #if the file doesn't exist, skip the index; it could happen because the measurement dates may not continuous (e.g. the weekend was skipped)
    else data_list[[i]] <- read.table(file_name, skip=no_skip)  #read the text file and saved into the item of the list
  }
  
  for(i in ini_no2:end_no2) {
    file_name <- paste(relative_path1,relative_path2,"/", as.character(ini_date), "Eco_", as.character(as.integer(i)), "_", plate_id, ".TXT",sep="")
    e <-try(read.table(file_name, skip=no_skip), silent=FALSE)   #error management
    if(class(e) == "try-error") next  #if the file doesn't exist, skip the index
    else data_list[[i]] <- read.table(file_name, skip=no_skip)
  } 
  data_list  #output (return value) of this function
}

#example
test_data<-load_ecoplate_data_demo01(relative_path2="2015_Feb", plate_id="T2") 
test_data[[10]]
is.data.frame(test_data[[10]])   #The element of the list is data frame
length(test_data)
test_data[[10]]$X1

#With different data-saving structure
load_ecoplate_soil_2014_Dec <- function(relative_path1="./text_file/", relative_path2="2014_Dec", ini_date01=20141216, end_date01=20141231, ini_date02=20150101, end_date02=20150114, plate_id="Eco_T3", no_skip=5)
{
  #prepare the indices for for loop  
  j1 = ini_date01
  p1 = end_date01 - ini_date01 + 1
  j2 = ini_date02   #change j value because the incubation was across years.
  p2 = end_date02 - ini_date02 + 1
  
  data_list <-list() #generate empty list
  
  for(i in 1:p1) {
    file_name <- paste(relative_path1,relative_path2,"/", as.character(j1), "/", as.character(as.integer(j1)), plate_id, ".TXT",sep="")
    data_list[[i]] <- read.table(file_name, skip=no_skip)
    j1=j1+1
  }
  
  for(i in (p1+1):(p1+p2)) {
    file_name <- paste(relative_path1,relative_path2,"/", as.character(j2), "/", as.character(as.integer(j2)), plate_id, ".TXT",sep="")
    data_list[[i]] <- read.table(file_name, skip=no_skip)
    j2=j2+1
  } 
  data_list  #output
}


#for loading microcosm data
load_ecoplate_microcosm_2014 <- function(relative_path="./text_file/microcosm_data/", file_name="Exp1day0.txt", OTU_name=otu_loss_name, no_samples=24, list="data")
{
  
  data_list<-list()  #generate empty list
  name_list<-data.frame()
  path_name <- paste(relative_path, file_name,sep="")
  
  #loading ecoplate sample names
  if(list=="sample") {
    for(i in 1: no_samples) {
      n1<-scan(path_name, what= "raw", sep="\t", skip = (1+12*(i-1)), nlines=1)
      name_list[i,1]<-n1[1]
    }
  
    #converting data in text file to sample ID (1 = control, others = composition with 19 spp)
    for(j in 1: 21) {
      for(i in 1: no_samples) {
        if(grepl(OTU_name[j], name_list[i,1])==TRUE) {
          name_list[i,2] = j;
        }
      }
    }
    colnames(name_list)<-c("sample_date", "composition")
  }
  
  
  #loading ecoplate color depth data
  if(list=="data") {
    for(i in 1: no_samples) {
      s1<-scan(path_name, what= "raw", sep="\t", skip = (3+12*(i-1)), nlines=8)
      s1d<-data.frame(matrix(s1, ncol=13, byrow=T))
      rownames(s1d)<-s1d[,1]
      s2<-s1d[,-1]
      colnames(s2) <-c("X1", "X2","X3","X4", "X5", "X6", "X7", "X8","X9", "X10","X11", "X12")
      for(j in 1:12) s2[,j] <- as.numeric(levels(s2[,j])[s2[,j]])  #Converting Factor to Numeric
      data_list[[i]]<-s2
    
    }
  }
  if(list=="data") {
    data_list
  }
  else {name_list}
}

#[2] Definitions of functions to data arrangement (integration, maximum, average, minimum, etc)

#The function for taking integration (=calculation of the area below the color development curve) even when there are gaps in measurement dates
###############list of parameters#######################
#data_eco: data source of ecopalte color patterns (list format)
#period: lenth of integration
#Output is a data frame
integ_ecoplate <- function(data_eco, period=length(data_eco))
{
  data_integ <- data.frame()  #prepare a data frame
  data_integ <- data_eco[[1]] #Note that "data_eco" should be a list; each item is the information of color development from each measurement date.
  if(period >= length(data_eco)) period = length(data_eco)  #period could be longer than the maximum length because other samples may have longer sampling period (length)
  for(i in 2:period) {
    
    #first check the NULL dates (lack in measurement)
    if(is.null(data_eco[[i]])) {
      for(j in (i+1):length(data_eco)) { #search for all NULL dates until non-Null date
        if(!is.null(data_eco[[j]])) {  #when encountring the non-Null date (j)
          for(k in i: (j-1)) data_eco[[k]] = ((j - k)*data_eco[[i-1]] + (k - i + 1)*data_eco[[j]])/(j - i + 1)  #use linear line of two non-NULL dates (i-1 & j)
          break   #skip further loop when non-zero data is found
          }#end of if data_eco[[j]]
      }#end of for j
    }#end of if data_eco[[i]]  
    data_integ <- data_integ + data_eco[[i]]  #each element of data_eco[[i]] was individually added.
  }#end of for loop i
  data_integ/period  #output, take the average value over the integration period, i.e., the normalization by integration period
}

#example
test_data<-load_ecoplate_data_demo01(relative_path2="2015_Feb", plate_id="T2") 
length(test_data)
test_integ<-integ_ecoplate(test_data)
test_integ
is.data.frame(test_integ)   #The element of the list is data frame


#The function for averaging triplicate
#The parameter is data frame
#Output is a vector
ave_ecoplate <- function(data_f){
  data_ave1<-(data_f$X1+data_f$X5+data_f$X9)/3.0  #take average
  data_ave2<-(data_f$X2+data_f$X6+data_f$X10)/3.0
  data_ave3<-(data_f$X3+data_f$X7+data_f$X11)/3.0
  data_ave4<-(data_f$X4+data_f$X8+data_f$X12)/3.0
  data_sum<-append(append(append(data_ave1,data_ave2),data_ave3),data_ave4)  #combine data
  data_sum_nor<-data_sum - data_sum[1] #normalizing by water well
  data_sum_nor #output
}

#The function for taking maximum of triplicate
#The parameter is data frame
#Output is a vector
max_ecoplate <- function(data_f){
  data_max<-max(data_f$X1[1],data_f$X5[1], data_f$X9[1])
  for(i in 2:8) data_max <-append(data_max, max(data_f$X1[i], data_f$X5[i], data_f$X9[i]))
  for(i in 1:8) data_max <-append(data_max, max(data_f$X2[i], data_f$X6[i], data_f$X10[i]))
  for(i in 1:8) data_max <-append(data_max, max(data_f$X3[i], data_f$X7[i], data_f$X11[i]))
  for(i in 1:8) data_max <-append(data_max, max(data_f$X4[i], data_f$X8[i], data_f$X12[i]))
  data_max_nor <- data_max - data_max[1]  #normalizing by water well
  data_max_nor  #output
}

#The function for takign minumu of triplicate
#The parameter is data frame
#Output is a vector
min_ecoplate <- function(data_f){
  data_min<-min(data_f$X1[1], data_f$X5[1], data_f$X9[1])
  for(i in 2:8) data_min <-append(data_min, min(data_f$X1[i], data_f$X5[i], data_f$X9[i]))
  for(i in 1:8) data_min <-append(data_min, min(data_f$X2[i], data_f$X6[i], data_f$X10[i]))
  for(i in 1:8) data_min <-append(data_min, min(data_f$X3[i], data_f$X7[i], data_f$X11[i]))
  for(i in 1:8) data_min <-append(data_min, min(data_f$X4[i], data_f$X8[i], data_f$X12[i]))
  data_min_nor <- data_min - data_min[1]  #normalizing by water well
  data_min_nor  #output
}

#For non-integrated dataset
#Functional definition using max, min, or average
#data_f: the list of ecoplate, noting that data_f[[i]] is still a list
#final_max=FALSE: this is the case when the final date measurement is used.
#final_max=TRUE: this is the case when the maximum value along the incubation period is used.
#Output is a vector
stat_non_integ_ecoplate <- function(data_f, method="max", final_max=FALSE)
{
  
  if(method == "max") {
    data_summary <- max_ecoplate(data_f[[1]])
    for(i in 2:length(data_f)) {
      if(is.null(data_f[[i]])) next;  #error management, skipping the non-measured dates
      data_summary<-rbind(data_summary, max_ecoplate(data_f[[i]]))
    }#end of for i
  }
  if(method == "min") {
    data_summary <- min_ecoplate(data_f[[1]])
    for(i in 2:length(data_f)) {
      if(is.null(data_f[[i]])) next;  #error management, skipping the non-measured dates
      data_summary<-rbind(data_summary, min_ecoplate(data_f[[i]]))
    }#end of for i
  }
  if(method == "average") {
    data_summary <- ave_ecoplate(data_f[[1]])
    for(i in 2:length(data_f)) {
      if(is.null(data_f[[i]])) next;  #error management, skipping the non-measured dates
      data_summary<-rbind(data_summary, ave_ecoplate(data_f[[i]]))
    }#end of for i
  }
  
  if(final_max==FALSE) data_summary[nrow(data_summary),] #output the final date data
  else{
    data_summary2 <-max(data_summary[,1])
    for(i in 2:32) {
      data_summary2 <- append(data_summary2, max(data_summary[,i]))
    }
    data_summary2 #output the maximum value for each substrate from independent date
  }#end of else
}

#example
stat_non_integ_ecoplate(test_data, method="max", final_max=FALSE)
stat_non_integ_ecoplate(test_data, method="max", final_max=TRUE)

###################Loading and arranging the data######################################
#[0] loading format information
#Date, treatment information from Xitou dataset
format_d<-read.csv("format_xitou_pattern.csv", header=T)

#Index information for data of microcosm experiments
otu_loss_name=c("control", "s1", "s5", "s201", "s203", "s204", "s207", "s208", "s210", "s211", "s212", "s213", "s216", "s217", "s218", "s219", "s223", "s225", "s254", "s258", "s259")

#metadata from microcosm experiments
format_microcosm<-read.csv("./COG_data.csv", header=T)

#[1]loading the data

################Loading and arranging all of the data from Xitou Ecoplate
#loading the data
{
  data_20142015_Xitou_Eco<-list()  #Prepare an empty list
  
  data_20142015_Xitou_Eco[[1]] <-load_ecoplate_soil_2014_Dec(plate_id="Eco_N1", end_date02=20150114)
  data_20142015_Xitou_Eco[[2]]<-load_ecoplate_soil_2014_Dec(plate_id="Eco_N2", end_date02=20150114)
  data_20142015_Xitou_Eco[[3]]<-load_ecoplate_soil_2014_Dec(plate_id="Eco_N3", end_date02=20150114)
  data_20142015_Xitou_Eco[[4]]<-load_ecoplate_soil_2014_Dec(plate_id="Eco_T1", end_date02=20150114)
  data_20142015_Xitou_Eco[[5]]<-load_ecoplate_soil_2014_Dec(plate_id="Eco_T2", end_date02=20150114)
  data_20142015_Xitou_Eco[[6]]<-load_ecoplate_soil_2014_Dec(plate_id="Eco_T3", end_date02=20150114)
  data_20142015_Xitou_Eco[[7]]<-load_ecoplate_data_demo01(relative_path2="2015_Feb", ini_date=20150228, plate_id="N1")
  data_20142015_Xitou_Eco[[8]]<-load_ecoplate_data_demo01(relative_path2="2015_Feb", ini_date=20150228, plate_id="N2")
  data_20142015_Xitou_Eco[[9]]<-load_ecoplate_data_demo01(relative_path2="2015_Feb", ini_date=20150228, plate_id="N3")
  data_20142015_Xitou_Eco[[10]]<-load_ecoplate_data_demo01(relative_path2="2015_Feb", ini_date=20150228, plate_id="T1")
  data_20142015_Xitou_Eco[[11]]<-load_ecoplate_data_demo01(relative_path2="2015_Feb", ini_date=20150228, plate_id="T2")
  data_20142015_Xitou_Eco[[12]]<-load_ecoplate_data_demo01(relative_path2="2015_Feb", ini_date=20150228, plate_id="T3")
  data_20142015_Xitou_Eco[[13]]<-load_ecoplate_data_demo01(relative_path2="2015_May", ini_date=20150507, end_no2=19, plate_id="N1")
  data_20142015_Xitou_Eco[[14]]<-load_ecoplate_data_demo01(relative_path2="2015_May", ini_date=20150507, end_no2=19, plate_id="N2")
  data_20142015_Xitou_Eco[[15]]<-load_ecoplate_data_demo01(relative_path2="2015_May", ini_date=20150507, end_no2=19, plate_id="N3")
  data_20142015_Xitou_Eco[[16]]<-load_ecoplate_data_demo01(relative_path2="2015_May", ini_date=20150507, end_no2=19, plate_id="T1")
  data_20142015_Xitou_Eco[[17]]<-load_ecoplate_data_demo01(relative_path2="2015_May", ini_date=20150507, end_no2=19, plate_id="T2")
  data_20142015_Xitou_Eco[[18]]<-load_ecoplate_data_demo01(relative_path2="2015_May", ini_date=20150507, end_no2=19, plate_id="T3")
  data_20142015_Xitou_Eco[[19]]<-load_ecoplate_data_demo01(relative_path2="2015_July", ini_date=20150703, plate_id="T1")
  data_20142015_Xitou_Eco[[20]]<-load_ecoplate_data_demo01(relative_path2="2015_July", ini_date=20150703, plate_id="T2")
  data_20142015_Xitou_Eco[[21]]<-load_ecoplate_data_demo01(relative_path2="2015_July", ini_date=20150703, plate_id="T3")
  data_20142015_Xitou_Eco[[22]]<-load_ecoplate_data_demo01(relative_path2="2015_July", ini_date=20150703, plate_id="N1")
  data_20142015_Xitou_Eco[[23]]<-load_ecoplate_data_demo01(relative_path2="2015_July", ini_date=20150703, plate_id="N2")
  data_20142015_Xitou_Eco[[24]]<-load_ecoplate_data_demo01(relative_path2="2015_July", ini_date=20150703, plate_id="N3")
}
#example
data_20142015_Xitou_Eco[[24]]
length(data_20142015_Xitou_Eco[[24]])

#Calculating integration (& averaging) of data
{
integ_20142015_Xitou_Eco<-list()
for(i in 1:24) integ_20142015_Xitou_Eco[[i]] <-integ_ecoplate(data_20142015_Xitou_Eco[[i]])
}
#example
integ_20142015_Xitou_Eco[[24]]
is.data.frame(integ_20142015_Xitou_Eco[[24]])

#Treatments on triplicates (average, maximum, or minimum)

#for Averaging the intergrated values
{
ave_integ_20142015_Xitou_Eco<-list()
for(i in 1:24) ave_integ_20142015_Xitou_Eco[[i]] <-ave_ecoplate(integ_20142015_Xitou_Eco[[i]])
}

#for Taking max of the intergrated values
{
max_integ_20142015_Xitou_Eco<-list()
for(i in 1:24) max_integ_20142015_Xitou_Eco[[i]] <-max_ecoplate(integ_20142015_Xitou_Eco[[i]])
}

#for Taking min of the intergrated values
{
min_integ_20142015_Xitou_Eco<-list()
for(i in 1:24) min_integ_20142015_Xitou_Eco[[i]] <-min_ecoplate(integ_20142015_Xitou_Eco[[i]])
}

#calculating the summary without integration
{
max_non_integ_final_20142015_Xitou_Eco<-list()   #for the case using the final day data only
min_non_integ_final_20142015_Xitou_Eco<-list()
ave_non_integ_final_20142015_Xitou_Eco<-list()
max_non_integ_max_20142015_Xitou_Eco<-list()  #for the case using the data at maximum
min_non_integ_max_20142015_Xitou_Eco<-list()
ave_non_integ_max_20142015_Xitou_Eco<-list()

for(i in 1:24) {  
  max_non_integ_final_20142015_Xitou_Eco[[i]] <-stat_non_integ_ecoplate(data_20142015_Xitou_Eco[[i]], method="max", final_max=FALSE)
  min_non_integ_final_20142015_Xitou_Eco[[i]] <-stat_non_integ_ecoplate(data_20142015_Xitou_Eco[[i]], method="min", final_max=FALSE)
  ave_non_integ_final_20142015_Xitou_Eco[[i]] <-stat_non_integ_ecoplate(data_20142015_Xitou_Eco[[i]], method="average", final_max=FALSE)
  max_non_integ_max_20142015_Xitou_Eco[[i]] <-stat_non_integ_ecoplate(data_20142015_Xitou_Eco[[i]], method="max", final_max=TRUE)
  min_non_integ_max_20142015_Xitou_Eco[[i]] <-stat_non_integ_ecoplate(data_20142015_Xitou_Eco[[i]], method="min", final_max=TRUE)
  ave_non_integ_max_20142015_Xitou_Eco[[i]] <-stat_non_integ_ecoplate(data_20142015_Xitou_Eco[[i]], method="average", final_max=TRUE)
}

}

#############################End of loading and arranging data from Xitou############################

###############################For microcosm data#######################################
#Data loading from microcosm experiments
{
  #color data
  data_microcosm01<-load_ecoplate_microcosm_2014(file_name="Exp1day0.txt", OTU_name=otu_loss_name, no_samples=24, list="data")
  #sample_ID
  name_microcosm01<-load_ecoplate_microcosm_2014(file_name="Exp1day0.txt", OTU_name=otu_loss_name, no_samples=24, list="sample")
  data_microcosm02<-load_ecoplate_microcosm_2014(file_name="Exp2day0.txt", OTU_name=otu_loss_name, no_samples=27, list="data")
  name_microcosm02<-load_ecoplate_microcosm_2014(file_name="Exp2day0.txt", OTU_name=otu_loss_name, no_samples=27, list="sample")
  data_microcosm03<-load_ecoplate_microcosm_2014(file_name="Exp3day0.txt", OTU_name=otu_loss_name, no_samples=28, list="data")
  name_microcosm03<-load_ecoplate_microcosm_2014(file_name="Exp3day0.txt", OTU_name=otu_loss_name, no_samples=28, list="sample")
  data_microcosm04<-load_ecoplate_microcosm_2014(file_name="Exp4day0.txt", OTU_name=otu_loss_name, no_samples=27, list="data")
  name_microcosm04<-load_ecoplate_microcosm_2014(file_name="Exp4day0.txt", OTU_name=otu_loss_name, no_samples=27, list="sample")
}

#calculating the summary without integration for microcosm experiments
{
max_non_integ_max_01_microcosm<-list()
min_non_integ_max_01_microcosm<-list()
ave_non_integ_max_01_microcosm<-list()
max_non_integ_max_02_microcosm<-list()
min_non_integ_max_02_microcosm<-list()
ave_non_integ_max_02_microcosm<-list()
max_non_integ_max_03_microcosm<-list()
min_non_integ_max_03_microcosm<-list()
ave_non_integ_max_03_microcosm<-list()
max_non_integ_max_04_microcosm<-list()
min_non_integ_max_04_microcosm<-list()
ave_non_integ_max_04_microcosm<-list()

for(i in 1:length(data_microcosm01)) {
  max_non_integ_max_01_microcosm[[i]] <-max_ecoplate(data_microcosm01[[i]])
  min_non_integ_max_01_microcosm[[i]] <-min_ecoplate(data_microcosm01[[i]])
  ave_non_integ_max_01_microcosm[[i]] <-ave_ecoplate(data_microcosm01[[i]])
}
for(i in 1:length(data_microcosm02)) {
  max_non_integ_max_02_microcosm[[i]] <-max_ecoplate(data_microcosm02[[i]])
  min_non_integ_max_02_microcosm[[i]] <-min_ecoplate(data_microcosm02[[i]])
  ave_non_integ_max_02_microcosm[[i]] <-ave_ecoplate(data_microcosm02[[i]])
}
for(i in 1:length(data_microcosm03)) {
  max_non_integ_max_03_microcosm[[i]] <-max_ecoplate(data_microcosm03[[i]])
  min_non_integ_max_03_microcosm[[i]] <-min_ecoplate(data_microcosm03[[i]])
  ave_non_integ_max_03_microcosm[[i]] <-ave_ecoplate(data_microcosm03[[i]])
}
for(i in 1:length(data_microcosm04)) {
  max_non_integ_max_04_microcosm[[i]] <-max_ecoplate(data_microcosm04[[i]])
  min_non_integ_max_04_microcosm[[i]] <-min_ecoplate(data_microcosm04[[i]])
  ave_non_integ_max_04_microcosm[[i]] <-ave_ecoplate(data_microcosm04[[i]])
}


#For the whole data######################
max_non_integ_max_total_microcosm<-append(max_non_integ_max_01_microcosm, max_non_integ_max_02_microcosm)
max_non_integ_max_total_microcosm<-append(max_non_integ_max_total_microcosm, max_non_integ_max_03_microcosm)
max_non_integ_max_total_microcosm<-append(max_non_integ_max_total_microcosm, max_non_integ_max_04_microcosm)

}
###############End of loading and arranging microcosms data############################################


##############PART3 Ecoplate analysis combined with chemical similarity###############################
library("picante")
library("devtools")  #for installing packages that are not uploaded to the common library depository
library("GUniFrac")
#library("ggplot2")
library("labdsv")
library("SYNCSA")
library("pforeach") #for palarell computing


#Script for loading text file of ecoplate for further calculations

#Sample codes for parallel computing
#i<-0
#pforeach(j=1:10)(
#  {i<<-i+1}
#)
#i
#ss
#x<-0
#y<-0
#cat(1,2)
#x<-pforeach(i=1:3000, .c=list)({
#  y=i*2
#  y2=i*3
#  z=cbind(y, y2)
#  z
#})
#z
############################################## i)Definitions of functions##############################################
{
##########Function, which calculate percentile(quantile)-based multifunctionality (MF) and with or without weighting by chemical_similarity, then conduct simple linear analysis to check if the multifunctionality can be explained by month and treatment####################
###############list of parameters#######################
#data: data source of ecopalte color patterns (list format)
#thres: the threshold for the quantile-based multifunctionality, 0.9 means the top 10% values are converted to 1 (presence of function) while others become 0 (absence of the function)
#format_data: metadata (sampling dates, treatments, etc)
#chem_n: list of chemical names
#chem_max: the similarity matrix for ecoplate substrates
#perm: number of permulations for statistical test
#n_cores: number of CPU cores used for parallel permutations
Chem_Disim_MF_xitou <- function(data=ave_integ_20142015_Xitou_Eco, thres=0.9, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, perm=999, n_cores=2) 
{
  #Loading all data (24 samples), some modification is necessary for the data size but in our case there are 24 samples in total.
  file_summary <- data[[1]]
  for(i in 2:length(data)) file_summary <- rbind(file_summary, data[[i]])  #data are combined in to a single file
  file_summary <-file_summary[,c(-1)] #exclusing water data (control, zero)
  row.names(file_summary)<-format_data$sample  #Add name (sample ID)
  colnames(file_summary)<-chem_n #Add name (substrate name)
  
  #Convert data into binary values by quantile (quantile-based multifunctionality)
  mf<-list()  #prepare the output list
  for(i in 1:31) mf[[i]] <-(file_summary[,i] > quantile(file_summary[,i], thres))  # return true if the color value qunatile is greater than threshold
  M_mf<-as.data.frame(mf)  #converted into dataframe
  row.names(M_mf)<-format_data$sample
  colnames(M_mf)<-chem_n
  M_mf[M_mf==TRUE]<-1  #converting T,F to 1 and 0; ts is now dataframe to represent the binary matrix (E_BC in Figure 3)
  
  #Calculating tree based on chemical similarity matrix
  clus.hier<-hclust(chem_matrix, method="average")
  tree<-as.phylo(clus.hier)  #convert to phytogenetic tree format; this process is necessary to use phylogenetic diversity (PD) function from picante for chemical similarity tree.
  
  #Using function in picante (phylogenetic community analysis)
  pd.result <- pd(M_mf, tree, include.root=T)  #Claculate chemical-similarity-weighted multifunctionality, and non-weighted (raw) multifunctionality
  colnames(pd.result)<-c("Chem_MF", "MF")   #add name
  pd.result<-cbind(format_data, pd.result)  #incorporate sample data (day & treatment)
  
  #Plot chemical-similarity-weighted multifunctionality: Normal (control) vs Treatment (trenching)
  plot(subset(pd.result, treatment=="N")$day, subset(pd.result, treatment=="N")$Chem_MF, ylim=c(0,6))
  par(new=T)
  plot(subset(pd.result, treatment=="T")$day, subset(pd.result, treatment=="T")$Chem_MF, col=3,ylim=c(0,6))
  
  
  ############Conduct a simple linear regression model#####################
  #MF is chemical simmilarity independent multifuncitonality (0-31)
  model_nw<-lm(MF~day*treatment, pd.result)
  #PD is chemical-similarity-weighted multifuncitonality 
  model_w<-lm(Chem_MF~day*treatment, pd.result)
  
  #Output1
  cat("Non-weighted MF can be explained by day, treatment, and their interaction?\n")
  print(summary(model_nw))
  cat("Weighted MF can be explained by day, treatment, and their interaction?\n")
  print(summary(model_w))
  
  
  ####Importantly, it is not sure if the explanation power (R2) of the linear regression on chemical-similarity weighted MF is statistically meaningfull or not due to the uncertainty of the chemical-similarity definition. Therefore, it is necessary to conduct permutation test through comparing R2 value from "model_w" and R2 values coming from the analysis using "randomly-generated chemical similarity tree".
  
  #To calculate R2 values from ANOVA results
  anv_w<-anova(model_w)
  R2_w=(anv_w$`Sum Sq`[1]+anv_w$`Sum Sq`[2]+anv_w$`Sum Sq`[3])/(anv_w$`Sum Sq`[1]+anv_w$`Sum Sq`[2]+anv_w$`Sum Sq`[3]+anv_w$`Sum Sq`[4])
  
  
  #######Checking the significance of the chemical information (permutation test)
  chem_dist.mat<-as.matrix(chem_matrix)   #format arrangement on data
  chem_dist.dist<-as.dist(chem_matrix)    #format arrangement on data
  size<-nrow(chem_dist.mat)
  leng<-length(chem_dist.dist)
  
  p_iU = 0  #Psudo-P value for chemical weight
  p_iL = 0  #Psudo-P value for chemical weight
  
  #############Use parallel calculations ".cores" specifies the number of CPU cores to be used.
  R2_temp<-pforeach(i = 1: perm, .c=list, .cores=n_cores) ({
    
  
    ###generating random tree, directly shufflting distance (similarity) matrix###############
    z<-sample(1:leng, leng, replace=FALSE) #re-sampling without replacement
    temp_d<-chem_dist.dist
    for(i in 1:leng) temp_d[i]<-chem_dist.dist[z[i]]
    #############################################################################
    
    temp_chem<-hclust(temp_d, method="average")    #clustering on randomly-generated similarity matrix
    temp.ph<-as.phylo(temp_chem)  #randomly-generated tree
    plot(temp.ph)
    
    ####Repeat the identical analysis using randomly-genreated tree
    pd.result_temp <- pd(M_mf, temp.ph, include.root=T)
    pd.result_temp<-cbind(format_data, pd.result_temp)
    model_w_temp<-lm(PD~day*treatment, pd.result_temp)
  
    #Calculate R2 value
    anv_w_temp<-anova(model_w_temp)
    R2_w_temp=(anv_w_temp$`Sum Sq`[1]+anv_w_temp$`Sum Sq`[2]+anv_w_temp$`Sum Sq`[3])/(anv_w_temp$`Sum Sq`[1]+anv_w_temp$`Sum Sq`[2]+anv_w_temp$`Sum Sq`[3]+anv_w_temp$`Sum Sq`[4])
    
    #if(R2_w_temp > R2_w) p_i = p_i + 1  #Counting the case when randomly-generated tree give a greater R2 value than actually estimated tree
  })
  
  #Count the cases when R2 value from random-tree is greater or smaller than R2_w in order to calculate Psuedo-P value of permutation
  for(i in 1: perm) if(R2_temp[[i]]>= R2_w) p_iU = p_iU + 1
  for(i in 1: perm) if(R2_temp[[i]]<= R2_w) p_iL = p_iL + 1
  
  #Output2 of permutation results
  cat("The probability with which the # of R2 values from permutated tree >= observed R2 value is more extreme under null hypothesis\n (P-Value for the chemical tree)\n")
  print((p_iU+1)/(perm+1))
  cat("The probability with which the # of R2 values from permutated tree <= observed R2 value is more extreme under null hypothesis\n (P-Value for the chemical tree)\n")
  print((p_iL+1)/(perm+1))
  
}

#Function, which returns cummulative R2 values of day, teatment, and their interaction by PERMANOVA, based on distance matrix
#Note that the results from this function and the related functions below of PERMANOVA were not shown in the manuscript.
###############list of parameters#######################
#data: data source of ecopalte color patterns (list format)
#format_data: metadata (sampling dates, treatments, etc)
#chem_n: list of chemical names
#chem_max: the similarity matrix for ecoplate substrates
#meth: the method to calculate the dissimilarity of color patterns
#perm: number of permulations for PERMANOVA
#perm2: number of premutations for chemical similarity tree
#n_cores: number of CPU cores used for parallel permutations
w_s_adonis_xitou <- function(data=max_integ_20142015_Xitou_Eco, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, meth="bray", perm=9999, perm2=99, n_cores=2) {

  ####IMPORTANT: The comments within this function are simplified.Check every comment on the function above "Chem_Disim_MF_xitou"##########   
  #First compile all of the 24 day's data
  file_summary <- data[[1]]
  for(i in 2:length(data)) file_summary <- rbind(file_summary, data[[i]])
  file_summary <-file_summary[,c(-1)] #excluding water data (control, zero)
  file_summary[file_summary <0]<-0 #convert negative values to zero, which is necessary to calculate the bray-curtis distance
  file_summary <- as.data.frame(file_summary) #converting to dataframe
  row.names(file_summary)<-format_data$sample #set names
  colnames(file_summary)<-chem_n              #set names   
  
  #For the analysis after conversion of color depth to presence/absence, by quantile-based binarization (use three thresholds)
  mf0.9<-list()  #for threshold 0.9
  for(i in 1:31) mf0.9[[i]] <-(file_summary[,i] > quantile(file_summary[,i], 0.9))
  mf0.9_summary<-as.data.frame(mf0.9)
  row.names(mf0.9_summary)<-format_data$sample
  colnames(mf0.9_summary)<-chem_n
  mf0.9_summary[mf0.9_summary==TRUE]<-1  #converting T,F to 1 and 0
  
  mf0.7<-list() #for threshold 0.7
  for(i in 1:31) mf0.7[[i]] <-(file_summary[,i] > quantile(file_summary[,i], 0.7))
  mf0.7_summary<-as.data.frame(mf0.7)
  row.names(mf0.7_summary)<-format_data$sample
  colnames(mf0.7_summary)<-chem_n
  mf0.7_summary[mf0.7_summary==TRUE]<-1  #converting T,F to 1 and 0mf0.7<-list()
  
  mf0.5<-list() #for threshold 0.5
  for(i in 1:31) mf0.5[[i]] <-(file_summary[,i] > quantile(file_summary[,i], 0.5))
  mf0.5_summary<-as.data.frame(mf0.5)
  row.names(mf0.5_summary)<-format_data$sample
  colnames(mf0.5_summary)<-chem_n
  mf0.5_summary[mf0.5_summary==TRUE]<-1  #converting T,F to 1 and 0mf0.5<-list()

  #Calculating default dissimilarity based on presence/absence (without chemical-similarity information) [baseline of the analysis]
  dist_binary_uw0.9<-vegdist(mf0.9_summary, method="euclidian")
  dist_binary_uw0.7<-vegdist(mf0.7_summary, method="euclidian")
  dist_binary_uw0.5<-vegdist(mf0.5_summary, method="euclidian")
  
  #######Calculating the default continuous dissimilarity (without chemical-similarity information) [baseline of the analysis]
  dist_cont_uw<-vegdist(file_summary, method=meth)  
  
  #Calculating tree based on chemical similarity matrix
  clus.hier<-hclust(chem_matrix, method="average")
  tree<-as.phylo(clus.hier)  #convert to phytogenetic tree
  
  
  #Calculating unifrac distance (chemical-similarity-weighted color pattern dissimilarity) on the binarized data by using {picante:unifrac}
  dist_binary_cw0.9<-unifrac(mf0.9_summary, tree)
  dist_binary_cw0.7<-unifrac(mf0.7_summary, tree)
  dist_binary_cw0.5<-unifrac(mf0.5_summary, tree)
  #error control
  dist_binary_cw0.9[is.nan(dist_binary_cw0.9)]<-0
  dist_binary_cw0.7[is.nan(dist_binary_cw0.7)]<-0
  dist_binary_cw0.5[is.nan(dist_binary_cw0.5)]<-0
  
  #Calculating generalized unifrac distance (chemical-similarity-weighted color patten dissimilarity) from continuous value using GUniFrac
  unif_w<-GUniFrac(file_summary, tree) #not weighted by color depth, but noting the difference from the value from the binarized data (we don't use this)
  dist_cont_GuniFrac<-as.dist(unif_w$unifracs[,,"d_1"])  #Weighted distance by exponent 1.0
  dist_cont_GuniFrac0.5<-as.dist(unif_w$unifracs[,,"d_0.5"])  #Weighted distance by exponent 0.5
  
  ######Calculating chemical-similarity-weighted dissimilarity using fuzzy function
  #1) Prepare the belonging (see the supporting information)
  weighting_tree<-belonging(chem_matrix, standardize = TRUE)
  
  #2) Since the order of substrates is different between file_summary and weighting_tree, this can be fitted by sorting
  file_summary_sorted<-file_summary[,sort(colnames(file_summary))]
  weighting_tree_sorted<-weighting_tree[sort(rownames(weighting_tree)),sort(colnames(weighting_tree))]
  
  #3) Fuzzy weighted color pattern matrix is calcualted
  file_summary_cwF <- as.matrix(file_summary_sorted)%*%as.matrix(weighting_tree_sorted)

  #4) Fuzzy-weighted dissimilarity is calculated on the matrix generated at the process 3) 
  dist_cont_cwF<-vegdist(file_summary_cwF, method=meth)


  
  ###########PERMANOVA for each of the different indices of color patten dissimilarity 
  #binarized, without chemical-similarity 
  diff_binary_uw0.9<-adonis(dist_binary_uw0.9 ~ format_data$day*format_data$treatment, permutations=perm)
  diff_binary_uw0.7<-adonis(dist_binary_uw0.7 ~ format_data$day*format_data$treatment, permutations=perm)
  diff_binary_uw0.5<-adonis(dist_binary_uw0.5 ~ format_data$day*format_data$treatment, permutations=perm)
  
  #binazied, with chemical-similarity
  diff_binary_cw0.9<-adonis(dist_binary_cw0.9 ~ format_data$day*format_data$treatment, permutations=perm)
  diff_binary_cw0.7<-adonis(dist_binary_cw0.7 ~ format_data$day*format_data$treatment, permutations=perm)
  diff_binary_cw0.5<-adonis(dist_binary_cw0.5 ~ format_data$day*format_data$treatment, permutations=perm)
  
  #continuous, without chemical-similarity
  diff_cont_uw<-adonis(dist_cont_uw ~ format_data$day*format_data$treatment, permutations=perm)
  #continuous, with chemical-similarity through generalized uniFrac distance
  diff_cont_GuniFrac<-adonis(dist_cont_GuniFrac ~ format_data$day*format_data$treatment, permutations=perm)
  diff_cont_GuniFrac0.5<-adonis(dist_cont_GuniFrac0.5 ~ format_data$day*format_data$treatment, permutations=perm)
  #continuous, with chemical-similarity through Fuzzy weighting
  diff_cont_cwF<-adonis(dist_cont_cwF ~ format_data$day*format_data$treatment, permutations=perm)
  
  #Output1, P values 
  cat("The P values of data, treatment and their interaction (binarized w/o CW 0.9):", diff_binary_uw0.9$aov.tab$Pr, "\n")
  cat("The P values of data, treatment and their interaction (binarized w/ CW 0.9):", diff_binary_cw0.9$aov.tab$Pr, "\n")
  cat("The P values of data, treatment and their interaction (binarized w/o CW 0.7):", diff_binary_uw0.7$aov.tab$Pr, "\n")
  cat("The P values of data, treatment and their interaction (binarized w/ CW 0.7):", diff_binary_cw0.7$aov.tab$Pr, "\n")
  cat("The P values of data, treatment and their interaction (binarized w/o CW 0.5):", diff_binary_uw0.5$aov.tab$Pr, "\n")
  cat("The P values of data, treatment and their interaction (binarized w/ CW 0.5):", diff_binary_cw0.5$aov.tab$Pr, "\n")
  cat("The P values of data, treatment and their interaction (continuous w/o CW):", diff_cont_uw$aov.tab$Pr, "\n")
  cat("The P values of data, treatment and their interaction (continuous w/ CW GuniFrac):", diff_cont_GuniFrac$aov.tab$Pr, "\n")
  cat("The P values of data, treatment and their interaction (continuous w/ CW GuniFrac 0.5):", diff_cont_GuniFrac0.5$aov.tab$Pr, "\n")
  cat("The P values of data, treatment and their interaction (continuous w/ CW fuzzy):", diff_cont_cwF$aov.tab$Pr, "\n")
  
  #Output2, R2 values
  cat("The R2 of data, treatment and their interaction (binarized w/o CW 0.9):", 1-diff_binary_uw0.9$aov.tab$R2[4], "\n")
  cat("The R2 of data, treatment and their interaction (binarized w/ CW 0.9):", 1-diff_binary_cw0.9$aov.tab$R2[4], "\n")
  cat("The R2 of data, treatment and their interaction (binarized w/o CW 0.7):", 1-diff_binary_uw0.7$aov.tab$R2[4], "\n")
  cat("The R2 of data, treatment and their interaction (binarized w/ CW 0.7):", 1-diff_binary_cw0.7$aov.tab$R2[4], "\n")
  cat("The R2 of data, treatment and their interaction (binarized w/o CW 0.5):", 1-diff_binary_uw0.5$aov.tab$R2[4], "\n")
  cat("The R2 of data, treatment and their interaction (binarized w/ CW 0.5):", 1-diff_binary_cw0.5$aov.tab$R2[4], "\n")
  cat("The R2 of data, treatment and their interaction (continuous w/o CW):", 1-diff_cont_uw$aov.tab$R2[4], "\n")
  cat("The R2 of data, treatment and their interaction (continuous w/ CW GuniFrac):", 1-diff_cont_GuniFrac$aov.tab$R2[4], "\n")
  cat("The R2 of data, treatment and their interaction (continuous w/ CW GuniFrac 0.5):", 1-diff_cont_GuniFrac0.5$aov.tab$R2[4], "\n")
  cat("The R2 of data, treatment and their interaction (continuous w/ CW fuzzy):", 1-diff_cont_cwF$aov.tab$R2[4], "\n")
  
  ####Importantly, it is not sure if the explanation power (R2) of PERMANOVA on chemical-similarity-weighted color pattens is statistically meaningfull or not due to the uncertainty of the chemical-similarity definition. Therefore, it is necessary to conduct permutation test through comparing R2 value from the above models and R2 values coming from the analysis using "randomly-generated chemical similarity tree".  
  
  #Preparation
  chem_dist.mat<-as.matrix(chem_matrix)
  chem_dist.dist<-as.dist(chem_matrix)
  size<-nrow(chem_dist.mat)
  leng<-length(chem_dist.dist)
  
  #Summarize the R2 values from the above models
  R2_binary_cw0.9 = 1-diff_binary_cw0.9$aov.tab$R2[4]
  R2_binary_cw0.7 = 1-diff_binary_cw0.7$aov.tab$R2[4]
  R2_binary_cw0.5 = 1-diff_binary_cw0.5$aov.tab$R2[4]
  R2_cont_cwF = 1-diff_cont_cwF$aov.tab$R2[4]
  R2_cont_GuniFrac = 1-diff_cont_GuniFrac$aov.tab$R2[4]
  R2_cont_GuniFrac0.5 = 1-diff_cont_GuniFrac0.5$aov.tab$R2[4]
  
  #For upper P values
  Up_i_binary_cw0.9 = 0  #Psudo-P value for chemical weight
  Up_i_binary_cw0.7 = 0  #Psudo-P value for chemical weight
  Up_i_binary_cw0.5 = 0  #Psudo-P value for chemical weight
  Up_i_cont_cwF = 0  #Psudo-P value for chemical weight
  Up_i_cont_GuniFrac = 0  #Psudo-P value for chemical weight
  Up_i_cont_GuniFrac0.5 = 0  #Psudo-P value for chemical weight
  
  #For lower P values
  Lp_i_binary_cw0.9 = 0  #Psudo-P value for chemical weight
  Lp_i_binary_cw0.7 = 0  #Psudo-P value for chemical weight
  Lp_i_binary_cw0.5 = 0  #Psudo-P value for chemical weight
  Lp_i_cont_cwF = 0  #Psudo-P value for chemical weight
  Lp_i_cont_GuniFrac = 0  #Psudo-P value for chemical weight
  Lp_i_cont_GuniFrac0.5 = 0  #Psudo-P value for chemical weight
  
  #parallel computing
  #perm2=5
  R2_temp<-pforeach(i = 1: perm2, .c=list, .cores=n_cores) ({
   
    #generating random tree, directly shufflting distance matrix##############
    z<-sample(1:leng, leng, replace=FALSE)
    temp_d<-chem_dist.dist
    for(i in 1:leng) temp_d[i]<-chem_dist.dist[z[i]]
    ###############################################################
    
    temp_chem<-hclust(temp_d, method="average")  #clustering on random-matrix
    temp.ph<-as.phylo(temp_chem) #generating tree
    plot(temp.ph)
    
    tree<-temp.ph
    chem_matrix<-temp_d
    
    #Calculating through binarized data by using {picante:unifrac}
    dist_binary_cw0.9<-unifrac(mf0.9_summary, tree)
    dist_binary_cw0.7<-unifrac(mf0.7_summary, tree)
    dist_binary_cw0.5<-unifrac(mf0.5_summary, tree)
    #error control
    dist_binary_cw0.9[is.nan(dist_binary_cw0.9)]<-0
    dist_binary_cw0.7[is.nan(dist_binary_cw0.7)]<-0
    dist_binary_cw0.5[is.nan(dist_binary_cw0.5)]<-0
    
    #Calculating from continuous value using GUniFrac
    unif_w<-GUniFrac(file_summary, tree)  #not used
    dist_cont_GuniFrac<-as.dist(unif_w$unifracs[,,"d_1"])  #Weighted distance
    dist_cont_GuniFrac0.5<-as.dist(unif_w$unifracs[,,"d_0.5"])  #Weighted distance
    
    #Calculating weighted dissimilarity using fuzzy function
    weighting_tree<-belonging(chem_matrix, standardize = TRUE)
    file_summary_sorted<-file_summary[,sort(colnames(file_summary))]
    weighting_tree_sorted<-weighting_tree[sort(rownames(weighting_tree)),sort(colnames(weighting_tree))]
    file_summary_cwF <- as.matrix(file_summary_sorted)%*%as.matrix(weighting_tree_sorted)
    dist_cont_cwF<-vegdist(file_summary_cwF, method=meth)
    
    #Calculating continuous dissimilarity without chemical-similarity info
    dist_cont_uw<-vegdist(file_summary, method=meth)
    
    #Calculating default dissimilarity based on presence/absence without chemical-similarity info
    dist_binary_uw0.9<-vegdist(mf0.9_summary, method="euclidian")
    dist_binary_uw0.7<-vegdist(mf0.7_summary, method="euclidian")
    dist_binary_uw0.5<-vegdist(mf0.5_summary, method="euclidian")
    
    #PERMANOVA, don't need permutation
    diff_binary_uw0.9<-adonis(dist_binary_uw0.9 ~ format_data$day*format_data$treatment, permutations=2)
    diff_binary_uw0.7<-adonis(dist_binary_uw0.7 ~ format_data$day*format_data$treatment, permutations=2)
    diff_binary_uw0.5<-adonis(dist_binary_uw0.5 ~ format_data$day*format_data$treatment, permutations=2)
    
    diff_binary_cw0.9<-adonis(dist_binary_cw0.9 ~ format_data$day*format_data$treatment, permutations=2)
    diff_binary_cw0.7<-adonis(dist_binary_cw0.7 ~ format_data$day*format_data$treatment, permutations=2)
    diff_binary_cw0.5<-adonis(dist_binary_cw0.5 ~ format_data$day*format_data$treatment, permutations=2)
    
    diff_cont_uw<-adonis(dist_cont_uw ~ format_data$day*format_data$treatment, permutations=2)
    diff_cont_cwF<-adonis(dist_cont_cwF ~ format_data$day*format_data$treatment, permutations=2)
    diff_cont_GuniFrac<-adonis(dist_cont_GuniFrac ~ format_data$day*format_data$treatment, permutations=2)
    diff_cont_GuniFrac0.5<-adonis(dist_cont_GuniFrac0.5 ~ format_data$day*format_data$treatment, permutations=2)
    
    #Calcualte R2 values from random tree
    R2_binary_cw0.9temp = 1-diff_binary_cw0.9$aov.tab$R2[4]
    R2_binary_cw0.7temp = 1-diff_binary_cw0.7$aov.tab$R2[4]
    R2_binary_cw0.5temp = 1-diff_binary_cw0.5$aov.tab$R2[4]
    R2_cont_cwFtemp = 1-diff_cont_cwF$aov.tab$R2[4]
    R2_cont_GuniFractemp = 1-diff_cont_GuniFrac$aov.tab$R2[4]
    R2_cont_GuniFrac0.5temp = 1-diff_cont_GuniFrac0.5$aov.tab$R2[4]
    
    cbind(R2_binary_cw0.9temp, R2_binary_cw0.7temp,R2_binary_cw0.5temp,R2_cont_cwFtemp, R2_cont_GuniFractemp, R2_cont_GuniFrac0.5temp)
  })
  
  #R2_temp
  
  #Counting the case when randomly-generated tree give a greater R2 value (or lower R2 value) than actually estimated tree
  for(i in 1:perm2) {
    if(R2_temp[[i]][1] >= R2_binary_cw0.9) Up_i_binary_cw0.9 = Up_i_binary_cw0.9 + 1 
    if(R2_temp[[i]][2] >= R2_binary_cw0.7) Up_i_binary_cw0.7 = Up_i_binary_cw0.7 + 1  
    if(R2_temp[[i]][3] >= R2_binary_cw0.5) Up_i_binary_cw0.5 = Up_i_binary_cw0.5 + 1  
    if(R2_temp[[i]][4] >= R2_cont_cwF) Up_i_cont_cwF = Up_i_cont_cwF + 1  
    if(R2_temp[[i]][5] >= R2_cont_GuniFrac) Up_i_cont_GuniFrac = Up_i_cont_GuniFrac + 1  
    if(R2_temp[[i]][6] >= R2_cont_GuniFrac0.5) Up_i_cont_GuniFrac0.5 = Up_i_cont_GuniFrac0.5 + 1  
    
    if(R2_temp[[i]][1] <= R2_binary_cw0.9) Lp_i_binary_cw0.9 = Lp_i_binary_cw0.9 + 1 
    if(R2_temp[[i]][2] <= R2_binary_cw0.7) Lp_i_binary_cw0.7 = Lp_i_binary_cw0.7 + 1  
    if(R2_temp[[i]][3] <= R2_binary_cw0.5) Lp_i_binary_cw0.5 = Lp_i_binary_cw0.5 + 1  
    if(R2_temp[[i]][4] <= R2_cont_cwF) Lp_i_cont_cwF = Lp_i_cont_cwF + 1  
    if(R2_temp[[i]][5] <= R2_cont_GuniFrac) Lp_i_cont_GuniFrac = Lp_i_cont_GuniFrac + 1  
    if(R2_temp[[i]][6] <= R2_cont_GuniFrac0.5) Lp_i_cont_GuniFrac0.5 = Lp_i_cont_GuniFrac0.5 + 1  
  }
  #Output3 
  cat("The probability with which the # of R2 values from permutated tree >= observed R2 value is more extreme under null hypothesis\n")
  cat("binary_cw0.9\t", (Up_i_binary_cw0.9+1)/(perm2+1), "\n") 
  cat("binary_cw0.7\t", (Up_i_binary_cw0.7+1)/(perm2+1), "\n") 
  cat("binary_cw0.5\t", (Up_i_binary_cw0.5+1)/(perm2+1), "\n") 
  cat("continuous GuniFrac weighted\t", (Up_i_cont_GuniFrac+1)/(perm2+1), "\n") 
  cat("continuous GuniFrac 0.5 weighted\t", (Up_i_cont_GuniFrac0.5+1)/(perm2+1), "\n")
  cat("continuous fuzzy weighted\t", (Up_i_cont_cwF+1)/(perm2+1), "\n") 
  cat("The probability with which the # of R2 values from permutated tree <= observed R2 value is more extreme under null hypothesis\n")
  cat("binary_cw0.9\t", (Lp_i_binary_cw0.9+1)/(perm2+1), "\n") 
  cat("binary_cw0.7\t", (Lp_i_binary_cw0.7+1)/(perm2+1), "\n") 
  cat("binary_cw0.5\t", (Lp_i_binary_cw0.5+1)/(perm2+1), "\n") 
  cat("continuous GuniFrac weighted\t", (Lp_i_cont_GuniFrac+1)/(perm2+1), "\n") 
  cat("continuous GuniFrac 0.5 weighted\t", (Lp_i_cont_GuniFrac0.5+1)/(perm2+1), "\n") 
  cat("continuous fuzzy weighted\t", (Lp_i_cont_cwF+1)/(perm2+1), "\n") 
}


#Function, which returns cummulative R2 values of day, teatment, and their interaction by distance-based RDA with or without weighing by chemical similarity
###############list of parameters#######################
#data: data source of ecopalte color patterns (list format)
#format_data: metadata (sampling dates, treatments, etc)
#chem_n: list of chemical names
#chem_max: the similarity matrix for ecoplate substrates
#meth: the method to calculate the dissimilarity of color patterns
#perm: number of permulations for PERMANOVA
#perm2: number of premutations for chemical similarity tree
#n_cores: number of CPU cores used for parallel permutations
W_s_dbRDA <- function(data=max_integ_20142015_Xitou_Eco, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, meth="bray", thres=0.9, perm=1999, perm2=999, n_cores=2) {
  
  ####IMPORTANT: The comments within this function are simplified.Check every comment on the function above "Chem_Disim_MF_xitou"########## 
  #First compile all of the 24 day's data
  file_summary <- data[[1]]
  for(i in 2:length(data)) file_summary <- rbind(file_summary, data[[i]])
  file_summary <-file_summary[,c(-1)] #exclusing water data (control, zero)
  file_summary[file_summary <0]<-0 #convert negative values to zero
  file_summary <- as.data.frame(file_summary) #converting to dataframe
  row.names(file_summary)<-format_data$sample
  colnames(file_summary)<-chem_n
  
  #For binarized (by color depth) analysis
  mf<-list()
  for(i in 1:31) mf[[i]] <-(file_summary[,i] > quantile(file_summary[,i], thres))
  mf_summary<-as.data.frame(mf)
  row.names(mf_summary)<-format_data$sample
  colnames(mf_summary)<-chem_n
  mf_summary[mf_summary==TRUE]<-1  #converting T,F to 1 and 0
  
  #Calculating tree based on chemical similarity matrix
  clus.hier<-hclust(chem_matrix, method="average")
  tree<-as.phylo(clus.hier)  #convert to phytogenetic tree
  
  #Calculating weighted similarity on binarized databy using {picante:unifrac}
  binary_cw_dist<-unifrac(mf_summary, tree)
  #error control
  binary_cw_dist[is.nan(binary_cw_dist)]<-0
  
  #Calculating chemically weighted dissimilarity on continous data using GUniFrac
  unif_cont_cw<-GUniFrac(file_summary, tree)
  cont_dist_cw<-as.dist(unif_cont_cw$unifracs[,,"d_1"])  #Weighted distance
  cont_dist_cw0.5<-as.dist(unif_cont_cw$unifracs[,,"d_0.5"])  #Weighted distance
  
  #Calculating weighted dissimilarity using fuzzy function
  weighting_tree<-belonging(chem_matrix, standardize = TRUE)
  
  #Since the order of substrates is different between file_summary and weighting_tree, this can be fitted by sorting
  file_summary_sorted<-file_summary[,sort(colnames(file_summary))]
  weighting_tree_sorted<-weighting_tree[sort(rownames(weighting_tree)),sort(colnames(weighting_tree))]
  
  #Fuzzy weighted
  file_summary_cwF <- as.matrix(file_summary_sorted)%*%as.matrix(weighting_tree_sorted)
  cont_dist_cwF<-vegdist(file_summary_cwF, method=meth)
  
  #Calculating dissimilarity on continuous data, w/o chemical weight 
  cont_uw_dist<-vegdist(file_summary, method=meth)
  
  #Calculating default dissimilarity based on presence/absence, w/o chemical weight 
  binary_uw_dist<-vegdist(mf_summary, method="euclidian")
  
  #Distance-based RDA (db-RDA) by vegan
  binary_cw_effect<-capscale(binary_cw_dist~format_data$day*format_data$treatment)
  cont_cw_effect<-capscale(cont_dist_cw~format_data$day*format_data$treatment)
  cont_cw0.5_effect<-capscale(cont_dist_cw0.5~format_data$day*format_data$treatment)
  cont_cwF_effect<-capscale(cont_dist_cwF~format_data$day*format_data$treatment)
  cont_uw_effect<-capscale(cont_uw_dist~format_data$day*format_data$treatment)
  binary_uw_effect<-capscale(binary_uw_dist~format_data$day*format_data$treatment)
  
  #Permutation test
  sig_binary_cw<-anova(binary_cw_effect, permutations=perm, by='terms')
  sig_cont_cw<-anova(cont_cw_effect, permutations=perm, by='terms')
  sig_cont_cw0.5<-anova(cont_cw0.5_effect, permutations=perm, by='terms')
  sig_cont_cwF<-anova(cont_cwF_effect, permutations=perm, by='terms')
  sig_cont_uw<-anova(cont_uw_effect, permutations=perm, by='terms')
  sig_binary_uw<-anova(binary_uw_effect, permutations=perm, by='terms')
  
  
  #Extract the constrained variances
  real_total_chi_binary_cw= binary_cw_effect$CCA$tot.chi + binary_cw_effect$CA$tot.chi
  real_total_chi_cont_cw= cont_cw_effect$CCA$tot.chi + cont_cw_effect$CA$tot.chi
  real_total_chi_cont_cw0.5= cont_cw0.5_effect$CCA$tot.chi + cont_cw0.5_effect$CA$tot.chi
  real_total_chi_cont_cwF= cont_cwF_effect$CCA$tot.chi + cont_cwF_effect$CA$tot.chi
  real_total_chi_cont_uw= cont_uw_effect$CCA$tot.chi + cont_uw_effect$CA$tot.chi
  real_total_chi_binary_uw= binary_uw_effect$CCA$tot.chi + binary_uw_effect$CA$tot.chi
  
  #OUTPUT 1: the explanation power of the mode (by the contrained variance)
  cat("The constrained variance by UW binary data by treatment and their interaction:", binary_uw_effect$CCA$tot.chi/real_total_chi_binary_uw, "\n")
  cat("The constrained variance by CW binary data by treatment and their interaction:", binary_cw_effect$CCA$tot.chi/real_total_chi_binary_cw, "\n")
  cat("The constrained variance by CW continuous data by treatment and their interaction:", cont_cw_effect$CCA$tot.chi/real_total_chi_cont_cw, "\n")
  cat("The constrained variance by CW0.5 continuous data by treatment and their interaction:", cont_cw0.5_effect$CCA$tot.chi/real_total_chi_cont_cw0.5, "\n")
  cat("The constrained variance by CW Fuzzy continuous data by treatment and their interaction:", cont_cwF_effect$CCA$tot.chi/real_total_chi_cont_cwF, "\n")
  cat("The constrained variance by UW continuous data by treatment and their interaction:", cont_uw_effect$CCA$tot.chi/real_total_chi_cont_uw, "\n")
  
  
  #OUTPUT2: The P values
  cat("The P values for day, treatment and ther interactions by UW binary data are", sig_binary_uw$`Pr(>F)`,"\n")
  cat("The P values for day, treatment and ther interactions by CW binary data are", sig_binary_cw$`Pr(>F)`,"\n")
  cat("The P values for day, treatment and ther interactions by CW continuous data are", sig_cont_cw$`Pr(>F)`,"\n")
  cat("The P values for day, treatment and ther interactions by CW0.5 continuous data are", sig_cont_cw0.5$`Pr(>F)`,"\n")
  cat("The P values for day, treatment and ther interactions by CW Fuzzy continuous data are", sig_cont_cwF$`Pr(>F)`,"\n")
  cat("The P values for day, treatment and ther interactions by UW continuous data are", sig_cont_uw$`Pr(>F)`,"\n")
  
  
  #Checking the significance of the chemical information
  chem_dist.mat<-as.matrix(chem_matrix)
  chem_dist.dist<-as.dist(chem_matrix)
  size<-nrow(chem_dist.mat)
  leng<-length(chem_dist.dist)
  
  
  Constrained_binary_cw = binary_cw_effect$CCA$tot.chi/real_total_chi_binary_cw
  Constrained_cont_cwF = cont_cwF_effect$CCA$tot.chi/real_total_chi_cont_cwF
  Constrained_cont_cw = cont_cw_effect$CCA$tot.chi/real_total_chi_cont_cw
  Constrained_cont_cw0.5 = cont_cw0.5_effect$CCA$tot.chi/real_total_chi_cont_cw0.5
  
  Up_i_binary_cw = 0  #Psudo-P value for chemical weight
  Up_i_cont_cwF = 0  #Psudo-P value for chemical weight
  Up_i_cont_cw = 0  #Psudo-P value for chemical weight
  Up_i_cont_cw0.5 = 0  #Psudo-P value for chemical weight
  
  Lp_i_binary_cw = 0  #Psudo-P value for chemical weight
  Lp_i_cont_cwF = 0  #Psudo-P value for chemical weight
  Lp_i_cont_cw = 0  #Psudo-P value for chemical weight
  Lp_i_cont_cw0.5 = 0  #Psudo-P value for chemical weight
  
  
  Constrained_temp<-pforeach(i = 1: perm2, .c=list, .cores=n_cores) ({
    #generating random tree, just shuffling the identity of substrates############
    # z<-sample(1:size, size, replace=FALSE)  #random sequence from 1 to size without overlapping
    # temp <- chem_dist.mat
    # temp_p <- temp[z,] # permute rows 
    # temp_p <- temp_p[,z] # permute columns  With this permutation, we can keep the 0 distance for the same substances
    # colnames(temp_p)<-chem_name
    # row.names(temp_p)<-chem_name
    # temp_d<-as.dist(temp_p)
    ###################################################
    
    #generating random tree, directly shufflting distance matrix##############
    z<-sample(1:leng, leng, replace=FALSE)
    temp_d<-chem_dist.dist
    for(i in 1:leng) temp_d[i]<-chem_dist.dist[z[i]]
    ###############################################################
    
    temp_chem<-hclust(temp_d, method="average")
    temp.ph<-as.phylo(temp_chem)
    plot(temp.ph)
    
    tree<-temp.ph
    chem_matrix<-temp_d
    
    #Calculating through binarized data by using {picante:unifrac}
    dist_binary_cw<-unifrac(mf_summary, tree)
    #error control
    dist_binary_cw[is.nan(dist_binary_cw)]<-0
    
    #Calculating chemically weighted dissimilarity on continous data using GUniFrac
    unif_cont_cw<-GUniFrac(file_summary, tree)
    cont_dist_cw<-as.dist(unif_cont_cw$unifracs[,,"d_1"])  #Weighted distance
    cont_dist_cw0.5<-as.dist(unif_cont_cw$unifracs[,,"d_0.5"])  #Weighted distance
    
    #Calculating weighted dissimilarity using fuzzy function
    weighting_tree<-belonging(chem_matrix, standardize = TRUE)
    
    #Since the order of substrates is different between file_summary and weighting_tree, this can be fitted by sorting
    file_summary_sorted<-file_summary[,sort(colnames(file_summary))]
    weighting_tree_sorted<-weighting_tree[sort(rownames(weighting_tree)),sort(colnames(weighting_tree))]
    
    #Fuzzy weighted
    file_summary_cwF <- as.matrix(file_summary_sorted)%*%as.matrix(weighting_tree_sorted)
    cont_dist_cwF<-vegdist(file_summary_cwF, method=meth)
    
    
    
    #Distance-based RDA (db-RDA) by vegan
    binary_cw_effect<-capscale(dist_binary_cw~format_data$day*format_data$treatment)
    cont_cw_effect<-capscale(cont_dist_cw~format_data$day*format_data$treatment)
    cont_cw0.5_effect<-capscale(cont_dist_cw0.5~format_data$day*format_data$treatment)
    cont_cwF_effect<-capscale(cont_dist_cwF~format_data$day*format_data$treatment)
    
    #Extract the constrained variances
    real_total_chi_binary_cw= binary_cw_effect$CCA$tot.chi + binary_cw_effect$CA$tot.chi
    real_total_chi_cont_cw= cont_cw_effect$CCA$tot.chi + cont_cw_effect$CA$tot.chi
    real_total_chi_cont_cw0.5= cont_cw0.5_effect$CCA$tot.chi + cont_cw0.5_effect$CA$tot.chi
    real_total_chi_cont_cwF= cont_cwF_effect$CCA$tot.chi + cont_cwF_effect$CA$tot.chi
    
    
    #Calcualte constrained values from random tree
    Constrained_binary_cw_temp = binary_cw_effect$CCA$tot.chi/real_total_chi_binary_cw
    Constrained_cont_cwF_temp = cont_cwF_effect$CCA$tot.chi/real_total_chi_cont_cwF
    Constrained_cont_cw_temp = cont_cw_effect$CCA$tot.chi/real_total_chi_cont_cw
    Constrained_cont_cw0.5_temp = cont_cw0.5_effect$CCA$tot.chi/real_total_chi_cont_cw0.5
    
    cbind(Constrained_binary_cw_temp, Constrained_cont_cwF_temp, Constrained_cont_cw_temp, Constrained_cont_cw0.5_temp)
  })
  
  #Counting the case when randomly-generated tree give a greater R2 value than actually estimated tree
  for(i in 1:perm2) {
    if(Constrained_temp[[i]][1] >= Constrained_binary_cw) Up_i_binary_cw = Up_i_binary_cw + 1
    if(Constrained_temp[[i]][2] >= Constrained_cont_cwF) Up_i_cont_cwF = Up_i_cont_cwF + 1  
    if(Constrained_temp[[i]][3] >= Constrained_cont_cw) Up_i_cont_cw = Up_i_cont_cw + 1  
    if(Constrained_temp[[i]][4] >= Constrained_cont_cw0.5) Up_i_cont_cw0.5 = Up_i_cont_cw0.5 + 1  
    
    if(Constrained_temp[[i]][1] <= Constrained_binary_cw) Lp_i_binary_cw = Lp_i_binary_cw + 1
    if(Constrained_temp[[i]][2] <= Constrained_cont_cwF) Lp_i_cont_cwF = Lp_i_cont_cwF + 1  
    if(Constrained_temp[[i]][3] <= Constrained_cont_cw) Lp_i_cont_cw = Lp_i_cont_cw + 1  
    if(Constrained_temp[[i]][4] <= Constrained_cont_cw0.5) Lp_i_cont_cw0.5 = Lp_i_cont_cw0.5 + 1  
    
  }
  #Output2 
  cat("The probability with which the # of the constrained variation from permutated trees >= observed value is more extreme under null hypothesis\n")
  cat("binary_cw\t", (Up_i_binary_cw+1)/(perm2+1), "\n") 
  cat("continuous fuzzy weighted\t", (Up_i_cont_cwF+1)/(perm2+1), "\n") 
  cat("continuous GuniFrac weighted\t", (Up_i_cont_cw+1)/(perm2+1), "\n") 
  cat("continuous GuniFrac 0.5 weighted\t", (Up_i_cont_cw0.5+1)/(perm2+1), "\n") 
  cat("The probability with which the # of the constrained variation from permutated trees <= observed value is more extreme under null hypothesis\n")
  cat("binary_cw\t", (Lp_i_binary_cw+1)/(perm2+1), "\n") 
  cat("continuous fuzzy weighted\t", (Lp_i_cont_cwF+1)/(perm2+1), "\n") 
  cat("continuous GuniFrac weighted\t", (Lp_i_cont_cw+1)/(perm2+1), "\n") 
  cat("continuous GuniFrac 0.5 weighted\t", (Lp_i_cont_cw0.5+1)/(perm2+1), "\n") 
  
  
}


#Function, which calculate quantile-based multifunctionality (MF) and with or without weighting by chemical_similarity, then simple linear analysis, on microcosm experiments, test the regression of MF by the degree of COG reduction
###############list of parameters#######################
#data01~data04: data sources of ecopalte color patterns (list format), which are separately saved
#format_data01~format_data04: data-specific metadata (sampling dates, treatments, etc)
#format_dataF: common metadata for data01~data04 
#chem_n: list of chemical names
#thres: the threshold for calculating quantile-based MF
#chem_max: the similarity matrix for ecoplate substrates
#perm: number of permulations for chemical similarity tree
#n_cores: number of CPU cores used for parallel permutations
#show_tree: Whether the permutated chemical similarity tree is shown or not
Chem_Disim_MF_microcosm <- function(data01=max_non_integ_max_01_microcosm,data02=max_non_integ_max_02_microcosm,data03=max_non_integ_max_03_microcosm, data04=max_non_integ_max_04_microcosm, format_data01=name_microcosm01,format_data02=name_microcosm02,format_data03=name_microcosm03,format_data04=name_microcosm04, format_dataF=format_microcosm, chem_n=chem_name, thres=0.1, chem_matrix=fp.dist_c, perm=999, ncores=2, show_tree=FALSE)
{
  ####IMPORTANT: The comments within this function are simplified.Check every comment on the function above "Chem_Disim_MF_xitou"########## 
  
  #Loading all data, preparing the lists to stock information
  data<-list()
  format_data<-list()
  pd.result<-list()
  ts<-list()
  data[[1]]<-data01
  data[[2]]<-data02
  data[[3]]<-data03
  data[[4]]<-data04
  format_data[[1]]<-format_data01
  format_data[[2]]<-format_data02
  format_data[[3]]<-format_data03
  format_data[[4]]<-format_data04
  
  #Calculating tree based on chemical similarity matrix
  clus.hier<-hclust(chem_matrix, method="average")
  tree<-as.phylo(clus.hier)  #convert to phytogenetic tree
  
  for(k in 1:4) {
    #Data from four different dates should be combined into the single list
    file_summary <- data[[k]][[1]]
    for(i in 2:length(data[[k]])) file_summary <- rbind(file_summary, data[[k]][[i]])
    file_summary <-file_summary[,c(-1)] #exclusing water data (control, zero)
    row.names(file_summary)<-format_data[[k]]$sample
    colnames(file_summary)<-chem_n
    file_summary
    
    #Convert data into binary values by quantile (quantile-based multifunctionality) for each experiment dates
    mf<-list()
    for(i in 1:31) mf[[i]] <-(file_summary[,i] > quantile(file_summary[,i], thres))
    ts[[k]]<-as.data.frame(mf)
    row.names(ts[[k]])<-format_data[[k]]$sample
    colnames(ts[[k]])<-chem_n
    ts[[k]][ts[[k]]==TRUE]<-1  #converting T,F to 1 and 0
    
    #Using function in picante (phylogenetic community analysis)
    pd.result[[k]] <- pd(ts[[k]], tree, include.root=T)
    pd.result[[k]]<-cbind(format_data[[k]], pd.result[[k]])   #Combined with sample name
    
    pd.result[[k]]<-cbind(pd.result[[k]], format_dataF$COG_reduction[pd.result[[k]]$composition]) #Combine with COG reduction information
    colnames(pd.result[[k]])<-c("sample_data", "composition", "Chem_MF", "MF", "COG_reduction")
  }#end of for k
  
  #The result should be combined into a single matrix
  pd_summary<-pd.result[[1]]
  for(k in 2:4) pd_summary<-rbind(pd_summary, pd.result[[k]])
  
  model_nw<-lm(MF~COG_reduction, pd_summary)   #liner model of explaining MF by COG reduction for data without chemical similarity
  #model_nw <-glm(MF~COG_reduction, pd_summary, family=poisson)
  #summary(model_nw)
 
  model_w<-lm(Chem_MF~COG_reduction, pd_summary) #the model for chemical-similarity-weighted data
  
  ###For graphical presentation
  #for regression lines
  xx <- c(0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400)
  #regression lines equations
  func_nw<-function(x) model_nw$coefficients[1] + model_nw$coefficients[2]*x
  func_w<-function(x) model_w$coefficients[1] + model_w$coefficients[2]*x
  
  plot(pd_summary$COG_reduction, pd_summary$MF, xlim=c(0, 1400), ylim=c(0, max(pd_summary$MF)))
  lines(xx, func_nw(xx), col=1, lty=1)
  par(new=T)
  plot(pd_summary$COG_reduction, pd_summary$Chem_MF, xlim=c(0, 1400),ylim=c(0, max(pd_summary$MF)),col=3, ylab="")
  lines(xx, func_w(xx), col=3, lty=1)
  
  
  #Output1
  cat("Non-weighted MF can be explained by reduction of COG?\n")
  print(summary(model_nw))  #show summary of the statistics
  cat("Weighted MF can be explained by reduction of COG?\n")
  print(summary(model_w))   #show summary of the statistics
  
  #To calculate R2 values from ANOVA results
  anv_w<-anova(model_w)
  R2_w=(anv_w$`Sum Sq`[1])/(anv_w$`Sum Sq`[1]+anv_w$`Sum Sq`[2])
  #R2_w: the explanation power of the linear model using chemical-similarity-information
  
  #Checking the significance of the chemical information
  chem_dist.mat<-as.matrix(chem_matrix)
  chem_dist.dist<-as.dist(chem_matrix)
  size<-nrow(chem_dist.mat)
  leng<-length(chem_dist.dist)
  
  p_iU = 0  #Psudo-P value for chemical weight
  p_iL = 0  #Psudo-P value for chemical weight
  
  #perm=10
  #show_tree=TRUE
  #for parallel calculation
  R2_temp<-pforeach(i = 1: perm, .c=list, .cores=ncores) ({
    
    #generating random tree, just shuffling the identity of substrates############
    # z<-sample(1:size, size, replace=FALSE)  #random sequence from 1 to size without overlapping
    # temp <- chem_dist.mat
    # temp_p <- temp[z,] # permute rows 
    # temp_p <- temp_p[,z] # permute columns  With this permutation, we can keep the 0 distance for the same substances
    # colnames(temp_p)<-chem_name
    # row.names(temp_p)<-chem_name
    # temp_d<-as.dist(temp_p)
    ###################################################
    
    #generating random tree, directly shufflting distance matrix##############
    z<-sample(1:leng, leng, replace=FALSE)
    temp_d<-chem_dist.dist
    for(i in 1:leng) temp_d[i]<-chem_dist.dist[z[i]]
    ###############################################################
    
    temp_chem<-hclust(temp_d, method="average")
    temp.ph<-as.phylo(temp_chem)
    if(show_tree==TRUE) plot(temp.ph)
    
    #Using function in picante (phylogenetic community analysis)
    for(k in 1:4) {
      pd.result[[k]] <- pd(ts[[k]], temp.ph, include.root=T)
      pd.result[[k]]<-cbind(format_data[[k]], pd.result[[k]])   #Combined with sample name
      
      pd.result[[k]]<-cbind(pd.result[[k]], format_dataF$COG_reduction[pd.result[[k]]$composition]) #Combine with COG reduction
      colnames(pd.result[[k]])<-c("sample_data", "composition", "Chem_MF", "MF", "COG_reduction")
    }#end of for k
    
    pd_summary_temp<-pd.result[[1]]
    for(k in 2:4) pd_summary_temp<-rbind(pd_summary_temp, pd.result[[k]])
    
    model_w_temp<-lm(Chem_MF~COG_reduction, pd_summary_temp)  #make a linear model for permutated tree
    
    anv_w_temp<-anova(model_w_temp)
    R2_w_temp=(anv_w_temp$`Sum Sq`[1])/(anv_w_temp$`Sum Sq`[1]+anv_w_temp$`Sum Sq`[2])
  }) #end of parallel calculation
  
  for(i in 1: perm)  {
    if(R2_temp[i] >= R2_w) p_iU = p_iU + 1  #Counting the case when randomly-generated tree give a greater R2 value than actually estimated tree
    if(R2_temp[i] <= R2_w) p_iL = p_iL + 1
  }
  
  
  #Output2 
  cat("The probability with which the # of R2 values from permutated tree >= observed R2 value is more extreme under null hypothesis\n (P-Value for the chemical tree)\n")
  print((p_iU+1)/(perm+1))
  cat("The probability with which the # of R2 values from permutated tree <= observed R2 value is more extreme under null hypothesis\n (P-Value for the chemical tree)\n")
  print((p_iL+1)/(perm+1))
  
  plot(pd_summary$COG_reduction, pd_summary$MF, xlim=c(0, 1400), ylim=c(0, max(pd_summary$MF)))
  lines(xx, func_nw(xx), col=1, lty=1)
  par(new=T)
  plot(pd_summary$COG_reduction, pd_summary$Chem_MF, xlim=c(0, 1400),ylim=c(0, max(pd_summary$MF)),col=3, ylab="")
  lines(xx, func_w(xx), col=3, lty=1)
  
}
}
##############################################End of definitons of functions###########################################


############################################ii)Analyses for microcosm experiments######################################
{
  
  #The case with taking averages of triplicate of the data at the final day: threshold values range from 0.1 to 0.9
  #test with small permutation number (9)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.1, chem_matrix=fp.dist_a, perm=9)
  
  #Using Tree-a with different threshold T from 0.1 to 0.9
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.1, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.2, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.3, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.4, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.5, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.6, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.7, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.8, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.9, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  
  #Using Tree-b with different threshold T from 0.1 to 0.9
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.1, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.2, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.3, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.4, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.5, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.6, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.7, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.8, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.9, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  
  #Using Tree-c with different threshold T from 0.1 to 0.9
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.1, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.2, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.3, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.4, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.5, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.6, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.7, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.8, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=ave_non_integ_max_01_microcosm,data02=ave_non_integ_max_02_microcosm,data03=ave_non_integ_max_03_microcosm, data04=ave_non_integ_max_04_microcosm, thres=0.9, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  
  #The case with taking maximum of triplicate: threshold values range from 0.1 to 0.9
  #Using Tree-a
  Chem_Disim_MF_microcosm(thres=0.1, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.2, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.3, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.4, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.5, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.6, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.7, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.8, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.9, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  
  #Using Tree-b
  Chem_Disim_MF_microcosm(thres=0.1, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.2, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.3, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.4, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.5, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.6, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.7, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.8, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.9, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  
  #Using Tree-c
  Chem_Disim_MF_microcosm(thres=0.1, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.2, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.3, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.4, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.5, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.6, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.7, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.8, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(thres=0.9, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  
  #The case with taking minimum of triplicate: threshold values range from 0.1 to 0.9
  #Using Tree-a
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.1, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.2, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.3, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.4, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.5, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.6, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.7, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.8, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.9, chem_matrix=fp.dist_a, perm=999, ncores = 20)
  
  #Using Tree-b
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.1, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.2, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.3, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.4, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.5, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.6, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.7, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.8, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.9, chem_matrix=fp.dist_b, perm=999, ncores = 20)
  
  #Using Tree-c
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.1, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.2, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.3, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.4, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.5, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.6, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.7, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.8, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  Chem_Disim_MF_microcosm(data01=min_non_integ_max_01_microcosm,data02=min_non_integ_max_02_microcosm,data03=min_non_integ_max_03_microcosm, data04=min_non_integ_max_04_microcosm, thres=0.9, chem_matrix=fp.dist_c, perm=999, ncores = 20)
  
  
}
########################################End of analyses for microcosm experiments######################################


#################################### iii)Analyses for Xitou data#######################################################
{
#PERMANOVA comparision on Xitou data
#Note that the results from this function and the related functions below of PERMANOVA were not shown in the manuscript.
  
{
#Taking integraton with the averages of the triplicate  
w_s_adonis_xitou(data=ave_integ_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_a, perm2=999, n_cores=20)
w_s_adonis_xitou(data=ave_integ_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_b, perm2=999, n_cores=20)
w_s_adonis_xitou(data=ave_integ_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_c, perm2=999, n_cores=20)

#Taking integraton with the maximum of the triplicate
w_s_adonis_xitou(data=max_integ_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_a, perm2=999, n_cores=20)
w_s_adonis_xitou(data=max_integ_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_b, perm2=999, n_cores=20)
w_s_adonis_xitou(data=max_integ_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_c, perm2=999, n_cores=20)

#Taking integraton with the minimum of the triplicate
w_s_adonis_xitou(data=min_integ_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_a, perm2=999, n_cores=20)
w_s_adonis_xitou(data=min_integ_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_b, perm2=999, n_cores=20)
w_s_adonis_xitou(data=min_integ_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_c, perm2=999, n_cores=20)

#Taking day at maximum with the averages of the triplicate
w_s_adonis_xitou(data=ave_non_integ_max_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_a, perm2=999, n_cores=20)
w_s_adonis_xitou(data=ave_non_integ_max_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_b, perm2=999, n_cores=20)
w_s_adonis_xitou(data=ave_non_integ_max_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_c, perm2=999, n_cores=20)

#Taking day at maximum with the maximum of the triplicate
w_s_adonis_xitou(data=max_non_integ_max_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_a, perm2=999, n_cores=20)
w_s_adonis_xitou(data=max_non_integ_max_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_b, perm2=999, n_cores=20)
w_s_adonis_xitou(data=max_non_integ_max_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_c, perm2=999, n_cores=20)

#Taking day at maximum with the minimum of the triplicate
w_s_adonis_xitou(data=min_non_integ_max_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_a, perm2=999, n_cores=20)
w_s_adonis_xitou(data=min_non_integ_max_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_b, perm2=999, n_cores=20)
w_s_adonis_xitou(data=min_non_integ_max_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_c, perm2=999, n_cores=20)

#Taking final day with the averages of the triplicate
w_s_adonis_xitou(data=ave_non_integ_final_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_a, perm2=999, n_cores=20)
w_s_adonis_xitou(data=ave_non_integ_final_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_b, perm2=999, n_cores=20)
w_s_adonis_xitou(data=ave_non_integ_final_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_c, perm2=999, n_cores=20)

#Taking final day with the maximum of the triplicate
w_s_adonis_xitou(data=max_non_integ_final_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_a, perm2=999, n_cores=20)
w_s_adonis_xitou(data=max_non_integ_final_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_b, perm2=999, n_cores=20)
w_s_adonis_xitou(data=max_non_integ_final_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_c, perm2=999, n_cores=20)

#Taking final day with the minimum of the triplicate
w_s_adonis_xitou(data=min_non_integ_final_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_a, perm2=999, n_cores=20)
w_s_adonis_xitou(data=min_non_integ_final_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_b, perm2=999, n_cores=20)
w_s_adonis_xitou(data=min_non_integ_final_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_c, perm2=999, n_cores=20)


}


#Linear model with/without chemical similarity weighting Multifunctionality (MF) on Xitou data
{

#Taking integration with the averages of triplicate, with the threshold T = 0.9, 0.7, and 0.5
Chem_Disim_MF_xitou(data=ave_integ_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_integ_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_integ_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_integ_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_integ_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20) 
Chem_Disim_MF_xitou(data=ave_integ_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_integ_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_integ_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)    
Chem_Disim_MF_xitou(data=ave_integ_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)  

#Taking integration with the maximum of triplicate, with the threshold T = 0.9, 0.7, and 0.5
Chem_Disim_MF_xitou(data=max_integ_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_integ_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_integ_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_integ_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_integ_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20) 
Chem_Disim_MF_xitou(data=max_integ_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_integ_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_integ_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)    
Chem_Disim_MF_xitou(data=max_integ_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)  

#Taking integration with the minimum of triplicate, with the threshold T = 0.9, 0.7, and 0.5
Chem_Disim_MF_xitou(data=min_integ_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_integ_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_integ_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_integ_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_integ_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20) 
Chem_Disim_MF_xitou(data=min_integ_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_integ_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_integ_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_integ_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)  

#Taking day at mamimum with the averages of triplicate, with the threshold T = 0.9, 0.7, and 0.5
Chem_Disim_MF_xitou(data=ave_non_integ_max_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_non_integ_max_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_non_integ_max_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)
Chem_Disim_MF_xitou(data=ave_non_integ_max_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_non_integ_max_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_non_integ_max_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)
Chem_Disim_MF_xitou(data=ave_non_integ_max_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_non_integ_max_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_non_integ_max_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)

#Taking day at mamimum with the maximum of triplicate, with the threshold T = 0.9, 0.7, and 0.5
Chem_Disim_MF_xitou(data=max_non_integ_max_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_non_integ_max_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_non_integ_max_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)
Chem_Disim_MF_xitou(data=max_non_integ_max_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_non_integ_max_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_non_integ_max_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)
Chem_Disim_MF_xitou(data=max_non_integ_max_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_non_integ_max_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_non_integ_max_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)

#Taking day at mamimum with the minimum of triplicate, with the threshold T = 0.9, 0.7, and 0.5
Chem_Disim_MF_xitou(data=min_non_integ_max_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_non_integ_max_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_non_integ_max_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)
Chem_Disim_MF_xitou(data=min_non_integ_max_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_non_integ_max_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_non_integ_max_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)
Chem_Disim_MF_xitou(data=min_non_integ_max_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_non_integ_max_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_non_integ_max_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)

#Taking final day with the averages of triplicate, with the threshold T = 0.9, 0.7, and 0.5
Chem_Disim_MF_xitou(data=ave_non_integ_final_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_non_integ_final_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_non_integ_final_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)
Chem_Disim_MF_xitou(data=ave_non_integ_final_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_non_integ_final_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_non_integ_final_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)
Chem_Disim_MF_xitou(data=ave_non_integ_final_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_non_integ_final_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=ave_non_integ_final_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)

#Taking final day with the maximum of triplicate, with the threshold T = 0.9, 0.7, and 0.5
Chem_Disim_MF_xitou(data=max_non_integ_final_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_non_integ_final_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_non_integ_final_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)
Chem_Disim_MF_xitou(data=max_non_integ_final_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_non_integ_final_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_non_integ_final_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)
Chem_Disim_MF_xitou(data=max_non_integ_final_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_non_integ_final_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=max_non_integ_final_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)

#Taking final day with the minimum of triplicate, with the threshold T = 0.9, 0.7, and 0.5
Chem_Disim_MF_xitou(data=min_non_integ_final_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_non_integ_final_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_non_integ_final_20142015_Xitou_Eco, thres=0.9, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)
Chem_Disim_MF_xitou(data=min_non_integ_final_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_non_integ_final_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_non_integ_final_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)
Chem_Disim_MF_xitou(data=min_non_integ_final_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_a,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_non_integ_final_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_b,perm=999, n_cores=20)  
Chem_Disim_MF_xitou(data=min_non_integ_final_20142015_Xitou_Eco, thres=0.5, format_data=format_d,chem_matrix=fp.dist_c,perm=999, n_cores=20)

}


#distance-based RDA analysis with/without chemical simiralirty weighting on xitou data
{
#Taking integration with the averages of triplicate, with the threshold T = 0.9, 0.7, and 0.5
W_s_dbRDA(ave_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.9, perm2=999, n_cores=20)
W_s_dbRDA(ave_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.7, perm2=999, n_cores=20)
W_s_dbRDA(ave_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.5, perm2=999, n_cores=20)
W_s_dbRDA(ave_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.9, perm2=999, n_cores=20)
W_s_dbRDA(ave_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.7, perm2=999, n_cores=20)
W_s_dbRDA(ave_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.5, perm2=999, n_cores=20)
W_s_dbRDA(ave_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.9, perm2=999, n_cores=20)
W_s_dbRDA(ave_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.7, perm2=999, n_cores=20)
W_s_dbRDA(ave_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.5, perm2=999, n_cores=20)

#Taking integration with the maximum of triplicate, with the threshold T = 0.9, 0.7, and 0.5
W_s_dbRDA(max_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.9, perm2=999, n_cores=20)
W_s_dbRDA(max_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.7, perm2=999, n_cores=20)
W_s_dbRDA(max_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.5, perm2=999, n_cores=20)
W_s_dbRDA(max_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.9, perm2=999, n_cores=20)
W_s_dbRDA(max_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.7, perm2=999, n_cores=20)
W_s_dbRDA(max_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.5, perm2=999, n_cores=20)
W_s_dbRDA(max_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.9, perm2=999, n_cores=20)
W_s_dbRDA(max_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.7, perm2=999, n_cores=20)
W_s_dbRDA(max_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.5, perm2=999, n_cores=20)

#Taking integration with the minimum of triplicate, with the threshold T = 0.9, 0.7, and 0.5
W_s_dbRDA(min_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.9, perm2=999, n_cores=20)
W_s_dbRDA(min_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.7, perm2=999, n_cores=20)
W_s_dbRDA(min_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.5, perm2=999, n_cores=20)
W_s_dbRDA(min_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.9, perm2=999, n_cores=20)
W_s_dbRDA(min_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.7, perm2=999, n_cores=20)
W_s_dbRDA(min_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.5, perm2=999, n_cores=20)
W_s_dbRDA(min_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.9, perm2=999, n_cores=20)
W_s_dbRDA(min_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.7, perm2=999, n_cores=20)
W_s_dbRDA(min_integ_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.5, perm2=999, n_cores=20)

#Taking day at maximum with the averages of triplicate, with the threshold T = 0.9, 0.7, and 0.5
W_s_dbRDA(ave_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.9, perm2=999)
W_s_dbRDA(ave_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.7, perm2=999)
W_s_dbRDA(ave_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.5, perm2=999)
W_s_dbRDA(ave_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.9, perm2=999)
W_s_dbRDA(ave_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.7, perm2=999)
W_s_dbRDA(ave_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.5, perm2=999)
W_s_dbRDA(ave_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.9, perm2=999)
W_s_dbRDA(ave_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.7, perm2=999)
W_s_dbRDA(ave_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.5, perm2=999)

#Taking day at maximum with the maximum of triplicate, with the threshold T = 0.9, 0.7, and 0.5
W_s_dbRDA(max_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.9, perm2=999)
W_s_dbRDA(max_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.7, perm2=999)
W_s_dbRDA(max_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.5, perm2=999)
W_s_dbRDA(max_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.9, perm2=999)
W_s_dbRDA(max_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.7, perm2=999)
W_s_dbRDA(max_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.5, perm2=999)
W_s_dbRDA(max_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.9, perm2=999)
W_s_dbRDA(max_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.7, perm2=999)
W_s_dbRDA(max_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.5, perm2=999)

#Taking day at maximum with the minimum of triplicate, with the threshold T = 0.9, 0.7, and 0.5
W_s_dbRDA(min_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.9, perm2=999)
W_s_dbRDA(min_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.7, perm2=999)
W_s_dbRDA(min_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.5, perm2=999)
W_s_dbRDA(min_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.9, perm2=999)
W_s_dbRDA(min_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.7, perm2=999)
W_s_dbRDA(min_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.5, perm2=999)
W_s_dbRDA(min_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.9, perm2=999)
W_s_dbRDA(min_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.7, perm2=999)
W_s_dbRDA(min_non_integ_max_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.5, perm2=999)

#Taking final day with the averages of triplicate, with the threshold T = 0.9, 0.7, and 0.5
W_s_dbRDA(ave_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.9, perm2=999)
W_s_dbRDA(ave_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.7, perm2=999)
W_s_dbRDA(ave_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.5, perm2=999)
W_s_dbRDA(ave_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.9, perm2=999)
W_s_dbRDA(ave_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.7, perm2=999)
W_s_dbRDA(ave_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.5, perm2=999)
W_s_dbRDA(ave_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.9, perm2=999)
W_s_dbRDA(ave_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.7, perm2=999)
W_s_dbRDA(ave_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.5, perm2=999)

#Taking final day with the maximum of triplicate, with the threshold T = 0.9, 0.7, and 0.5
W_s_dbRDA(max_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.9, perm2=999)
W_s_dbRDA(max_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.7, perm2=999)
W_s_dbRDA(max_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.5, perm2=999)
W_s_dbRDA(max_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.9, perm2=999)
W_s_dbRDA(max_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.7, perm2=999)
W_s_dbRDA(max_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.5, perm2=999)
W_s_dbRDA(max_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.9, perm2=999)
W_s_dbRDA(max_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.7, perm2=999)
W_s_dbRDA(max_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.5, perm2=999)

#Taking final day with the minimum of triplicate, with the threshold T = 0.9, 0.7, and 0.5
W_s_dbRDA(min_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.9, perm2=999)
W_s_dbRDA(min_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.7, perm2=999)
W_s_dbRDA(min_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_a, thres=0.5, perm2=999)
W_s_dbRDA(min_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.9, perm2=999)
W_s_dbRDA(min_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.7, perm2=999)
W_s_dbRDA(min_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_b, thres=0.5, perm2=999)
W_s_dbRDA(min_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.9, perm2=999)
W_s_dbRDA(min_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.7, perm2=999)
W_s_dbRDA(min_non_integ_final_20142015_Xitou_Eco, chem_matrix=fp.dist_c, thres=0.5, perm2=999)

}
}
########################################End of analyses for Xitou data#################################################

########################################End of PART3###################################################################

###########################PART4 Additional analysis for SI############################################################

#######################i) Checking the distribution of statistical values from randomly-generated similarity tree################
###########This part is for Fig.S6

#The modified function from "Chem_Disim_MF_xitou"
Chem_Disim_MF_xitou_hist <- function(data=ave_integ_20142015_Xitou_Eco, thres=0.9, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, perm=9, n_cores=1) 
{
  #Loading all data (24 samples), some modification is necessary for the data size but in our case there are 24 samples in total.
  file_summary <- data[[1]]
  for(i in 2:length(data)) file_summary <- rbind(file_summary, data[[i]])  #data are combined in to a single file
  file_summary <-file_summary[,c(-1)] #exclusing water data (control, zero)
  row.names(file_summary)<-format_data$sample  #Add name (sample ID)
  colnames(file_summary)<-chem_n #Add name (substrate name)
  
  #Convert data into binary values by quantile (quantile-based multifunctionality)
  mf<-list()  #prepare the output list
  for(i in 1:31) mf[[i]] <-(file_summary[,i] > quantile(file_summary[,i], thres))  # return true if the color value qunatile is greater than threshold
  M_mf<-as.data.frame(mf)  #converted into dataframe
  row.names(M_mf)<-format_data$sample
  colnames(M_mf)<-chem_n
  M_mf[M_mf==TRUE]<-1  #converting T,F to 1 and 0; ts is now dataframe to represent the binary matrix (E_BC in Figure 3)
  
  #Calculating tree based on chemical similarity matrix
  clus.hier<-hclust(chem_matrix, method="average")
  tree<-as.phylo(clus.hier)  #convert to phytogenetic tree format; this process is necessary to use phylogenetic diversity (PD) function from picante for chemical similarity tree.
  
  #Using function in picante (phylogenetic community analysis)
  pd.result <- pd(M_mf, tree, include.root=T)  #Claculate chemical-similarity-weighted multifunctionality, and non-weighted (raw) multifunctionality
  colnames(pd.result)<-c("Chem_MF", "MF")   #add name
  pd.result<-cbind(format_data, pd.result)  #incorporate sample data (day & treatment)
  
  #Plot chemical-similarity-weighted multifunctionality: Normal (control) vs Treatment (trenching)
  plot(subset(pd.result, treatment=="N")$day, subset(pd.result, treatment=="N")$Chem_MF, ylim=c(0,6))
  par(new=T)
  plot(subset(pd.result, treatment=="T")$day, subset(pd.result, treatment=="T")$Chem_MF, col=3,ylim=c(0,6))
  
  
  ############Conduct a simple linear regression model#####################
  #MF is chemical simmilarity independent multifuncitonality (0-31)
  model_nw<-lm(MF~day*treatment, pd.result)
  #PD is chemical-similarity-weighted multifuncitonality 
  model_w<-lm(Chem_MF~day*treatment, pd.result)
  
  #Output1
  #cat("Non-weighted MF can be explained by day, treatment, and their interaction?\n")
  #print(summary(model_nw))
  #cat("Weighted MF can be explained by day, treatment, and their interaction?\n")
  #print(summary(model_w))
  
  
  ####Importantly, it is not sure if the explanation power (R2) of the linear regression on chemical-similarity weighted MF is statistically meaningfull or not due to the uncertainty of the chemical-similarity definition. Therefore, it is necessary to conduct permutation test through comparing R2 value from "model_w" and R2 values coming from the analysis using "randomly-generated chemical similarity tree".
  
  #To calculate R2 values from ANOVA results
  anv_w<-anova(model_w)
  R2_w=(anv_w$`Sum Sq`[1]+anv_w$`Sum Sq`[2]+anv_w$`Sum Sq`[3])/(anv_w$`Sum Sq`[1]+anv_w$`Sum Sq`[2]+anv_w$`Sum Sq`[3]+anv_w$`Sum Sq`[4])
  
  
  #######Checking the significance of the chemical information (permutation test)
  chem_dist.mat<-as.matrix(chem_matrix)   #format arrangement on data
  chem_dist.dist<-as.dist(chem_matrix)    #format arrangement on data
  size<-nrow(chem_dist.mat)
  leng<-length(chem_dist.dist)
  
  p_iU = 0  #Psudo-P value for chemical weight
  p_iL = 0  #Psudo-P value for chemical weight
  
  #############Use parallel calculations ".cores" specifies the number of CPU cores to be used.
  R2_temp<-pforeach(i = 1: perm, .c=list, .cores=n_cores) ({
    
    ###generating random tree, just shuffling the identity of substrates (we decided not to use this)############
    # z<-sample(1:size, size, replace=FALSE)  #random sequence from 1 to size without overlapping
    # temp <- chem_dist.mat
    # temp_p <- temp[z,] # permute rows 
    # temp_p <- temp_p[,z] # permute columns  With this permutation, we can keep the 0 distance for the same substances
    # colnames(temp_p)<-chem_name
    # row.names(temp_p)<-chem_name
    # temp_d<-as.dist(temp_p)
    #############################################################################################################
    
    ###generating random tree, directly shufflting distance (similarity) matrix###############
    z<-sample(1:leng, leng, replace=FALSE) #re-sampling without replacement
    temp_d<-chem_dist.dist
    for(i in 1:leng) temp_d[i]<-chem_dist.dist[z[i]]
    #############################################################################
    
    temp_chem<-hclust(temp_d, method="average")    #clustering on randomly-generated similarity matrix
    temp.ph<-as.phylo(temp_chem)  #randomly-generated tree
    plot(temp.ph)
    
    ####Repeat the identical analysis using randomly-genreated tree
    pd.result_temp <- pd(M_mf, temp.ph, include.root=T)
    pd.result_temp<-cbind(format_data, pd.result_temp)
    model_w_temp<-lm(PD~day*treatment, pd.result_temp)
    
    #Calculate R2 value
    anv_w_temp<-anova(model_w_temp)
    R2_w_temp=(anv_w_temp$`Sum Sq`[1]+anv_w_temp$`Sum Sq`[2]+anv_w_temp$`Sum Sq`[3])/(anv_w_temp$`Sum Sq`[1]+anv_w_temp$`Sum Sq`[2]+anv_w_temp$`Sum Sq`[3]+anv_w_temp$`Sum Sq`[4])
    
    #if(R2_w_temp > R2_w) p_i = p_i + 1  #Counting the case when randomly-generated tree give a greater R2 value than actually estimated tree
  })
  
  R2_temp
  #Count the cases when R2 value from random-tree is greater or smaller than R2_w in order to calculate Psuedo-P value of permutation
  #for(i in 1: perm) if(R2_temp[[i]]>= R2_w) p_iU = p_iU + 1
  #for(i in 1: perm) if(R2_temp[[i]]<= R2_w) p_iL = p_iL + 1
  
  #Output2 of permutation results
  #cat("The probability with which the # of R2 values from permutated tree >= observed R2 value is more extreme under null hypothesis\n (P-Value for the chemical tree)\n")
  #print((p_iU+1)/(perm+1))
  #cat("The probability with which the # of R2 values from permutated tree <= observed R2 value is more extreme under null hypothesis\n (P-Value for the chemical tree)\n")
  #print((p_iL+1)/(perm+1))
  
}

#for tree a. b. and c, respectively
random_dist_a<-Chem_Disim_MF_xitou_hist(data=min_non_integ_max_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_a,perm=999,n_cores=20)  
random_dist_b<-Chem_Disim_MF_xitou_hist(data=min_non_integ_max_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_b,perm=999,n_cores=20)  
random_dist_c<-Chem_Disim_MF_xitou_hist(data=min_non_integ_max_20142015_Xitou_Eco, thres=0.7, format_data=format_d,chem_matrix=fp.dist_c,perm=999,n_cores=20)  

#Output is a list
is.list(random_dist_a)
length(random_dist_a)
#Then, need to convert to vector
hist_a<-vector()
hist_a[1]
for(i in 1: length(random_dist_a)) hist_a[i]<-random_dist_a[[i]]
hist_b<-vector()
hist_b[1]
for(i in 1: length(random_dist_b)) hist_b[i]<-random_dist_b[[i]]
hist_c<-vector()
hist_c[1]
for(i in 1: length(random_dist_c)) hist_c[i]<-random_dist_c[[i]]

#Plot the results
b_break =c(0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25)
hist(hist_a, col = "#ff00ff40", border = "#ff00ff", breaks = b_break, xlim=c(0.0, 0.25), freq=FALSE)
hist(hist_b, col = "#0000ff40", border = "#0000ff", breaks = b_break, xlim=c(0.0, 0.25), freq=FALSE)
hist(hist_c, col = "#00ff0040", border = "#00ff00", breaks = b_break, xlim=c(0.0, 0.25), freq=FALSE)
##########End of Checkign the distibution of statistical values#####################


###################ii)Comparing the chemical dissimilarity matrix and in situ color development dissimiarlity matrix#############

#Function, which calculate the pairwise correlation between all substrate color developments in ecoplates and the chemical-similarity matrix. This is for comparing the performance of three chemical similarity trees for Fig. S5
###############list of parameters#######################
#data: data from ecoplates 
#meth: method for calculating correlation
#chem_n: list of chemical names
#chem_dist: distance matrix of chemical similarity
#perm: number of permutation for Mantel correlation between ecoplate patterns and chemical-similarity trees
distance_substrate <- function(data=max_integ_20142015_Xitou_Eco, meth="pearson", chem_n=chem_name, chem_dist, perm=19999, tree=TRUE) {
  
  #First compile all of the 24 day's data
  file_summary <- data[[1]]
  for(i in 2:length(data)) file_summary <- rbind(file_summary, data[[i]])
  file_summary <-file_summary[,c(-1)] #exclusing water data (control, zero)
  file_summary[file_summary <0]<-0 #convert negative values to zero
  file_summary <- as.data.frame(file_summary) #converting to dataframe
  colnames(file_summary)<-chem_n
  
  cor_substrate <- cor(file_summary, method=meth) #calcualte all pair-wise correlation
  
  cor_substrate[cor_substrate <0]<-0   #neglecting negative correlation (assumption for similarity based on ecoplate)
  dist_substrate<-1-cor_substrate  #converting from high correlation to shorter distance (smaller dissimilarity)
  insitu_dist<-as.dist(dist_substrate) #converting to distance matrix format in R
  clus.hier_insitu<-hclust(insitu_dist, method="average")  #clustering
  tree_insitu.ph<-as.phylo(clus.hier_insitu) #convert to phylogenetic tree format
  if(tree==TRUE) plot(tree_insitu.ph)  #plot the tree
  else plot(insitu_dist, chem_dist)
  mantel(insitu_dist, chem_dist, permutations=perm)  #Matel correlation test on chemical-similarity matrix and and ecoplate correlation matrix
  
}


#Applying to Xitou best data
distance_substrate(data=max_integ_20142015_Xitou_Eco, meth="pearson", chem_n=chem_name, chem_dist=fp.dist_a, perm=19999) 
distance_substrate(data=max_integ_20142015_Xitou_Eco, meth="pearson", chem_n=chem_name, chem_dist=fp.dist_b, perm=19999) 
distance_substrate(data=max_integ_20142015_Xitou_Eco, meth="pearson", chem_n=chem_name, chem_dist=fp.dist_c, perm=19999) 

#Applying to microcosm best data
distance_substrate(data=max_non_integ_max_total_microcosm, meth="pearson", chem_n=chem_name, chem_dist=fp.dist_a, perm=19999) 
distance_substrate(data=max_non_integ_max_total_microcosm, meth="pearson", chem_n=chem_name, chem_dist=fp.dist_b, perm=19999) 
distance_substrate(data=max_non_integ_max_total_microcosm, meth="pearson", chem_n=chem_name, chem_dist=fp.dist_c, perm=19999) 

#Combine all ecoplates together
total_xitou_microcosm<-append(max_integ_20142015_Xitou_Eco, max_non_integ_max_total_microcosm)
distance_substrate(data=total_xitou_microcosm, meth="pearson", chem_n=chem_name, chem_dist=fp.dist_c, perm=19999) 

#For trees drawing only
distance_substrate(data=max_integ_20142015_Xitou_Eco, meth="pearson", chem_n=chem_name, chem_dist=fp.dist_a, perm=99) 
distance_substrate(data=max_non_integ_max_total_microcosm, meth="pearson", chem_n=chem_name, chem_dist=fp.dist_a, perm=99) 
total_xitou_microcosm<-append(max_integ_20142015_Xitou_Eco, max_non_integ_max_total_microcosm)
distance_substrate(data=total_xitou_microcosm, meth="pearson", chem_n=chem_name, chem_dist=fp.dist_c, perm=99) 


##################End of chemical dissimilarity matrix and insitu color development dissimilarity matrix######################



#####################iii)Sub_analysis to investigate the effect of integration period on the explanation power of the statistical model for specific dataset (the best datasets) for Figure S4#########################


#########Preparation of data
#1) Calculating integration (& averaging) of data with different periods for 1-30 days
#2) In total, 24 samples exist for Xitou data, we need to stock 24 x 30 different integrated data into list 
#3) List format is basically one-dimensional vector (list[[i]]), we prepare the index matrix with the size 24x30.
{
  #How long period should be prepared for the whole sampling?
  sample_no=24
  for(i in 1:sample_no) print(length(data_20142015_Xitou_Eco[[i]]))
  #-->30 day
  max_period=30
  
  #preparing list and index matrix
  integ_20142015_Xitou_Eco_sub<-list()
  i_index<-matrix(0, nrow=sample_no, ncol=max_period)  #the index matrix
  for(i in 1:sample_no) for(j in 1:max_period) i_index[i,j] <-(i-1)*max_period + j #assign integer to the element of i_index
  
  #Calcualting average (=integration/period) for each intergration period
  for(j in 1:max_period) {
    for(i in 1:sample_no) {
      integ_20142015_Xitou_Eco_sub[[i_index[i,j]]] <-integ_ecoplate(data_20142015_Xitou_Eco[[i]], period=j)
    }
  }
  
  #test
  integ_20142015_Xitou_Eco_sub[[i_index[24,1]]]
  integ_20142015_Xitou_Eco_sub[[i_index[24,25]]]
  
  #Treatments on triplicates
  #for Taking max of the intergrated values
  #j is for different integration period (1-30)
  #We choose this setting because this is the one of the best explanation power R2 (Fig.6)
  max_integ_20142015_Xitou_Eco_sub<-list()
  for(j in 1:max_period) {
    for(i in 1:sample_no) {
      max_integ_20142015_Xitou_Eco_sub[[i_index[i,j]]] <-max_ecoplate(integ_20142015_Xitou_Eco_sub[[i_index[i,j]]])
    }
  }
  
  #test, each element is vector
  is.vector(max_integ_20142015_Xitou_Eco_sub[1])
  max_integ_20142015_Xitou_Eco_sub[[i_index[24,25]]]
  max_integ_20142015_Xitou_Eco_sub[[i_index[24,30]]]


#This is the extended function for the function above ("Chem_Disim_MF_xitou") to investigate the impact of integration period on the statistical power (R2)

#The parameter, data, should be a list
#Return type is the vector  
Chem_Disim_MF_xitou_for_integ_period <- function(data=max_integ_20142015_Xitou_Eco_sub, thres=0.9, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, period = 10)
{
  #preparation of the indices
  max_period=30
  sample_no=24
  i_index<-matrix(0, nrow=sample_no, ncol=max_period)  #prepare an empty matrix with value 0 for all elements
  for(i in 1:sample_no) for(j in 1:max_period) i_index[i,j] <-(i-1)*max_period + j
  
  #First compile all of the 24 day's data
  file_summary <- data[[i_index[1, period]]]
  for(i in 2:sample_no) file_summary <- rbind(file_summary, data[[i_index[i, period]]])
  file_summary <-file_summary[,c(-1)] #exclusing water data (control, zero)
  row.names(file_summary)<-format_data$sample
  colnames(file_summary)<-chem_n
  
  #Convert data into binary values by quantile (quantile-based multifunctionality)
  mf<-list()
  for(i in 1:31) mf[[i]] <-(file_summary[,i] > quantile(file_summary[,i], thres))
  ts<-as.data.frame(mf)
  row.names(ts)<-format_data$sample
  colnames(ts)<-chem_n
  ts[ts==TRUE]<-1  #converting T,F to 1 and 0
  
  #Calculating tree based on chemical similarity matrix
  clus.hier<-hclust(chem_matrix, method="average")
  tree<-as.phylo(clus.hier)  #convert to phytogenetic tree
  
  #Using function in picante (phylogenetic community analysis)
  pd.result <- pd(ts, tree, include.root=T)
  colnames(pd.result)<-c("Chem_MF", "MF")
  pd.result<-cbind(format_data, pd.result)
  
  #In fact, MF is chemical simmilarity independent multifuncitonality (0-31)
  model_nw<-lm(MF~day*treatment, pd.result)
  #PD is chemical-similarity-weighter multifuncitonality 
  model_w<-lm(Chem_MF~day*treatment, pd.result)
  
  #To calculate R2 values from ANOVA results
  anv_nw<-anova(model_nw)
  R2_nw=(anv_nw$`Sum Sq`[1]+anv_nw$`Sum Sq`[2]+anv_nw$`Sum Sq`[3])/(anv_nw$`Sum Sq`[1]+anv_nw$`Sum Sq`[2]+anv_nw$`Sum Sq`[3]+anv_nw$`Sum Sq`[4])
  
  anv_w<-anova(model_w)
  R2_w=(anv_w$`Sum Sq`[1]+anv_w$`Sum Sq`[2]+anv_w$`Sum Sq`[3])/(anv_w$`Sum Sq`[1]+anv_w$`Sum Sq`[2]+anv_w$`Sum Sq`[3]+anv_w$`Sum Sq`[4])
  
  x<-1:2 #just prepare output vector with dummy value
  
  x[1] = R2_nw  #chemical-information-Non-Weighted result
  x[2] = R2_w   #chemical-information-weighted result
  
  return(x)
  
  
}

#This is the extended function for the function above ("w_s_adonis_xitou") to investigate the impact of integration period on the statistical power (R2)
w_s_adonis_xitou_for_integ_period <- function(data=max_integ_20142015_Xitou_Eco_sub, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, meth="bray", period=10) {
  
  
  #preparation of the indices
  max_period=30
  sample_no=24
  i_index<-matrix(0, nrow=sample_no, ncol=max_period)
  for(i in 1:sample_no) for(j in 1:max_period) i_index[i,j] <-(i-1)*max_period + j
  
  #First compile all of the 24 day's data
  file_summary <- data[[i_index[1, period]]]
  for(i in 2:sample_no) file_summary <- rbind(file_summary, data[[i_index[i, period]]])
  file_summary <-file_summary[,c(-1)] #exclusing water data (control, zero)
  file_summary[file_summary <0]<-0 #convert negative values to zero
  file_summary <- as.data.frame(file_summary) #converting to dataframe
  row.names(file_summary)<-format_data$sample
  colnames(file_summary)<-chem_n
  
  #For non-weighted (by color depth) analysis
  mf0.9<-list()  #for threshold 0.9
  for(i in 1:31) mf0.9[[i]] <-(file_summary[,i] > quantile(file_summary[,i], 0.9))
  mf0.9_summary<-as.data.frame(mf0.9)
  row.names(mf0.9_summary)<-format_data$sample
  colnames(mf0.9_summary)<-chem_n
  mf0.9_summary[mf0.9_summary==TRUE]<-1  #converting T,F to 1 and 0
  
  mf0.7<-list() #for threshold 0.7
  for(i in 1:31) mf0.7[[i]] <-(file_summary[,i] > quantile(file_summary[,i], 0.7))
  mf0.7_summary<-as.data.frame(mf0.7)
  row.names(mf0.7_summary)<-format_data$sample
  colnames(mf0.7_summary)<-chem_n
  mf0.7_summary[mf0.7_summary==TRUE]<-1  #converting T,F to 1 and 0
  
  mf0.5<-list() #for threshold 0.5
  for(i in 1:31) mf0.5[[i]] <-(file_summary[,i] > quantile(file_summary[,i], 0.5))
  mf0.5_summary<-as.data.frame(mf0.5)
  row.names(mf0.5_summary)<-format_data$sample
  colnames(mf0.5_summary)<-chem_n
  mf0.5_summary[mf0.5_summary==TRUE]<-1  #converting T,F to 1 and 0mf0.5<-list()

  #Calculating tree based on chemical similarity matrix
  clus.hier<-hclust(chem_matrix, method="average")
  tree<-as.phylo(clus.hier)  #convert to phytogenetic tree
  
  
  #Calculating through binarized data by using {picante:unifrac}
  dist_binary_cw0.9<-unifrac(mf0.9_summary, tree)
  dist_binary_cw0.7<-unifrac(mf0.7_summary, tree)
  dist_binary_cw0.5<-unifrac(mf0.5_summary, tree)
  #error control
  dist_binary_cw0.9[is.nan(dist_binary_cw0.9)]<-0
  dist_binary_cw0.7[is.nan(dist_binary_cw0.7)]<-0
  dist_binary_cw0.5[is.nan(dist_binary_cw0.5)]<-0
  
  #Calculating from continuous value using GUniFrac
  unif_w<-GUniFrac(file_summary, tree)
  dist_cont_GuniFrac<-as.dist(unif_w$unifracs[,,"d_1"])  #Weighted distance
  dist_cont_GuniFrac0.5<-as.dist(unif_w$unifracs[,,"d_0.5"])  #Weighted distance
  
  #Calculating weighted dissimilarity using fuzzy function
  weighting_tree<-belonging(chem_matrix, standardize = TRUE)
  
  #Since the order of substrates is different between file_summary and weighting_tree, this can be fitted by sorting
  file_summary_sorted<-file_summary[,sort(colnames(file_summary))]
  weighting_tree_sorted<-weighting_tree[sort(rownames(weighting_tree)),sort(colnames(weighting_tree))]
  
  #Fuzzy weighted
  file_summary_cwF <- as.matrix(file_summary_sorted)%*%as.matrix(weighting_tree_sorted)
  
  dist_cont_cwF<-vegdist(file_summary_cwF, method=meth)
  #Calculating continuous dissimilarity 
  dist_cont_uw<-vegdist(file_summary, method=meth)
  
  #Calculating default dissimilarity based on presence/absence 
  dist_binary_uw0.9<-vegdist(mf0.9_summary, method="euclidian")
  dist_binary_uw0.7<-vegdist(mf0.7_summary, method="euclidian")
  dist_binary_uw0.5<-vegdist(mf0.5_summary, method="euclidian")
  
  
  #PERMANOVA, without permutation (just set as 2)
  diff_binary_uw0.9<-adonis(dist_binary_uw0.9 ~ format_data$day*format_data$treatment, permutations=2)
  diff_binary_uw0.7<-adonis(dist_binary_uw0.7 ~ format_data$day*format_data$treatment, permutations=2)
  diff_binary_uw0.5<-adonis(dist_binary_uw0.5 ~ format_data$day*format_data$treatment, permutations=2)
  
  diff_binary_cw0.9<-adonis(dist_binary_cw0.9 ~ format_data$day*format_data$treatment, permutations=2)
  diff_binary_cw0.7<-adonis(dist_binary_cw0.7 ~ format_data$day*format_data$treatment, permutations=2)
  diff_binary_cw0.5<-adonis(dist_binary_cw0.5 ~ format_data$day*format_data$treatment, permutations=2)
  
  diff_cont_uw<-adonis(dist_cont_uw ~ format_data$day*format_data$treatment, permutations=2)
  diff_cont_cwF<-adonis(dist_cont_cwF ~ format_data$day*format_data$treatment, permutations=2)
  
  diff_cont_GuniFrac<-adonis(dist_cont_GuniFrac ~ format_data$day*format_data$treatment, permutations=2)
  diff_cont_GuniFrac0.5<-adonis(dist_cont_GuniFrac0.5 ~ format_data$day*format_data$treatment, permutations=2)
  
  x<-1:10 #just prepare output vector with dummy value
  
  x[1] = 1-diff_binary_uw0.9$aov.tab$R2[4]
  x[2] = 1-diff_binary_uw0.7$aov.tab$R2[4]
  x[3] = 1-diff_binary_uw0.5$aov.tab$R2[4]
  x[4] = 1-diff_binary_cw0.9$aov.tab$R2[4]
  x[5] = 1-diff_binary_cw0.7$aov.tab$R2[4]
  x[6] = 1-diff_binary_cw0.5$aov.tab$R2[4]
  x[7] = 1-diff_cont_uw$aov.tab$R2[4]
  x[8] = 1-diff_cont_cwF$aov.tab$R2[4]
  x[9] = 1-diff_cont_GuniFrac$aov.tab$R2[4]
  x[10] = 1-diff_cont_GuniFrac0.5$aov.tab$R2[4]
  
  return(x)
}

#This is the extended function for the function above ("W_s_dbRDA") to investigate the impact of integration period on the statistical power (R2)
W_s_dbRDA_for_integ_period <- function(data=max_integ_20142015_Xitou_Eco_sub, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, meth="bray", thres=0.9, period=10) {
  
  #preparation of the indices
  max_period=30
  sample_no=24
  i_index<-matrix(0, nrow=sample_no, ncol=max_period)
  for(i in 1:sample_no) for(j in 1:max_period) i_index[i,j] <-(i-1)*max_period + j
  
  #First compile all of the 24 day's data
  file_summary <- data[[i_index[1, period]]]
  for(i in 2:sample_no) file_summary <- rbind(file_summary, data[[i_index[i, period]]])
  file_summary <-file_summary[,c(-1)] #exclusing water data (control, zero)
  file_summary[file_summary <0]<-0 #convert negative values to zero
  file_summary <- as.data.frame(file_summary) #converting to dataframe
  row.names(file_summary)<-format_data$sample
  colnames(file_summary)<-chem_n
  
  #For binarized (by color depth) analysis
  mf<-list()
  for(i in 1:31) mf[[i]] <-(file_summary[,i] > quantile(file_summary[,i], thres))
  mf_summary<-as.data.frame(mf)
  row.names(mf_summary)<-format_data$sample
  colnames(mf_summary)<-chem_n
  mf_summary[mf_summary==TRUE]<-1  #converting T,F to 1 and 0
  
  #Calculating tree based on chemical similarity matrix
  clus.hier<-hclust(chem_matrix, method="average")
  tree<-as.phylo(clus.hier)  #convert to phytogenetic tree
  
  #Calculating weighted similarity on binarized databy using {picante:unifrac}
  binary_cw_dist<-unifrac(mf_summary, tree)
  #error control
  binary_cw_dist[is.nan(binary_cw_dist)]<-0
  
  #Calculating chemically weighted dissimilarity on continous data using GUniFrac
  unif_cont_cw<-GUniFrac(file_summary, tree)
  cont_dist_cw<-as.dist(unif_cont_cw$unifracs[,,"d_1"])  #Weighted distance
  cont_dist_cw0.5<-as.dist(unif_cont_cw$unifracs[,,"d_0.5"])  #Weighted distance
  
  #Calculating weighted dissimilarity using fuzzy function
  weighting_tree<-belonging(chem_matrix, standardize = TRUE)
  
  #Since the order of substrates is different between file_summary and weighting_tree, this can be fitted by sorting
  file_summary_sorted<-file_summary[,sort(colnames(file_summary))]
  weighting_tree_sorted<-weighting_tree[sort(rownames(weighting_tree)),sort(colnames(weighting_tree))]
  
  #Fuzzy weighted
  file_summary_cwF <- as.matrix(file_summary_sorted)%*%as.matrix(weighting_tree_sorted)
  cont_dist_cwF<-vegdist(file_summary_cwF, method=meth)
  
  #Calculating dissimilarity on continuous data, w/o chemical weight 
  cont_uw_dist<-vegdist(file_summary, method=meth)
  
  #Calculating default dissimilarity based on presence/absence, w/o chemical weight 
  binary_uw_dist<-vegdist(mf_summary, method="euclidian")
  
  #Distance-based RDA (db-RDA) by vegan
  binary_cw_effect<-capscale(binary_cw_dist~format_data$day*format_data$treatment)
  cont_cw_effect<-capscale(cont_dist_cw~format_data$day*format_data$treatment)
  cont_cw0.5_effect<-capscale(cont_dist_cw0.5~format_data$day*format_data$treatment)
  cont_cwF_effect<-capscale(cont_dist_cwF~format_data$day*format_data$treatment)
  cont_uw_effect<-capscale(cont_uw_dist~format_data$day*format_data$treatment)
  binary_uw_effect<-capscale(binary_uw_dist~format_data$day*format_data$treatment)
  
  #Extract the constrained variances
  real_total_chi_binary_cw= binary_cw_effect$CCA$tot.chi + binary_cw_effect$CA$tot.chi
  real_total_chi_cont_cw= cont_cw_effect$CCA$tot.chi + cont_cw_effect$CA$tot.chi
  real_total_chi_cont_cw0.5= cont_cw0.5_effect$CCA$tot.chi + cont_cw0.5_effect$CA$tot.chi
  real_total_chi_cont_cwF= cont_cwF_effect$CCA$tot.chi + cont_cwF_effect$CA$tot.chi
  real_total_chi_cont_uw= cont_uw_effect$CCA$tot.chi + cont_uw_effect$CA$tot.chi
  real_total_chi_binary_uw= binary_uw_effect$CCA$tot.chi + binary_uw_effect$CA$tot.chi
  
  x<-1:6
  x[1] <- binary_uw_effect$CCA$tot.chi/real_total_chi_binary_uw
  x[2] <- binary_cw_effect$CCA$tot.chi/real_total_chi_binary_cw
  x[3] <- cont_cw_effect$CCA$tot.chi/real_total_chi_cont_cw
  x[4] <- cont_cw0.5_effect$CCA$tot.chi/real_total_chi_cont_cw0.5
  x[5] <- cont_cwF_effect$CCA$tot.chi/real_total_chi_cont_cwF
  x[6] <- cont_uw_effect$CCA$tot.chi/real_total_chi_cont_uw
  
  return(x)
}




#For checking R2 by linear model for each measure
#test
Chem_Disim_MF_xitou_for_integ_period(data=max_integ_20142015_Xitou_Eco_sub, thres=0.9, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, period = 10)
#matrix for output with different threshold value
Var_lm0.9<-matrix(0, nrow=2, ncol=max_period)
Var_lm0.7<-matrix(0, nrow=2, ncol=max_period)
Var_lm0.5<-matrix(0, nrow=2, ncol=max_period)

#calculate R2 value from the integration period=1 to 30 with the threshold=0.9,0.7,0.5
for(j in 1:30) {
  x<-Chem_Disim_MF_xitou_for_integ_period(data=max_integ_20142015_Xitou_Eco_sub, thres=0.9, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, period=j)
  for(i in 1: 2) Var_lm0.9[i,j] <-x[i]
  x<-Chem_Disim_MF_xitou_for_integ_period(data=max_integ_20142015_Xitou_Eco_sub, thres=0.7, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, period=j)
  for(i in 1: 2) Var_lm0.7[i,j] <-x[i]
  x<-Chem_Disim_MF_xitou_for_integ_period(data=max_integ_20142015_Xitou_Eco_sub, thres=0.5, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, period=j)
  for(i in 1: 2) Var_lm0.5[i,j] <-x[i]
}

#Plot the results
plot(0,0, type="n", xlim=c(0,30), ylim=c(0,1), xlab="integration period", ylab="R2")  #draw the frame of the graph
for(i in 1:2) {
  points(1:30, Var_lm0.9[i,], pch=i, col=i)
  lines(1:30, Var_lm0.9[i,], lty=i, col=i)
  points(1:30, Var_lm0.7[i,], pch=i, col=i+2)
  lines(1:30, Var_lm0.7[i,], lty=i, col=i+2)
  points(1:30, Var_lm0.5[i,], pch=i, col=i+4)
  lines(1:30, Var_lm0.5[i,], lty=i, col=i+4)
}
labels <- c("binary_nw0.9", "binary_cw0.9", "binary_nw0.7", "binary_cw0.7", "binary_nw0.5", "binary_cw0.5")
legend("bottomright", legend=labels, col=1:6, pch=1:6, lty=1:6)
title("linear model for binarized multifunctionality")

#For checking R2 by PERMANOVA for each measure
#test
w_s_adonis_xitou_for_integ_period(data=max_integ_20142015_Xitou_Eco_sub, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, meth="bray", period=15)

#matrix for output
Var_adonis<-matrix(0, nrow=10, ncol=max_period)

for(j in 1:30) {
  x<-w_s_adonis_xitou_for_integ_period(data=max_integ_20142015_Xitou_Eco_sub, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, meth="bray", period=j)
  for(i in 1: 10) Var_adonis[i,j] <-x[i]
}

#Plot the results
plot(0,0, type="n", xlim=c(0,30), ylim=c(0,1), xlab="integration period", ylab="R2")
for(i in 1:10) {
  points(1:30, Var_adonis[i,], pch=i, col=i)
  lines(1:30, Var_adonis[i,], lty=i, col=i)
}
labels <- c("binary_uw0.9", "binary_uw0.7", "binary_uw0.5", "binary_cw0.9", "binary_cw0.7", "binary_cw0.5", "cont_uw", "cont_cwF", "cont_cwGuniFrac", "cont_cwGuniFrac0.5")
legend("bottomright", legend=labels, col=1:10, pch=1:10, lty=1:10)
title("PERMANOVA")

#For checking contrained variation by db-RDA for each measure, which is used for Fig.S6
#test
W_s_dbRDA_for_integ_period(data=max_integ_20142015_Xitou_Eco_sub, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, meth="bray", thres=0.9, period=4)

Var_dbRDA<-matrix(0, nrow=6, ncol=max_period)  #Prepare the matrix to stock results; each result for each integration period is a vector with 6 values

for(j in 1:30) {
  x<-W_s_dbRDA_for_integ_period(data=max_integ_20142015_Xitou_Eco_sub, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, meth="bray", thres=0.8, period=j)
  for(i in 1: 6) Var_dbRDA[i,j] <-x[i]
}

plot(0,0, type="n", xlim=c(0,30), ylim=c(0,1), xlab="integration period", ylab="R2")
for(i in 1:6) {
  points(1:30, Var_dbRDA[i,], pch=i, col=i)
  lines(1:30, Var_dbRDA[i,], lty=i, col=i)
}
labels <- c("binary_uw", "binary_cw",  "cont_cwGuniFrac", "cont_cwGuniFrac0.5", "cont_cwF", "cont_uw")
legend("bottomright", legend=labels, col=1:6, pch=1:6, lty=1:6)
title("RDA")


}


############################iv) TEST COMPARE PERMANOVA and db-RDA##########################################

w_s_adonis_xitou_test <- function(data=max_integ_20142015_Xitou_Eco, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, meth="bray", thres=0.9, perm=199) {
  
  ####IMPORTANT: The comments within this function are simplified.Check every comment on the function above "Chem_Disim_MF_xitou"##########   
  #First compile all of the 24 day's data
  file_summary <- data[[1]]
  for(i in 2:length(data)) file_summary <- rbind(file_summary, data[[i]])
  file_summary <-file_summary[,c(-1)] #excluding water data (control, zero)
  file_summary[file_summary <0]<-0 #convert negative values to zero, which is necessary to calculate the bray-curtis distance
  file_summary <- as.data.frame(file_summary) #converting to dataframe
  row.names(file_summary)<-format_data$sample #set names
  colnames(file_summary)<-chem_n              #set names   
  
  #For the analysis after conversion of color depth to presence/absence, by quantile-based binarization (use three thresholds)
  mf_given<-list()  #for given threshold 
  for(i in 1:31) mf_given[[i]] <-(file_summary[,i] > quantile(file_summary[,i], thres))
  mf_summary<-as.data.frame(mf_given)
  row.names(mf_summary)<-format_data$sample
  colnames(mf_summary)<-chem_n
  mf_summary[mf_summary==TRUE]<-1  #converting T,F to 1 and 0

  #Calculating default dissimilarity based on presence/absence (without chemical-similarity information) [baseline of the analysis]
  dist_binary_uw<-vegdist(mf_summary, method="euclidian")

  #######Calculating the default continuous dissimilarity (without chemical-similarity information) [baseline of the analysis]
  dist_cont_uw<-vegdist(file_summary, method=meth)  
  

  
  ###########PERMANOVA for each of the different indices of color patten dissimilarity 
  #binarized, without chemical-similarity 
  diff_binary_uw<-adonis(dist_binary_uw ~ format_data$day*format_data$treatment, permutations=perm)

  #continuous, without chemical-similarity
  diff_cont_uw<-adonis(dist_cont_uw ~ format_data$day*format_data$treatment, permutations=perm)

  #Output1, P values 
  #cat("The P values of data, treatment and their interaction (binarized w/o CW):", diff_binary_uw$aov.tab$Pr, "\n")
  #cat("The P values of data, treatment and their interaction (continuous w/o CW):", diff_cont_uw$aov.tab$Pr, "\n")

  #Output2, R2 values
  #cat("The R2 of data, treatment and their interaction (binarized w/o CW):", 1-diff_binary_uw$aov.tab$R2[4], "\n")
  #cat("The R2 of data, treatment and their interaction (continuous w/o CW):", 1-diff_cont_uw$aov.tab$R2[4], "\n")

  #Display the summary results
  #diff_binary_uw
  diff_cont_uw
  
}

W_s_dbRDA_test <- function(data=max_integ_20142015_Xitou_Eco, format_data=format_d, chem_n=chem_name, chem_matrix=fp.dist_c, meth="bray", thres=0.9, perm=199) {
  
  ####IMPORTANT: The comments within this function are simplified.Check every comment on the function above "Chem_Disim_MF_xitou"########## 
  #First compile all of the 24 day's data
  file_summary <- data[[1]]
  for(i in 2:length(data)) file_summary <- rbind(file_summary, data[[i]])
  file_summary <-file_summary[,c(-1)] #exclusing water data (control, zero)
  file_summary[file_summary <0]<-0 #convert negative values to zero
  file_summary <- as.data.frame(file_summary) #converting to dataframe
  row.names(file_summary)<-format_data$sample
  colnames(file_summary)<-chem_n
  
  #For binarized (by color depth) analysis
  mf<-list()
  for(i in 1:31) mf[[i]] <-(file_summary[,i] > quantile(file_summary[,i], thres))
  mf_summary<-as.data.frame(mf)
  row.names(mf_summary)<-format_data$sample
  colnames(mf_summary)<-chem_n
  mf_summary[mf_summary==TRUE]<-1  #converting T,F to 1 and 0
  
  
  
  #Calculating dissimilarity on continuous data, w/o chemical weight 
  cont_uw_dist<-vegdist(file_summary, method=meth)
  
  #Calculating default dissimilarity based on presence/absence, w/o chemical weight 
  binary_uw_dist<-vegdist(mf_summary, method="euclidian")
  
  #Distance-based RDA (db-RDA) by vegan
  cont_uw_effect<-capscale(cont_uw_dist~format_data$day*format_data$treatment)
  binary_uw_effect<-capscale(binary_uw_dist~format_data$day*format_data$treatment)
  
  #Permutation test
  sig_cont_uw<-anova(cont_uw_effect, permutations=perm, by='terms')
  sig_binary_uw<-anova(binary_uw_effect, permutations=perm, by='terms')
  
  
  #Extract the constrained variances
  real_total_chi_cont_uw= cont_uw_effect$CCA$tot.chi + cont_uw_effect$CA$tot.chi
  real_total_chi_binary_uw= binary_uw_effect$CCA$tot.chi + binary_uw_effect$CA$tot.chi
  
  #OUTPUT 1: the explanation power of the mode (by the contrained variance)
  #cat("The constrained variance by UW binary data by treatment and their interaction:", binary_uw_effect$CCA$tot.chi/real_total_chi_binary_uw, "\n")
  #cat("The constrained variance by UW continuous data by treatment and their interaction:", cont_uw_effect$CCA$tot.chi/real_total_chi_cont_uw, "\n")
  
  
  #OUTPUT2: The P values
  #cat("The P values for day, treatment and ther interactions by UW binary data are", sig_binary_uw$`Pr(>F)`,"\n")
  #cat("The P values for day, treatment and ther interactions by UW continuous data are", sig_cont_uw$`Pr(>F)`,"\n")
  
  
  #Display summary results
  #binary_uw_effect
  cont_uw_effect
}


w_s_adonis_xitou_test(data=ave_integ_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_a, thres=0.9, perm=999)
W_s_dbRDA_test(data=ave_integ_20142015_Xitou_Eco, format_data=format_d,chem_matrix=fp.dist_a, thres=0.9, perm=999)


##############End of PART4############################################################################

