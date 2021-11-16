


################################################################################################## 
#          Project: R script for selecting potential halogenated persistent organic pollutants  
#          from high resolution mass spectrometry data. The script is described in:
#
#          "Compositional space: A guide for environmental chemists on the identification 
#          of persistent and bioaccumulative organics using mass spectrometry" by Xianming Zhang, 
#          Robert A. Di Lorenzo, Paul A. Helm, Eric J. Reiner, Philip H. Howard, Derek C. G. Muir, 
#	   John G. Sled, Karl J. Jobst (https://doi.org/10.1016/j.envint.2019.05.002)
#
# Program name: HaloSelect.R
#
# Program developed by: Xianming Zhang, Robert A. Di Lorenzo and Karl J. Jobst
#
# Program developed on an Apple-Darwin15.6.0 platform with R version 3.5.0; R packages: openxlsx 
# version 4.0.17;plotly version 4.7.1 
#
#
# Purpose: To select halogenated organics from (LC or GC) high resolution mass spectrometry data 
#	    on the basis of isotope peak clustering and ratios.
#         
#
# Revision History  :
#==============================================================
# Date (YYYYMMDD)   Author   Version   Revised Content
#--------------------------------------------------------------
# 20180725          XZ       1         Program created 
#
#
#==============================================================
# Instructions:
#
# Following installation of R and requisite packages (openxlsx, plotly):
#
# (1) Create a folder containing the input file, which can either be in the .xlsx or .csv format.
#     The input file should contain colunms named 'id','mz','rt',and 'intb', which correspond to 
#     identity numbers of each peak, accurate mass (m/z), retention time, and intensity.
#
# (2) Save this script to the same folder and specify the folder location, input filename 
#     (e.g.filename <- "2585.xlsx"), and unit of retention time (e.g. rt_unit <- "s" or "min").
#
# (3) Run the script.
#
# (4) After the script is finished running, there will be a "data.csv" file containing the 
#     isotopic cluster information, i.e., annotation of A, A+1(Ap1), A+2(Ap2) and A-2(Am2) peaks.
#     The "output.csv" file reports the subset of peaks that satisfy the criteria for the 
#     presence of Cl/Br. The user may also design their own compositional space filter and apply 
#     it to the "data.csv" file.
#
#################################################################################################

#Import data from .csv or .xlsx file

if(length(grep('.xlsx',filename))==1){
  library(openxlsx)
  data<-read.xlsx(filename,
                  sheet=1,
                  startRow=1, #cols=1:10,
                  colNames=TRUE,skipEmptyCols=TRUE)
  if (sum(colnames(data)=='id')!=1 | sum(colnames(data)=='rt')!=1 | sum(colnames(data)=='mz')!=1 | sum(colnames(data)=='intb')!=1) {
    stop('Check the columns (and names) of the input file ')
  }
} else if (length(grep('.csv',filename))==1){
  data_temp<-read.csv(filename)
  if (sum(colnames(data)=='id')!=1 | sum(colnames(data)=='rt')!=1 | sum(colnames(data)=='mz')!=1 | sum(colnames(data)=='intb')!=1) {
    stop('Check the columns (and names) of the input file ')
  }  
} else {
  stop('Not the right type of input file')
}

data<-data[,c('id','mz','rt','intb')]

if (rt_unit=='s'){
  data[,'rt']<-data[,'rt']/60 # change rt from second to minute
  #data[,c('rt', 'rtmin','rtmax')]<-data[,c('rt', 'rtmin','rtmax')]/60 # change rt from second to minute
}

data<-data[data$rt>1,]

data$cluster <- NA
data<-data[order(data$rt),]
m<-1

for (n in 1:dim(data)[1]) {
  if (is.na(data$cluster[n])){
    #critera<-abs(data$rt[n]-data$rt)<0.03 & abs(data$rtmin[n]-data$rtmin)<0.03 & abs(data$rtmax[n]-data$rtmax)<0.03
    #critera<-abs(data$rt[n]-data$rt)<0.02 & abs(data$rtmax[n]-data$rtmin[n]-(data$rtmax-data$rtmin))<0.03 
    critera<-abs(data$rt[n]-data$rt)<0.03  
    data$cluster[critera]<-m
    m<-m+1
  }
}

data<-data[order(data$cluster),]

# Masses and tolerances for isotopic clustering

p1H<-1.00783
p1C <- 1.00335
p2ClBr <- 1.997
p1Si <- 0.99956
pCpHwindow<-0.00224
window <- 0.005

data$Mp1H <- NA
data$Mp1C <- NA
data$Mp1Si <- NA
data$Mp2ClBr <- NA

# Create dataframe for isotopic clusters

ClBrSiCluster<- data.frame(M= rep(NA,dim(data)[1]),Mp2= NA,Mp4= NA,Mp6= NA,Mp8= NA,Mp10= NA,Mp12= NA,Mp14= NA)
intensity<- data.frame(intensityM= rep(0,dim(data)[1]),intensityMp2= 0,intensityMp4= 0,intensityMp6= 0,intensityMp8= 0,intensityMp10= 0,intensityMp12= 0,intensityMp14= 0)

data$Ap2 <- 0
data$Am2 <- 0
data$ApSi <- NA
data$intensityApSi <- 0
data$ApC <- NA
data$intensityApC <- 0
data$ApH <- NA
data$intensityApH <- 0

flag <- data.frame("ChemNo"= 1:dim(data)[1],
                   "PeakMpH"=FALSE,
                   "isotopePeakC"=FALSE,
                   "isotopePeakSi"=FALSE,
                   "isotopePeakClBr"=FALSE,
                   "Ap2GTApSi"=FALSE,
                   "Ap2GTApC"=FALSE,
                   "ClBrMonoIsotope"=FALSE,
                   "ClBrNonMonoIsotopes"=FALSE)


#Flag isotope peaks ----------->
pb <- txtProgressBar(min = 0, max = dim(data)[1], style = 3)
print('find MpH and C, Si, Cl, Br isotope peaks')

for (n in 1:dim(data)[1]) {
  isMp1H <- (data$cluster - data$cluster[n] == 0) & (data$mz-data$mz[n] < p1H+pCpHwindow) & (data$mz-data$mz[n] > p1H-pCpHwindow) & (data$mz-data$mz[n] != 0)
  tryCatch(data$Mp1H[n] <- data$mz[isMp1H], error = function(e) NULL) #XZ note: in the mass list, there can be multiple mz within the criteria above and result in warning message "number of items to replace is not a multiple of replacement length
  flag$PeakMpH[which(isMp1H)] <- TRUE
  
  isMp1C <- (data$cluster - data$cluster[n] == 0) & (data$mz-data$mz[n] < p1C+pCpHwindow) & (data$mz-data$mz[n] > p1C-pCpHwindow) & (data$mz-data$mz[n] != 0)
  tryCatch(data$Mp1C[n] <- data$mz[isMp1C], error = function(e) NULL) #XZ note: in the mass list, there can be multiple mz within the criteria above and result in warning message "number of items to replace is not a multiple of replacement length
  flag$isotopePeakC[which(isMp1C)] <- TRUE
  
  isMp1Si <- (data$cluster - data$cluster[n] == 0) & (data$mz-data$mz[n] < p1Si+window) & (data$mz-data$mz[n] > p1Si-window) & (data$mz-data$mz[n] != 0)
  tryCatch(data$Mp1Si[n] <- data$mz[isMp1Si], error = function(e) NULL)
  flag$isotopePeakSi[which(isMp1Si)] <- TRUE
  
  isMp2ClBr <- (data$cluster - data$cluster[n] == 0) & (data$mz-data$mz[n] < p2ClBr+window) & (data$mz-data$mz[n] > p2ClBr-window) & (data$mz-data$mz[n] != 0)
  tryCatch(data$Mp2ClBr[n] <- data$mz[isMp2ClBr], error = function(e) NULL)
  flag$isotopePeakClBr[which(isMp2ClBr)] <- TRUE
  
  # update progress bar
  setTxtProgressBar(pb, n)
}

flag$isotopePeak <- flag$PeakMpH |flag$isotopePeakC | flag$isotopePeakSi | flag$isotopePeakClBr
ClBrSiCluster$M[!flag$isotopePeak] <- data$mz[!flag$isotopePeak]
intensity$intensityM[!flag$isotopePeak] <- data$intb[!flag$isotopePeak]
#<-----------Flag isotope peaks 

#Filter Cl/Br compounds ---------->
pb <- txtProgressBar(min = 0, max = dim(data)[1], style = 3)#initialize progress bar
print('Filter Cl/Br Compounds')

for (n in 1:dim(data)[1]) {
  if (!is.na(ClBrSiCluster$M[n])) {
    
    for (m in 2:8){
      isMpx <- (data$cluster - data$cluster[n] == 0) & (data$mz-ClBrSiCluster[n,m-1] < p2ClBr+window) & (data$mz-ClBrSiCluster[n,m-1] > p2ClBr-window) 
      tryCatch(ClBrSiCluster[n,m] <- data$mz[which(isMpx)], error = function(e) NULL)
      tryCatch(intensity[n,m] <- data$intb[which(isMpx)], error = function(e) NULL)
    }
    
    maxposition <-which.max(intensity[n,!is.na(intensity[n,])])
    if(maxposition >1){
      if (intensity[n,maxposition-1]/intensity[n,maxposition]>0.95 & intensity[n,maxposition-1]/intensity[n,maxposition]<1.05) {
        maxposition<-maxposition-1 
      }  
    }
    
    tryCatch(data$Ap2[n] <- 100*intensity[n,maxposition+1]/intensity[n,maxposition], error = function(e) NULL)
    tryCatch(data$Am2[n] <- 100*intensity[n,maxposition-1]/intensity[n,maxposition], error = function(e) NULL)
    
    isApSi <- data$cluster - data$cluster[n] == 0 & data$mz-ClBrSiCluster[n,maxposition] < p1Si+window & data$mz-ClBrSiCluster[n,maxposition] >  p1Si-window 
    tryCatch(data$ApSi[n] <- data$mz[which(isApSi)], error = function(e) NULL)
    tryCatch(data$intensityApSi[n] <- data$intb[which(isApSi)], error = function(e) NULL)
    
    isApC <- data$cluster - data$cluster[n] == 0 & data$mz-ClBrSiCluster[n,maxposition] < p1C+0.001 & data$mz-ClBrSiCluster[n,maxposition] > p1C-0.001 
    tryCatch(data$ApC[n] <- data$mz[which(isApC)], error = function(e) NULL)
    tryCatch(data$intensityApC[n] <- data$intb[which(isApC)], error = function(e) NULL)
    
    isApH <- data$cluster - data$cluster[n] == 0 & data$mz-ClBrSiCluster[n,maxposition] < p1H+window & data$mz-ClBrSiCluster[n,maxposition] > p1H-window 
    tryCatch(data$ApH[n] <- data$mz[which(isApH)], error = function(e) NULL)
    tryCatch(data$intensityApH[n] <- data$intb[which(isApH)], error = function(e) NULL)
    
    tryCatch(flag$Ap2GTApSi[n] <- 100*data$intensityApSi[n]/intensity[n,maxposition] < data$Ap2[n], error = function(e) NULL)
    tryCatch(flag$Ap2GTApC[n]<- 100*data$intensityApC[n]/intensity[n,maxposition] < data$Ap2[n], error = function(e) NULL)
    
    flag$ClBrMonoIsotope[n] <- isTRUE(flag$Ap2GTApSi[n]) & isTRUE(flag$Ap2GTApC[n]) &  #495.7552 rt3.58 has A+1C/H,  
      #((data$Am2[n]==0) & data$Ap2[n] >30) | (data$Am2[n]>25 & data$Ap2[n] > 85-data$Am2[n] & data$Ap2[n]>0)&
      #((data$Am2[n]==0 & data$Ap2[n] >25) | (data$Ap2[n] > 134-2.28*data$Am2[n] & data$Ap2[n]>0)) &
      ((data$Am2[n]==0 & data$Ap2[n] >25) | (data$Am2[n]>15 & data$Ap2[n]>0)) &
      #!(data$mz[n]>700 & data$rt[n]<6) &
      (100*data$intensityApSi[n]/apply(intensity[n,], 1, max)<0.0643*data$mz[n]+6.8) #&
    #(100*data$intensityApSi[n]/intensity[n,maxposition]<0.0643*data$mz[n]+6.8)
    #(data$intensityApSi[n]<intensity[n,maxposition+1])
    
    if(flag$ClBrMonoIsotope[n]){
      #print('label other isotopic peaks of ClBr compounds')
      for (m in 2:8){
        isMpx <- (data$cluster - data$cluster[n] == 0) & (data$mz-ClBrSiCluster[n,m-1] < p2ClBr+window) & (data$mz-ClBrSiCluster[n,m-1] > p2ClBr-window) 
        tryCatch(flag$ClBrNonMonoIsotopes[isMpx] <- TRUE, error = function(e) NULL)
      } 
    }
    
  }
  
  # update progress bar
  setTxtProgressBar(pb, n)
}
#<----------Filter Cl/Br compounds

#Output filtered Cl/Br compounds ----------->
output<-cbind(data[flag$ClBrMonoIsotope,'id'],ClBrSiCluster[flag$ClBrMonoIsotope,],intensity[flag$ClBrMonoIsotope,])
colnames(output)[1]<-'id(InOriginalList)'
write.csv(output, file = paste0('filteredClBrCompoundsFrom',gsub('.csv','',gsub('.xlsx','',filename)),'.csv'),row.names=FALSE)
#<----------