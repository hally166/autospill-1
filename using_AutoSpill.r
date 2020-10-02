#This is a little messy and a work in progress.  
#This is me practically using autospill and exporting the result in an fcs file for use on a flow cytometer.

###RUNNING AUTOSPILL
#devtools::install_github("carlosproca/autospill")
devtools::install_github("hally166/autospill-1") #I cant get Carlos's version to work on windows

#library( autospill )
library( autospillCH )
library( flowCore )

# set parameters

asp <- get.autospill.param( "final.step" )

# read flow controls

control.dir <- "C:/Users/chall/Desktop/CompensationFolder"
control.def.file <- "C:/Users/chall/Desktop/CompensationFolder/fcs_control.csv" #for DiVA use the filter, laser, and PMY details


flow.control <- read.flow.control( control.dir, control.def.file, asp )

# gate events before calculating spillover

flow.gate <- gate.flow.data( flow.control, asp )


# get initial spillover matrices from untransformed and transformed data

marker.spillover.unco.untr <- get.marker.spillover( TRUE, flow.gate,
                                                    flow.control, asp )
marker.spillover.unco.tran <- get.marker.spillover( FALSE, flow.gate,
                                                    flow.control, asp )


# refine spillover matrix iteratively

refine.spillover.result <- refine.spillover( marker.spillover.unco.untr,
                                             marker.spillover.unco.tran, flow.gate, flow.control, asp )

###CHECKING THE RESULTS AND EXPORTING THE COMPENSATION TO AN FCS FILE FOR USE ON A SORTER
#I hate scientific numbers
options(scipen = 4)

#function to get the nemas of the parameters
autovect_verbose<- function(ff){
  c<- data.frame(ff@parameters@data)
  d<- grep("FSC|SSC|Time", c$name, invert = TRUE, value = TRUE)
  return(unname(d))
}

#load the files to apply the compensation to
myfiles <- list.files(path="C:/Users/chall/Desktop/Kerstin", pattern=".fcs$")
fs <- flowCore::read.flowSet(myfiles, path="C:/Users/chall/Desktop/Kerstin")

#compensate files - internal
comp <- fsApply(fs, function(x) spillover(x)[[1]], simplify=FALSE) #change the [[3]] dependent on the spillover keyword
fs_comp <- compensate(fs, comp)

#compensate files - autospill - this is all about converting the csv into a compensation matrix for use in flowCore
comp_AS <- read.csv("C:/Users/chall/Desktop/Kerstin/table_compensation/autospill_compensation.csv",check.names = FALSE,header=TRUE)
x<-autovect_verbose(fs[[1]])
colnames(comp_AS)[1]<-"fluor"
x<-c("fluor",x)
comp_AS <- comp_AS[,x]
comp_AS<-comp_AS[order(factor(comp_AS$fluor,levels=x[2:ncol(comp_AS)])),]
comp_AS[sapply(comp_AS, is.numeric)] <- abs(comp_AS[sapply(comp_AS, is.numeric)])
comp_AS <- comp_AS[,2:ncol(comp_AS)]
comp_AS <- do.call(cbind, lapply(comp_AS, as.numeric))
fs_comp_AS <- compensate(fs, comp_AS)

#transform files - this function chooses all the parameters, except FSC SSC and Time, you can add more if you wish
x<-autovect_verbose(fs_comp[[1]])
tf<-estimateLogicle(fs_comp[[1]], channels =  x)
fs_trans<-transform(fs_comp,tf)

x<-autovect_verbose(fs_comp_AS[[1]])
tf<-estimateLogicle(fs_comp_AS[[1]], channels =  x)
fs_trans_AS<-transform(fs_comp_AS,tf)

#export the compensation as a FCS file for import into your favorite sorter software
fFrame <- fs[[8]]
keyword(fFrame)$SPILL<-comp_AS
write.FCS(fFrame,"autospill.fcs")

#plot the results - examples
autoplot(fs_trans[6:11],x="FSC-A",y="SSC-A", bins = 256)
autoplot(fs_trans[6:11],x="670/30 YG-C-A",y="670/30 R-C-A", bins = 256)
autoplot(fs_trans[6:11],x="582/15 YG-E-A",y="670/30 R-C-A", bins = 256)
autoplot(fs_trans[6:11],x="530/30 B-B-A",y="SSC-A", bins = 256)
autoplot(fs_trans[6:11],x="450/50 V-F-A",y="SSC-A", bins = 256)

autoplot(fs_trans_AS,x="FSC-A",y="SSC-A", bins = 256)
autoplot(fs_trans_AS,x="670/30 YG-C-A",y="670/30 R-C-A", bins = 256)
autoplot(fs_trans_AS,x="582/15 YG-E-A",y="670/30 R-C-A", bins = 256)
autoplot(fs_trans_AS,x="SSC-A",y="530/30 B-B-A", bins = 256)
autoplot(fs_trans_AS,x="SSC-A",y="450/50 V-F-A", bins = 256)

autoplot(fs_trans[6:11],x="FSC-A",y="SSC-A", bins = 256)
autoplot(fs_trans_AS[6:11],x="FSC-A",y="SSC-A", bins = 256)
autoplot(fs_trans[6:11],x="670/30 YG-C-A",y="670/30 R-C-A", bins = 256)
autoplot(fs_trans_AS[6:11],x="670/30 YG-C-A",y="670/30 R-C-A", bins = 256)
autoplot(fs_trans[6:11],x="582/15 YG-E-A",y="670/30 R-C-A", bins = 256)
autoplot(fs_trans_AS[6:11],x="582/15 YG-E-A",y="670/30 R-C-A", bins = 256)
autoplot(fs_trans[6:11],x="SSC-A",y="530/30 B-B-A", bins = 256)
autoplot(fs_trans_AS[6:11],x="SSC-A",y="530/30 B-B-A", bins = 256)
autoplot(fs_trans[6:11],x="SSC-A",y="450/50 V-F-A", bins = 256)
autoplot(fs_trans_AS[6:11],x="SSC-A",y="450/50 V-F-A", bins = 256)
