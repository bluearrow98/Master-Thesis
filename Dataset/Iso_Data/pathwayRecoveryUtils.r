divideIsoData <- function(path){
    isodata <- read.csv("ExpressionData.csv", header=TRUE, stringsAsFactors=FALSE)

    isodata <- isodata[,2:ncol(isodata)]

    isodata_chloropl <- isodata[,c("DXPS1","DXPS2.cla1.","DXPS3","DXR","MCT","CMK","MECPS","HDS","HDR","IPPI1","GPPS","GGPPS2","GGPPS6","GGPPS8","GGPPS10","GGPPS11","GGPPS12","PPDS1","PPDS2mt")]

    isodata_mitochon<- isodata[,c("UPPS1","DPPS2","GGPPS1mt","GGPPS5","GGPPS9")]

    isodata_cytopl <- isodata[,c("AACT1","AACT2","HMGS","HMGR1","HMGR2","MK","MPDC1","MPDC2","IPPI2","FPPS1","FPPS2","DPPS1","DPPS3","GGPPS3","GGPPS4")]

    if (path == 1) {
        return(isodata_chloropl)
    }else if (path == 2){
        return(isodata_mitochon)
    }else{
        return(isodata_cytopl)
    }
}


setPathwayProb <- function(data) {

signal <- list()

# Number of variables
p <- dim(data)[2]

# Data scaling
X <- as.matrix(data)
X <- scale(X)


# Construct a Volatility matrix
C <- diag(p)


signal$p <- p
signal$n <- dim(X)[1]
signal$C <- C
signal$X <- X

return(signal)
}