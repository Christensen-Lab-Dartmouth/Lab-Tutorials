# Original custom plotting function from RPMM authors
# New function
plotMethByClass = function(BETA, CLASS, sep="red",
                           col=colorRampPalette(c("yellow","black","blue"),space="Lab")(64)){
  
  nColor = length(col)
  nCpG = dim(BETA)[2]
  nSubj = dim(BETA)[1]
  index = split(1:nSubj, CLASS)
  nClass = length(index)
  plot(c(0,nCpG), c(0,nSubj), type="n", xlab="", ylab="",
       xaxt="n", yaxt="n", bty="n")
  
  ordCpG = hclust(dist(t(BETA)), method="ward")$ord
  BETA = BETA[,ordCpG]
  
  k=0
  for(i in 1:nClass){
    ii = index[[i]]
    nii = length(ii)
    ord = hclust(dist(BETA[ii,]), method="ward")$ord
    for(j in 1:nii){
      colori = ceiling(BETA[ii[ord[j]],]*nColor)
      rect(1:nCpG-1,k,1:nCpG,k+1,col=col[colori],density=-1,border=NA)
      k = k+1
    }
    if(i<nClass) lines(c(0,nCpG), c(k,k), col=sep)
  }
  nn = cumsum(c(0,unlist(lapply(index,length))))
  axp = 0.5*(nn[-1]+nn[-nClass-1])
  axis(2, axp, names(index), las=2)
}

## NOT RUN
# plotMethByClass(autoB,RCCclass1)
## END NOT RUN
