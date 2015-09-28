CLOSER<-function(Input, codePath, output, sampleName="closer){

	   # CLOSE-R main functions
	   # calls all the sub functions inside of subFunc.R

       setwd(output)
       library("DPpackage")
       library(grid)
       library(ggplot2)

       R_scripts = paste(codePath, list.files(path=(codePath), pattern=".R$"), sep="/")

       for (i in 1:length(R_scripts)){
	       source(R_scripts[i])
       }
       set.seed(123)


       ############ segmentation preparation
       segs = segPrep(Input=Input)
       # output: LAF_median and LRR_mean for each segment


       ############ Dirichlet Process

       DPcluster = runDP(x=segs[,1], y=segs[,2], disp.param = 0.45, max.iter = 30, tolerance = .001)
       # cluster segments based on LAF_median and LRR_mean
       segs.cluster = as.data.frame(cbind(DPcluster[[2]], segs))
       colnames(segs.cluster)<-c("cluster_id", "LAF_median","LRR_mean","length","size","chromosome","start","end")
       Ncluster = DPcluster[[3]]

       ############ calculate minCNR and majCNR

       CNR = getASCN(LAF=segs[,1], R=2^(segs[,2]), scale=1)
       CNR.mat = as.data.frame(cbind(CNR[[1]], CNR[[2]], factor(segs.cluster[,1]), segs.cluster[,5],segs.cluster[,6]))
       CNR.mat[,3] =factor(CNR.mat[,3])
       colnames(CNR.mat) = c("minCNR", "majCNR", "cluster2", "size", "chr")     

       ############ find density peaks

       CNR.fit = findPeaks(y=CNR[[1]])
       fit.mat = as.data.frame(CNR.fit[[1]])
       peak.mat = as.data.frame(CNR.fit[[2]])

       purity = getPurity(peak.mat = peak.mat)
       clonality = getClonality(peak.mat = peak.mat)


       ############ plot CNR clusters

       plotCNR(sampleName =sampleName, segs.cluster = segs.cluster, CNR.mat =CNR.mat, fit.mat =fit.mat, peak.mat = peak.mat, purity = purity, clonality = clonality)

       ############ estimate copy number status

       centroid = getCentroid(segs.cluster = segs.cluster, Ncluster = Ncluster)
       CNstatus = getCNstatus(centroid, CNR.mat, segs, segs.cluster, Ncluster)

       write.table(CNstatus, paste(sampleName,  ".CNstatus.txt", sep=""), col.names=colnames(CNstatus), row.names=F, quote=F)
       # output: .CNstatus.txt
}


