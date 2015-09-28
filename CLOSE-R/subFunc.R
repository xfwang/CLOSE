# all sub-functions in CLOSER

# convert BAF to LAF values
BAFtoLAF<-function(BAF.vec){
      
      # BAF.vec: a vector of BAF values (numeric)

      LAF.vec <- BAF.vec
	  BAF.folding.index<-(1:length(BAF.vec))[BAF.vec>0.5]
         if (length(BAF.folding.index)>1){
         LAF.vec[BAF.folding.index]<-1-LAF.vec[BAF.folding.index]
      }
      return(LAF.vec)
}


# find the density peaks of minor allele copy number ratio (minCNR)
findPeaks<-function(y){
   
   # y: minCNR for all segments; the first element of the list returned by getASCN()
   # return: a list with two elements (fit_mat , peak_mat)
           # 1) fit_mat: Minor_Allelic_CN, density for each segments
           # 2) peak_mat: peak_value, peak_y for each identified peak

   # DPdensity parameters
   mcmc <- list(nburn = 1000, nsave = 1000, nskip = 1, ndisplay = 100)
   prior1 <- list(alpha = 1, m1 = rep(0, 1), psiinv1 = diag(0.5, 1), nu1 = 4, tau1 = 1,  tau2 = 100)
   prior2 <- list(alpha=1,m1=rep(0,1),psiinv2=solve(diag(0.25,1)), nu1=1,nu2=4,tau1=1,tau2=1)
   prior3 <- list(alpha=1,m1=rep(0,1),psiinv2=solve(diag(0.25,1)), nu1=1,nu2=10,tau1=1,tau2=10)


   fit <- DPdensity(y = y, prior = prior2, mcmc = mcmc, state = state, status = TRUE)
   peak_pos<-(1:length(fit$dens))[peaks(density=fit$dens)]
   peak_value<-fit$x1[peak_pos]
   peak_y<-fit$dens[peak_pos]
   peak_mat<-as.data.frame(cbind(peak_value, peak_y))
   if ((max(peak_y)-min(peak_y)) > 0.8*(max(fit$dens)-min(fit$dens))){
        fit <- DPdensity(y = y, prior = prior3, mcmc = mcmc, state = state, status = TRUE)
        peak_pos<-(1:length(fit$dens))[peaks(density=fit$dens)]
        peak_value<-fit$x1[peak_pos]
        peak_y<-fit$dens[peak_pos]
        peak_mat<-as.data.frame(cbind(peak_value, peak_y))
        if ((max(peak_y)-min(peak_y)) > 0.8*(max(fit$dens)-min(fit$dens))){
            fit <- DPdensity(y = y, prior = prior1, mcmc = mcmc, state = state, status = TRUE)
            peak_pos<-(1:length(fit$dens))[peaks(density=fit$dens)]
            peak_value<-fit$x1[peak_pos]
            peak_y<-fit$dens[peak_pos]
            peak_mat<-as.data.frame(cbind(peak_value, peak_y))
        }
   }
  
  fit_mat<-as.data.frame(cbind(fit$x1, fit$dens))
  colnames(fit_mat)<-c("Minor_Allelic_CN","density")
  colnames(peak_mat)<-c("peak_value", "peak_y")
  return(list(fit_mat, peak_mat))
}


# calualte minCNR (minor copy number ratio), majCNR (major copy number ratio)
getASCN<-function(LAF, R, scale=1){

	  # parameters:
	  # 1) LAF: vector of LAF values
      # 2) ratio of read depth (tumor/normal) for each segments: usualy 2^LRR
      # 3) scale: scale factor for library size
      
      # return: list with 2 elements (minCNR and majCNR for each segment)

      minCNR<-2*LAF*R/scale
      majCNR<-2*R*(1-LAF)/scale
      return(list(minCNR, majCNR))
}


# find centroid of the CNR clsuters
getCentroid<-function(segs.cluster, Ncluster){
	   
       # parameters:
       # 1) segs.cluster 
       # 2) Ncluster: number of clusters; third element of the list returned from runDP()

       # return: centroid of each cluster

       # 1. get the mean values of LAF and LRR for each cluster as the centroid
	   LAF_centroid_vec<-c()
       logR_centroid_vec<-c()
       for (i in 1:Ncluster){
          points<-segs.cluster[segs.cluster[,1]==i,]
          LAF_centroid_vec<-c(LAF_centroid_vec,mean(points[,2]))
          logR_centroid_vec<-c(logR_centroid_vec, mean(points[,3]))
        }
      
       # 2. calculate minCNR and majCNR for centroid
       centroid_min_majCNR<-getASCN(LAF=LAF_centroid_vec, R=2^(logR_centroid_vec), scale=1)
       centroid<-as.data.frame(cbind(centroid_min_majCNR[[1]], centroid_min_majCNR[[2]], rep(0, length(Ncluster) ), rep(16, length(Ncluster)) ))
       colnames(centroid)<-c("minCNR","majCNR", "chr", "shape")
       return(centroid)
}


# get clonality from number of peaks
getClonality<-function(peak.mat){
	# paramters:
	# 1) peak.mat: peak_value, peak_y for each identified peak; the second element of the list returned from findPeaks()

    # return: a single value of clonality, numeric

	clonality<-dim(peak.mat)[1]-2
    if (clonality<=0){
       clonality<-0
    }
    return(clonality)
}


# estimate the copy number status: normal, gainm loss, LOH, ambigous
getCNstatus<-function(centroid, CNR.mat, segs, segs.cluster, Ncluster){

    # 1. get the status for each segments based on cluster information
    # 1.1 estimate the copy number status of each cluster centroid
	normal_index<-(1:dim(centroid)[1])[centroid[,1]>=0.75 & centroid[,1]<=1.25 & centroid[,2]>=0.5 & centroid[,2]<=1.5]
    loss_loh_index<-(1:dim(centroid)[1])[centroid[,1]<=0.3 ]
    neutral_loh_index<-(1:dim(centroid)[1])[centroid[,1]<=0.3 & centroid[,2]>=1.8 & centroid[,2]<=2.2]
    gain_index<-(1:dim(centroid)[1])[centroid[,1]>=0.75 & centroid[,2]>=1.5]
    ambigu_index<-(1:dim(centroid)[1])[-c(normal_index, loss_loh_index, gain_index)]

    cluster_status<-rep("Normal", rep(dim(centroid)[1]))
    cluster_status[loss_loh_index]<-"Loss"
    cluster_status[neutral_loh_index]<-"CN Neutral LOH"
    cluster_status[gain_index]<-"Gain"
    cluster_status[ambigu_index]<-"Ambiguous"

    # 1.2 assign the CN status of cluster centroid to all segments in that cluster
    seg_status_cluster<-rep("N", dim(segs)[1])
    for (i in 1:Ncluster){
      cluster_seg<-(1:dim(segs.cluster)[1])[segs.cluster[,1]==i]
      seg_status_cluster[cluster_seg]<-rep(cluster_status[i], length(cluster_seg))
    }
    

    # 2. assign status to each segments using segments infomation, not cluster
    seg_normal_index<-(1:dim(CNR.mat)[1])[CNR.mat[,1]>=0.75 & CNR.mat[,1]<=1.25 & CNR.mat[,2]>=0.5 & CNR.mat[,2]<=1.5]
    seg_loss_loh_index<-(1:dim(CNR.mat)[1])[CNR.mat[,1]<=0.3]
    seg_neutral_loh_index<-(1:dim(CNR.mat)[1])[CNR.mat[,1]<=0.3 & CNR.mat[,2]>=1.8 & CNR.mat[,2]<=2.2]
    seg_gain_index<-(1:dim(CNR.mat)[1])[CNR.mat[,1]>=0.75 & CNR.mat[,2]>=1.5]
    seg_ambigu_index<-(1:dim(CNR.mat)[1])[-c(seg_normal_index, seg_loss_loh_index, seg_gain_index)]

    seg_status<-rep("Normal", rep(dim(CNR.mat)[1]))
    seg_status[seg_loss_loh_index]<-rep("Loss", length(seg_loss_loh_index))
    seg_status[seg_neutral_loh_index]<-rep("CN Neutral LOH", length(seg_neutral_loh_index))   
    seg_status[seg_gain_index]<-rep("Gain", length(seg_gain_index))
    seg_status[seg_ambigu_index]<-rep("Ambiguous", length(seg_ambigu_index))


    f<-cbind(paste("chr",segs[,5], sep=""), segs[,6:7],segs[,1:2], CNR.mat[,1:2], segs.cluster[,1], seg_status_cluster, seg_status)
    f2<-f[order(segs[,5], f[,2]),]
    colnames(f2)<-c("chr","start","end","LAF_median", "logR_mean", "minCNR", "majCNR", "cluster", "status_cluster", "status_seg")
    return(f2)
}


# get purity from peak values of minor allele CN ratio
getPurity<-function(peak.mat){
    
    peak.value <- peak.mat[,1]
	return(peak.value[length(peak.value)]-peak.value[1])
}


# find peaks from density
peaks <- function(density,span=3){

   z <- embed(density, span)
   s <- span%/%2
   v <- max.col(z) == 1 + s
   result <- c(rep(FALSE,s),v)
   result <- result[1:(length(result)-s)]
   result
}


# plot clusters based on minCNR and majCNR 
plotCNR<-function(sampleName, segs.cluster, CNR.mat, fit.mat, peak.mat, purity, clonality){

	    # plot clusters based on minCNR and majCNR 
	    pdf(paste(sampleName, ".plotCNR.pdf", sep=""), width=14, height=8)

        print(ggplot(segs.cluster, aes_string(x = "LAF_median", y = "LRR_mean", color = "cluster_id", label="chromosome")) +geom_point(aes_string(size="size"), alpha=0.3)+scale_size_continuous(range = c(5,15))+geom_text(aes_string(color="cluster_id"))+xlim(0,0.5)+ylim(-2,2)+theme_bw() +theme(panel.grid.major = element_blank(),panel.border = element_rect(colour = "black"),legend.position="none", plot.title=element_text(face="bold",size=20)))

        print(ggplot(CNR.mat, aes_string(x = "minCNR", y = "majCNR", color = "cluster2", label="chr")) +geom_point(aes_string(size="size"), alpha=0.3)+scale_size_continuous(range = c(5,15))+geom_text(aes_string(color="cluster2"))+xlim(min(0, fit.mat[,1]),max(CNR.mat[,1], fit.mat[,1]))+ylim(0,min(6, max(CNR.mat[,2])))+theme_bw() +theme(panel.grid.major = element_blank(),panel.border = element_rect(colour = "black"),legend.position="none", plot.title=element_text(face="bold", size=20))+geom_hline(aes_string(yintercept=1), linetype="dashed")+geom_vline(aes_string(xintercept=1), linetype="dashed"))

        print(ggplot(fit.mat, aes_string(x="Minor_Allelic_CN", y="density"))+geom_line( color="darkred")+xlim(min(0, fit.mat[,1]),max(CNR.mat[,1], fit.mat[,1]))+ylim(0,max(fit.mat[,2])+2)+geom_text(data=peak.mat, x=peak.mat[,1], y=peak.mat[,2]+0.2, label=paste("peak ",signif(peak.mat[,1],3), sep=""))+ geom_text(data=NULL, x=max(fit.mat[,1])*0.6, y=(max(fit.mat[,2])+2)*0.9, label=paste("Purity parameter = ",signif(purity,3), "\n", "Clonality parameter = ", clonality, sep="" ))+theme_bw() +ggtitle("Dirchlet Process Fit of Minor Allelic CN")+theme(panel.grid.major = element_blank(),panel.border = element_rect(colour = "black"),legend.position="none", plot.title=element_text(face="bold", size=14)))

        #print(ggplot(CNR.mat, aes_string(x = "minCNR", y = "majCNR", color = "cluster3", label="chr")) +geom_point(aes_string(size="size"), alpha=0.3)+scale_size_continuous(range = c(5,15))+geom_text(aes_string(color="cluster3"))+xlim(min(0, fit.mat[,1]),max(CNR.mat[,1], fit.mat[,1]))+ylim(0,min(max(CNR.mat[,2])))+theme_bw() +theme(panel.grid.major = element_blank(),panel.border = element_rect(colour = "black"),legend.position="none", plot.title=element_text(face="bold", size=20))+geom_hline(aes_string(yintercept=1), linetype="dashed")+geom_vline(aes_string(xintercept=1), linetype="dashed"))

        dev.off()
 
}


# plot copy number status estimated from getCNstatus chr by chr
plotCNstatus.chr<-function(CNstatus, BAF, LRR, sampleName){
   
     # paramters:
     # 1) CNstatus: output from getCNstatus()
     # 2) BAF matrix: three columns (chr, position, BAF values)
     # 3) LRR matrix: three columns (chr, position, log2(read depth ratio) values)
     # 3) sampleName: string; output prefix

     chr.list<-as.data.frame(table(CNstatus[,1]))[,1]

     for (j in 1:length(chr.list)){
            LRR.chr<-LRR[LRR[,1]==as.character(chr.list[j]),]
            BAF.chr<-BAF[BAF[,1]==as.character(chr.list[j]),]
            seg.chr<-CNstatus[CNstatus[,1]==as.character(chr.list[j]),]
            seg.LRR.mat<-data.frame(start=seg.chr[,2], end=seg.chr[,3],LRR=seg.chr[,5],type=rep("LRR", dim(seg.chr)[1]))
            seg.BAF.mat<-data.frame(start=seg.chr[,2], end=seg.chr[,3],BAF=seg.chr[,4],type=rep("BAF", dim(seg.chr)[1]))


            if (dim(LRR.chr)[1]>0 && dim(BAF.chr)[1]>0){
                     pdf(paste(sampleName,".",chr.list[j], ".pdf", sep=""), width = 8, height=5.5)

                    status<-rep(0, dim(seg.chr)[1]+5)
                    status[(1:dim(seg.chr)[1])[seg.chr[,9]=="Gain"]]<-1
                    status[(1:dim(seg.chr)[1])[seg.chr[,9]=="Loss" | seg.chr[,9]=="CN Neutral LOH" ]]<--1
                    status[length(status)-1]<-1.5
                    status[length(status)]<- -1.5
                    status.label<-c(as.character(seg.chr[,9]), "Normal","CN Neutral LOH","Ambiguous","Gain","Loss")

                             combined<-data.frame(Position=c(as.numeric(BAF.chr[,2]), as.numeric(LRR.chr[,2]) ,  as.numeric(seg.chr[,2]),c(-100000001:-100000005)) , start.value=c(as.numeric(BAF.chr[,3]), as.numeric(LRR.chr[,3]) ,  as.numeric(status)), end=c(as.numeric(BAF.chr[,2]), as.numeric(LRR.chr[,2]) , as.numeric(seg.chr[,3]),c(-100000001:-100000005)), end.value=c(as.numeric(BAF.chr[,3]), as.numeric(LRR.chr[,3]) , as.numeric(status)), status.label=c(rep(NA,dim(BAF.chr)[1]), rep(NA, dim(LRR.chr)[1]), status.label), type=c(rep("BAF", dim(BAF.chr)[1]), rep("LRR", dim(LRR.chr)[1]), rep("CN Status", 5+dim(seg.chr)[1])))
                    combined$type<-factor(combined$type, levels=c("CN Status", "LRR","BAF")) # fix the order of facet panels
                    combined$status.label=factor(combined$status.label, levels=c("Normal","Gain","Loss","CN Neutral LOH","Ambiguous"))

                    g0<-ggplot(subset(combined,type=="CN Status"), aes(Position, start.value))+ coord_cartesian(xlim=c(0-1000000,  max(combined[,3])+1000000))+theme_bw() +theme(strip.text.y = element_text(size = 10, face="bold"), axis.title=element_text(face="bold", size="10"), axis.text.y=element_text(face="bold", size="12"))+geom_segment(aes(x = Position,y = start.value,xend = end,yend = end.value,colour = status.label),size=3) +labs(x = "", y = "")+scale_y_continuous(breaks=c(-1,0,1), labels=c("L", "N", "G")) 
                    p0<-ggplotGrob(g0)

                    g<-ggplot(combined, aes(Position, start.value))+ coord_cartesian(xlim=c(0-1000000,  max(combined[,3])+1000000))+theme_bw() +theme(strip.text.y = element_text(size = 14, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"), axis.title=element_text(face="bold", size="14"), legend.position="top",  axis.text.y=element_text(face="bold", size="12"), legend.key = element_blank())+labs(x = "Position", y = "")+facet_grid(type ~ . , scales = "free")

                    g1<-g+geom_segment(data=subset(combined,type=="CN Status"),aes(x = Position,y = start.value,xend = end,yend = end.value,colour = status.label),size=3) + scale_color_manual("Copy Number Status: ", values = c("black","red", "dodgerblue2", "green","grey"  ))
                    g2<-g1+ geom_point(data=subset(combined,type=="LRR"),color="darkgrey", size=1)+geom_segment(data=seg.LRR.mat,aes(x = start,y = LRR,xend = end,yend = LRR),size=1, color="blue")
                    g3<-g2+ geom_point(data=subset(combined,type=="BAF"),color="darkgrey", size=1) +geom_segment(data=seg.BAF.mat,aes(x = start,y = BAF,xend = end,yend = BAF),size=1, color="blue")

                    gg_table <- ggplot_gtable(ggplot_build(g3))
                    gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"
                    gg_table[["grobs"]][[2]]<-p0[["grobs"]][[2]]
                    grid.draw(gg_table)
                    dev.off()                        
            }
                            
     }

}


# plot LRR, BAF, and CNstatus for whole genome
plotCNstatus.WG<-function(CNstatus, BAF, LRR, sampleName){
	    # paramters are the same as plotCNstatus.chr

	    ##### prepare matrix for plot
	    hg19.length<-c(249250621,243199373,198022430, 191154276, 180915260,171115067,159138663,146364022, 141213431,135534747,135006516,133851895, 115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895, 51304566)
        cum.length<-cumsum(hg19.length)
        length(cum.length)
        cum.length<-c(0,cum.length)
        

        # remove chrX and chrY,chrM
        auto.chr<-paste("chr",1:22, sep="")
        LRR.auto<-subset(LRR, LRR[,1]%in%auto.chr)
        BAF.auto<-subset(BAF, BAF[,1]%in%auto.chr)

       
        chr<-apply(LRR.auto,1, function(x) as.numeric(unlist(strsplit(as.character(x[1]), split="r"))[2]))
        LRR.chr<-cbind(chr, LRR.auto[,-1])
        LRR.chr<-LRR.chr[order(LRR.chr[,1],LRR.chr[,2]),]
        row.names(LRR.chr)<-c(1:dim(LRR.chr)[1])
        new.pos<-unlist(apply(LRR.chr, 1, function(x) {as.numeric(x[2])+cum.length[as.numeric(x[1])]}))	
        LRR.plot.mat<-data.frame(pos=as.numeric(new.pos),log2R=LRR.auto[,3])


        chr<-apply(BAF.auto,1, function(x) as.numeric(unlist(strsplit(as.character(x[1]), split="r"))[2]))
        BAF.chr<-cbind(chr, BAF.auto[,2:3])
        BAF.chr<-BAF.chr[order(BAF.chr[,1],BAF.chr[,2]),]
        row.names(BAF.chr)<-c(1:dim(BAF.chr)[1])
        new.pos<-unlist(apply(BAF.chr, 1, function(x) {as.numeric(x[2])+cum.length[as.numeric(x[1])]}))	
        BAF.plot.mat<-data.frame(pos=as.numeric(new.pos),BAF=BAF.auto[,3])

  
        chr<-apply(CNstatus,1, function(x) as.numeric(unlist(strsplit(as.character(x[1]), split="r"))[2]))
        CNstatus.chr<-cbind(chr, CNstatus[,-1])
        CNstatus.chr<-CNstatus.chr[order(CNstatus.chr[,1],CNstatus.chr[,2]),]
        row.names(CNstatus.chr)<-c(1:dim(CNstatus.chr)[1])

        new.start<-unlist(apply(CNstatus.chr, 1, function(x) {as.numeric(x[2])+cum.length[as.numeric(x[1])]}))
        new.end<-unlist(apply(CNstatus.chr, 1, function(x) {as.numeric(x[3])+cum.length[as.numeric(x[1])]}))

        status<-rep(0, dim(CNstatus)[1]+5)
        status[(1:dim(CNstatus.chr)[1])[CNstatus.chr[,9]=="Gain"]]<-1
        status[(1:dim(CNstatus.chr)[1])[CNstatus.chr[,9]=="Loss" | CNstatus.chr[,9]=="CN Neutral LOH" ]]<--1
        status[length(status)-1]<-1.5
        status[length(status)]<- -1.5
        status.label<-c(as.character(CNstatus.chr[,9]), "Normal","CN Neutral LOH","Ambiguous","Gain","Loss")

        CNS.plot.mat<-data.frame(start=c(as.numeric(new.start),c(-1:-5)), end=c(as.numeric(new.end),c(-1:-5)), status.start=as.numeric(status), status.end= as.numeric(status), status.label=status.label)	

  
        combined<-data.frame(start=c(as.numeric(BAF.plot.mat[,1]), as.numeric(LRR.plot.mat[,1]) , as.numeric(CNS.plot.mat[,1])) , start.value=c(as.numeric(BAF.auto[,3]), as.numeric(LRR.auto[,3]) , as.numeric(CNS.plot.mat[,3])), end=c(as.numeric(BAF.plot.mat[,1]), as.numeric(LRR.plot.mat[,1]) , as.numeric(CNS.plot.mat[,2])), end.value=c(as.numeric(BAF.auto[,3]), as.numeric(LRR.auto[,3]), as.numeric(CNS.plot.mat[,4])), status.label=c(rep(NA,dim(BAF.plot.mat)[1]), rep(NA, dim(LRR.plot.mat)[1]), status.label ), type=c(rep("BAF", dim(BAF.plot.mat)[1]), rep("LRR", dim(LRR.plot.mat)[1]), rep("CN Status", dim(CNS.plot.mat)[1])))
        combined$type<-factor(combined$type, levels=c("CN Status", "LRR","BAF")) # fix the order of facet panels
        combined$status.label=factor(combined$status.label, levels=c("Normal","Gain","Loss","CN Neutral LOH","Ambiguous"))


        combined<-subset(combined, start.value>=-2)
        combined<-subset(combined, start.value<=2)

	    ##### plot CN status, BAF, and LRR

        pdf(paste(sampleName,".","CNstatus.WG.pdf", sep=""), width = 8, height=5.5)
        g0<-ggplot(subset(combined,type=="CN Status"), aes(start, start.value))+ coord_cartesian(xlim=c(0, 2881033286))+theme_bw() +theme(strip.text.y = element_text(size = 14, face="bold"), axis.title=element_text(face="bold", size="14"), axis.text.y=element_text(face="bold", size="12"))+geom_segment(aes(x = start,y = start.value,xend = end,yend = end.value,colour = status.label),size=3) +labs(x = "", y = "")+scale_y_continuous(breaks=c(-1,0,1), labels=c("L", "N", "G"))
        p0<-ggplotGrob(g0)

        g<-ggplot(combined, aes(start, start.value))+ coord_cartesian(xlim=c(0, 2881033286))+theme_bw() +theme(strip.text.y = element_text(size = 14, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"), axis.title=element_text(face="bold", size="14"), legend.position="top", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(face="bold", size="12"), legend.key = element_blank())+geom_vline(xintercept = cum.length,colour="darkgrey", size=0.2)+facet_grid(type ~ . , scales = "free")

        g1<-g+geom_segment(data=subset(combined,type=="CN Status"),aes(x = start,y = start.value,xend = end,yend = end.value,colour = status.label),size=3) + scale_color_manual("Copy Number Status: ", values = c("black","red", "dodgerblue2", "green","grey"  ))+ labs(x = "", y = "")
        g2<-g1+ geom_point(data=subset(combined,type=="LRR"),color="grey", size=0.5) + labs(x = "", y = "")
        g3<-g2+ geom_point(data=subset(combined,type=="BAF"),color="grey", size=0.5) + labs(x = "", y = "")

        for (i in 1:22) {	
	        g3<-g3+annotation_custom(grob = textGrob(i,gp=gpar(fontsize=12,fontface=2), hjust=1, rot=90),  xmin = cum.length[i], xmax = cum.length[i+1], ymin = -7, ymax = -7)
        }

       gg_table <- ggplot_gtable(ggplot_build(g3))
       gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"
       gg_table[["grobs"]][[2]]<-p0[["grobs"]][[2]]
       grid.draw(gg_table)
       dev.off()
}



# cluster segments based on LAF_median and LRR_mean using dirichlet process   
runDP <- function(x, y, disp.param = 0.45, max.iter = 100, tolerance = .001){
    
    # parameters:
    # 1) x: x variable (LAF for each segment)
    # 2) y: y variable (LRR for each segment)
    
    # return: 
    # modified from http://statistical-research.com/dirichlet-process-infinite-mixture-models-and-clustering/

    myData <- cbind( x, y )
    n <- dim(myData)[1]
    var.range<-c(0,0)
    var.range[1]<-max(myData[,1])-min(myData[,1])
    var.range[2]<-max(myData[,2])-min(myData[,2])
    small.scale<-which(var.range==min(var.range))
  
    myData[,small.scale]<-myData[,small.scale]*4
    pick.clusters <- rep(1, n)
    k <- 1
  
    mu <- matrix( apply(myData,2,mean), nrow=1, ncol=ncol(myData) )
  
    is.converged <- FALSE
    iteration <- 0
  
    ss.old <- Inf
    ss.curr <- Inf
    itr.cluster<-matrix(NA, nrow=n, ncol=0)
    while ( !is.converged & iteration < max.iter ) { 
        # Iterate until convergence
         iteration <- iteration + 1
         for( i in 1:n ) { 
         # Iterate over each observation and measure the distance each observation' from its mean center's squared distance from its mean
         distances <- rep(NA, k)
      
             for( j in 1:k ){
                 distances[j] <- sum( (myData[i, ] - mu[j, ])^2 ) # Distance formula.
              }
      
             if( min(distances) > disp.param ) { # If the dispersion parameter is still less than the minimum distance then create a new cluster
                 k <- k + 1
                 pick.clusters[i] <- k
                 mu <- rbind(mu, myData[i, ])
              } else {
                 pick.clusters[i] <- which(distances == min(distances))
              } 
          }
    
          ####### Calculate new cluster means

          for( j in 1:k ) {
               if( length(pick.clusters == j) > 0 ) {
                    mu[j, ] <- apply(subset(myData,pick.clusters == j), 2, mean)
               }
          }
    
          ####### Test for convergence
          
          ss.cur<-0
          for (i in 1:n){
          	ss.cur<-ss.curr+sum((myData[i,]-mu[pick.clusters[i],])^2)
          }
          ss.diff<-ss.old-ss.curr
          ss.old<-ss.curr
          if (!is.nan(ss.diff) & ss.diff < tolerance){
                 is.converged<-TRUE
          }



          itr.cluster<-cbind(itr.cluster,pick.clusters)
    }
  
    centers <- data.frame(mu)
    ret.val <- list("centers" = centers, "cluster" = factor(pick.clusters),
                  "k" = k, "iterations" = iteration, "itr_cluster" = itr.cluster)
  
    return(ret.val)
}


runFalcon<-function(input, sampleName, threshold=0.15){
       # function that runs Falcon
       # threshold: for a given segments,if for both alleles the abs(copy numbers  -1) < threshold, this segment will not be considered as CNV

       library(falcon)

       AN = as.numeric(input[,5])
       BN = as.numeric(input[,6])
       AT = as.numeric(input[,7])
       BT = as.numeric(input[,8])
       pos = as.numeric(input[,2])

       rdep = median(AT+BT)/median(AN+BN)

       result = c()

       for (i in 1:22){
           name = paste("chr",i,sep="")
           ids = which(input[,1]==name)
           cn = getASCN(AT[ids],BT[ids],AN[ids],BN[ids],rdep=rdep)
           #save(cn,file=paste(sampleName,"_chr",i,".Rdata",sep=""))
           png(paste(sampleName,"_chr",i,".png",sep=""),width=1200,height=750,pointsize=23)
           view(AT[ids],BT[ids],AN[ids],BN[ids],cn$tauhat,cn$ascn,pos[ids],rdep=rdep)
           dev.off()

           tauhat0 = c(1,cn$tauhat,length(ids)+1)
           m = length(tauhat0)
           start = tauhat0[1:(m-1)]
           end = tauhat0[2:m]-1
           temp = cbind(rep(i,m-1),ids[start], ids[end], pos[ids[start]], pos[ids[end]], round(t(cn$ascn),digits=3))
           result = rbind(result, temp)
       }

       removeid = c()

       for (j in 1:dim(result)[1]){
           a = result[j,6]
           b = result[j,7]
           if (abs(a-1)<threshold && abs(b-1)<threshold){
                removeid = c(removeid, j)
           }
       }

 colnames(result) = c("chr", "SNP.start", "SNP.end", "pos.start", "pos.end", "minor.CN", "major.CN")
 #return(result[-removeid,])
 return(result)
 # write.table(result[-removeid,], file=paste(samplename,"_ascn_table.txt",sep=""), quote=F, row.names=F, sep="\t")

}



# prepare segmentations for CLOSE-R
segPrep<-function(Input){
    
    # paramters:
    # 1) Input:the Input file of CLOSER 

    # return: LAF_median and LRR_mean for each segment

    size_vec = region_length_vec = c()
    for (i in 1:dim(Input)[1]){
            region_length<-as.numeric(Input[i,3])-as.numeric(Input[i,2])+1
            region_length_vec<-c(region_length_vec, region_length)
            size<-as.numeric(region_length)/1000000
            if (size<1) {size<-1}
            if (size>8) {size<-8}
            size_vec<-c(size_vec,size)
    }

    mat<-cbind(Input[,4:5], region_length_vec,  size_vec, Input[,1:3] )
    return(mat)

}




VCFprep<-function(Normal.vcf, Tumor.vcf,  filter = FALSE){

       # calculate BAF and LRR from the .vcf files of Normal and Tumor sample

       # paramters:
       # 1) Normal.vcf and Tumor.vcf are the path to the two .vcf files
       # 2) filter: logical. Whether to filter the .vcf files based on FILTER (PASS)
       
       library(VariantAnnotation)

       N.vcf = readVcf(Normal.vcf, "hg19") 
       if (!("AO"%in%names(geno(N.vcf)) && "RO"%in%names(geno(N.vcf)))) {
            #N.valid=F
                 stop("Normal vcf failed vcf foramt check: No RO/AO information!\n")
            }else{
                  cat ("Normal.vcf passed vcf foramt check\n")
            }
       T.vcf = readVcf(Tumor.vcf, "hg19")
       if (!("AO"%in%names(geno(T.vcf)) && "RO"%in%names(geno(T.vcf)))) {
            #T.valid=F
            stop("Tumor vcf failed vcf foramt check: No RO/AO information!\n")
            }else{
                  cat ("Tumor.vcf passed vcf foramt check\n")
            }
       

       T.FILTER = filt(T.vcf)
       N.FILTER = filt(N.vcf)


       if (filter) {
            T.pass = which(T.FILTER=="PASS")
            N.pass = which(N.FILTER=="PASS")
       } else {
            T.pass = 1:length(T.FILTER)
            N.pass = 1:length(N.FILTER)
       }
       
       T.FILTER = T.FILTER[T.pass]
       N.FILTER = N.FILTER[N.pass]

	   T.POS = start(ranges(rowRanges(T.vcf)))[T.pass]
	   N.POS = start(ranges(rowRanges(N.vcf)))[N.pass]

       T.CHR = as.character(seqnames(rowRanges(T.vcf)))[T.pass]
       N.CHR = as.character(seqnames(rowRanges(N.vcf)))[N.pass]


       #T.GT = t(sapply(geno(T.vcf)$GT[T.pass,], function(x) unlist(strsplit(as.character(x), split="/"))))
       #N.GT = t(sapply(geno(N.vcf)$GT[N.pass,], function(x) unlist(strsplit(as.character(x), split="/"))))

       T.GT= as.vector(geno(T.vcf)$GT[T.pass,])
       N.GT= as.vector(geno(N.vcf)$GT[N.pass,])

       T.REF = as.vector(ref(T.vcf))[T.pass]
       N.REF = as.vector(ref(N.vcf))[N.pass]

       T.ALT = unstrsplit(CharacterList(alt(T.vcf)), sep = ",")[T.pass]
       N.ALT = unstrsplit(CharacterList(alt(N.vcf)), sep = ",")[N.pass]

       T.AO = unstrsplit(CharacterList(geno(T.vcf)$AO), sep=",")[T.pass]
       N.AO = unstrsplit(CharacterList(geno(N.vcf)$AO), sep=",")[N.pass]

       T.RO = as.vector(geno(T.vcf)$RO)[T.pass]
       N.RO = as.vector(geno(N.vcf)$RO)[N.pass]

       T.mat = cbind(T.CHR, T.POS, T.REF, T.ALT, T.GT, T.AO, T.RO)
       N.mat = cbind(N.CHR, N.POS, N.REF, N.ALT, N.GT, N.AO, N.RO)

       T.DP = as.vector(geno(T.vcf)$DP)[T.pass]
       N.DP = as.vector(geno(N.vcf)$DP)[N.pass]

       normalize.factor = sum(as.numeric(T.DP))/(sum(as.numeric(N.DP)))
       ### match Normal and Tumor positions

       T.index = N.index = c()


       for (i in 1:22){
            name = paste("chr", i, sep="")
            T.id = which(T.CHR==name)
            N.id = which(N.CHR==name)
            T.pos = T.POS[T.id]
            N.pos = N.POS[N.id]
            comp = match(T.pos, N.pos)
            T.match = (1:length(T.id))[-which(is.na(comp))]
            N.match = comp[-which(is.na(comp))]
            T.index = c(T.index,T.id[T.match])
            N.index = c(N.index,N.id[N.match])

       }

       T.match.mat = T.mat[T.index,]
       N.match.mat = N.mat[N.index,]
       
       ### get the heterozygous SNPs in normal sample
       N.GT.mat = matrix(as.numeric(sapply(as.list(N.match.mat[,5]), function(x) {strsplit(as.character(x), split="/")[[1]]})), ncol=2, byrow=T)

       dif = N.GT.mat[,2] - N.GT.mat[,1]
       het = which(dif>0)
       
       N.matrix.het = N.match.mat[het,]
       T.matrix.het = T.match.mat[het,]
       mydata = cbind(N.matrix.het[,c(1,2,3)], T.matrix.het[,3], N.matrix.het[,4], T.matrix.het[,4], N.matrix.het[,5],T.matrix.het[,5],N.matrix.het[,6],T.matrix.het[,6], N.matrix.het[,7], T.matrix.het[,7])
       colnames(mydata) = c("chr", "pos", "N.ref", "T.ref", "N.alt", "T.alt", "N.GT", "T.GT", "N.AO", "T.AO", "N.RO","T.RO")

       ### remove the SNPs where Tumor and normal sample has different ref or alt alleles
       ref.dif = which(as.character(mydata[,3])!=as.character(mydata[,4]))
       alt.dif = which(as.character(mydata[,5])!=as.character(mydata[,6])) # ref.dif is uauslly a subset of alt.dif

       mydata1 = mydata[-unique(c(ref.dif, alt.dif)),]

       ### remove the SNPs that has more than 2 alternative alleles in normal sample
       mydata2 = mydata1[which(as.character(mydata1[,7])=="0/1" | as.character(mydata1[,7])=="0/2" | as.character(mydata1[,7])=="1/2"),]
       #N.GT.mat = matrix(as.numeric(sapply(as.list(mydata2[,7]), function(x) {strsplit(as.character(x), split="/")[[1]]})), ncol=2, byrow=T)
       #T.GT.mat = matrix(as.numeric(sapply(as.list(mydata2[,8]), function(x) {strsplit(as.character(x), split="/")[[1]]})), ncol=2, byrow=T)
 
       mydata2 = mydata2[!apply(mydata2,1,function(x) {max(unlist(strsplit(x[7], split="/"))) !=  max(unlist(strsplit(x[8], split="/")))}),]
       ### get number of reads for the two alleles A and B
       ### AN: number of reads for allele A
       ### BN: number of reads for allele B

       reads = matrix(0,dim(mydata2)[1],4)
       for (i in 1:dim(mydata2)[1]) {
            #cat("i is ", i, "\n")
            if (as.character(mydata2[i,7])=="0/1"){
                  reads[i,]=c(as.numeric(mydata2[i,11]), as.numeric(mydata2[i,9]), as.numeric(mydata2[i,12]), as.numeric(mydata2[i,10]))
            }else {
                  temp1 = as.numeric(unlist(strsplit(as.character(mydata2[i,9]), split=",")))
                  temp2 = as.numeric(unlist(strsplit(as.character(mydata2[i,10]),split=",")))
                
                  if (mydata2[i,7]=="1/2") {
                        reads[i,] = c(temp1, temp2)
                  } else {
                        reads[i,] = c(as.numeric(mydata2[,11]),temp1[2], as.numeric(mydata2[,12]),temp2[2])
                  }
             }
       }

       
       Afreq.N = round(reads[,1]/(reads[,1]+reads[,2]),3)
       Afreq.T = round(reads[,3]/(reads[,3]+reads[,4]),3)
       LRR = log(((reads[,3]+reads[,4])/normalize.factor)/(reads[,1]+reads[,2]),2)
       output = cbind(mydata2[,c(1,2,3,5)], reads, Afreq.N, Afreq.T, LRR)
       colnames(output) = c("chr", "pos", "ref", "alt", "AN", "BN", "AT", "BT", "Afreq.N", "Afreq.T", "LRR")

       return(output)
}




# compare copy number estimated by CLOSER to TCGA resutls
compareToTCGA<-function(CNstatus, TCGA,sampleName){
	  # parameters:
	  # 1) CNstatus: output from getCNstatus()
      # 2) TCGA: matrix with 4 columns 
            #(1) chr: integer
            #(2) start: numeric
            #(3) end: numeric
            #(4) seg_mean from TCGA data
      # 3) sampleName: output prefix


	  	hg19.length<-c(249250621,243199373,198022430, 191154276, 180915260,171115067,159138663,146364022, 141213431,135534747,135006516,133851895, 115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895, 51304566)
        cum.length<-cumsum(hg19.length)
        length(cum.length)
        cum.length<-c(0,cum.length)

        new.start<-unlist(apply(TCGA, 1, function(x) {as.numeric(x[2])+cum.length[as.numeric(x[1])]}))
        new.end<-unlist(apply(TCGA, 1, function(x) {as.numeric(x[3])+cum.length[as.numeric(x[1])]}))

        col<-rep("Normal", dim(TCGA)[1])
        col[TCGA[,4]< -0.05]<-"Loss"
        col[TCGA[,4]> 0.05]<-"Gain"
        seg_mean<-TCGA[,4]
        for (i in 1:length(seg_mean)){
            if (TCGA[i,4]>2){
               seg_mean[i]=2
             }else if (TCGA[i,4]< -2){
              seg_mean[i]=-2
             }
        }

        TCGA.mat<-data.frame(start=new.start, end=new.end, segment_mean.start=rep(0,length(seg_mean)), segment_mean.end=as.numeric(seg_mean),col=factor(col))


        chr<-apply(CNstatus,1, function(x) as.numeric(unlist(strsplit(as.character(x[1]), split="r"))[2]))
        CNstatus.chr<-cbind(chr, CNstatus[,-1])
        CNstatus.chr<-CNstatus.chr[order(CNstatus.chr[,1],CNstatus.chr[,2]),]
        row.names(CNstatus.chr)<-c(1:dim(CNstatus.chr)[1])

        new.start<-unlist(apply(CNstatus.chr, 1, function(x) {as.numeric(x[2])+cum.length[as.numeric(x[1])]}))
        new.end<-unlist(apply(CNstatus.chr, 1, function(x) {as.numeric(x[3])+cum.length[as.numeric(x[1])]}))

        status<-rep(0, dim(CNstatus)[1]+5)
        status[(1:dim(CNstatus.chr)[1])[CNstatus.chr[,9]=="Gain"]]<-1
        status[(1:dim(CNstatus.chr)[1])[CNstatus.chr[,9]=="Loss" | CNstatus.chr[,9]=="CN Neutral LOH" ]]<--1
        status[length(status)-1]<-1.5
        status[length(status)]<- -1.5
        status.label<-c(as.character(CNstatus.chr[,9]), "Normal","CN Neutral LOH","Ambiguous","Gain","Loss")

        CNS.plot.mat<-data.frame(start=c(as.numeric(new.start),c(-1:-5)), end=c(as.numeric(new.end),c(-1:-5)), status.start=as.numeric(status), status.end= as.numeric(status), status.label=status.label)	


        combined<-data.frame(start=c(as.numeric(CNS.plot.mat[,1]), as.numeric(TCGA.mat[,1])) , start.value=c(as.numeric(CNS.plot.mat[,3]), as.numeric(TCGA.mat[,3])), end=c(as.numeric(CNS.plot.mat[,2]), as.numeric(TCGA.mat[,2])), end.value=c(as.numeric(CNS.plot.mat[,4]), as.numeric(TCGA.mat[,4])), status.label=c(status.label,col), type=c( rep("CLOSE-R", dim(CNS.plot.mat)[1]), rep("SNP Array", dim(TCGA.mat)[1])))
        combined$type<-factor(combined$type, levels=c("CLOSE-R", "SNP Array")) # fix the order of facet panels
        combined$status.label=factor(combined$status.label, levels=c("Normal","Gain","Loss","CN Neutral LOH","Ambiguous"))



        pdf(paste(sampleName,".","compareToTCGA.pdf", sep=""), width = 8, height=5)
        g0<-ggplot(subset(combined,type=="CLOSE-R"), aes(start, start.value,color=status.label))+ coord_cartesian(xlim=c(0, 2881033286))+theme_bw() +theme(strip.text.y = element_text(size = 14, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"), axis.title=element_text(face="bold", size="14"), legend.position="top", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(face="bold", size="12"), legend.key = element_blank())+geom_vline(xintercept = cum.length,colour="darkgrey", size=0.2)+geom_segment(aes(x = start,y = start.value,xend = end,yend = end.value,colour = status.label),size=3) +scale_color_manual("Copy Number Status: ", values = c("black","red", "dodgerblue2", "green","grey"  ))+labs(x = "", y = "")+scale_y_continuous(breaks=c(-1,0,1), labels=c("L", "N", "G"))
        p0<-ggplotGrob(g0)

        g_0<-ggplot(subset(combined,type=="SNP Array"), aes(start, start.value,color=status.label,fill=status.label))+ coord_cartesian(ylim=c(-2,2),xlim=c(0, 2881033286))+theme_bw() +theme(strip.text.y = element_text(size = 14, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"), axis.title=element_text(face="bold", size="14"), legend.position="top", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(face="bold", size="12"), legend.key = element_blank())+geom_vline(xintercept = cum.length,colour="darkgrey", size=0.2)+geom_segment(aes(x = start,y = start.value,xend = end,yend = end.value,colour = status.label),size=3) +scale_color_manual("Copy Number Status: ", values = c("black","red", "dodgerblue2", "green","grey"  ))+labs(x = "", y = "")+geom_rect( aes(xmin = start,xmax=end, ymin = 0, ymax=end.value, color=status.label,fill=status.label))
        p_0<-ggplotGrob(g_0)


        g<-ggplot(combined, aes(start, start.value, color=status.label,fill=status.label))+ coord_cartesian(xlim=c(0, 2881033286))+theme_bw() +theme(strip.text.y = element_text(size = 14, face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_rect(colour = "black"), axis.title=element_text(face="bold", size="14"), legend.position="top", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(face="bold", size="12"), legend.key = element_blank())+geom_vline(xintercept = cum.length,colour="darkgrey", size=0.2)+scale_color_manual("Copy Number Status: ", values = c("black","red", "dodgerblue2", "green","grey"  ))+scale_fill_manual("Copy Number Status: ", values = c("black","red", "dodgerblue2", "green","grey"  ))+ labs(x = "", y = "")+facet_grid(type ~ . , scales = "free")

        g1<-g+geom_segment(data=subset(combined,type=="CLOSE-R"),aes(x = start,y = start.value,xend = end,yend = end.value, color=status.label,fill=status.label),size=3) 
        g2<-g1+geom_rect(data=subset(combined,type=="SNP Array"), aes(xmin = start,xmax=end, ymin = 0, ymax=end.value, color=status.label,fill=status.label))
        for (i in 1:22) {	
	        g2<-g2+annotation_custom(grob = textGrob(i,gp=gpar(fontsize=12,fontface=2), hjust=1, rot=90),  xmin = cum.length[i], xmax = cum.length[i+1], ymin = -2.2, ymax = -2.2)
        }

       gg_table <- ggplot_gtable(ggplot_build(g2))
       gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"
       gg_table[["grobs"]][[2]]<-p0[["grobs"]][[2]]
       gg_table[["grobs"]][[3]]<-p_0[["grobs"]][[2]]

       grid.draw(gg_table)
       dev.off()

}








