#Global MAP estimation
#xuefeng.wang@stonybrook.edu; xuefeng.wang@yale.edu


##function to estimate copy number given purity and ploidy
getBC<-function(x, y, la, p){
   #x:LRR; y: LAF; la: ploidy* ; p: purity
   #B: minor allele;  C; total allele
   C<-(2^x*la)+2*((1-p)/p)*(2^x-1)
   B<-C*y+(2*y-1)*((1-p)/p)
  return(c(B,C))
 }
getBA<-function(x, y, la, p,bias=0){
   #x:LRR; y: LAF; la: ploidy* ; p: purity
   #B: minor allele;  C; total allele
   y<-ifelse(abs(y-0.5)<bias,0.5, y)
   C<-(2^x*la)+2*((1-p)/p)*(2^x-1)
   B<-C*y+(2*y-1)*((1-p)/p)
   A=C-B
  return(c(B,A))
 }
#getBC(x=0,y=1/4,la=4,p=1)
#getBC(x=0,y=1/4,la=4,p=1)


##function to calculate mean.LRR and mean.LAF (given purity and ploidy)
getLL<-function(B, C, la, p, bias=0.04){
  #B: la: ploidy ; p: purity; B: B allele;  C; total allele
  mu.R<- log(((C*p+2*(1-p))/(la*p+2*(1-p))),2)
  #mu.R<- log(((C*p+2*(1-p))/(2*la),2)
  mu.L<- (B*p+(1-p))/(C*p+2*(1-p)) #become LAF!
  mu.L<-min(mu.L,1-mu.L)
  if (mu.L==0.5) mu.L=mu.L-bias #control for mapping bias (can be determined by checking all samples )
  return(c(mu.R,mu.L))
}
#getLL(B=2,C=2,la=2,p=1,bias=0)


##function to generate different B allele and total allele # combinations
BC.expand<-function(x=7) {
  #x is number of maxium ploidy consider
  amatrix=NULL
  for (i in 1:x) {
     C=i; B=seq(0,floor(C/2)) #only minor
     amatrix[[i]]=cbind(B,C)
  }
  data.out<-do.call(rbind,amatrix)
  row.names(data.out)<-paste(data.out[,1],data.out[,2],sep="/")
  return(data.out)
}

##Canonical point calculation (generate LRR/LAF points)
##output: each row is purity/ploidy comb, each col is allele combo
Cpoints<-function(p_grid,la_grid) {
    #p, la are purity/ploidy grid vector p<-seq(0,1,1) la<-seq(0.1,2.5,0.5)
  PP.matrix<-data.matrix(expand.grid(p_grid,la_grid)); colnames(PP.matrix)<-c("purity","ploidy")
  BC.matrix<-BC.expand(x=5); colnames(BC.matrix)<-c("B","total")
  LRR_LAF.array<-array(dim=c(nrow(PP.matrix),nrow(BC.matrix),2))
  for (i in 1:nrow(PP.matrix)) {
    for (j in 1:nrow(BC.matrix)) {
       #LRR_LAF.array[i,j,]<-ifelse(BC.matrix[j,2]>PP.matrix[i,2], NA, getLL(B=BC.matrix[j,1],C=BC.matrix[j,2],la=PP.matrix[i,2],p=PP.matrix[i,1]))  #if model loss only
       LRR_LAF.array[i,j,]<-getLL(B=BC.matrix[j,1],C=BC.matrix[j,2],la=PP.matrix[i,2],p=PP.matrix[i,1])
    }
  }
  name1=paste(PP.matrix[,1],PP.matrix[,2],sep="/")
  name2=paste(BC.matrix[,1],BC.matrix[,2],sep="/")
  return(list(cLL=LRR_LAF.array, name1=name1, name2=name2))
}
#p_grid<-seq(0,1,0.1); la_grid<-seq(0.1,2.5,0.1)  #la is relative average ploidy
#Cpoints(p_grid,la_grid)

#Canonical line calcualtion (given ploidy)
Clines.demo<-function(p_grid=seq(0,1,0.02), la) {
 #p_grid are purity/ploidy grid vector
 ## p_grid<-seq(0,1,0.02); la<-4
 PP.matrix<-data.matrix(expand.grid(p_grid,la)); colnames(PP.matrix)<-c("purity","ploidy")
 BC.matrix<-BC.expand(x=6); colnames(BC.matrix)<-c("B","total")
 LRR_LAF.array<-array(dim=c(nrow(PP.matrix),nrow(BC.matrix),2))
 MM.array<-array(dim=c(nrow(PP.matrix),nrow(BC.matrix),2))
  for (i in 1:nrow(PP.matrix)) {
    for (j in 1:nrow(BC.matrix)) {
       #LRR_LAF.array[i,j,]<-ifelse(BC.matrix[j,2]>PP.matrix[i,2], NA, getLL(B=BC.matrix[j,1],C=BC.matrix[j,2],la=PP.matrix[i,2],p=PP.matrix[i,1]))  #if model loss only
       LRR_LAF.array[i,j,]<-getLL(B=BC.matrix[j,1],C=BC.matrix[j,2],la=PP.matrix[i,2],p=PP.matrix[i,1],bias=0)
       MM.array[i,j,]<-getBA(x=LRR_LAF.array[i,j,1],y=LRR_LAF.array[i,j,2],la=PP.matrix[i,2],p=1)

    }
  }
  name1=paste(PP.matrix[,1],PP.matrix[,2],sep="/")
  name2=paste(BC.matrix[,1],BC.matrix[,2],sep="/")

  #par(mfrow=c(2,2))
  pdf("4.pdf",width=10, height=5); layout(matrix(1:2, 1, 2,byrow=T))

  #g= 1: 0/1" 2 "0/2" 3 "1/2" 4 "0/3" 5 "1/3" 6 "0/4" 7 "1/4" 8 "2/4"
  for (g in 1:15) {
    xdata=LRR_LAF.array[,g,2]; ydata=LRR_LAF.array[,g,1]
    if (g==1) {plot(xdata,ydata, xlim=c(0,0.5),ylim=c(-2,2),type="l",col=g,xlab="LAF",ylab="LRR"); text(mean(xdata),mean(ydata), name2[g])} else {
    points (xdata,ydata, xlim=c(0,0.5),ylim=c(-2,2),type="l",col=g); text(mean(xdata),mean(ydata),name2[g])}
  }

  for (g in 1:15) {
    ydata=MM.array[,g,2]; xdata=MM.array[,g,1]
    if (g==1) {plot(xdata,ydata, xlim=c(0,3),ylim=c(0,3),type="l",col=g,xlab="minor allele CN",ylab="major allele CN"); text(mean(xdata),mean(ydata), name2[g])} else {
    points (xdata,ydata, xlim=c(0,3),ylim=c(0,3),type="l",col=g); text(mean(xdata),mean(ydata),name2[g])}
  }
  dev.off()
}
#Clines.demo(la=2); Clines.demo(la=4)
Clines.demo(la=4)

##Function to find the mode of a distribution
#From http://stackoverflow.com/questions/16255622/peak-of-the-kernel-density-estimation
dmode <- function(x) {
      den <- density(x, kernel=c("gaussian"))
        ( den$x[den$y==max(den$y)] )
    }
n.modes <- function(x) {
       den <- density(x, kernel=c("gaussian"))
       den.s <- smooth.spline(den$x, den$y, all.knots=TRUE, spar=0.8)
         s.0 <- predict(den.s, den.s$x, deriv=0)
         s.1 <- predict(den.s, den.s$x, deriv=1)
       s.derv <- data.frame(s0=s.0$y, s1=s.1$y)
       nmodes <- length(rle(den.sign <- sign(s.derv$s1))$values)/2
       if ((nmodes > 10) == TRUE) { nmodes <- 10 }
          if (is.na(nmodes) == TRUE) { nmodes <- 0 }
       ( nmodes )
}
peaks <- function(density,span=3){

   # find peaks from density

   z <- embed(density, span)
   s <- span%/%2
   v <- max.col(z) == 1 + s
   result <- c(rep(FALSE,s),v)
   result <- result[1:(length(result)-s)]
   result
}
##x <- runif(1000,0,100); plot(density(x)); abline(v=dmode(x))

##laplace denstiy and weight function
getLAFwt1<-function(laf){
     laplace <-function (y,m=0,s=4){
      #f(y) = exp(-abs(y-m)/s)/(2*s)
     d=exp(-abs(y-m)/s)/(2*s)
     return(d)
    }
   wt<-laplace(laf*10,s=1.5)*10+1  #change s for different wts
   return(wt)
}
#laf<-seq(0,0.5,0.01)
#plot(laf,getLAFwt(laf),ylim=c(0,6))


##main function to estimate global purity and ploidy
#p_grid<-seq(0.1,1,0.1); la_grid<-seq(1,6,1)  #la is ploidy
#CanP<-Cpoints(p_grid,la_grid)
fitPP<-function(OLL,slength,nor.rm=TRUE,nor.wt=FALSE, LAF.wt=c(1,6),LRR.center=TRUE){
   #input oLL is the observed LRR_LAF matrix, each per column, l
   #slength: length of each segments
   #
   ##filter noninformative segments
   mydata<-cbind(OLL,slength)
   mydata<-mydata[mydata[,1]>-2&mydata[,1]<2,]  #remove all LRR that is smaller than -2 and larger than 2
   if (LRR.center==TRUE) {
         mydata[,1]<-scale(mydata[,1],center=median(mydata[,1]),scale=FALSE) #centralized around LRR median
   }
   #remove normal segments
   #mydata<-cbind(mydata, (abs(mydata[,1]-0)<0.05)&(abs(mydata[,2]-0.5)<0.05))
   #mydata<-mydata[mydata[,4]==0,]
   if (nor.wt==TRUE) {

        LAF.peak<-dmode(mydata[,2])
        foo1<<-LAF.peak
        mydata<-cbind(mydata, (abs(mydata[,1]-0)<0.2)&(abs(mydata[,2]-LAF.peak)<0.05))
        mydata[mydata[,4]==1,3]=0.5*(mydata[mydata[,4]==1,3])

   }
   if (nor.rm==TRUE) {
        #foo1<<-LAF.peak
        LAF.peak<-dmode(mydata[,2])
        mydata<-cbind(mydata, (abs(mydata[,1]-0)<0.2)&(abs(mydata[,2]-LAF.peak)<0.05))
        mydata<-mydata[mydata[,4]==0,]
   }
   foo2<<-mydata
   OLL<-mydata[,1:2];   slength<-mydata[,3]

   if (LAF.wt[1]==1) {
      slength[OLL[,2]<0.25]=LAF.wt[2]*(slength[OLL[,2]<0.25])
   }
   if (LAF.wt[1]==2) {
      b=LAF.wt[2]; a=(1-b)*2
      slength=(a*OLL[,2]+b)*slength
   }

   sum.dist=NULL
   for (i in 1: nrow(CanP$cLL[,,1])) {
     # i=20
     LL.coor<-cbind(CanP$cLL[i,,1], CanP$cLL[i,,2])
     min.dist<-NULL; copy<-NULL
     for (j in 1:nrow(OLL)) {
        dist<-apply(LL.coor,1,function(x)dist(rbind(x,OLL[j,])))
        min.dist[j]<-min(dist,na.rm=TRUE)
        temp<-data.frame(CanP$name2,dist)
        copy[j]<-as.character(temp[order(temp[,2]),][1,1])
     }
     sum.dist[i]<-sum(min.dist*slength);
     #sum.dist[i]<-sum(min.dist);
   }
   out<-data.frame(CanP$name1,sum.dist)
   out<-out[order(out[,2]),]
   return(head(out,10))
}
#fitPP(OLL,slength)


CPdemo<-function(p=c(0.95,0.95),la=c(2,4)){
  p_grid<-p; la_grid<-la  #la is relative average ploidy
  foo=Cpoints(p_grid,la_grid)
  par(mfrow=c(2,2))
  plot(cbind(foo$cLL[,,2][1,], foo$cLL[,,1][1,]),xlim=c(0,0.5),ylim=c(-2,2))
  plot(cbind(foo$cLL[,,2][4,], foo$cLL[,,1][4,]),xlim=c(0,0.5),ylim=c(-2,2))
  #getBA(x=data0[i,5],y=data0[i,4],la=4,p=0.92)
}
#CPdemo(p=c(0.95,0.95),la=c(2,4))

CPdemo<-function(p=c(0.8,1),la=c(2,2)){
  p_grid<-p; la_grid<-la  #la is relative average ploidy
  foo=Cpoints(p_grid,la_grid)
  par(mfrow=c(2,2))
  plot(cbind(foo$cLL[,,2][1,], foo$cLL[,,1][1,]),xlim=c(0,0.5),ylim=c(-2,2))
  plot(cbind(foo$cLL[,,2][4,], foo$cLL[,,1][4,]),xlim=c(0,0.5),ylim=c(-2,2))

  data1<-cbind(foo$cLL[,,1][1,], foo$cLL[,,2][1,])
  data2<-cbind(foo$cLL[,,1][4,], foo$cLL[,,2][4,])
  BA<-matrix(nrow=nrow(data1),ncol=2)
  for (i in 1:nrow(data1)) {
    BA[i,]=getBA(x=data1[i,1],y=data1[i,2],la=la[1],p=1)
  }
  plot(BA,xlim=c(0,3),ylim=c(0,5))

  BA<-matrix(nrow=nrow(data1),ncol=2)
  for (i in 1:nrow(data1)) {
    BA[i,]=getBA(x=data2[i,1],y=data2[i,2],la=2,p=1)
  }
  plot(BA,xlim=c(0,3),ylim=c(0,5))
}
#CPdemo()
#This function should be replaced by get.MCpoints(p_grid=0.95,la=2)



################# Testing Data ##################
#data<-read.table("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/CNV_status/CNcalls/yurolTM.0.45.CNV.output.txt",head=T)
data0<-read.table("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/CNV_status/CNcalls/yurolTM.0.45.CNV.output.txt",head=T)
data<-data0[data0$logR_mean>-2&data0$logR_mean<2,]#remove outlier
OLL<-data.matrix(data[,c(5,4)])
slength<-(data$end-data$start)/1000000
##
p_grid<-seq(0.2,1,0.1); la_grid<-seq(1,5,1)  #la is ploidy
CanP<-Cpoints(p_grid,la_grid)
#fitPP(OLL,slength,LAF.wt=c(1,5))
fitPP(OLL,slength,LAF.wt=c(1,4))



##get ploidy estimates of 76 samples
sub<-read.csv("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/sample_path.csv",head=T)[1:76,1:4]
tnames<-as.character(sub$Tumor)
results<-NULL
for (i in 1:length(tnames)) {
    fname<-paste("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/CNV_status/CNcalls/",tnames[i],".0.45.CNV.output.txt",sep="")
    data0<-read.table(fname,head=T)
    data<-data0[data0$logR_mean>-2&data0$logR_mean<2,]
    OLL<-data.matrix(data[,c(5,4)])
    slength<-(data$end-data$start)/1000000
    ##
    p_grid<-seq(0.1,1,0.01); la_grid<-seq(1,5,1)  #la is ploidy
    CanP<-Cpoints(p_grid,la_grid)
    results[i]<-as.character(fitPP(OLL,slength,LAF.wt=c(1,5))[1,1])
}
OUT<-cbind(tnames,results)
save(OUT,file="76ploidy1.Rdata")

sub<-read.csv("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/sample_path.csv",head=T)[1:76,1:4]
tnames<-as.character(sub$Tumor)
results<-NULL
for (i in 1:length(tnames)) {
    fname<-paste("~/Dropbox/xiaoqing_share/exome_CNV_MEGAfolder/CNV_status/CNcalls/",tnames[i],".0.45.CNV.output.txt",sep="")
    data0<-read.table(fname,head=T)
    data<-data0[data0$logR_mean>-2&data0$logR_mean<2,]
    OLL<-data.matrix(data[,c(5,4)])
    slength<-(data$end-data$start)/1000000
    ##
    p_grid<-seq(0.1,1,0.01); la_grid<-seq(1,5,1)  #la is ploidy
    CanP<-Cpoints(p_grid,la_grid)
    results[i]<-as.character(fitPP(OLL,slength,LAF.wt=c(1,4))[1,1])
}
OUT<-cbind(tnames,results)
save(OUT,file="76ploidy2.Rdata")






###
par(mfrow=c(2,2))
BA<-matrix(nrow=nrow(data0),ncol=2)
for (i in 1:nrow(data0)) {
  #i=1
  BA[i,]=getBA(x=data0[i,5],y=data0[i,4],la=4,p=0.92)
}
data1<-cbind(data0,BA)
names(data1)[11:12]<-c("min","maj")
plot(data1$min,data1$maj)
for (i in 1:nrow(data0)) {
  #i=1
  BA[i,]=getBA(x=data0[i,5],y=data0[i,4],la=2,p=0.9)
}
data1<-cbind(data0,BA)
names(data1)[11:12]<-c("min","maj")
plot(data1$min,data1$maj)



################# lab ##################

p_grid<-c(0.95, 0.95); la_grid<-c(2,4)  #la is relative average ploidy
foo=Cpoints(p_grid,la_grid)
par(mfrow=c(2,2))
plot(cbind(foo$cLL[,,2][1,], foo$cLL[,,1][1,]),xlim=c(0,0.5),ylim=c(-2,2))
plot(cbind(foo$cLL[,,2][4,], foo$cLL[,,1][4,]),xlim=c(0,0.5),ylim=c(-2,2))




