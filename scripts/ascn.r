
normalvcf = "yurolN.vcf"
tumorvcf = "yurolTM.vcf"
samplename = "yurol"

####### Data Preprocessing

mydata = scan(normalvcf, what="character")
N = length(mydata)
id = which(mydata[1:1000]=="#CHROM")
mymatrix0 = matrix(mydata[id:N], (N-id+1)/10, byrow=T)
N0 = dim(mymatrix0)[1]-1
mymatrix = mymatrix0[2:(N0+1),]
colnames(mymatrix) = mymatrix0[1,]

N.matrix = mymatrix

mydata = scan(tumorvcf, what="character")
N = length(mydata)
id = which(mydata[1:1000]=="#CHROM")
mymatrix0 = matrix(mydata[id:N], (N-id+1)/10, byrow=T)
N0 = dim(mymatrix0)[1]-1
mymatrix = mymatrix0[2:(N0+1),]
colnames(mymatrix) = mymatrix0[1,]

TM.matrix = mymatrix

N.pass = which(N.matrix[,7]=="PASS")
TM.pass = which(TM.matrix[,7]=="PASS")

N.matrix.pass = N.matrix[N.pass,]
TM.matrix.pass = TM.matrix[TM.pass,]


### match Normal and Tumor positions

index1 = index2 = c()

for (i in 1:22){
  name = paste("chr",i,sep="")
  id1 = which(N.matrix.pass[,1]==name)
  id2 = which(TM.matrix.pass[,1]==name)
  pos1 = N.matrix.pass[id1,2]
  pos2 = TM.matrix.pass[id2,2]
  temp = match(pos1, pos2)
  NA.id = which(is.na(temp))
  match1 = (1:length(id1))[-NA.id]
  match2 = temp[-NA.id]
  index1 = c(index1, id1[match1])
  index2 = c(index2, id2[match2])
}

## check whether the ref and alt alleles are the same
ref1 = N.matrix.pass[id1,4]
alt1 = N.matrix.pass[id1,5]
ref2 = TM.matrix.pass[id2,4]
alt2 = TM.matrix.pass[id2,5]
for (j in 1:length(match1)){
  if (ref1[match1[j]]!=ref2[match2[j]] || alt1[match1[j]]!=alt2[match2[j]]) print(j)
}


N.matrix.match = N.matrix.pass[index1,]
TM.matrix.match = TM.matrix.pass[index2,]

N0 = dim(N.matrix.match)[1]

N.gt = matrix(0,N0,2)
for (i in 1:N0){
  temp = N.matrix.match[i,10]
  N.gt[i,] = as.numeric(strsplit(strsplit(temp,":")[[1]][1],"/")[[1]])
}

dif = N.gt[,2] - N.gt[,1]
het = which(dif>0)

N.matrix.het = N.matrix.match[het,]
TM.matrix.het = TM.matrix.match[het,]

# ref.dif = which(N.matrix.het[,4]!=TM.matrix.het[,4])
# alt.dif = which(N.matrix.het[,5]!=TM.matrix.het[,5])

mydata = cbind(N.matrix.het[,c(1,2,4)], TM.matrix.het[,4], N.matrix.het[,5], TM.matrix.het[,5], N.matrix.het[,10], TM.matrix.het[,10])

colnames(mydata) = c("chr", "pos", "ref.N", "ref.TM", "alt.N", "alt.TM", "N", "TM")

mydata = data.frame(mydata)

ref.dif = which(as.character(mydata$ref.N) != as.character(mydata$ref.TM))
alt.dif = which(as.character(mydata$alt.N) != as.character(mydata$alt.TM))  # ref.dif turns out to be a subset of alt.dif

## remove ref.dif, alt.dif
mydata1 = mydata[-unique(c(ref.dif,alt.dif)),]
N1 = dim(mydata1)[1]
N.gt = TM.gt = matrix(0,N1,2)
for (i in 1:N1){
  N.gt[i,] = as.numeric(strsplit(strsplit(as.character(mydata1[i,7]),":")[[1]][1],"/")[[1]])
  TM.gt[i,] = as.numeric(strsplit(strsplit(as.character(mydata1[i,8]),":")[[1]][1],"/")[[1]])
}

## get reads
reads = matrix(0,N1,4)

for (i in 1:N1){
  if (N.gt[i,2]==2 && N.gt[i,1]==1){
    temp1 = as.numeric(strsplit(strsplit(as.character(mydata1[i,7]),":")[[1]][2],",")[[1]])
    temp2 = as.numeric(strsplit(strsplit(as.character(mydata1[i,8]),":")[[1]][2],",")[[1]])
    reads[i,] = c(temp1, temp2)
  }else if (N.gt[i,2]==2 && N.gt[i,1]==0){
    temp1 = as.numeric(strsplit(strsplit(as.character(mydata1[i,7]),":")[[1]][2],",")[[1]])
    temp2 = as.numeric(strsplit(strsplit(as.character(mydata1[i,8]),":")[[1]][2],",")[[1]])
    temp3 = strsplit(as.character(mydata1[i,7]),":")[[1]]
    temp4 = strsplit(as.character(mydata1[i,8]),":")[[1]]
    reads[i,] = c(as.numeric(temp3[7]), temp1[2], as.numeric(temp4[7]), temp2[2])
  }else{
    temp1 = strsplit(as.character(mydata1[i,7]),":")[[1]]
    temp2 = strsplit(as.character(mydata1[i,8]),":")[[1]]
    reads[i,] = as.numeric(c(temp1[7],temp1[2],temp2[7],temp2[2]))
  }
}

reads = cbind(mydata1[1:4],reads)
colnames(reads) = c("chr","pos","ref","alt","AN","BN","AT","BT")
reads = data.frame(reads)

AN = reads$AN
BN = reads$BN
AT = reads$AT
BT = reads$BT
pos = reads$pos

output = cbind(reads, round(AN/(AN+BN),digits=3), round(AT/(AT+BT),digits=3))
colnames(output) = c("chr","pos","ref","alt","AN","BN","AT","BT", "Afreq.N", "Afreq.TM")
write.table(output, file=paste(samplename,"_summary.txt",sep=""), quote=F, row.names=F, sep="\t")

#### analyzing using falcon
library(falcon)
rdep = median(AT+BT)/median(AN+BN)

hardthres = function(v, low=0.9, high=1.1){
  n = length(v)
  for (i in 1:n){ if (v[i]>low && v[i]<high) v[i] = 1 }
  v
}

threshold = 0.20

for (i in 1:22){
  name = paste("chr",i,sep="")
  ids = which(reads$chr==name)
  cn = getASCN(AT[ids],BT[ids],AN[ids],BN[ids],rdep=rdep)
  save(cn,file=paste(samplename,"_chr",i,".Rdata",sep=""))
  # ascn = hardthres(cn$ascn,1-threshold,1+threshold)
  png(paste(samplename,"_chr",i,".png",sep=""),width=1200,height=750,pointsize=23)
  view(AT[ids],BT[ids],AN[ids],BN[ids],cn$tauhat,cn$ascn,pos[ids],rdep=rdep)
  dev.off()
}


result = c()
for (i in 1:22){
  name = paste("chr",i,sep="")
  ids = which(reads$chr==name)
  load(paste(samplename,"_chr",i,".Rdata",sep=""))
  tauhat0 = c(1,cn$tauhat,length(ids)+1)
  m = length(tauhat0)
  start = tauhat0[1:(m-1)]
  end = tauhat0[2:m]-1
  temp = cbind(rep(i,m-1),ids[start], ids[end], pos[ids[start]], pos[ids[end]], round(t(cn$ascn),digits=3))
  result = rbind(result, temp)
}
removeid = c()
for (i in 1:dim(result)[1]){
  a = result[i,6]
  b = result[i,7]
  if (abs(a-1)<0.15 && abs(b-1)<0.15){
    removeid = c(removeid,i)
  }
}
colnames(result) = c("chr", "SNP.start", "SNP.end", "pos.start", "pos.end", "minor.CN", "major.CN")
write.table(result[-removeid,], file=paste(samplename,"_ascn_table.txt",sep=""), quote=F, row.names=F, sep="\t")

