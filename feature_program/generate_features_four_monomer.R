# get the current working directory
filepath <- getwd()
#pdbname <- "2h8n_A"
# read four pdbnames
pdbpath <- paste(filepath,"/data/pdb/",sep="")
setwd(pdbpath)
pdbname <- dir()
pdbname1 <- pdbname[1]
pdbname2 <- pdbname[2]
pdbname3 <- pdbname[3]
pdbname4 <- pdbname[4]

##############################################################################
# two functions
pdbread <- function(filename){
  maxnum <- 200000
  posall <- c()
  protein <- file(filename,open="rt")
  
  for(i in 1:200000){
    lines <- readLines(protein,1)
    if(length(lines)==0) 
      break 
    else if(substring(lines,1,4)=="ATOM" & substring(lines,15,15)=="A"){
      pos <- c(substring(lines,22,22),as.numeric(substring(lines,23,26)),as.numeric(substring(lines,31,38)),as.numeric(substring(lines,39,46)),as.numeric(substring(lines,47,54)))
      posall <- rbind(posall,pos)
    }
  }
  close(protein)
  posall <- data.frame(posall)
  return(posall)
}

k_nearest <- function(k,chain){
  m <- length(chain[,1])
  distance <- matrix(0,m,m)
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      distance[i,j] <- sqrt(sum((chain[i,3:5]-chain[j,3:5])^2))
      distance[j,i] <- distance[i,j]
    }
  }
  
  dd <- t(apply(distance,1,sort,decreasing=F))
  cc <- t(apply(distance,1,order,decreasing=F))
  
  distance_no <- cc[,2:(k+1)]
  chain_distance <- dd[,2:(k+1)]
  return(list(one=distance_no, two=chain_distance))
}
##############################################################################
# monomer1
filename1 <- paste(pdbpath,pdbname1,sep="")
pos1 <- pdbread(filename1)
posfile1 <- paste(filepath,"/",pdbname1,".txt",sep="")
write.table(pos1,file=posfile1,row.names = FALSE,col.names = FALSE, quote=FALSE)

chain1<-read.table(posfile1)
file.remove(posfile1)

d1 <- k_nearest(3,chain1)
dd1 <- d1$one
cc1 <- d1$two

fpath1 <- paste(filepath,"/results/features_by_chain/",substring(pdbname1,1,nchar(pdbname1)-4),"_output.txt",sep="")
feature1 <- read.table(fpath1)
m <- length(feature1[,1])
f_nearest1 <- data.frame()
for(i in 1:m){
  f_near1 <- cbind(feature1[dd1[i,1],5:13],feature1[dd1[i,2],5:13],feature1[dd1[i,3],5:13])
  f_nearest1 <- rbind(f_nearest1,f_near1)
}
f1 <- cbind(feature1,f_nearest1,cc1)

ffile1 <- paste(filepath,"/",substring(pdbname1,1,nchar(pdbname1)-4),"_feature.txt",sep="")
write.table(f1,file=ffile1,row.names = FALSE,col.names = FALSE, quote=FALSE)
##############################################################################
# monomer2
filename2 <- paste(pdbpath,pdbname2,sep="")
pos2 <- pdbread(filename2)
posfile2 <- paste(filepath,"/",pdbname2,".txt",sep="")
write.table(pos2,file=posfile2,row.names = FALSE,col.names = FALSE, quote=FALSE)

chain2<-read.table(posfile2)
file.remove(posfile2)

d2 <- k_nearest(3,chain2)
dd2 <- d2$one
cc2 <- d2$two

fpath2 <- paste(filepath,"/results/features_by_chain/",substring(pdbname2,1,nchar(pdbname2)-4),"_output.txt",sep="")
feature2 <- read.table(fpath2)
m <- length(feature2[,1])
f_nearest2 <- data.frame()
for(i in 1:m){
  f_near2 <- cbind(feature2[dd2[i,1],5:13],feature2[dd2[i,2],5:13],feature2[dd2[i,3],5:13])
  f_nearest2 <- rbind(f_nearest2,f_near2)
}
f2 <- cbind(feature2,f_nearest2,cc2)

ffile2 <- paste(filepath,"/",substring(pdbname2,1,nchar(pdbname2)-4),"_feature.txt",sep="")
write.table(f2,file=ffile2,row.names = FALSE,col.names = FALSE, quote=FALSE)
##############################################################################
# monomer3
filename3 <- paste(pdbpath,pdbname3,sep="")
pos3 <- pdbread(filename3)
posfile3 <- paste(filepath,"/",pdbname3,".txt",sep="")
write.table(pos3,file=posfile3,row.names = FALSE,col.names = FALSE, quote=FALSE)

chain3<-read.table(posfile3)
file.remove(posfile3)

d3 <- k_nearest(3,chain3)
dd3 <- d3$one
cc3 <- d3$two

fpath3 <- paste(filepath,"/results/features_by_chain/",substring(pdbname3,1,nchar(pdbname3)-4),"_output.txt",sep="")
feature3 <- read.table(fpath3)
m <- length(feature3[,1])
f_nearest3 <- data.frame()
for(i in 1:m){
  f_near3 <- cbind(feature3[dd3[i,1],5:13],feature3[dd3[i,2],5:13],feature3[dd3[i,3],5:13])
  f_nearest3 <- rbind(f_nearest3,f_near3)
}
f3 <- cbind(feature3,f_nearest3,cc3)

ffile3 <- paste(filepath,"/",substring(pdbname3,1,nchar(pdbname3)-4),"_feature.txt",sep="")
write.table(f3,file=ffile3,row.names = FALSE,col.names = FALSE, quote=FALSE)
##############################################################################
# monomer4
filename4 <- paste(pdbpath,pdbname4,sep="")
pos4 <- pdbread(filename4)
posfile4 <- paste(filepath,"/",pdbname4,".txt",sep="")
write.table(pos4,file=posfile4,row.names = FALSE,col.names = FALSE, quote=FALSE)

chain4<-read.table(posfile4)
file.remove(posfile4)

d4 <- k_nearest(3,chain4)
dd4 <- d4$one
cc4 <- d4$two

fpath4 <- paste(filepath,"/results/features_by_chain/",substring(pdbname4,1,nchar(pdbname4)-4),"_output.txt",sep="")
feature4 <- read.table(fpath4)
m <- length(feature4[,1])
f_nearest4 <- data.frame()
for(i in 1:m){
  f_near4 <- cbind(feature4[dd4[i,1],5:13],feature4[dd4[i,2],5:13],feature4[dd4[i,3],5:13])
  f_nearest4 <- rbind(f_nearest4,f_near4)
}
f4 <- cbind(feature4,f_nearest4,cc4)

ffile4 <- paste(filepath,"/",substring(pdbname4,1,nchar(pdbname4)-4),"_feature.txt",sep="")
write.table(f4,file=ffile4,row.names = FALSE,col.names = FALSE, quote=FALSE)




