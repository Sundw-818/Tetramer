
# give filepath and pdbname, the function pdbread can return the position of CA in the pdbfile

filepath <- getwd()
#filepath <- "/Volumes/MyPassport/科研/四体蛋白质界面残基对预测/renew/data"
pdbname <- "2h8n"

filename <- paste(filepath,"/",pdbname,".pdb",sep="")

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

posall <- pdbread(filename)
chains <- table(posall$X1)
chains <- as.data.frame(chains)
# chains include chain name and the number of residues in this chain

# divide four chains and save
pos1 <- posall[which(posall$X1==chains$Var1[1]),]
pos2 <- posall[which(posall$X1==chains$Var1[2]),]
pos3 <- posall[which(posall$X1==chains$Var1[3]),]
pos4 <- posall[which(posall$X1==chains$Var1[4]),]

posfile1 <- paste(filepath,"/",pdbname,"_",as.character(chains[1,1]),".txt",sep="")
posfile2 <- paste(filepath,"/",pdbname,"_",as.character(chains[2,1]),".txt",sep="")
posfile3 <- paste(filepath,"/",pdbname,"_",as.character(chains[3,1]),".txt",sep="")
posfile4 <- paste(filepath,"/",pdbname,"_",as.character(chains[4,1]),".txt",sep="")
write.table(pos1,file=posfile1,row.names = FALSE,col.names = FALSE, quote=FALSE)
write.table(pos2,file=posfile2,row.names = FALSE,col.names = FALSE, quote=FALSE)
write.table(pos3,file=posfile3,row.names = FALSE,col.names = FALSE, quote=FALSE)
write.table(pos4,file=posfile4,row.names = FALSE,col.names = FALSE, quote=FALSE)

# read the position of 4 chains save before
chain1<-read.table(posfile1)
chain2<-read.table(posfile2)
chain3<-read.table(posfile3)
chain4<-read.table(posfile4)

# use the position of pdb to find the k nearest residues of each residue
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

d1 <- k_nearest(5,chain1)
dd1 <- d1$one
cc1 <- d1$two

fpath1 <- paste(filepath,"/",pdbname,"_",chains[1,1],"_output.txt",sep="")
feature1 <- read.table(fpath1)
m <- length(feature1[,1])
f_nearest1 <- data.frame()
for(i in 1:m){
  f_near1 <- cbind(feature1[dd1[i,1],5:13],feature1[dd1[i,2],5:13],feature1[dd1[i,3],5:13],feature1[dd1[i,4],5:13],feature1[dd1[i,5],5:13])
  f_nearest1 <- rbind(f_nearest1,f_near1)
}
f1 <- cbind(feature1,f_nearest1,cc1)

d2 <- k_nearest(5,chain2)
dd2 <- d2$one
cc2 <- d2$two

fpath2 <- paste(filepath,"/",pdbname,"_",chains[2,1],"_output.txt",sep="")
feature2 <- read.table(fpath2)
m <- length(feature2[,1])
f_nearest2 <- data.frame()
for(i in 1:m){
  f_near2 <- cbind(feature2[dd2[i,1],5:13],feature2[dd2[i,2],5:13],feature2[dd2[i,3],5:13],feature2[dd2[i,4],5:13],feature2[dd2[i,5],5:13])
  f_nearest2 <- rbind(f_nearest2,f_near2)
}
f2 <- cbind(feature2,f_nearest2,cc2)

d3 <- k_nearest(5,chain3)
dd3 <- d3$one
cc3 <- d3$two

fpath3 <- paste(filepath,"/",pdbname,"_",chains[3,1],"_output.txt",sep="")
feature3 <- read.table(fpath3)
m <- length(feature3[,1])
f_nearest3 <- data.frame()
for(i in 1:m){
  f_near3 <- cbind(feature3[dd3[i,1],5:13],feature3[dd3[i,2],5:13],feature3[dd3[i,3],5:13],feature3[dd3[i,4],5:13],feature3[dd3[i,5],5:13])
  f_nearest3 <- rbind(f_nearest3,f_near3)
}
f3 <- cbind(feature3,f_nearest3,cc3)

d4 <- k_nearest(5,chain4)
dd4 <- d4$one
cc4 <- d4$two

fpath4 <- paste(filepath,"/",pdbname,"_",chains[4,1],"_output.txt",sep="")
feature4 <- read.table(fpath4)
m <- length(feature4[,1])
f_nearest4 <- data.frame()
for(i in 1:m){
  f_near4 <- cbind(feature4[dd4[i,1],5:13],feature4[dd4[i,2],5:13],feature4[dd4[i,3],5:13],feature4[dd4[i,4],5:13],feature4[dd4[i,5],5:13])
  f_nearest4 <- rbind(f_nearest4,f_near4)
}
f4 <- cbind(feature4,f_nearest4,cc4)

ffile1 <- paste(filepath,"/",pdbname,"_",chains[1,1],"_feature.txt",sep="")
ffile2 <- paste(filepath,"/",pdbname,"_",chains[2,1],"_feature.txt",sep="")
ffile3 <- paste(filepath,"/",pdbname,"_",chains[3,1],"_feature.txt",sep="")
ffile4 <- paste(filepath,"/",pdbname,"_",chains[4,1],"_feature.txt",sep="")
write.table(f1,file=ffile1,row.names = FALSE,col.names = FALSE, quote=FALSE)
write.table(f2,file=ffile2,row.names = FALSE,col.names = FALSE, quote=FALSE)
write.table(f3,file=ffile3,row.names = FALSE,col.names = FALSE, quote=FALSE)
write.table(f4,file=ffile4,row.names = FALSE,col.names = FALSE, quote=FALSE)