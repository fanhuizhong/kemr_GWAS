args <- commandArgs(TRUE)
gene = args[1]
result = args[2]
c = 0
phenotype = '/home/wuq/projects/kmerEvol/data_version3/chr_number.csv'
phe = as.matrix(read.csv('/home/wuq/projects/kmerEvol/data_version3/chr_number.csv',header=F))
distance = as.matrix(read.csv('/home/wuq/projects/kmerEvol/data_version3/distance_chr.csv',header=F))
distance = distance[,-1]
distance =matrix(as.numeric(distance),dim(distance),dimnames=dimnames(distance))
y = as.numeric(phe[,2])

a = as.matrix(read.csv(gene,header=T))
a = t(a)
aa = cbind(y,a)
bb = c(48,98,89,102,84,62,42,100,32,31,94,63,93,11,13,92,70,50,66,19,57,81,79,28,37,95,61,8,59,33,25,53,86,23,68,55,54,14,101,5,21,16,47,73,4,10,34,75,41,74,27,76,56,72,15,88,96,78,20,52,35,40,18,77,6,49,99,36,44,91,39,1,87,71,58,17,46,12,9,60,30,7,24,65,90,51,22,29,45,3,69,80,82,26,38,67,97,64,2,85,43,83)
hh =cbind(bb,aa)
aa = hh[order(hh[,1]),]
aa = aa[,-1]
aa = aa[-c(34,45,56,61,65,70),]
mm = aa[,1]
nn = cbind(mm,distance)
dat=as.data.frame(nn)
wy=glm(mm~distance,data=dat)
ys=matrix(wy$residuals,,1)
aa = cbind(ys,aa[,c(2:ncol(aa))])

dd = c()
for (i in 1:nrow(aa)){
    bb = aa[i,]
    aa1 = subset(bb,bb!=0&bb!=-1)
    if (length(aa1)/length(bb)<0.90){
    dd = rbind(dd,i)
    }
    }
if (length(dd)==0){
ee = aa
}else{
ee = as.matrix(aa[-dd,])
}
ff = c()
for (ii in 1:nrow(ee)){
    bb1 = ee[ii,-1]
    bb1[which(bb1==0)] = NA
    bb1[which(bb1==-1)] = NA
    aa2 = subset(bb1,bb1!=0&bb1!=-1)
    bb1[is.na(bb1)] = mean(aa2)
    ff = rbind(ff,c(ee[ii,1],bb1))  
    }
rownames(ff) = rownames(ee)
yy1 = ff[,1]
m = ff[,-1]
pca=prcomp(m)
dcor=cor(m)
deig=eigen(dcor)

eSNP=predict(pca)
eSNP=matrix(as.numeric(eSNP),dim(eSNP),dimnames=dimnames(eSNP))
eSNP=as.matrix(eSNP)
p = matrix(0,1,ncol(eSNP))
for(l in 1:ncol(eSNP)){
mod <- glm(yy1~ eSNP[,l], family=gaussian)
p[l] <- summary(mod)$coefficients[2,4]
}
z=-2*sum(log(p))
pz=1-pchisq(z,2*ncol(eSNP))
st=cbind(gene,z,pz)
write.table(st,result,quote=F,row.names=F,col.names=F,sep="\t")
