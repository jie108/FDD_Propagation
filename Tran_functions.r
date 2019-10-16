#################################
 ################## auxiliary functions: specify transition prob. and get absorbing probabilities
 Tran.prob<-function(FDD.in, state.dir){
## parameter: FDD.in: 27 by 26 FDD prob. matrix: each row is one voxel, each column is one direction; row sum equals to one
## state.dir : state + direction => state correspondence matrix
## return:  tran: 125 by 125 transition probability matrix among the 125 states

tran.in=matrix(0,27,125)  ## transition prob. from the 27 in-states
 for (i in 1:27){
     index.c=state.dir[i,]
     tran.in[i, index.c]=FDD.in[i,]
 }

 rownames(tran.in)=paste("in",1:27,sep="")
 colnames(tran.in)=c(paste("in",1:27,sep=""), paste("out",1:98,sep=""))


 tran.out=matrix(0,98,125) ##transition prob. from the 98 out-states: note they are absorbing states
  for (i in 1:98){
   tran.out[i,i+27]=1

  }

 rownames(tran.out)=paste("out",1:98,sep="")
 colnames(tran.out)=c(paste("in",1:27,sep=""), paste("out",1:98,sep=""))

  tran=rbind(tran.in, tran.out)

  return(tran)
 }

 Tran.abs<-function(tran, max.step=100, thre=1e-6){##get absorbing prob. starting from each state
 ## parameter: tran, 125 by 125 transition prob. matrix
 ## max.step: maximum step used to get absorbing prob.
 ## thre: threshold for convergence

 tran.old=tran

 count=1
 conv=99
 while(count<=max.step && conv>thre){

 tran.new=tran.old%*%tran
 count=count+1
 conv=max(abs(tran.new-tran.old))
 tran.old=tran.new
 }

 tran.abs=tran.new
 return(list("tran.abs"=tran.abs,"conv"=conv))

 }


 ### symmetrize FDD and make the total prob. equals to one:
 FDD.symm<-function(FDD.v, dir.pair){
 ##parameter: FDD.v -- an FDD vector on 26 directions; dir.pair: 13 pair of  directions 
 ##return: symmetrized and normalized FDD vector 
  result=FDD.v
  
  for (i in 1:13){
    index1=dir.pair[i,1]
    index2=dir.pair[i,2]
    prob.c=0.5*(FDD.v[index1]+FDD.v[index2])
    result[index1]= prob.c
    result[index2]= prob.c  
  }
  
  result=result/sum(result)
  
  return(result)
 
 
 
 }
 
 FDD.symm.check<-function(FDD.v,dir.pair){
 ##check symmetry of FDD 
 result=TRUE
  
  for (i in 1:13){
    index1=dir.pair[i,1]
    index2=dir.pair[i,2]
    result=result*(FDD.v[index1]==FDD.v[index2])
  }
  
  return(result)
 
 
 }
 
 FDD.thre<-function(FDD.v, thre=1/26){
 ### set small FDD prob. to zero and renormalize to have total prob. 1. 
 ### This will increase directionality at higher-level super voxels
 ###FDD.v: FDD flow prob. ther: thresholding value, default = 1/26: prob. of a direction when uniform FDD is used (for voxles without coherent flow directions)
 ### return: thresholded and renormalized FDD
 
 result=FDD.v
 result[FDD.v<thre]=0
 if(sum(result)>0)
 result=result/sum(result)
 
 return(result)
 }
 
 
 FDD.thre_new<-function(FDD.v, thre=0.1){ ## threshold all prob. < thre * max
 ### set small FDD prob. to zero and renormalize to have total prob. 1. 
 ### This will increase directionality at higher-level super voxels
 ###FDD.v: FDD flow prob. ther: thresholding value, default = 1/26: prob. of a direction when uniform FDD is used (for voxles without coherent flow directions)
 ### return: thresholded and renormalized FDD
 
 result=FDD.v
 thre.c=thre*max(FDD.v)
 result[FDD.v<thre.c]=0
 result=result/sum(result)
 
 }
 
 FDD.thre_gap<-function(FDD.v, thre=1){ ## threshold at the thre(=k) biggest gap
 ### set small FDD prob. to zero and renormalize to have total prob. 1. 
 ### This will increase directionality at higher-level super voxels
 ###FDD.v: FDD flow prob. ther: thresholding value, default = 1/26: prob. of a direction when uniform FDD is used (for voxles without coherent flow directions)
 ### return: thresholded and renormalized FDD
 
 k=thre
 result=FDD.v
 
 temp=sort(unique(FDD.v))
 temp.d=diff(temp)
 temp.index= order(temp.d, decreasing=TRUE)  ##find the kth biggest gap 
 index.c=temp.index[k]   ##find the kth biggest gap 
 
 thre.c=temp[index.c+1]
 result[FDD.v<thre.c]=0
 result=result/sum(result)
 
 }
 
 #### calculate FDD from a single peak FOD 
 FOD.to.FDD<-function(FOD, grid.0, dir.v.norm){
 ## parameters: FOD: single peak FOD: n by 1 vector 
 ## grid.0:  the corresponding grid points, 3 by n matrix; dir.v.norm: normalized direction vectors (26)
 ## return: FDD: 26 by 1 vector 
 
 ##  group the grid  points to each direction
 group.id=numeric(ncol(grid.0))
for (i in 1:ncol(grid.0)){
 grid.c=grid.0[,i]
 inner.c=dir.v.norm%*%matrix(grid.c)
 group.id[i]=which.max(inner.c)
} 

temp=table(group.id)  
group.prop=temp/sum(temp)  ## proportion of grid points in each group

### FDD.z: FDD of the FOD aligned with z-direction, resized by group size (some of the 26 directions cover more areas than others)
FDD=numeric(26)
for (i in 1:26){
FDD[i]= sum(FOD[group.id==i])/group.prop[i]
}
FDD=FDD/sum(FDD)

return(FDD)
 
 }
 
 
#### calculate sharing probability for a grid
 Prob.share<-function(grid.0, dir.v.norm, k=2, thre=0){
 ## grid.0:  the corresponding grid points, 3 by n matrix; dir.v.norm: normalized direction vectors (26)
 ## k: sharing between the k -closest directions 
 ## thre: cutoff value to thresholding small prob. ; should be no larger than 1/26
 ## return: FDD: 26 by 1 vector 
 
 ## for each grid point, get the sharing scheme to the 26 FDD directions 
 #dist.0=matrix(0,ncol(grid.0), 26)
 grid.prob=matrix(0,ncol(grid.0), 26)
 
for (i in 1:ncol(grid.0)){
 grid.c=grid.0[,i]
 
 inner.c=dir.v.norm%*%matrix(grid.c)
 dist.c=acos(inner.c)
 #dist.0[i,]=dist.c

 index.c=which.min(dist.c)[1]
 dist.min=min(dist.c)
  if (dist.min==0){       ## no ambuguity
   grid.prob[i,index.c]=1  
   }else{
   order.c=order(dist.c)
   grid.prob[i,order.c[1:k]]=1/dist.c[order.c[1:k]]/sum(1/dist.c[order.c[1:k]])
   grid.prob[i,grid.prob[i,]<thre]=0
  
   }
  }

return(grid.prob)
 
 }
 
 
  
 #### calculate FDD from an FOD 
 FOD.to.FDD_new<-function(FOD, grid.0, dir.v.norm, k=2, thre=0){
 ## parameters: FOD: single peak FOD: n by 1 vector 
 ## grid.0:  the corresponding grid points, 3 by n matrix; dir.v.norm: normalized direction vectors (26)
 ## k: sharing between the k -closest directions 
 ## thre: cutoff value to thresholding small prob. ; should be no larger than 1/26
 ## return: FDD: 26 by 1 vector 
 
 ## for each grid point, get the sharing scheme to the 26 FDD directions 
 #dist.0=matrix(0,ncol(grid.0), 26)
 grid.prob=matrix(0,ncol(grid.0), 26)
 
for (i in 1:ncol(grid.0)){
 grid.c=grid.0[,i]
 
 inner.c=dir.v.norm%*%matrix(grid.c)
 dist.c=acos(inner.c)
 #dist.0[i,]=dist.c

 index.c=which.min(dist.c)[1]
 dist.min=min(dist.c)
  if (dist.min==0){       ## no ambuguity
   grid.prob[i,index.c]=1  
   }else{
   order.c=order(dist.c)
   grid.prob[i,order.c[1:k]]=1/dist.c[order.c[1:k]]/sum(1/dist.c[order.c[1:k]])
   grid.prob[i,grid.prob[i,]<thre]=0
  
   }
  }

## FOD to FDD
FOD.u=FOD/sum(FOD)
FDD=t(matrix(FOD.u, 1, ncol(grid.0))%*%grid.prob)


return(FDD)
 
 }
 
 
 ### orthogonal transformation to transform a unit-vector to (0,0,1)  
 Orth.tran<-function(theta,phi){
 ## parameters:  theta: polar angle, [0,pi); phi: azimuthal angle, [0,2*pi) 
 ## return: the   orthogonal matrix to transform (theta, phi) to (0,0,1)  
 
 phi.tran=matrix(0,3,3)            ## in x-y plane, rotate to x-axis
 phi.tran[1:2,1:2]=matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)),2,2)
 phi.tran[3,3]=1
 theta.tran=matrix(0,3,3)          ## in x-z plane, rotate to z-axis
 theta.tran[c(1,3),c(1,3)]=matrix(c(cos(theta), sin(theta),-sin(theta),cos(theta)),2,2)
 theta.tran[2,2]=1
 
 A=theta.tran%*%phi.tran           ##orthorgonal matrix to transform (theta, phi) to (0,0) [i.e.,z-axis]
 
 return(A)
 
 }