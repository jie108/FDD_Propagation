

### process FOD estimates from matlab
## these estiamtes have already been aligned such that, the true fibers across all replicates have the same direction
## we want to rotate them to a specific direction, and then calcualte the FDD
## these (noisy) FDD will then be used to build a level 2 supervoxel
### no prob. sharing


load("state.Rdata")
load("super1.Rdata")
source("../Tran_functions.r")

##true FDD
FDD1=FDD.in[1,,1]  ## with crossing
FDD2=FDD.in[2,,1]  ## x22
FDD3=FDD.in[8,,1]  ## x122

FDD.x22.tr=FDD2    ##x-22
FDD.x122.tr=FDD3  ## x-122
FDD.xy22.tr=FDD1  ##crossing

###
library(R.matlab)
#fod.est=readMat("1fib_N81_ave.MAT")

fod.est=readMat("1fib_N41_ave.MAT")
fod.est=fod.est$fod.classo.tr  ## get classo estimate
dim(fod.est)
N_rep=dim(fod.est)[1]

grid.ea=readMat("equal_angle_grid.MAT")
grid.ea=grid.ea[[1]]          ## 3 by 2562
dim(grid.ea)

############
 ###   22.5 degree with x
theta=pi/2  ## polar angle, [0,pi)
 phi=pi/8  ## azimuthal angle, [0,2*pi)
 x0=sin(theta)*cos(phi)
 y0=sin(theta)*sin(phi)
 z0=cos(theta)
 dir.x22=c(x0,y0,z0)

 A=Orth.tran(theta,phi)

 ### apply t(A) to the grid
 grid.0=t(A)%*%grid.ea

######### FDD
FDD.x22=matrix(0,N_rep,26)        ##22.5 with x-
colnames(FDD.x22)=paste("dir",1:26)
#fac.s=5  ##stretching factor to sharpen the fod estimation
fac.s=1
for (i in 1:N_rep){

  FDD.x22[i,]=FOD.to.FDD_new(fod.est[i,]^fac.s, grid.0, dir.v.norm,k=2)

  FDD.x22[i,]=FDD.thre_gap(FDD.x22[i,], thre=1)  ##thresholding at the biggest gap 
  FDD.x22[i,]=FDD.symm(FDD.x22[i,], dir.pair)        ##check symmetry: true
  ## FDD.symm.check(FDD.x22,dir.pair)              ##check symmetry: true
}

 index.dir.z0=which(dir.v[,3]==0)
  apply(FDD.x22,1,sum)

#### mean and sd
FDD.x22.ave=apply(FDD.x22,2,mean)
FDD.x22.bias=FDD.x22.tr-FDD.x22.ave
FDD.x22.sd=apply(FDD.x22,2,sd)

round(cbind(FDD.x22.bias, FDD.x22.sd)*100,2)[index.dir.z0,] ##bias, sd
cbind(dir.v, round(FDD.x22.ave*100,2),round(FDD.x22.tr*100,2))   ##compare with truth

 temp=cbind(abs(FDD.x22.bias), FDD.x22.sd, FDD.x22.tr)
temp1=temp[FDD.x22.tr>0,]
temp2=temp[FDD.x22.tr==0,]

round(apply(temp1,2,mean),3)    ##0.01, 0.141, 0.25          
round(apply(temp1/temp1[,3],2, mean),3)    ##0.041, 0.563
round(apply(temp2,2,mean),3)     ##0.002, 0.007



 ######### 122.5 with x-
 ###   122.5 degree with x
theta=pi/2  ## polar angle, [0,pi)
 phi=pi/8+pi/2  ## azimuthal angle, [0,2*pi)
 x0=sin(theta)*cos(phi)
 y0=sin(theta)*sin(phi)
 z0=cos(theta)
 dir.x122=c(x0,y0,z0)

 A=Orth.tran(theta,phi)

 ### apply t(A) to the grid
 grid.0=t(A)%*%grid.ea


 ##
 FDD.x122=matrix(0,N_rep,26)        ##122.5 with x-
colnames(FDD.x122)=paste("dir",1:26)

for (i in 1:N_rep){

  FDD.x122[i,]=FOD.to.FDD_new(fod.est[i,], grid.0, dir.v.norm,k=2)
  FDD.x122[i,]=FDD.thre_gap(FDD.x122[i,], thre=1)  ##thresholding at the biggest gap 
  FDD.x122[i,]=FDD.symm(FDD.x122[i,], dir.pair)        ##check symmetry: true
  ## FDD.symm.check(FDD.x22,dir.pair)              ##check symmetry: true
}

  apply(FDD.x122,1,sum)

## mean and sd
FDD.x122.ave=apply(FDD.x122,2,mean)

FDD.x122.bias=FDD.x122.tr-FDD.x122.ave
FDD.x122.sd=apply(FDD.x122,2,sd)

round(cbind(FDD.x122.bias, FDD.x122.sd)*100,2)[index.dir.z0,]
cbind(dir.v, round(FDD.x122.ave*100,2),round(FDD.x122.tr*100,2))



####
####   x22 and x122 crossing
fod.est.xy=readMat("sep90_N41_ave.MAT")
fod.est.xy=fod.est.xy$fod.classo.tr          ##classo estimates
dim(fod.est.xy)
N_rep=dim(fod.est.xy)[1]

 ###   22.5 degree with x
theta=pi/2  ## polar angle, [0,pi)
 phi=pi/8  ## azimuthal angle, [0,2*pi)
 x0=sin(theta)*cos(phi)
 y0=sin(theta)*sin(phi)
 z0=cos(theta)
 dir.x22=c(x0,y0,z0)

 temp1=matrix(c(1,0,0,0,0,1,0,1,0),3,3)
 temp2=matrix(c(cos(phi),sin(phi),0,-sin(phi),cos(phi),0,0,0,1),3,3)
 A=t(temp2%*%temp1)
 A%*%dir.x22
 A%*%dir.x122
 
 ### apply t(A) to the grid
 grid.0=t(A)%*%grid.ea


 ##
 FDD.xy22=matrix(0,N_rep,26)        ##22 and 122.5 with x-
colnames(FDD.xy22)=paste("dir",1:26)

for (i in 1:N_rep){

  FDD.xy22[i,]=FOD.to.FDD_new(fod.est.xy[i,], grid.0, dir.v.norm,k=2)
 
   FDD.xy22[i,]=FDD.thre_gap(FDD.xy22[i,], thre=1)  ##thresholding at the biggest gap 
  FDD.xy22[i,]=FDD.symm(FDD.xy22[i,], dir.pair)        ##check symmetry: true
  ## FDD.symm.check(FDD.x22,dir.pair)              ##check symmetry: true
}

 
apply(FDD.xy22,1,sum)
  
## mean and sd 
FDD.xy22.ave=apply(FDD.xy22,2,mean)
FDD.xy22.bias=FDD.xy22.tr-FDD.xy22.ave
FDD.xy22.sd=apply(FDD.xy22,2,sd)

round(cbind(FDD.xy22.bias, FDD.xy22.sd, FDD.xy22.tr)*100,2)[index.dir.z0,]
cbind(dir.v, round(FDD.xy22.ave*100,2),round(FDD.xy22.tr*100,2))


temp=cbind(abs(FDD.xy22.bias), FDD.xy22.sd, FDD.xy22.tr)
temp1=temp[FDD.xy22.tr>0,]
temp2=temp[FDD.xy22.tr==0,]

round(apply(temp1,2,mean),3)    ##0.015, 0.135, 0.125          
round(apply(temp1/temp1[,3],2, mean),3)    ##0.122, 1.08
round(apply(temp2,2,mean),3)     ##0.005, 0.012

########
save(FDD.x22, FDD.x122, FDD.xy22,file="FOD.x22.N41.Rdata")
  
  
  ###############results: FDD bias and sd
cbind(dir.v, round(FDD.x22.ave*100,2),round(FDD.x22.tr*100,2))   ##one fiber, x-22
cbind(dir.v, round(FDD.x122.ave*100,2),round(FDD.x122.tr*100,2)) ##one fiber, x-122
cbind(dir.v, round(FDD.xy22.ave*100,2),round(FDD.xy22.tr*100,2)) ## two fibers 

round(cbind(FDD.x22.bias, FDD.x22.sd)*100,2)[index.dir.z0,]
round(cbind(FDD.x122.bias, FDD.x122.sd)*100,2)[index.dir.z0,]
round(cbind(FDD.xy22.bias, FDD.xy22.sd)*100,2)[index.dir.z0,]

 