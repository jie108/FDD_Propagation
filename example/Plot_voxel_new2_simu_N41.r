### simulation based on super 2 design
## input FDD are from estimated FOD based on N=81 directions and bvalue=1000, SNR=20, stored in FOD.x22.Rdata


load("state.Rdata")
load("super2.Rdata")
load("FOD.x22.N41.Rdata")
source("../Tran_functions.r")
dim(FDD.x22)
dim(FDD.x122)
dim(FDD.xy22)
rep.num=dim(FDD.x22)[1]

############## Part I: plot in x-y plane, one level 2 supervoxel
## x-y coordinates of the center of each voxel
## x-y coordinates of the center of each voxel
coor.xy<-matrix(0,81,2)
count=1
 for (i in -4:4){
  for (j in -4:4){
  coor.xy[count,]=c(i,j)
  count=count+1
  }
 }

 coor.bound<-matrix(0,81,4)  ## x-limits and y-limits 

count=1
 for (i in -4:4){
  for (j in -4:4){
  coor.xy[count,]=c(i,j)
  coor.bound[count, c(1,2)]=c(i-0.5,i+0.5)  ## x-limit
  coor.bound[count, c(3,4)]=c(j-0.5,j+0.5)  ##y-limit
  
  count=count+1
  }
 }
 
 
###############
##plot and center of the voxels
 #dev.new(width=5, height=5)

 coor.s1=matrix(0,9,2)  ##coordinates of the nine level 1 super voxels
 rownames(coor.s1)=paste("sup",1:9)
 index.s1=numeric(9)  ## index of centering voxels of level 1 super voxels

count=1
count.s1=1
for (i in -4:4){
  for (j in -4:4){
    col.c=3
     if (is.element(i,c(-3,0,3)) && is.element(j,c(-3,0,3))){
    col.c=2 
    coor.s1[count.s1,]=coor.xy[count,]
    index.s1[count.s1]=count
     count.s1=count.s1+1
     } 
    #oints(coor.xy[count,1],coor.xy[count,2], pch=16, cex=1, col=col.c)
    count=count+1
   }
 }  
 
### indices of level 0 voxel in each level 1 super voxel
 index.01<-matrix(0,9,9)
 rownames(index.01)=paste("sup",1:9)
 colnames(index.01)=paste("vol",1:9)
 

 for (i in 1:9){
  index.s1.c=index.s1[i]
  coor.c=coor.xy[index.s1.c,]
  temp=which(abs(coor.xy[,1]-coor.c[1])<=1&abs(coor.xy[,2]-coor.c[2])<=1)
  index.01[i,]=temp   
   }



 #############################
dim(fib.ind)
 
 ########## assign x-y plane FOD on a grid on the circle
fdd.u=numeric(26)+1/26
dir.v.t=rbind(c(0,0,0),dir.v)
  
result=NULL

N_rep=100 ##100 replicates 

for (repl in 1:N_rep){### 
print(repl)
set.seed(repl*100+51)  ##set random number seed for randomly choosing FDDs 

FDD.xy.c=array(0,dim=c(9,26,9))  ## build FDD for each of the nine voxels on x-y plane
    
  for (i in 1:9){      ##supvoxle i
   for (j in 1:9){    ##voxel j
   index.c=index.01[i,j]
   fib.c=  fib.ind[index.c,]
   
   FDD.c=numeric(26)
   if(fib.c[4]==1){ ## crossing fibers chhose from FDD.xy22
   temp=sample(1:rep.num,1)
   FDD.c=FDD.xy22[temp,]
   }else{
    if(fib.c[2]==1){##one fiber, choose from FDD.x22
    temp=sample(1:rep.num,1)
    FDD.c=FDD.x22[temp,]
     }
     
     if(fib.c[3]==1){##one fiber, choose from FDD.x122
     temp=sample(1:rep.num,1)
    FDD.c=FDD.x122[temp,]
      } 
    }##end else 
   
    if(sum(FDD.c)==0&&(is.element(index.c, index.s1))){## no fiber, but a centering voxel 
    
   FDD.c=fdd.u ##use uniform FDD 
    }
  
   FDD.xy.c[j,,i]=FDD.c   
   #print(FDD.symm.check(FDD.c,dir.pair))   
   
        
   } ## end of i 
  } ## end of j
  
  
 #apply(FDD.xy.c, c(1,3),sum)                      ##some are zero: voxels with no fiber  
 #index.dir.z0=which(dir.v[,3]==0)    ##check whether all flow is within x-y directions, no z-direction flow , yes
 #FDD.xy.temp=FDD.xy.c[,index.dir.z0,]
 #apply(FDD.xy.temp,c(1,3),sum)
 

############### set FDD  in 3D: for all three z slides, use the same flow pattern

FDD.in.c=array(0,dim=c(27,26,9))  ## FDD for each of the 27 voxels in 3D: for all three z slides, use the same flow pattern
 

  for (i in 1:9){  ##super voxel i
   for (j in 1:9){ ## voxel j
     coor.c= coor.xy[index.01[i,j],] ## current voxel
     center.c=coor.xy[index.s1[i],]  ##current center
     # points(coor.c[1],coor.c[2], col=sum(FDD.xy[j,,i])+1)
     
     index.c=which((dir.v.t[,1]==(coor.c[1]-center.c[1]))&(dir.v.t[,2]==(coor.c[2]-center.c[2])))
    
      for (k in index.c){
      FDD.in.c[k,,i]=FDD.xy.c[j,,i]
        #  print(FDD.symm.check(FDD.c,dir.pair))   
      

      }## end of k
    
   } ## end of i 
  } ## end of j
 
# apply(FDD.in.c, c(1,3),sum) 
 
 
 #############
 ######### get FDD for level 1 super voxel
 FDD.out1.c=matrix(0,26,9)

 rownames(FDD.out1.c)=paste("dir",1:26)
 colnames(FDD.out1.c)=paste("sup",1:9)
 
 prob.out1.c=matrix(0,98,9)                  ##98 absorbing state prob
 rownames(prob.out1.c)=paste("out-state",1:98)
 colnames(prob.out1.c)=paste("sup",1:9)
 
 
 for (k in 1:9){
## specify transition matrix: use state.dir and FDD.in
FDD.in.cc=FDD.in.c[,,k]
tran=Tran.prob(FDD.in.cc, state.dir)

 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 tran.abs=temp[[1]]
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 #all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]
 
 ###renormlize to sum 1: due to the voxels with no fibers
  prob.abs.out=prob.abs.out/sum(prob.abs.out)
 prob.out1.c[,k]=prob.abs.out
 #cbind(state.out, round(prob.abs.out*100,2))

 
### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton

FDD.out.cc=matrix(prob.abs.out,1,98)%*%group.out
FDD.out.cc=matrix(FDD.out.cc,26,1)
FDD.out1.c[,k]=FDD.out.cc
 # print(FDD.symm.check(FDD.out.cc,dir.pair)) ##not symmetric anymore due to voxels with zero fdd  

 } ### end of k loop

###################################################################
################################################################### 
 ###################
# apply(FDD.out1.c,2,sum)
# FDD.out1.temp=FDD.out1.c[index.dir.z0,]            ##check whether no z-direction flow
# apply(FDD.out1.temp,2,sum)
# cbind(dir.v, round(FDD.out1.c*100,2))[index.dir.z0,]
 ### 
 
 #######################
 # build level 2 voxel FDD
 ## FDD.in for all 27 level 1 sv
  dimnam=NULL
dimnam[[1]]= paste("vox",1:27)
dimnam[[3]]=paste("sup",1)
FDD.in1.c=array(0,dim=c(27,26,1), dimnames=dimnam)  ## FDD for each of the 27 voxels in 3D: for all three z slides, use the same flow pattern
 

  for (i in 1:1){  ##level 2 super voxel i
   for (j in 1:9){ ## level 1 super voxel j
     coor.c= coor.xy[index.s1[j],] ## current voxel
     
     index.c=which(dir.v.t[,1]==(coor.c[1]/3)&dir.v.t[,2]==(coor.c[2]/3))
    
      for (k in index.c){
      FDD.in1.c[k,,i]=FDD.out1.c[,j]
          
      }## end of k
    
   } ## end of i 
  } ## end of j
 
 apply(FDD.in1.c, c(1,3),sum) 
apply(FDD.in1.c[,index.dir.z0,], 1,sum)   ##no z-direction flow
  
 ####
 ######### get FDD for level 2 super voxel
 FDD.out2.c=matrix(0,26,1)

 rownames(FDD.out2.c)=paste("dir",1:26)
 colnames(FDD.out2.c)=paste("sup",1:1)
 
 prob.out2.c=matrix(0,98,1)                  ##98 absorbing state prob
 rownames(prob.out2.c)=paste("out-state",1:98)
 colnames(prob.out2.c)=paste("sup",1:1)
 
 
 for (k in 1:1){
## specify transition matrix: use state.dir and FDD.in
FDD.in.cc=FDD.in1.c[,,k]
tran=Tran.prob(FDD.in.cc, state.dir)

 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 tran.abs=temp[[1]]
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 #all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]
 
 ###renormlize to sum 1: due to the voxels with no fibers
  prob.abs.out=prob.abs.out/sum(prob.abs.out)
 prob.out2.c[,k]=prob.abs.out
 #cbind(state.out, round(prob.abs.out*100,2))

 
### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton

FDD.out.cc=matrix(prob.abs.out,1,98)%*%group.out
FDD.out.cc=matrix(FDD.out.cc,26,1)
FDD.out2.c[,k]=FDD.out.cc
  #print(FDD.symm.check(FDD.out.cc,dir.pair)) ##not symmetric anymore   

 } ### end of k loop
 
 ###################
 #apply(FDD.out2.c,2,sum)
 #FDD.out2.temp=FDD.out2.c[index.dir.z0,]            ##check whether no z-direction flow
 #sum(FDD.out2.temp)
 #cbind(dir.v, round(FDD.out2.c*100,2))[index.dir.z0,]
 
 ###
#dir.out.z0=which(state.out[,3]==0)
#prob.out2.temp=prob.out2.c[dir.out.z0,]
#apply(prob.out2.c,2,sum)
#sum(prob.out2.temp)

#cbind(state.out, round(prob.out2.c*100,2))[dir.out.z0,]
 
 ###record result.
 result.c=list("FDD.in0"=FDD.in.c, "FDD.out1"=FDD.out1.c, "prob.out1"=prob.out1.c, "FDD.in1"=FDD.in1.c, "FDD.out2"=FDD.out2.c, "prob.out2"=prob.out2.c)
  result[[repl]]=result.c
 } ##end of rep-loop
 
 
###################### 
##########compare with truth
FDD.out2.tr=FDD.out2
 prob.out2.tr=prob.out2
cbind(dir.v, round(FDD.out2*100,2))[index.dir.z0,]
cbind(state.out, round(prob.out2*100,2))[dir.out.z0,]
     
FDD.out2.re=matrix(0,N_rep,26)  ##estiamted results
prob.out2.re=matrix(0,N_rep,98)
prob.out2.thre=matrix(0,N_rep,98) ##apply thresholding 
 for (repl in 1:N_rep){
 result.c=result[[repl]]
 FDD.out2.re[repl,]=result.c$FDD.out2
 prob.out2.re[repl,]=result.c$prob.out2
 prob.out2.thre[repl,]=FDD.thre(result.c$prob.out2, thre=1/98)
 
 }
  
FDD.out2.ave=apply(FDD.out2.re,2, mean)      
prob.out2.ave=apply(prob.out2.re,2, mean)      
prob.out2.thre.ave=apply(prob.out2.thre,2, mean)      

cbind(dir.v, round(FDD.out2.ave*100,2))[index.dir.z0,]
cbind(state.out, round(prob.out2.ave*100,2))[dir.out.z0,]
cbind(state.out, round(prob.out2.thre.ave*100,2))[dir.out.z0,]
    
   
FDD.out2.bias=FDD.out2.tr-FDD.out2.ave
prob.out2.bias=prob.out2.tr-prob.out2.ave
prob.out2.thre.bias=prob.out2.tr-prob.out2.thre.ave


FDD.out2.sd=apply(FDD.out2.re, 2,sd)
prob.out2.sd=apply(prob.out2.re, 2,sd)
prob.out2.thre.sd=apply(prob.out2.thre, 2,sd)
  

round(cbind(FDD.out2.bias, FDD.out2.sd)*100,2)[index.dir.z0,]
round(cbind(prob.out2.bias, prob.out2.sd)*100,2)[dir.out.z0,]
round(cbind(prob.out2.thre.bias, prob.out2.thre.sd)*100,2)[dir.out.z0,]


 cbind(state.out, round(prob.out2.tr*100,2),round(prob.out2.ave*100,2),round(prob.out2.thre.ave*100,2))
 round(cbind(prob.out2.tr, prob.out2.thre.bias, prob.out2.thre.sd)*100,2)[dir.out.z0,]


####
 save.image(file="super2_simu_N41.Rdata")     
 
##################
##### plot the out state with prob .
state.out.z0=state.out[dir.out.z0,]
norm.temp=sqrt(apply(state.out.z0^2,1,sum))
state.out.u= state.out.z0[,1:2]/matrix(norm.temp,nrow(state.out.z0),2)
len.a=15  ## baseline arrow length
wid.a=45  ##baseline arraow width

par(mfrow=c(2,1))

##true out-prob

prob.out.c=prob.out2.tr[dir.out.z0,]
plot(0,0,xlab="", ylab="",xlim=c(-2,2), ylim=c(-2,2), type='n', axes=FALSE, asp=1)

 for (i in 1:nrow(state.out.u)){
   
  if(prob.out.c[i]==0){
  arrow.c=state.out.u[i,]*len.a*0.15
   lwd.c =wid.a*0.07
   col.c=gray(0.4)
    arr.len.c=0
    lty.c=3
  }else{
  arrow.c=state.out.u[i,]*prob.out.c[i]*len.a
  lwd.c=wid.a*prob.out.c[i]
  #col.c=gray(prob.out.c[i])
  col.c=1
  arr.len.c=0.1
  lty.c=1
  }
  
  
  arrows(0,0, arrow.c[1],arrow.c[2],length=arr.len.c, angle=20, lwd=lwd.c, col=col.c, lty=lty.c)

 }
 
 ##mean estimated with 2-sd bands
  
 #prob.out.c=prob.out2.ave[dir.out.z0]
prob.out.c=prob.out2.thre.ave[dir.out.z0]
sum(prob.out.c) ## on 37% prob. in x-y plane
plot(0,0,xlab="", ylab="",xlim=c(-2,2), ylim=c(-2,2), type='n', axes=FALSE, asp=1)

 for (i in 1:nrow(state.out.u)){
   
  if(prob.out.c[i]<0){
  arrow.c=state.out.u[i,]*len.a*0.15
   lwd.c =wid.a*0.07
   col.c=gray(0.4)
    arr.len.c=0
    lty.c=3
  }else{
  arrow.c=state.out.u[i,]*prob.out.c[i]*len.a
  lwd.c=wid.a*prob.out.c[i]
   #lwd.c=wid.a*0.07
  #col.c=gray(prob.out.c[i])
  col.c=1
  arr.len.c=0.1
  lty.c=1
  }
  
  
  arrows(0,0, arrow.c[1],arrow.c[2],length=arr.len.c, angle=20, lwd=lwd.c, col=col.c, lty=lty.c)

   #lwd.c=wid.a*(prob.out.c[i]+2*prob.out2.sd[i])
    #arrows(0,0, arrow.c[1],arrow.c[2],length=arr.len.c, angle=20, lwd=lwd.c, col=gray(0.6), lty=lty.c)

 }


 par(mfrow=c(1,1))
 
 #### results

> round(cbind(FDD.out2.bias, FDD.out2.sd)*100,2)[index.dir.z0,]
#        sup 1 FDD.out2.sd
# dir 2   3.27        4.06
# dir 5   4.29        2.83
# dir 8   1.24        4.00
# dir 11  2.10        3.79
# dir 16  2.00        4.26
# dir 19  0.40        4.49
# dir 22  4.43        3.06
# dir 25  3.39        4.21
# > round(cbind(prob.out2.bias, prob.out2.sd)*100,2)[dir.out.z0,]
#              sup 1 prob.out2.sd
# out-state 3   1.11         1.62
# out-state 6   5.79         4.80
# out-state 9   2.58         2.17
# out-state 12  0.00         0.00
# out-state 17  1.72         2.46
# out-state 36  0.00         0.00
# out-state 45  0.66         1.94
# out-state 48  3.51         5.21
# out-state 51  1.24         3.15
# out-state 60  2.66         6.25
# out-state 71  0.00         0.00
# out-state 76  0.47         1.74
# out-state 79  0.00         0.00
# out-state 82  2.56         2.72
# out-state 91  6.13         4.52
# out-state 96  1.05         1.83



temp=cbind(abs(prob.out2.bias), prob.out2.sd, prob.out2.tr)
temp1=temp[prob.out2.tr>0,]
temp2=temp[prob.out2.tr==0,]

round(apply(temp1,2,mean),3)    ##0.025, 0.032, 0.083          
round(apply(temp1/temp1[,3],2, mean),3)    ##0.305, 0.474
round(apply(temp2,2,mean),3)     ##0.003, 0.003

