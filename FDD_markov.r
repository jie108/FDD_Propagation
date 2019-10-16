### calculate 3D- FDD using a Markov transition process

################ Part I: states, direction vectors, neighboring states
### generate the 3D-coordinate of all states
### there are two types of states: 27 in-state corresponding to the 27 voxels of a supervoxel
### 98 out-state, corresponding to the outward neighboring voxels of the 26 in-state voxel

### 26 direction/displacement vectors
dir.v=matrix(0,26,3)
count=1
for (i in c(-1,0,1)){
      for (j in c(-1,0,1)){
       for (k in c(-1,0,1)){
         cur.v=c(i,j,k)
         if(all(cur.v==0)==FALSE){
          dir.v[count,]=cur.v
          count=count+1
         }
  }
 }
}

colnames(dir.v)=c("x","y","z")
rownames(dir.v)=paste("V",1:26, sep="")
##which two directions are paired (about the center)
dir.pair=matrix(0,13,2)
dir.pair[,1]=1:13
dir.pair[,2]=27-(1:13)


### find the coordinates of all  states
state.out1=NULL
 for (l in 1:26){ ## get the out-state (absorbing states): note that their distance to (0,0,0) must be greater than sqrt{3}
  cur.s=dir.v[l,]
    
    for (i in c(-1,0,1)){
      for (j in c(-1,0,1)){
       for (k in c(-1,0,1)){
         cur.v=c(i,j,k)
         if(all(cur.v==0)==FALSE){
          temp=cur.s+cur.v
           if(sum(temp^2)>3)
            {
             state.out1=rbind(state.out1, temp)
            }
          
         }
    }
  }
 }
}
##get unique out-states
state.out=unique(state.out1)        ##98 by 3
colnames(state.out)=c("x","y","z")
rownames(state.out)=paste("out",1:98,sep="")

## get in-state
 state.in=rbind(c(0,0,0),dir.v)                      ##27 in states 
 colnames(state.in)=c("x","y","z")
 rownames(state.in)=paste("in",1:27,sep="")

## combine all states together and give each an index from 1:125
state.all=rbind(state.in, state.out)
state.all.index=cbind(state.all, 1:(nrow(state.all)))
colnames(state.all.index)=  c("x","y","z","index")    ##27+98=125 states in total 
 
 
#### for each in-state, find the state (index) each of the direction vectors points to 
##### for each out-state, and each direction vector, whether it is pointed to from an in-state along this direction
state.dir=matrix(0,27,26)
out.dir=matrix(0,98,26)

 for (l in 1:27){
   cur.s=state.in[l,]
    
    count=1
     for (i in c(-1,0,1)){
      for (j in c(-1,0,1)){
       for (k in c(-1,0,1)){
         cur.v=c(i,j,k)
         if(all(cur.v==0)==FALSE){
    
           cur.out=cur.s+cur.v
           temp=apply(state.all, 1, identical, cur.out)
           index.c=which(temp)
           state.dir[l,count]=index.c  
            if(index.c>27){
             out.dir[index.c-27,count]=1
            }
           
          count=count+1
         }
       }
    }
  }
 
}

rownames(state.dir)=paste("in",1:27,sep="")
colnames(state.dir)=paste("V",1:26, sep="")
rownames(out.dir)=paste("out",1:98,sep="")
colnames(out.dir)=paste("V",1:26, sep="")
tt=apply(out.dir,1,sum)
table(tt)
#tt
# 1  2  3  4  6  9 
# 8 24 12 24 24  6 
cbind(state.out,tt)

#############################################################
#### grouping the 98 out-states into 26 directions 
#############################################################
#### for each out-state: calculate its directional vector from (0,0,0)
angle.out=matrix(0,98,26)  ##angle between the directional vector (from 0,0,0) of an out-state and each of the 26 directions 
temp=sqrt(apply(dir.v^2,1,sum))
dir.v.norm=dir.v/matrix(temp,26,3,byrow=FALSE)
 
 for (l in 1:98){
  cur.s=state.out[l,]
  cur.v=cur.s/sqrt(sum(cur.s^2))
  cur.inner=dir.v.norm%*%matrix(cur.v,3,1)
  cur.inner=sign(cur.inner)*(abs(cur.inner)-(1e-6))  ##away a little bit from the boundary 1 and -1 to avoid numerical errors 
  angle.out[l,]=acos(cur.inner)
 }
 
 rownames(angle.out)=paste("out",1:98,sep="")
 colnames(angle.out)=paste("V",1:26,sep="")

angle.out.sort=t(apply(angle.out,1,sort))     ## for each out-state, sort its angles-with-directions profile 
angle.profile=unique(angle.out.sort)      ##unique angles-with-directions profiles: 6 unique profiles  
round(angle.profile/pi*180)

#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25] [,26]
#out1    0   35   35   35   55   55   55   71   71    71    90    90    90    90    90    90   109   109   109   125   125   125   145   145   145   180
#out2   16   19   45   45   48   48   55   71   76    76    79    79    90    90   101   101   104   104   109   125   132   132   135   135   161   164
#out3    0   35   35   45   45   60   60   60   60    90    90    90    90    90    90    90    90   120   120   120   120   135   135   145   145   180
#out5   19   30   30   35   55   62   62   66   66    73    73    90    90    90    90   107   107   114   114   118   118   125   145   150   150   161
#out6   18   27   39   39   51   51   63   72   72    72    75    75    90    90   105   105   108   108   108   117   129   129   141   141   153   162
#out9    0   45   45   45   45   55   55   55   55    90    90    90    90    90    90    90    90   125   125   125   125   135   135   135   135   180

### classify the 98 out-states by its angles-with-directions profile : 6 classes in total 
angle.status=numeric(98)   ##record angle profile of each out-state

  for (l in 1:98){
     cur.ap=angle.out.sort[l,]
     temp=NULL
   for (j in 1:6){
     temp[j]=sum(abs(cur.ap-angle.profile[j,]))
    }
   angle.status[l]=which.min(temp)
  
  }

  state.out.status=cbind(state.out, angle.status)
  temp.order=order(angle.status)
  state.out.status[temp.order,]

table(angle.status)
 #angle.status
 #1  2  3  4  5  6 
 #8 24 12 24 24  6 

  ##profile 1:    8 out-states: three 2, e.g., (2,2,2): consistent with a  3d-corner direction vector; no ambiguous in classification
  ##profile 2: 24 out-states: two 2, one 1, e.g., (2,2,1): in between two direction vectors (a 3d-corner and a 2d-corner, e.g., (1,1,1) and (1,1,0)); 55%-45% share
  ##profile 3: 12 out-states: two 2, one 0, e.g., (2,2,0): consistent with a 2d-corner direction vector;  no ambiguous in classification
  ## profile 4: 24 out-states: one 2, two 1, e.g., (2,1,1): closest (19 degree) to a 3d-corner direction (e.g., (1,1,1)), also close (30 degree) to two 2-d corner direction (e.g, (1,1,0),(1,0,1)); 45%-27.5%-27.5% share?
  ## profile 5: 24 out-states: one 2, one 1, one 0, e.g., (2,1,0):  closest (18 degree) to a 2-d corner direction (e.g., (1,1,0)), then (27degree) to a vertical direction (e.g., (1,1,0)); 60%-40% share?
  ## profile 6: 6 out-state: one 2, two 0, e.g., (2,0,0): consistent with a vertical direciton (e.g., (1,0,0)); no ambiguous
   
##
################## grouping scheme of the 98 out-states to 26 directions according to their profiles 
 group.out=matrix(0,98,26) ## what percentage of the absorbing prob. of an out-state should go to each of the 26 direction
  
  for (l in 1:98){
    temp=order(angle.out[l,])
    if(angle.status[l]==1 || angle.status[l]==3 || angle.status[l]==6){ ##profiles 1,3,6, classify to the closest direction
     group.out[l,temp[1]]=1
    }
    
    if(angle.status[l]==2){ ##55%-45%% share between the two cloest directions 
         
         group.out[l,temp[1]]=0.55
         group.out[l,temp[2]]=0.45
    }
  
   if(angle.status[l]==4){ ##45%-27.5-27.5% share between the threecloest directions 
         
         group.out[l,temp[1]]=0.45
         group.out[l,temp[2]]=0.275
         group.out[l,temp[3]]=0.275
    }
  
       if(angle.status[l]==5){ ##60%-40% share between the two cloest directions 
         
         group.out[l,temp[1]]=0.60
         group.out[l,temp[2]]=0.40
    }
  
  }
  
 
 rownames(group.out)=paste("out",1:98,sep="")
 colnames(group.out)=paste("V",1:26,sep="")
 apply(group.out,1,sum)         ##row sum equals to 1
 
 
 
 

##########
save.image(file="state.Rdata")


#######################
###################### Part II: get transition probability,    and FDD examples 
### specify transition probabilities based on FDD at each of the 27 in-state voxels 
## FDD.in: a 27 by 26 matrix of FDD:  symmetric about (0,0,0)
## each row corresponds to an in-state, the columns are the FDD values for the corresponding direction, >=0 and sum =1


 ################################## Part I : level 0 voxels flow ==> level 1 voxel
######
### Example 1:  FDD.in: prob. 1/6 at V14 V13, +/- (0,0,1), V11, V16 +/-(0,1,0), V5, V22 +/-(1,0,0), for each voxel 
### in terms of the supervoxel, more flow along the vertical directions
load("state.Rdata")
source("Tran_functions.r")

FDD.in=matrix(0,27,26)   
FDD.in[,c(5, 11, 13, 14,16,22)]=1/6    
                                      
apply(FDD.in, 1,sum) ##check: row sum=1
all(FDD.in>=0)       ## check: all nonnegative

## specify transition matrix: use state.dir and FDD.in
tran=Tran.prob(FDD.in, state.dir)
dim(tran)
apply(tran, 1,sum)       ## row sum =1: tran[i,j] is the transition prob. from state i to state j
 
 ############ 
 #tran2=tran%*%tran   ### two-step transition probabilities 
 
 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 temp[[2]]
 tran.abs=temp[[1]]
 
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]

 ## indices of out-states having nonzero absorbing prob. (starting from the center)
index.out=which(prob.abs.out!=0)
length(index.out) 
cbind(state.out[index.out,], prob.abs.out[index.out])


### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton

FDD.out=matrix(prob.abs.out,1,98)%*%group.out
FDD.out=matrix(FDD.out,26,1)
FDD.out=FDD.symm(FDD.out, dir.pair)

sum(FDD.out)       ## total prob. 1
cbind(dir.v, round(FDD.out*100,2))   ## more flow along the vertical directions (8%), and uniform along other directions (2~3%): 
 
save.image(file="example1.Rdata")
 
 
 
######
### Example 2:  FDD.in: prob. 1/4 at V11, V16 +/-(0,1,0), V5, V22 +/-(1,0,0), for each voxel; no z-direction flow; i.e., a x-y plane cross at each voxel 
### in terms of the supervoxel, the flow is rather uniform  within x-y plane (more along vertical directions), but no z-direction 
load("state.Rdata")
source("Tran_functions.r")

FDD.in=matrix(0,27,26)   
FDD.in[,c(5, 11, 16,22)]=1/4  
apply(FDD.in, 1,sum) ##check: row sum=1
all(FDD.in>=0)       ## check: all nonnegative

## specify transition matrix: use state.dir and FDD.in
tran=Tran.prob(FDD.in, state.dir)
dim(tran)
apply(tran, 1,sum)       ## row sum =1: tran[i,j] is the transition prob. from state i to state j
 

 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 temp[[2]]
 tran.abs=temp[[1]]
 
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]

 ## indices of out-states having nonzero absorbing prob. (starting from the center)
index.out=which(prob.abs.out!=0)
length(index.out) 
cbind(state.out[index.out,], prob.abs.out[index.out])


### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton

FDD.out=matrix(prob.abs.out,1,98)%*%group.out
FDD.out=matrix(FDD.out,26,1)
FDD.out=FDD.symm(FDD.out, dir.pair)
sum(FDD.out)       ## total prob. 1
cbind(dir.v, round(FDD.out*100,2))   ## no z-direction flow, 4 vertical dirction: 17.5%; 4 corner directions: 7.5%
 
save.image(file="example2.Rdata")
 
 
######
### Example 3:  FDD.in: prob. 1/6 on each vertical directions in the center voxel (voxel 1); flow continued to the neighboring voxels with same direciton 

load("state.Rdata")
source("Tran_functions.r")

FDD.in=matrix(0,27,26)   
FDD.in[1,c(5, 22, 11,16, 13,14)]=1/6

FDD.in[6,c(5,22)]=1/2
FDD.in[23,c(5,22)]=1/2

FDD.in[12,c(11,16)]=1/2
FDD.in[17,c(11,16)]=1/2

FDD.in[14,c(13,14)]=1/2
FDD.in[15,c(13,14)]=1/2

  
apply(FDD.in, 1,sum) ##check: row sum=1
all(FDD.in>=0)       ## check: all nonnegative

## specify transition matrix: use state.dir and FDD.in
tran=Tran.prob(FDD.in, state.dir)
dim(tran)
apply(tran, 1,sum)       ## row sum =1: tran[i,j] is the transition prob. from state i to state j
 

 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 temp[[2]]
 tran.abs=temp[[1]]
 
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]

 ## indices of out-states having nonzero absorbing prob. (starting from the center)
index.out=which(prob.abs.out!=0)
length(index.out) 
cbind(state.out[index.out,], prob.abs.out[index.out])


### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton

FDD.out=matrix(prob.abs.out,1,98)%*%group.out
FDD.out=matrix(FDD.out,26,1)
FDD.out=FDD.symm(FDD.out, dir.pair)
sum(FDD.out)       ## total prob. 1
cbind(dir.v, round(FDD.out*100,2))   ## 6 vertical dirctions: 16.7% each, consistent with expectation
 
save.image(file="example3.Rdata")
 
 ################
 ######
### Example 4a:  FDD.in: prob. 1/2 at V11, V16 +/-(0,1,0), for each voxel; no x- z-direction flow, only y-direction flow at each voxel.  
### in terms of the supervoxel, the flow should also be only along y-direction
load("state.Rdata")
source("Tran_functions.r")

FDD.in=matrix(0,27,26)   
FDD.in[,c(11, 16)]=1/2  
apply(FDD.in, 1,sum) ##check: row sum=1
all(FDD.in>=0)       ## check: all nonnegative

## specify transition matrix: use state.dir and FDD.in
tran=Tran.prob(FDD.in, state.dir)
dim(tran)
apply(tran, 1,sum)       ## row sum =1: tran[i,j] is the transition prob. from state i to state j
 

 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 temp[[2]]
 tran.abs=temp[[1]]
 
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]

 ## indices of out-states having nonzero absorbing prob. (starting from the center)
index.out=which(prob.abs.out!=0)
length(index.out) 
cbind(state.out[index.out,], prob.abs.out[index.out])


### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton

FDD.out=matrix(prob.abs.out,1,98)%*%group.out
FDD.out=matrix(FDD.out,26,1)
FDD.out=FDD.symm(FDD.out, dir.pair)
sum(FDD.out)       ## total prob. 1
cbind(dir.v, round(FDD.out*100,2))   ## starting from center,, no x- z-direction flow, only y-direction flow; consistent with expectation 
 
save.image(file="example4a.Rdata")
 
######################
### Example 4b:  FDD.in: prob. 1/2 at V5, V22 +/-(1,0,0), for each voxel; no y- z-direction flow, only x-direction flow at each voxel.  
### in terms of the supervoxel, the flow should also be only along y-direction
load("state.Rdata")
source("Tran_functions.r")

FDD.in=matrix(0,27,26)   
FDD.in[,c(5, 22)]=1/2  
apply(FDD.in, 1,sum) ##check: row sum=1
all(FDD.in>=0)       ## check: all nonnegative

## specify transition matrix: use state.dir and FDD.in
tran=Tran.prob(FDD.in, state.dir)
dim(tran)
apply(tran, 1,sum)       ## row sum =1: tran[i,j] is the transition prob. from state i to state j
 

 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 temp[[2]]
 tran.abs=temp[[1]]
 
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]

 ## indices of out-states having nonzero absorbing prob. (starting from the center)
index.out=which(prob.abs.out!=0)
length(index.out) 
cbind(state.out[index.out,], prob.abs.out[index.out])


### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton

FDD.out=matrix(prob.abs.out,1,98)%*%group.out
FDD.out=matrix(FDD.out,26,1)
FDD.out=FDD.symm(FDD.out, dir.pair)
sum(FDD.out)       ## total prob. 1
cbind(dir.v, round(FDD.out*100,2))   ## starting from center, no y- z-direction flow, only x-direction flow; consistent with expectation 
 
save.image(file="example4b.Rdata")
  
 
 
 ####################  Part II: level 1 supervoxel flows ==> level 2 supervoxel
 ### example I: 9 level one supervoxels; at each z-slides: the middle one is built from example 2; the four side ones around the middle one are from example4a, example 4b;
 ### the four corner ones have uniform (isotropic) flow; no z-direction flow
 
load("state.Rdata")
source("Tran_functions.r")

FDD.in.l1=matrix(1/26,27,26)       ##level 1 voxels FDD 
load("example2.Rdata")
index=c(1,14,15)
for (i in index){
FDD.in.l1[i,] =FDD.out         ##center voxels at each of the three z-slides
}

load("example4a.Rdata")
index=c(12,17,13,18,11,16)
for (i in index){
FDD.in.l1[i,] =FDD.out         ##two y-direction voxels  at each of the three z-slides
}  


load("example4b.Rdata")
index=c(6,23,7,24,5,22)
for (i in index){
FDD.in.l1[i,] =FDD.out         ##two x-direction voxels  at each of the three z-slides
}  
##

all(apply(FDD.in.l1,1, FDD.symm.check, dir.pair=dir.pair))  ##check symmetry
apply(FDD.in.l1, 1,sum) ##check: row sum=1
all(FDD.in.l1>=0)       ## check: all nonnegative

## specify transition matrix: use state.dir and FDD.in
tran=Tran.prob(FDD.in.l1, state.dir)
dim(tran)
apply(tran, 1,sum)       ## row sum =1: tran[i,j] is the transition prob. from state i to state j
 

 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 temp[[2]]
 tran.abs=temp[[1]]
 
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]

 ## indices of out-states having nonzero absorbing prob. (starting from the center)
index.out=which(prob.abs.out!=0)
length(index.out) 
cbind(state.out[index.out,], prob.abs.out[index.out])


### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton

FDD.out.l1=matrix(prob.abs.out,1,98)%*%group.out
FDD.out.l1=matrix(FDD.out.l1,26,1)
FDD.out.l1=FDD.symm(FDD.out.l1,dir.pair)
sum(FDD.out.l1)       ## total prob. 1
cbind(dir.v, round(FDD.out.l1*100,2))   ## starting from center,, no z-direction flow, most prob. (17.72*4) on the four x-y plane vertical directions; consistent with expectation 
 
FDD.out.thre=FDD.thre(FDD.out.l1, thre=1/26)    ##thresholding and renomarlize to get rid of random flow of water / incoherent minuate fibers 
#FDD.out.thre=FDD.symm(FDD.out.thre,dir.pair)
cbind(dir.v, round(FDD.out.thre*100,2)) 

save.image(file="example_super1.Rdata")                                      
                                       
 
####################
 ### example II: 9 level one supervoxels; at the middle z-slide (z=0):  x-y vertical cross flow passing the center level 1 super voxel 
 ###  also z-direction vertical flow passing the center level 1 super voxel
 ### So FDD of the center level 1 supervoxel is from example 1. In the middle z-slide, the four neighboring level 1 supervoxels have FDD from example4ab.   
 
 ### the four corner ones have uniform (isotropic) flow; no z-direction flow
 
load("state.Rdata")
source("Tran_functions.r")
FDD.in.l1=matrix(1/26,27,26)       ##level 1 voxels FDD 

###
load("example1.Rdata")
index=c(1)
for (i in index){
FDD.in.l1[i,] =FDD.out         ##center voxels at each of the three z-slides
}

load("example4a.Rdata")
index=c(12,17)
for (i in index){
FDD.in.l1[i,] =FDD.out         ##two y-direction voxels  at each of the three z-slides
}  


load("example4b.Rdata")
index=c(6,23)
for (i in index){
FDD.in.l1[i,] =FDD.out         ##two x-direction voxels  at each of the three z-slides
}  
##

index=c(14,15)
for (i in index){
FDD.in.l1[i,c(13,14)] =0.5         ##center voxels at z= +/- 1, vertical z-direction flow 
}  

all(apply(FDD.in.l1,1, FDD.symm.check, dir.pair=dir.pair))  ##check symmetry
apply(FDD.in.l1, 1,sum) ##check: row sum=1
all(FDD.in.l1>=0)       ## check: all nonnegative

## specify transition matrix: use state.dir and FDD.in
tran=Tran.prob(FDD.in.l1, state.dir)
dim(tran)
apply(tran, 1,sum)       ## row sum =1: tran[i,j] is the transition prob. from state i to state j
 

 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 temp[[2]]
 tran.abs=temp[[1]]
 
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]

 ## indices of out-states having nonzero absorbing prob. (starting from the center)
index.out=which(prob.abs.out!=0)
length(index.out) 
cbind(state.out[index.out,], prob.abs.out[index.out])


### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton

FDD.out.l1=matrix(prob.abs.out,1,98)%*%group.out
FDD.out.l1=matrix(FDD.out.l1,26,1)
FDD.out.l1=FDD.symm(FDD.out.l1,dir.pair)
sum(FDD.out.l1)       ## total prob. 1
cbind(dir.v, round(FDD.out.l1*100,2))   ## starting from center,  9% for each of the six vertical direction, ~2% prob. for each of the other 20 directions

FDD.out.thre=FDD.thre(FDD.out.l1, thre=1/26)
#FDD.out.thre=FDD.symm(FDD.out.thre,dir.pair)
cbind(dir.v, round(FDD.out.thre*100,2)) 

 
save.image(file="example_super2.Rdata")                                      
 
 ############################ 
 ############################Part III: level 0 voxles have FDD from a smooth FOD (having prob. on all directions), rather than Dirac functions.
 ### Example 1:  FDD.xyz for each voxel 
load("state.Rdata")
load("FDD.Rdata")
source("Tran_functions.r")

cbind(dir.v, round(FDD.xyz*100,2))

FDD.in=matrix(0,27,26)   
FDD.in= matrix(FDD.xyz,27,26, byrow=TRUE)   
                                      
apply(FDD.in, 1,sum) ##check: row sum=1
all(FDD.in>=0)       ## check: all nonnegative
FDD.symm.check(FDD.in,dir.pair)

## specify transition matrix: use state.dir and FDD.in
tran=Tran.prob(FDD.in, state.dir)
dim(tran)
apply(tran, 1,sum)       ## row sum =1: tran[i,j] is the transition prob. from state i to state j
 
 ############ 
 #tran2=tran%*%tran   ### two-step transition probabilities 
 
 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 temp[[2]]
 tran.abs=temp[[1]]
 
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]

 ## indices of out-states having nonzero absorbing prob. (starting from the center)
index.out=which(prob.abs.out!=0)
length(index.out) 
cbind(state.out[index.out,], prob.abs.out[index.out])


### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton

FDD.out=matrix(prob.abs.out,1,98)%*%group.out
FDD.out=matrix(FDD.out,26,1)
FDD.out=FDD.symm(FDD.out, dir.pair)

sum(FDD.out)       ## total prob. 1
cbind(dir.v, round(FDD.out*100,2))   ## the flow is rather uniform, a bit more flow along the six vertical directions; so pre-thresholding of FFD.in may be needed
 
save.image(file="example.fod1.Rdata")
 
 
 
 ######
### Example 2:  FDD.xyfor each voxel; no z-direction flow; i.e., a x-y plane cross at each voxel 
### in terms of the supervoxel, the flow is rather uniform  within x-y plane (more along vertical directions), but no z-direction 
load("state.Rdata")
load("FDD.Rdata")
source("Tran_functions.r")

FDD.in=matrix(0,27,26)   
FDD.in=matrix(FDD.xy,27,26,byrow=TRUE)  
apply(FDD.in, 1,sum) ##check: row sum=1
all(FDD.in>=0)       ## check: all nonnegative

## specify transition matrix: use state.dir and FDD.in
tran=Tran.prob(FDD.in, state.dir)
dim(tran)
apply(tran, 1,sum)       ## row sum =1: tran[i,j] is the transition prob. from state i to state j
 

 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 temp[[2]]
 tran.abs=temp[[1]]
 
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]

 ## indices of out-states having nonzero absorbing prob. (starting from the center)
index.out=which(prob.abs.out!=0)
length(index.out) 
cbind(state.out[index.out,], prob.abs.out[index.out])


### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton

FDD.out=matrix(prob.abs.out,1,98)%*%group.out
FDD.out=matrix(FDD.out,26,1)
FDD.out=FDD.symm(FDD.out, dir.pair)
sum(FDD.out)       ## total prob. 1
cbind(dir.v, round(FDD.out*100,2))   ## 4 x-y dirctions: 7.5%; small prob. on other directions.
cbind(dir.v, FDD.thre(FDD.out, thre=1/26))
 
save.image(file="example.fod2.Rdata")
 
 

 #####
### Example 3:  FDD.x22 and FDD.x122 (90 degree cross, 45 degree with x-)  for each voxel; no z-direction flow; i.e., a x-y plane cross at each voxel 
### in terms of the supervoxel, the flow is rather uniform  within x-y plane , but no z-direction 
load("state.Rdata")
load("FDD_new.Rdata")
source("Tran_functions.r")

cbind(dir.v, round(FDD.xy*100,2))
FDD.in=matrix(0,27,26)   
FDD.in=matrix(FDD.xy,27,26,byrow=TRUE)  
apply(FDD.in, 1,sum) ##check: row sum=1
all(FDD.in>=0)       ## check: all nonnegative

## specify transition matrix: use state.dir and FDD.in
tran=Tran.prob(FDD.in, state.dir)
dim(tran)
apply(tran, 1,sum)       ## row sum =1: tran[i,j] is the transition prob. from state i to state j
 

 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 temp[[2]]
 tran.abs=temp[[1]]
 
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]

 ## indices of out-states having nonzero absorbing prob. (starting from the center)
index.out=which(prob.abs.out!=0)
length(index.out) 
cbind(state.out[index.out,], prob.abs.out[index.out])


### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton
cbind(state.out, round(prob.abs.out*100,2))

FDD.out=matrix(prob.abs.out,1,98)%*%group.out
FDD.out=matrix(FDD.out,26,1)
FDD.out=FDD.symm(FDD.out, dir.pair)
sum(FDD.out)       ## total prob. 1
cbind(dir.v, round(FDD.out*100,2))   ## (+/- 1,=/- 1,0): 10.5%; (+/-1 1,0,0) and (0,+/- 1,0) 14%
cbind(dir.v, round(FDD.thre(FDD.out, thre=1/26)*100,2))
 
save.image(file="example.fod3.Rdata")


#####
### Example 4:  FDD.x22 and FDD.x122 (90 degree cross, 45 degree with x-)  at centering 5 voxels; rthe four corner one only has one direction ; no z-direction flow; i.e., a x-y plane cross at each voxel 
### in terms of the supervoxel, the flow is rather uniform  within x-y plane , but no z-direction 
load("state.Rdata")
load("FDD_new.Rdata")
source("Tran_functions.r")

cbind(dir.v, round(FDD.xy*100,2))
FDD.in=matrix(0,27,26)   
FDD.in=matrix(FDD.xy,27,26,byrow=TRUE) 
index.x=c(1,2,3,24,25,26)+1
for (i in index.x){
FDD.in[i,] =FDD.x22
}

index.y=c(7,8,9,18,19,20)+1
for (i in index.y){
FDD.in[i,] =FDD.x122
}




apply(FDD.in, 1,sum) ##check: row sum=1
all(FDD.in>=0)       ## check: all nonnegative

## specify transition matrix: use state.dir and FDD.in
tran=Tran.prob(FDD.in, state.dir)
dim(tran)
apply(tran, 1,sum)       ## row sum =1: tran[i,j] is the transition prob. from state i to state j
 

 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 temp[[2]]
 tran.abs=temp[[1]]
 
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]

 ## indices of out-states having nonzero absorbing prob. (starting from the center)
index.out=which(prob.abs.out!=0)
length(index.out) 
cbind(state.out[index.out,], prob.abs.out[index.out])


### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton
cbind(state.out, round(prob.abs.out*100,2))

FDD.out=matrix(prob.abs.out,1,98)%*%group.out
FDD.out=matrix(FDD.out,26,1)
FDD.out=FDD.symm(FDD.out, dir.pair)
sum(FDD.out)       ## total prob. 1
cbind(dir.v, round(FDD.out*100,2))   ## (+/- 1,=/- 1,0): 10.5%; (+/-1 1,0,0) and (0,+/- 1,0) 14%
cbind(dir.v, round(FDD.thre(FDD.out, thre=1/26)*100,2))
 
save.image(file="example.fod4.Rdata")
   

###########
#####
### Example 5:  FDD.x22 and FDD.x122 (90 degree cross, 45 degree with x-)  at the centering voxel; the four corner one only has no fiber, so uniform flow; 
### other voxels have either FDD.x22 or FDD.x122
##no z-direction flow; i.e., a x-y plane cross at each voxel 
### in terms of the supervoxel, the flow is rather uniform  within x-y plane , but no z-direction 
load("state.Rdata")
load("FDD_new.Rdata")
source("Tran_functions.r")

cbind(dir.v, round(FDD.xy*100,2))
FDD.in=matrix(1/26,27,26)   
FDD.in[1,]=FDD.xy  ##centering voxel
index.x=c(4,5,6,21,22,23)+1
for (i in index.x){
FDD.in[i,] =FDD.x22
}

index.y=c(10,11,12,15,16,17)+1
for (i in index.y){
FDD.in[i,] =FDD.x122
}




apply(FDD.in, 1,sum) ##check: row sum=1
all(FDD.in>=0)       ## check: all nonnegative

## specify transition matrix: use state.dir and FDD.in
tran=Tran.prob(FDD.in, state.dir)
dim(tran)
apply(tran, 1,sum)       ## row sum =1: tran[i,j] is the transition prob. from state i to state j
 

 ## get absorbing prob.
 temp=Tran.abs(tran, max.step=100, thre=1e-6)
 temp[[2]]
 tran.abs=temp[[1]]
 
 
 ## absorbing prob. for each state starting from the center voxel(state1)
 prob.abs=tran.abs[1,]
 all(abs(prob.abs[1:27])<1e-6)    ##for in-state, should be zero
 prob.abs.out=prob.abs[28:125]

 ## indices of out-states having nonzero absorbing prob. (starting from the center)
index.out=which(prob.abs.out>1/26)
length(index.out) 
cbind(state.out[index.out,], prob.abs.out[index.out])


### grouping the 98 out-state absorbing probabilities " prob.abs.out" into 26 directions "FDD" of the higher level supervoxel,  accroding to group.out classificaiton
cbind(state.out, round(prob.abs.out*100,2))

FDD.out=matrix(prob.abs.out,1,98)%*%group.out
FDD.out=matrix(FDD.out,26,1)
FDD.out=FDD.symm(FDD.out, dir.pair)
sum(FDD.out)       ## total prob. 1
cbind(dir.v, round(FDD.out*100,2))   ## (+/- 1,=/- 1,0): 7.5%; (+/-1 1,0,0) and (0,+/- 1,0) 17.5%
cbind(dir.v, round(FDD.thre(FDD.out, thre=1/26)*100,2))
 
save.image(file="example.fod5.Rdata")
   

                
                

                


 
 
 

