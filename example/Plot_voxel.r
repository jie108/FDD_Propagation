### plot a level2 super voxels with fiber bundles

load("state.Rdata")
############## Part I: plot in x-y plane, one level 2 supervoxel
## x-y coordinates of the center of each voxel

load("super2.Rdata")
coor.xy<-matrix(0,81,2)
count=1
 for (i in -4:4){
  for (j in -4:4){
  coor.xy[count,]=c(i,j)
  count=count+1
  }
 }


###############
##plot center of the voxels
 #dev.new(width=5, height=5)
 ##par(mfrow=c(2,1))
plot(x=0,y=0, xlab="", ylab="",xlim=c(-5,5), ylim=c(-5,5), type='n', axes=FALSE, asp=1)         ##asp=1, makes x-, y- same scale
 ## add voxel center
count=1
for (i in -4:4){
  for (j in -4:4){
    col.c=3
     if (is.element(i,c(-3,0,3)) && is.element(j,c(-3,0,3))){
    col.c=2 
     } 
    points(coor.xy[count,1],coor.xy[count,2], pch=16, cex=2, col=col.c)
    count=count+1
   }
 }  
 
 
## add voxel border lines
xlm=seq(-4.5,4.5,0.1)
ylm=seq(-4.5,4.5,0.1)
for (i in -5:4){
x.c=rep(i+0.5, length(ylm))
y.c=rep(i+0.5, length(xlm))

col.c=3
lty.c=2
lwd.c=4
if(is.element(i,c(-5,-2,1,4)))
{
col.c=2
lty.c=1
lwd.c=5
}

points(x.c, ylm, type='l', lty=lty.c, lwd=lwd.c, col=col.c)
points(xlm, y.c, type='l', lty=lty.c, lwd=lwd.c, col=col.c)

}

## add fibers:  
 fib1=pi/8       ## fiber 1 directions     , relative to x-axis
 fib2=pi/8+pi/2  # fiber 2 directions      , relative to x-axis
 wd1=3           ##fiber 1 width : 
 wd2=3           ## fiber 2 width: parallel moving the fiber along the orthgonal direction by wd, both up and down
 vox1=c(0,0)      ## starting voxel for fib1
 vox2=c(0,0)      ##starting voxel for fib2  
  l.t=11      ## length of fiber lines
 mag=1    ## magnify the spacing between fibers 
 
 ##fiber bundle 1
 b1=tan(fib1)  ##slope of fiber 1
  for (i in -wd1:wd1){
    for (j in 0:0){
     vox.c=vox1+c(i,j)*mag
     a1.c=vox.c[2]-b1*vox.c[1] ##intercept of fiber 1
     line.c=a1.c+b1*xlm
     len.c=sqrt((xlm[1]-xlm[length(xlm)])^2+(line.c[1]-line.c[length(xlm)])^2)
     fac.c=l.t/len.c  ## factor 
     
     points(fac.c*xlm,fac.c*line.c, type='l', lwd=4, col=gray(0.4))
    } 
  }
 
 

 ##fiber bundle 2 
 b2=tan(fib2)  ##slope of fiber 1
  for (i in 0:0){
    for (j in -wd2:wd2){
     vox.c=vox2+c(i,j)*mag
     a2.c=vox.c[2]-b2*vox.c[1] ##intercept of fiber 1
     line.c=(ylm-a2.c)/b2
     len.c=sqrt((ylm[1]-ylm[length(ylm)])^2+(line.c[1]-line.c[length(ylm)])^2)
     fac.c=l.t/len.c  ## factor 
     
     points(fac.c*line.c, fac.c*ylm, type='l', lwd=4, col=gray(0.4))
    } 
  }

##################
##### plot the out state with prob .
state.out.z0=state.out[dir.out.z0,]
norm.temp=sqrt(apply(state.out.z0^2,1,sum))
state.out.u= state.out.z0[,1:2]/matrix(norm.temp,nrow(state.out.z0),2)
len.a=15  ## baseline arrow length
wid.a=55  ##baseline arraow width

#prob.out.c=prob.out2.temp
prop.out.c=prob.out2.prop[dir.out.z0]

plot(0,0,xlab="", ylab="",xlim=c(-2,2), ylim=c(-2,2), type='n', axes=FALSE, asp=1)

 for (i in 1:nrow(state.out.u)){
   
  if(prob.out.c[i]==0){
  arrow.c=state.out.u[i,]*len.a*0.15
   lwd.c =wid.a*0.1
   col.c=gray(0.1)
    arr.len.c=0
    lty.c=3
  }else{
  arrow.c=state.out.u[i,]*prob.out.c[i]*len.a
  lwd.c=wid.a*prob.out.c[i]
  #col.c=gray(prob.out.c[i])
  col.c=1
  arr.len.c=0.15
  lty.c=1
  }
  
  
  arrows(0,0, arrow.c[1],arrow.c[2],length=arr.len.c, angle=20, lwd=lwd.c, col=col.c, lty=lty.c)

 }

 par(mfrow=c(1,1))


###################split screen
#split.screen( figs = c( 1, 2 ) ) ##screen 1, 2
#split.screen( figs = c( 3, 1 ), screen = 1 )   ##screen 3, 4,5
#split.screen( figs = c( 2, 3 ), screen = 2 )   ## screen 6,7,8,9,10,11
  
  
############################
###################### part II:        three level 1 supervoxels, with different fiber width passing through it. 
## x-y coordinates of the center of each voxel
coor1.xy<-matrix(0,9,2)
count=1
 for (i in -1:1){
  for (j in -1:1){
  coor1.xy[count,]=c(i,j)
  count=count+1
  }
 }

##################################
###############
xlm=seq(-1.5,1.5,0.1)
ylm=seq(-1.5,1.5,0.1)
 
 fib1=pi/8       ## fiber 1 directions     , relative to x-axis
 fib2=pi/8+pi/2  # fiber 2 directions      , relative to x-axis
 vox1=c(0,0)      ## starting voxel for fib1
 vox2=c(0,0)      ##starting voxel for fib2  
 l.t=c(4,4,3.5)      ## length of fiber lines
 l.num=c(1,2,3)  ## number of lines
 mag=c(0.5,0.7,0.9)    ## magnify the spacing between fibers 

################################    
par(mfrow=c(2,3))
for (k in 1:3){           ##screen 6,8,10
#screen(2*k+4)
##plot center of the voxels
 #dev.new(width=5, height=5)
plot(x=0,y=0, xlab="", ylab="",xlim=c(-2,2), ylim=c(-2,2), type='n', axes=FALSE, asp=1)         ##asp=1, makes x-, y- same scale
points(coor1.xy[,1],coor1.xy[,2], pch=16, cex=2, col=3)
points(0,0,pch=16, cex=2, col=2)

 
 
## add voxel border lines

for (i in -2:1){
x.c=rep(i+0.5, length(ylm))
y.c=rep(i+0.5, length(xlm))

col.c=3
lty.c=2
lwd.c=4

if(is.element(i,c(-2,1)))
{
col.c=2
lty.c=1
lwd.c=5
}

points(x.c, ylm, type='l', lty=lty.c, lwd=lwd.c, col=col.c)
points(xlm, y.c, type='l', lty=lty.c, lwd=lwd.c, col=col.c)

}

 
## add fibers:  

 wd1=l.num[k]           ##fiber 1 width : 
 wd2=l.num[k]           ## fiber 2 width: parallel moving the fiber along the orthgonal direction by wd, both up and down

 ##fiber bundle 1
 b1=tan(fib1)  ##slope of fiber 1
  for (i in -wd1:wd1){
    for (j in 0:0){
     vox.c=vox1+c(i,j)*mag[k]
     a1.c=vox.c[2]-b1*vox.c[1] ##intercept of fiber 1
     line.c=a1.c+b1*xlm
     len.c=sqrt((xlm[1]-xlm[length(xlm)])^2+(line.c[1]-line.c[length(xlm)])^2)
     fac.c=l.t[k]/len.c  ## factor 
     
     points(fac.c*xlm,fac.c*line.c, type='l', lwd=4, col=gray(0.4))
  
    } 
  }
 
 

 ##fiber bundle 2 
 b2=tan(fib2)  ##slope of fiber 1
  for (i in 0:0){
    for (j in -wd2:wd2){
     vox.c=vox2+c(i,j)*mag[k]
     a2.c=vox.c[2]-b2*vox.c[1] ##intercept of fiber 1
     line.c=(ylm-a2.c)/b2
     len.c=sqrt((ylm[1]-ylm[length(ylm)])^2+(line.c[1]-line.c[length(ylm)])^2)
     fac.c=l.t[k]/len.c  ## factor 
     
     points(fac.c*line.c, fac.c*ylm, type='l', lwd=4, col=gray(0.4))
    } 
  }
 
 }  ##end of k loop
 
######### add level 1 out-ward prob.
 load("super1.Rdata")

##### plot the out state with prob .
state.out.z0=state.out[dir.out.z0,]
norm.temp=sqrt(apply(state.out.z0^2,1,sum))
state.out.u= state.out.z0[,1:2]/matrix(norm.temp,nrow(state.out.z0),2)
len.a=15  ## baseline arrow length
wid.a=55  ##baseline arraow width


##par(mfrow=c(1,3))
for (k in 1:3){
##screen(2*k+5)

#prob.out.c=prob.out.temp[,k]
prob.out.c=prob.out.prop[dir.out.z0,k]

plot(0,0,xlab="", ylab="",xlim=c(-2,2), ylim=c(-2,2), type='n', axes=FALSE, asp=1)

 for (i in 1:nrow(state.out.u)){
   if(prob.out.c[i]==0){
    arrow.c=state.out.u[i,]*len.a*0.1
    ##arrow.c=state.out.u[i,]*len.a*0.2
    lwd.c=wid.a*0.08
    lty.c=3
    col.c=gray(0.1)
    arr.len.c=0
   }else{
  arrow.c=state.out.u[i,]*prob.out.c[i]*len.a
  
  ##arrow.c=state.out.u[i,]*len.a*0.2
  lwd.c=wid.a*prob.out.c[i]
  #col.c=gray(prob.out.c[i])
  lty.c=1
   
 col.c=1
 arr.len.c=0.1
  }

  arrows(0,0, arrow.c[1],arrow.c[2],length=arr.len.c, angle=20, lwd=lwd.c, col=col.c, lty=lty.c)

 }
} ##  end of k

##par(mfrow=c(1,1))


###################
############ draw plots for three different FDD.in
 load("super1.Rdata")
FDD1=FDD.in[1,index.dir.z0,1]  ## with crossing
FDD2=FDD.in[2,index.dir.z0,1]  ## with crossing
FDD3=FDD.in[8,index.dir.z0,1]  ## with crossing
FDD.u=cbind(FDD1, FDD2, FDD3)

#######   fiber information
fib1=pi/8       ## fiber 1 directions     , relative to x-axis
 fib2=pi/8+pi/2  # fiber 2 directions      , relative to x-axis
 vox1=c(0,0)      ## starting voxel for fib1
 vox2=c(0,0)      ##starting voxel for fib2  
 l.t=3      ## length of fiber lines
 l.num=c(1,1,1)  ## number of lines
 mag=c(0.6,0.6,0.6)    ## magnify the spacing between fibers 
 
b1=tan(fib1)  ##slope of fiber 1
 b2=tan(fib2)  ##slope of fiber 2

##### 
dir.z0=dir.v.norm[index.dir.z0,1:2]
len.a=c(7,5,5)  ## baseline arrow length
wid.a=c(35,30,30)  ##baseline arraow width

lim=1.5
xlm.c=seq(-lim,lim,0.01)
ylm.c=seq(-lim,lim,0.01)

########
par(mfrow=c(3,1))
for (k in 1:3){
prob.out.c=FDD.u[,k]
plot(0,0,xlab="", ylab="",xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), type='n', axes=FALSE, asp=1)

 for (i in 1:nrow(dir.z0)){
 
  if(prob.out.c[i]==0){
  arrow.c=dir.z0[i,]*len.a[k]*0.25
  #arrow.c=dir.z0[i,]*len.a[k]*0.5
  lwd.c=wid.a*0.2
  #col.c=gray(prob.out.c[i])
  col.c=gray(0.1)
  arr.len.c=0
  lty.c=3 
 }else{
  arrow.c=dir.z0[i,]*prob.out.c[i]*len.a[k]
  #arrow.c=dir.z0[i,]*len.a[k]*0.5
  lwd.c=wid.a*prob.out.c[i]
  ##col.c=gray(prob.out.c[i])
  col.c=1
  arr.len.c=0.2
  lty.c=1
  }
  
  arrows(0,0, arrow.c[1],arrow.c[2],length=arr.len.c, angle=20, lwd=lwd.c, col=col.c, lty=lty.c)

 }
 
 ## add voxle border lines
 x.c=rep(-lim, length(ylm.c))
y.c=rep(-lim, length(xlm.c))
points(x.c, ylm.c, type='l', lty=2, lwd=4, col=3)
points(xlm.c, y.c, type='l', lty=2, lwd=4, col=3)
 
  x.c=rep(lim, length(ylm.c))
y.c=rep(lim, length(xlm.c))
points(x.c, ylm.c, type='l', lty=2, lwd=4, col=3)
points(xlm.c, y.c, type='l', lty=2, lwd=4, col=3)

### add fibers

## add fibers:  

 wd1=l.num[k]           ##fiber 1 width : 
 wd2=l.num[k]           ## fiber 2 width: parallel moving the fiber along the orthgonal direction by wd, both up and down

 if(k==1 || k==2){
 ##fiber bundle 1
 
  for (i in -wd1:wd1){
    for (j in 0:0){
     vox.c=vox1+c(i,j)*mag[k]
     a1.c=vox.c[2]-b1*vox.c[1] ##intercept of fiber 1
     line.c=a1.c+b1*xlm.c
     len.c=sqrt((xlm.c[1]-xlm.c[length(xlm.c)])^2+(line.c[1]-line.c[length(xlm.c)])^2)
     fac.c=l.t/len.c  ## factor 
     
     points(fac.c*xlm.c,fac.c*line.c, type='l', lwd=4, lty=1,col=gray(0.4))
  
    } 
  }
 
} 


if(k==1 || k==3){
 ##fiber bundle 2 

  for (i in 0:0){
    for (j in -wd2:wd2){
     vox.c=vox2+c(i,j)*mag[k]
     a2.c=vox.c[2]-b2*vox.c[1] ##intercept of fiber 1
     line.c=(ylm.c-a2.c)/b2
     len.c=sqrt((ylm.c[1]-ylm.c[length(ylm.c)])^2+(line.c[1]-line.c[length(ylm.c)])^2)
     fac.c=l.t/len.c  ## factor 
     
     points(fac.c*line.c, fac.c*ylm.c, type='l', lwd=4, lty=1, col=gray(0.4))
    } 
  }
  
}  
  
}## end k loop

par(mfrow=c(1,1))


#############################
################ draw FDD.x45
######## 
###   45 degree with x

theta=pi/2  ## polar angle, [0,pi)
 phi=pi/4  ## azimuthal angle, [0,2*pi)
 x0=sin(theta)*cos(phi)
 y0=sin(theta)*sin(phi)
 z0=cos(theta)
 dir.x45=c(x0,y0,z0)
 
 A=Orth.tran(theta,phi) 
 ### apply t(A) to the grid
 grid.0=t(A)%*%grid.ea
 
 ######### FDD
  FDD.x45=FOD.to.FDD_new(FOD.z, grid.0, dir.v.norm)      
  FDD.x45=FDD.thre(FDD.x45,thre=1/26) 
  FDD.x45=FDD.symm(FDD.x45, dir.pair)        ##check symmetry: true
   FDD.symm.check(FDD.x45,dir.pair)              ##check symmetry: true
  sum(FDD.x45)
  cbind(dir.v, round(FDD.x45*100,2))
  
  
  ### draw FDD
  vox1=c(0,0)      ## starting voxel for fib1
 vox2=c(0,0)      ##starting voxel for fib2  
 l.t=3      ## length of fiber lines
 l.num=c(1,1,1)  ## number of lines
 mag=c(0.6,0.6,0.6)    ## magnify the spacing between fibers 
 
 fib1=phi       ## fiber 1 directions     , relative to x-axis
 #fib2=pi/8+pi/2  # fiber 2 directions      , relative to x-axis
b1=tan(fib1)  ##slope of fiber 1
# b2=tan(fib2)  ##slope of fiber 2

##### 
dir.z0=dir.v.norm[index.dir.z0,1:2]
len.a=c(7,5,3)  ## baseline arrow length
wid.a=c(35,30,30)  ##baseline arraow width


lim=1.5
xlm.c=seq(-lim,lim,0.01)
ylm.c=seq(-lim,lim,0.01)

########
prob.out.c=FDD.x45[index.dir.z0]
plot(0,0,xlab="", ylab="",xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), type='n', axes=FALSE, asp=1)

k=3
 for (i in 1:nrow(dir.z0)){
 
  if(prob.out.c[i]==0){
  arrow.c=dir.z0[i,]*len.a[k]*0.4
  #arrow.c=dir.z0[i,]*len.a[k]*0.5
  lwd.c=wid.a[k]*0.25
  #col.c=gray(prob.out.c[i])
  col.c=gray(0.1)
  arr.len.c=0
  lty.c=3 
 }else{
  arrow.c=dir.z0[i,]*prob.out.c[i]*len.a[k]
  #arrow.c=dir.z0[i,]*len.a[k]*0.5
  lwd.c=wid.a[k]*prob.out.c[i]
  ##col.c=gray(prob.out.c[i])
  col.c=1
  arr.len.c=0.2
  lty.c=1
  }
  
  arrows(0,0, arrow.c[1],arrow.c[2],length=arr.len.c, angle=20, lwd=lwd.c, col=col.c, lty=lty.c)

 }
 
 ## add voxle border lines
 x.c=rep(-lim, length(ylm.c))
y.c=rep(-lim, length(xlm.c))
points(x.c, ylm.c, type='l', lty=2, lwd=4, col=3)
points(xlm.c, y.c, type='l', lty=2, lwd=4, col=3)
 
  x.c=rep(lim, length(ylm.c))
y.c=rep(lim, length(xlm.c))
points(x.c, ylm.c, type='l', lty=2, lwd=4, col=3)
points(xlm.c, y.c, type='l', lty=2, lwd=4, col=3)

### add fibers

## add fibers:  

 wd1=l.num[k]           ##fiber 1 width : 
 wd2=l.num[k]           ## fiber 2 width: parallel moving the fiber along the orthgonal direction by wd, both up and down

 
 ##fiber bundle 1
 
  for (i in -wd1:wd1){
    for (j in 0:0){
     vox.c=vox1+c(i,j)*mag[k]
     a1.c=vox.c[2]-b1*vox.c[1] ##intercept of fiber 1
     line.c=a1.c+b1*xlm.c
     len.c=sqrt((xlm.c[1]-xlm.c[length(xlm.c)])^2+(line.c[1]-line.c[length(xlm.c)])^2)
     fac.c=l.t/len.c  ## factor 
     
     points(fac.c*xlm.c,fac.c*line.c, type='l', lwd=4, lty=1,col=gray(0.4))
  
    } 
  }
 


