#Latest 
library(rPython)
library(softImpute)
library(NMF)
library(ggplot2)
library(RcppOctave)

#.CallOctave("solver_sNuclearBP.m",)
###########FUNCTIONS#############

#1. Generating and drawing random matrix

graphx <- function(ad, x2s,y2s,xbar1s, xbar2s,ybar1s, ybar1ps, MaxNumColorss,warmcolors, coldcolors, randchances) {

  x <- ifelse(ad==0, ad, NA)
  par(fig=c(0,x2s,0,y2s), new=FALSE)
  image(x=t(x[nrow(x):1,]), axes=FALSE, zlim=c(min(x,na.rm=T),max(x,na.rm=T)), col="#38424E", frame.plot=TRUE) #38424E C0D8D8 F6F6F6 A9B0B3
  
  par(fig=c(0,x2s,0,y2s), new=TRUE)
  grid(nx = ncol(ad), ny = nrow(ad), col ="black", lty = "solid", lwd = par("lwd"), equilogs = TRUE) #rgb(red=244/256, green=244/256, blue=244/256, alpha=0.2)
  
  par(fig=c(0,x2s,0,y2s), new=TRUE)
  x <- ifelse(ad>0, ad, NA)
  image(x=t(x[nrow(x):1,]), axes=FALSE, zlim=c(min(x,na.rm=T),max(x,na.rm=T)), col=warmcolors, font.main = 4)#, main="X = Axon Matrix . Dendrites Matrix = Connection Matrix"
  
  x <- ifelse(ad<0, ad, NA) 
  par(fig=c(0,x2s,0,y2s), new=TRUE)
  image(x=t(x[nrow(x):1,]), axes=FALSE, zlim=c(min(x,na.rm=T),max(x,na.rm=T)), col=coldcolors)

  rand <- matrix(runif(n=length(ad), min=0 , max=1),ncol=ncol(ad))
  x <-ifelse(ad!=0 & rand>randchances, ad, NA) #
  par(fig=c(0,x2s,0,y2s), new=TRUE)
  image(x=t(x[nrow(x):1,]), axes=FALSE, zlim=c(min(x,na.rm=T),max(x,na.rm=T)), col="#FAF0FE")#FF6666 FAF0FE BD99E1 F4F7F7 E3F6F3 C4EDE4
  
  x <-ifelse(ad!=0 & rand>randchances, NA, ad) #
  
  write.csv(x,"x.csv")

  
  par(fig=c(xbar1s,xbar2s,ybar1s,y2s),new=TRUE)
  y <- seq(min(ad[ad>0]),max(ad),length.out=MaxNumColorss)
  image(y=y,z=t(2:MaxNumColorss), axes=FALSE, frame.plot=TRUE, ann=FALSE, col=warmcolors[2:MaxNumColorss])
  
  par(fig=c(xbar1s,xbar2s,0,ybar1ps),new=TRUE)
  y <- seq(min(ad),max(ad[ad<0]),length.out=MaxNumColorss)
  image(y=y,z=t(2:MaxNumColorss), axes=FALSE, frame.plot=TRUE, ann=FALSE, col=coldcolors[2:MaxNumColorss])
  x <- x[-(c(1)),]
  return (x) 
}

#2. Drawing solved data

graph <- function(ad, file, x2s,y2s,xbar1s, xbar2s, ybar1s, ybar1ps, MaxNumColorss,warmcolors, coldcolors) {
  xs <- as.matrix(read.csv(file, sep=","))
  xs <- xs[,-c(1)]
  x <- xs
  
  x <- ifelse(ad==0, ad, NA)
  par(fig=c(0,x2s,0,y2s), new=FALSE)
  image(x=t(x[nrow(x):1,]), axes=FALSE, zlim=c(min(x,na.rm=T),max(x,na.rm=T)), col="#38424E", frame.plot=TRUE) #38424E C0D8D8 F6F6F6
  
  par(fig=c(0,x2s,0,y2s), new=TRUE)
  grid(nx = ncol(ad), ny = nrow(ad), col ="black", lty = "solid", lwd = par("lwd"), equilogs = TRUE) #rgb(red=244/256, green=244/256, blue=244/256, alpha=0.2)
  
  par(fig=c(0,x2s,0,y2s), new=TRUE)
  x <- ifelse(xs>0 & ad!=0, xs, NA)
  image(x=t(x[nrow(x):1,]), axes=FALSE, zlim=c(min(x,na.rm=T),max(x,na.rm=T)), col=warmcolors, font.main = 4)#, main="X = Axon Matrix . Dendrites Matrix = Connection Matrix"
  
  x <-ifelse(xs<0 & ad!=0, xs, NA)
  par(fig=c(0,x2s,0,y2s), new=TRUE)
  image(x=t(x[nrow(x):1,]), axes=FALSE, zlim=c(min(x,na.rm=T),max(x,na.rm=T)), col=coldcolors)
  
  
  par(fig=c(xbar1s,xbar2s,ybar1s,y2s),new=TRUE)
  y <- seq(min(ad[ad>0]),max(ad),length.out=MaxNumColorss)
  image(y=y,z=t(2:MaxNumColorss), axes=FALSE, frame.plot=TRUE, ann=FALSE, col=warmcolors[2:MaxNumColorss])
  
  par(fig=c(xbar1s,xbar2s,0,ybar1ps),new=TRUE)
  y <- seq(min(ad),max(ad[ad<0]),length.out=MaxNumColorss)
  image(y=y,z=t(2:MaxNumColorss), axes=FALSE, frame.plot=TRUE, ann=FALSE, col=coldcolors[2:MaxNumColorss])
  
  xs <-xs[-c(1),]

  return (xs)
}

#3. SVD of softimpute
softimpute <- function(number) {
  set.seed(number)
  bob <- as.matrix(read.csv("x.csv", sep=","))
  bob <- bob[,-c(1)]
  fits=softImpute(bob,trace=TRUE,type="svd")
  fita=softImpute(bob,trace=TRUE)
  fits2=softImpute(bob,rank.max=121,lambda=10,trace=TRUE,type="svd")
  fita2=softImpute(bob,rank.max=121,lambda=10,trace=TRUE)
  fits2$d
  write.csv(complete(bob,fits2),"x2.csv")
  
}

#4. ALS of softimpute
#softimpute2 <- function(number) {
  #set.seed(number)
  #bob <- as.matrix(read.csv("x.csv", sep=","))
  #bob <- bob[,-c(1)]
  #fits=softImpute(bob,trace=TRUE,type="als")
  #fita=softImpute(bob,trace=TRUE)
  #fits2=softImpute(bob,rank.max=3,lambda=1.9,trace=TRUE,type="als")
  #fita2=softImpute(bob,rank.max=3,lambda=1.9,trace=TRUE)
  #fits2$d
  #completion <- complete(bob,fits2)
  #write.csv(completion,"x3.csv")
#}

#5. Generating successrate/ Note: xna = surrogate; xa = original; x = testing
successrate <- function(xna, xa, x){
  sd_pos <- sd(xa[xa>0],na.rm=TRUE)
  sd_neg <- sd(xa[xa<0],na.rm=TRUE)
  poscount <- ifelse(((xa[is.na(xna) & xa>0] - sd_pos) < x[is.na(xna) & xa>0]) & (x[is.na(xna) & a>0] < (xa[is.na(xna) & xa>0] + sd_pos)), 1, 0)
  negcount <- ifelse( ((xa[is.na(xna) & xa<0] - sd_neg) < x[is.na(xna) & xa<0]) & (x[is.na(xna) & a<0] < (xa[is.na(xna) & xa<0] + sd_neg)), 1, 0)
  return (round (100*(sum(poscount) + sum(negcount))/sum(is.na(xna)),2))

}
######################################################################################################################
######################################################################################################################
######################################################################################################################
################       MAIN         #################

pdf(file = "myplot.pdf")

axons     <-   as.matrix(read.csv("Axons.csv", sep=","))
dendrites <- t(as.matrix(read.csv("Dendrites.csv", sep=","))) #tranpose to change row and col
x.conNum  <- axons %*% dendrites 
write.csv(x.conNum,"xconnum.csv")

axons[axons>0] <- runif(n=length(axons[axons>0]), min= 2.3, max= 5)
axons[axons<0] <- runif(n=length(axons[axons<0]), min=-2/3, max=-0.25)
dendrites[dendrites>0] <- runif(n=length(dendrites[dendrites>0]), min=2 , max=3)
x.ad <- axons %*% dendrites #dot product
write.csv(x.ad,"xad.csv")

#Initializing values
x2=0.85
y2=0.85
xbar1=0.68
xbar2=0.88
ybar1=0.294
ybar1p=0.556
MaxNumColors <- 2048
warmcolor=colorRampPalette(c("#D53D0C", "yellow"))(MaxNumColors)
coldcolor=colorRampPalette(c("#C6CAD3","#096F9F"))(MaxNumColors)

################Creating Original matrix################
a=graph(x.ad, "xad.csv", x2, y2,xbar1,xbar2,ybar1,ybar1p,MaxNumColors,warmcolor, coldcolor)

###############Creating SVD seednumber#################
seednumber <- Sys.time() 
seednumber <- gsub(':','',seednumber) 
seednumber <- gsub('-','',seednumber) 
seednumber <- gsub(' ','',seednumber) 
seednumber <- gsub('2016','',seednumber) 
seednumber <- gsub('0','',seednumber)

#x <- graphx(x.ad, x2,y2,xbar1, xbar2, ybar1, ybar1p, MaxNumColors,warmcolor, coldcolor, 0.5) #
#softimpute(seednumber) #
#b=graph(x.ad, "x2.csv", x2, y2,xbar1,xbar2,ybar1, ybar1p,MaxNumColors,warmcolor, coldcolor)
#successrate(x, a, b)

#########Creating "g", the list of lists/ gx = x values, gy = y values #########
gx <- c()
gy <- c()

g = list() 
for ( i in 1:10 ) {
  g[i] <- append(g[i],list())
} 

#Gaining successrates for different randoms and/or noise
count=10 #number of loops 
for (j in 1:count) {
  for (i in 1:10) {
    
    rand <- matrix(runif(n=length(x.ad), min=0 , max=1),ncol=ncol(x.ad))
    pmatrix <- matrix(runif(n=length(x.ad), min=min(x.ad[x.ad>0],na.rm=TRUE), max=max(x.ad[x.ad>0],na.rm=TRUE)))
    posmatrix <- ifelse(rand>i/10.0 & x.ad>0, pmatrix, 0)
    
    rand <- matrix(runif(n=length(x.ad), min=0 , max=1),ncol=ncol(x.ad))
    nmatrix <- matrix(runif(n=length(x.ad), min=min(x.ad[x.ad<0],na.rm=TRUE), max=max(x.ad[x.ad<0],na.rm=TRUE)))
    negmatrix <- ifelse(rand>i/10.0 & x.ad<0, nmatrix, 0)
    
    xrandom <- posmatrix + negmatrix 
    xrandom <- ifelse(xrandom==0 & x.ad!=0,x.ad,xrandom)
    
    x <- graphx(xrandom, x2,y2,xbar1, xbar2, ybar1,ybar1p, MaxNumColors,warmcolor, coldcolor, 0.5) #Note: 0.1 = 90% random
    seednumber <- Sys.time() 
    seednumber <- gsub(':','',seednumber) 
    seednumber <- gsub('-','',seednumber) 
    seednumber <- gsub(' ','',seednumber) 
    seednumber <- gsub('2016','',seednumber) 
    seednumber <- gsub('0','',seednumber)
    softimpute(seednumber)
    b=graph(x.ad, "x2.csv", x2, y2,xbar1,xbar2,ybar1, ybar1p,MaxNumColors,warmcolor, coldcolor)
    gx <- c(gx, 100-(i*10.0))
    
    s <-successrate(x, a, b)
    gy <- c(gy, s) 
    g[[i]] <- append(g[[i]],s) 
    
  }
}

#Compiling and saving datas
data.raw2 <- readRDS("compiled.Rda")
data.raw <- data.frame(group = rep(c('99','90','80','70','60','50','40','30','20','10'),each=count),data = c(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]],g[[6]],g[[7]],g[[8]],g[[9]],g[[10]]) )
data.raw <- rbind(data.raw, data.raw2)
#data.raw <- append(data.raw, data.raw2)
saveRDS(data.raw,file="compiled.Rda")

gx2 <- readRDS("compiledx.Rda")
gy2 <- readRDS("compiledy.Rda")


#gx <- rbind(gx2, gx)
#gy <- rbind(gy2, gy)
gx <- append(gx2, gx)
gy <- append(gy2, gy)

saveRDS(gx, file="compiledx.Rda")
saveRDS(gy, file="compiledy.Rda")

#Obtaining sd, mean, etc.
data.summary <- data.frame(
  Noise=levels(data.raw$group),
  Success_Rate=tapply(data.raw$data, data.raw$group, mean),
  n=tapply(data.raw$data, data.raw$group, length),
  sd=tapply(data.raw$data, data.raw$group, sd)
)

#Graphing scatterplot with linear regresson 
par(fig=c(0,x2,0,y2), new=FALSE)

plot(gx, gy, pch = 5, cex = 1, col = "blue", main = "SVD- Noise vs. Success Rate (50% random).", xlab = "Noise (%)", ylab = "Success Rate (%)")
abline(lm(gy ~ gx))

#Error bars 
par(fig=c(0,x2,0,y2), new=FALSE)
ggplot(data.summary, aes(x = Noise, y = Success_Rate)) +  
  geom_point(shape=1) +
  geom_abline() +
  geom_errorbar(aes(ymin=Success_Rate-sd, ymax=Success_Rate+sd)) +
  ggtitle("Bar plot with standard deviation as error bars")


dev.off()



#########################################Random in surrogate matrix

#gx2 <- c()
#gy2 <- c()
#foo = seq(1, 100, by=2.5)
#for (i in 1:99){
 # x <- graphx(x.ad, x2,y2,xbar1, xbar2, ybar1,ybar1p, MaxNumColors,warmcolor, coldcolor, i/1.0)
#  softimpute(seednumber)
 # b=graph(x.ad, "x2.csv", x2, y2,xbar1,xbar2,ybar1, ybar1p,MaxNumColors,warmcolor, coldcolor)
#  gx2 <- c(gx2, 100.0-(i*1.0))
 # gy2 <- c(gy2, successrate(x, a, b))
#}
#par(fig=c(0,x2,0,y2), new=FALSE)
#plot(gx2, gy2, pch = 5, cex = 1.3, col = "blue", main = "SVD- Random vs. Success Rate.", xlab = "Random (%)", ylab = "Success Rate (%)")
#abline(lm(gy2 ~ gx2))


#dev.off()


###########################################Others methods than SVD
#seednumber <- Sys.time() 
#seednumber <- gsub(':','',seednumber) 
#seednumber <-gsub('-','',seednumber) 
#seednumber <-gsub(' ','',seednumber) 
#seednumber <- gsub('2016','',seednumber) 
#seednumber <- gsub('0','',seednumber)

#softimpute(seednumber)
#b=graph(x.ad, "x2.csv", x2, y2,xbar1,xbar2,ybar1p,MaxNumColors,warmcolor, coldcolor)

#softimpute2(seednumber)
#c=graph(x.ad, "x3.csv", x2, y2,xbar1,xbar2,ybar1p,MaxNumColors,warmcolor, coldcolor)

####################PYTHON using numpy and pandas: MEAN
#python.load("Numpy.py")
#d=graph(x.ad, "x4.csv", x2, y2,xbar1,xbar2,ybar1p,MaxNumColors,warmcolor, coldcolor)

####################PYTHON using numpy and pandas: MEDIAN
#python.load("Numpy2.py")
#e=graph(x.ad, "x5.csv", x2, y2,xbar1,xbar2,ybar1p,MaxNumColors,warmcolor, coldcolor)

############
########PYTHON using numpy and pandas: MOST FREQUENT
#python.load("Numpy3.py")
#f=graph(x.ad, "x6.csv", x2, y2,xbar1,xbar2,ybar1p,MaxNumColors,warmcolor, coldcolor)

##########################################Calc

#b2<-mean((a[is.na(x)]-b[is.na(x)])^2)
#c2<-mean((a[is.na(x)]-c[is.na(x)])^2)
#d2<-mean((a[is.na(x)]-d[is.na(x)])^2)
#e2<-mean((a[is.na(x)]-e[is.na(x)])^2)
#f2<-mean((a[is.na(x)]-f[is.na(x)])^2)

#paste("[Average of SI SVD:",b2,sep=" ") 
#paste("[Average of SI ALS:",c2,sep=" ") 
#paste("[Average of PANDA Mean:",d2,sep=" ") 
#paste("[Average of PANDA Median:",e2,sep=" ") 
#paste("[Average of PANDA MostFrequent:",f2,sep=" ") 

#############################Analyzing SVD


#successrate(x, a, b)


###########

#dev.off()  

#Next tasks:
#1. error bars on success rate
#2. change success rate criteria - standard deviation 



