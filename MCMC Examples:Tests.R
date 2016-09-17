source('/Users/reillyvillanueva/Documents/Thesis/Data/MCMC functions.R')

##### Find all calendar years with C14 around given value #######
singInit(5000)

#### Find vector of initial calendar years for given C14 values #####
multInit(rep(5000,5))

###### mcmc09 : Single Carbon Datew/ Uniform Prior ####
    # Example with year 30,000
    # 50,000 iterations
    carbTest <- mcmc09(30000, 50000)
    # plots it
    plot(carbTest$arr[calC$Year > 34000 & calC$Year < 36000] ~
             calC$Year[calC$Year>34000 & calC$Year<36000],type="l")
    plot(post(30000)[calC$Year > 34000&calC$Year<36000] ~ 
             calC$Year[calC$Year>34000&calC$Year<36000],type="l")
    
    lines(density(carbTest$arr),col="red")


##### lagMCMC : 2 carbon dates w/ known lag b/w them #####
## Example with 1st carbon date 10,060
# 2nd carbon date 10,090
# intitial value of 11,600 for first carbon date
# 20 year lag between, from 1st carbon
# with 50,000 iterations
lagTest <- lagMCMC(10060,10090,11600,20,50000)
# turn mcmc output into distribution
lag <- postMCMC(lagTest$states)
# plot
plot(lagTest$states,type="l")


##### Old Code, Better Code Later 
# # Example of unknownLag
# ## First carbon date is 16,646
# ## initial date 19,700
# ## Second carbon is 33,120
# ## initial date 37,600
# # 50,000 iterations
# hamburger <- unknownLag(16645,33120,19700,37600,50000)
# # turn mcmc output from first date into distribution 
# ketchup <- postMCMC(hamburger$states)
# # turn mcmc output from second date into distribution
# mustard <- postMCMC(hamburger$states2)
# # turn distributions into density plot
# mayo <- kde2d(hamburger$states,hamburger$states2,h=c(1,1),n=25)
# # look at the density plot, (cc is a color pallete made earlier)
# image(mayo,col=cc)
# # contour(mayo) #uncomment for contour plot

###### genUnknownLag : Go forward one iteration ######
genUnknownLag(25000,29650)

####### twoTrees : Procedure for 2 Trees #######
b1<-c(1201,1200,1189,1193)
l1<-c(23,15,3,1)
b2<-c(1259,1263,1272,1274)
l2<-c(27,12,4,1)
lobster <- lagCar(b2,l2)
plot(lobster[calC$Year<1500 & calC$Year>1000]~calC$Year[calC$Year<1500 & calC$Year>1000], type="l")
tray <- twoTrees(b1,l1,b2,l2,30000,c(1130,1265))
contour(tray[[1]], axes=FALSE)
axis(1,at=c(0,1), labels=c(tray[[2]],tray[[3]]))
axis(2,at=c(0,1), labels=c(tray[[2]],tray[[3]]))

######## unknownLagMult  : no known lag values #### 
b1<-c(1201,1200,1189,1193)
unknTest <- unknownLagMult(b1, rep(7000,4), iterations = 100000)
aa <- postMCMC(unknTest[,1])
plot(aa[1:4000], type='l')
ab <- postMCMC(unknTest[,2])
plot(ab[1:4000], type='l')
ac <- postMCMC(unknTest[,3])
plot(ac[1:4000], type='l')
ad <- postMCMC(unknTest[,4])
plot(ad[1:4000], type='l')

###### multCarb : Procedure for multiple carbon dates #####
carbonDates <- c(2460,2530,2490)
treeNum <-c(1,2,3)
treeLag <- c(0,0,0)
initial <- c(2675,2685,2695)
treeTest <- multCarb(carbonDates,treeNum,treeLag,10000,initial)

## Another multCarb Example
# begin test of n-tree procedure:
carbonDates <- c(b1,b2,1270,1250,1258,1265,1270)
treeNames <- c(rep(1,4),rep(2,4),rep(3,5))
treeLags <- c(l1,l2,49,23,9,6,1)
initialGuesses <- c(1120,1190,1200)

testRun <- multCarb(carbonDates,treeNames,treeLags,50000,initialGuesses)
tree1 <- postMCMC(testRun[,1])
tree2 <- postMCMC(testRun[,2])
tree3 <- postMCMC(testRun[,3])
plot(tree1[calC$Year<1300&calC$Year>1000]~calC$Year[calC$Year<1300&calC$Year>1000],type="l",main='Tree 1')
plot(tree2[calC$Year<1250&calC$Year>1150]~calC$Year[calC$Year<1250&calC$Year>1150],type="l",main='Tree 2')
plot(tree3[calC$Year<1250&calC$Year>1150]~calC$Year[calC$Year<1250&calC$Year>1150],type="l",main='Tree 3')
lines(tree2,col="red")
plot(post(1274)[calC$Year<1250&calC$Year>1150]~calC$Year[calC$Year<1250&calC$Year>1150],type="l")
plot(post(1270)[calC$Year<1250&calC$Year>1150]~calC$Year[calC$Year<1250&calC$Year>1150],type="l")


runVar1 <- rep(0,50000)
for (i in 1:50000){
    runMean1[i] <- mean(testRun[1:i,1])
    runMean2[i] <- mean(testRun[1:i,2])
    runMean3[i] <- mean(testRun[1:i,3])
}

for (i in 1:50000){
    runVar1[i] <- var(testRun[1:i,1])
    runVar2[i] <- var(testRun[1:i,2])
    runVar3[i] <- var(testRun[1:i,3])
}

plot(runMean1~seq(1,50000))
plot(runMean2~seq(1,50000))
plot(runMean3~seq(1,50000))

plot(runVar1~seq(1,50000))
plot(runVar2~seq(1,50000))
plot(runVar3~seq(1,50000))




