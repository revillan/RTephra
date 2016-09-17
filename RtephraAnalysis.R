############
# This is a replication of the data analysis in "Geochemical characterization and dating 
# R tephra, a postglacial marker bed in Mount Rainier National Park, Washington, USA"
# 2016, Samolcyzk, Wallace, Cubley, Osborn, Clark
#
# Specifically, replicating Fig. 4
############

library(ggplot2)
library(nlme)

## Read in Calibration data
IntCal13 <- read.table("http://www.radiocarbon.org/IntCal13%20files/intcal13.14c",
                       skip=8,sep=",")[,c(1,2,5)]

names(IntCal13) <- c("Year","V2","V5")

# Linear Interpolation of 
CalC <- approx(IntCal13$Year, IntCal13$V2, n=50000, xout=c(1:50000))
CalS <- approx(IntCal13$Year, IntCal13$V5, n=50000, xout=c(1:50000))

Reimer13 <- data.frame("Year"=CalC$x, "C14 Age"=CalC$y, "Sigma"=CalS$y)


## posterior density function for Mean year
## Numerical integration in 1 dimension
post <- function (carb, sigma, data = Reimer13){
    year <- 1:50000
    prior <- rep(50000^-1, 50000)
    poster <- dnorm(carb,data$C14.Age, sigma) * prior
    poster/sum(poster)
}

##### Turn single MCMC output into Distribution ####
# create posterior dist'n from MCMC output
postMCMC <- function(vector, data = Reimer13){
    out <- rep(0 , length(data[[1]]) )
    for (i in (min(vector)):max(vector)) {
        out[i] <- length(vector[which(vector==i)])
    }
    out/sum(out)
}

###### Pick a single initial year given Carbon Data ####
singInit <- function(carb, data = Reimer13){
    # take all calendar years in ball centered on carb w/ radius 5
    data$Year[which( carb - 150 < data$C14.Age & data$C14.Age < carb + 150 )]
}


#### Pick multiple initial years where year i < i+1 ######
multInit <- function(carb){
    # initilize vector
    inits <- rep( 0, length(carb) )
    for (i in seq_along(carb) ) {
        # take smallest year for first element
        if (i == 1){
            inits[i] <- min( singInit( carb[i] ) )
            # take largest year for last element    
        } else if (i == length(carb) ) {
            inits[i] <- max( singInit( carb[i] ) )
            # otherwise take smallest year greater than previous year
        } else { 
            years <- singInit(carb[i])
            inits[i] <- min( years[ which(years > inits[i-1] + 15) ])
        }
    }
    # return the vector of years
    inits
}



# Given carb datem chooses
# carb is a carbon date
# init is the calendar year we're starting with
# lowBound and upBound are the years we have to stay between
# prior is the prior distribution
genUnknownLag <- function(carb,init,sigma, prior, lowBound=0,upBound=50000, data = Reimer13){
    reject <- 0
    carbpost <- dnorm(carb,data$C14.Age[init],sigma)*prior[init]
    # draw new year
    newYear1 <- rnorm(1,init,40)
    # make sure it's between appropriate bounds
#     print(lowBound)
#     print(upBound)
#     print(init)
#     print(newYear1)
    while (newYear1 > (upBound-2) | newYear1 < (lowBound+2) ) {
        newYear1 <- rnorm(1, init, 40)
    }
    newYear1 <- round(newYear1)
    newPost1 <- dnorm(carb,data$C14.Age[newYear1],sigma)*prior[newYear1]
    # see if newYear1 has a higher probability than inti
    if (newPost1/carbpost >= 1) {
        init <- newYear1
    } else { 
            draw <- runif(1,0,1)
            #accept/ reject
            if (newPost1/carbpost >= draw) {
                init <- newYear1
            } else reject <- 1
    }
    df <- c(init,reject)
    df
}


####### Unknown Lag for multiple densities #######
# carb: vector of c14 values least to greatest
# sigma: vector of error values
unknownLagMult <- function(carb, sigma, inits, iterations = 50000, prior = rep(50000^-1,50000)) {
    reject <- 0
    len <- length(carb)
    results <- matrix(0, nrow = iterations, ncol = len)
    results[1,] <- inits
    bounds <- rep(0,len-1) 
    for (i in 2:iterations){
        for (j in 1:(len-1)){
            bounds[j] <- ( results[i-1,j+1] + results[i-1,j] ) / 2
        }
        for (j in 1:len){
            if (j == 1 ){
                k <- genUnknownLag( carb[j], init = results[i-1,j], sigma[j],
                                              upBound = bounds[j], prior = prior)
            } else if (j == len){
                k <- genUnknownLag(carb[j], init = results[i-1,j], sigma[j],
                                              lowBound = bounds[j-1],prior = prior)
            } else {
                k <- genUnknownLag(carb[j], init = results[i-1,j], sigma[j],
                                              lowBound = bounds[j-1],  upBound = bounds[j], 
                                              prior = prior)
            }
            reject <- reject + k[2]
            results[i,j] <- k[1]
        }
    }
    print( reject / (5 * (iterations-1)) )
    results
}




######## Analyzing the data for Rtephra ########
carbons <- c(8905,8920,8760,8890,8990)
sigmas <- c(20,60,80,40,60)
initValues <- multInit(carbons)

set.seed(96)
ptm <- proc.time()
stratRes <- unknownLagMult(carbons, sigmas, initValues, iterations = 5000000)
proc.time() - ptm

p <- seq(15000, 5000000, by = 15)



acf2 <- acf(stratRes[,5], plot = FALSE)
acf22 <- with(acf2, data.frame(lag, acf))
uo <- ggplot(data = acf22, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_bw() +
    labs(x = "Lag", y = "ACF")

acf2 <- acf(stratRes[p,5], plot = FALSE)
acf22 <- with(acf2, data.frame(lag, acf))
u <- ggplot(data = acf22, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + theme_bw() +
    labs(x = "Lag", y = "ACF")

pdf("stratLagFive.pdf", width=6, height = 2)
grid.arrange(uo, u, ncol=2, top = "Sample 5")
dev.off()


acfPlot(stratRes[,1], labels = FALSE)
title(main = "Sample 1", xlab = "Lag", ylab= "ACF")
acfPlot(stratRes[p,1], labels = FALSE)
title( "Sample 1", xlab = "Lag", ylab= "ACF")
acfPlot(stratRes[,2],labels = FALSE)
title(main = "Sample 2", xlab = "Lag", ylab= "ACF")
acfPlot(stratRes[p,2],labels = FALSE)
title(main = "Sample 2", xlab = "Lag", ylab= "ACF")
acfPlot(stratRes[,3], labels = FALSE)
title(main = "Sample 3", xlab = "Lag", ylab= "ACF")
acfPlot(stratRes[p,3], labels = FALSE)
title(main = "Sample 3", xlab = "Lag", ylab= "ACF")
acfPlot(stratRes[,4], labels = FALSE)
title(main = "Sample 4", xlab = "Lag", ylab= "ACF")
acfPlot(stratRes[p,4], labels = FALSE)
title(main = "Sample 4", xlab = "Lag", ylab= "ACF")
acfPlot(stratRes[,5], labels = FALSE)
title(main = "Sample 5", xlab = "Lag", ylab= "ACF")
acfPlot(stratRes[p,5], labels = FALSE)
title(main = "Sample 5", xlab = "Lag", ylab= "ACF")
dev.off()

# Turn chains into densities
a <- postMCMC( stratRes[p,1] )
b <- postMCMC( stratRes[p,2] )
c <- postMCMC( stratRes[p,3] )
d <- postMCMC( stratRes[p,4] )
f <- postMCMC( stratRes[p,5] )


# The num. int. densities for each individual date
aa <- post( carbons[1], sigmas[1] )
bb <- post( carbons[2], sigmas[2] )
cc <- post( carbons[3], sigmas[3] )
dd <- post( carbons[4], sigmas[4] )
ff <- post( carbons[5], sigmas[5] )

q <- 50000

#throw it all into a data.frame
dataF <- data.frame( Data = c(a,b,c,d,f), Factor = c( rep("1", q), rep("2",q),
                rep("3",q), rep("4",q), rep("5",q)), Indiv = c(aa,bb,cc,dd,ff), 
                pindex = rep(1:q, 5) ) 

##### Stacked results plot with both MCMC and Num Int
pdf("Fig4Arep.pdf", height = 4, width = 6.5)
ggplot(data = dataF) + geom_area(aes (x = pindex, y = Data), fill = "peru") + 
    geom_line(aes(x = pindex, y = Indiv)) + facet_grid(Factor ~ .) + 
    xlim(9600,10600) + ylab("Density") + xlab("Modeled Age (BP)") #+ theme_bw()
dev.off()


##### Running Mean Plots ######
runlen <- seq(1,20000)

q <- 5000000
runMean <- data.frame(mc=as.vector(stratRes), Factor = c( rep("1", q), rep("2",q),
                        rep("3",q), rep("4",q), rep("5",q)), pindex = rep(1:q, 5) )

pdf("stackrunmean.pdf", height = 3.5, width=6)
ggplot(data = runMean) + geom_point(aes(y = cumsum(mc)/seq(along=pindex), x = pindex), size = 0.1) + 
    facet_grid(Factor ~ .) + xlim(c(1,20000)) + theme_bw() + 
     labs(x="Iterations",  y="Running Mean (Years BP)")
dev.off()

arun <- qplot(y= cumsum(stratRes[runlen,1])/seq(along=stratRes[runlen,1]),
      x = seq(along=stratRes[runlen,1]), geom="line") + ylab("Running Mean") + 
    xlab("Iterations") + theme_bw() + geom_line(y = mean(stratRes[,1]), col="peru")

brun <- qplot(y = cumsum(stratRes[runlen,2])/seq(along=stratRes[runlen,2]), 
      x =seq(along=stratRes[runlen,2]), geom="line") + ylab("Running Mean") + 
    xlab("Iterations") + theme_bw() + geom_line(y = mean(stratRes[,2], alpha = 0.6), col="peru") + 
    ylim(c(10000,10035))

pdf("runMeanB.pdf",width = 6, height = 1.75)
brun
dev.off()

crun <- qplot(y = cumsum(stratRes[runlen,3])/seq(along=stratRes[runlen,3]), 
      x = seq(along=stratRes[runlen,3]), geom="line") + ylab("Running Mean") +
    xlab("Iterations") + theme_bw() + geom_line(y = mean(stratRes[,3]), col="peru")

drun <- qplot(y = cumsum(stratRes[runlen,4])/seq(along=stratRes[runlen,4]), 
      x = seq(along=stratRes[runlen,4]), geom="line") + ylab("Running Mean") +
    xlab("Iterations") + theme_bw() + geom_line(y = mean(stratRes[,4]), col="peru")

frun <- qplot(y = cumsum(stratRes[runlen,5])/seq(along=stratRes[runlen,5]), 
      x = seq(along=stratRes[runlen,5]), geom="line") + ylab("Running Mean") +
      xlab("Iterations") + theme_bw() + geom_line(y = mean(stratRes[,5]), col="peru")

pdf("stratrun.pdf",height = 5, width = 6)
grid.arrange(arun,brun,crun,drun,frun,ncol=1)
dev.off()

hdpRegion <- function(density, CI = .95){
    uniroot( function(density,CI,h) {sum(density[which(density>h)]) - CI},
             c(0,1), density = density, CI = CI)$root
}


a95 <- hdpRegion(a, CI=.954)
aconf <- qplot(x = 1:50000, y = a, geom="line") + xlim(c(9850,10050)) + theme_bw() + 
    geom_ribbon(aes(ymax = ifelse(a > a95, a,0), 
                    x = c(1:50000)), ymin = 0, fill="peru", alpha=0.6) +
    xlab("Age (Years BP)") + ylab("Density") + labs(title="Sample 1") + 
    geom_vline(xintercept=mean(stratRes[15000:5000000,1])) + 
    geom_vline(xintercept=median(stratRes[15000:5000000,1]), linetype="dashed") 


b95 <- hdpRegion(b, CI=.954)
bconf <- qplot(x = 1:50000, y = b, geom="line") + xlim(c(9900,10150)) + theme_bw() + 
    geom_ribbon(aes(ymax = ifelse(b > b95, b,0), 
                    x = c(1:50000)), ymin = 0, fill="peru", alpha=0.6) +
    xlab("Age (Years BP)") + ylab("Density") + labs(title="Sample 2") +
    geom_vline(xintercept=mean(stratRes[15000:5000000,2])) + 
    geom_vline(xintercept=median(stratRes[15000:5000000,2]),linetype="dashed") 

c95 <- hdpRegion(c, CI=.954)
cconf <- qplot(x = 1:50000, y = c, geom="line") + xlim(c(10000,10150)) + theme_bw() + 
    geom_ribbon(aes(ymax = ifelse(c > c95, c,0), 
                    x = c(1:50000)), ymin = 0, fill="peru", alpha=0.6) +
    xlab("Age (Years BP)") + ylab("Density") + labs(title="Sample 3") + 
    geom_vline(xintercept=mean(stratRes[15000:5000000,3])) + 
    geom_vline(xintercept=median(stratRes[15000:5000000,3]), linetype="dashed") 

d95 <- hdpRegion(d, CI=.954)
dconf <- qplot(x = 1:50000, y = d, geom="line") + xlim(c(10075,10200)) + theme_bw() + 
    geom_ribbon(aes(ymax = ifelse(d > d95, d,0), 
                    x = c(1:50000)), ymin = 0, fill="peru", alpha=0.6) +
    xlab("Age (Years BP)") + ylab("Density") + labs(title="Sample 4") + 
    geom_vline(xintercept=mean(stratRes[15000:5000000,4])) + 
    geom_vline(xintercept=median(stratRes[15000:5000000,4]), linetype="dashed") 

f95 <- hdpRegion(f, CI=.954)
fconf <- qplot(x = 1:50000, y = f, geom="line") + xlim(c(10100,10300)) + theme_bw() + 
    geom_ribbon(aes(ymax = ifelse(f > f95, f,0), 
                    x = c(1:50000)), ymin = 0, fill="peru", alpha=0.6) +
    xlab("Age (Years BP)") + ylab("Density") + labs(title="Sample 5") +     
    geom_vline(xintercept=mean(stratRes[15000:5000000,5])) +
    geom_vline(xintercept=median(stratRes[15000:5000000,5]), linetype="dashed") 

pdf("stratconf.pdf", height=4.5,width=6)
grid.arrange(aconf,bconf,cconf,dconf,fconf,ncol=2)
dev.off()





