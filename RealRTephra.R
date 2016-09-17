library(ggplot2)



## Read in Calibration data
IntCal13 <- read.table("http://www.radiocarbon.org/IntCal13%20files/intcal13.14c",
                       skip=8,sep=",")[,c(1,2,5)]

# Linear Interpolation of 
CalC <- approx(IntCal13$V1, IntCal13$V2, n=50000, xout=c(1:50000))
CalS <- approx(IntCal13$V1, IntCal13$V5, n=50000, xout=c(1:50000))

Reimer13 <- data.frame("Year"=CalC$x, "C14 Age"=CalC$y, "Sigma"=CalS$y)


# posterior density function for Mean year
# Numerical integration in 1 dimension
post <- function (carb, sigma ,data = Reimer13){
    year <- 1:50000
    prior <- rep(50000^-1, 50000)
    poster <- dnorm(carb, data$C14.Age, sigma) * prior
    poster/sum(poster)
}


##### Turn single MCMC output into Distribution ####
# create posterior dist'n from MCMC output
postMCMC <- function(vector){
    out <- seq(1:50000)
    for(i in 1:50000){
        out[i] <- length(vector[which(vector==i)])
    }
    out/sum(out)
}

####### Main MCMC Functionss ###
# Given (single) carb date and previous calendar year,
# chooses the next calendar year value
# carb is a carbon date
# init is the calendar year we're starting with
# lowBound and upBound are the years we have to stay between (order constraints)
# prior is the prior distribution
# var is the variance of the std dev of the proposal density
genUnknownLag <- function(carb,init,sigma, prior, sd = 25,lowBound=0, upBound=50000, data = Reimer13){
    reject <- 0

    carbpost <- dnorm(carb,data$C14.Age[init],sigma)*prior[init]
    # draw new year
    newYear1 <- rnorm(1,init, sd)
    # make sure it's between appropriate bounds
    #    print(lowBound)
       #  print(upBound)
         #print(init)
    while (newYear1 > (upBound-2) | newYear1 < (lowBound+2) ) {
        newYear1 <- rnorm(1, init, sd)
    }
    newYear1 <- round(newYear1)

    newPost1 <- dnorm(carb,data$C14.Age[newYear1],sigma)*prior[newYear1]
    # see if newYear1 has a higher probability than init
    if (newPost1/carbpost >= 1) {
        init <- newYear1
    } else { 
        draw <- runif(1,0,1)
        #accept/ reject
        if (newPost1/carbpost >= draw) {
            init <- newYear1
        } else reject <- 1
    }
    # returns a vector of the next calendar year and 0/1 for whether or not the update was accepted
    df <- c(init,reject)
    df
}


# Keeps track of the Markov Chain
# carb = vector of carbon dates
# sigma = vector of sigma values associated with the carbon dates
# inits = vector of initial calendar year guesses associated with carb dates
# iterations = length of the markov chain
# This function is written assuming we're dating the R Tephra, so nothing about the input 
# specifies we're dating a sample without a carbon date
unknownLagR <- function(carb, sigma, inits, iterations = 50000, prior = rep(50000^-1,50000)) {
    reject <- 0
    len <- length(carb)
    # so we're storing the results of MCMC in a matrix where every column corresponds to a sample
    # this matrix won't hold the (undated) tephra sample
    results <- matrix(0, nrow = iterations, ncol = len) 
    # first row is our initial calendar year guesses
    results[1,] <- inits
    # keep a vector of boundaries between samples
    bounds <- rep(0,len-1) 
    # separate vector for the R Tephra
    tephra <- rep(0,iterations)
    for (i in 2:iterations){    # loops over iterations
        print(i)                    ## this line prints out the current iteration, could be commented
        for (j in 1:(len-1)){   
            # calculate the boundaries for the dated sample for the current iteration
            # bounds[1] is the boundary b/w samples 1 and 2
            bounds[j] <- ( results[i-1,j+1] + results[i-1,j] ) / 2
        }
        
        # Reset of the boundaries to be the R Tephra calendar year guess
        bounds[2] <- round(runif(1, results[i-1,2]+2, results[i-1,3]-2))
        # and then save it to the tephra vector
        tephra[i-1] <- bounds [2]

        # loop over the carbon dated samples to find next iteration calendar year values
        for (j in 1:len){ 
            if (j == 1 ){
                k <- genUnknownLag( carb[j], init = results[i-1,j], sigma[j],
                                              upBound = bounds[j], prior = prior)
            } else if (j == 5){
                k <- genUnknownLag(carb[j], init = results[i-1,j], sigma[j],
                                              lowBound = bounds[j-1],prior = prior)
            } else {
                k <- genUnknownLag(carb[j], init = results[i-1,j], sigma[j],
                                              lowBound = bounds[j-1],  upBound = bounds[j], 
                                              prior = prior)
            }
            
            # update results matrix and rejection count
            reject <- reject + k[2]
            results[i,j] <- k[1]
 
            
        }
    }
    print( reject / (5 * (iterations-1)) ) # prints out the rejection rate for the MC
    # return a data.frame containing the results matrix and the R tephra vector
    out <- data.frame(results,tephra)
    out
}


### Calculate Highest Posterior Density region ###
hdpRegion <- function(density, CI = .95){
    uniroot( function(density,CI,h) {sum(density[which(density>h)]) - CI},
             c(0,1), density = density, CI = CI)$root
}



################ R Tephra Analysis ################

# The radiocarbon data
carbons <- c(8905,8920,8760,8890,8990)
sigmas <- c(20,60,80,40,60)

iters <- 50000 

set.seed(96)
ptm <- proc.time()
rTephra <- unknownLagR(carbons, sigmas, c(9704,9739,9755,9771,10250), iterations = iters)
proc.time() - ptm


tephra <- postMCMC(rTephra$tephra)

# this thins the chain so we're not looking at all the values
p <- seq(10000, iters, by = 5)

# Turn chains into densities
a <- postMCMC( rTephra[p,1] )
b <- postMCMC( rTephra[p,2] )
c <- postMCMC( rTephra[p,3] )
d <- postMCMC( rTephra[p,4] )
f <- postMCMC( rTephra[p,5] )

aa <- post( carbons[1], sigmas[1] )
bb <- post( carbons[2], sigmas[2] )
cc <- post( carbons[3], sigmas[3] )
dd <- post( carbons[4], sigmas[4] )
ff <- post( carbons[5], sigmas[5] )
gg <- rep(NA, 50000)


q <- 50000

#throw it all into a data.frame
dataF <- data.frame( Data = c(a,b,tephra,c,d,f), Factor = c( rep("1", q), rep("2",q), rep("R Tephra",q) ,
                     rep("3",q), rep("4",q), rep("5",q)), Indiv = c(aa,bb,gg,cc,dd,ff), 
                     pindex = rep(1:q, 6) ) 

# Recode sample names as a factor variable, to get the order correct
dataF$Factor <- factor(dataF$Factor,levels=c("1","2","R Tephra","3","4","5"))

# pdf("tephStack.pdf", width = 6.25, height = 4.5)
ggplot(data = dataF) + geom_area(aes (x = pindex, y = Data), fill = "peru") + 
    geom_line(aes(x = pindex, y = Indiv)) + facet_grid(Factor ~ .) + 
    xlim(9700,10400) + ylab("Density") + xlab("Modeled Age (BP)") #+ theme_bw()
#dev.off()



# R Tephra 95% confidence region
teph95 <- hdpRegion(tephra, CI=.95)
tephconf <- qplot(x = 1:50000, y = tephra, geom="line") + xlim(c(9950,10150)) +# theme_bw() + 
    geom_ribbon(aes(ymax = ifelse(tephra > teph95, tephra,0), 
                    x = c(1:50000)), ymin = 0, fill="peru", alpha=0.6) +
    xlab("Age (Years BP)") + ylab("Density")  + 
    geom_vline(xintercept=mean(rTephra[p,6])) + 
    geom_vline(xintercept=median(rTephra[p,6]), linetype="dashed") 

#pdf("conTeph.pdf", height = 2.5, width = 5)
    tephconf
#dev.off()


