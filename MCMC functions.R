setwd("/Users/reillyvillanueva/Documents/Thesis/Data")

##### Read in calibration data + make data frame #####


    # read the data into R, but only the columns we need
    reimer09Data <- read.table("intcal09.14c", header=TRUE, skip = 9,sep=",")[c(1:3521),c(1,2,5)]
    
    #linear interpolation of data wrt C14 age
    CalC <- approx(reimer09Data$CAL.BP, reimer09Data$X14C.age, n=50000, xout=c(1:50000))
    plot(CalC, ylab="C-14 Age", xlab="Years BP",type="l")
    
    
    #linear interpolation wrt Sigma
    CalS <- approx(reimer09Data$CAL.BP, reimer09Data$Sigma, n=50000, xout=c(1:50000))
    plot(CalS, ylab="Sigma", xlab="Years BP",type="l")
    
    #put the values I care about into a data frame
    calC <- data.frame("Year"=CalC$x, "C14 Age"=CalC$y, "Sigma"=CalS$y)


##### Turn single MCMC output into Distribution ####
# create posterior dist'n from MCMC output
postMCMC <- function(vector){
    out <- seq(1:50000)
    for(i in 1:50000){
        out[i] <- length(vector[which(vector==i)])
    }
    out/sum(out)
}

##### Finding out how much space between consecutive years ####
# crawfish <- rep(0,49999)
# for (i in seq_along(crawfish)){
#     crawfish[i] <- calC$C14.Age[i] - calC$C14.Age[i+1]
# }
# min(crawfish)   # -8.8
# max(crawfish)    # 6.2
# mean(crawfish)    # -0.93


###### Pick a single initial year given Carbon Data ####
singInit <- function(carb){
    # take all calendar years in ball centered on carb w/ radius 5
    calC$Year[which( carb - 5 < calC$C14.Age & calC$C14.Age < carb + 5 )]
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
            inits[i] <- min(years[ which(years > inits[i-1]) ] )
        }
    }
    # return the vector of years
    inits
}

##### Single Carbon Date w/ Uniform Prior #####
# MCMC with uniform prior
mcmc09 <- function(carb, iterations){
    reject <- 0
    init <- multInit(carb)
    #initial vector
    arr <- rep(0,iterations)
    # loop through iterations
    for (i in 1:iterations) {
        #create density for initial guess
        postOrg <- dnorm(carb,calC$C14.Age[init],calC$Sigma[init]) # *PRIOR
        # pick a new year from a normal dist around the inital guess
        newYear <- rnorm(1, init, 100)
        # make sure the new year is positive integer & in the support
        # is this really the best way to do this?????
        while (newYear > 500000 | newYear < 0){
            newYear <- rnorm(1, init, 100)
        }
        newYear <- round(newYear)
        # create density for proposed new year
        postNew <- dnorm(carb, calC$C14.Age[newYear],calC$Sigma[newYear])  # *PRIOR
        # random draw to accept/reject with
        ###### use if statement to test if greater than 1 first ######
        if (postNew/postOrg >= 1) {
            init <- newYear
        } else {
            #accept/ reject
            draw <- runif(1,0,1)
            if (postNew/postOrg >= draw){
                init <- newYear
            } else reject = reject + 1
        }
        # put the year into the appropriate part of vector
        arr[i] <- init
    }
    # output vector
    acceptanceRate <- (iterations-reject)/iterations
    output <- data.frame(arr,acceptanceRate)
    output
}


###### 2 carbon dates w/ known lag b/w them #####
# MCMC Code for 2 Carbon Dates with known lag
# init is an initial guess for carb year 
lagMCMC <- function(carb, carb2, init, lag, iterations){
    reject <- 0
    #prior <- 50000^-1
    #initial vector
    states <- rep(0,iterations)
    # loop through iterations
    for (i in 1:iterations) {
        # create density for initial guess
        # assumes lag is positive value by putting lag-1 zeros to end of carb2's density
        postOrg <- dnorm(carb,calC$C14.Age[init],calC$Sigma[init])*
            dnorm(carb2,calC$C14[init+lag],calC$Sigma[init+lag])
        # pick a new year from a normal dist around the inital guess
        newYear <- rnorm(1, init, 100)
        # make sure the new year is positive integer & in the support
        while (newYear > 500000-lag | newYear < 0-lag){
            newYear <- rnorm(1, init, 100)
        }
        newYear <- round(newYear)
        # create density for proposed new year
        postNew <- dnorm(carb,calC$C14.Age[newYear],calC$Sigma[newYear])*
            dnorm(carb2,calC$C14[newYear+lag],calC$Sigma[newYear+lag])
        # random draw to accept/reject with
        ###### use if statement to test if greater than 1 first ######
        if (postNew/postOrg >= 1) {
            init <- newYear
        } else {
            #accept/ reject
            draw <- runif(1,0,1)
            if (postNew/postOrg >= draw){
                init <- newYear
            } else reject = reject + 1
        }
        # put the year into the appropriate part of vector
        states[i] <- init
    }
    # output vector
    acceptanceRate <- (iterations-reject)/iterations
    output <- data.frame(states)
    print(acceptanceRate)
    output
}

##### Meta-Comment : This really feels like old code, there's cleaner code later 
#####                that does the same thing in half the lines
#####                it's cool.
#####                There are some tests of this function anyway
# MCMC function for two carbon dates with unknown lag time
# Assuming carb < carb2
# unknownLag <- function(carb,carb2,init1,init2,iterations){
#     states <- rep(0,iterations)
#     states2 <- rep(0,iterations)
#     reject <- 0
#     reject2 <- 0
#     for(i in 1:iterations){
#         carbpost <- dnorm(carb,calC$C14.Age[init1],calC$Sigma[init1])
#         newYear1 <- rnorm(1,init1,100)
#         while(newYear1 > init2 | newYear1 < 0){
#             newYear1 <- rnorm(1, init1, 100)
#         }
#         newYear1 <- round(newYear1)
#         newPost1 <- dnorm(carb,calC$C14.Age[newYear1],calC$Sigma[newYear1])
#         if (newPost1/carbpost >= 1) {
#             init1 <- newYear1
#         } else {
#             #accept/ reject
#             draw <- runif(1,0,1)
#             if (newPost1/carbpost >= draw){
#                 init1 <- newYear1
#             } else reject = reject + 1
#         }
#         states[i] <- init1
#         ### Start the carb2 part
#         carb2post <- dnorm(carb2,calC$C14.Age[init2],calC$Sigma[init2])
#         newYear2 <- rnorm(1,init2,100)
#         while(newYear2 > 50000 | newYear2 < init1){
#             newYear2 <- rnorm(1, init2, 100)
#         }
#         newYear2 <- round(newYear2)
#         newPost2 <- dnorm(carb2,calC$C14.Age[newYear2],calC$Sigma[newYear2])
#         if (newPost2/carb2post >= 1) {
#             init2 <- newYear2
#         } else {
#             #accept/ reject
#             draw <- runif(1,0,1)
#             if (newPost2/carb2post >= draw){
#                 init2 <- newYear2
#             } else reject2 = reject2 + 1
#         }
#         states2[i] <- init2
#     }
#     state <- data.frame(states,states2)
#     state
# }


###### Go forward one iteration ######
# Given carb date, moves forward one iteration
# carb is a carbon date
# init is the calendar year we're starting with
# lowBound and upBound are the years we have to stay between
# prior is the prior distribution, defaults to uniform

genUnknownLag <- function(carb,init, sigma, lowBound=0,upBound=50000, prior = rep(50000^-1,50000)){
    carbpost <- dnorm(carb,calC$C14.Age[init],sigma)*prior[init]
    # draw new year
    newYear1 <- rnorm(1,init,10)
    # make sure it's between appropriate bounds
    while (newYear1 > upBound | newYear1 < lowBound){
        newYear1 <- rnorm(1, init, 100)
    }
    newYear1 <- round(newYear1)
    newPost1 <- dnorm(carb,calC$C14.Age[newYear1],sigma)*prior[newYear1]
    # see if newYear1 has a higher probability than init
    if ( identical(newPost1, numeric(0)) ) {
        newPost1 <- 0
    } else if (is.na(newPost1) ) {
        newPost1 <- 0
    }
    print(carbpost)
    print(newPost1)
    if (newPost1/carbpost >= 1) {
        init <- newYear1
    } else {
        #accept/ reject
        draw <- runif(1, 0, 1)
        if (newPost1/carbpost >= draw){
            init <- newYear1
        }
    }
    init
}


######## Procedure for 2 Trees ############
# Takes two trees, two lag values
# iterations & init
# outputs 
twoTrees <- function(tree1,lag1,tree2,lag2,iterations,init, prior){
    priors <- matrix(nrow=50000,ncol=2)
    priors[,1] <- lagCar(tree1,lag1)
    priors[,2] <- lagCar(tree2,lag2)
    # Have priors. Numerical Integration b/w trees
    min <- min(calC$Year[which(priors[,1]>1E-5)])
    max <- max(calC$Year[which(priors[,2]>1E-5)])
    
    # initialize a matrix to hold the vectors
    posts <- matrix(0, nrow=max-min+1, ncol=max-min+1)
    # loop over all max(tree1) values
    prior <- (.5*(max-min)^2)^-1
    for (i in min:(max-1)){
        #calculation of posterior for current carb1
        inorm <- dnorm(tree1[which.max(tree1)],calC$C14.Age[i],calC$Sigma[i])
        # find appropriate number of zeros (spacing out the matrix)
        zeros <- rep(0,i-min+1)
        # calculate posterior for carb1=i and carb2= carb1+1:max and put into matrix
        posts[i-min+1,] <- c(zeros,inorm*
                                 dnorm(tree2[which.max(tree2)],
                                       calC$C14.Age[(i+1):max], 
                                       calC$Sigma[(i+1):max]) * prior) 
    }
    # return the posterior matrix
    postMaT <- list(posts,min,max)
}


####### Unknown Lag for multiple densities #######
# carb: vector of c14 values least to greatest
# sigma: vector of error values
unknownLagMult <- function(carb, sigma, iterations = 50000, prior = rep(50000^-1,5000)) {
    len <- length(carb)
    results <- matrix(0, nrow = iterations, ncol = len)
    results[1,] <- multInit(carb)
    for (i in 2:iterations){
        for (j in 1:len){
            if (j == 1 ){
                results[i,j] <- genUnknownLag(carb[j], init = results[i-1,j], sigma[j],
                                upBound = results[i-1,j+1], prior = prior)
            } else if (j == len){
                results[i,j] <- genUnknownLag(carb[j], init = results[i-1,j], sigma[j],
                                lowBound = results[i-1,j-1],prior = prior)
            } else {
                results[i,j] <- genUnknownLag(carb[j], init = results[i-1,j], sigma[j],
                                lowBound = results[i-1,j-1],  upBound = results[i-1,j+1], 
                                prior = prior)
            }
        }
    }
    results
}


######## Procedure for multiple carbon dates #####
# n trees
# carb is a vector of ALL carbon dates in the sample from least to greatest within the ring
# tree is a vector giving which tree each date is associated with
# lag is a vector giving the lag from the oldest ring. All lag values are positive, decreasing
# Should basically be a data frame so ith entry of each refers to same piece of data
# init is a vector of initial values for last date in each tree
multCarb <- function(carb,tree,lag,iterations=50000,init){
    # see how many trees there are
    nameTrees <- unique(tree)
    # nubmber of unique trees
    numTrees <- length(nameTrees)
    # initialize vectors to hold indices of beginning and end of each tree
    beginTrees <- rep(0,numTrees)
    endTrees <- rep(0,numTrees)
    # loop to fill those vectors
    for (i in 1:numTrees){
        beginTrees[i] <- min(which(tree==nameTrees[i]))
        endTrees[i] <- max(which(tree==nameTrees[i]))
    }
    # make a matrix to put the posteriors in
    priors <- matrix(0, nrow=50000, ncol=numTrees)
    # loop to get the posterior dist'n of the known lag trees
    for (i in (1:numTrees)){
        start <- beginTrees[i]
        end <- endTrees[i]
        priors[,i] <- lagCar(carb[start:end],lag[start:end])
    }
    
    # Great. Now we only have to work with one date for each tree
    # Use the new distributions to do unknown lag posterior estimation
    # Going to use these lags as a prior (weighting the dnorm results)
    
    # Initialize a matrix to put all the MCMC output
    post <- matrix(0,ncol=numTrees,nrow=iterations+1)
    # put in initial values in first row
    post[1,] <- init
    # set up boundary conditions on the first and last columns
    for (i in (2:(iterations+1))){
        for (j in (1:numTrees)){
            # Arguements for genUnknownLag: carb, init, lowBound, upBound, prior
            # special cases for first and last tree:
            carbIndex <- endTrees[j]
            # middle tree case
            if (j != 1 & j != numTrees){
                post[i,j] <- genUnknownLag(carb[carbIndex],post[i-1,j],post[i,j-1],post[i-1,j+1],priors[,j])
            } else if (j==1){
                post[i,j] <- genUnknownLag(carb[carbIndex],post[i-1,j],0,post[i-1,j+1],priors[,j])
            } else post[i,j] <- genUnknownLag(carb[carbIndex],post[i-1,j],post[i,j-1],50000,priors[,j])
        }
    }
    # return the matrix
    post
}