
getwd()
setwd("/Users/reillyvillanueva/Documents/thesis/data")
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

####### Making Pretty Colors #######
# with RColorBrewer
# because red + yellow is gross
ccc <- colorRampPalette(rev(brewer.pal(9,'RdPu')))
cc <- ccc(44)

# 2004 --------------------------------------------------------------------
# read the data into R, but only the columns we need
reimer04Data <- read.csv("Reimer.dat")[,c("YearBP","C14age","Sigma")]

#linear interpolation of data wrt C14 age
calCurveC <- approx(reimer04Data$YearBP, reimer04Data$C14age, n=26001)
plot(calCurveC, ylab="C-14 Age", xlab="Years BP",type="l")

#linear interpolation wrt Sigma
calCurveS <- approx(reimer04Data$YearBP, reimer04Data$Sigma, n=26001)
plot(calCurveS, ylab="Sigma", xlab="Year BP",type="l")

#put the values I care about into a data frame
calCurve <- data.frame("Year"=calCurveC$x, "C14 Age"=calCurveC$y, "Sigma"=calCurveS$y)


# posterior density function of Mean year
# 2004
posterior <- function(carb){
    year <- 0:26000
    prior <- rep(26001^-1,26001)
    post <- dnorm(carb,calCurve$C14.Age,calCurve$Sigma)*prior
    post/sum(post)
}

# These two are identical for Carbon Age = 9500
#2004
plot(posterior(9500)[calCurve$Year<13000&calCurve$Year>9000] ~ 
         calCurve$Year[calCurve$Year<13000&calCurve$Year>9000],type="l")

#These two are pretty different. Borderline year just isn't really part of 2004 cal data
#2004
plot(posterior(22000)[calCurve$Year>9000] ~ 
         calCurve$Year[calCurve$Year>9000],type="l")

# MCMC with uniform prior
mcmc04 <- function(carb, iterations, init){
    # take initial year to index value by adding 1 (since yearBP starts at 0)
    init <- init + 1
    # initial vector
    arr <- rep(0,iterations)
    # loop through iterations
    for (i in 0:iterations) {
        #create density for initial guess
        postOrg <- dnorm(carb,calCurve$C14.Age[init],calCurve$Sigma[init])
        # pick a new year from a normal dist around the inital guess
        # is this sigma value ok
        newYear <- rnorm(1, init, 100)
        # make sure the new year is positive integer & in the support
        newYear <- round(abs(newYear))
        if(newYear > 26000){
            newYear <- 26000
        }
        # create density for proposed new year
        postNew <- dnorm(carb, calCurve$C14.Age[newYear+1],calCurve$Sigma[newYear+1])
        # random draw to accept/reject with
        draw <- runif(1,0,1)
        # accept/ reject
        if (postNew/postOrg >= draw){
            init <- newYear
        }
        # put the year into the appropriate part of vector
        arr[i] <- init + 1
    }
    # output vector
    arr
}

## Example for C14 = 350
plate <- mcmc04(350,1000000,375)
mean <- mean(plate)
plot(posterior(350)[calCurve$Year < 600] ~ calCurve$Year[calCurve$Year < 600],type="l")
abline(v=mean,col="red")

## Example for C14 = 20828
saucer <- mcmc04(20828,1000000, 20000)
mean <- mean(saucer[saucer!=26000])
plot(posterior(20828)[calCurve$Year > 20000] ~ calCurve$Year[calCurve$Year>20000],type="l")
abline(v=mean, col="red")

## Example for C14 = 6500
coffee <- mcmc04(6500,1000000,7428)
mean <- mean(coffee)
plot(posterior(6500)[calCurve$Year < 7500 & calCurve$Year > 7400] ~ 
         calCurve$Year[calCurve$Year < 7500 & calCurve$Year > 7400],type="l")
abline(v=mean, col="red")

## Example for C14 = 11683, this is a curvy part!
mug <- mcmc04(11683,1000000, 13510)
mean <- mean(mug)
stdDev <- sd(mug)
stdDev
plot(mug~c(1:1000000))
plot(posterior(11683)[calCurve$Year < 17000 & calCurve$Year > 13000] ~ 
         calCurve$Year[calCurve$Year < 17000 & calCurve$Year > 13000],type="l")
abline(v=mean,col="red")

# 2009 --------------------------------------------------------------------

#### Read in Calibration Data, Interpolate, Data.Frame ####

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

#posterior density function of Mean year
post <- function (carb){
    year <- 1:50000
    prior <- rep(50000^-1, 50000)
    poster <- dnorm(carb,calC$C14.Age, calC$Sigma)*prior
    poster/sum(poster)
}


# These two are identical for Carbon Age = 9500
plot(post(9500)[calC$Year<13000&calC$Year>9000] ~ 
         calC$Year[calC$Year<13000&calC$Year>9000],type="l")

#These two are pretty different. Borderline year just isn't really part of 2004 cal data
plot(post(22000)[calC$Year<40000&calC$Year>9000] ~ 
         calC$Year[calC$Year<40000&calC$Year>9000],type="l")

#just another example year: 20828
plot(posterior(20828)[calCurve$Year > 20000] ~ calCurve$Year[calCurve$Year>20000],type="l")

# MCMC with uniform prior
mcmc09 <- function(carb, iterations, init){
    reject <- 0
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

# Example with year 30,000
# 50,000 iterations
# and 35,000 as initial year of mcmc
carbTest <- mcmc09(30000, 50000, 35000)
# plots it
plot(carbTest$arr[calC$Year > 34000 & calC$Year < 36000] ~
         calC$Year[calC$Year>34000 & calC$Year<36000],type="l")
plot(post(30000)[calC$Year > 34000&calC$Year<36000] ~ 
         calC$Year[calC$Year>34000&calC$Year<36000],type="l")

lines(density(carbTest$arr),col="red")

# 2 carbon dates, 1 lag value
# Lag assumed to be positive, so carb < carb2
lagCarb <- function(carb,carb2,lag){
    prior <- rep(50000^-1, 50000) 
    ageL <- calC$C14.Age[lag:length(calC$C14.Age)]
    sigL <- calC$Sigma[lag:length(calC$Sigma)]
    poster <- dnorm(carb,calC$C14.Age, calC$Sigma)*c(dnorm(carb2,ageL,sigL),rep(0,lag-1))*prior
    poster/sum(poster)
}

# 1st carbon date: 10,060
# 2nd carbon date: 10,090
# 20 year lag
lagTest <- lagCarb(10060,10090,20)
# plot it
plot(lagTest[calC$Year>11000 & calC$Year<12000]
     ~ calC$Year[calC$Year>11000 & calC$Year<12000],type="l")
lines(post(10060), col="red")
lines(post(10090), col='blue')




## >2 carbon dates, carb ordered newest to oldest
## lag is always from first element of carb
lagCar <- function(carb,lag){
    #make lag same length as carb
    lag <- c(0,lag)
    prior <- rep(50000^-1,50000)
    len <- length(lag)
    # initialize posterior
    posterior <- dnorm(carb[1],calC$C14.Age,calC$Sigma)

    for (i in 2:len){
        #calculate zeros to ass to vector to make length = 50000
        lagAdd <- rep(0,lag[i]-1)
        #make C14 and Sigma vectors with appropriate vector length and correct subsets of cal Curve
        cAge <- c( calC$C14.Age[lag[i]:length(calC$C14.Age)], lagAdd)
        sig <- c( calC$Sigma[lag[i]:length(calC$Sigma)], lagAdd)
        # use subsetted C14 and Sigma values to calculate posterior for particular carbon age.
        # multiply it by current posterior value
        posterior <- posterior*dnorm(carb[i], cAge,sig)
    }
    posterior <- posterior*prior
    posterior/sum(posterior)
}

# carbon dates 10,060 and 10,090
# 20 year lag
multTest <- lagCar(c(10060,10090),20)
plot(multTest[calC$Year < 12000 & calC$Year > 11000] 
     ~ calC$Year[calC$Year < 12000 & calC$Year > 11000],type="l")

## Plot the result of lagCar vs. lagCarb
# they match
plot(lagTest[calC$Year>11000 & calC$Year<12000]
     ~ calC$Year[calC$Year>11000 & calC$Year<12000],type="l")
lines(multTest, col="red")

## Example with more carbon dates
# Lots of carbon dates
lagCar(c(6040,14000,16583,25522,38924),c(8960,11342,20400,34000))

# MCMC Code for 2 Carbon Dates with known lag
#init is an initial guess for carb year 
lagMCMC <- function(carb, carb2, init, lag, iterations){
    reject <- 0
    #prior <- 50000^-1
    #initial vector
    states <- rep(0,iterations)
    # loop through iterations
    for (i in 1:iterations) {
        #create density for initial guess
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

#create posterior from MCMC output
postMCMC <- function(vector){
    out <- seq(1:50000)
    for(i in 1:50000){
        out[i] <- length(vector[which(vector==i)])
    }
    out/sum(out)
}


## Example with 1st carbon date 10,060
# 2nd carbon date 10,090
# intitial value of 11,600 for first carbon date
# 20 year lag between, from 1st carbon
# with 50,000 iterations
lagTest <- lagMCMC(10060,10090,11600,20,50000)
# turn mcmc output into distribution
lag <- postMCMC(lagTest$states)
# plot
plot(lagCarb(10060,10090,20)[calC$Year>11600 & calC$Year<11700]
     ~ calC$Year[calC$Year>11600 & calC$Year<11700],type="l")
lines(lag,col="red")
plot(lag~calC$Year,type='l')

plot(lagTest$states,type="l")


##### Meta-Comment : This really feels like old code, there's cleaner code later 
#####                that does the same thing in half the lines
#####                it's cool.
#####                There are some tests of this function anyway
# MCMC function for two carbon dates with unknown lag time
# Assuming carb < carb2
unknownLag <- function(carb,carb2,init1,init2,iterations){
    states <- rep(0,iterations)
    states2 <- rep(0,iterations)
    reject <- 0
    reject2 <- 0
    for(i in 1:iterations){
        carbpost <- dnorm(carb,calC$C14.Age[init1],calC$Sigma[init1])
        newYear1 <- rnorm(1,init1,100)
        while(newYear1 > init2 | newYear1 < 0){
            newYear1 <- rnorm(1, init1, 100)
        }
        newYear1 <- round(newYear1)
        newPost1 <- dnorm(carb,calC$C14.Age[newYear1],calC$Sigma[newYear1])
        if (newPost1/carbpost >= 1) {
            init1 <- newYear1
        } else {
            #accept/ reject
            draw <- runif(1,0,1)
            if (newPost1/carbpost >= draw){
                init1 <- newYear1
            } else reject = reject + 1
        }
        states[i] <- init1
        ### Start the carb2 part
        carb2post <- dnorm(carb2,calC$C14.Age[init2],calC$Sigma[init2])
        newYear2 <- rnorm(1,init2,100)
        while(newYear2 > 50000 | newYear2 < init1){
            newYear2 <- rnorm(1, init2, 100)
        }
        newYear2 <- round(newYear2)
        newPost2 <- dnorm(carb2,calC$C14.Age[newYear2],calC$Sigma[newYear2])
        if (newPost2/carb2post >= 1) {
            init2 <- newYear2
        } else {
            #accept/ reject
            draw <- runif(1,0,1)
            if (newPost2/carb2post >= draw){
                init2 <- newYear2
            } else reject2 = reject2 + 1
        }
        states2[i] <- init2
    }
    state <- data.frame(states,states2)
    state
}

# Example of unknownLag
## First carbon date is 16,646
    ## initial date 19,700
## Second carbon is 33,120
    ## initial date 37,600
# 50,000 iterations
hamburger <- unknownLag(16645,33120,19700,37600,50000)
# turn mcmc output from first date into distribution 
ketchup <- postMCMC(hamburger$states)
# turn mcmc output from second date into distribution
mustard <- postMCMC(hamburger$states2)
# turn distributions into density plot
mayo <- kde2d(hamburger$states,hamburger$states2,h=c(1,1),n=25)
# look at the density plot, (cc is a color pallete made earlier)
image(mayo,col=cc)
# contour(mayo) #uncomment for contour plot


# Another Example of unknownLag
## First carbon is 31,400
    # initial date is 35,900
## Second carbon date is 31,350
    # initial is 36,000
# 100,000 iterations
picnic <- unknownLag(31400,31350,35900,36000,100000)
# turn mcmc output from first date into distribution
melon <- postMCMC(picnic$states)
# turn mcmc output from second date into distribution
salad <- postMCMC(picnic$states2)
#### To interpolate or not to interpolate.....????
# turns into density plot
potato <- kde2d(picnic$states,picnic$states2,c(1,1),n=25)
# look at density plot
image(potato, col=cc)
# contour(potato) #uncomment for contour plot


######### Numerical Integration with Unknown Lag Time ########
# carb1 < carb2
# min < max, the range over which we integrate over
# Prior: Given carb1, uniform on carb2
uLag <- function(carb1,carb2,min,max){
    # initialize a matrix to hold the vectors
    posts <- matrix(0, nrow=max-min, ncol=max-min)
    # loop over all carb1 values
    for (i in min:(max-1)){
        #calculation of posterior for current carb1
        inorm <- dnorm(carb1,calC$C14.Age[i],calC$Sigma[i])
        # find appropriate number of zeros (spacing out the matrix)
        zeros <- rep(0,i-min)
        prior <- (i-length((i+1):max))^-1
        # calculate posterior for carb1=i and carb2= carb1+1:max and put into matrix
        posts[i-min+1,] <- c(zeros, 
                inorm*dnorm(carb2,calC$C14.Age[(i+1):max],calC$Sigma[(i+1):max])*prior) 
    }
    # return the posterior matrix
    posts
}

# carb1 < carb2
# Conditional uniform prior on carb1
# Chooses appropriate integration bounds for you
uuLag <- function(carb1,carb2){
    #finding bounds by finding where carb1 has positive probabilty
    posterior1 <- post(carb1)
    posterior2 <- post(carb2)
    min <- min(calC$Year[which(posterior1>1E-5)])
    max <- max(calC$Year[which(posterior2>1E-5)])
    # initialize a matrix to hold the vectors
    posts <- matrix(0, nrow=max-min, ncol=max-min)
    # loop over all carb1 values
    for (i in min:(max-1)){
        #calculation of posterior for current carb1
        inorm <- dnorm(carb1,calC$C14.Age[i],calC$Sigma[i])
        # find appropriate number of zeros (spacing out the matrix)
        zeros <- rep(0,i-min)
        prior <- (i-length((i+1):max))^-1
        # calculate posterior for carb1=i and carb2= carb1+1:max and put into matrix
        posts[i-min+1,] <- c(zeros, 
                             inorm*dnorm(carb2,calC$C14.Age[(i+1):max],calC$Sigma[(i+1):max])*prior) 
    }
    # return the posterior matrix
    posts
}


# carb1 < carb2
# chooses appropriate integration bounds
# Like uuLag, but uniform on carb2 instead of 1
# 1E-5 is arbitrary, but seemed to work
qLag <- function(carb1,carb2){
    #finding bounds by finding where carb1 has positive probabilty
    posterior1 <- post(carb1)
    posterior2 <- post(carb2)
    min <- min(calC$Year[which(posterior1>1E-5)])
    max <- max(calC$Year[which(posterior2>1E-5)])
    # initialize a matrix to hold the vectors
    posts <- matrix(0, nrow=max-min, ncol=max-min)
    # loop over all carb1 values
    for (i in (min+1):max){
        #calculation of posterior for current carb1
        inorm <- dnorm(carb2,calC$C14.Age[i],calC$Sigma[i])
        # find appropriate number of zeros (spacing out the matrix)
        zeros <- rep(0,max-i)
        prior <- (i-length(i:max))^-1
        # calculate posterior for carb1=i and carb2= carb1+1:max and put into matrix
        posts[,i-min] <- c(inorm*dnorm(carb1,calC$C14.Age[min:(i-1)],calC$Sigma[min:(i-1)])*prior,
                           zeros) 
    }
    # return the posterior matrix
    posts
}



# Same as uLag, but uniform prior 
# Prior: every point in the triangle given same priority
sLag <- function(carb1,carb2,min,max){
    # initialize a matrix to hold the vectors
    posts <- matrix(0, nrow=max-min, ncol=max-min)
    # loop over all carb1 values
    prior <- (.5*(max-min)^2)^-1
    for (i in min:(max-1)){
        #calculation of posterior for current carb1
        inorm <- dnorm(carb1,calC$C14.Age[i],calC$Sigma[i])
        # find appropriate number of zeros (spacing out the matrix)
        zeros <- rep(0,i-min)
        # calculate posterior for carb1=i and carb2= carb1+1:max and put into matrix
        posts[i-min+1,] <- c(zeros, 
                             inorm*dnorm(carb2,calC$C14.Age[(i+1):max],calC$Sigma[(i+1):max])*prior) 
    }
    # return the posterior matrix
    posts
}

# chooses integration bounds
# carb1 < carb2
# boundary determined by function rather than given as input
ssLag <- function(carb1,carb2){
    posterior1 <- post(carb1)
    posterior2 <- post(carb2)
    min <- min(calC$Year[which(posterior1>1E-5)])
    max <- max(calC$Year[which(posterior2>1E-5)])
    # initialize a matrix to hold the vectors
    posts <- matrix(0, nrow=max-min, ncol=max-min)
    # loop over all carb1 values
    prior <- (.5*(max-min)^2)^-1
    for (i in min:(max-1)){
        #calculation of posterior for current carb1
        inorm <- dnorm(carb1,calC$C14.Age[i],calC$Sigma[i])
        # find appropriate number of zeros (spacing out the matrix)
        zeros <- rep(0,i-min+1)
        # calculate posterior for carb1=i and carb2= carb1+1:max and put into matrix
        posts[i-min+1,] <- c(zeros, 
                             inorm*dnorm(carb2,calC$C14.Age[(i+1):max],calC$Sigma[(i+1):max])*prior) 
    }
    # return the posterior matrix
    posts
}

basket <- uLag(31400,31350,35600,35700)
image(basket, col = grey(seq(0, 1, length = 256)))
bask <- uLag(31400,31350,35000,36000)
image(bask, col = grey(seq(0, 1, length = 256)), useRaster=TRUE)
bas <- uLag(15000,15050,18000,18600)
image(bas,col=cc,axes=FALSE,useRaster=TRUE,xlab="carb1",ylab="carb2",main="uniform given one carbon date")
axis(1,at=c(0,1), labels=c(18000,18600))
axis(2,at=c(0,1), labels=c(18000,18600))
contour(bas, axes=FALSE)
bass <- sLag(15000,15050,18000,18600)
image(bass,col=cc,axes=FALSE,useRaster=TRUE,xlab="carb1",ylab="carb2",main="Uniform over whole space")
axis(1,at=c(0,1), labels=c(18000,18600))
axis(2,at=c(0,1), labels=c(18000,18600))
contour(bass, axes=FALSE)


# These don't take integration bounds
# Works better
# numbers are carbon dates
vas <- uuLag(15000,15050)
image(vas,col=cc,axes=FALSE,useRaster=TRUE,xlab="carb1",ylab="carb2",main="uniform given one carbon1")
contour(vas)
vass <- ssLag(15000,15050)
image(vass,col=cc,axes=FALSE,useRaster=TRUE,xlab="carb1",ylab="carb2",main="Uniform over whole space")
contour(vass)
tass <- qLag(15000,15050)
contour(tass)

# Contour plot of one prior vs the other
difference <- vas-vass
contour(difference)


# Conditioning on carb1 vs carb2:
diff <- vas-tass
contour(diff)
contour(tass-vass)


# So highest probability point is the same in both
wmvas <- which.max(vas)
c(row(vas)[wmvas], col(vas)[wmvas])
wmvass <- which.max(vass)
c(row(vass)[wmvass], col(vass)[wmvass])

# and a matrix of differences!
differences <- vas-vass
contour(differences)
# Looks like the contours for the individual posteriors, so don't think uniform matters that much...


#comparison to MCMC 
mcbas <- unknownLag(15000,15050,18100,18150,25000)
mcmcbas <- kde2d(mcbas$states,mcbas$states2)
image(mcmcbas,col=cc, main="MCMC")



ind1<- post(15000)
ind2 <- post(16000)
plot(post(15000)[calC$Year<19000&calC$Year>18000] ~ 
         calC$Year[calC$Year<19000&calC$Year>18000],type="l")
lines(post(15050),col="red")
plot(post(16000)[calC$Year<20000&calC$Year>17000] ~ 
         calC$Year[calC$Year<20000&calC$Year>17000],type="l")


########## Groups of Trees With Unknown Lag Between #########

# Given carb datem chooses
# carb is a carbon date
# init is the calendar year we're starting with
# lowBound and upBound are the years we have to stay between
# prior is the prior distribution
genUnknownLag <- function(carb,init,lowBound=0,upBound=50000,prior){
    carbpost <- dnorm(carb,calC$C14.Age[init],calC$Sigma[init])*prior[init]
    # draw new year
    newYear1 <- rnorm(1,init,100)
    # make sure it's between appropriate bounds
    while (newYear1 > upBound | newYear1 < lowBound){
        newYear1 <- rnorm(1, init, 100)
    }
    newYear1 <- round(newYear1)
    newPost1 <- dnorm(carb,calC$C14.Age[newYear1],calC$Sigma[newYear1])*prior[newYear1]
    # see if newYear1 has a higher probability than inti
    if (newPost1/carbpost >= 1) {
        init <- newYear1
    } else
        #accept/ reject
        draw <- runif(1,0,1)
        if (newPost1/carbpost >= draw){
            init <- newYear1
        }
    init
}


## >2 carbon dates, carb ordered youngest to oldest
# carb is a vector of carbon dates
## lag is always from last element of carb
# Returns posterior dist'n of last (largest) carbon date. all others inferred
# because we care about when the tree died.
# numerical integration with known lag
lagCar <- function(carb,lag){
    #prior <- rep(50000^-1,50000)
    len <- length(lag)
    # initialize posterior
    # calculate zeros to vector to make length = 50000
    lagAdd <- rep(0,lag[1:(len)]-1)
    # make C14 and Sigma vectors with appropriate vector length and correct subsets of cal Curve
    cAge <- c(calC$C14.Age[lag[1:(len)]:length(calC$C14.Age)], lagAdd)
    sig <- c(calC$Sigma[lag[1:(len)]:length(calC$Sigma)], lagAdd)
    # use subsetted C14 and Sigma values to calculate posterior for particular carbon age.
    # multiply it by current posterior value
    posterior <- posterior*dnorm(carb[1:(len-1)], cAge,sig)
    posterior/sum(posterior)
}

######## Procedure for 2 Trees ############
twoTrees <- function(tree1,lag1,tree2,lag2,iterations,init){
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
                            dnorm(tree2[which.max(tree2)],calC$C14.Age[(i+1):max],
                                  calC$Sigma[(i+1):max])*prior) 
    }
    # return the posterior matrix
    postMaT <- list(posts,min,max)
}

b1<-c(1201,1200,1189,1193)
l1<-c(23,15,3,1)
b2<-c(1259,1263,1272,1274)
l2<-c(27,12,4,1)
lobster <- lagCar(b2,l2)
plot(lobster[calC$Year<1500 & calC$Year>1000]~calC$Year[calC$Year<1500 & calC$Year>1000],
     type="l")
tray <- twoTrees(b1,l1,b2,l2,30000,c(1130,1265))
contour(tray[[1]], axes=FALSE)
axis(1,at=c(0,1), labels=c(tray[[2]],tray[[3]]))
axis(2,at=c(0,1), labels=c(tray[[2]],tray[[3]]))


############# The Actual Procedure ###########
# n trees
# carb is a vector of ALL carbon dates in the sample from least to greatest within the ring
# tree is a vector giving which tree each date is associated with
# lag is a vector giving the lag from the oldest ring. All lag values are positive, decreasing
# Should basically be a data frame so ith entry of each refers to same piece of data
# init is a vector of initial values for last date in each tree
multCarb <- function(carb,tree,lag,iterations, init){
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

carbonDates <- c(2460,2530,2490)
treeNum <-c(1,2,3)
treeLag <- c(0,0,0)
initial <- c(2675,2685,2695)
treeTest <- multCarb(carbonDates,treeNum,treeLag,10000,initial)

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




###### n-Tree procedure using numerical integration: #######
# n trees
# carb is a vector of ALL carbon dates in the sample from least to greatest within the ring
# tree is a vector giving which tree each date is associated with
# lag is a vector giving the lag from the oldest ring. All lag values are positive, decreasing
# Should basically be a data frame so ith entry of each refers to same piece of data
# init is a vector of initial values for last date in each tree
multCarbNumInt <- function(carb,tree,lag){
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
    
    # Now to integrate between the trees:
    # initialize an array to hold the data:
    posterior1 <- post(carb[endTrees[1]])
    posterior2 <- post(carb[endTrees[numTrees]])
    min <- min(calC$Year[which(priors[,1]>1E-5)])
    max <- max(calC$Year[which(priors[,numTrees]>1E-5)])
    leng <- rep(max-min,3)
    post.array <- array(data=0,dim=length(leng))
    tree3 <- endTrees[3]
    for (j in min:(max-1)){  # loop over last tree
        zer <- rep(0,j-min+1)
        jay <- c(zer,dnorm(carb[tree3],calC$C14.Age[(j+1):max],calC$Sigma[(j+1):max])*priors[(min:max),3])
        for (i in min:(max-1)){   # loop over first tree, vectorize 2nd tree
            #calculation of posterior for current carb1
            inorm <- dnorm(carb[endTrees[1]],calC$C14.Age[i],calC$Sigma[i])*priors[carb[endTrees[1]],2]
            # find appropriate number of zeros (spacing out the array)
            zeros <- rep(0,i-min+1)
            # calculate posterior for carb1=i and carb2= carb1+1:max and put into matrix
            post.array[i-min+1, ,j-min+1] <- c(zeros, 
                                 inorm*dnorm(carb[endTrees[2]],calC$C14.Age[(i+1):max],calC$Sigma[(i+1):max])
                                 *priors[(min:max),1])*jay 
        }
    }
}


carbonDates <- c(b1,b2,1270,1250,1258,1265,1270)
treeNames <- c(rep(1,4),rep(2,4),rep(3,5))
treeLags <- c(l1,l2,49,23,9,6,1)
numTest <- multCarbNumInt(carbonDates,treeNames,treeLags)


