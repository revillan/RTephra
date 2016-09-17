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

#just another example year: 20828
plot(posterior(20828)[calCurve$Year > 20000] ~ calCurve$Year[calCurve$Year>20000],type="l")

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