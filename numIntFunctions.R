############
# All of the functions take uncalibrated carbon dates as input and through 
# numerical integration calculate the posterior density for the calendar year the 
# samples originate from.
# 
# They all rely a calibration curve (Reimer13) to find likely calendar years.
#
############

## Read in Calibration data
IntCal13 <- read.table("http://www.radiocarbon.org/IntCal13%20files/intcal13.14c",
                       skip=8,sep=",")[,c(1,2,5)]

# Linear Interpolation of 
CalC <- approx(IntCal13$V1, IntCal13$V2, n=50000, xout=c(1:50000))
CalS <- approx(IntCal13$V1, IntCal13$V5, n=50000, xout=c(1:50000))

Reimer13 <- data.frame("Year"=CalC$x, "C14 Age"=CalC$y, "Sigma"=CalS$y)


# posterior density function for Mean year
# Numerical integration in 1 dimension
post <- function (carb, data = Reimer13){
    year <- 1:50000
    prior <- rep(50000^-1, 50000)
    poster <- dnorm(carb,data$C14.Age, data$Sigma) * prior
    poster/sum(poster)
}

# chooses integration bounds
# carb1 < carb2

numInt2D <- function(carb1,carb2, data = Reimer13, prior = (.5*(max-min)^2)^-1, 
                  cutoff = 1E-5){
    # find the posterior densities of the carbon dates by themselves
    posterior1 <- post(carb1)
    posterior2 <- post(carb2)
    # find the bounds of integration by seeing where the posteriors for the 
    # boundary samples are not 0
    min <- min(data$Year[which(posterior1 > cutoff)])
    max <- max(data$Year[which(posterior2 > cutoff)])
    # initialize a matrix to hold the vectors
    posts <- matrix(0, nrow=max-min, ncol=max-min)
    # loop over all carb1 values
    for (i in (min+1):max){
        #calculation of posterior for current carb1
        inorm <- dnorm(carb1,data$C14.Age[i],data$Sigma[i])
        # find appropriate number of zeros (spacing out the matrix)
        zeros <- rep(0,max-i)
        # calculate posterior for carb1=i and carb2= carb1+1:max 
        # and put into matrix
        posts[,i-min] <- c(zeros, 
                           inorm*dnorm(carb2,data$C14.Age[min:(i-1)],
                                       data$Sigma[min:(i-1)]) * prior) 
    }
    # return the posterior matrix
    posts
}

test2 <- numInt2D(15000,15050)
which(test2 == max(test2), arr.ind=TRUE)
dim(test2)
contour(test2[37500:38777,4300:5350])




### 3D Array
# 3 carb values where carb1 < carb2 < carb3
numInt3D <- function(carb1,carb2,carb3, data=Reimer13, cutoff=1E-5){
    # calculate poster for oldest and youngest samples 
    posterior1 <- post(carb1)
    posterior3 <- post(carb3)
    # find the bounds of integration by seeing where the posteriors for the 
    # boundary samples are not 0
    min <- min(data$Year[which(posterior1 > cutoff)])
    max <- max(data$Year[which(posterior3 > cutoff)])
    len <- max - min
    # initialize empty array with each dimension of length len
    posts <- array(data = 0, rep(len, 3))
    # looping over carb1 values
    for (j in (min+1):max-2){
        jnorm <- dnorm(carb1,data$C14.Age[j], data$Sigma[j])
          # looping over carb2 values with carb2 > carb1
          for (i in (j+1):max-1){
            inorm <- dnorm(carb2, data$C14.Age[i], data$Sigma[j])
            # looping over carb3 > carb2
            for (k in (i+1):max){
                knorm <- dnorm(carb3, data$C14.Age[k], data$Sigma[k])
                # inserting values into array
                posts[j-min, i-min, k-min] <- jnorm * inorm * knorm
            }
        }
    }
    posts
}


test3 <- numInt3D(17000,17025,17050)
which(test3 == max(test3), arr.ind=TRUE)
contour(test3[38, ,])

# Test where cutoff value is smaller
test3.1 <- numInt3D(15000,15025,15050, cutoff = 1E-6)
dim(test3.1)

# Test with greater distance between samples
## This call crashes R for me. 
test3.2 <- numInt3D(15000,17000,17050)
