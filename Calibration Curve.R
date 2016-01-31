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


#### Single carbon date into distribution ####
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



###### ggplot2 ######

######  Raw Cal-Curve Data #####3
# error bar aesthetics 
errorLimits <- aes(ymin = X14C.age - Sigma,
                   ymax = X14C.age + Sigma,
                   )

# whole raw plot SLOWWWWWW    UGLYYYY
rei <- ggplot( reimer09Data, aes( x = CAL.BP, y = X14C.age) )
rei + geom_point(color="salmon",size=0.3) +
      geom_errorbar(errorLimits, size=0.3, color="blue") +
      scale_x_log10() +
      scale_y_log10() 

# Zommed-in slice of cal curve
plot_slice <- reimer09Data %>% filter(22000 < CAL.BP & CAL.BP < 24000)

# make ggplot for slice
slice <- ggplot(plot_slice, aes( x= CAL.BP, y = X14C.age) )


# add layers to look pretty
slice + geom_point(color="salmon",size=2) +
  #  coord_trans(xtrans = "log10", ytrans = "log10") +
    geom_errorbar(errorLimits, size=0.3, color="blue")

####### Interpolated Cal Curve #####

errorLimits <- aes(ymin = C14.Age - Sigma,
                   ymax = C14.Age + Sigma,
)

# whole interpolate 
mer <- ggplot(calC, aes(x = Year, y = C14.Age))
mer + geom_line()

# filter whole into part
mer_slice <- calC %>% filter( 23000 < Year & Year < 24000 )

# ggplot of partial
slice2 <- ggplot( mer_slice, aes(x = Year, y = C14.Age))
slice2 + geom_line() + 
         geom_errorbar(errorLimits, color = "turquoise", alpha=0.5)
        


