#### Read in Calibration Data, Interpolate, Data.Frame ####
library(ggplot2)
library(dplyr)
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
errorLimits <- aes(ymin = C14.Age - Sigma,
                   ymax = C14.Age + Sigma,
                   )

# whole raw plot SLOWWWWWW    UGLYYYY
rei <- ggplot( reimer09Data, aes( x = CAL.BP, y = X14C.age) )
rei + geom_point(color="salmon",size=0.3) +
      geom_errorbar(errorLimits, size=0.3, color="blue") +
      scale_x_log10() +
      scale_y_log10() 

# Zommed-in slice of cal curve
plot_slice <- Reimer13 %>% filter(10000 < Year & Year < 10500)

# make ggplot for slice
slice <- ggplot(plot_slice, aes( x= Year, y = C14.Age) )


# add layers to look pretty
slice + geom_point(size=2) +
  #  coord_trans(xtrans = "log10", ytrans = "log10") +
    geom_errorbar(errorLimits, size=0.3, color="peru")

####### Interpolated Cal Curve #####

errorLimits <- aes(ymin = C14.Age - Sigma,
                   ymax = C14.Age + Sigma,
)

# whole interpolate 
mer <- ggplot(Reimer13, aes(x = Year, y = C14.Age))
mer + geom_line()

# filter whole into part
mer_slice <- Reimer13 %>% filter( 9075 < Year & Year < 9250 )

# ggplot of partial
pdf("calseg.pdf", width=6, height = 3.5)
slice2 <- ggplot( mer_slice, aes(x = Year, y = C14.Age))
slice2 + geom_line() + 
         geom_ribbon(errorLimits, alpha=0.3) +
        theme_bw() + xlab("Year (Cal BP)") + ylab("Conventional Carbon Date")
dev.off()

## pic of cal curve in Intro
pdf("fullcal.pdf", width=6,height = 4)
ggplot(Reimer13, aes(x = Year, y = C14.Age)) + geom_line(size = 1.5)  +# theme_bw() + 
    geom_abline() + ylim(0,50000) + labs(x="Calendar Date (Cal BP)", y = "Conventional Carbon Date")
dev.off()



Combo <- left_join(IntCal13,Reimer13, by = "Year")

Combo_slice <- Combo %>% filter(Year < 8000 & Year > 7800)

pdf("interpol.pdf",width=6,height = 3.2)
ggplot(data = Combo_slice) + geom_line(aes(x=Year, y=C14.Age)) + 
    geom_point(aes(x=Year, y = V2)) +# theme_bw() + 
    labs(x= "Year (Cal BP)", y = "Conventional Carbon Date")
dev.off()


pdf("partcal.pdf", width = 6, height = 3)
ggplot(Combo_slice, aes(x = Year, y = C14.Age)) + geom_line(size = 1.5) +
    geom_hline(aes(yintercept = 7040), size = 0.5) +
    labs(x="Calendar Date (Cal BP)", y = "Conventional Carbon Date")# + theme_bw() 
dev.off()
