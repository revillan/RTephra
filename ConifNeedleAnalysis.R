library(gridExtra)


mcmcSing <- function(carb, sigma ,iterations, data = Reimer13){
    reject <- 0
    init <- sample(singInit(carb), 1 )
    #initial vector
    arr <- rep(0,iterations)
    # loop through iterations
    for (i in 1:iterations) {
        #create density for initial guess
        postOrg <- dnorm(carb,data$C14.Age[init],sigma) # *PRIOR
        # pick a new year from a normal dist around the inital guess
        newYear <- rnorm(1, init, 40)
        # make sure the new year is positive integer & in the support
        while (newYear > 500000 | newYear < 0){
            newYear <- rnorm(1, init, 40)
        }
        newYear <- round(newYear)
        # create density for proposed new year
        postNew <- dnorm(carb, data$C14.Age[newYear],sigma)  # *PRIOR
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
    print(acceptanceRate)
    arr
}

needleNum <- post(8905, 20)
set.seed(64)
needleMCMC <- mcmcSing(8905,20,2000000)

keep <- seq(10000,2000000, by = 5)

needleMDens <- postMCMC(needleMCMC[keep])
needle <- data.frame(needleMDens, needleNum,"index"=c(1:50000))



pdf("num&mcmc.pdf",width=6,height=4)
ggplot(data = needle) + 
    geom_line(aes( x = c(1:50000), y = needleNum)) +
    geom_line(aes(x = c(1:50000), y = needleMDens ),col = "peru", alpha = 0.85) +
    xlim(c(9850,10200)) + #theme_bw() + 
    xlab("Age (Cal BP)") + ylab("Density")
dev.off()

p <- hdpRegion(needleNum)
numNeed <- ggplot(data=needle) + 
    geom_ribbon(aes(ymax = ifelse(needleNum > p,needleNum,0), x = index),
                    ymin = 0, fill="peru", alpha= 0.6) + 
    geom_line(aes( x = index, y = needleNum)) +
    xlim(c(9900,10200)) + theme_bw() +
    xlab("Age (Cal BP)") + ylab("Density")
    
q <- hdpRegion(needleMDens)
mcmcNeed <- ggplot(data=needle) + 
    geom_ribbon(aes(ymax = ifelse(needleMDens > q,needleMDens,0),
                     x = index), ymin = 0, fill="peru", alpha = 0.6) + 
    geom_line(aes( x = index, y = needleMDens)) + 
    xlim(c(9900,10200)) + theme_bw() +
    xlab("Age (Cal BP)") + ylab("Density")

pdf("sidebyside.pdf",height=2,width=6)
grid.arrange(numNeed,mcmcNeed, ncol=2)
dev.off()

which(needleNum > p) # (9923,10069) U (10117,10171)
which(needleMDens > q) # (9923,10069) U (10117,10170)

pdf("needlerunmean.pdf",width=6,height=2)
qplot(x=c(1:50000),y=cumsum(needleMCMC[1:50000])/seq(along=needleMCMC[1:50000]),
      geom="line") + geom_line(aes(y=mean(needleMCMC)), col="tan3") + 
    theme_bw() + xlab("Iterations") + ylab("Running Mean") + ylim(c(10025,10060))
dev.off()

needleNum[max(needleNum)]

########## Credible Interval Plot #########

cred <- post(9000,60)
pdf("credintexam.pdf", width = 6, height = 3)
qplot(y = cred, x = c(1:50000),geom="line") + xlim(c(9750,10350)) +
    geom_ribbon(aes(ymax = ifelse(cred > hdpRegion(cred),cred,0), 
    x = c(1:50000)), ymin = 0, fill="peru", alpha=0.6) +
    xlab("Age (Years BP)") + ylab("Density") +# theme_bw() + 
    geom_line(y = hdpRegion(cred), col="peru")
dev.off()
