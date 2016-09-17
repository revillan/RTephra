###### Test of lagCar, >2 carbon dates w/ known lag ####
knownLagTest <- lagCar( c(15000,15060), 40 )
plot(knownLagTest[15000:20000]~calC$Year[15000:20000], type= "l")
