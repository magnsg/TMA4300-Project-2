load("rain.rda")#Set working directory correctly
plot(rain)
sapply(rain,class)

#a
plot(rain$day,rain$n.rain)

print(rain$n.years[61])

#b

prob_y_t <- function(y_t,t){
  dbinom(y_t,size = (rain$n.years[t]), prob = rain$n.rain[t]/rain$n.years[t])
}

plot(prob_y_t(c(0:30),1))
plot(prob_y_t(20,rain$day))

