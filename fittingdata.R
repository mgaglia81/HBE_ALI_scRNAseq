## nls function

library(minpack.lm)

# this code was designed for use with HEK-Blue IFN-lambda reporter cells (Invivogen)
# defining x and y coordinates : for x coordinates use the log transformed standard concentrations of IL-28

# example standard: x <- c(6, 5, 4, 3, 2, 1, 0, -1, -2) is 1000000, 100000, 10000, 1000, 100, 10, 1, 0.1, 0.001 fg/ul
x <- c(6, 5, 4, 3, 2, 1, 0, -1, -2)

#y should change relative to each plate - enter average absorbance values for each standard 

y <- c(abs1,	abs2,	abs3,	....)

# output to be present as PNG file 
#png(file ="nls.png") 

# Taking the model to get fitted 
# I() clarifies that + and - are mathematical

m <- nlsLM(y~ a+((b-a)/(1+2.718282^(-d*(x-c)))), start = list(a = 0.4, b=2, c=1.5, d=2)) 

#a is min, b is max, c is inflection, d is gradient

# plot the graph 
plot(x, y, col.lab ="red",  
     col.axis ="red",
     xlab = "Standard (log transformed picograms)",
     ylab = "Absorbance",
     main = "Non-linear least squares fit for IFN standard") 

# plot the graph with new fitting line 
# or regression line 
lines(x, predict(m)) 

# saving the file 
#dev.off() 

# print minimum residual or error value 
print(sum(resid(m)^2)) 
SSR <- sum(resid(m)^2)
GOF <- cor(y, predict(m)) 

plot(x, y, col.lab ="red",  
     col.axis ="red",
     xlab = "Standard (log transformed picograms)",
     ylab = "Absorbance",
     main = "Non-linear least squares fit",
     col.sub = "blue",
     sub = bquote("SSR" == .(SSR) ~ "GOF" == .(GOF)),
     lines(x, predict(m)))

summary(m)

#will vary for each plate 
#x=((log(( (b-a) /(y-a) )-1))/(-d))+c
#y=a+(b-a)/(1+exp(-lambda*(x-c)))
#plot(a+((b-a)/(1+2.718282^(-d*(x-c)))))

q = ((log(( (b-a) /(v-a) )-1))/(-d))+c

#enter values from summary(m) to make the nls equation

a = 0.105703 #example
b = 1.350947 #example
c = 3.145969 #example
d = 6.262297 #example

q = ((log(( (b-a) /(v-a) )-1))/(-d))+c

#enter the absorbance values for each sample
v = c(absA, absB, absC, ...)

#q will output the concentration for each sample
q
