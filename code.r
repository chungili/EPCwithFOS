# Read the raw data
dt = read.csv("uci-secom.csv", header = F)
# There are 592 variables in the dataset.
# In the CSV file, the first column represents time stamps, 
# and the subsequent columns represent the quality characteristics.

# function for data pre-processing
prepro = function(x){
  #filters out the missing values
  xfull = x[!is.na(x)]
  end = floor(length(xfull)*0.9)
  xtrim = xfull[1:end]
  # filters out and select only unique values
  xfinal = unique(xtrim)
  return(xfinal)
}

# The 3nd, 12th, 216th and 221st variables are selected and 
# then named as x1, x11, x215, and x220, respectively.

x1 = prepro(dt$V3)
x11 = prepro(dt$V12)
x215 = prepro(dt$V216)
x220 = prepro(dt$V221)

### FOS
FOS = function(x, pn = 0.9, alpha = 0.0027){
  Rx = sort(x)
  n = length(x)
  n1 = n + 1
  Xu = function(u){
    if (0<u && u<=1/n1) {
      Rx[1]+(Rx[2]-Rx[1])*log((n1)*u)
    } else if (1/n1<u && u<n/n1) {
      (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
    } else {
      Rx[n]-(Rx[n]-Rx[n-1])*log(n1*(1-u))
    }
  }
  findkl = function(r, n, alpha = alpha, pn = pn ) pbeta(alpha/2, r, n-r+1)-(1+pn)/2
  findku = function(r, n, alpha = alpha, pn = pn ) pbeta(1-alpha/2, r, n-r+1)-(1-pn)/2
  kl <- uniroot(findkl, c(0, n1), n = n, alpha = alpha, pn = pn)$root
  ku <- uniroot(findku, c(0, n1), n = n, alpha = alpha, pn = pn)$root
  ur = kl/n1
  us = ku/n1
  LCL = Xu(ur)
  UCL = Xu(us)
  return(list(UCL=UCL,LCL=LCL))
}

### BH
BH = function(x, pn = 0.9, alpha = 0.0027, B = 5000){
  n = length(x)
  q.00135 = q.99865 = q.5 = c()
  for (j in 1:B) {
    boot.samp = sample(x, n, replace = TRUE)
    Rx = sort(boot.samp)
    n1 = n+1
    limit11 = function(u, n){
      if (0<u && u<=1/n1) {
        Rx[1]+(Rx[2]-Rx[1])*log((n+1)*u)
      } else if (1/n1<u && u<n/n1) {
        (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
      } else {
        Rx[n]-(Rx[n]-Rx[n-1])*log(n1*(1-u))
      }
    }
    
    q.00135[j] = limit11(u = alpha/2, n)
    q.99865[j] = limit11(u = 1-alpha/2, n)
    q.5[j] = limit11(u = 0.5, n) 
  }

  #Two-Sided limits, Hutson (2002) 
  
  Rx = sort(q.00135)
  n = B
  n1 = n+1
  
  limit1 = function(u){
    if (0<u && u<=1/n1) {
      Rx[1]+(Rx[2]-Rx[1])*log((n1)*u)
    } else if (1/n1<u && u<n/n1) {
      (1-(n1*u-floor(n1*u)))*Rx[floor(n1*u)]+(n1*u-floor(n1*u))*Rx[floor(n1*u)+1]
    } else {
      Rx[n]-(Rx[n]-Rx[n-1])*log(n1*(1-u))
    }
  }
  
  LCL = limit1(1-(1+pn)/2)
  
  Rx = sort(q.99865)

  UCL = limit1((1+pn)/2)
  
  return(list(UCL = UCL, LCL = LCL))
}

### Goedhart et al. (2020)
GSD = function(x, pn = 0.9, alpha = 0.0027){
  library(tolerance)
  out <- nptol.int(x = x, alpha = (1-pn), 
                   P = (1-alpha), side = 2,
                   method = "YM", upper = NULL, lower = NULL)
  out1 = data.frame(out)
  LCL = out1$X2.sided.lower
  UCL = out1$X2.sided.upper
  return(list(UCL = UCL, LCL = LCL))
}
