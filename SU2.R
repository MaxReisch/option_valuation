#http://linanqiu.github.io/2015/10/26/Binomial-European-Option-Pricing-in-R/


#Calculate up&Down

delta = 3/12

u = exp(sigma*sqrt(delta))
d = exp(-sigma*sqrt(delta))



# Build binomial tree

stock_tree = function(S, sigma, delta, N) {
  tree = matrix(0, nrow=N+1, ncol=N+1)
  u = exp(sigma*sqrt(delta))
  d = exp(-sigma*sqrt(delta))
  for (i in 1:(N+1)) {
    for (j in 1:i) {
      tree[i,j] = S * u^(j-1) * d^((i-1)-(j-1))
    }
  }
  return(tree)
}


tree =stock_tree(S,sigma,delta,N = 1)

# Calculate the riskneutral probability

q_prob = function(r, delta, sigma) {
  u = exp(sigma*sqrt(delta))
  d = exp(-sigma*sqrt(delta))
  return((exp(r*delta) - d)/(u-d))
}



#Backcalculation

value_binomial_option = function(tree, sigma, delta, r, barrier) {
  q = q_prob(r, delta, sigma)
  option_tree = matrix(0, nrow=nrow(tree), ncol=ncol(tree))
  for (i in nrow(tree):1) {
    for(j in 1:i) {
      if(tree[i,j] >= barrier) {
        option_tree[i, j] = 0
      } else {
        option_tree[i, j] = (barrier -tree[i,j])*0.01
      }
      
    }
    
  }
    
  for (i in (nrow(tree)-1):1) {
    for(j in 1:i) {
      option_tree[i, j] = ((1-q)*option_tree[i+1,j] + q*option_tree[i+1,j+1])/exp(r*delta)
    }
  }
  return(option_tree)
}

# Write all together 
binomial_certificate = function(S, sigma, r, t, T, barrier,N){
  q = q_prob(r=r, delta=delta, sigma=sigma)
  tree = stock_tree(S=S, sigma=sigma, delta=delta, N=N)
  certificate = value_binomial_option(tree = tree, sigma= sigma, delta=delta, r=r, barrier=barrier)
  return(list(q=q, stock = tree, certificate = certificate, price= certificate[1,1]))
}

# Market data

S = 12661
sigma = 0.1539
t = 95
T = 365
r = - 0.00035  #EURIBOR 3M
barrier = 13000


binomial_certificate(S = S, sigma = sigma, r = r, t = t, T = T, barrier = barrier, N = 1)
binomial_certificate(S = S, sigma = sigma, r = r, t = t, T = T, barrier = barrier, N = 10)
binomial_certificate(S = S, sigma = sigma, r = r, t = t, T = T, barrier = barrier, N = 100)
binomial_certificate(S = S, sigma = sigma, r = r, t = t, T = T, barrier = barrier, N = 1000)

Tm = T-t


black_scholes <- function(S, k, Tm, r, sigma){
  values <- c(1)
  d1 <- (log(S/k) + (r+(sigma^2)/2)*(Tm))/(sigma*sqrt(Tm))
  d2 <- (log(S/k) + (r-(sigma^2)/2)*(Tm))/(sigma*sqrt(Tm))
  
  values <- (k*exp(-r*Tm)*pnorm(-d2)- S*pnorm(-d1))*0.01
  
  values
}

black_scholes(S = S, k = barrier, Tm = t, r = r, sigma = sigma)
