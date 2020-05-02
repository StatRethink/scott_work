# RETHINKING STATS 

# SECTION 2.3.2.1

ans = dbinom(6, size=9, prob=0.5)
print(ans)


# SECTION 2.4.3 - GRID APPROXIMATION

# define grid
p_grid <- seq( from=0 , to=1 , length.out=20 )

# define prior
prior <- rep( 1 , 20 )
#prior <- ifelse( p_grid < 0.5 , 0 , 1 )
#prior <- exp( -5*abs( p_grid - 0.5 ) )

# compute likelihood at each value in grid
likelihood <- dbinom( 6 , size=9 , prob=p_grid )

# compute product of likelihood and prior
unstd.posterior <- likelihood * prior

# standardize the posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

# plot results

plot( p_grid , posterior , type="b" ,
      xlab="probability of water" , ylab="posterior probability" )
mtext( "20 points" )


# SECTION 2.4.4 - QUADRATIC APPROXIMATION

library(rethinking)
globe.qa <- quap(
  alist(
    W ~ dbinom( W+L, p),    # binomial likelihood
    p ~ dunif(0,1)          # uniform prior
  ),
  data=list(W=6, L=3))

# display summary of quadratic approximation
precis( globe.qa)
  
# analytical calculation
W <- 6
L <- 3
curve( dbeta(x, W+1, L+1 ), from=0, to=1)
# quadratic approximation
curve( dnorm( x, 0.67, 0.16 ), lty=2, add=TRUE )


# SECTION 2.4.5 - MARKOV CHAIN MONTE CARLO

n_samples <- 1000
p <- rep( NA, n_samples)
p[1] <- 0.5
W <- 6
L <- 3
for ( i in 2:n_samples){
  p_new <- rnorm( 1, p[i-1], 0.1 )
  if ( p_new < 0 ) p_new <- abs( p_new )
  if ( p_new > 1 ) p_new <- 2-p_new
  q0 <- dbinom( W, W+L, p[i-1] )
  q1 <- dbinom( W, W+L, p_new )
  p[i] <- ifelse( runif(1) < q1/q0, p_new, p[i-1] )
}

dens( p, xlim=c(0,1) )
curve( dbeta( x, W+1, L+1 ), lty=2, add=TRUE)


# SECTION 2.6 - PRACTICE

# QUESTION 2H1

Pr_A = 0.5
Pr_B = 0.5
Pr_TA = 0.1
Pr_TB = 0.2

Pr_T = Pr_A * Pr_TA + Pr_B * Pr_TB

Pr_A = Pr_TA * Pr_A / Pr_T
Pr_B = Pr_TB * Pr_B / Pr_T

Pr_T2 = Pr_A * Pr_TA + Pr_B * Pr_TB

  # alternatively...

p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep(1,1000)
likelihood <- p_grid * dbinom( 1, size=1 , prob=0.1 ) + ( 1-p_grid) * dbinom( 1, size=1, prob=0.2)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)

plot( p_grid , posterior , type="b" ,
      xlab="probability of success" , ylab="posterior probability" )
mtext( "Twins" )

samples <- sample( p_grid , size=1e4 , replace=TRUE , prob=posterior )
plot( samples ) 


  # alternatively ...

n <- 100000
species <- rbinom( n, size=1, prob=0.5 )
offspringA <- ( 1 + rbinom( n, size=1, prob=0.1 ) )
offspringB <- ( 1 + rbinom( n, size=1, prob=0.2 ) )
offspring <- (1 - species ) * offspringA + species * offspringB
results <- species*10 + offspring
table(results)/n






# CHAPTER 3

Pr_Positive_Vampire <- 0.95  #True positive rate of test
Pr_Positive_Mortal <- 0.01   #False positive rate of test
Pr_Vampire <- 0.001          #Only 1 in 1000 people are vampires
Pr_Positive <- Pr_Positive_Vampire * Pr_Vampire + Pr_Positive_Mortal * (1 - Pr_Vampire)
( Pr_Vampire_Positive <- Pr_Positive_Vampire * Pr_Vampire / Pr_Positive)


# SECTION 3.1 Sampling from a grid-approximate posterior

p_grid <- seq( from=0, to=1, length.out=1000)
prob_p <- rep(1, 1000)
prob_data <- dbinom(6, size=9, prob=p_grid)
posterior <- prob_data * prob_p
posterior <- posterior / sum(posterior)

plot( p_grid , posterior , type="b" ,
      xlab="probability of success" , ylab="posterior probability" )
mtext( "Example" )

samples <- sample( p_grid, prob=posterior, size=10000, replace=TRUE)
plot( samples )     # samples are most common around 0.6

library(rethinking)
dens( samples )

# SECTION 3.2.1

sum( posterior[ p_grid < 0.5 ] )

sum( samples < 0.5 ) / 10000
sum( samples > 0.5 & samples < 0.75) / 10000

# SECTION 3.2.2

quantile( samples, 0.8)
quantile( samples, c(0.1, 0.9) )

p_grid <- seq( from=0, to=1, length.out=1000 )
prior <- rep(1, 1000)
likelihood <- dbinom( 3, size=3, prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
samples <- sample( p_grid, size=10000, replace=TRUE, prob=posterior )

PI( samples, prob=0.5 )

HPDI( samples, prob=0.5 )


# SECTION 3.2.3

p_grid[ which.max(posterior)]
chainmode( samples, adj=0.01 )
(mean(samples))
(median(samples))

sum( posterior * abs( 0.5 - p_grid ) )

loss <- sapply( p_grid, function(d) sum( posterior * abs( d - p_grid ) ) )
p_grid[ which.min(loss) ]
median(samples)

loss <- sapply( p_grid, function(d) sum( posterior * ( d - p_grid )^2 ) )
p_grid[ which.min(loss) ]
mean(samples)



# SECTION 3.3.1

dbinom( 0:2, size=2, prob=0.7 )
rbinom( 1, size=2, prob=0.7 )
rbinom( 10, size=2, prob=0.7 )
dummy_w <- rbinom(100000, size=2, prob=0.7 )
table(dummy_w)/100000

dummy_w <- rbinom( 100000, size=9, prob=0.7 )
simplehist( dummy_w, xlab="dummy water count" )

dummy_w <- rbinom( 100000, size=100, prob=0.7 )
simplehist( dummy_w, xlab="dummy water count" )


# SECTION 3.3.2

w <- rbinom( 10000, size=9, prob=0.6)
simplehist(w)

w <- rbinom( 10000, size=9, prob=samples)
simplehist(w)


# SECTION 3.5 PRACTICE

# QUESTION 3H1

birth1 <- c(1,0,0,0,1,1,0,1,0,1,0,0,1,1,0,1,1,0,0,0,1,0,0,0,1,0,
            0,0,0,1,1,1,0,1,0,1,1,1,0,1,0,1,1,0,1,0,0,1,1,0,1,0,0,0,0,0,0,0,
            1,1,0,1,0,0,1,0,0,0,1,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,0,1,0,1,1,0,
            1,0,1,1,1,0,1,1,1,1)
birth2 <- c(0,1,0,1,0,1,1,1,0,0,1,1,1,1,1,0,0,1,1,1,0,0,1,1,1,0,
            1,1,1,0,1,1,1,0,1,0,0,1,1,1,1,0,0,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,
            1,1,1,0,1,1,0,1,1,0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,1,0,0,1,0,0,1,1,
            0,0,0,1,1,1,0,0,0,0)

library(rethinking)
data(homeworkch3)

sum(birth1) + sum(birth2)

p_grid <- seq( from=0, to=1, length.out=1000)
prob_p <- rep(1, 1000)

n=200
likelihood <- dbinom(sum(birth1) + sum(birth2), size=n, prob=p_grid)
posterior <- likelihood * prob_p
posterior <- posterior / sum(posterior)

plot( p_grid , posterior , type="b" ,
      xlab="probability of success" , ylab="posterior probability" )
mtext( "Example" )
p_grid[ which.max(posterior) ]


samples <- sample( p_grid, prob=posterior, size=10000, replace=TRUE)
plot( samples )     


p_grid[ which.max(posterior)]
chainmode( samples, adj=0.01 )
(mean(samples))
(median(samples))

loss <- sapply( p_grid, function(d) sum( posterior * abs( d - p_grid ) ) )
p_grid[ which.min(loss) ]
median(samples)

loss <- sapply( p_grid, function(d) sum( posterior * ( d - p_grid )^2 ) )
p_grid[ which.min(loss) ]
mean(samples)

plot(loss)
