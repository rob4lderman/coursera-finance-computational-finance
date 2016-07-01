
# Coursera: Intro to Computational Finance: Week 4: Matrix Review and Time Series Concepts

* [Matrix Review: Portfolio Math](#mrpm)
* [MultiVariate Normal Distribution](#mvnd)
* [Time Series Concepts](#tsc)
* [Auto-correlation and Auto-covariance](#acaac)
* [White Noise Process](#wnp)
* [Non-Stationary Processes](#nsp)
* [Moving Average Processes](#map)
* [Auto-Regressive Processes](#arp)


## <a name="mrpm"></a>Matrix Review: Portfolio Math


    R_p = x1*R1 + x2*R2 + x3*R3 ....

    R_p ~ N(mu_p, sigma_p^2)

    E[R_p] = x1*mu1 + x2*mu2 + x3*mu3 ....

    VaR(R_p) = sigma_p^2 = x1^2*sigma_1^2 +
                           x2^2*sigma_2^2 +
                           x3^2*sigma_3^2 + 
                           ... +
                           2 * x1 * x2 * sigma_12 + 
                           2 * x1 * x3 * sigma_13  +  
                           ...


    R = [ R1 ]      # Random vector = vector of random variables
        [ R2 ]
        [ R3 ]

    MU = [ mu1 ]
         [ mu2 ]
         [ mu3 ]

    x = [ x1 ]      # weights
        [ x2 ]
        [ x3 ]

    covariance matrix: n x n matrix, where n= number of R's

    Sigma = [ sigma_1^2     sigma_12        sigma_13   ]
            [ sigma_21      sigma_2^2       sigma_23   ]
            [ sigma_31      sigma_23        sigma_3^2  ]


    E[R_p] = x' * MU

    Var(R_p) = x' * Sigma * x

    R: t(x) %*% Sigma %*% x

    R_p = x' * R

    R: crossprod(x',R)



### Covariance Between Two Portfolios:


    x = [ x1 ]      # weights
        [ x2 ]
        [ x3 ]

    y = [ y1 ]      # weights
        [ y2 ]
        [ y3 ]

    R_p,x = x'*R
    R_p,y = y'*R


    Covariance between portfolios:

    cov(R_p,x, R_p,y) = x' * Sigma * y
                    
                      = y' * Sigma * x

    Note that R_p is a random variable: linear combo of individual asset returns


## <a name="mvnd"></a>MultiVariate Normal Distribution


* The PDF of the BiVariate Normal Distribution...
    * which is a big nasty formula
* can be re-written in matrix algebra


Define:

    Z = [ X ]   x = [ x ]   mu = [ mu_x ]  
        [ Y ]       [ y ]        [ mu_y ]


    covariance matrix: Sigma = [ sigma_x^2  sigma_xy  ] 
                               [ sigma_xy   sigma_y^2 ]


The BiVariate Normal Distribution PDF becomes:


                1
    f(x) = ------------------------ * exp( -1/2 (x - mu)' * Sigma^-1 * (x - mu) )
             2pi * det(Sigma)^1/2


    det(Sigma) = sigma_x^2 * sigma_y^2 - sigma_xy^2
               
                 sigma_xy^2 = sigma_x^2 * sigma_y^2 * rho_xy^2

               = sigma_x^2 * sigma_y^2 * (1 - rho_xy^2)

    Z ~ N(mu, Sigma)


    (x - mu)' * Sigma^-1 * (x - mu)     # quadratic form: vector-matrix' * square-matrix * vector-matrix
                                        # formula for an ellipse

## <a name="tsc"></a>Time Series Concepts

* Time Series Process 
    * Stochastic (random) process
    * sequence of random variables indexed by time
* e.g stock prices
* e.g. momentum trading
* "auto-correlation"

.

    Infinite series:
    { ..., Y1, Y2, ... Yt, Yt+1 ... } = {Yt}_t=-INF..INF

    Observed series, length T
    { y1, y2, ... y_T } = {y_t}_t=1..T


### Stationary processes

* {Yt} is stationary if ALL aspects of behavior are UNCHANGED by shifts in time
    * Y may have correlation across time
    * but all Y has the exact same correlation across all time
* A stochastic process is **strictly stationary** if..
    * for any given finite r
    * and any set of subscripts t1, t2 ... tr
    * the joint distribution of (Yt, Yt1, Yt2, ... Ytr)
    * depends ONLY on t1 - t, t2 - t, tr - t
    * but does NOT depend on the specific value of t
    * E.g. the distribution of (Y1, Y5) is the same as (Y12, Y16)
    * Yt has the same mean, variance, and moments for all t
    * any transformation of a strictly stationary process...
        * e.g. {g(Yt)} = {Yt^2}
        * is also strictly stationary


## <a name="accac"></a>Auto-correlation and Auto-covariance

Covariance (Weakly) Stationary Processes: {Yt}

    E[Yt] = mu,                 for all t

    Var(Yt) = sigma^2,          for all t

    Cov(Yt, Yt-j) = gamma_j     depends on j, not t


* `Cov(Yt, Yt-j) = gamma_j` is called the **j-lag auto-covariance**
* auto-covariance measures DIRECTION 
    * but not STRENGTH 
    * of linear association
* auto-correlation measures both...
    * DIRECTION and STRENGTH 
    * of linear association

.

                                cov(Yt,Yt-j)             gamma_j
    corr(Yt, Yt-j) = rho_j = ------------------------ = ---------
                              sqrt(Var(Yt)*Var(Yt-j)     sigma^2


    By stationarity: var(Yt) = Var(Yt-j) = sigma^2

* Autocorrelation function (ACF)
    * Plot of `rho_j` vs j
    * j is the lag = number of periods between random variables
    * Yt, Yt-j
    * R: ARMAacf



## <a name="wnp"></a>White Noise Process

* **Gaussian White Noise**, {Yt}
    * represents random draws from N(0, sigma^2)
    * iid = independent and identically distributed

.

    Yt ~ iid N(0, sigma^2)

    Yt is independent of Ys for all t != s

    cov(Yt, Yt-s) = 0       for all t != s



* **Independent White Noise**
    * same as Guassian, but can use different dist than normal dist
        * e.g. student's t distribution
    * E[Yt] = 0
    * Var(Yt) = sigma^2
    * Yt independent of Ys for all t != s
* **Weak White Noise Process**
    * Yt ~ W N(0, sigma^2)
    * E[Yt] = 0
    * Var(Yt) = sigma^2
    * cov(Yt,Ys) = 0, for t != s
    * does NOT assume that Yt and Ys are INDEPENDENT
    * only assume Yt and Ys are UN-CORRELATED
        * INDEPENDENCE implies cov = 0
        * however, cov = 0 does NOT imply INDEPENDENCE
    * e.g. monthly return data tends to be uncorrelated
        * however the square of the return data tends to be correlated !


## <a name="nsp"></a>Non-Stationary Processes


* Stationary: common structure over time
* Non-stationary: something is different over time
    * "something" could be the mean, variance, or whatever
   

#### Example: deterministically trending process:

    Yt = beta_0 + beta_1 * t + e_t

    error_t ~ W N(0, sigma_e^2)

    E[Yt] = beta_0 + beta_1 * t     depends on t

Note: simple "detrending" transformation yields a stationary process:

    Xt = Yt - beta_0 - beta_1 * t = e_t

* subtract out the time-dependent part
* e.g. using a mean=0 expected return
    * unrealistic for long time intervals
    * more realistic for short time intervals
        * e.g. 1 day or 1 hour
    * can use mean=0 expected return for short time intervals
    * this can be useful, simplifies the math



#### Example: Random Walk

* Random Walk
    * named after the "drunken sailor" walking home from the bar
* non-stationary process
* the VARIANCE of the processes is what changes over time
    * gets larger over time 
    * because it becomes harder to predict further into future

.

    Yt = Yt-1 + e_t

    e_t ~ W N(0, sigma_e^2)

    Y0 is fixed.  Apply recursive substitution:

    Yt = Y0 + SUM_j=1..t [ e_j ]
    

    E[Yt] = Y0 + E[ SUM_j e_j ]

          = Y0  


    Var(Yt) = Var(Y0 + SUM_j e_j)

            = Y0 + Var(e_1) + Var(e_2) + ... + Var(e_t)

              Assume all covariances = 0

    Var(Yt) = sigma_e^2 * t     # depends on t



De-trending:

    delta_Yt = Yt - Yt-1 = e_t


## <a name="map"></a>Moving Average Processes

* MA Process
* Yt is correlated with Yt-1
    * but NOT correlated with t-2, t-3, ...
* Is it covariant stationary?
    * mean is constant over time
    * var is constant over time
    * time-dependence only depends on the size of the lag
    * MA(1) model IS covariant stationary
        * mean is constant over time
        * var is constant over time
        * time-dependence only depends on size of lag
        

.

    MA(1) Model

    Yt = mu + e_t + theta * e_t-1

    e_t ~ N(0, sigma_e^2)


    E[Yt] = mu + E[e_t] + theta * E[e_t]

          = mu          # constant over time

    
    Var(Yt) = Var(mu + e_t + theta * e_t-1)

              Covar(e_t, e_t-1) = 0

            = Var(e_t) + Var(theta * e_t-1)

            = Var(e_t) + theta^2 * Var(e_t-1)

            = sigma_e^2 * (1 + theta^2)     # constant over time


    Cov(Yt, Yt-1) = Cov(mu + e_t + theta * e_t-1, mu + e_t-1 + theta * e_t-2)

                  = E[ (Yt - mu) * (Yt-1 - mu) ]

                  = E[ (e_t + theta * e_t-1) * (e_t-1 + theta * e_t-2) ]

                    Cov(e_t, e_t-1) = 0

                    Cov(e_t-1, e_t-1) != 0

                    Cov(e_t-1, e_t-2) = 0

                  = E[ e_t * e_t-1 +
                       e_t * theta * e_t-2 + 
                       theta * e_t-1^2 +
                       theta^2 * e_t-1 * e_t-2 ]

                  = E[ e_t * e_t-1 ] + 
                    theta * E[ e_t * e_t-2 ] + 
                    theta * E[e_t-1^2] + 
                    theta^2 * E[e_t-1 * e_t-2]

                     Cov(e_t, e_t-1) = E[e_t * e_t-1] - E[e_t] * E[e_t-1]

                                       E[e_t] = E[e_t-1] = 0

                                     = E[e_t * e_t-1]

                                     = 0     

                  = 0 + theta * 0 + theta * Var(e_t-1) + theta^2 * 0

    Cov(Yt, Yt-1) = theta * sigma_e^2 
    
                  = gamma_1


             gamma_1      theta * sigma_e^2
    rho_1 = ---------  = -------------------------
             sigma^2       sigma_e^2(1 + theta^2)

                              theta
                       = ---------------
                          (1 + theta^2)

                       
                       = direction and strength of 1-period time dependence

                       = depends on time lag, not on t

                       = covariant stationary !!

    Note: max(rho_1) = 1/2


    
    gamma_2 = 2-period time covariance 
    
            = Cov(Yt, Yt-2) 

            = E[ (Yt - mu) * (Yt-2 - mu) ]

            = E[ (e_t + theta * e_t-1) * (e_t-2 + theta * e_t-3) ]

            = 0     # cuz there's no e_t-1^2 term

    In general:

    gamma_j = Cov(Yt, Yt-j) = 0,  for all j > 1 in our MA(1) model


### Example: MA(1) model for overlapping CC monthly returns

    r_t = 1-month cc return ~ N(mu,sigma^2)

    Consider a MONTHLY time series of 2-MONTH cc returns

    r_t(2) = r_t + r_t-1
    r_t-1(2) = r_t-1 + r_t-2
    r_t-2(2) = r_t-2 + r_t-3

    The two-month returns OVERLAP.

    Cov(r_t(2), r_t-1(2)) != 0

    rho_1 = 0.5 (can be shown)

    The stochastic process {r_t(2)} follows a MA(1) process


* Note that there's no time dependence in the underlying data, rt
* However, we've created time dependence...
    * by the way we've structured the data (overlapping)


## <a name="arp"></a>Auto-Regressive Processes

    AR(1) Model (mean-adjusted form)

    Yt - mu = phi * (Yt-1 - mu) + e_t

        e_t ~ iid N(0,sigma_e^2)

    -1 < phi < 1: AR(1) model is covariance stationary

    if phi=1, this becomes the Random Walk model (non-stationary)

    geometric decay: cor(Yt, Yt-1) > cor(Yt, Yt-2)
    correlation decreases over time

    Note: if phi < 0, correlations can change sign 
          e.g cor(Yt, Yt-1) < 0
              cor(Yt, Yt-2) > 0


* **Ergodicity**: time-dependence dies out as time interval increases
* e.g. Auto-regressive process
    * time-dependence decays geometrically
* Covariant stationary + ergodicity processes shows **Mean-Reversion**
    * the process tends to revert toward the mean
    * **Speed of mean-reversion** in Auto-Regressive processes
        * depends on phi
        * if phi is close to 1...
            * takes longer to revert to mean
            * close to random walk (which is NON-mean-reverting)
        * if phi is close to 0
            * takes shorter to revert to mean
            * close to white noise
        * if phi > 1...
            * "explosive" process
            * quickly trends toward infinity



#### Properties of AR(1) model:

    E[Yt] = mu
    
    Var(Yt) = sigma^2 

            = sigma_e^2 / (1 - phi^2)

    Cov(Yt, Yt-1) = gamma_1 
                
                  = sigma^2 * phi

    Cor(Yt, Yt-1) = rho_1 
    
                  = gamma_1 / sigma^2 
                  
                  = phi

    Cov(Yt, Yt-j) = gamma_j 
    
                  = sigma^2 * phi^j

    Cor(Yt, Yt-j) = rho_j 
    
                  = gamma_j / sigma^2 
                  
                  = phi^j

    Since -1 < phi < 1,

    lim_j..INF [ rho_j ] = phi^j = 0

* AR model is good for
    * interest rates
    * growth rates
        * GDP
        * money, velocity
        * real wages
        * unemployment rate
    * if there is time-dependence in the data,
        * then you can do forecasting
* not good for
    * stock returns



