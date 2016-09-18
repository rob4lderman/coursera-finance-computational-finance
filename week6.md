
# Coursera: Intro to Computational Finance: Week 6: Constant Expected Return Model

* [Constant Expected Return Model](#cerm)
* [Regression Model Representation (CER Model)](#rmr)
* [Monte Carlo Simulation](#mcs)
* [Monte Carlo Simulation: Multivariate Returns](#mcsmr)
* [Random Walk Model (RW)](#rwm)
* [Estimating Parameters of CER model](#epocerm)
* [Properties of Estimators](#poe)
* [Mean Squared Error](#mse)
* [Standard Errors](#se)
* [Asymptotic Properties of Estimators](#apoe)
* [Central Limit Theorem](#clt)
* [Confidence Intervals](#ci)
* [Stylized Facts for the Estimation of CER Model Parameters](#sfeocermp)
* [Monte Carlo Simulation to Evaluate CER Model Parameters](#mcstecermp)
* [Value-at-Risk in the CER Model](#varcerm)


## <a name="cerm"></a>Constant Expected Return Model

* Assumptions
    * normal distribution of CC monthly returns
    * covariance stationarity
        * mean, var constant over time
        * "constant expected return" = constant mean
    * recall (week1): multi-period CC return is equal to sum of single period CC returns
        * we like ADDITIVE (from week 1)!
        * sum of two normally distributed random variables...
            * is ALSO a random variable
    * recall (week2): simple returns are not suited for normal distribution
        * because it allows for impossible values
            * e.g. `R < -1`
* Note: we use monthly returns because...
    * they most closely resemble a normal distribution 
    * (as opposed to daily returns).
    * TODO: WHY?!
        * daily returns gives us a larger sample size
        * CER Model estimators become more accurate
* Essentially we are treating returns as...
    * random draws from a multi-variate normal distribution

.

    r_it = cc return on asset i in month t

    i = 1 .. N assets
    t = 1 .. T months

    r_it ~ N(mu_i, sigma_i^2), for all i and t

    mu_i = E[r_it]                  (constant over time - covariant stationarity)

    sigma_i^2 = var(r_it)           (constant over time)

    sigma_ij = cov(r_it, r_jt)      (constant over time)

    rho_ij = cor(r_it, r_jt)        (constant over time)


## <a name="rmr"></a>Regression Model Representation (CER Model)

* The CER Model can be alternative represented as a regression model
* the monthly return, `r_it`...
    * is the mean return 
    * plus some error term
        * error is gaussian white noise
        * with sd equal to the sd of the asset returns
    * the error terms for two assets are correlated
        * but not correlated over time

.

    r_it = mu_i + e_it

    e_it ~ iid N(0, sigma_i^2)      

    cov(e_it, e_jt) = sigma_ij      : covariance between assets

    cor(e_it, e_jt) = rho_ij        

    # errors are NOT correlated over time
    cov(e_it, e_js) = 0             : for t != s, for all i,j
    cov(e_it, e_is) = 0             : an asset's errors are NOT correlated over time


    
* Interpretation of `e_it`:
    * random news
    * news affecting asset i...
        * might be correlated with news affecting asset j
    * news is UNcorrelated over time

.

    e_it = r_it - mu_i

    unexpected = actual -  expected
    news         return    return

    No news,    e_it = 0,   r_it = mu_it
    Good news,  e_it > 0,   r_it > mu_it
    Bad news,   e_it < 0,   r_it < mu_it


### CER Regression Model with Standardized News Shocks

    r_it = mu_i + e_it

    e_it = sigma_i * z_it

    z_it ~ iid N(0,1)

    cov(z_it, z_jt) = cor(z_it, z_jt) 
    
                    = rho_ij

    cov(z_it, z_js) = 0     # uncorrelated over time

    z_it = "standardized" random news shock (+ for good news, - for bad news)

    sigma_i = magnitude of effect that standardized news shock has on asset returns


* This representation is useful for going from...
    * a model with normally distributed returns
    * to a model with NON-normally distributed returns
    * just plug in the appropraite distribution for `z_it`


### CER Model in Matrix Notation

    Nx1 vectors for returns, expected returns, and errors:

    r_t = [ r_1t, r_2t, ... r_Nt ]'

    mu = [ mu_1, mu_2, ... mu_N]'

    e_t = [ e_1t, e_2t, ... e_Nt]'


    NxN symmetric covariance matrix:

    Sigma = [ sigma_1^2  sigma_12   ...  sigma_1N  ]
            [ sigma_12   sigma_2^2  ...  sigma_2N  ]
            [ ...        ...             ...       ]
            [ sigma_1N   sigma_2N   ...  sigma_N^2 ]


    r_t = mu + e_t

    e_t ~ GWN(0,Sigma)

    This implies:  r_t ~ iid N(mu, Sigma)


## <a name="mcs"></a>Monte Carlo Simulation

* Use random number generator
    * to create simulated values
    * for a model
    * values pulled from a random variable distribution
* Run many simulations
    * to find the "typical" behavior
* Useful for:
    * reality check on the model
        * do simulated data look like real data?
    * create "what if?" scenarios
    * study properties of statistics 
        * computed from proposed model

Example: simulate monthly returns on MSFT from CER Model

    # Specify model parameters based on sample statistics:
    mu_i = 0.03         # monthly expected return
    sigma_i = 0.11      # monthly sd

    # r_it = mu_i + e_it, t=1,2,...100    (simulate 100 iterations)

    # e_it = iid N(0, sigma_i^2)

    # Simulation requires generating random numbers from a normal
    # distribution in order to determine e_it.

    set.seed(111)
    e_it = rnorm(100, mean=0, sd=sigma_i)
    r_it = mu_i + e_it

    # TODO: Four-panel plot
    hist(r_it, probability=T)
    boxplot(r_it)
    plot(density(r_it))
    qqnorm(r_it)
    qqplot(r_it)


## <a name="mcsmr"></a>Monte Carlo Simulation: Multivariate Returns

* Simulate using random values 
    * from multi-variate normal distribution

.

    mu = [ mu_sbux  ]
         [ mu_msft  ]
         [ mu_sp500 ]

    Sigma = [ sigma_sbux^2      sigma_sbux_msft  sigma_sbux_sp500  ]
            [ sigma_msft_sbux   sigma_msft^2     sigma_msft_sp500  ]
            [ sigma_sp500_sbux  sigma_sp500_msft sigma_sp500^2     ]


    r_t = mu + e_t

    e_it ~ iid N(0, sigma_i^2)

    cov(e_it, e_jt) = sigma_ij

    e_t = rmvnorm(100, sd=Sigma )
    


## <a name="rwm"></a>Random Walk Model (RW)

* CER Model for cc returns is EQUIVALENT to...
    * the random walk (RW) model for log stock prices

.

    r_t = ln(Pt/Pt-1) 

        = ln(Pt) - ln(Pt-1)

    ln(Pt) = ln(Pt-1) + r_t

           = ln(P0) + SUM_s=1..t r_s

            r_s = mu + e_s

           = ln(P0) + t*mu + SUM_s=1..t e_s


* Interpretation: Log price at time=t...
    * equals initial price, `ln(P0)`
    * plus expected growth in prices, `E[ln Pt] = t*mu`
    * plus accumulation of "news", `SUM_s=1..t e_s`

.

    Pt = exp(ln(Pt))

       = exp(ln(P0) + t*mu + SUM_s=1..t e_s)

       = P0 * exp(t*mu) * exp(SUM_s=1..t e_s)

       = initial  * expected   * unexpected
         price      growth in    growth in
                    price        price

    # Plot: PO * exp(t*mu)  - expected growth (linear)
    # Plot: exp(SUM_s=1..t e_s) - random growth (non-linear)
    # Plot: plot1 + plot2

* Note: prices are non-stationary
    * Pt+1 depends on Pt
* returns are stationary
* random-walk variance increases over time
    * `var( SUM_s=1..t e_s ) = sum of variances`


## <a name="epocerm"></a>Estimating Parameters of CER model


Parameters:

    mu_i            expected return for asset i
    sigma_i^2       variance of asset i
    sigma_ij        covariance of assets i and j
    rho_ij          correlation of assets i and j

* Estimate using observed sample of historical data
* **Estimator**: a method/formula for generating an estimate
    * a random variable
* **Estimate**: one particular value generated by the estimator
    * one outcome from the random variable
* **Plug-in principle**: estimate model paramters using sample statistics
    * easy and intuitive way to think about estimators
    * for CER Model, same as least-squares estimate from regression model
* **Maximum Likelihood**: another method for estimating parameters
    * more general technique
    * for CER Model, maximum likelihood is same as plug-in principle

.

    # returns.mat: matrix of SBUX, MSFT, SP500
    muhat.vals = apply(returns.mat, 2, mean)

    sigma2hat.vals = apply(returns.mat, 2, var)

    sigmahat.vals = apply(returns.mat, 2, sd)

    cov.mat = var(returns.mat)
    cor.mat = cor(returns.mat)

    covhat.vals = cov.mat[ lower.tri(cov.mat) ]
    rhohat.vals = cor.mat[ lower.tri(cor.mat) ]

    names(covhat.vals) = names(rhohat.vals) = c("sbux,msft", "sbux,sp500", "msft,sp500")


## <a name="poe"></a>Properties of Estimators

    theta = parameters to be estimated
    ^theta = estimator of theta from random sample

* `^theta` is a random variable
* pdf of `^theta` depends on pdf of random variables in random sample
* props of `^theta` can be derived analytically
    * using probability theory
    * or estimated using central Limit Theorem
        * avgs of random variables
    * or estimated using monte carlo simulation
        * to get distribution of underlying random variables

#### Bias

    bias(^theta, theta) = E[^theta] - theta

    ^theta is UN-biased if E[^theta] = theta

#### Precision

    SE(^theta) = standard error of ^theta
            
               = sqrt(var(^theta))

               = sqrt( E[(^theta - E[^theta])^2]

               = sigma_^theta

* If an estimator `^theta` is unbiased...
* its distribution is centered around "the truth", `theta`
* its "spread" is equal to its standard error


### Properties of Estimators for CER Model

* mean is unbiased
* variance is unbiased
* covariance is unbiased
* sd is BIASED
* correlation is BIASED
* BIAS tends to be small
    * goes to 0 as sample size T gets large
* use bootstrapping to estimate size of BIAS


## <a name="mse"></a>Mean Squared Error
    
MSE = Square the diff between the estimated value and the true value

    MSE(^theta, theta) = E[ (^theta -theta)^2 ]

                       = Var(^theta) + bias(^theta,theta)^2


## <a name="se"></a>Standard Errors

* we want the estimation error to be small
* SE quantifies the estimation error
* SE is the SD of the distribution of the estimator
    * tells you how far from "the truth" the estimate typically is

.

    SE(^mu_i) = sigma_i / sqrt(T)

    SE(^sigma_i^2) ~= sigma_i^2 / sqrt(T/2)

    SE(^sigma_i) ~= sigma_i / sqrt(2*T)
    # Note: SE(^sigma_i) is smaller than SE(^mu) by factor of sqrt(2)
    # i.e. we estimate volatility much better than mean

    SE(^sigma_ij):  no easy formula

    SE(^rho_ij) ~= (1 - rho_ij^2) / sqrt(T)


* Large SE means imprecise estimate
* precision increases with sample size
    * `SE -> 0 as T -> INF`
* `^sigma_i` is generally more precise than..
    * `^mu_i`
    * `^rho_ij`
* SE formulas are approximations for:
    * `^sigma_i`
    * `^rho_ij`
    * monte carlo simulations and bootstrapping can be used to get better approximations
* SE formulas depend on **unknown values of parameters**
    * are not practically useful
    * e.g. we don't know true `sigma_i`
    * e.g. we don't know true `rho_ij`
* Practically useful formulas compute ESTIMATED SE
    * replace unknows with estimates
* Note: there are errors to the estimate of the SE
    * SE of SE
    * infinite recursion
    * so we don't usually care about error in the SE estimate

.

    se.muhat = sigmahat.vals / sqrt(nobs)
    se.sigma2hat = sigma2hat.vals / sqrt(nobs/2)
    se.sigmahat = sigmahat.vals / sqrt(2*nobs)

    se.rhohat = (1-rhohat.vals^2)/sqrt(nobs)


* Is the SE big?
    * compare magnitude of SE to the estimate




## <a name="apoe"></a>Asymptotic Properties of Estimators

* **Asymptotic properties**
    * properties that result as the sample size goes to INF
* **Consistency** 
    * estimate converges on truth as sample size goes to INF
    * `bias(^theta,theta) = 0, as T -> INF`
    * `SE(^theta) = 0, as T -> INF`
    * In CER Model, all estimators are consistent
        * `^mu_i` is normally distributed
            * `^mu_i ~ N( mu_i, sigma_i^2/T )`
        * exact distributions of `^sigma_i`, `^sigma_ij`, `^rho_ij` not normal
            * however as sample size T goes to INF
            * distributions approach normal
            * due to Central Limit Theorem
        


### <a name="clt"></a>Central Limit Theorem

    Let X1, ... XT be iid random variables, with:

        E[Xt] = mu

        var(Xt) = sigma^2

    Then,

        X' - mu        X' - mu
       --------- = --------------- ~ N(0,1) as T -> INF
        SE(X')      sigma/sqrt(T)

        (note: X' is "X-bar", not "X-transpose")

    Equivalently,

        X' ~ N(mu, SE(X')^2) ~ N(mu, sigma^2/T)

        for large enough T

    "X' is asymptotically normally distributed with mean mu and variance sigma^2/T."


### <a name="ci"></a>Confidence Intervals

* A range of estimated values...
    * that we can say with a certain probability
    * contains the true value 
    * i.e. probability that the interval covers the parameter
        * NOT the probability of the parameter being in the interval
        * the interval is what's RANDOM
        * the parameter is NOT random
            * it has some fixed (unknown) value
* An approximate 95% confidence interval for `theta`:
    * `^theta +/- 2 * SE(^theta)`
    * the range will include the true value, `theta`
    * with a probability approximately equal to 0.95
    * i.e, `Pr( ^theta-2*SE(^theta) <= theta <= ^theta+2*SE(^theta) ) ~= 0.95`
* An approximate 99% confidence interval for `theta`:
    * `^theta +/- 3 * SE(^theta)`
    * the range will include the true value, `theta`
    * with a probability approximately equal to 0.99
    * i.e, `Pr( ^theta-3*SE(^theta) <= theta <= ^theta+3*SE(^theta) ) ~= 0.99`

95% conf intervals for mu, sigma, rho:

    mu.lower = muhat.vals - 2 * se.muhat
    mu.upper = muhat.vals + 2 * se.muhat
    mu.width = mu.upper - mu.lower

    sigma.lower = sigmahat.vals - 2 * se.sigmahat
    sigma.upper = sigmahat.vals + 2 * se.sigmahat
    sigma.width = sigma.upper - sigma.lower

    rho.lower = rhohat.vals - 2 * se.rhohat
    rho.upper = rhohat.vals + 2 * se.rhohat
    rho.width = rho.upper - rho.lower


## <a name="sfeocermp"></a>Stylized Facts for the Estimation of CER Model Parameters

* The mean is not estimated very precisely
    * Large SE relative to size of mean estimates
* Standard Deviations and correlations are estimated more precisely than mean


## <a name="mcstecermp"></a>Monte Carlo Simulation to Evaluate CER Model Parameters

* Evaluate:
    * bias
    * SE
    * Confidence INterval coverage
* Procedure:
    * create many simulated samples from CER model
        * specify true values for mean, variance, etc
    * compute parameter estimates from each simulated sample
    * compute mean and SD (SE) of estimates over simulated samples
    * compute 95% conf interval
    * count the number of intervals that cover the true parameter

#### Evaluating mean, var, sd: 

    # true values
    mu = 0.05
    sd = 0.10
    n.obs = 100
    n.sim = 1000

    sim.means = rep(0,n.sim)
    sim.vars = rep(0,n.sim)
    sim.sds = rep(0,n.sim)

    # generate simulated values
    set.seed(111)
    for (sim in 1:n.sim) {

        # generate simulated sample: n.obs randoms from normal dist
        sim.ret = rnorm(n.obs, mean=mu, sd=sd)

        # compute sample mean, var, sd
        sim.means[sim] = mean(sim.ret)
        sim.vars[sim] = var(sim.ret)
        sim.sds[sim] = sqrt(sim.vars[sim])
    }

    # estimate of bias
    mean(sim.means) - mu
    mean(sim.vars) - sd^2
    mean(sim.sds) - sd

    
    sd(sim.means)       # monte-carlo standard error estimate
    sd/sqrt(n.obs)      # analytic (true) SE estimate

    sd(sim.vars)
    sd^2/sqrt(n.obs/2)

    sd(sim.sds)
    sd/sqrt(2*n.obs)


#### Evaulation of 95% Confidence Interval coverage:

    mu.lower = rep(0,n.sim)
    mu.upper = rep(0,n.sim)

    # generate simulated values
    set.seed(111)
    for (sim in 1:n.sim) {

        # generate simulated sample: n.obs randoms from normal dist
        sim.ret = rnorm(n.obs, mean=mu, sd=sd)

        # compute 95% conf int
        mu.hat = mean(sim.ret)
        se.mu.hat = sd(sim.ret) / sqrt(n.obs)

        mu.lower[sim] = mu.hat - 2 * se.mu.hat
        mu.upper[sim] = mu.hat + 2 * se.mu.hat
    }

    # determine how many conf intervals actually covered the true value
    in.interval = (mu >= mu.lower) & (mu <= mu.upper)
    sum(in.interval) / n.sims
    # [1] 0.934


#### Evaluate distribution of correlations


    
    n.obs = 100
    n.sim = 1000

    sim.corrs = matrix(0,n.sim,3)   
    colnames(sim.corrs) = c("sbux,msft", "sbux,sp500", "msft,sp500")

    set.seed(111)
    for (sim in 1:n.sim) {
        sim.ret = rmvnorm(n.obs, mean=muhat.vals, cov=cov.mat)
        cor.mat = cor(sim.ret)
        sim.corrs[sim,] = cor.mat[lower.tri(cor.mat)]
    }


    # compare monte carlo results to true values for rhohat
    apply(sim.corrs,2,mean)

    # compare to analytic SE values for rhohat
    apply(sim.corrs,2,sd)




## <a name="varcerm"></a>Value-at-Risk in the CER Model

    In the CER Model:

        r_it ~ iid N(mu_i, sigma_i^2)

        r_it = mu_i + sigma_i * z_it

        z_it ~ iid N(0,1)

    The alpha quantile may be expressed as:

        q_alpha(r) = mu_i + sigma_i * q_alpha(z)

        q_alhpa(z) = standard normal quantile

    Estimated quantiles:

        ^q_alpha(r) = ^mu_i + ^sigma_i * q_alpha(z)

    Estimated V-a-R:

        ^VaR_alpha = (exp(^q_alpha(r)) - 1) * W0

            ^q_alpha(r) = ^mu_i + ^sigma_i * q_alpha(z)

            W0 = initial investment

R:

    qhat.05 = muhat.vals + sigmahat.vals * qnorm(0.05)
    
    W0 = 100000
    VaRhat.05 = (exp(qhat.05)-1) * W0


### Standard Error for Value-at-Risk estimate

We can compute `SE(^q_alpha(r))` using:

    ^q_alpha(r) = ^mu_i + ^sigma_i * q_alpha(z)

    var(^q_alpha(r)) = var(^mu_i) + q_alpha(z)^2 * var(^sigma_i) + 2 * cov(^mu_i,^sigma_i)

        cov(^mu_i,^sigma_i) = 0     # this can be shown analytically

    SE(^q_alpha(r)) = sqrt( var(^mu_i) + q_alpha(z)^2 * var(sigma_i) )

However, computing `SE(^VaR_alpha)` is not straightforward:

    var(^VaR_alpha) = var( (exp(^q_alpha(r)) - 1) * W0 )

* Computing `var(^VaR_alpha)` is hard.  
    * assumes large sample size
* There's an easier way: **bootstrapping!!**
