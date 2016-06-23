


# Coursera: Intro to Computational Financing: Week 7: Bootstrapping and Hypothesis Testing

* [Bootstrapping](#boot)
* [Bootstrapping in the CER Model](#bootcer)
* [Bootstrapping Procedure](#bootproc)
* [Bootstrap 95% Confidence Intervals](#bootci)
* [Performing Bootstrap in R:](#bootr)
* [Bootstrapping Value-at-Risk](#bootvar)
* [Hypothesis testing](#hyptest)
* [Power of Test](#power)
* [Hypothesis Testing in CER Model](#hyptestcer)
* [Chi-square random variable and distribution](#chisq)
* [Student's t random variable and distribution](#studentst)
* [Example: Hypothesis Test for Specific Coefficient Value](#hyptestvalue)
* [Example: Test for Normal Distribution](#hyptestnormal)
* [Example: Test for No Autocorrelation](#hyptestautocor)
* [Example: Diagnostics for Constant Parameters (Rolling Avgs)](#hyptestparams)
* [Summary of Hypothesis Testing in CER Model](#summaryhyptest)


## <a name="boot"></a>Bootstrapping

* Bootstrapping: Simulate many values to estimate stats like SE
    * resampling methods
    * fewer assumptions
        * does not rely on normal distributions
    * greater accuracy
        * does not rely on large sample sizes
        * in contrast to CLT
    * generality
        * same method applies to wide variety of stats
* avoids having to do messy analytical mathematics
    * which are often approximations
    * that rely on large samples and CLT


### <a name="bootcer"></a>Bootstrapping in the CER Model

    CER Model:

        r_t = mu + e_t

        e_t ~ iid N(0, sigma^2)


    Observed sample: {r_1, r_2, ... r_T}

    Goals: Compute using bootstrap: 
        1. ^SE(^mu) 
        2. ^SE(^sigma)
        3. 95% confidence intervals for mu and sigma

    Recall: The analytical formulas are:

        ^SE(^mu) = ^sigma/sqrt(T)

        ^SE(^sigma) ~= ^sigma/sqrt(2T)

        95% conf int ~= ^theta +/- 2 * ^SE(^theta)

                ^theta = ^mu, ^sigma


### <a name="bootproc"></a>Bootstrapping Procedure

"Non-parametric Bootstrapping": no parameters, no models

1. Resample (w/ replacement)
    * create B bootstrap samples from the original observed data
    * Each bootstrap sample has T observations
        * same size as original sample
    * very similar to monte carlo method
        * monte carlo: sample randomly from an assumed distribution
        * bootstrap: sample from observed data
    * resampling w/ replacement OK if observations are uncorrelated
    * if observations are correlated...
        * correlation must be preserved in bootstrap sample
        * e.g. if adjacent observations are correlated,
            * then select adjacent observations when bootstrapping
2. Calculate bootstrap distribution for statistic of interest, `^theta`
    * for each bootstrap sample, compute `^theta*`
    * there will be B values of `^theta*`
3. Use the bootstrap distribution
    * gives info about the shape, center, and spread...
    * of the unknown pdf `f(^theta)`



    Bootstrap samples:

        {r*_11, r*_12, ... r*_1T} = 1st bootstrap sample
            ...
        {r*_B1, r*_B2, ... r*_BT} = B'th bootstrap sample


    B bootstrap values of the estimated parameter, ^theta:

        {theta*_1, ^theta*_2, ... ^theta*_B}


    Bootstrap estimate of bias:  bootstrap mean - estimate

        E[^theta*] =    1/B * SUM_j=1..B ^theta*_j 

        bias_boot(^theta, theta) = E[^theta*] - ^theta


    Bootstrap estimate of ^SE(^theta): sample SD of bootstrap values

                         1
        Var(^theta*) = -----  SUM_j=1..B (^theta*_j - E[^theta*])^2
                        B-1


        ^SE_boot(^theta) = sqrt( Var(^theta*) )




### <a name="bootci"></a>Bootstrap 95% Confidence Intervals

1. If bootstrap distributions is symmetric and looks normal, use:
    * `^theta +/- 2 * ^SE_boot(^theta)`
2. If bootstrap distribution is NOT symmetric and looks NON-normal, use:
    * `[q*_.025, q*_.975]`
    * `q*_.025` = 2.5% quantile from bootstrap distribution
    * `q*_.975` = 97.5% quantile from bootstrap distribution


## <a name="bootr"></a>Performing Bootstrap in R:

Two options:

1. brute force
2. R package boot


### Brute force

Perform steps by hand:

1. sample w/ replacement
    * using `sample()`
2. compute statistic from B bootstrap samples
    * using for loops and `apply()`
3. compute bootstrap bias and SE
    * using `mean()` and `sd()` functions


Recall: Estimated Standard Errors:

    se.muhat = sigmahat.vals/sqrt(nobs)
    rbind(muhat.vals, se.muhat)

    se.sigma2hat = sigma2hat.vals / sqrt(nobs/2)
    rbind(sigma2hat.vals, se.sigma2hat)

    se.sigmahat = sigmahat.vals/sqrt(2*nobs)
    rbind(sigmahat.vals, se.sigmahat)


Brute force Bootstrap:

    B = 999     # choose odd number for quantiles for conf int
    muhat.boot = rep(0,B)
    nobs = length(MSFT)

    # generate bootstrap samples via resampling
    for (i in 1:B) {
        boot.data = sample(MSFT, nobs, replace=T)
        muhat.boot[i] = mean(boot.data)
    }

    # bootstrap bias
    mean(muhat.boot) - muhat.MSFT

    # bootstrap SE
    sd(muhat.boot)

    # analytic SE
    sigmahat.MSFT / sqrt(length(MSFT))


    # plot it
    par(mfrow=c(1,2))

    hist(muhat.boot, col="slateblue1")
    abline(v=muhat.MSFT, col="white", lwd=2)
    qqnorm(muhat.boot)
    qqline(muhat.boot)

    par(mfrow=c(1,1))


### R package boot    


Implements variety of bootstrapping functions

    boot()          # bootstrap user-supplied function
    boot.ci()       # bootstrap conf int
    

#### Example: bootstrapping sample mean

    # function for bootstrapping sample mean
    mean.boot = function(x, idx) {
        # arguments:
        # x:    data to be resampled
        # idx:  vector of scrambled indices created by boot() function
        #
        # return: mean value computed using resampled data

        ans = mean(x[idx])
        ans
    }


    MSFT.mean.boot = boot(MSFT, 
                          statistic = mean.boot,
                          R=999)
    class(MSFT.mean.boot)
    # [1] "boot"

    MSFT.mean.boot
    # ..
    #
    # Bootstrap Statistics:
    #       original            bias        std. error
    # t1* {sample-mean}     {bootstrap      {bootstrap
    #                        estimate        estimate
    #                        of bias}        of SE}

    plot(MSFT.mean.boot)
    # plots histogram and qqplot against the normal

    # compare bootstrap std.error vs. analytical estimate
    se.muhat.MSFT = sigmahat.MSFT / sqrt(length(MSFT))


#### Bootstrap Confidence Interval:

    boot.ci(MSFT.mean.boot, 
            conf=0.95, 
            type=c("norm", "perc"))     # "norm": use normal
                                        # "perc": use quantiles


#### Example: Bootstrapping sample SD

    sd.boot = function(x, idx) {
        sd(x[idx])  # sample sd
    }

    MSFT.sd.boot = boot(MSFT,
                        statistic = sd.boot,
                        R=999)

    # compare boot SE with analytic SE based on CLT:
    MSFT.sd.boot

    se.sigmahat.MSFT = sigmahat.MSFT / sqrt(2*length(MSFT))

    # plot it
    plot(MSFT.sd.boot)
    

## <a name="bootvar"></a>Bootstrapping Value-at-Risk

VaR in CER Model:

    ^VaR_.05 = (exp(^q_.05) - 1) * W0

    ^q_.05 = ^mu + ^sigma * q_.05(z)

Bootstrapping can be used to compute:

    ^SE(^VaR_.05)

as well as confidence intervals.


#### Example: Bootstrapping NOrmal VaR


    valueAtRisk.boot = function(x, idx, p=0.05, w=100000) {
        # Arguments:
        #   x:      data to be resampled
        #   idx:    vector of scrambled indices created by boot() function
        #   p:      probability value for VaR calculation
        #   w:      value of initial investment
        #
        # Returns: Value-at-Risk computed using resampled data
        #   

        q = mean(x[idx]) + sd(x[idx]) * qnorm(p)
        # or: q = qnorm(p, mean=mean(x[idx]), sd=sd(x[idx]))
        VaR = (exp(q) - 1 ) * w
        VaR
    }

    MSFT.VaR.boot = boot(MSFT, 
                         statistic = valueAtRisk.boot,
                         R=999,
                         p=0.05,
                         w=100000)

    plot(MSFT.VaR.boot)

    boot.ci(boot.out = MSFT.VaR.boot,
            conf = 0.95,
            type=c("norm", "perc"))



## <a name="hyptest"></a>Hypothesis testing

Procedure:

1. Specify hypothesis to be tested
    * H0: null hypothesis vs.
    * H1: alternative hypothesis
    * Some Finance tests:
        * is return data normally distributed?
        * are CER model parameters constant over time?
            * e.g. mean
            * correlations
2. Specify significance level of test
    * level = Pr(Reject H0|H0 is true)
    * level = Pr(Type 1 Error)
        * Type I error: detect effect (H1) when no effect exists (H0)
        * Type II error: detect NO effect (H0) when an effect exists (H1)
    * common levels: 
        * 5% `(p<0.05)`
        * 1% `(p<0.01)`
3. Construct test statistic, T, from observed data
4. Use test statistic T to evaluate data evidence regarding H0
    * |T| is big : evidence against H0
    * |T| is small : evidence in favor of H0 
    * **p-values**: probability of observing test statistic T

R: 

    pt(T, df=x)             = probabiliity of observing t-statistic <= T
    1 - pt(abs(T), df = x)  = probability of observing t-statistic > abs(T)

    # two-sided: 
    2 * (1 - pt(abs(T), df=x)


### <a name="power"></a>Power of Test

    Power = 1 - Pr(Type II Error)

          = Pr(Reject H0|H0 is false)

* Goal: construct test with high power
* Problem: Impossible to simultaneously have..
    * level ~= 0
    * power ~= 1
    * As `level -> 0`, `power -> 0`


## <a name="hyptestcer"></a>Hypothesis Testing in CER Model

Example: Test for specific value.

    H0: mu_i = *mu_i        H1: mu_i !=  *mu_i
    H0: sigma_i = *sigma_i  H1: sigma_i !=  *sigma_i
    H0: rho_ij = *rho_ij    H1: rho_ij !=  *rho_ij

Example: Test for sign

    H0: mu_i = 0            H1: mu_i != 0
    H0: rho_ij = 0          H1: rho_ij != 0

Example: Test for normal distribution

    H0: r_it ~ iid N(mu_i, sigma_i^2)
    H1: r_it ~ NOT normal

Example: Test for no autocorrelation

    H0: p_j = corr(r_it, r_it-j) = 0,   j>1
    H1: p_j = corr(r_it, r_it-j) != 0,  for some j

Example: Test of constant parameters

    H0: mu_i, sigma_i, rho_ij are constant over entire sample
    H1: mu_i, sigma_i, rho_ij changes in some sub-sample



## <a name="chisq"></a>Chi-square random variable and distribution

    Let Z1 ... Zq be iid N(0,1) random variables

        X = Z1^2 + Z2^2 + ... Zq^2

    Then,

        X ~ chi^2(q)

        q = degrees of freedom (d.f.)

    Properties:

        X > 0   (squared terms)

        E[X] = q

        chi^2(q) -> Normal as q -> INF


    R:
    rchisq
    dchisq
    pchisq
    qchisq



## <a name="studentst"></a>Student's t random variable and distribution

    Let:
        Z ~ N(0,1)

        X ~ chi^2(q)

        Z and X are independent

    Then,

                Z
        T =  --------   ~ t_q
             sqrt(X/q)

        q = degrees of freedome


    Properties of t_q distribution:

        E[T] = 0

        skew(T) = 0

                  3q - 6
        kurt(T) = -------- ,  q > 4
                   q - 4

        T -> N(0,1) as q -> INF  (q >= 60)


    "t-statitistic": follows a Student's T distribution

* `t_T-1` is similar to Normal but with fatter tails
* d.f. = sample size - number of estimated parameters



### <a name="hyptestvalue"></a>Example: Hypothesis Test for Specific Coefficient Value


    H0: mu_i = *mu_i
    H1: mu_i != *mu_i

    1. Test statistic:

                            ^mu_i - *mu_i
        t[mu_i = *mu_i] = ------------------
                              ^SE(^mu_i)

        E.g. t = 2 means ^mu_i is two standard deviations (^SE(^mu_i))
                   away from *mu_i.


#### Distribution of t-statistic under H0

    Under the assumptions of the CER Model, and H0: mu_i = *mu_i

        t-statistic ~ t_T-1

        t-statistic follows a t distribution with T-1 degrees of freedom


    2. Set significance level:

            Pr(Type I error) = 5%

        Two-sided test, so critical values are at q_.025 and q_.975

    3. Decision rule

        Reject H0 if abs(t) > q_.975


#### Useful Rule of Thumb:

    If T >= 60, then q_.975 ~= 2, so:

        Reject H0 if abs(t) > 2


R:

    nobs = nrow(returns.z)
    muhat.vals = apply(returns.z,2,mean)
    sigmahat.vals = apply(returns.z,2,sd)
    se.muhat = sigmahat.vals/sqrt(nobs)

    t.stats = muhat.vals/se.muhat
    abs(t.stats)

    # p-values:
    2 * (1 - pt(abs(t.stats), df=nobs-1))

    
    t.test.msft = t.test(returns.z[,"msft"],
                         alternative="two.sided",
                         mu=0,
                         conf.level=0.95)
    class(t.test.msft)
    # [1] "htest"



#### Relationship between Hypothesis Test and Confidence Intervals

    Hypothesis Test Rule of Thumb:

        Reject H0 at 5% if abs(t) > 2

    Confidence Interval:

        Reject H0 at 5% if hypothesized value does NOT lie in 95% confidence interval



## <a name="hyptestnormal"></a>Example: Test for Normal Distribution


    H0: r_t ~ iid N(mu,sigma^2)

    H1: r_t ~ not normal


1\. Test statistic:  Jarque-Bera statistic

          T               (^kurt - 3)^2
    JB = --- * (^skew^2 + -------------- )
          6                     4
        
    R: tseries :: jarque.bera.test 



Intuition: 

    If      
        r_t ~ iid N(mu,sigma^2), 
    then    
        ^skew(r_t) ~= 0 and ^kurt(r_t) ~= 3, 
    and
        JB ~= 0


    ELSE 
        if r_t is NOT normally distributed,
    then
        ^skew(r_t) != 0 and/or ^kurt(r_t) != 3, 
    and
        JB >> 0
        
 
Distribution of JB under H0:

    if HO is true, then

        JB ~ chi^2(2)
   

2\. Set significance level and determine critical value:

    Pr(Type I error) = 5%

    Pr(chi^2(2) > cv) = 0.05

        cv = qchisq(0.95, df=2) ~= 6

3\. Decision Rule:

    Reject H0 at 5% level if JB > 6


4\. P-Value of test

    significance level at which test is just rejected

        Pr(chi^2(2) > JB)

R:

    sbux.skew = skewness(returns.z[,"sbux"])
    sbux.ekurt = kurtosis(returns.z[,"sbux"])

    JB = nobs / 6 * (sbux.skew^2 + sbux.ekurt^2 / 4)
    # [1] 24.34  -> reject H0

    p.value = 1 - pchisq(JB, df=2)

    library(tseries)
    jarque.bera.test(returns.z[,"sbux"])


## <a name="hyptestautocor"></a>Example: Test for No Autocorrelation


Recall the jth lag autocorrelation:

    rho_j = cor(r_t, r_t-j)

          = cov(r_t, r_t-j) / var(r_t)

Hypothesis:

    H0: rho_j = 0, for all j = 1, ... q

    H1: rho_j != 0, for some j

1\. Estimate rho\_j using sample autocorrelation:

    cov(r_t, r_t-j) = 1/T * SUM_t=j+1..T (r_t - ^mu) * (r_t-j - ^mu)

    var(r_t) = 1/T * SUM_t=1..T (r_t - ^mu)^2

    rho_j = cov(r_t, r_t-j) / var(r_t)
               

Recall: Based on CLT:

    ^rho_ij ~ N(rho_ij, SE(^rho_ij)^2)

                  1 - rho_ij^2
    SE(^rho_ij) = -------------
                    sqrt(T)

Result: Under H0, if T is large then:

    ^rho_j ~ N(0, 1/T), for all j >= 1

    SE(^rho_j) = 1/sqrt(T)


2\. Test Statistic

    t = ^rho_j / SE(^rho_j)  = sqrt(T) * ^rho_j

    95% conf interval: ^rho_j +/- 2 * 1/sqrt(T)

3\. Decision Rule

    Reject H0 at 5% level if abs(t) > 2

    ^rho_j > 2/sqrt(T)   or ^rho_j < -2/sqrt(T)

Note: the dotted blue lines on the ACF graph are +/- 2/sqrt(T)

    acf(returns.ts[,"sbux"])


## <a name="hyptestparams"></a>Example: Diagnostics for Constant Parameters

    H0: mu_i / sigma_i / rho_ij are constant over time
    H1: mu_i / sigma_i / rho_ij are NOT constant over time

* Formal test statistics are available but require advanced stats
    * see R package strucchange
* Informal graphical diagnostics:
    * Rolling estimates


### Rolling Means

Compute estimate of `mu_i` over rolling windows of length n < T


    ^mu_it(n) = 1/n * SUM_j=0..n-1 r_it-j

    R:  zoo :: rollapply

    H0: if mu_i is constant, then ^mu_it(n) should stay fairly constant


    roll.muhat = rollapply(returns.z[,"sbux"], 
                           width=24,
                           FUN=mean,
                           align="right")   # which time index (right-endpoint) to apply to rolling mean


### Rolling Variances and Standard Deviations

                      1
    ^sigma_it^2(n) = --- * SUM_j=0..n-1 (r_it-j - ^mu_it(n))^2 
                     n-1

    ^sigma_it(n) = sqrt( ^sigma_it^2(n) )


    roll.sigma2hat = rollapply(returns.z[,"sbux"],
                               width=24,
                               FUN=var,
                               align="right" )

    roll.sigmahat = rollapply(returns.z[,"sbux"],
                              width=24,
                              FUN=sd,
                              align="right" )



### Rolling Covariances and Correlations

                      1
    ^sigma_ij,t(n) = --- * SUM_k=0..n-1 (r_it-k - ^mu_it(n)) * (r_jt-k - ^mu_jt(n))
                     n-1

                        ^sigma_ij,t(n)
    ^rho_ij,t(n) = -----------------------------
                    ^sigma_it(n) * ^sigma_jt(n)


    rhohat = function(x) {
        # x: 2 column matrix of returns
        # return: correlation coefficient (rho) of the two cols of returns
        cor(x)[1,2]     # cor() returns cor matrix. [1,2] returns the correlation coefficient 
    }

    roll.rhohat = rollapply(returns.z[,c("sp500","sbux")],
                            width=24,
                            FUN=rhohat,
                            by.column=F,    # apply FUN to both columns at the same time, not each col individually
                            align="right")

* Correlations tend to vary widely over time
* especially during financial crises
    * assets tend to become highly correlated during crises
    * you thought your portfolio was diversified...
        * with uncorrelated assets...
    * and then a crisis hits
        * and all your "uncorrelated" assets become correlated
        * and your portfolio tanks
* TODO: Need to analyze asset correlations during crises
    * to find good stabilizers for crisis times



## <a name="summaryhyptest"></a>Summary of Hypothesis Testing in CER Model

* Can often reject hypothesis that monthly returns are "NOT" normally distributed
    * prof said "NOT"
    * slide omitted "NOT"
* Typically cannot reject hypothesis that monthly returns are "NOT" uncorrelated over time
    * doesn't appear to be time-dependence in the returns themselves
* Rolling window estimates show that mean, sd, and correlations are NOT constant over time
    * RED FLAG!
    * implications for portfolio theory

