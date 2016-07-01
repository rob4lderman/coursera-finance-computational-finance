

# Coursera: Intro to Computational Finance: Week 10: Portfolio Risk Budgeting, Single Index Model

* [Portfolio Risk Budgeting](#portriskbudget)
* [Risk Decomposition using Euler's Theorem](#riskdecomp)
* [Marginal Contribution to Risk (MCR)](#mcr)
* [Beta](#beta)
* [Sharpe Single Index Model](#simodel)
* [Statistical Properties of SI Model (Unconditional)](#simodelstats)
* [Statistical Properties of SI Model (Conditional on `R_Mt = r_Mt`)](#simodelstatscond)
* [Decomposition of Total Variance](#decompvar)
* [SI Model and Portfolios](#simodelport)
* [Estimating the SI Model](#estimatesimodel)
* [Least Squares Estimates](#lse)
* [SI Model Using Matrix Algebra](#simodelmatrix)
* [SI Model for 4 Asset Portfolio](#simodel4port)
* [Estimating SI Model Covariance Matrix](#simodelcovmat)
* [Hypothesis Testing in SI Model](#simodelhyptest)

## <a name="portriskbudget"></a> Portfolio Risk Budgeting

* getting more popular
* in terms of how to structure a portfolio
* how much risk is driven by a specific security(s)
* Additively decompose portfolio risk..
    * into contributions from individual assets
    * "risk budgets"
* risk measures:
    * var and sd are easily decomposed
    * other measures, e.g. VaR, not so easy
* for risk measures that are...
    * homogeneous functions of degree k=1 in the portfolio weights
        * portfolio return
        * portfolio expected return
        * portfolio sd
        * NOT portfolio var
        * VaR
        * expected short fall
        * expected tail loss
    * **Euler's theorem** provides a general method for decomposition


### Homogenous function of degree k: 

    f(alpha*x) = alpha^k * f(x)

    f(c * x1, c * x2, c * x3, ...) = c * f(x1, x2, x3, ...)
    
    f(x1,x2) = x1 + x2              : Homogeneous
    f(x1,x2) = x1^2 + x2^2          : Not homogeneous
    f(x1,x2) = sqrt(x1^2 + x2^2)    : Homogeneous


### Euler's Theorem

    Let f(x1, x2, ... xn) = f(x) 

    f(x) is continuous, differentiable, and homogeneous of degree one.

    Then f(x) can be additively decomposed as follows:

                 d-f(x)        d-f(x)
    f(x) = x1 * -------- + x2 -------- + ...
                 d-x1           d-x2
                            
                 d-f(x)
         = x' * -------
                  d-x


    d-f(x)   [ d-f(x)/d-x1 ]
    ------ = [ d-f(x)/d-x2 ]
     d-x     [     ...     ]
    (nx1)    [ d-f(x)/d-xn ]


## <a name="riskdecomp"></a>Risk Decomposition using Euler's Theorem

Let RMp(x) denote a portfolio risk measure that is a homogeneous function
of degree one in the portfolio weight vector x.

    RMp(x) = sigma_p(x) 
           
           = sqrt(x' * Sigma * x)

Euler's Theorem gives risk decomposition:

                   d-RMp(x)        d-RMp(x)
    RMp(x) = x1 * ---------- + x2 ---------- + ...
                    d-x1            d-x2

                   d-RMp(x)
           = x' * ----------
                    d-x

                 d-RMp(x)
    MCR_i(RM) = ----------      : Marginal Contribution to Risk (MCR) from asset i
                  d-x_i

                  delta-RMp(x)
              ~= -------------  : change in Risk per change in asset i
                   delta-x_i      (holding all other assets constant - 
                                   which isn't possible for portfolio weights)
                 

    CR_i(RM) = x_i * MCR_i(RM)  : contribution of risk (CR) = weighted MCR


    RMp(x) = CR_1(RM) + CR_2(RM) + ... + CR_n(RM)


                 CR_i(RM)
    PCR_i(RM) = ---------
                  RMp(x)

### Risk Decomposition of Portfolio SD


    RMp(x) = sigma_p(x) 
    
           = sqrt(x' * Sigma * x)

    Euler decomposition:

                           d-sigma_p(x)
        sigma_p(x) = x' * ----------
                             d-x


        d-sigma_p(x)    d(x' * Sigma * x)^1/2  
       ------------- = ----------------------- 
          d-x                   d-x            

                        1
                     = --- (x' * Sigma * x)^-1/2 * 2 * Sigma * x
                        2

                            Sigma * x
                     = ----------------------
                        (x' * Sigma * x)^1/2

                            Sigma * x
                     = ----------------------
                            sigma_p(x)

    The partial derivative of sigma_p(x) wrt x_i is the ith row of the above.
   
    R: PerformanceAnalytics :: StdDev performs SD decomposition
    R: PerformanceAnalytics :: Var    performs Value-at-Risk decomposition


### <a name="mcr"></a>Marginal Contribution to Risk (MCR)


              delta-sigma_p
    MCR_i ~= ---------------     : "Change in sigma_p per change in asset i weight"
                delta-x_i

    delta-sigma_p ~= MCR_i * delta-x_i

* However, changing one asset weight 
* requires changing one or more other asset weights

Example: changing one other asset weight:

    delta-x_i = - delta-x_j

    delta-sigma_p ~= MCR_i * delta-x_i + MCR_j * delta-x_j

                  ~= (MCR_i - MCR_j) * delta-x_i


## <a name="beta"></a>Beta

              cov(R_i, R_p(x))
    Beta_i = ------------------
                 var(R_p(x))

    Beta_i = Regression coefficient when regressing R_i vs R_p(x)

           : Measure of asset contribution to portfolio risk


    MCR_i = Beta_i * sigma_p(x)

    CR_i = x_i * Beta_i * sigma_p(x)

    PCR_i = x_i * Beta_i


#### When `Beta_i > 1`

    MCR_i > sigma_p(x)

    CR_i > x_i * sigma_p(x)

    PCR_i > x_i 

    Increasing allocation will increase volatility of portfolio


#### When `Beta_i < 1`

    MCR_i < sigma_p(x)

    CR_i < x_i * sigma_p(x)

    PCR_i < x_i 

    Increasing allocation will decrease volatility of portfolio

### Beta as a Measure of Portfolio Risk

* Asset specific risk can be diversified away via portfolios
* what remains is "portfolio risk"
    * measured by Beta
* riskiness of asset should be judged in a portfolio context
* beta measures the portfolio risk of an asset

### Beta and Risk Return Tradeoff

* If Beta is an appropriate measure of (undiversifiable) asset risk,
* Then the asset's expected return should depend on Beta
    * CAPM!



## <a name="simodel"></a>Sharpe Single Index Model

* precursor to CAPM
* extension of CER Model
    * to allow for common source of risk for all assets
    * CER Model is "special case" of Single Index Model
        * where `Beta_i` = 0 for all assets
        * `alpha_i` = `mu_i`

.

    R_it = alpha_i + Beta_i * R_Mt + err_it

    i = 1, ... N assets
    t = 1, ... T time

    alpha_i asset return when R_Mt = 0, constant over time
    Beta_i: cov(R_it, R_Mt)/var(R_Mt),  constant over time
          : contribution of asset i to market risk
          : magnitude of asset move in response to market move

    R_Mt: return on diversified market index portfolio
    err_it: random error term unrelated to R_Mt



#### Assumptions:

    cov(R_Mt, err_is) = 0,      for all t,s

    cov(err_is, err_js) = 0,    for all i,j, t,s

    err_it ~ iid N(0, sigma_e,i^2)

    R_M,t ~ iid N(mu_M, sigma_M^2)

* error terms are unrelated to the market-wide return
* error terms are firm-specific
    * i.e. uncorrelated
    * different than CER Model
* error terms are normally distributed
* market returns are normally distributed


#### Interpretation of error term, `err_it`

    err_it = R_it - alpha_i - Beta_i * R_Mt

* Return on market index `R_Mt` captures common "market-wide" news
* `Beta_i` measures sensitivity of asset i to "market-wide" news
* Random error term `err_it` captures "firm-specific" news
    * unrelated to market-wide news
* Returns are correlated only through their common exposure to..
    * common market-wide news
    * captured by `Beta_i`
    * i.e. Returns otherwise are NOT correlated


## <a name="simodelstats"></a>Statistical Properties of SI Model (Unconditional)

    mu_i = E[R_it] 
    
         = alpha_i + Beta_i * mu_M

           (firm-specific) + (market movements)


    sigma_i^2 = var(R_it) 
    
              = Beta_i^2 * sigma_M^2 + sigma_e,i^2

                   (market news)    +   (firm-specific news)

    sigma_ij = cov(R_it, R_jt) 
    
             = sigma_M^2 * Beta_i * Beta_j

    R_it ~ N(mu_i, sigma_i^2) 
    
         ~ N(alpha_i + Beta_i * mu_M, Beta_i^2 * sigma_M^2 + sigma_e,i^2


Covariances:

    sigma_ij = 0: Beta_i or Beta_j is 0
    sigma_ij > 0: Beta_i, Beta_j have the same sign
    sigma_ij < 0: Beta_i, Beta_j have different signs



## <a name="simodelstatscond"></a>Statistical Properties of SI Model (Conditional on `R_Mt = r_Mt`)

    E[R_it|R_Mt=r_Mt] = alpha_i * Beta_i * r_Mt

    sigma_i^2 = var(R_it|R_Mt=r_Mt) 
    
              = sigma_e,i^2

    sigma_ij = 0

    R_it|R_Mt=r_Mt ~ N(alpha_i + Beta_i * r_Mt, sigma_e,i^2)


## <a name="decompvar"></a>Decomposition of Total Variance


    sigma_i^2 = var(R_it) 
    
              = Beta_i^2 * sigma_M^2 + sigma_e,i^2

                 (market-variance)   + (non-market-variance)


        Beta_i^2 * sigma_M^2   sigma_e,i^2
    1 = -------------------- + -------------
           sigma_i^2            sigma_i^2

      =     R_i^2 + 1 - R_i^2


    where,

             Beta_i^2 * sigma_M^2 
    R_i^2 =  --------------------  = proportion of market variance for asset i
                sigma_i^2         

    1 - R_i^2 = proportion of non-market-variance for asset i

    similar to R^2 in regression analysis


#### Sharpe's Rule of Thumb

* A typical stock has `R_i^2 = 0.30`
    * i.e. 30% of asset's total variance 
    * is due to/explained by market variance
    * THAT's NOT MUCH!


### SI Model Return Covariance Matrix
    

            [ sigma_1^2     sigma_12    sigma_13  ]
    Sigma = [ sigma_12      sigma_2^2   sigma_23  ]   
            [ sigma_13      sigma_23    sigma_3^2 ]

            [ Beta_1^2*sigma_M^2 + sigma_e,1^2      sigma_M^2*Beta_1*Beta_2             sigma_M^2*Beta_1*Beta_3          ]
          = [ sigma_M^2*Beta_1*Beta_2               Beta_2^2*sigma_M^2 + sigma_e,2^2    sigma_M^2*Beta_2*Beta_3          ]
            [ sigma_M^2*Beta_1*Beta_3               sigma_M^2*Beta_1*Beta_3             Beta_3^2*sigma_M^2 + sigma_e,3^2 ]


                        [ Beta_1^2      Beta_1*Beta_2   Beta_1*Beta_3 ]   [ sigma_e,1^2     0           0        ]
          = sigma_M^2 * [ Beta_1*Beta_2 Beta_2^2        Beta_2*Beta_3 ] + [    0         sigma_e,2^2    0        ]
                        [ Beta_1*Beta_3 Beta_2*Beta_3   Beta_3^2      ]   [    0            0        sigma_e,3^2 ]

            [ Beta_1 ]
    Let B = [ Beta_2 ]
            [ Beta_3 ]

            [ sigma_e,1^2     0           0        ]
    Let D = [    0         sigma_e,2^2    0        ]
            [    0            0        sigma_e,3^2 ]            

    Sigma = sigma_M^2 * Beta * Beta'     +     D

             (covariance due to market)  + (asset-specific variances)


## <a name="simodelport"></a>SI Model and Portfolios

    R_1t = alpha_1 + Beta_1 * R_Mt + err_1t

    R_2t = alpha_2 + Beta_2 * R_Mt + err_2t

    x1 + x2 = 1


#### Portfolio Return

alphas and betas and error terms are additive:

    R_p,t = x1 * R_1t + x2 * R_2t

          = alpha_p + Beta_p * R_Mt + err_p,t


    alpha_p = x1 * alpha_1 + x2 * alpha_2
    Beta_p = x1 * Beta_1 + x2 * Beta_2
    err_pt = x1 * err_1t + x2 * err_2t


### SI Model with Large Portfolios

For large N, the error terms cancel out:

    R_p,t = alpha' + Beta' * R_Mt + err_t

          = alpha' + Beta' * R_Mt 

    var(R_p,t) = Beta'^2 * var(R_Mt)

    R_p^2 ~= 1

Implications:

* all non-market variance (error term) is diversified away
* magnitude of portfolio variance 
    * is porportional to market variance
    * determined by portfolio Beta
* Approximately 100% of portfolio variance
    * is due to market variance
* Sharpe and others argue that Beta risk should drive asset risk premiums
    * i.e. asset expected return
    * because asset-specific risk can be diversified away
        * in a large portfolio
        * leaving only market risk
        * and the asset's sensitivity to market risk: BETA!
    * therefore only BETA should influence asset risk premium/expected return
        * high BETA should have high returns
        * low BETA should have low returns
    * basis for CAPM
    * some empirical support for this model


## <a name="estimatesimodel"></a>Estimating the SI Model


    R_it = alpha_i + Beta_i * R_Mt + err_it

    i = 1, ... N assets
    t = 1, ... T time

    err_it ~ iid N(0, sigma_e,i^2)

    R_M,t ~ iid N(mu_M, sigma_M^2)


    E[R_it] = mu_i      

            = alpha_i + Beta_i * mu_M


    Var(R_it) = sigma_i^2 
    
              = Beta_i^2 * sigma_M^2 + sigma_e,i^2


    alpha_i = mu_i - Beta_i * mu_M

              cov(R_it, R_Mt)
    Beta_i = ------------------
                 var(R_Mt)

Main parameters to estimate:

    alpha_i

    Beta_i

    sigma_e,i^2

Three ways to estimate parameters:

1. **Plug in Principle**
    * estimate model parameters using sample statistics
2. **Linear Least Squares**
    * fit a regression model using Least Squares
    * turns out to be exactly the same as plug-in principle in this case
3. **Maximum Likelihood Function**
    * turns out to be exactly the same as plug-in principle in this case


#### Plug-in Principle:

        ^alpha_i = ^mu_i - ^Beta_i * ^mu_M

        ^err_it = R_it - ^alpha_i - ^Beta_i * ^mu_M

                        1
        sigma_e,i^2 = ----- SUM_t=1..T ^err_t^2
                       T-2

                  T-2: two degrees of freedom
                       two estimated parameters: ^alpha_i, ^Beta_i


## <a name="lse"></a>Least Squares Estimates

The SI model is in the same format as a linear regression model 

    R_it = alpha_i + Beta_i * R_Mt + err_it

R:

    msft.fit = lm(formula = msft ~ sp500, data = si.df)

    confint(msft.fit, level=0.95)
    coef(msft.fit)
    residuals(msft.fit)
    fitted(msft.fit)        # the regression line


### Statistical Properties of Least Squares Estimates


The SE estimates are reported in the `lm` output.


                          ^sigma_e,i          
    ^SE(^alpha_i) = ------------------------  * sqrt( 1/T * SUM_t=1..T R_Mt^2 )
                      sqrt(T * ^sigma_M^2)

                           ^sigma_e,i      
    ^SE(^Beta_i)  =  ----------------------
                       sqrt(T * ^sigma_M^2)


For large enough T, CLT tells us:

    ^alpha_i ~ N(alpha_i, ^SE(^alpha_i)^2)

    ^Beta_i ~ N(Beta_i, ^SE(^Beta_i)^2)


Approx 95% conf intervals:

    ^alpha_i +/- 2 * ^SE(^alpha_i)

    ^Beta_i +/- 2 * ^SE(^Beta_i)



* SE goes down with smaller error terms `^sigma_e,i`
* SE for Beta is smaller with larger `^sigma_M^2`
    * greater market variance makes it easier to estimate Beta
    * where Beta is the asset's sensitivity to market variance
* SE goes to 0 as T goes to INF
    * alpha and Beta are **consistent estimators**
* SE for the following can be computed via bootstrap (no analytical solution):
    * `sigma_e,i^2`
    * `sigma_e,i`
    * `R-squared`


## <a name="simodelmatrix"></a>SI Model Using Matrix Algebra


    R_it = alpha_i + Beta_i * R_Mt + err_it

    Observations: t = 1,...,T


    [ R_i1 ]              [ 1 ]             [ R_M1 ]   [ err_i1 ]
    [ R_i2 ]              [ 1 ]             [ R_M2 ]   [ err_i2 ]
    [ ...  ]  = alpha_i * [ 1 ] + Beta_i *  [ ...  ] + [ ...    ] 
    [ R_iT ]              [ 1 ]             [ R_MT ]   [ err_iT ]


      R_i     = alpha_i *   1   + Beta_i *    R_M    +    err_i

                             [ alpha_i ]
              = [ 1  R_M ] * [ Beta_i  ]  + err_i


              = X * gamma_i + err _ 


           X = [ 1  R_M ]


           gamma_i = [ alpha_i ]
                     [ Beta_i  ]


First order conditions (set partial derviatives wrt alpha and Beta):

    [ 1'*R_i   ]     [ 1'*1     1'*R_M   ]   [ ^alpha_i ]
    [ R_M'*R_i ]  =  [ 1'*R_M   R_M'*R_M ] * [ ^Beta_i  ]

        X'*R_i    =  X' * X * ^gamma_i

        ^gamma_i  = (X'*X)^-1 * X' * R_i

        ^gamma_i: the least squares estimates of alpha_i and Beta_i:


## <a name="simodel4port"></a>SI Model for 4 Asset Portfolio

    # returns: equally weighted portfolio
    port = (si.df$sbux + 
            si.df$msft + 
            si.df$nord + 
            si.df$boeing) / 4

    new.data = data.frame(si.df, port)
    port.fit = lm(port ~ sp500, data=new.data)
    summary(port.fit)

    beta.sbux = coef(lm(sbux~sp500, data=si.df))[2]
    beta.msft = coef(lm(msft~sp500, data=si.df))[2]
    beta.nord = coef(lm(nord~sp500, data=si.df))[2]
    beta.boeing = coef(lm(boeing~sp500, data=si.df))[2]

    (beta.sbux + beta.msft + beta.nord + beta.boeing)/4

    coef(port.fit)[2]


* Portfolio Beta should be weighted avg of asset betas
* Portfolio Beta tends to be close to one
* Portfolio "fits" data better than individual assets
    * Portfolio approaches market returns as #assets increase
    * R-squared is higher
    * Beta estimated more preceisly
        * lower SE(Beta)
    * residuals are lower
        * lower `sigma_e`
        * diversification effect


## <a name="simodelcovmat"></a>Estimating SI Model Covariance Matrix

* In the CER Model, we estimated covariances...
    * by looking at sample covariances
* In the SI Model, covariances are estimated...
    * from their Betas
    * i.e. their common exposure to the market

Recall, in the SI Model:

    sigma_ij = cov(R_it, R_jt) 
            
             = sigma_M^2 * Beta_i * Beta_j

            [ sigma_1^2     sigma_12    sigma_13  ]
    Sigma = [ sigma_12      sigma_2^2   sigma_23  ]   
            [ sigma_13      sigma_23    sigma_3^2 ]

                        [ Beta_1^2      Beta_1*Beta_2   Beta_1*Beta_3 ]   [ sigma_e,1^2     0           0        ]
          = sigma_M^2 * [ Beta_1*Beta_2 Beta_2^2        Beta_2*Beta_3 ] + [    0         sigma_e,2^2    0        ]
                        [ Beta_1*Beta_3 Beta_2*Beta_3   Beta_3^2      ]   [    0            0        sigma_e,3^2 ]

            [ Beta_1 ]
    Let B = [ Beta_2 ]
            [ Beta_3 ]

            [ sigma_e,1^2     0           0        ]
    Let D = [    0         sigma_e,2^2    0        ]
            [    0            0        sigma_e,3^2 ]            

    Sigma = sigma_M^2 * Beta * Beta'     +     D

             (covariance due to market)  + (asset-specific variances)


Estimate Sigma using plug-in principle: 

    ^Sigma = ^sigma_M^2 * ^Beta * ^Beta' + ^D

R:  
    beta.vec <- c( coef(lm(sbux~sp500, data=si.df))[2],
                   coef(lm(msft~sp500, data=si.df))[2],
                   coef(lm(nord~sp500, data=si.df))[2],
                   coef(lm(boeing~sp500, data=si.df))[2] )

    D.mat <- diag(c( sig2e.sbux,
                     sig2e.msft,
                     sig2e.nord,
                     sig2e.boeing ) )

    cov.market <- sig2.sp500 * beta.vec %*% t(beta.vec)


    cov.si <- cov.market + D.mat

    # compare with sample covariances
    print(cov.hat, digits=4)
    print(cov.si, digits=4)


## <a name="simodelhyptest"></a>Hypothesis Testing in SI Model

* Test:
    * returns are normally distributed
    * error terms are normally distributed
    * params are constant over time
    * beta is signficantly different from 0 or specified value

Hypothesis: Beta is 0:

              ^Beta_i - 0
    t.stat = -------------
              ^SE(^Beta_i)



