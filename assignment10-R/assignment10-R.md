
# Coursera: Intro to Computational Finance: Assignment 10: Sharpe Single Index Model

Note: the course doesn't actually have an assignment 10.  Just playing around with the SI Model.


```r
    library(tseries)    # get.hist.quote
    library(zoo)        # coredata
    library(quadprog)   # solve.QP

    #
    # Load price data from yahoo
    #
    SBUX_prices <- get.hist.quote(instrument="sbux", 
                                  start="2001-01-01",
                                  end="2015-12-31", 
                                  quote="AdjClose",
                                  provider="yahoo", 
                                  origin="1970-01-01",
                                  compression="m", 
                                  retclass="zoo", 
                                  quiet = TRUE)
    MSFT_prices <- get.hist.quote(instrument="msft", 
                                  start="2001-01-01",
                                  end="2015-12-31", 
                                  quote="AdjClose",
                                  provider="yahoo", 
                                  origin="1970-01-01",
                                  compression="m", 
                                  retclass="zoo", 
                                  quiet = TRUE)
    IBM_prices <-  get.hist.quote(instrument="ibm", 
                                  start="2001-01-01",
                                  end="2015-12-31", 
                                  quote="AdjClose",
                                  provider="yahoo", 
                                  origin="1970-01-01",
                                  compression="m", 
                                  retclass="zoo", 
                                  quiet = TRUE)
    SP500_prices <- get.hist.quote(instrument="^gspc", 
                                   start="2001-01-01",
                                   end="2015-12-31", 
                                   quote="AdjClose",
                                   provider="yahoo", 
                                   origin="1970-01-01",
                                   compression="m", 
                                   retclass="zoo", 
                                   quiet = TRUE)

    #
    # Compute simple returns, means, sd, cov
    # Portfolio theory assumes simple returns (as opposed to cc returns)
    # 
    all_prices <- merge(SP500_prices, IBM_prices, MSFT_prices, SBUX_prices)

    # diff: computes pt1 - pt0
    # lag: shifts all_prices by k=-1 (so that pt-1 -> pt)
    simple_returns <- diff(all_prices) / lag(all_prices,k=-1)
    simple_returns.mat <- coredata(simple_returns)

    asset_names <- c("SP500", "IBM", "MSFT", "SBUX")
    colnames(simple_returns.mat) <- asset_names

    dim(simple_returns.mat)
```

```
## [1] 179   4
```

```r
    rbind(head(simple_returns.mat),
          tail(simple_returns.mat))
```

```
##                SP500          IBM        MSFT         SBUX
##        -0.0922907358 -0.107019201 -0.03377691 -0.046309019
##        -0.0642047105 -0.037237183 -0.07309321 -0.108922854
##         0.0768143618  0.197130342  0.23885717 -0.088071066
##         0.0050901871 -0.027833836  0.02110698  0.008785605
##        -0.0250353891  0.015205715  0.05521826  0.178278600
##        -0.0107401501 -0.073039677 -0.09328760 -0.215652054
## [174,]  0.0197420297 -0.004119073  0.05775761  0.080380521
## [175,] -0.0625808182 -0.079463536 -0.06194892 -0.052980449
## [176,] -0.0264428316 -0.019744405  0.01700364  0.038932604
## [177,]  0.0829831178 -0.033731148  0.18933580  0.100809283
## [178,]  0.0005048693  0.004624335  0.03944405 -0.015682161
## [179,] -0.0175301852 -0.012910646  0.02079115 -0.022153421
```

```r
    mu.vec <- apply(simple_returns.mat, 2, mean)
    sigma.vec <- apply(simple_returns.mat, 2, sd)

    rbind(mu.vec, sigma.vec)
```

```
##                 SP500         IBM        MSFT       SBUX
## mu.vec    0.003205399 0.004668487 0.007989258 0.01686569
## sigma.vec 0.043358161 0.067195801 0.073275423 0.08565586
```

### 1. Compute SI Model estimates (alpha,beta,`sigma_e,i^2`) for each security from sample statistics, using SP500 as the market index.

#### Review: SI Model:

    R_it = alpha_i + beta_i * R_Mt + err_it

    i = 1..N assets
    t = 1..T time

    return-for-asset-i-at-time-t = alpha-for-asset-i + 
                                   beta-for-asset-i * market-return-at-time-t + 
                                   error-term-for-asset-i-at-time-t

               cov(R_it, R_Mt)
    ^beta_i = ---------------
                var(R_Mt)

    ^alpha_i = E[R_i] - beta_i * E[R_Mt]
        
             = ^mu_i - beta_i * mu_M

    ^err_it = R_it - ^alpha_i - ^Beta_i * ^mu_M

                    1
    sigma_e,i^2 = ----- SUM_t=1..T ^err_t^2
                   T-2


#### R:


```r
    beta.vec <- apply(simple_returns.mat, 2, function(R_i) { cov(R_i, simple_returns.mat[,"SP500"]) / sigma.vec["SP500"]^2 })   

    alpha.vec <- mu.vec - beta.vec * mu.vec["SP500"]

    simple_returns.df <- as.data.frame(simple_returns.mat)
    errors.mat <- mapply(function(R_i,alpha_i,beta_i) { R_i - alpha_i - beta_i * simple_returns.df$SP500  }, 
                         simple_returns.df, 
                         alpha.vec, 
                         beta.vec)

    # above mapply is same as:
    # errors.mat2 <- matrix(0, 
    #                       nrow=nrow( simple_returns.mat ),
    #                       ncol=ncol( simple_returns.mat ))

    # colnames(errors.mat2) <- asset_names

    # for (i in seq_along(asset_names)) {
    #     errors.mat2[,i] = simple_returns.mat[,i] - alpha.vec[i] - beta.vec[i] * simple_returns.mat[,"SP500"]
    # }

    sigma_e.vec <- apply(errors.mat, 2, sd)

    rbind(alpha.vec,
          beta.vec,
          sigma_e.vec)
```

```
##                    SP500        IBM        MSFT       SBUX
## alpha.vec   4.336809e-19 0.00159085 0.004532139 0.01382926
## beta.vec    1.000000e+00 0.96014148 1.078530043 0.94728661
## sigma_e.vec 7.109311e-18 0.05274677 0.056413669 0.07516627
```

### 1. Compute SI Model estimates (alpha,beta,`sigma_e,i^2`) for each security from linear regression model.



```r
    linear.models <- mapply( function(R_i) { lm(R_i ~ simple_returns.df$SP500) },
                             simple_returns.df,
                             SIMPLIFY=F)

    alpha.lm.vec <- sapply( linear.models, function(lm_i) { coef(lm_i)[1] } )
    names(alpha.lm.vec) <- asset_names

    beta.lm.vec <- sapply( linear.models, function(lm_i) { coef(lm_i)[2] } )
    names(beta.lm.vec) <- asset_names

    errors.lm.mat <- sapply( linear.models, function(lm_i) { residuals(lm_i) } )
    sigma_e.lm.vec <- apply(errors.lm.mat, 2, sd)

    rbind(alpha.lm.vec,
          beta.lm.vec,
          sigma_e.lm.vec)
```

```
##                       SP500        IBM        MSFT       SBUX
## alpha.lm.vec   1.037275e-18 0.00159085 0.004532139 0.01382926
## beta.lm.vec    1.000000e+00 0.96014148 1.078530043 0.94728661
## sigma_e.lm.vec 3.364181e-18 0.05274677 0.056413669 0.07516627
```


