
# Coursera: Intro to Computational Finance: Week 5: Descriptive Statistics

## Covariant Stationarity

* Looking at asset returns
* Asset returns are realization of covariant stationary time series
* Gaussian White Noise is the benchmark
* think critically about time-dependence
    * MA: moving average processes
    * AR: auto-regressive processes
    * both of which are covariant stationary
        * time-dependence only depends on lag
        * NOT on absolute time
* S&P500: top 500 US corps in term of market cap
    * S&P500 index: each corp weighted by market cap
    * requires constant re-balancing as prices change
    * TODO: is there a way to arbitrage this?
* What are the stylized facts about...
    * an asset's return
    * compared to the market return (S&P500)
* Looking at MSFT returns, 1998 - 2008
    * more volatility between 2000 - 2002
    * less volatility between 2003 - 2006
    * is it covariant stationary?
        * NO! 
        * the variance/volatility appears to change over time
    * mean appears constant over time
* MSFT and S&P500 returns appear to be correlated
    * beta!

**STYLIZED FACT:** Stocks with higher volatility tend to have higher expected return.
Compensation for the higher volatility.  But not always true.


## Descriptive Sample Statistics

* Typically we don't know the model/distribution...
    * of the underlying population data
    * that generates the random sample we observe
* Descriptive sample statistics
    * attempt to model the underlying distribution
    * to determine the unknown underlying PDF
    * and capture observed dependencies in the data
    
    
    Observed sample:  { X1 = x1, ... XT = xt } = {xt}t=1..T
    
    Observations from a stochastic (random) process.


### Comparing asset returns against Gaussian White Noise


    # MSFT: contains MSFT monthly returns
    set.seed(123)
    gwn = rnorm( length(MSFT), mean=mean(MSFT), sd=std(MSFT) )

    hist(MSFT, main="Histogram of MSFT Monthly CC Returns", col="stateblue1")

    hist(gwn, main="Histogram of Gaussian White Noise", col="stateblue1")

* Histogram: crude estimate of PDF
* Mean: center-of-balance of histogram
* Variance: Spread of data in histogram
* Skewness: Asymmetry of distribution
* Excess Kurtosis: thickness of tail data
* MSFT returns have "fatter tails" (excess kurtosis) than GWN
* S&P500 returns are negative-skewed
    * long left tail
* MSFT returns has higher spread (greater variance) than S&P500


#### Smoothed Histogram:

    # R: density()

    MSFT.density = density(MSFT)
    plot(MSFT.density, 
         type="l", 
         xlab="monthly return",
         ylab="probability density estimate",
         main="Smoothed histogram of MSFT cc monthly returns",
         col="orange",
         lwd=2)

    # overlay smoothed histogram
    hist(MSFT, 
         main="Histogram of MSFT Monthly CC Returns", 
         col="stateblue1",
         probability=T)     # converts data to density data

    points(MSFT.density,
           type="l",
           col="orange",
           lwd=2)


### Empirical Quantiles/Percentiles

* Sample quantiles are known as Empirical Quantiles
* Quartiles:
    * `q_.25`: first quartile
    * `q_.50`: second quartile (median)
    * `q_.75`: third quartile
    * `q_.75 - q_.25`: interquartile range (IQR)

.

    quantile(MSFT)
    quantile(MSFT, probs=c(0.01, 0.05))

    # compare to standard normal quantiles
    qnorm(p=c(0.01, 0.05), mean=mean(MSFT), sd=sd(MSFT))

Useful R functions for sample statistics:

    sort
    min
    max
    range
    quantile
    median
    IQR
    summary

### Historical Value-at-Risk

* uses empirical quantiles from data
    * as opposed to standard normal quantiles

. 

    VaR = W0 * q_alpha(R)               # simple returns

    VaR = W0 * (exp(q_alpha(R)) - 1)    # CC returns
    


### Sample Statistics

"Plug-in Principle" plug-in sample statistic for underlying population statistic

    Sample Mean = ^mu_x 
                = 1/T * SUM_t=1..T xt

    ^mu_x: the hat ^ indicates the value is an ESTIMATE from a sample

    Sample Variance = s_x^2 
                    = ^sigma_x^2 
                    = 1/(T-1) * SUM_t=1..T (xt - ^mu_x)^2

    Sample Standard Deviation = sqrt(s_x^2)

    Sample Skewness = ^skew 
                    = 1/(T-1) * SUM_t=1..T (xt - ^mu_x)^3 / s_x^3

    Sample Kurtosis = ^kurt 
                    = 1/(T-1) * SUM_t=1..T (xt - ^mu_x)^4 / s_x^4

    Sample Excess Kurtosis = ^kurt - 3


    # R functions:

    mean()      base pkg                sample mean
    colMeans()  base
    var()       stats                   sample variance
    sd()        stats                   sample standard deviation
    skewness()  PerformanceAnalytics    sample skewness
    kurtosis()  PerformanceAnalytics    sample EXCESS kurtosis (sample kurtosis - 3)

    # use apply() to apply functions over rows/cols of matrix or dataframe.
    # MSFT.SP500.mat: matrix contains two columns: returns for (1) MSFT and (2) SP500
    apply(MSFT.SP500.mat, 2, mean)      # 2 = apply to columns
    apply(MSFT.SP500.mat, 2, sd)      
    apply(MSFT.SP500.mat, 2, skewness)     
    apply(MSFT.SP500.mat, 2, kurtosis)     


## Empirical CDF


    Recall: CDF = F(x) = Pr(X <= x)

    Empirical CDF of a random sample:

    ^F(x) = 1/n * (#xi <= x)

             # of xi values <= x
          = ---------------------
                sample size

    R: ecdf

    # for gaussian white noise:
    n1 = length(gwn)
    plot(sort(gwn),
         (1:n1)/n1,
         type="s",
         ylim=c(0,1),
         ylab="#xi <= x")

    # TODO: plot empirical CDF vs normal CDF


## QQ Plot: Quantile-Quantile Plot

* QQ Plot
    * Compares empirical quantiles (y-axis) 
    * vs standard normal quantiles (x-axis)
    * actually, can be compared against any appropriate distribution

R:

    qqnorm()
    qqline()

    par(mfrow=c(2,2))
    qqnorm(gwn)
    qqline(gwn)
    qqnorm(MSFT)
    qqline(MSFT)
    qqnorm(SP500)
    qqline(SP500)
    par(mfrow=c(1,1))


## Outliers

* Outliers can really mess with sample statistics
    * except quantiles, which aren't affected as much
* "Robust Statistics": handle outliers well

.

    Moderate Outlier:

    IQR = q_.75 - q_.25: 50% of range

    ^q_.75 + 1.5 * IQR < x < ^q_.75 + 3 * IQR
    ^q_.25 - 3 * IQR   < x < ^q_.25 + 1.5 * IQR

    ^q_.75 + 1.5 * IQR = 75% of range higher than 75% quantile

    Extreme Outlier:

    x > ^q_.75 + 3 * IQR
    x < ^q_.25 - 3 * IQR

    ^q_.75 + 3 * IQR = 150% of range higher than 75% quantile


## Graphical Measures
    
### BoxPlots

* BoxPlots are robust to outliers
* uses median, quartiles to draw the box
* useful for comparing distributions/outliers of multiple distributions

.

    boxplot(MSFT, 
            outchar=T, 
            main="Boxplot of MSFT monthly CC returns",
            ylab="monthly cc return")


    boxplot(gwn,
            MSFT,
            SP500,
            names=c("gwn","MSFT","SP500"),
            outchar=T,
            main="Compare return distributions",
            ylab="monthly cc return")


### Four Graph Summary

1. histogram of monthly cc returns
2. boxplot of monthly cc returns
3. an empirical smoothed histogram
4. a QQ plot compared to normal distribution

.

    par(mfrow=c(2,2))

    hist(MSFT,
         probability=T,
         main="MSFT monthly cc returns",
         ylab="cc return")

    boxplot(MSFT,
            outchar=T,
            ylab="cc return")

    plot(MSFT.density,
         type="l",
         xlab="cc return",
         ylab="density estimate",
         main="Smoothed density")

    qqnorm(MSFT)
    qqline(MSFT)

    par(mfrow=c(1,1))


## BiVariate Sample Statistics

    Sample Covariance = s_xy 
                      = ^sigma_xy 
                      = 1/(T-1) SUM_t=1..T (xt - ^mu_x) * (yt - ^mu_y)

    var(cbind(MSFT,SP500,gwn))

    Sample Correlation = ^rho_xy 
                       = s_xy / s_x / s_y

    cor(cbind(MSFT,SP500,gwn))


### Scatterplot
    
    # plot SP500 vs MSFT returns
    plot(x=MSFT,
         y=SP500,
         main="Monthly CC returns on MSFT and SP500")

    # add lines for means
    abline(h=mean(SP500))
    abline(h=mean(MSFT))


    # pair-wise scatterplots between all datasets
    pairs(cbind(gwn,MSFT,SP500))


## Time-Series Descriptive Statistics

    Sample Autocovariance = ^gamma_j
                    
                          = 1/(T-1) SUM_t=j+1..T (xt - ^mu_x)(xt-j - ^mu_x)

                            j is time lag

    Sample Autocorrelation = ^rho_j
            
                           = ^gamma_j / ^sigma_x^2

    Sample Autocorrelation Function (SACF): Plot ^rho_j vs j

                    
## Summary of Stylized Facts for Monthly CC Returns

* Returns appear to be approximately normally distributed
    * some noticeable negative skewness
    * and excess kurtosis
* Many assets are contemporaneously correlated
* Returns are approximately UNcorrelated over time
    * no serial correlation


## Stylized Facts for DAILY Returns

* Daily returns are NOT normally distributed
    * empirical distributions have fatter tails
    * more outliers
* Daily returns are not auto-correlated
* Daily returns are NOT independent over time
* **Volatility Clustering**
    * volatility measures tend be auto-correlated
    * volatility = variance = squared returns
* Daily ABSOLUTE value of returns tend to be auto-correlated
    * similar to volatility/variance (squared returns)
* See Engle's GARCH model



