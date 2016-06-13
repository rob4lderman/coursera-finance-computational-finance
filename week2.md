
# Coursera: Intro to Computational Finance and Financial Econometrics: Week 2: Probability Review


## Univariate Random Variables

* A random variable X can take on a given set of values
    * "set of values" = Sample Space, S_x
* where the likelihood of of X being any given value
    * in the sample space, S_x
    * is determined by the random variable's probability distribution function
* Future Price of a stock (e.g. 1 month from now)...
    * is a random variable
    * we're trying to predict what the future price will be
* Example, MSFT:
    * X = price of MSFT in one month
        * S_x = { R: 0 < X <= M }
        * M is some large number
    * X = simple return on a one month investment
        * S_x = { R: -1 <= X < M }
        * M is some large number
    * X = 1 if price goes up; X = 0 if price goes down
        * S_x = {0,1}

### Discrete Random Variables

* can take on a finite number of n different values
    * x1, x2, x3, ... xn
    * X (capital X): the random variable 
    * x1, x2 (little x): specific values of the random variable
* **pdf** is function such that p(x) = P(X = x)

PDF must satisfy:

1. p(x) >= 0 for all x E S_x; 
    * p(x) = 0 for all x !E S_x
2. SUM_x E S_x [ p(x) ] = 1
3. p(x) <= 1 for all x E S_x


### Bernouli Distribution

* Consider two mutually exclusive events
    * "success" (X=1) and "failure" (X=0)
* P(X=1) = pi, where 0 < pi < 1
* P(X=0) = 1 - pi
* pdf = p(x) = pi^x * (1 - pi)^(1-x)
    * x = {0,1}


### Continuous Random Variable

* can take on any real value
* pdf = f(x) such that
    * Pr(X E A) = integral_A [ f(x) * dx ]
    * Pr(X E A) = "probability X is in the range A"
    * Pr(X E A) = "area under the probability curve over interval A"

PDF must satisfy:

1. f(x) >= 0
2. integral_INF [ f(x) * dx ] = 1


### Uniform Distribution

* X ~ U[a,b]
    * "~" = "is distributed as"
    * X is distributed uniformly over the interval a,b

.

                   1/(b-a), for a <= x <= b
    PDF: f(x) = {  0,       otherwise


    f(x) >= 0, provided b > a



## Cumulative Distribution Function (CDF)

CDF = F(x) = Pr(X <= x)

Features of the CDF:

* F(-INF) = 0
* F(INF) = 1
* Pr(X > x) = 1 - F(x)
* Pr(x1 < X <= x2) = F(x2) - F(x1)
* the derivative of the CDF is the PDF
    * if X is continuous
* Continuously increasing function
    * the derivative is never negative
* CDF is a plot of probability (y-axis) vs S_x (x-axis)
    * y-axis is P(X <= x), where x E S_x

Note: For a continuous RV:

* Pr(X <= x) = Pr(X < x)
* Pr(X = x) = 0



## Quantiles of a Distrubution

* X is a rv with continuous CDF F(x)
* The alpha * 100% quantile of F(x)...
    * is the value `q_alpha`, such that
    * `F(q_alpha) = alpha = Pr(X <= q_alpha)`
    * `q_alpha` is the quantile
    * `q_alpha` E S_x
* i.e. the area under the PDF to the left of `q_alpha`
* if F(x) has an inverse, F^-1(x), then:
    * `q_alpha = F^-1(alpha)`
    * F^-1(x) aka the "quantile function"
* alhpa = 0.5 = 50% quantile = "median"


Quantiles are often used in financial analysis to tell us the probability of loss

* X represents rate of return
* with 1% probability, how much could we lose?
    * ANSWER = 0.01 quantile
    * which is the rate of return with 1% probability of occurring
        * technically, 1% probability that X <= 0.01 quantile rate of return



## Standard Normal Distribution

    
                       1                  1
    f(x) = phi(x) = ----------- * exp( - --- * x^2 )
                     sqrt(2*pi)           2


    Pr( -1 <= x <= 1 ) ~= 0.67
    Pr( -2 <= x <= 2 ) ~= 0.95
    Pr( -3 <= x <= 3 ) ~= 0.99


* In finance, monthly returns tend to approximate a normal distribution
* but not all the time...
    * so we'll be looking at other distributions as well for modeling returns

R functions:

* pnorm(x): Pr(X <= x)
* qnorm(alpha): F^-1(alpha)
* dnorm(x): PDF(x) 


## Expected Value, Variance, Skewness, Kurtosis

* Expected Value: Center of mass
    * probability-weighted avg 
* Variance and Standard Deviation: dispersion about mean
    * dispersion of returns = risk
* Skewness: symmetry about mean
    * some assets may be skewed toward losses (or gains)
    * `Skew(x) > 0` = positive skewness
        * more values above the mean
        * long right tail
    * `Skew(x) < 0` = negative skewness
        * more values below the mean
        * long left tail
* Kurtosis: tail thickness
    * extreme gains/losses
    * Normal distribution implies that extreme events are unlikely
        * but we know from experience that they happen
        * more often than predicted by normal distribution
    * large kurtosis = more extreme values
    * small kurtosis = fewer extreme values
    * Kurt(normal) = 3
        * benchmark value
        * Kurt(X) > 3: fat tails
        * Kurt(X) < 3: thin tails
        * "excess kurtosis": kurt more or less than 3

.

    Discrete: E[X] = mu_x = SUM_Sx [ x * p(x) ]

    Continuous: E[X] = integral_Sx [ x * pdf(x) ]


    Var(X) =  sigma_x^2 = E[ (X - E[X])^2 ]

                        = E[X^2] = E[X]^2



                    X - mu_x
    Skew(X) = E[ ( ---------- )^3 ]
                    sigma_x


                        X - mu_x
    Kurtosis(X) = E[ ( ---------- )^4 ]
                        sigma_x


    


## General Normal Distribution

    
                        
                       1                            1        x - mu_x
    f(x) = phi(x) = ----------------------* exp( - --- * ( ------------ )^2 )
                     sqrt(2*pi*sigma_x^2)           2        sigma_x



    Quantile:
    
    q_alpha = mu_x + sigma_x * PHI^-1(alpha) = mu_x + sigma_x * z_alpha

    where PHI^-1(alpha) = "quantile" function for the standard normal



* The normal distribution may not be appropriate for simple returns
    * allows for impossible values (e.g. R < -1)
* More appropriate for cc returns
    * `r = ln(1 + R)`
    * `R = e^rt - 1`
    * r can take on values < -1
        * R still: -1 <= R <= 1
    * in fact can take on any value between -INF and +INF


## The Log-Normal Distribution

    X ~ N(mu_x, sigma_x^2), -INF < X < INF

    Y = exp(X) ~ lognormal(mu_x, sigma_x^2), 0 < Y < INF

    E[Y] = exp(mu_x + sigma_x^2/2)

    Var(Y) = exp(2*mu_x + sigma_x^2)( exp(sigma_x^2) - 1 )


Example: continuously compounded returns

    r ~ N(0.05, 0.50^2)

    1 + R ~ lognormal(0.05, 0.50^2)

R:

* rlnorm
* plnorm
* qlnorm
* dlnorm


## Student's t distribution

* similar to normal
* extra parameter: degrees of freedom
    * controls the thickness of the tails
    * smaller degrees of freedom = thicker tails
    * as degrees of freedom goes to INF
        * t-distribution approaches standard normal
    * if df close to 4...
        * kurtosis is large and the tails are thick
    * if df < 4...
        * kurtosis is INF
* R:
    * dt
    * pt
    * rt
    * qt



.

    df = degrees of freedom

    X ~ t_df

                gamma(df + 1)                     x^2         df+1
    PDF(x) =  --------------------------- * (1 + ----- )^( - ------ )
               sqrt(df*pi) * gamma(df/2)          df           2

              -INF < x < INF, df > 0


    gamma(z) = integral_0-INF [ t^(z-1) * e^-t * dt ]


    E[X] = 0,                   df > 1

    Var(X) =  df / (df-2),      df > 2

    skew(X) = 0,                df > 3

    kurt(X) = 6 / (df-4) - 3,   df > 4


## Linear Functions of Random Variables


* Let X be a random variable
* `Y = a * X + b`

.

    E[Y] = E[X] * a + b

    Var(Y) = a^2 * Var(X)

    sigma_y = a * sigma_x


* Let `X ~ N(mu_x, sigma_x^2)`
* and `Y = aX + b`
* then `Y ~ N(mu_y, sigma_y^2)`



## Value at Risk

* Value-at-risk: how much money can you lose, with a certain probability?
* Example: a $10K investement in MSFT, for 1 month
    * R = simple monthly return
    * `R ~ N(0.05, 0.10^2)`
* Calculate how much we can lost with a specified probability, alpha


Questions to answer:

1. What's the probability distribution of end-of-month wealth, W1?
    * W1 = $10K ( 1 + R )
    * Linear function of R...
    * So W1 is also normally-distributed random variable
2. What is Pr(W1 < 9000)?
3. What value of R produces W1 = 9000?
4. What's the monthly value-at-risk (VaR), with 5% probability?
    * i.e. how much can we lose if `R <= q_.05`
    * `VaR_alpha = $W0 * q_alpha(R)`
    * Note: people often assume R is normally distributed
    * however, this is not often the case in practice
    * this leads to UNDERestimating the actual VaR
    * using a t-distribution yields a more accurate estimate


1\. What's the probability distribution of end-of-month wealth, W1?
W1 is a Linear Function of R.  R is normally distributed.  Therefore
W1 is also normally distributed.

    
    E[W1] = $10K * (1 + E[R])
          = 10,500

    Var(W1) = $10K^2 * Var(R)
            = 1,000,000

    W1 ~ N(10500, 1000^2)
            
2\. What is Pr(W1 < 9000)?

    pnorm(9000, mean=105000, sd=1000)
    # [1] = 0.067


3\. What value of R produces W1 = 9000?

    FV = PV * (1 + R)
    1 + R = PV/FV
    R = PV/FV - 1
    R = -0.10

Notice that -0.10 is the 6.7% quantile of the distribution of R:

    q_.067 = Pr(R < -0.10) = 0.067



4\. What's the monthly value-at-risk (VaR), with 5% probability?

    # q_.05(R) = Pr(R < q_.05(R)) 
    qnorm(0.05, mean=0.05, sd=0.10)
    # [1] -0.1144854

    FV = PV * (1 + R),  R = q_.05(R) = -0.1144

    VaR = PV - PV * (1 + R)
    VaR = PV * R

    5% VaR = 10000 * -0.1144 
           = $1144
    
    so we could lose $1144 OR MORE, with 5% probability



## Value-at-Risk for continuously compounded returns

    r = ln(1 + R)
    R = e^r - 1

    r ~ N(mu_r, sigma_r^2)
    W0 = initial investment

1. compute alpha quantile or normal distribution
    * `q_alpha(r) = mu_r + sigma_r * z_alpha`
2. convert alpha quantile for r into alpha quantile for R
    * `q_alpha(R) = e^q_alpha(r) - 1`
3. compute VaR_alpha
    * `VaR_alpha = $W0 * q_alpha(R)`

.

    W0 = 10000
    r ~ N(0.05, 0.10^2)
    
    alpha = 0.05
    1. q_.05(r) = mu_r + sigma_r * z_.05
                = -0.114

    2. q_.05(R) = e^q_.05(r) - 1
                = -0.108

    3. VaR_.05 = W0 * q_.05(R)
               = 1080


