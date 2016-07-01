
# Coursera: Intro to Computational Finance and Financial Econometrics: Week 3: Probability Review Continued

* [Location Scale Model](#lsm)
* [BiVariate Probability Distribution](#bvpd)
* [Conditional probability](#cp)
* [BiVariate Continuous Distributions](#bvcd)
* [Covariance and Correlation](#cac)
* [General BiVariate Normal Distribution](#gbvnd)
* [Linear Combination of Two Random Variables](#lcotrv)
* [Example: Portfolio Analysis](#epa)
* [Example: Multi-period continously compounded returns](#empccr)


## <a name="lsm"></a>Location Scale Model

* For converting from standardized RV to non-standardized RV
* location = mu
* scale = sigma
* PRO: we can deduce quantiles of non-standard from standard

.

    X ~ N(mu, sigma^2)
    
    Z = (X - mu) / sigma ~ N(0,1)

    X = mu + sigma * Z


### Quantiles of Normal Distribution


    Pr(Z <= Z_alpha) = alpha

    Pr( mu + sigma * Z <=  mu + sigma * Z_alpha ) = alpha

    q(x)_alpha = mu + sigma * Z_alpha   

    Pr( X <= q(x)_alpha ) 


Example: Review **Value-at-Risk**: 


    VaR(alpha) = q(R)_alpha * PV

    # convert quantile from standard normal quantile
    q(R)_alpha = mu + sigma * Z_alpha 


## <a name="bvpd"></a>BiVariate Probability Distribution


* Joint probability distribution
    * Given RV's X,Y
    * `S_xy` = all possible combinations of X and y
    * p(x,y) = probability that X=x and Y=y at the same time
        * "joint probability"
* Marginal Distributions
    * From joint probability, can deduce individual probability distributions
    
.

    P(X = x1) = SUM_all-y p(x1,y)

    P(Y = y1) = SUM_all-x p(x,y1)


### <a name="cp"></a>Conditional probability

* If we know Y, what does this tell us about X?
* If X and Y are INDEPENDENT...
    * then knowing Y tells us NOTHING about X
* If X and Y are DEPENDENT...
    * then knowing Y tells us something about X


Example:

                       Pr(X=0, Y=0)
    Pr(X=0|Y=0 ) = -------------------
                         Pr(Y=0)

                                joint probabilty
    Conditional probability = --------------------
                               marginal probabilty


* Conditial probability distributions have their own mean and variance
* mean/variance of the RV given a certain condition
    * the "certain condition" must be on a DEPENDENT RV
    * otherwise the condition doesn't have any affect
* conditional variance is typically LESS than UN-conditional variance
    * because the condition typically...
    * removes some of the uncertainty 
    * which reduces the variation of the DEPENDENT variable

.

    Conditional mean = E[X|Y=y] = SUM_sx [ x * Pr(X=x|Y=y) ]        


    Conditional var = Var(X|Y=y) = SUM_sx [ (x - mu_x|Y=y)^2 * Pr(X=x|Y=y) ] 

    
### Independent Random Variables


    p(x,y) = p(x) * p(y)

    p(x|y) = p(x)   # knowledge of y does not affect p(x)

    p(y|x) = p(y)   # knowledge of x does not affect p(y)

    # derivation
    p(x|y) = p(x,y) / p(y)
        
           = p(x) * p(y) / p(y)
        
           = p(x)
                

## <a name="bvcd"></a>BiVariate Continuous Distributions


* The joint pdf of X and Y...
* is a non-negative function, f(x,y) 
    * 3-dimensional surface
* computed using a double integral
* across the range of X and the range of Y

.

    integral_INF [ integral_INF [ f(x,y) * dx * dy ] = 1

    Pr(x1 <= X <= x2, y1 <= Y <= y2) = integral_x1-x2 [ integral_y1-y2 [ f(x,y) * dx * dy ] ]

                                     = VOLUME under multi-dimensional probablity surface


Marginal Distributions:


    f(x) = integral_INF [ f(x,y) * dy ]     # integrate over Y


    f(y) = integral_INF [ f(x,y) * dx ]     # integrate over X


Conditional Distribution:

    f(x|y) = f(x,y) / f(y)

    f(y|x) = f(x,y) / f(x)


Conditional Mean:

    E[X|Y=y] = integral [ x * p(x|y) * dx ]

    E[Y|X=x] = integral [ y * p(y|x) * dy ]


Conditional Variance:

    Var(X|Y=y) = integral [ (x - mu_x|y)^2 * p(x|y) * dx ]

    Var(Y|X=x) = integral [ (y - mu_y|x)^2 * p(y|x) * dy ]


Independence:

    f(x|y) = f(x)

    f(y|x) = f(y)

    f(x,y) = f(x) * f(y)

* Independence propsitions above is useful in practice... 
* because it provides an easy way to construct the joint pdf
* from the product of the marginal pdf's

R: 

    mvtnorm # multi-variate normal distributions
    pmvnorm(lower=c(-1,-1), upper=c(1,1))


## <a name="cac"></a>Covariance and Correlation


* **Covariance**: measures direction
* but NOT strength
* of LINEAR relationship between two RV's
* does NOT say anything about any other type of relationship
    * e.g. quadratic relationship
    * difficult to measure non-linear relationship
* Covariance can be POSITIVE or NEGATIVE
    * unlike Variance, which can only be POSITIVE
* For a Graphical Description of Covariance...
    * Plot:
        * y-axis: `y-mu_y`
        * x-axis: `x-mu_x`

.

    sigma_xy = E[(X-mu_x)(Y-mu_y)]

             = SUM_s_xy (x - mu_x) * (y - mu_y) * p(x,y)


Noteworthy Properties of Covariance:


    Cov(X,Y) = Cov(Y,X)     # symmetric relationship

    Cov(aX,bY) = a * b * Cov(X,Y) = a * b * sigma_xy

    Cov(X,X) = Var(X)

    X,Y independent: Cov(X,Y) = 0

    Cov(X,Y) = 0 DOES NOT IMPLY X,Y independent

    Cov(X,Y) = E[X*Y] - E[X]*E[Y]

    recall: Var(X) = E[X*X] - E[X]*E[X]



* **Correlation**: measures direction 
* AND strength
* of LINEAR relationship between two RV's

.

                            Cov(X,Y)
    rho_xy = Cor(X,Y) = --------------------
                          sigma_x * sigma_y

                 sigma_xy
           = ------------------
             sigma_x * sigma_y


           = "scaled variance"

           = "Pearson correlation"

           

Noteworthy properties of correlation:


    Bounded: -1 <= rho_xy <= 1

    rho_xy = 1  if Y = aX + b and a > 0

    rho_xy = -1 if Y = aX + b and a < 0

    rho_xy = 0 if and only if sigma_xy = 0

    rho_xy = 0 DOES NOT IMPLY X and Y are independent

    rho_xy = 0 IMPLIES independence IFF X and Y are normally distributed
    

* Above is known as "Pearson Correlation"
* There are other formulas for computing correlation
    * e.g. "rank correlation"
    * which ranks the data points
    * and uses the rank in the correlation computation


### <a name="gbvnd"></a>General BiVariate Normal Distribution

* Big long equation
* depends on 5 parms:
    * `mu_x`
    * `mu_y`
    * `sigma_x`
    * `sigma_y`
    * `rho_xy`
* distribution is a 3-d liberty bell
    * like standard normal bivariate distribution
    * except centered on `mu_x, mu_y`
    * higher `rho_xy` results in "smushed" bell
    * i.e. narrower distribution
    * `rho_xy = 1` is a 2-d distribution


## <a name="lcotrv"></a>Linear Combination of Two Random Variables


* Define new random variable Z
* as a linear combination of X and Y
* with joint pdf p(x,y)

.

    Z = aX + bY


Mean/Expected Value:

    mu_z = E[Z] = a * E[X] + b * E[Y]

    Derivation:

    E[Z] = E[aX + bY]

         = SUM_x SUM_y (ax + by) * p(x,y)

         = SUM_x SUM_y ax * p(x,y) + SUM_x SUM_y by * p(x,y)

         = a * SUM_x SUM_y x * p(x,y) + b * SUM_x SUM_y y * p(x,y)

         = a * SUM_x x * SUM_y p(x,y) + b * SUM_y y * SUM_x p(x,y)

            SUM_y p(x,y) = p(x)     # marginal distribution

         = a * SUM_x x * p(x) + b * SUM_y * p(y)

            SUM_x x * p(x) = E[X]

         = a * E[X] + b * E[Y]


Variance:

    sigma_z^2 = Var(Z) = Var(aX + bY)

              = a^2 * Var(X) + b^2 * Var(Y) + 2*a*b*Cov(X,Y)

              = a^2 * sigma_x^2 + b^2 * sigma_y^2 + 2*a*b*sigma_xy



    Derivation:

    Var(Z) = Var(aX + bY)

           = E[ (Z - mu_z)^2 ]

           = E[ (aX + bY - (a*mu_x + b*mu_y))^2 ]

           = E[ (a(X - mu_x) + b(Y - mu_y))^2 ]

           = E[ a^2(X - mu_x)^2 + b^2(Y - mu_y)^2 + 2ab(X-mu_x)(Y-mu_y) ]

           = a^2 * E[(X-mu_x)^2] + b^2 * E[(Y-mu_y)^2] + 2ab * E[(X-mu_x)(Y-mu_y)]

           = a^2 * Var(X) + b^2 * Var(Y) + 2ab * Cov(X,Y)


Normality:

     If X is normal, and Y is normal, then Z is normal.

   

## <a name="epa"></a>Example: Portfolio Analysis


    R_a = return on asset A

    R_b = return on asset B

    sigma_ab = Cov(R_a, R_b)

                                 sigma_ab
    rho_ab = Cor(R_a, R_b) = -------------------
                             sigma_a * sigma_b


    x_a = share of wealth in asset A
    x_b = share of wealth in asset B
    
    x_a + x_b = 1 (all wealth invested)

    R_p = x_a * R_a + x_b * R_b 
        = portfolio return
        = linear combination of two RV's
        

    E[R_p] = x_a * E[R_a] + x_b * E[R_b]

    Var(R_p) = x_a^2 * Var(R_a) +
               x_b^2 * Var(R_b) +
               2 * x_a * x_b + Cov(R_a, R_b)

    stdev(R_p) = sqrt(Var(R_p))



Generalize to N RV's:

    Z = a1*X1 + a2*X2 + ... aN*XN = SUM_i=1..N [ ai * Xi ]

    E[Z] = SUM_i=1..N [ ai * E[Xi] ]


    Var(Z) = a1^2*Var(X1) + a2^2*Var(X2) + ...
             + all pairwise covariance terms

    N variance terms
    N * (N-1) covariance terms
    Way more covariance terms: biggest factor for portfolio variance


### Portfolio calculations using Matrix Algebra


Setup covariance matrix:

            xA          xB         xC
    xA  Var(xA)     Cov(xA,xB)  Cov(xA,xC)
    xB  Cov(xB,xA)  Var(xB)     Cov(xB,xC)
    xC  Cov(xC,xA)  Cov(xC,xB)  Var(xC)

More later...



## <a name="empccr"></a>Example: Multi-period continously compounded returns

* CC monthly returns
* normally distributed
* returns are UN-correlated with each other
* every month is a random draw from the normal distribution

.

    r = ln(1 + R) = monthly CC return

    r_t ~ N(mu, sigma^2) for all t

    Cov(r_t, r_s) = 0 for all t != s


Annual return:

    r_t(12) = SUM_j=0..11 r_t-j

           = r_t + r_t-1 + ... + r_t-11      # CC returns are additive

           = linear combination of UN-correlated monthly returns

    E[r_t(12)] = SUM_j=0..11 [ E[r_t-j] ]

               = 12 * mu

               = 12 * monthly return


    Var(r_t(12)) = Var( SUM_j=0..11 [ r_t-j ] )

                 = SUM_j=0..1 [ Var(r_t-j) ]

                 = 12 * sigma^2         # all covariance terms drop out, bc Cov = 0

    stdev(r_t(12)) = sqrt(12) * sigma   # "square-root of time rule"


    r_t(12) ~ N(12*mu, 12*sigma^2)


    # Quantiles scaled up from standard normal quantile, Z:

    q_alpha(r_t(12)) = 12*mu + sqrt(12)*sigma * q_alpha(Z)



* Square root of time rule, depends on:
    * same distribution every period
    * returns are uncorrelated across periods
    * used often in practice



