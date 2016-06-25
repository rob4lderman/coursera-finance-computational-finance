
# Coursera: Intro to Computational Finance: Week 8: Portfolio Theory

* [Portfolio Theory](#porttheory)
* [Portfolio Value-at-Risk](#portvar)
* [Portfolio Frontier](#portfront)
* [Global Minimum Variance Portfolio](#gmvport)
* [Portfolio Theory with a Risk-Free Asset](#portriskfree)
* [Portfolio Frontier: Risk-free Asset + Risky Asset, Sharpe Ratio](#portriskfreerisky)
* [Tangency Portfolio: Risk-Free Asset + 2 Risky Assets](#portriskfreetworisky)
* [Efficient Portfolio Examples](#effport)
* [Matrix Algebra: Portfolio w/ More Than 2 Assets](#matrixport)
* [Covariance Between Two Portfolios](#portcov)
* [Matrix Algebra: Computing Global Minimum Variance Portfolio](#matrixgmvport)
* [Matrix Algebra: Finding Efficient Portfolios](#matrixeffport)

## <a name="porttheory"></a>Portfolio Theory

Investment in two Risky Assets

    RA: SIMPLE return for asset A
    RB: SIMPLE return for asset B
    W0: initial wealth


* Note: We use SIMPLE returns...
    * because then the portfolio return can be computed as...
    * the weighted avg of the simple returns
* This is not true if we use CC returns
    * however some analysts use CC returns anyway
    * not exact but close enough
* CER Model, however, applies to CC returns
    * NOT simple returns
    * so we assume CER model applies to simple returns as well

#### Assumptions:

    RA and RB are described by the CER model

    Ri ~ iid N(mu_i, sigma_i^2)

    cov(RA,RB) = sigma_AB

    cor(RA,RB) = rho_AB

    xA: share of wealth in A
    xB: share of wealth in B

    Long position: xA, xB > 0

    Short position: xA < 0 or xB < 0

    Allocate all wealth between assets: xA + xB = 1


* Investors like high E[Ri]
* Investors like low Var(Ri)
* Investment horizon is one period
    * e.g. one month
    * one year

#### Portfolio Return: 

    Rp = xA * RA + xB * RB

    E[Rp] = mu_p = xA * mu_A + xB * mu_B


#### Portfolio Distribution:

    Var(Rp) = sigma_p^2 = xA^2 * sigma_A^2 +
                          xB^2 * sigma_B^2 +
                          2 * xA * xB * sigma_AB

    Rp ~ iid N(mu_p, sigma_p^2)


* Note: we are assuming all asset returns are normally distributed
* Portfolio return, Rp, is linear combination of assets
* Therefore, Rp is also normally distributed
* it's possible to use different distributions
    * e.g. student's t


#### End of Period Wealth:

    W1 = (1 + Rp) * W0

       = W0 * (1 + xA*RA + xB * RB)

    W1 ~ N(W0*(1+mu_p), sigma_p^2 * W0^2)


## Example: Portfolio Calculations

Note: Asset A has higher risk+return than B

    mu_A = 0.175
    mu_B = 0.055

    sigma_A = sqrt(0.067)
    # [1] 0.2588436
    sigma_B = sqrt(0.013)
    # [1] 0.1140175

    sigma_AB = -0.004875

    rho_AB = sigma_AB / sigma_A / sigma_B
    
#### Example: equally weighted portfolio

    xA = xB = 0.5
    mu_p = xA * mu_A + xB * mu_B

    sigma2_p = xA^2 * sigma_A^2 +
               xB^2 * sigma_B^2 +
               2 * xA * xB * sigma_AB
    sigma_p = sqrt(sigma2_p)

    weighted_sigma_p <- xA * sigma_A + xB * sigma_B

    c("mu_p"=mu_p, "sigma_p"=sigma_p, "weighted_sigma_p"=weighted_sigma_p)
    #        mu_p          sigma_p weighted_sigma_p 
    #   0.1150000        0.1325236        0.1864306

* Note: `sigma_p` is LESS THAN weighted avg of sigma's
* This is because the assets are not perfectly correlated
* This represents the **Diversification Benefit**


#### Example: Long-short portfolio:
    
    # short B (sell now) and use proceeds to invest in A
    xA = 1.5
    xB = -0.5

    mu_p = xA * mu_A + xB * mu_B

    sigma2_p = xA^2 * sigma_A^2 +
               xB^2 * sigma_B^2 +
               2 * xA * xB * sigma_AB
    sigma_p = sqrt(sigma2_p)

    weighted_sigma_p <- xA * sigma_A + xB * sigma_B

    c("mu_p"=mu_p, "sigma_p"=sigma_p, "weighted_sigma_p"=weighted_sigma_p)
    #        mu_p          sigma_p weighted_sigma_p 
    #   0.2350000        0.4016373        0.3312566

* Note: portfolio has higher return and SD than asset A


## <a name="portvar"></a>Portfolio Value-at-Risk

Since Rp is normally distributed, VaR is straighforward:

    VaR_p,alpha = q_p,alpha * W0

                = (mu_p + sigma_p * q_alpha(z)) * W0


#### R:

    W0 = 1e5

    VaR_A.05 = qnorm(0.05, mean=mu_A, sd=sigma_A) * W0
    VaR_B.05 = (mu_B + sigma_B * qnorm(0.05)) * W0

    # For equally weighted portfolio:
    xA = 0.5
    xB = 0.5

    mu_p = xA * mu_A + xB * mu_B
    sigma2_p = xA^2 * sigma_A^2 +
               xB^2 * sigma_B^2 +
               2 * xA * xB * sigma_AB
    sigma_p = sqrt(sigma2_p)
    
    VaR_p.05 = (mu_p + sigma_p * qnorm(0.05)) * W0

    # Weighted avg of VaR_A and VaR_B
    VaR_AB.05 = xA * VaR_A.05 + xB * VaR_B.05

    c("VaR_A.05"=VaR_A.05, 
      "VaR_B.05"=VaR_B.05, 
      "VaR_p.05"=VaR_p.05, 
      "VaR_AB.05"=VaR_AB.05)
    #  VaR_A.05  VaR_B.05  VaR_p.05 VaR_AB.05 
    # -25075.98 -13254.22 -10298.19 -19165.10 

* Note: Portfolio VaR is NOT a weighted avg of asset VaR
    * UNLESS `rho_AB=1`
* Portfolio VaR is actually LESS THAN both asset's VaR
    * reflects diversification benefit



## <a name="portfront"></a>Portfolio Frontier

* Vary portfolio weights
* compute mean and variance
* plot mean (y-axis) against variance (x-axis)
* shape of frontier depends on correlations
* if correlation = -1...
    * there exists a portfolio with var=0
* if correlation - 1...
    * there is NO DIVERSIFICATION BENEFIT
* Diversification Benefit even when correlation is positive

#### R:

    xA = seq(from=-0.4, to=1.4, by=0.1)
    xB = 1 - xA
    mu_p = xA * mu_A + xB * mu_B
    sigma2_p = xA^2 * sigma_A^2 +
               xB^2 * sigma_B^2 +
               2 * xA * xB * sigma_AB
    sigma_p = sqrt(sigma2_p)

    plot(x=sigma_p, 
         y=mu_p, 
         type="b",
         pch=16,
         ylim=c(0,max(mu_p)),
         xlim=c(0,max(sigma_p)),
         main="Portfolio Frontier",
         xlab=expression(sigma[p]),
         ylab=expression(mu[p]),
         col=c(rep("red",6), rep("green",13)))  # red are "dominated" portfolios
    text(x=sigma_A, 
         y=mu_A, 
         labels="Asset A",
         pos=4)
    text(x=sigma_B, 
         y=mu_B, 
         labels="Asset B",
         pos=4)
    text(x=0.1*max(sigma_p),
         y=0.9*max(mu_p),
         labels=c(expression(rho[AB]),paste("=",rho_AB)),
         pos=c(2,4))


* Long-short portfolios are beyond the data points marked Asset A and Asset B
* Green portfolios are **efficient portfolios**
    * Portfolios with the highest expected return 
    * for a given level of risk
* Red portfolios are "dominated"/inefficient portfolios
* all portfolios "above" the min-variance portfolio are "efficient"
* shape of efficient frontier depends on correlation coefficient
    * narrower parabolas for lower/negative correlations
        * corr = -1: 
            * parabola becomes a triangle
            * min-variance is 0
    * wider parabolas for higher/positive correlations
        * corr = 1:
            * a straight line
            * no diversification benefit
* Investor risk preferences:
    * risk averse will hold portfolios near the min variance 
    * risk tolerant will hold portfolios with higher expected return


## <a name="gmvport"></a>Global Minimum Variance Portfolio


    min(xA,xB) sigma_p^2 = xA^2 * sigma_A^2 + 
                           xB^2 * sigma_B^2 + 
                           2 * xA * xB * sigma_AB

    constraint: xA + xB = 1


### Review: Constrained Optimization

             y = f(x,z) = x^2 + z^2

    min(x,z) y = f(x,z)

    s.t. x + z = 1


Solving Methods:

* Analytical
    * Substitution
        * substitute one var in terms of the other
    * Lagrange Multipliers
        * augment function with extra term 
        * extra term imposes constraints
* Numerical
    * Excel: Solver
    * R: quadprog :: solve.QP


### Method: Substitution

    z = 1 - x

    y = f(x, x-1) = x^2 + (1-x)^2


First order conditions:
        
         d
    0 = --- (x^2 + (1 - x)^2) 
         dx

      = 2x + 2 * (1 - x) * -1

      = 4x - 2

      x = 0.5
      z = 0.5


### Method: Lagrange Multipliers

1\. Put constraints in homogenous form

    x + z = 1  =>  x + z - 1 = 0

2\. Form Lagrangian function

    L(x,z,lambda) = x^2 + z^2 + lambda(x + z - 1)

    lambda = Lagrange multiplier

3\. Minimize Lagrangian function

    min(x,z,lambda) L(x,z,lambda)


First order conditions:

         dL
    0 = ---- = 2x + lambda
         dx

         dL
    0 = ---- = 2z + lambda
         dz

            dL
    0 = ---------- = x + z - 1
         d-lambda


Three linear equations, three unknowns:


    2x = -lambda

    2z = -lambda

    x = z

    2z - 1 = 0

    z = 0.5

    x = 0.5


## <a name="portriskfree"></a>Portfolio Theory with a Risk-Free Asset

* Risk-free asset:
    * asset with fixed, known rate of return over investment horizon
    * usually use US t-bill (under 1yr) 
    * or t-bond (over 1yr)
    * t-bill / t-bond rate is only NOMINALLY risk free
        * not real
        * not accounting for inflation 

#### Properties of risk-free asset:

    Rf = return

    E[Rf] = rf = constant

    Var(Rf) = 0

    cov(Rf,Ri) = 0,     Ri: any other asset


#### Portfolios of Risk Asset and Risk-Free Asset

    xf = share of risk-free

    xB = share of wealth in asset B

    xf + xB = 1

#### Portfolio Return:

    Rp = xf * rf + xB * RB

       = (1 - xB) * rf + xB * RB

       = rf + xB * (RB - rf)        # (RB - rf): excess return on asset B

#### Portfolio excess return:

    Rp - rf = xB * (RB - rf)


#### Portfolio Distribution:

    mu_p = E[Rp] = rf + xB * (mu_B - rf)

    sigma_p^2 = Var(Rp) = xB^2 * sigma_B^2      

    sigma_p = xB * sigma_B

    Rp ~ N(mu_p, sigma_p^2)


* **Risk Premium**: `mu_b - rf`
    * excess return over risk free
    * premium to compensate for additional risk
* **Leveraged Investment**
    * borrow at T-bill rate to buy more of asset B
        * in reality, can't borrow at t-bill rate
    * increases expected return and risk


## <a name="portriskfreerisky"></a>Portfolio Frontier: Risk-free Asset + Risky Asset, Sharpe Ratio

    sigma_p = xB * sigma_B

    xB = sigma_p / sigma_B

    mu_p = rf + xB * (mu_B - rf)

         = rf + sigma_p / sigma_B * (mu_B - rf)

         = rf + (mu_B - rf)/sigma_B  * sigma_p


     mu_B - rf
    ----------- = SR_B = Asset B Sharpe Ratio
      sigma_B

                = excess expected return per unit of risk


* **Sharpe ratio** is commonly used to rank assets
* Higher sharpe ratios are preferred over lower Sharpe ratio
* because it's a coefficient of expected return per unit of risk
    * the bigger the coefficient,
    * the more return per increased unit of risk
    * "bang for your buck"

R:

    rf = 0.03

    xA = seq(from=0, to=1.4, by=0.1)

    mu_p.A = rf + xA * (mu_A - rf)
    sigma_p.A = xA * sigma_A
    sharpe.A = (mu_A - rf) / sigma_A
    # [1] 0.5601839

    xB = seq(from=0, to=1.4, by=0.1)

    mu_p.B = rf + xB * (mu_B - rf)
    sigma_p.B = xB * sigma_B

    sharpe.B = (mu_B - rf) / sigma_B
    # [1] 0.2192645


    plot(x=sigma_p.A, 
         y=mu_p.A, 
         type="b",
         pch=16,
         ylim=c(0,max(mu_p.A, mu_p.B)),
         xlim=c(0,max(sigma_p)),
         main="Portfolios of T-Bills + 1 Risky Asset",
         xlab=expression(sigma[p]),
         ylab=expression(mu[p]),
         col="green")  
    text(x=sigma_A, 
         y=mu_A, 
         labels="Asset A",
         pos=4)
    lines(x=sigma_p.B,
          y=mu_p.B,
          type="b",
          pch=16,
          col="red")
    text(x=sigma_B, 
         y=mu_B, 
         labels="Asset B",
         pos=4) 

* All portfolios of T-Bills + Asset A dominate T-Bills + Asset B
* Higher Sharpe Ratio
    * which is why sharpe ratio is often used to rank assets


## <a name="portriskfreetworisky"></a>Tangency Portfolio: Risk-Free Asset + 2 Risky Assets

* Collapse 2 Risky Assets into single portfolio
* so the problem becomes T-Bills + portfolio (single risky asset)
* max Sharpe Ratio for t-Bills + portfolio is...
    * the line from rf to the tangency of the efficient frontier
    * "tangency portfolio"
* T-Bills + Tangency portfolio dominates all other portfolios


### Mutual Fund Separation Theorem

* Efficient portfolios are combinations of two portfolios (mutual funds):
    * T-Bill portfolio
    * Tangency Portfolio 
        * portfolio of assets A and B
        * with the maximum Sharpe Ratio
* Implication: all investors hold assets A and B
    * according to proportions in the tangency portfolio
    * regardless of risk preference


### Computing the tangency portfolio


                        mu_p - rf
    max(xA,xB) SR_p = -------------
                         sigma_p

    mu_p = xA * mu_A + xB * mu_B
         = x' * mu

    sigma_p^2 = xA^2 * sigma_A^2 + 
                xB^2 * sigma_B^2 + 
                2 * xA * xB * sigma_AB
              = x' * Sigma * x

    s.t. xA + xB = 1

This can be solved numerically or analytically.

    L(x,lambda) = (x' * mu - rf) * (x' * Sigma * x)^-1/2 + lambda(x' * 1 - 1)

    This is long and tedious to solve, but eventually you get to:
        
           Sigma^-1 * (mu - rf * 1)
    x = -----------------------------
         1' * Sigma^-1 * (mu -rf * 1)



#### R Example:

    xA.tan = 0.46
    xB.tan = 0.54

    mu_p.tan = xA.tan * mu_A + xB.tan * mu_B
    # [1] 0.1102

    sigma2_p.tan = xA.tan^2 * sigma_A^2 + 
                   xB.tan^2 * sigma_B^2 + 
                   2 * xA.tan * xB.tan * sigma_AB

    sigma_p.tan = sqrt(sigma2_p.tan)
    # [1] 0.124684


Efficient portfolios have the following characteristics:

    mu_p.eff = rf + x_p.tan( mu_p.tan - rf)

    sigma_p.eff = x_p.tan * sigma_p.tan

    # x_p.tan = share of wealth invested in tangency portfolio
    #           (the rest in T-Bills)


Plot everything:

    x_p.tan <- seq(from=0, to=1.4, by=0.1)
    mu_p.eff <- rf + x_p.tan * (mu_p.tan - rf)
    sigma_p.eff <- x_p.tan * sigma_p.tan
    

    plot(x=sigma_p, 
         y=mu_p, 
         type="b",
         pch=16,
         ylim=c(0,max(mu_p)),
         xlim=c(0,max(sigma_p)),
         main="Portfolio Frontier",
         xlab=expression(sigma[p]),
         ylab=expression(mu[p]),
         col="black")
    text(x=sigma_A, 
         y=mu_A, 
         labels="Asset A",
         pos=4)
    text(x=sigma_B, 
         y=mu_B, 
         labels="Asset B",
         pos=4)
    lines(x=sigma_p.A,
          y=mu_p.A,
          type="b",
          pch=16,
          col="blue")
    lines(x=sigma_p.B,
          y=mu_p.B,
          type="b",
          pch=16,
          col="red")
    lines(x=sigma_p.eff,
          y=mu_p.eff,
          type="b",
          pch=16,
          col="green")
    text(x=sigma_p.tan,
         y=mu_p.tan,
         labels="Tangency",
         pos=4)


## <a name="effport"></a>Efficient Portfolio Examples

Find the efficient portfolio with same risk as asset B:

    mu_p.eff = rf + x_p.tan( mu_p.tan - rf)

    sigma_p.eff = x_p.tan * sigma_p.tan

    # x_p.tan = share of wealth invested in tangency portfolio
    #           (the rest in T-Bills)

    # Set sigma_p.eff equal to sigma_B and solve for x_p.tan:

    sigma_B = x_p.tan * sigma_p.tan

    x_p.tan = sigma_B / sigma_p.tan
    # [1] 0.9144521

    xf = 1 - x_p.tan
    # [1] 0.08554793

    
Find the efficient portfolio with same expected return as asset B:

    mu_p.eff = rf + x_p.tan( mu_p.tan - rf)

    sigma_p.eff = x_p.tan * sigma_p.tan

    # Set mu_p.eff equal to mu_B and solve for x_p.tan:

    mu_B = rf + x_p.tan( mu_p.tan - rf)

    x_p.tan = (mu_B - rf) / (mu_p.tan - rf)
    # [1] 0.3117207

    xf = 1 - x_p.tan
    # [1] 0.6882793
    


## <a name="matrixport"></a>Matrix Algebra: Portfolio w/ More Than 2 Assets


    RA, RB, RC

    xA + xB + xC = 1

    Rp = xA * RA + xB * RB + xC * RC

    mu_p = xA * mu_A + xB * mu_B + xC * mu_C

    Cov(Ri,Rj) = sigma_ij

In Matrix Notation:

        [ xA ]
    x = [ xB ]
        [ xC ]
        
        [ RA ]
    R = [ RB ]
        [ Rc ]

         [ mu_A ]
    mu = [ mu_B ]
         [ mu_C ]

            [ sigma_A^2     sigma_AB    sigma_AC    ]
    Sigma = [ sigma_BA      sigma_B^2   sigma_BC    ]
            [ sigma_CA      sigma_CB    sigma_C^2   ]

                Note: sigma_BA = sigma_AB


Portfolio Weights Sum To One:

    x' * 1 = 1

Portfolio Return:

    x' * R 

Portfolio Expected Return:

    x' * mu

    R: t(x) %*% mu

       ..or..

       crossprod(x,mu)


Portfolio Variance:

    x' * Sigma * x

    R: t(x) %*% Sigma %*% x


Portfolio Distribution:

    R_p,x ~ N(mu_p,x, sigma_p,x^2)


### <a name="portcov"></a>Covariance Between Two Portfolios

Two portfolios of different weights: x and y

        [ xA ]
    x = [ xB ]
        [ xC ]

        [ yA ]
    y = [ yB ]
        [ yC ]
   
Portfolio Returns:

    x' * R
    y' * R

Covariance:

    cov(R_p,x, R_p,y) = x' * Sigma * y
                
                      = y' * Sigma * x


### Matrix: Example: R


    asset.names <- c("MSFT", "NORD", "SBUX")
    mu.vec <- c(0.0427, 0.0015, 0.0285)
    names(mu.vec) <- asset.names

    sigma.mat <- matrix(c(0.0100, 0.0018, 0.0011,
                          0.0018, 0.0109, 0.0026,
                          0.0011, 0.0026, 0.0199),
                        nrow=3,
                        ncol=3)
    dimnames(sigma.mat) <- list(asset.names, asset.names)

Plot risk-return characteristics of the assets:

    sigma.vec <- sqrt(diag(sigma.mat))
    plot(x=sigma.vec,
         y=mu.vec,
         pch=16,
         main="Risk/Return Characeteristics",
         xlim=c(0,max(sigma.vec) * 1.5),
         ylim=c(0,max(mu.vec) * 1.5),
         xlab=expression(sigma[p]),
         ylab=expression(mu[p]),
         col="slateblue1")

    for (i in seq_along(asset.names)) {
        text(x=sigma.vec[i],
             y=mu.vec[i],
             labels=asset.names[i],
             pos=4)
    }


Equally Weighted Portfolio:

    x.vec <- rep(1,3)/3
    names(x.vec) <- asset.names

    mu_p.x <- t(x.vec) %*% mu.vec
    # or: mu_p.x <- crossprod(x.vec, mu.vec)

    sigma2_p.x <- t(x.vec) %*% sigma.mat %*% x.vec

    sigma_p.x <- sqrt(sigma2_p.x)

    sum(x.vec)
    mu_p.x
    sigma_p.x


Long-short Portfolio:

    y.vec <- c(0.8, 0.4, -0.2)
    names(y.vec) <- asset.names

    mu_p.y <- t(y.vec) %*% mu.vec
    # or: mu_p.y <- crossprod(y.vec, mu.vec)

    sigma2_p.y <- t(y.vec) %*% sigma.mat %*% y.vec

    sigma_p.y <- sqrt(sigma2_p.y)

    sum(y.vec)
    mu_p.y
    sigma_p.y


Add portfolios to plot:

    points(x=c(sigma_p.x, sigma_p.y),
           y=c(mu_p.x, mu_p.y),
           pch=16,
           col="black")
    text(x=sigma_p.x,
         y=mu_p.x,
         labels="EQUAL-WEIGHT",
         pos=4)
    text(x=sigma_p.y,
         y=mu_p.y,
         labels="LONG-SHORT",
         pos=2)


Covariance and correlation between portfolios:

    sigma_p.xy = t(x.vec) %*% sigma.mat %*% y.vec

    rho_xy = sigma_p.xy / sigma_p.x / sigma_p.y


* TODO: Generate random portfolios and plot them
* Note: portfolios mean/var plot with 2 assets...
    * lies on line
    * portfolio frontier curve
* Portfolios mean/var plot with 3+ assets...
    * fills the space
    * outlines by the frontier curve


### Matrix Calculus

     d
    --- x' * A * x =  2 * A * x
     dx

     d
    --- x' * y = y
     dx


## <a name="matrixgmvport"></a>Matrix Algebra: Computing Global Minimum Variance Portfolio


         [ mA ]
     m = [ mB ]
         [ mC ]

     min(mA,mB,mC) sigma2_p,m = m' * Sigma * m

     s.t: m' * 1 = 1

1. Analytical solution via matrix algebra
    * Lagrangian
2. Numerical solution 
    * excel solver

### Analytical solution: Lagrangian 


    L(m,lambda) = m' * Sigma * m + lambda(m' * 1 - 1)

First Order conditions:

         d-L(m,lambda)     d                 d
    0 = -------------- =  ---- m'*Sigma*m + --- lambda(m'*1 - 1)
             dm            dm                dm

      = 2*Sigma*m + lambda*1


         d-L(m,lambda)        d                     d
    0 = -------------- =  --------  m'*Sigma*m + --------- lambda(m'*1 - 1)
           d-lambda       d-lambda                d-lambda

      = m'*1 - 1


Write First Order Conditions in Matrix Form:

    [ 2*Sigma   1 ]     [ m      ]   [ 0 ]  (3x1)
    [ 1'        0 ]  *  [ lambda ] = [ 1 ]  (1x1)


Linear System:

    Am * zm = b

         [ 2*Sigma   1 ]
    Am = [ 1'        0 ]

         [ m      ]
    zm = [ lambda ]

        [ 0 ]
    b = [ 1 ]

    zm = Am^-1 * b

R:

    top.mat <- cbind(2 * sigma.mat, rep(1,3))
    bot.mat = c(rep(1,3), 0)
    Am.mat <- rbind(top.mat, bot.mat)

    b.vec = c(rep(0,3), 1)

    zm.mat <- solve(Am.mat) %*% b.vec
    m.vec <- zm.mat[1:3,1]

    mu_p.min <- t(m.vec) %*% mu.vec
    sigma2_p.min <- t(m.vec) %*% sigma.mat %*% m.vec
    sigma_p.min <- sqrt(sigma2_p.min)

Plot it:

    points(x=sigma_p.min,
           y=mu_p.min,
           pch=16,
           col="black")
    text(x=sigma_p.min,
         y=mu_p.min,
         labels="GLOBAL MIN",
         pos=2)



## <a name="matrixeffport"></a>Matrix Algebra: Finding Efficient Portfolios


Two ways to formulate the problem:

1. Find portfolio with highest expected return for a given level of risk
2. Find portfolio with smallest risk for a given level of expected return
    * easier than #1
    * constrained optimization w/ 2 constraints

Two methods for solving:

1. Analytically 
    * matrix algebra
    * lagrangian
2. Numerically
    * solver in Excel


Formulation 2\.:

    min(xA,xB,xC) sigma2_p,x = x' * Sigma * x

    s.t:
        
        mu_p,x = x' * mu = mu_p,t = target return

        x'*1 = 1

IN homogenous form:

    0 = x'*mu - mu_p,t

    0 = x'*1 - 1 

Lagrangian:

    L(x,lambda1, lambda2) = x'*Sigma*x + 
                            lambda1 * [ x'*mu - mu_p,t ] +
                            lambda2 * [ x'*1 - 1 ]

* Note: when you take the derivative wrt lambda1...
    * you force the first constraint to hold
* when you take the derivative wrt lambda2...
    * you force the second constraint to hold

First Order Conditions:

         d-L
    0 = ----- = 2*Sigma*x + lambda1*mu + lambda2*1
         dx

            d-L
    0 = ----------- = x'*mu - mu_p,t
         d-lambda1

            d-L
    0 = ----------- = x'*1 - 1
         d-lambda2

Five Linear Equations, with Five Unknowns:

    (xA,xB,xC,lambda1,lambda2)


In Matrix Notation:

    [ 2*Sigma  mu  1 ]   [    x    ]    [   0    ]
    [ mu'      0   0 ] * [ lambda1 ] =  [ mu_p,t ]
    [ 1'       0   0 ]   [ lambda2 ]    [   1    ]

            Ax         *     zx      =      b0

    zx = Ax^1 * b0


#### Example: Find Efficient Portfolio with same mean as MSFT:

    top.mat <- cbind(2*sigma.mat, mu.vec, rep(1,3))
    mid.vec <- c(mu.vec, 0, 0)
    bot.vec <- c(rep(1,3), 0, 0)

    Ax.mat <- rbind(top.mat, mid.vec, bot.vec)

    b_msft.vec <- c(rep(0,3), mu.vec["MSFT"], 1)

    zx.mat = solve(Ax.mat) %*% b_msft.vec
    x.vec <- zx.mat[1:3,1]

    mu_p.eff1 <- t(x.vec) %*% mu.vec
    sigma2_p.eff1 <- t(x.vec) %*% sigma.mat %*% x.vec
    sigma_p.eff1 <- sqrt(sigma2_p.eff1)


Plot it:


    points(x=sigma_p.eff1,
           y=mu_p.eff1,
           pch=16,
           col="black")
    text(x=sigma_p.eff1,
         y=mu_p.eff1,
         labels="EFF1",
         pos=2)


#### Example: Find Efficient Portfolio with same mean as SBUX:
    
    b_sbux.vec <- c(rep(0,3), mu.vec["SBUX"], 1)

    zx.mat = solve(Ax.mat) %*% b_sbux.vec
    x.vec <- zx.mat[1:3,1]

    mu_p.eff2 <- t(x.vec) %*% mu.vec
    sigma2_p.eff2 <- t(x.vec) %*% sigma.mat %*% x.vec
    sigma_p.eff2 <- sqrt(sigma2_p.eff2)


Plot it:


    points(x=sigma_p.eff2,
           y=mu_p.eff2,
           pch=16,
           col="black")
    text(x=sigma_p.eff2,
         y=mu_p.eff2,
         labels="EFF2",
         pos=2)





