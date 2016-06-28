
# Coursera: Intro to Computational Finance: Week 9: Portfolio Theory Computations

* [Computing the Portfolio Frontier](#portfront)
* [Strategy for Plotting the Portfolio Frontier](#plotportfront)
* [Computing the Tangency Portfolio](#tanport)
* [Mutual Fund Separation Theorem (again)](#mutfund)
* [Portfolio Value-at-Risk](#portvar)
* [Portfolio Theory with NO Short Sales](#portnoshort)
* [Using solve.QP in R:](#solveqp)
* [Global Minimum Variance Portfolio w/ No Short Sales](#gmvportnoshort)
* [Finding the Efficient Frontier with No Short Sales](#efffrontnoshort)
* [R Functions / Packages for Portfolio Analysis](#rfuncs)
* [Statistical Analysis of Efficient Portfolios](#statsport)
* [Bootstrapping Efficient Portfolios](#bootport)
* [Rolling Efficient Portfolios Over Time](#rollport)



## <a name="portfront"></a>Computing the Portfolio Frontier

* We could find the portfolio frontier...
    * by using lagrangian multipliers
    * and a set of expected returns 
        * from the global min var portfolio on up
* However, there are easier ways:
    * The portfolio frontier can be represented as..
    * **convex combinations** of any two frontier portfolios

Let x be a frontier portfolio that solves:

    min(x) sigma_p,x^2 = x' * Sigma * x

    s.t. mu_p,x = x' * mu = mu_p,0

    x' * 1 = 1

Let y be another frontier portfolio:

    min(y) sigma_p,y^2 = y' * Sigma * y

    s.t. mu_p,y = y' * mu = mu_p,1

    y' * 1 = 1

Let alpha be any constant, then z is also a frontier portfolio:

    z = alpha * x + (1 - alpha) * y

    mu_p,z = z' * mu 
    
           = alpha * mu_p,x + (1 - alpha) * mu_p,y

    sigma_p,z^2 = z' * Sigma * z

                = alpha^2 * sigma_p,x^2 + (1 - alpha)^2 * sigma_p,y^2 + 2 * alpha * (1 - alpha) * sigma_x,y

    sigma_x,y = cov(R_p,x, R_p,y) 
        
             = x' * Sigma * y

* if alpha=1, z = x
* if alpha=0, z = y
* if alpha=0.5, then z is halfway between x and y
* if alpha=2, then z is a long-short combo of x and y
* note that portfolio x and y are essentially represented as single assets
* so the portfolio computations involving just two assets shows up again
    * can be applied to portfolios comprised of multiple assets

   
#### Example: Find efficient portfolio with expected return 0.05

    mu_p,z = alpha * mu_p,x + (1 - alpha) * mu_p,y = 0.05


             mu_p,z - mu_p,y
    alpha = ------------------ 
             mu_p,x - mu_p,y

    z = alpha * x + (1 - alpha) * y


## <a name="plotportfront"></a>Strategy for Plotting the Portfolio Frontier


1. 1st portfolio: start with global min var portfolio
2. 2nd portfolio: find the frontier portfolio with
    * expected return equal to max asset in the portfolio
3. grid of alpha values -1 to 1
    * compute all z frontier portfolios
4. plot means and variances to form the efficient frontier


R:

    # m.vec:        global min var portfolio weights
    # mu_p.m:       global min var portfolio mean
    # sigma_p.m:    global min var portfolio sd
    #
    # x.vec:        frontier portfolio #2 weights
    # mu_p.x:       frontier portfolio #2 mean
    # sigma_p.x:    frontier portfolio #2 sd
    #
    # sigma.mat     covariance matrix
    


    
    # cov between m and x:
    sigma_mx = t(m.vec) %*% sigma.mat %*% x.vec

    # set up alpha vector
    alpha.vec = seq(from=-1, to=1, by=0.1)
    alpha.n = length(alpha.vec)

    # output for frontier portfolios' weights, means, variances
    z.mat = matrix(0, nrow=alpha.n, ncol=3)     # weights
    mu_p.z.vec = rep(0, alpha.n)                      # means
    sigma2_p.z.vec = rep(0, alpha.n)                  # variances

    for (i in 1:alpha.n) {
        alpha = alpha.vec[i]

        z = alpha * m.vec + (1 - alpha) * x.vec
        mu_p.z.vec[i] = alpha * mu_p.m + (1 - alpha) * mu_p.x
        sigma_p.z.vec[i] = alpha^2 * sigma_p.m^2 + 
                           (1 - alpha) * sigma_p.x^2 +
                           2 * alpha * (1 - alpha) * sigma_mx
    }




## <a name="tanport"></a>Computing the Tangency Portfolio

* The tangency portfolio is relevent
* once you add a RISK-FREE asset
* The new portfolio is a combo of
    * the tangency portfolio of risky assets
    * the risk-free asset
* tangency portfolio has MAX sharpe slope of all portfolios
    * for a given risk-free rate
    * sharpe slope/tangency portfolio will be different
        * for different risk-free rates

Analytical solution for tangency portfolio:

    L(x,lambda) = (x' * mu - rf) * (x' * Sigma * x)^-1/2 + lambda(x' * 1 - 1)

    This is long and tedious to solve, but eventually you get to:
        
           Sigma^-1 %*% (mu - rf * 1)
    x = -----------------------------------
         1' %*% Sigma^-1 %*% (mu - rf * 1)

#### Remarks

* if risk-free rate is LESS THAN global min var portfolio "rate" (expected return)
    * then the tangency portfolio / sharpe slope is POSITIVE
* if risk-free rate EQUALS global min var portfolio rate
    * then there is NO tangency portfolio
* if risk-free rate is GREATER THAN global min var portfolio rate
    * then the sharpe slope is NEGATIVE
    * optimal returns involve SHORTING the tangency portfolio


### <a name="mutfund"></a>Mutual Fund Separation Theorem (again)

* Every efficient portfolio can be represented as...
* a convex combination of..
    * tangency portfolio
    * risk-free asset
* Given a desired level of risk or return,
    * you can compute the relative weights
    * of the tangency portfolio and the risk-free asset
        

Efficient Portfolios:

    xt = tangency portfolio weight
    xf = risk-free asset weight

    xt + xf = 1
    xf = 1 - xt

    mu_p,e = rf + xt * (mu_p,t - rf)

    mu_p,t = t' * mu

    sigma_p,e = xt * sigma_p,t

    sigma_p,t = sqrt( t' * Sigma * t )


### <a name="portvar"></a>Portfolio Value-at-Risk

    VaR_alpha = W0 * q_alpha(R_p)

    R_p: portfolio simple returns (a random variable with a distribution)

    q_alpha(R_p): alpha percentile of R_p distribution

    W0: initial investment

If returns `R_p` are normally distributed, then:

    R_p ~ iid N(mu_p, sigma_p^2)

    q_alpha(R_p) = mu_p + sigma_p * q_alpha(z)

For a given portfolio, x:

    mu_p,x = x' * mu

    sigma_p,x = sqrt( x' * Sigma * x )



## <a name="portnoshort"></a>Portfolio Theory with NO Short Sales


* exchanges may restrict short sales under certain conditions
* some institutions are prevented from short selling
* some accounts don't allow short sales
    * e.g. retirement accounts
* short selling often requires substantial credit requirements


#### Markowitz Formulation of the optimization problem:

    min(x) sigma_p,x^2 = x' * Sigma * x

    mu_p,x = x' * mu 
           
           = mu_p,0

    s.t:
        x' * 1 = 1

        x_i >= 0

* No analytic solution.
* Lagrangian multipliers / constrained optimization...
    * only works for equality constraints
    * not inequality constraints
* only Numeric solution
    * Excel: Solver
    * R: quadprog :: solve.QP
* efficient frontier cannot be constructed from any two efficient portfolios
    * must use brute force for every target expected return
* There may be no solution for some target expected returns
* frontier for no-short-sales lies INSIDE short-sales frontier


#### Quadratic Programming problems are of the form:

    min(x) 1/2 x' * D * x - d' * x

    A_eq'  * x  = b_eq,       for l equality constraints
    A_neq' * x >= b_neq,    for m inequality constraints

    D:      nxn matrix
    x,d:    nx1 vectors
    A_neq': mxn matrix
    b_neq:  mx1 vector
    A_eq':  lxn matrix
    b_eq:   lx1 vector


#### Applying to Portfolio Optimization

    D = 2*Sigma

    d = (0,...,0)'

m=2 equality constraints:

    x' * mu = mu_p,0

    x' * 1 = 1

    A_eq'  * x  = b_eq,       for l equality constraints

            [ mu' ]
    A_eq' = [ 1'  ]
    (2xn)

            [ mu_p,0 ]
    b_eq =  [ 1      ]
    (2x1)

l=n inequality constraints:

    x_i >= 0, i=1,...,n

    A_neq' * x >= b_neq,    for m inequality constraints

    A_neq' = I_n 
    (nxn)

    b_neq = (0,...,0)'
    (nx1)


### <a name="solveqp"></a>Using solve.QP in R:

* solve.QP assumes...
    * the inequality constraints
    * and equality constraints
    * are combined into a single matrix A'
    * and a single vector b

For the portfolio problem, we have:

          [ A_eq'  ]
    A' =  [ A_neq' ]

          [ b_eq  ]
    b  =  [ b_neq ]

         [ mu' ]            [ mu_1  1  1  0  0 ]
    A' = [ 1'  ]        A = [ mu_2  1  0  1  0 ]
         [ I_n ]            [ mu_n  1  0  0  1 ]

         [ mu_p,0 ]     # equality constraint
    b  = [   1    ]     # equality constraint
         [   0    ]     # inequality constraints - one for each asset weight

    qp.out <- solve.QP(Dmat=D.mat, 
                       dvec=d.vec, 
                       Amat=A.mat,
                       bvec=b.vec,
                       meq = 2)       # number of equality constraints
         

## <a name="gmvportnoshort"></a>Global Minimum Variance Portfolio w/ No Short Sales

Optimization problem:

    min(m) sigma_p,m^2 = m' * Sigma * m

    s.t:
        m' * 1 = 1

        m_i >= 0


    A_eq' = 1'      # only m=1 equality constraints
    (1xn)

    b_eq = 1
    (1x1)

    A_neq' = I_n    # l=n inequality constraints
    (nxn)

    b_neq = (0,...,0)'
    (nx1)

For solve.QP:

          [ A_eq'  ]
    A' =  [ A_neq' ]

          [ b_eq  ]
    b  =  [ b_neq ]


         [ 1'  ]
    A' = [ I_n ]

         [ 1 ]
    b  = [ 0 ]


R:

    D.mat = 2 * sigma.mat
    d.vec = rep(0,3)

    A.mat = cbind(rep(1,3), diag(3))   # A' comes out correct 
    b.vec = c(1, rep(0,3))

    args(solve.QP)
    qp.out <- solve.QP(D.mat, 
                       d.vec, 
                       A.mat,
                       b.vec,
                       meq = 1)       # number of equality constraints
    class(qp.out)
    names(qp.out)
    qp.out$solution     # portfolio weights
    qp.out$value        # portfolio variance



## <a name="efffrontnoshort"></a>Finding the Efficient Frontier with No Short Sales

* Start with global min var portfolio
* setup grid of expected returns between
    * global min var portfolio 
    * max individual asset return
* Solve the Markowitz algorithm with no short sales
    * for each target expected return in the grid
* Test to see if feasible solution exists
    * for target expected return GREATER THAN max individual asset return
    * how is that possible?

R:


    #       [ A_eq'  ]
    # A' =  [ A_neq' ]
    #
    #       [ b_eq  ]
    # b  =  [ b_neq ]
    #
    #      [ mu' ]            [ mu_1  1  1  0  0 ]
    # A' = [ 1'  ]        A = [ mu_2  1  0  1  0 ]
    #      [ I_n ]            [ mu_n  1  0  0  1 ]
    #
    #      [ mu_p,0 ]
    # b  = [   1    ]
    #      [   0    ]
    #
    # mu_p.gmv.ns:  mean return for global min var portfolio with no short sales
    # sigma.mat:    covariance matrix
    # mu.vec:       asset means

    k = 10          # number of portfolios to compute
    mu.vals = seq(from=mu_p.gmv.ns, to=max(mu.vec), length.out=k)

    N <- length(mu.vec)

    A.mat.t <- rbind( t(mu.vec),
                      rep(1,N),
                      diag(N) )
    A.mat <- t(A.mat.t)

    D.mat <- 2 * sigma.mat
    d.vec = rep(0,N)

    x.mat <- matrix(0, nrow=N, ncol=k)  # to be filled with portfolio weights
    colnames(x.mat) <- names(mu.vec)

    for (i in seq_along(mu.vals)) {
        mu.target <- mu.vals[i]
        b.vec <- c(mu.target, 1, rep(0,N))

        qp.out <- solve.QP(Dmat=D.mat,
                           dvec=d.vec,
                           Amat=A.mat,
                           bvec=b.vec,
                           meq=2)       # 2 equality constraints

        x.mat[i,] <- qp.out$solution
    }




## <a name="rfuncs"></a>R Functions / Packages for Portfolio Analysis


* tseries :: portfolio.optim()
* Rmetrics :: fPortfolio 
    * extensive collection of functions
* quadprog :: solve.QP
* PortfolioAnalytics  (PerformanceAnalytics)
* Eric Zivot's (prof) functions
* MY STUFF: [my.portfolio.r](my.portfolio.r)

Zivot's functions:

    my.portfolio = getPortfolio(er=mu.vec, 
                                cov.mat=sigma.mat,
                                weights=x.vec)

    attributes(my.portfolio)
    plot(my.portfolio)
    

    gmv.portfolio = globalMin.portfolio(er=mu.vec,
                                        cov.mat=sigma.mat)


    target.return = mu.vec["MSFT"]
    eff.portfolio = efficient.portfolio(er=mu.vec,
                                        cov.mat=sigma.mat,
                                        target.return)


    tan.portfolio = tangency.portfolio(er=mu.vec,
                                       cov.mat=sigma.mat,
                                       risk.free=rf)


    ef = efficient.frontier(er=mu.vec,
                            cov.mat=sigma.mat,
                            alpha.min=-2,
                            alpha.max=1.5,
                            nport=20)
    attributes(ef)
    plot(ef, plot.assets=T)


## <a name="statsport"></a>Statistical Analysis of Efficient Portfolios

The CER Model and Efficient Portfolios:

    i = asset, 1 .. N
    t = month, 1 .. T
    R_it = return on asset i in month t

    R_it  ~ iid N(mu_i, sigma_i^2)

    cov(R_it, R_jt) = sigma_ij

CER Parameters are estimated from historical data:

    mu_i     = ^mu_i
    sigma_i  = ^sigma_i
    sigma_ij = ^sigma_ij

* The estimates are RANDOM VARIABLES
    * and therefore are subject to STANDARD ERROR.
* Sharpe ratios and efficient portfolios...
    * are functions of the estimates
    * therefore they are also RANDOM VARIABLES
    * and are subject to STANDARD ERROR
* No easy analytical formulas for SE
* Use BOOTSTRAPPING!


#### Estimate of Sharpe Ratio

            ^mu_i - rf
    ^SR_i = -----------
             ^sigma_i

    SE(^SR_i): no easy formula


#### Estimate of Global Minimum Variance Portfolio

            ^Sigma^-1 * 1
    ^m = -------------------
          1' * ^Sigma^-1 * 1

    ^Sigma: estimated covariance matrix
    ^m:     estimated weights of the GMV portfolio
    
    SE(^m): SE(^m_i):   each element in m has its own SE.
    SE(^m_i): no easy formula

#### Estimate of GMV portfolio mean and variance

    ^mu_p,^m = ^m' * ^mu

    ^sigma_p,^m = sqrt( ^m' * ^Sigma * m )

    SE(^mu_p,^m):    no easy formula
    SE(^sigma_p,^m): no easy formula


## <a name="bootport"></a>Bootstrapping Efficient Portfolios

* Portfolio statistics to bootstrap:
    * weights
    * expected returns
    * standard deviations
* mean estimates of asset returns have large SE
    * therefore portfolios which depend on means also have large SE
* Things to bootstrap:
    * GMV portfolio
    * Any eff portfolio
    * Efficient frontier
    * Tangency portfolio

Strategy (e.g. bootstrapping GMV portfolio):

1. Re-sample from historical asset return data
    * to compute bootstrap distributions 
    * for mu.vec and Sigma.mat
2. Compute GMV portfolio for each bootstrap re-sample
3. analyze distributions of 
    * portfolio weights
    * portfolio mean and sd


## <a name="rollport"></a>Rolling Efficient Portfolios Over Time

* CER Model estimates are not constant over time
* Rolling estimates of mean, sd, and covariance change over time
* Things to compute from rolling estimates:
    * GMV portfolio
    * Any eff portfolio
    * Efficient frontier
    * Tangency portfolio

