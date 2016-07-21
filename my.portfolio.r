#
# Functions for Portfolio Analysis 
#
#
# gmv.weights.vec <- gmvPortfolio( mu.vec, Sigma.mat )
# tan.weights.vec <- tanPortfolio( mu.vec, Sigma.mat, risk.free.rate )
# eff.weights.vec <- effPortfolio( mu.vec, Sigma.mat, target.return )
# eff.weights.mat <- effFrontier(mu.vec, Sigma.mat, alpha.vec=seq(from=-1,to=2,by=0.1), gmv.portfolio.weights=NULL, eff.portfolio.weights=NULL) 
#
# gmv.weights.vec <- gmvPortfolio.noShort( mu.vec, Sigma.mat )
# tan.weights.vec <- tanPortfolio.noShort( mu.vec, Sigma.mat, risk.free.rate )
# eff.weights.vec <- effPortfolio.noShort( mu.vec, Sigma.mat, target.return )
# eff.weights.mat <- effFrontier.noShort(mu.vec, Sigma.mat, target.return.vec )
#
# portfolioWeightsBarplot(x.mat, ...)
#


gmvPortfolio <- function( mu.vec, Sigma.mat) {
    # Parameters:
    #   mu.vec: expected returns for all assets
    #   Sigma.mat: covariance matrix
    #
    # Return Value:
    #   x.vec: vector of weights for each asset in the global-min-var portfolio
    #
    #
    # Constrained optimization: Lagrangian Multipliers
    #
    # A.mat %*% z.vec = b.vec
    # z.vec = A.mat^-1 %*% b.vec
    #
    # A.mat = [ 2*Sigma  1 ]
    #         [ 1'       0 ]
    #
    # b.vec = [ 0 ]
    #         [ 0 ]
    #         [ 1 ]
    #
    # z.vec = [ x.vec  ]
    #         [ lambda ]
    N <- length(mu.vec)

    A.mat.top = cbind( 2*Sigma.mat, rep(1,N) )
    A.mat = rbind( A.mat.top, c(rep(1,N),0) )

    b.vec = c(rep(0,N),1)

    z.vec = solve(A.mat) %*% b.vec
    x.vec = z.vec[1:N,1]
    x.vec
}


gmvPortfolio.noShort <- function(mu.vec, Sigma.mat) {
    # Parameters:
    #   mu.vec:     expected returns for all assets
    #   Sigma.mat:  covariance matrix
    #
    # Return Value:
    #   x.vec: vector of weights for each asset in the global-min-var portfolio
    #
    # Quadratic Programming problem with inequality constraints.
    # Numerical solutions only! 
    # Makes use of solve.QP
    #
    #       [ A_eq'  ]
    # A' =  [ A_neq' ]
    #
    #       [ b_eq  ]
    # b  =  [ b_neq ]
    #
    #                         [ 1  1  0  0 ]
    # A' = [ 1'  ]        A = [ 1  0  1  0 ]
    #      [ I_n ]            [ 1  0  0  1 ]
    #
    # b  = [   1    ]
    #      [   0    ]
    #      [   0    ]
    #      [   0    ]
    #
    # D = 2*Sigma
    #
    # d = (0,...,0)'
    #

    N <- length(mu.vec)
    D.mat <- 2 * Sigma.mat
    d.vec = rep(0,N)

    A.mat.t <- rbind( rep(1,N),
                      diag(N) )
    A.mat <- t(A.mat.t)

    b.vec <- c(1, rep(0,N))

    qp.out <- solve.QP(Dmat=D.mat,
                       dvec=d.vec,
                       Amat=A.mat,
                       bvec=b.vec,
                       meq=1)       # 1 equality constraints (x'*1=1)
    x.vec <- qp.out$solution
    names(x.vec) <- names(mu.vec)
    x.vec
}



tanPortfolio <- function(mu.vec, Sigma.mat, risk.free.rate) {
    # Parameters:
    #   mu.vec: asset expected returns
    #   Sigma.mat: asset covariance matrix
    #   risk.free.rate: risk-free rate
    #
    # Return Value:
    #   t.vec: tangency portfolio weights
    #
    # Constrained optimization, maximize sharpe ratio:
    #
    #               Sigma.mat^-1 %*% (mu.vec - risk.free.rate * 1)
    #     t.vec = --------------------------------------------
    #               1' %*% Sigma.mat^-1 %*% (mu.vec - risk.free.rate * 1)

    rf <- risk.free.rate
    one.vec <- rep(1,length(mu.vec))
    Sigma.mat.inv <- solve(Sigma.mat)

    tan.top <- Sigma.mat.inv %*% (mu.vec - rf * one.vec)
    tan.bot <- t(one.vec) %*% Sigma.mat.inv %*% (mu.vec - rf * one.vec) 

    t.vec =  tan.top[,1] / tan.bot
    t.vec
}


tanPortfolio.noShort <- function(mu.vec, Sigma.mat, risk.free.rate) {
    # Parameters:
    #   mu.vec: asset expected returns
    #   Sigma.mat: asset covariance matrix
    #   risk.free.rate: risk-free rate
    #
    # Return Value:
    #   t.vec: tangency portfolio weights
    #
    # Constrained optimization, maximize sharpe ratio:
    #
    #                 mu_p - rf
    # max(x) SR_p = -------------
    #                  sigma_p
    # 
    # mu_p = x' * mu
    # 
    # sigma_p^2 = x' * Sigma * x
    # 
    # s.t:
    #       x' * 1 = 1
    #
    #       x_i >= 0
    #
    #         
    #       [ A_eq'  ]
    # A' =  [ A_neq' ]
    #
    #       [ b_eq  ]
    # b  =  [ b_neq ]
    #
    #                                  [ mu_1-rf  1  0  0 ]
    # A' = [ (mu_p - rf)' ]        A = [ mu_2-rf  0  1  0 ]
    #      [ I_n          ]            [ mu_n-rf  0  0  1 ]
    #
    #      
    # b  = [   1    ]       # equality constraint       (x'*1  = 1)
    #      [   0    ]       # inequality constraints    (x_i >= 0 )
    #
    # I don't understand this solution...
    #

    N <- length(mu.vec)
    D.mat <- 2 * Sigma.mat
    d.vec <- rep(0, N)

    er.excess <- mu.vec - risk.free.rate
    A.mat <- cbind(er.excess, diag(1,N))
    b.vec <- c(1, rep(0,N))

    qp.out <- solve.QP(Dmat=D.mat,
                       dvec=d.vec,
                       Amat=A.mat,
                       bvec=b.vec,
                       meq=1)
    t.vec <- qp.out$solution / sum(qp.out$solution)
    names(t.vec) <- names(mu.vec)
    t.vec
}




effPortfolio <- function(mu.vec, Sigma.mat, target.return) {
    # Parameters:
    #   mu.vec: asset expected returns
    #   Sigma.mat: asset covariance matrix
    #   target.return: target expected return of the efficient portfolio
    #
    # Return Value:
    #   x.vec: asset weights for the efficient portfolio
    #
    # Constrained optimization: minimize portfolio variance subject to constraints:
    #   1) asset weights must sum to 1 
    #   2) portfolio expected return must equal target return.
    #
    # min(xA,xB,xC) sigma2_p,x = x' * Sigma * x
    #
    #           s.t:
    # 
    #       mu_p,x = x' * mu = mu_p,t = target return
    #
    #       x'*1 = 1
    #
    # L(x,lambda1, lambda2) = x'*Sigma*x + 
    #                         lambda1 * [ x'*mu - mu_p,t ] +
    #                         lambda2 * [ x'*1 - 1 ]
    #
    # Take partial derivatives wrt x, lambda1, lambda2. 
    # End up with 5 equations, 5 unknowns.
    # 
    # [ 2*Sigma  mu  1 ]   [    x    ]    [   0    ]
    # [ mu'      0   0 ] * [ lambda1 ] =  [ mu_p,t ]
    # [ 1'       0   0 ]   [ lambda2 ]    [   1    ]
    # 
    #         Ax         *     zx      =      b0
    # 
    # zx = Ax^1 * b0
    #
    # 
    # TODO: PERF: separately solve for A.mat.inv and allow user to pass in

    N <- length(mu.vec)

    A.mat.top <- cbind( 2*Sigma.mat, mu.vec, rep(1,N) )
    A.mat.mid <- c( mu.vec, 0, 0 )
    A.mat.bot <- c( rep(1,N), 0, 0 )

    A.mat <- rbind(A.mat.top, 
                   A.mat.mid,
                   A.mat.bot)
    A.mat.inv <- solve(A.mat)

    b0.vec <- c(rep(0,N), target.return, 1)

    z.vec <- A.mat.inv %*% b0.vec
    x.vec <- z.vec[1:N,1]
    x.vec
}



effPortfolio.noShort <- function(mu.vec, Sigma.mat, target.return) {
    # Parameters:
    #   mu.vec: asset expected returns
    #   Sigma.mat: asset covariance matrix
    #   target.return: target expected return of the efficient portfolio
    #
    # Return Value:
    #   x.vec: asset weights for the efficient portfolio, no shorting allowed
    #
    # Quadratic Programming problem with inequality constraints.
    # Numerical solutions only! 
    # Makes use of solve.QP
    #
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
    #      [ mu_p,0 ]       # first equality constraint (x'*mu = mu_p,0)
    # b  = [   1    ]       # 2nd equality constraint   (x'*1  = 1)
    #      [   0    ]       # inequality constraints    (x_i >= 0 )
    #
    # D = 2*Sigma
    #
    # d = (0,...,0)'
    #

    N <- length(mu.vec)
    D.mat <- 2 * Sigma.mat
    d.vec = rep(0,N)

    A.mat.t <- rbind( t(mu.vec),
                      rep(1,N),
                      diag(N) )
    A.mat <- t(A.mat.t)

    b.vec <- c(target.return, 1, rep(0,N))

    qp.out <- solve.QP(Dmat=D.mat,
                       dvec=d.vec,
                       Amat=A.mat,
                       bvec=b.vec,
                       meq=2)       # 2 equality constraints
    x.vec <- qp.out$solution
    names(x.vec) <- names(mu.vec)
    x.vec
}




effFrontier <- function(mu.vec, 
                        Sigma.mat, 
                        alpha.vec=seq(from=-1,to=2,by=0.1),
                        gmv.portfolio.weights=NULL,
                        eff.portfolio.weights=NULL) {
    # Parameters:
    #   mu.vec: asset expected returns
    #   Sigma.mat: asset covariance matrix
    #   alpha.vec: sequence of alpha values (eff frontier = convex combo of two frontier portfolios)
    #   gmv.portfolio.weights: weights for the global min var portfolio. can be null
    #   eff.portfolio.weights: weights for another eff portfolio. can be null
    #
    # Return Value:
    #   x.mat: length(alpha.vec) X 3 matrix of asset weights for the efficient portfolio
    #
    # All efficient frontier portfolios can be represented as a convex ("sum-to-one") 
    # combination of two other eff frontier portfolios.
    #

    if (is.null(gmv.portfolio.weights)) {
        gmv.portfolio.weights = gmvPortfolio(mu.vec, Sigma.mat)
    }

    if (is.null(eff.portfolio.weights)) {
        eff.portfolio.weights = effPortfolio(mu.vec, Sigma.mat, target.return=max(mu.vec))
    }

    x.mat = sapply( as.array(alpha.vec),
                    function(alpha) { alpha * gmv.portfolio.weights + (1 - alpha) * eff.portfolio.weights }
                  )
    x.mat <- t(x.mat)
    x.mat

    #
    # Computing means and sigmas:
    #
    # eff.frontier.means <- apply(x.mat,
    #                             1,
    #                             function(x.vec) { t(x.vec) %*% mu.vec })
    #
    # eff.frontier.sigmas <- apply(x.mat,
    #                              1,
    #                              function(x.vec) { sqrt(t(x.vec) %*% Sigma.mat %*% x.vec) })
    #

}



effFrontier.noShort <- function(mu.vec, Sigma.mat, target.return.vec) {
    # Parameters:
    #   mu.vec: asset expected returns
    #   Sigma.mat: asset covariance matrix
    #   target.return.vec: sequence of target returns. Typically between gmv portfolio return and max individual asset return
    #
    # Return Value:
    #   x.mat: length(alpha.vec) X 3 matrix of asset weights for the efficient portfolios, no shorting allowed
    #
    # Quadratic Programming problem with inequality constraints.
    # Numerical solutions only! 
    # Makes use of solve.QP

    x.mat = sapply( as.array(target.return.vec),
                    function(target.return) { effPortfolio.noShort(mu.vec, Sigma.mat, target.return) }
                  )
    x.mat <- t(x.mat)
    x.mat
}




portfolioWeightsBarplot <- function(x.mat, ...) {
    # Parameters:
    #   x.mat: N x m matrix of portfolio weights
    #          N = number of portfolios
    #          m = number of weights (assets) per portfolio
    # 
    # Generates barplot of portfolio weights
    #

    test1 <- test2 <- t(x.mat)
    test1[test1<0] <- 0
    test2[test2>0] <- 0
    myrange <- c(min(colSums(test2)),max(colSums(test1)))
    barplot(test1,ylim=myrange,legend.text=T,...)
    barplot(test2,add=TRUE,ylim=rev(myrange),...)
}



my.portfolio.test <- function() {


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

    #
    # Compute simple returns, means, sd, cov
    # Portfolio theory assumes simple returns (as opposed to cc returns)
    # 
    all_prices <- merge(IBM_prices, MSFT_prices, SBUX_prices)

    # diff: computes pt1 - pt0
    # lag: shifts all_prices by k=-1 (so that pt-1 -> pt)
    simple_returns <- diff(all_prices) / lag(all_prices,k=-1)
    simple_returns <- coredata(simple_returns)

    asset_names <- c("IBM", "MSFT", "SBUX")
    colnames(simple_returns) <- asset_names

    mu.vec <- apply(simple_returns, 2, mean)
    sigma.vec <- apply(simple_returns, 2, sd)
    sigma.mat <- cov(simple_returns)


    par(mfrow=c(3,1))

    #
    # Mean-Var Plot: Plot the assets 
    #
    plot(x=sigma.vec,
         y=mu.vec,
         pch=16, 
         ylim=c(0, max(mu.vec) * 1.5), 
         xlim=c(0, max(sigma.vec) *1.5), 
         xlab=expression(sigma[p]), 
         ylab=expression(mu[p]))

    for (i in seq_along(asset_names)) {
        text(x=sigma.vec[i], 
             y=mu.vec[i], 
             labels=asset_names[i], 
             pos=4)
    }
    # ..OR..:
    #
    # text(x=sigma.vec,
    #      y=mu.vec,
    #      labels=names(sigma.vec),
    #      pos=4)

    #
    # Compute global min var portfolio and plot it.
    #
    x_p.gmv <- gmvPortfolio(mu.vec=mu.vec,
                            Sigma.mat=sigma.mat)
    mu_p.gmv <- t(x_p.gmv) %*% mu.vec
    sigma_p.gmv <- sqrt( t(x_p.gmv) %*% sigma.mat %*% x_p.gmv )

    points(x=sigma_p.gmv,
           y=mu_p.gmv,
           pch=16,
           col="orange")
    text(x=sigma_p.gmv,
         y=mu_p.gmv,
         labels="GMV",
         pos=2)


    #
    # Compute efficient portfolio w/ same return as individual asset returns, 
    # and plot them
    #
    for (i in seq_along(asset_names)) {

        x_p.eff <- effPortfolio(mu.vec=mu.vec, 
                                 Sigma.mat=sigma.mat, 
                                 target.return=mu.vec[i])
        mu_p.eff <- t(x_p.eff) %*% mu.vec
        sigma_p.eff <- sqrt( t(x_p.eff) %*% sigma.mat %*% x_p.eff)

        # x_p.eff
        # mu_p.eff
        # sigma_p.eff

        points(x=sigma_p.eff,
               y=mu_p.eff,
               pch=16,
               col="green")
        text(x=sigma_p.eff,
             y=mu_p.eff,
             labels=paste("EFF-",asset_names[i]),
             pos=2)
    }


    #
    # Now we want to plot the efficient frontier
    #
    eff.frontier.weights <- effFrontier(mu.vec=mu.vec,
                                        Sigma.mat=sigma.mat,
                                        gmv.portfolio.weights=x_p.gmv)
    eff.frontier.means <- apply(eff.frontier.weights,
                                1,
                                function(x.vec) { t(x.vec) %*% mu.vec })

    eff.frontier.sigmas <- apply(eff.frontier.weights,
                                 1,
                                 function(x.vec) { sqrt(t(x.vec) %*% sigma.mat %*% x.vec) })
    
    points(x=eff.frontier.sigmas,
           y=eff.frontier.means,
           type="b",
           pch=16,
           col="green")


    #
    # Highlight the tangency portfolio
    # 
    rf <- 0.03/12   # monthly
    x_p.tan <- tanPortfolio(mu.vec=mu.vec,
                            Sigma.mat=sigma.mat,
                            risk.free.rate=rf)
    mu_p.tan <- t(x_p.tan) %*% mu.vec
    sigma_p.tan <- sqrt( t(x_p.tan) %*% sigma.mat %*% x_p.tan )

    points(x=c(0, sigma_p.tan),
           y=c(rf, mu_p.tan),
           pch=16,
           cex=2,
           col="orange")
    text(x=c(0, sigma_p.tan),
         y=c(rf, mu_p.tan),
         labels=c(paste("rf=",rf),"TAN"),
         pos=4)


    #
    # Compute and plot tangency + risk-free portfolio combinations
    #
    tan.portfolio.alpha <- seq(from=0,to=2,by=0.1)
    tan.frontier.means <- rf + tan.portfolio.alpha * ( mu_p.tan - rf )
    tan.frontier.sigmas <- tan.portfolio.alpha * sigma_p.tan

    points(x=tan.frontier.sigmas,
           y=tan.frontier.means,
           pch=16,
           col="orange")


    #
    # NO SHORTING! Global Min Var Portfolio
    #
    x_p.gmv.ns <- gmvPortfolio.noShort(mu.vec, sigma.mat)
    mu_p.gmv.ns <- t(x_p.gmv.ns) %*% mu.vec
    sigma_p.gmv.ns <- sqrt(  t(x_p.gmv.ns) %*% sigma.mat %*% x_p.gmv.ns )

    points(x=sigma_p.gmv.ns,
           y=mu_p.gmv.ns,
           pch=16,
           col="red")
    text(x=sigma_p.gmv.ns,
         y=mu_p.gmv.ns,
         labels="GMV-NS",
         pos=4)

    #
    # NO SHORTING! Highlight the tangency portfolio
    # 
    rf <- 0.03/12   # monthly
    x_p.tan.ns <- tanPortfolio.noShort(mu.vec=mu.vec,
                                    Sigma.mat=sigma.mat,
                                    risk.free.rate=rf)
    mu_p.tan.ns <- t(x_p.tan.ns) %*% mu.vec
    sigma_p.tan.ns <- sqrt( t(x_p.tan.ns) %*% sigma.mat %*% x_p.tan.ns )

    points(x=c(0, sigma_p.tan.ns),
           y=c(rf, mu_p.tan.ns),
           pch=16,
           cex=2,
           col="red")
    text(x=c(0, sigma_p.tan.ns),
         y=c(rf, mu_p.tan.ns),
         labels=c(paste("rf=",rf),"TAN-NS"),
         pos=4)


    #
    # NO SHORTING! Eff frontier from GMV portfolio to MAX(mu.vec)
    # 
    target.return.vec <- seq(from=mu_p.gmv.ns, to=max(mu.vec), length.out=10)

    eff.frontier.ns.weights <- effFrontier.noShort( mu.vec=mu.vec,
                                                    Sigma.mat=sigma.mat,
                                                    target.return.vec=target.return.vec)

    eff.frontier.ns.means = apply(eff.frontier.ns.weights,
                                  1,
                                  function(x.vec) { t(x.vec) %*% mu.vec } )

    eff.frontier.ns.sigmas = apply(eff.frontier.ns.weights,
                                   1,
                                   function(x.vec) { sqrt( t(x.vec) %*% sigma.mat %*% x.vec ) } )
    points(x=eff.frontier.ns.sigmas,
           y=eff.frontier.ns.means,
           type="b",
           pch=16,
           col="red")


    #
    # Plot efficient frontier portfolio weights in a stacked bar plot.
    #
    portfolioWeightsBarplot( eff.frontier.weights, beside=T)


    #
    # Plot efficient tangent + risk-free portfolio weights in a stacked bar plot.
    #
    tan.frontier.weights <- t(sapply(tan.portfolio.alpha,function(alpha) { alpha * x_p.tan} ))

    # add risk-free weight
    tan.frontier.weights = cbind(1-tan.portfolio.alpha, tan.frontier.weights  )
    colnames(tan.frontier.weights)[1] = "RF"

    portfolioWeightsBarplot( tan.frontier.weights, beside=T)

    par(mfrow=c(1,1))

}



