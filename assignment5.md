
# Coursera: Intro to Computational Finance: Assignment 5: Linear Combinations of RVs, Portfolios, Time-Dependence Models 


1\. Let: 

* E[X]=2
* E[Y]=1
* Var[X]=3
* Var[Y]=2.5
* Cov(X,Y)=.9

What is E[.4X+.6Y]?

    # Linear Combination of two random variables
    # Z = aX + bY
    # E[Z] = a * E[X] + b * E[Y]

    E_x = 2
    E_y = 1
    E_z = 0.4 * E_x + 0.6 * E_y
    # [1] 1.4


2\. What is Var[.4X+.6Y]

    # Var(Z) = Var(aX + bY)
    #        = a^2 * Var(X) + b^2 * Var(Y) + 2*a*b*Cov(X,Y)
    a = 0.4
    b = 0.6
    sigma_xy = 0.9
    var_x = 3
    var_y = 2.5

    var_z = a^2 * var_x + b^2 * var_y + 2 * a * b * sigma_xy
    # [1] 1.812

3\. Suppose X and Y are returns on two assets, and w is your portfolio weight
in asset X with (1−w) being the portfolio weight in asset Y. What value of w
minimizes the variance of your portfolio?

[http://faculty.washington.edu/ezivot/econ424/portfolioTheoryMatrix.pdf](http://faculty.washington.edu/ezivot/econ424/portfolioTheoryMatrix.pdf)

    # To solve for global min variance:
    # Sigma = [ sigma_x^2 sigma_xy  ]
    #         [ sigma_xy  sigma_y^2 ]

    # Sigma = matrix( c(var_x, sigma_xy, sigma_xy, var_y), 
    #                 nrow=2,
    #                 ncol=2,
    #                 byrow=T )
    Sigma = rbind( c(var_x, sigma_xy),
                   c(sigma_xy, var_y) )

    # W.mat = [ w1 ]
    #         [ w2 ]

    # A.mat = [ 2*Sigma  1 ]
    #         [ 1'       0 ]
    A.mat.top = cbind( 2*Sigma, c(1,1) )
    A.mat = rbind( A.mat.top, c(1,1,0) )

    # b = [ 0 ]
    #     [ 0 ]
    #     [ 1 ]
    b = c(0,0,1)

    # z = [ W.mat  ]
    #     [ lambda ]

    # A.mat %*% z = b
    # z = A.mat^-1 %*% b

    z = solve(A.mat) %*% b
    #            [,1]
    # [1,]  0.4324324     <---- w
    # [2,]  0.5675676     <---- 1-w
    # [3,] -3.6162162
    

4\. What is the variance of the portfolio with the weights you derived in the
previous question?
    
    w = 0.4324324
    
    # var_p = w^2 * var_x + (1 - w)^2 * var_y + 2 * w * (1 - w) * sigma_xy
    x = c(w, 1-w)

    var_p = t(x) %*% Sigma %*% x
    #          [,1]
    # [1,] 1.808108


5\. What is the expected value of the portfolio with the weights you derived in
the previous question?

    R.mat = c(E_x, E_y) 

    R_p = t(x) %*% R.mat
    #          [,1]
    # [1,] 1.432432


6\. For the following questions, consider the AR(1) model:

    Y_t = 10 + .6 * Y_t−1 + e_t 

    e_t ~ N(0,2^2)

Is the process is covariance stationary? YES.

7\. What is the mean of this process? 

    Mean-adjusted:
    Yt - mu = phi * (Yt-1 - mu) + e_t
         Yt = mu + phi * Y_t-1 - phi * mu = 10 + phi * Y_t-1
              mu - phi * mu = 10
              mu = 10 / (1 - phi)
              mu = 10 / 0.4
              mu = 25 

8\. What is the variance of the process?

    # Var(Yt) = sigma^2 = sigma_e^2 / (1 - phi^2)
    sigma_e = 2
    phi = 0.6

    var_yt = sigma_e^2 / (1-phi^2)
    # [1] 6.25


