
# Coursera: Intro to Computational Finance: Assignment 3: BiVariate Computations


1\.  Consider the following joint distribution of X and Y: 

    X/Y	1	2	3
    1	0.1	0.2	0
    2	0.1	0	0.2
    3	0	0.1	0.3

Find the marginal distributions of X and Y.
Using the table above, compute E[X]. 


    # p(x) = SUM_all_y p(x,y)
    # p(x=1) = SUM [ p(x=1,y=1) + p(x=1,y=2) + p(x=1,y=3)

    joint_pdf = matrix(c(0.1,0.1,0, 0.2,0,0.1, 0,0.2,0.3), nrow=3, ncol=3)
    pdf_x = c( sum(joint_pdf[1,]), sum(joint_pdf[2,]), sum(joint_pdf[3,]) )
    pdf_y = c( sum(joint_pdf[,1]), sum(joint_pdf[,2]), sum(joint_pdf[,3]) )

    ev = function(pdf) {
        sum = 0
        for (i in seq_along(pdf)) {
            sum = sum + pdf[i] * i
        }
        sum
    }

    var1 = function(pdf) {
        mu = ev(pdf)
        sum = 0
        for (i in seq_along(pdf)) {
            sum = sum + (i - mu)^2 * pdf[i]
        }
        sum
    }


    e_x = ev(pdf_x) 
    # [1] 2.1

2\. Compute E[Y]:

    e_y = ev(pdf_y)
    # [1] 2.3

3\. Compute Var(X): 

    var_x = var1(pdf_x)
    # [1] 0.69

4\. Compute Var(Y): 

    var_y = var1(pdf_y)
    # [1] 0.61

5\. Compute stdev(X): 

    sd_x = sqrt(var_x)
    # [1] 0.8306624

6\. Compute stdev(Y): 

    sd_y = sqrt(var_y)
    # [1] 0.781025

7\. Compute Cov(X,Y):

    # Cov(X,Y) = sigma_xy = E[(X-mu_x)(Y-mu_y)]
    #          = SUM_s_xy (x - mu_x) * (y - mu) * p(x,y)

    covar1 = function(joint_pdf, mu_x, mu_y) {
        sum = 0
        for (x in 1:nrow(joint_pdf)) {
            for (y in 1:ncol(joint_pdf)) {
                sum = sum + (x - mu_x) * (y - mu_y) * joint_pdf[x,y]
            }
        }
        sum
    }

    var_xy = covar1(joint_pdf, e_x, e_y)
    # [1] 0.37

8\. Compute Cor(X,Y):

    rho_xy = var_xy / sd_x / sd_y
    # [1] 0.5703117

9\. Are X and Y independent?  NO.


10\. Let r denote the continuously compounded monthly return on Microsoft stock
and let W0 denote initial wealth to be invested over the month. Assume that r∼
iid N(0.04,(0.09)2) and that W0 = $100,000.

Determine the 1% and 5% value-at-risk (VaR) over the year on the
investment. Hint: to answer this question, you must determine the normal
distribution that applies to the annual (12 month) continuously compounded
return. This was done as an example in class.

    # r ~ N(0.04, 0.09^2)
    mu_r = 0.04
    sigma_r = 0.09

    # r_a = 12 * r
    mu_r_a = 12 * mu_r
    sigma_r_a = sqrt(12) * sigma_r

    W0 = 1e5

    # convert cc return quantile to simple return quantile
    # R = e^r - 1
    # q_.01(R) = e^(q_.01(r)) - 1
    #
    # V-a-R(1%) = (exp( q_0.01(r_a) ) - 1) * W0
    var_01 = W0 * ( exp( qnorm(0.01, mean=mu_r_a, sd=sigma_r_a) ) -1 ) 
    # [1] -21751.73
    
    # V-a-R(5%) = q_0.05(r_a) * W0
    var_05 = W0 * ( exp( qnorm(0.05, mean=mu_r_a, sd=sigma_r_a) ) -1 )
    # [1] -3228.205

