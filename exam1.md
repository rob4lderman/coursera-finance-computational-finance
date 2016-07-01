

# Coursera: Intro to Computational Finance: Mid-term Exam


### Questions 1 - 14 depend on the following information:

Consider a 1-month investment in two assets: 

* the Vanguard S&P 500 index (VFINX)
* the Vanguard Emerging Markets Stock Index (VEIEX). 
    
Suppose you buy one share of the S&P 500 fund and one share of the emerging
markets fund at the end of September, 2010 

    PVFINX_t-1=105.06
    PVEIEX_t-1=28.64, 

and then sell these shares at the end of October, 2010 

    PVFINX_t=109.04, 
    PVEIEX_t=29.49


### 1. What is the simple 1-month return for VFINX?

    # R = P_t/P_t-1 - 1

    P0 = 105.06
    P1 = 109.04

    P1/P0 - 1
    # [1] 0.03788311

### 2. What is the simple 1-month return for VEIEX?

    P0 = 28.64
    P1 = 29.49

    P1/P0 - 1
    # [1] 0.02967877

### 3. What is the continuously compounded (cc) 1-month return for VFINX?

    # r = ln(1 + R)

    log( 1 + 0.03788311)
    # 0.03718317

### 4. What is the continuously compounded (cc) 1-month return for VEIEX?


    log( 1 + 0.02967877)
    # [1] 0.02924688

### 5. Assume that you get the same monthly return for VFINX from Q1 every month for the next year. What is the annualized simple return? 


    # 1 + R_a = (1 + R_m)^12

    (1 + 0.03788311)^12 
    # [1] 1.562361

### 6. Assume that you get the same monthly return for VEIEX from Q2 every month for the next year. What is the annualized simple return? 


    ( 1 + 0.02967877)^12
    # [1] 1.420434


### 7. Assume that you get the same monthly return for VFINX from Q1 every month for the next year. What is the annualized cc return ? 


    # r_a = 12 * r_m
    # r_m = log(1 + R_m)

    log(1.562361)
    # [1] 0.4461981

    12 * log( 1 + 0.03788311)
    # [1] 0.4461981

### 8. Assume that you get the same monthly return for VEIEX from Q2 every month for the next year. What is the annualized cc return?

    log(1.420434)
    # [1] 0.3509625

### 9. Assume that you get the same monthly return for VFINX from Q1 every month for the next year. Approximately how much will $10,000 invested in VFINX be worth after 1 year? 

    # FV = PV * (1 + R)^12

    PV=10000

    PV * (1 + 0.03788311)^12
    # [1] 15623.61

### 10. Assume that you get the same monthly return for VEIEX from Q2 every month for the next year. How much will $10,000 invested in VEIEX be worth after 1 year? 


    PV=10000
    PV * (1 + 0.02967877)^12
    # [1] 14204.34


At the end of September, you have $10,000 to invest in VFINX and VEIEX over the
next month. Suppose you purchase $2,000 worth of VFINX and the remainder in
VEIEX.

### 11. What is the portfolio weight in VFINX?
    
    20%

### 12. What is the portfolio weight in VEIEX?

    80%


### 13. What is the 1-month simple portfolio return?


    w_vfinx = 0.2
    R_vfinx = 0.03788311

    w_veiex = 0.8
    R_veiex = 0.02967877

    R_p = w_vfinx * R_vfinx + w_veiex * R_veiex
    # [1] 0.03131964

### 14. What is the 1-month continuously compounded portfolio return?

    log(1 + R_p)
    # [1] 0.03083918



### For questions 15 - 20:

Let `r_VFINX` and `r_VEIEX` denote the monthly continuously compounded returns on VFINX and VEIEX and suppose that 

    r_VFINX ~ iid N(0.001,(0.05)^2),
    r_VEIEX ~ iid N(0.01,(0.09)^2).


### 17. Let W0=$100,000 be the initial wealth invested in each asset. Compute the 5% monthly Value-at-Risk for VFINX. 

    # VaR_.05 = W0 * q_.05(R)
    # q_.05(R) = exp(q_.05(r)) - 1
    # q_.05_r = mu + sigma * q_.05(z)

    mu_vfinx = 0.001
    sigma_vfinx = 0.05

    W0 = 100000
    r_05 = mu_vfinx + sigma_vfinx * qnorm(0.05)
    R_05 = exp(r_05) - 1

    W0 * R_05
    # [1] -7803.008

### 18. Let W0=$100,000 be the initial wealth invested in each asset. Compute the 5% monthly Value-at-Risk for VEIEX.

    mu_veiex = 0.01
    sigma_veiex = 0.09

    W0 = 100000
    r_05 = mu_veiex + sigma_veiex * qnorm(0.05)
    R_05 = exp(r_05) - 1

    W0 * R_05
    # [1] -12893.34

### 19. Let W0 = $100,000 be the initial wealth invested in each asset. Compute the 5% annual Value-at-Risk for VFINX.

    # annualized cc returns:
    # mu_a = 12 * mu_m
    # sigma_a = sqrt(12) * sigma_m
    
    mu_vfinx = 0.001 * 12
    sigma_vfinx = 0.05 * sqrt(12)

    W0 = 100000
    r_05 = mu_vfinx + sigma_vfinx * qnorm(0.05)
    R_05 = exp(r_05) - 1

    W0 * R_05
    # [1] -23882.88

### 20. Let W0 = $100,000 be the initial wealth invested in each asset. Compute the 5% annual Value-at-Risk for VEIEX.

    mu_veiex = 0.01 * 12
    sigma_veiex = 0.09 * sqrt(12)

    W0 = 100000
    r_05 = mu_veiex + sigma_veiex * qnorm(0.05)
    R_05 = exp(r_05) - 1

    W0 * R_05
    # [1] -32484.61

### 21. Let {Yt} represent a stochastic process. Under which of the following conditions is {Yt} covariance stationary? Check all that apply.

    E[Y_t] = mu_t               : NO - mu depends on t
    E[Y_t] = mu                 : YES - mu constant over time
    var(Y_t) = sigma_t^2        : NO - sigma_t^2 depends on t
    var(Y_t) = sigma^2          : YES - sigma_2 is constant over time
    cov(Y_t,Y_t−j) = gamma_t    : NO - cov depends on t
    cov(Y_t,Y_t−j) = gamma_j    : YES - cov depends only on the lag j


### For the following questions, consider the AR(1) model:

    Y_t = 3 + .45 * Y_t-1 + err_t 

        err_t ~ N (0,1.5^2)

### 26. The process is covariance stationary?

    Yes.  Phi = 0.45,  -1 < phi < 1: covariant stationary

### 27. What is the mean of the process?

    # Mean adjusted:
    # Y_t - mu = phi * (Yt-1 - mu) + e_t
    #
    # Y_t = mu + phi * Y_t-1 - mu * phi + e_t
    #
    #     = mu * (1 - phi) + phi * Y_t-1 + e_t
    #
    #   3 = mu * (1 - phi)
    #
    # mu = 3 / (1 - phi)

    3/(1-0.45)
    # [1] 5.454545

### 28. What is the variance of the process?


    # Var(Y_t) = sigma_e^2 / (1 - phi^2)

    1.5^2 / (1 - 0.45^2)
    # [1] 2.821317

### 29. What is the covariance of (Yt,Yt-1)?

    # Cov(Yt, Yt-1) = sigma^2 * phi

    2.821317 * 0.45
    # [1] 1.269593

### 30.What is the correlation of (Yt,Yt-1)? 

    # Cor(Yt, Yt-1) = phi
    0.45





