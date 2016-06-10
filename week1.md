

# Coursera: Intro to Computational Finance and Financial Econometrics

Course page: [https://class.coursera.org/compfinance-009](https://class.coursera.org/compfinance-009)  
Eric Zivot   


# Coursera: Intro to Computational Finance and Financial Econometrics: Week 1: Asset Return Calculations


## Time Value of Money


Future Value: (compounding annually)
    
    FV_1 = PV * (1 + r)
    FV_2 = PV * (1 + r) * (1 + r)

    FV_n = PV * (1 + r)^n

    
Present Value:

    PV = FV_n / (1 + r)^n


Compounded Annual Return:

    r =  ( FV_n / PV )^(1/n) - 1


Investment Horizon:

    n = ln(FV_n/PV) / ln(1 + r)

    E.g: how long does it take for my money to double?

    2 * PV = PV * (1 + r)^n

    n = ln(2PV/PV) / ln(1+r)
      = ln(2) / ln(1+r) 
      ~ 0.7 / r                 (if r is small)


Compounding m times per year:

    FV = PV * (1 + r/m)^(m*n)


Compounding continuously:

    FV = lim_m->inf PV * (1 + r/m)^m*n
       = PV * e^(r*n)


Effective Annual Rate:

    PV * (1 + r/m)^m*n = PV * (1 + R)^n

    R = (1 + r/m)^m - 1


Effective Annual Rate, Continuous Compounding:

    PV * e^r*n = PV * (1 + R)^n

    R = e^r - 1



## Asset Return Calculations


Holding period return:

    P_0 = price today
    P_t = price at time t

                             P_t - P_0 
    holding period return =  -----------
                               P_0

                          = % delta.P

Simple Returns:


    P_t = price at end of month t
         on asset that pays no dividends
    P_t-1 = price at end of month t-1

           P_t - P_t-1 
    R_t = -------------
             P_t-1

           P_t 
    R_t = ----- - 1
          P_t-1

                P_t 
    1 + R_t = ------
               P_t-1


Multi-period returns: Simple two-month return


                  P_t 
    1 + R_t(2) = ------
                 P_t-2



Multi-period returns: Relationship of two-month return to one-month returns

               P_t  
    R_t(2) = ------ - 1 
              P_t-2    

               P_t      P_t-1
           = ------- * ------- - 1
              P_t-1     P_t-2


           = (1 + R_t) * (1 + R_t-1) - 1

           = a.k.a geometric average of the one-month returns


    1 + R_t(2) = (1 + R_t) * (1 + R_t-1)

    R_t(2) = R_t + R_t-1 + R_t * R_t-1

           = R_t + R_t-1   (if R_t and R_t-1 are small)



Multi-period returns: General case, k-month return:


               P_t  
    R_t(k) = ------ - 1 
              P_t-k    


    1 + R_t(k) = (1 + R_t) * (1 + R_t-1) * ... * (1 + R_t-k+1)

               = PRODUCT_j=0..k-1 [ 1 + R_t-j ]




Multi-period returns: per-period compounded return over k periods:

                  [    P_t    ]
    (1 + R_t)^k = [ --------- ]
                  [  P_t-k    ]


                  [    P_t    ]^(1/k)
          R_t   = [ --------- ]       - 1
                  [  P_t-k    ]




## Portfolio Returns


Simple Portfolio Return (`R_p`) is the weighted average of individual asset returns.


    R_p = w1 * R_1 + w2 * R_2 + ...

        = SUM_i=1..n [ w_i * R_i ]


## Adjusting for Dividends


    D_t = dividend payment between months t-1 and t


                 P_t + D_t - P_t-1
    total.R_t = ------------------
                      P_t-1

                 P_t - P_t-1     D_t
              = ------------- + -----
                   P_t-1         P_t-1


              = capital gain  + dividend yield


                      P_t + D_t
    1 + total.R_t   = ------------
                         P_t-1



## Adjusting for Inflation 
                   
* REAL rate of return
* vs. NOMINAL rate of return
* deflate NOMINAL price P_t 
    * by an index of the general price level
    * CPI_t

.
                 P_t
    real.P_t = -------
                CPI_t

                
                 real.P_t - real.P_t-1 
    real.R_t  = -----------------------
                      real.P_t-1         

                  P_t        CPI_t-1
              = -------  * ----------- - 1
                 P_t-1        CPI_t


Or, define inflation as:

             CPI_t - CPI_t-1
    pi_t = --------------------
                 CPI_t-1

                1 + R_t
    real.R_t = ----------  - 1
                1 + pi_t


    real.R_t ~=  R_t - pi_t   (for small R and pi)



## Annualizing Returns


Example: **assume** same monthly return, `R_m`, for 12 months:


    compound annual gross return = 1 + R_a = (1 + R_m)^12

    compound annual net return   =     R_a = (1 + R_m)^12 - 1



Annualized return over k years:


                  [    P_t    ]
    (1 + R_t)^k = [ --------- ]
                  [   P_t-k   ]


                  [    P_t    ]^(1/k)
          R_t   = [ --------- ]       - 1
                  [   P_t-k   ]



## Continuously Compounded Returns


    R_t = simple return
    r_t = continuously compounded return
    r_t is always smaller than R_t

    r_t = ln(1 + R_t) 
        
        = ln( P_t / P_t-1 )     # simple return: 1+R_t = P_t / P_t-1

        ~= R_t  (if R_t is small)


    R_t = e^r*t - 1


    e^(r_t) = e^(ln(1+R_t))

            = e^(ln(P_t/P_t-1))

            = P_t / P_t-1

    P_t-1 * e^(r_t) = P_t


    r_t = ln( P_t / P_t-1 )

        = ln(P_t) - ln(P_t-1)

        = diff in log prices


* Compounded growth rates are typically plotted on a log-scale
* slope is equal to the growth rate


Multi-period continuously compounded return:

    r_t(2) = continuously compounded growth rate between t-2 and t

    r_t(2) = ln(1 + R_t(2))

           = ln( P_t / P_t-2 )

           = ln(P_t) - ln(P_t-2)



    P_t-2 * e^(r_t(2)) = P_t


Multi-period cc return, relationship to 1-period cc return:

                 P_t     P_t-1
    r_t(2) = ln( ----  * ------ )
                 P_t-1   P_t-2 


           = ln(P_t/P_t-1) + ln(P_t-1/P_t-2)

           = r_t + r_t-1

           = sum of 1-period cc returns

           WE LIKE ADDITIVE!!

           sum of two normally distributed random variables
           is also a random variable.

           NOT TRUE under multiplication.


Mult-period cc returns over k periods: General case:


    r_t(k) = ln(1 + R_t(k))


           = ln( P_t / P_t-k )


           = r_t + r_t-1 + ... + r_t-k+1



## Natural Logs and Exponentials


    ln(0) = -inf
    ln(1) = 0

    e^-inf = 0
    e^0 = 1
    e^1 = 2.7183

    ln(e^x) = x
    e^ln(x) = x

    d.ln(x)/dx = 1/x
    d.e^x/dx = e^x

    ln(x*y) = ln(x) + ln(y)
    ln(x/y) = ln(x) - ln(y)

    ln(x^y) = y * ln(x)
    
    e^x * e^y = e^(x+y)
    e^x * e^-y = e^(x-y)

    (e^x)^y = e^(x*y)



## Continously Compounded Portfolio Returns 


    R_p = w1 * r1 + w2 * r2 + ... 
    
    r_p = ln(1 + R_p)

        = ln(1 + w1 * r1 + w2 * r2 + ... )

    Portfolio returns are NOT additive.

    r_p ~= R_p   (if R_p is small)

    R_p > r_p    (if R_p is big)


## Continously Compunded Return, Adjusting for Inflation


    
    r_t = nominal
    pi_t = continously compounded cpi inflation rate

    real.r_t = ln(1 + real.R_t)

                     P_t        CPI_t-1
    1 + real.R_t = -------  * ----------- 
                    P_t-1        CPI_t


                     P_t             CPI_t-1
    real.r_t = ln( ------- ) + ln( -----------  )
                    P_t-1             CPI_t

             = ln(P_t) - ln(P_t-1) + ln(CPI_t-1) - ln(CPI_t)

             = r_t                 - pi_t

    real cc rate is nominal cc rate minus cc inflation rate

    WE LIKE ADDITIVE!



## Growth of $1 Graphs

* the thing about returns is...
    * if you lose 50% in one period
    * you need a 100% return in the next period to get back to whole
    * a +50% move does not cancel out a -50% move in the previous period
* so to better understand return performance...
    * people often use "Growth of $1" graphs
* given a series of daily (or monthly) returns
* compute the performance of $1 invest
    * after period 1: (1 + R_1)
    * after period 2: (1 + R_1) * (1 + R_2)
    * ...
    * after period n: (1 + R_1) * (1 + R_2) * ... * (1 + R_n)
