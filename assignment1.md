
# Coursera: Intro to Computational Finance: Assignment 1: Computing Asset Returns

Consider the following (actual) monthly adjusted closing price data for
Starbucks stock over the period December 2004 through December 2005 


    End of Month Price Data for Starbucks Stock

    December, 2004	$31.18
    January, 2005	$27.00
    February, 2005	$25.91
    March, 2005	    $25.83
    April, 2005	    $24.76
    May, 2005	    $27.40
    June, 2005	    $25.83
    July, 2005	    $26.27
    August, 2005	$24.51
    September, 2005	$25.05
    October, 2005	$28.28
    November, 2005	$30.45
    December, 2005	$30.51


1\. Using the data in the table, what is the simple monthly return between the
end of December 2004 and the end of January 2005?

    # simple return: R = P_t - P_t-1 / P_t-1
    pt0 = 31.18
    pt1 = 27.00
    R = (pt1 - pt0) / pt0
    # [1] -0.1340603

2\. If you invested $10,000 in Starbucks at the end of December 2004, how much
would the investment be worth at the end of January 2005?

    # FV = PV * (1 + R)
    PV = 10000
    R = -0.1340603
    FV = PV * (1 + R)
    # [1] 8659.397
   
3\. Using the data in the table, what is the continuously compounded monthly
return between December 2004 and January 2005?

    # continously compounded return: r_t = ln(1 + R_t) 
    R = -0.1340603
    r = log(1 + R)
    # [1] -0.14394

4\. Assuming that the simple monthly return you computed in Question 1 is the
same for 12 months, what is the annual return with monthly compounding?

    # compound annual gross return = 1 + R_a = (1 + R_m)^12
    # compound annual net return   =     R_a = (1 + R_m)^12 - 1
    R_m = -0.1340603
    R_a = (1 + R_m)^12 - 1
    # [1] -0.8222327

5\. Assuming that the continuously compounded monthly return you computed in
Question 3 is the same for 12 months, what is the continuously compounded
annual return? 

    # multi-period continously compounding returns are simply added together
    r_m = -0.14394
    r_a = r_m * 12
    # [1] -1.72728

6\. Using the data in the table, compute the actual simple annual return
between December 2004 and December 2005.

    # simple return: P_t - P_t-1 / P_t-1
    # December, 2004	$31.18
    # December, 2005	$30.51

    pt0 = 31.18
    pt1 = 30.51
    R = (pt1 - pt0) / pt0
    # [1] -0.02148813
    
7\. If you invested $10,000 in Starbucks at the end of December 2004, how much
would the investment be worth at the end of December 2005?

    # FV = PV * (1 + R)
    PV = 10000
    R = -0.02148813
    FV = PV * (1 + R)
    # [1] 9785.119

8\. Using the data in the table, compute the actual annual continuously
compounded return between December 2004 and December 2005.

    # continously compounded return: r_t = ln(1 + R_t) 
    R = -0.02148813
    r = log(1 + R)
    # [1] -0.02172236


