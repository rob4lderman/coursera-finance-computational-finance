

# Coursera: Intro to Computational Finance: Assignment #4: Matrix Algebra


Questions 1 - 9 refer to the following matrices and vectors: 

        [ 1 4 7 ]       [ 4 4 0 ]       [ 1 ]       [ 5 ]
    A=  [ 2 4 8 ]  B =  [ 5 9 1 ]   x = [ 2 ]   y = [ 2 ]
        [ 6 1 3 ]       [ 2 2 5 ]       [ 3 ]       [ 7 ]


1\. Compute the transpose of A.

    A = matrix( c(1,4,7,2,4,8,6,1,3), nrow=3, ncol=3, byrow=T )
    t(A)
    #      [,1] [,2] [,3]
    # [1,]    1    2    6
    # [2,]    4    4    1
    # [3,]    7    8    3


2\. Compute the transpose of B.

    B = matrix( c(4,4,0,5,9,1,2,2,5), nrow=3, ncol=3, byrow=T )
    t(B)
    #      [,1] [,2] [,3]
    # [1,]    4    5    2
    # [2,]    4    9    2
    # [3,]    0    1    5


3\. Compute the transpose of x.

    x = c(1,2,3)
    t(x)
    #      [,1] [,2] [,3]
    # [1,]    1    2    3


4\. Compute the transpose of y.

    y = c(5,2,7)
    t(y)
    #      [,1] [,2] [,3]
    # [1,]    5    2    7

5\. Compute A + B

    A + B
    #      [,1] [,2] [,3]
    # [1,]    5    8    7
    # [2,]    7   13    9
    # [3,]    8    3    8

6\. Compute A - B

    A - B
    #     [,1] [,2] [,3]
    # [1,]   -3    0    7
    # [2,]   -3   -5    7
    # [3,]    4   -1   -2

7\. Compute 2 * A

    2*A
    #      [,1] [,2] [,3]
    # [1,]    2    8   14
    # [2,]    4    8   16
    # [3,]   12    2    6

8\. Compute A x

    A %*% x
    #      [,1]
    # [1,]   30
    # [2,]   34
    # [3,]   17

9\. compute y' A x

    t(y) %*% A %*% x
    #      [,1]
    # [1,]  337

10\. Consider the system of equations: 

    z1 + z2 = 1

    2z1 + 4z2 = 2

Write the system using matrix notation as Az=b and solve for z.

    A = [ 1  1 ]
        [ 2  4 ]

    z = [ z1 ]
        [ z2 ]

    b = [ 1 ]
        [ 2 ]
        
    A = matrix(c(1,1,2,4), nrow=2, ncol=2, byrow=T)
    b = c(1,2)

    z = solve(A) %*% b
    #      [,1]
    # [1,]    1
    # [2,]    0


11\. Consider creating a portfolio of three assets denoted A, B and C. Assume the following information: 

            [ 0.01 ]            [ 0.10  0.30  0.10  ]
    mu =    [ 0.04 ]    Sigma = [ 0.30  0.15  -0.20 ]
            [ 0.02 ]            [ 0.10  -0.20 0.08  ]

Compute the expected return for an equally weighted portfolio (i.e., xA=xB=xC=1/3).

    x = c(1/3, 1/3, 1/3)
    mu = c(0.01, 0.04, 0.02)

    er_p = t(x) %*% mu
    #            [,1]
    # [1,] 0.02333333


12\. Continuing from the previous question, what is the variance for an equally weighted portfolio?

    Sigma = matrix(c(0.10, 0.30, 0.10,
                     0.30, 0.15, -0.20,
                     0.10, -0.20, 0.08), 
                   nrow=3, 
                   ncol=3,
                   byrow=T)

    var_p = t(x) %*% Sigma %*% x
    #            [,1]
    # [1,] 0.08111111



