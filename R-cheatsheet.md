

# R Cheatsheet


Main website: [http://cran.r-project.org/](http://cran.r-project.org/)


* [Basics](#basics)
* [Vectors and Matrices and Data Frames](#vectors)
* [Built-in commands / variables](#builtin)
* [Reading data](#readingdata)
* [Extract/Transform/Load (ETL) data](#etl)
* [Exploring data](#exploringdata)
* [Plotting: Base](#plotbase)
* [Plotting: Lattice](#plotlattice)
* [Plotting: ggplot2](#plotggplot2)
* [Plotting: devices](#plotdevices)
* [Knitr](#knitr)



## <a name="basics"></a>Basics

    x <- 1
    x = 1
    msg <- "hello"
    msg = "hello"

    print(x)
    # [1] 1       # x is a vector with one element

    x <- 1:20   # range
    x <- 2L
    is.na(x)            # missing value
    is.nan(x)           # not-a-number
    sum( is.na(x) )     # num of missing values



Class / Cast / Coersion / Types

    attributes(x)
    class(x)
    unclass(x)
    as.numeric(x)  
    as.logical(x)
    as.character(x)


Loops:

    # 1:length(<vector>)
    for (i in seq_along(<vector>)) {
    }

    # 1:len
    for (i in seq_len(<len>)) {
    }




## <a name="vectors"></a>Vectors and Matrices and Data Frames


* Vectors:
    * All elements are coerced to the same "least-common-denominator" type
* Lists:
    * Like vectors, but all elements may be different types
    * Elements are indexed by double-brackets [[1]]
* Matrices
    * Vectors with a dimension attribute
* Data Frame
    * For holding tabular data
    * Each element in the data frame is a list 
        * List of lists
    * All lists have the same length
    * Each element (list) is a column (vector) in the table
    * Every row has a name (often just the row index)
    * Every column has a name 
* Subsetting: by single bracket, `v[1]`
    * returns an object of the same type (except for matrix)
        * i.e, subset a vector, get a vector back. 
        * subset a list, get a list back
    * can be used to select a range
        * `v[0:5]`
    * can subset on column name
        * e.g  `df["col-name"]`
        * returns a list/data-frame with a single column.
        * ncol=1 ("col-name")
        * length=1 
            * 1 entry in the list
            * the "col-name" column
            * which itself is a vector

* Subsetting: by double brackets, `[[]]`
    * select a single element from a list / data frame
    * no ranges
    * returned object may not be a list / data frame
    * e.g df[["col-name"]]
         * returns the col-name value 
            * typically a vector
         * ncol= NULL/undefined (vectors do not have columns)
         * length=6 (or whatever the length of the vector is)
    * same as subsetting by $name

* Subsetting: by dollar, `df${col-name}`
    * like `[[]]`
* Subset FAQ: [http://www.ats.ucla.edu/stat/r/faq/subset_R.htm](http://www.ats.ucla.edu/stat/r/faq/subset_R.htm)


#### Vectors:

    x <- c(0.5, 0.6)
    x <- c(T,F)     
    x <- c(2, T)            # will coerce vector to common type (numeric in this case)
    x <- list(1, "a", T)
    x <- factor( c("yes", "no", "yes") , levels=c("yes", "no"))
    x <- numeric(5)     # x is a vector of length 5
    x <- factor( c("yes", "no", "yes") )

    x <- c("a", "b", "c", "c", "a")     # c = "combine"
    x[1]
    # [1] "a"

    x[1:4]
    # [1] "a" "b" "c" "c"

    x[ x > "a" ]
    # [1] "b" "c" "c"

    u <- x > "a"
    u
    # [1] FALSE TRUE TRUE TRUE FALSE

    x[u]
    # [1] "b" "c" "c"
    
    x <- list(foo=1:4, bar=0.6)
    x
    # $foo
    # [1] 1 2 3 4
    #  
    # $bar
    # [1] 0.6
    
    x[1]
    # $foo
    # [1] 1 2 3 4

    x[[1]]
    # [1] 1 2 3 4

    x$bar
    # [1] 0.6
    
    x[["bar"]]
    # [1] 0.6

    x["bar"]
    # $bar
    # [1] 0.6

    names(x) <- c("foo", "bar", "norf")
    x <- list(a=1, b=2, c=3)


#### Matrices: 

    x <- matrix( nrow=2, ncol=3 )
    x <- matrix(1:6, nrow=2, ncol=3 )           # cols are filled first
    x <- matrix(1:6, nrow=2, ncol=3, byrow=T )  # rows are filled first

    cbind(x, y)     # bind vectors x and y as columns of new matrix
    rbind(x, y)     # bind vectors x and y as rows of new matrix
    table(x)

    dim(x)   
    dimnames(m) <- list( c("a", "b"), c("c", "d") )     # for matrix

    # convert vector to col matrix
    y <- c(1,2,3)
    x = as.matrix(y)            # row=3, col=1, column vector matrix

    
    t(x)                        # transpose
    A %*% B                     # matrix multiplication
    solve(A)                    # matrix inverse AX = b (b = I)

    t(x) %*% y                  # cross product (x1 * y1 + x2 * y2 + ...)
    crossprod(x,y)

    x <- matrix(1:6, 2, 3)
    #      [,1] [,2] [,3]
    # [1,]    1    3    5
    # [2,]    2    4    6

    x[1, 2]
    # [1] 3

    x[1, ]
    # [1] 1 3 5

    x[1, , drop=FALSE]
    #      [,1] [,2] [,3]
    # [1,]    1    3    5

    # covariance and correlation matrices
    cov.mat = var(returns.mat)
    cor.mat = cor(returns.mat)

    # lower.tri 
    covhat.vals = cov.mat[ lower.tri(cov.mat) ]
    rhohat.vals = cor.mat[ lower.tri(cor.mat) ]

    diag(Sigma.mat)                 # extract diagonal
    diag(Sigma.mat) <- c(1,2,3)     # set diagonal


#### Data Frames:

    df <- read.table( ... )
    df <- read.csv( "file.csv" )
    data.matrix(df) // convert to matrix
    df <- data.frame( foo = 1:4, bar=c(T,F,T, T))

    summary(df)
    str(df)
    head(df)
    nrow(df)
    ncol(df)

    df[,1]
    df[1:2,"var2"]

    table(X$var2, useNA="ifany") // useNA to count the NA values

    sbux_prices_df <- sbux_df[, "Adj.Close", drop=FALSE]    # drop=F - preserve df dimensions
    rownames(sbux_prices_df) <- sbux_df$Date                # assign row names
    head(sbux_prices_df)
    price_1 <- sbux_prices_df["3/1/1994",1]
    price_2 <- sbux_prices_df["3/1/1995",1]

    colSums(is.na(X))
    all(colSums(is.na(X))==0)
    X$var1 %in% c("1","2")

    xtabs(var1 ~ var2 + var3,data=df)

    object.size

    library(Hmisc)






#### Data Tables: 

* type of data frame
* faster and more efficient than data frames
* [http://cran.r-project.org/web/packages/data.table/data.table.pdf](http://cran.r-project.org/web/packages/data.table/data.table.pdf)
* [http://blog.datacamp.com/data-table-cheat-sheet/](http://blog.datacamp.com/data-table-cheat-sheet/)

.

    tables()    # list all tables in memory
    DT[2,]
    DT[DT$y=="a",]
    DT[c(2,3)]                      # 2nd and 3rd rows
    DT[, list(mean(x), sum(z)) ]    # x and z are vars in the table.  reports mean of x and sum of z
    DT[, w:=z^2]                    # add new variable (named "w") to data frame
                                    # Note: does NOT create a copy of the data frame

    DT[ , m := { tmp <- (x+z); log2(tmp)} ]
    DT[ , a := x > 0 ]
    DT[ , b := mean(x+w),by=a ],    # new col b, take mean(x+w) GROUPED BY col a
    DT[ , .N, by=x]                 # count elements GROUPED BY x

    setkey()
    merge()


    install.packages("data.table")
    library(data.table)
    
    DF <- read.csv(...)
    DT <- data.table(DF)
    tables()
    
    DT[,mean(pwgtp15),by=SEX] // mean col val (pwgtp15) grouped by SEX col value
    SEX V1
    1: 1 99.80667
    2: 2 96.66534



## <a name="builtin"></a>Built-in commands / variables

    
Installing Packages:

    available.packages
    install.packages("slidify")
    library(slidify)            # not in quotes
    search()
    sessionInfo

Logic:

    any()
    all( df$var1 < 3 )


Maths:

    Inf             
    exp()           # e
    log()           # natural log (base e)
    log2()
    log10()
    logb(x, base=2)

    x %*% y         # matrix cross-product multiplication

    svd
    prcomp
    unzip


Working with Strings:

    tolower
    toupper
    strsplit
    sub         # substitution
    gsub        # global substitution
    grep
    grepl
    nchar
    substr
    paste
    paste0
    str_trim


Working with Dates:

    mutate( BGN_DATE = as.Date(BGN_DATE, format="%m/%d/%Y %H:%M:%S") ) %>%
    filter( BGN_DATE >= as.Date("2000-01-01") ) %>%
    mutate( YEAR = as.integer(strftime(BGN_DATE, format="%Y")) )
    asPOSIXlt
    asPOSIXct
    
    date() // just a string
    Sys.Date() // Date object
    format(d,"%m %d")
    as.Date(x, "<format>")

    df$Date.date = as.Date(df$Date, "%Y-%m-%d") 

    weekdays(d)
    months(d)
    julian(d)
    
    library(lubridate)
    ymd("20140108")
    mdy("..")
    dmy
    ymd_hms
    ymd_smh
    Sys.timezone
    wday(d,label=T)



apply Functions:

    # always returns a list
    lapply(f, mean)             # apply function 'mean' to every element of f

    sapply(x, mean)  
    # like lapply but "simplifies" result
    # i.e if lapply returns a list of vectors with length 1, sapply returns a single vector
    # if lapply returns a list of vectors all of the same length, sapply returns a matrix

    # apply(matrix, dim-to-retain, function)
    rowSums <- apply(x, 1, sum)
    colSums <- apply(x, 2, sum)
    rowMeans <- apply(x, 1, mean)
    colMeans <- apply(x, 2, mean)

    # mapply
    # "multi-variate" apply
    # takes multiple lists, passes elements from lists to function
    # sorta like "zip" + apply function
    x = c(5,3,0,4)
    y = c(4,4,1,3)
    mapply(function(x,y) { x+y }, x, y)


    # tapply
    # pass subsets of vector to function specified by factors
    # split + lapply/sapply
    

    # split
    # split a vector according to a factor variable
    # always returns a list



Statistics:


    rnorm   # normally-distributed random values 
    dnorm   # probability density at a point
    pnorm   # cumulative distribution 
    qnorm   # quantile ("inverse" of pnorm, pnorm^1)

    # Distributions: Normal, Binomial, Exponential, Gamma, Poisson, ...

    rpois   
    runif
    rbinom
    rexp


    set.seed(1)
    sample(1:10, 4)  // sample 4 random values from the vector 1:10


    q <- quantile(df$var1, probs=c(0.5, 0.75, 0.9))
    q <- quantile(df$var1, probs=seq(0,1,by=0.2))

    # use cut with quantile for categorizing continuous variables
    c <- cut(m$Rank, breaks=q, labels=1:5)
    cut2



Machine Learning:


    Hierarchical Clustering / K-means clustering
    
    dist(data.frame)
    - gives pair-wise distances between all points
    
    hclust(dist)
    - produces dendrogram
    
    heatmap(matrix)
    
    k-means clustering
    - initial guess where the cluster "centroids" are
    - cluster points based on centroids
    - recalculate centroids
    - repeat
    
    df <- data.frame(x, y)
    kmeansObj <- kmeans(df, centers=3)
    par(mar=repo(0.2,4))
    plot(x,y, col=kmeansObj$cluster, pch=19, cex=2)
    points(kmeansObj$centers, col=1:3, pch=3, cex=3, lwd=3)



Profiling and Debug:

    system.time()   # times an operation
    Rprof           # keeps track of call stack and how much time is spent in each function
        $by.total
        $by.self
    summaryRprof
    debug()



## <a name="readingdata"></a>Reading data


    read.table 
    write.table

    read.csv

    readLines           # read/write text files
    writeLines 

    source 
    dump                # textual format

    dget 
    dput                # textual format

    load 
    save

    unserialize 
    serialize

    download.file
    setInternet2(use=TRUE)      // for windows https unsupported url scheme errors


    # fixed-width cols:
    fw <- read.fwf("getdata-wksst8110.for", widths=c(-1,9,-5,4,4,-5,4,4,-5,4,4,-5,4,4), header=F, skip=4)
    str(fw)
  
    file
    gzfile
    bzfile
    url
    con <- file("foo.txt", "r")
    ?connection


Reading from databases:

    library(sqldf)
    sqldf
    
    mysql
    dbconn <- dbConnect(MySQL(), user="", host="", db="")
    result <- dbGetQuery(dbconn, "show databases;")
    result <- dbListTables(dbconn)
    result <- dbListFields(dbconn, "tableName")
    result <- dbReadTable(dbconn, "tableName")
    query <- dbSendQuery(dbconn, "select * from blah") // doesn't return results; use fetch()
    results <- fetch(query, n=10)
    dbClearResult(query)
    dbDisconnect(dbconn)


* hdf5 - hierarchical data format (for large datasets)
    * groups containing datasets + metadata
    * optimize reading data from disc in R

.

    source("http://bioconductor.org/biocLite.R")
    biocLite("rhdf5")
    library(rhdf5)
    h5createFile("example.h5")
    h5createGroup("example.h5", "foo")
    h5createGroup("example.h5", "bar")
    h5ls("example.h5")
    h5Write(A, "example.h5", "foo/A")
    h5Write(c(12,13,14), "example.h5", "foo/A", index=list(1:3,1)) // fill in sub-section 
    A <- h5Read("example.h5", "foo/A")

* other Dbs:
    * RPostgreSQL
    * RODBC
    * RMongo


Reading from the web:

    con <- url("http://...")
    readLines(con)
    close(con)
    htmlTreeParse(con) // for XML
    library(httr)
    con = GET(url)
    data = content(con, as="text")
    parsed = htmlParse(data, asText=T)
    xpathSApply(parsed, "//title")
    con = GET(url, authenticate("user", "pass) )
    goog = handle("http://google.com")
    GET(handle=goog, path="/")
    GET(handle=goog, path="search")


Web APIs:

    myapp = oauth_app("twitter", key="your-consumer-key", secret="your-consumer-secret")
    sig = sign_oauth1.0(myapp, token="your token", token_secret="your token secret")
    h = GET("https://...", sig)
    json1 = content(h)
    json2 = jsonlite::fromJSON(toJSON(json1))


Images:

    jpeg
    readbitmap
    png
    EBImage (bioconductor)


GIS (geographic info systems):

    rdgal
    rgeos
    raster

Music:

    tuneR
    seewave


XML:

    install.packages("XML")
    library(XML)
    doc <- xmlTreeParse("restaurants.xml", useInternal=T)
    root <- xmlRoot(doc)
    xmlName(root)
    names(root)             # names of child elements
    root[[1]]               # first child element
    root[[1]][[1]]          # first grandchild element
    
    xmlValue(root[[1]][[1]][[1]])   # recursively concats all values from all child nodes
    
    xmlValue(root[[1]][[1]])
    [1] "41021206Frankford2NORTHEASTERN"
    xmlSApply(root[[1]][[1]],xmlValue)
    name zipcode neighborhood councildistrict policedistrict location_1 
    "410" "21206" "Frankford" "2" "NORTHEASTERN" "" 

    # retrieve all <zipcode>{value}</zipcode> values for all nodes under root (recursively) 
    xpathSApply(root,"//zipcode",xmlValue) 


Excel xlsx:

First need to setup JAVA_HOME and PATH to avoid these errors:

    Error : .onLoad failed in loadNamespace() for 'rJava', details:
    call: inDL(x, as.logical(local), as.logical(now), ...)
    error: unable to load shared object 'C:/MyPrograms/R/R-3.1.2/library/rJava/libs/x64/rJava.dll':
    LoadLibrary failure: The specified module could not be found. <-- accompanied by pop-up window complaining about not finding jvm.dll
    Error: package or namespace load failed for ‘rJava’
    
    Error : .onLoad failed in loadNamespace() for 'rJava', details:
    call: fun(libname, pkgname)
    error: JAVA_HOME cannot be determined from the Registry
    Error: package or namespace load failed for ‘rJava’
    
    64-bit version
    > version
    _ 
    platform x86_64-w64-mingw32 
    arch x86_64 
    os mingw32 
    system x86_64, mingw32 
    status 
    major 3 
    minor 1.2 
    year 2014 
    month 10 
    day 31 
    svn rev 66913 
    language R 
    version.string R version 3.1.2 (2014-10-31)
    nickname Pumpkin Helmet 
    
    Sys.setenv(JAVA_HOME="C:\\fox\\java\\jre")
    p <- Sys.getenv("PATH")
    p1 <- paste("C:\\fox\\java\\jre\\bin\\j9vm;", p, sep="")
    Sys.setenv(PATH=p1)


    install.packages("xlsx")
    library(xlsx)
    
    data <- read.xlsx("getdata-data-DATA.gov_NGAP.xlsx",
                      sheetIndex=1,
                      rowIndex=c(18:23),
                      colIndex=c(7:15))
    sum(dat$Zip*dat$Ext,na.rm=T) 
    
    data <- read.xlsx("ThreeMysterySecurities.xlsx",
                      sheetIndex=1,
                      rowIndex=4:244,
                      colIndex=c(1:5,7:9))


## <a name="etl"></a>Extract/Transform/Load (ETL) data


[https://github.com/natapone/RepData_PeerAssessment1/blob/master/PA1_template.md](https://github.com/natapone/RepData_PeerAssessment1/blob/master/PA1_template.md)

    subset
    aggregate
    which
    dplyr 

    df <- read.table( ... )      # df is a data frame
    df <- read.csv( ... )


Subset rows in data frame:

    x <- df[!is.na(df$Ozone) & df$Ozone > 31 & df$Temp > 90,]

    df[(df$var1 <= 1 & df$var3 > 11), ]
    df[(df$var1 <= 1 | df$var3 > 11), ]
    
    df[(which(df$var1 > 8),]        # no NAs




Dealing with NA:

    x <- c(1, NA, 2, NA, 4)
    bad <- is.na(x)
    x[!bad]

    complete.cases(x, y)  // elements where both x and y are not NA


Sorting:

    df[order(df$var1),] 
    df[order(df$var1, df$var3),] 



RESHAPING data frames:

    melt(X, id=c("var1", "var2", "var3", measure.vars=c("var4", "var5"))
    dcast(Xmelt, var1 ~ variable)
    dcast(Xmelt, var1 ~ variable, mean)

    # or, split + apply + unlist == sapply
    ddply



dplyr: for transforming data frames   
[https://github.com/tom43214/RepData_PeerAssessment1/blob/master/PA1_template.md](https://github.com/tom43214/RepData_PeerAssessment1/blob/master/PA1_template.md)   


    library(dplyr)

    # %>% pipeline variable for chaining functions
    imputed <- activity %>%
               group_by(interval) %>%
               mutate(steps = ifelse(is.na(steps), mean(steps, na.rm = TRUE), steps))
    
    arrange
    filter
    select
    mutate
    rename
    summarize
    
    select(df, var2:var4)               # select subset of columns
    select(df, -(var2:var4))            # select subset of columns by exclusion
    filter(df, var1 > 30 & var2 > 80)
    arrange(df, var3)
    rename(df, newvar = oldvar, ...)
    mutate(df, newvar = var2 - mean(var2, na.rm=T))
    df2 <- mutate(df, newvar = factor( 1 * (var2 > 80), labels = c("cold", "hot")))
    gdf2 <- group_by(df2, newvar)
    summarize(gdf2, v1 = mean(var1), v2 = max(var2), v3 = median(var3), count = n(), uniques = n_distinct() )
    mutate(df, year= as.POSIXlt(date)$year + 1900)
    
    tbl_df
    View


tidyr: for cleaning/transforming data frames

    library(tidyr)

    # gather: key/ value are new column names in the resulting "tidy" dataset
    # the column names are put as values in the "key" column, and the column values are put in the 'value' column
    gather(df, key, value, columns...)      

    spread()  # complement of gather

    # separate the values in a column (by the sep regex) into new columns col1 and col2
    separate(df, col, into=c("col1", "col2"), sep=)   

    
    # Merge
    merge(x,y,by,by.x,by.y,all)
    merge(X,Y, by.x="y_id", by.y="id", all=T)
    
    intersect(names(X),names(Y))
    join
    join_all    # works on "id" col onl



## <a name="exploringdata"></a>Exploring data

    str
    summary
    nrow
    ncol

    # pair-wise scatterplots between all datasets
    pairs(cbind(gwn,MSFT,SP500))



## <a name="plotbase"></a>Plotting: Base

Plotting systems in R (note: can't be mixed)

* Base - build piece by piece with multiple commands
* Lattice - build all at once
* ggplot2 - combo

Base plotting system:

* boxplot
    * shows median, 25%, 75%
    * varwidth - shows size of sample
* barplot
    * bar chart
* density
    * like histogram, shows % of observations across values
* hist(x)
    * histogram
    * shows distribution
* plot(x, y)
    * scatterplot X vs Y
    * par - for setting parameters
    * mfrow - multiple plots on same device
    * lines
    * points
    * text
    * title
    * mtext
    * axis
    * example(points)
    * abline
    * legend


## <a name="plotlattice"></a>Plotting: Lattice

* Packages: 
    * graphics
    * grDevices
    * lattice
    * grid
* xyplot: scatterplot
* bwplot - box and whiskers plot (boxplot)
* histogram
* stripplot - like boxplot but with actual points
* dotplot - plot dots on "violin strings"
* splom - scatterplot matrix. like pairs in base plotting system
* levelplot - image data
* contourplot - image data
* trellis.par.set
* llines

.

    library(lattice)
    library(datasets)
    airquality <- transform(airquality, Month = factor(Month))

    # xyplot( y ~ x, data = <data frame containing y and x> )
    xyplot( Ozone ~ Wind | Month, data = airquality, layout=c(5,1) )
    
    p <- xyplot( y ~ x, data = dd)
    print(p)

    xyplot( y ~ x | f, function(x, y, ...) {
          panel.xyplot(x, y, ...)
          panel.lmline(x, y, col="red")
    })



## <a name="plotggplot2"></a>Plotting: ggplot2



* ggplot: "grammar of graphics"
* quick ref: [http://sape.inf.usi.ch/quick-reference/ggplot2/geom](http://sape.inf.usi.ch/quick-reference/ggplot2/geom)
* basic components of ggplot:
    * a data frame
    * aesthetic mappings - x,y, color, shape, size
    * geoms - geometrical objects like points, lines, bars, 
    * facets - for conditional plots
    * stats - statistical transforms like binning (bar charts), quantiles, smoothing
    * scales - applied to aesthetics
    * coordinate system - e.g cartesian (x-y) vs polar (pie)
    * factors are important for indicating subsets of data

.

    library(ggplot2)
    g <- ggplot(dd, aes(x, y)) # aesthetics
    summary(g)
    
    p <- g + geom_point() + geom_smooth(method="lm") + facet_grid(. ~ f)
    print(p)
    
    g + geom_point(color="blue")
    g + geom_point(aes(color=f)) # color by data
    g + geom_line() + ylim(-3,3) # will exclude outliers
    g + geom_line() + coord_cartesion(ylim = c(-3,3)) # will NOT exclude outliers
    
    facet_wrap()
    xlab()
    ylab()
    labs()
    ggtitle()
    theme(legend.position="none")
    theme_gray()
    theme_bw()


    ggplot(data=merge.sum, aes(x=Year, y=Emissions/1000)) + 
        geom_line(aes(group=1, col=Emissions)) + 
        geom_point(aes(size=2, col=Emissions)) + 
        ggtitle(expression('Total Emissions of PM'[2.5])) + 
        ylab(expression(paste('PM', ''[2.5], ' in kilotons'))) + 
        geom_text(aes(label=round(Emissions/1000,digits=2), size=2, hjust=1.5, vjust=1.5)) + 
        theme(legend.position='none') + scale_colour_gradient(low='black', high='red')


    qplot(x, y, data=dd)
    qplot(displ, hwy, data=mpg)
    qplot(displ, hwy, data=mpg, color=drv)
    qplot(displ, hwy, data=mpg, geom=c("point", "smooth"))
    qplot(displ, hwy, data=mpg, geom=c("point", "smooth"), method="lm")      # linear regression
    qplot(hwy, data=mpg, fill=drv)      # histogram
    qplot(displ, hwy, data=mpg, facets=. ~ drv)          # facets = row x col.  use "." for 1 row/col
    qplot(displ, hwy, data=mpg, facets= drv ~ ., binwidth=2)   # 1 col (".").  rows determined by drv variable.

    qplot(Date.date, Adj.Close, data=df, geom=c("point", "line"))



    smoothScatter

    qqplot      # plot quantiles

    matplot (spaghetti)

    library(maps)
    image       # heatmap (a 2D histogram)



Colors:

    rgb( r,g,b, alpha=)
    library(colorspace)
    
    grDevices
    
    colorRamp
    colorRampPalette
    colors()

* RColorBrewer
    * color palettes for plotting specific kinds of data patterns
    * 3 kinds of palettes:
    * sequential
    * for data from low to high
    * divergent
    * for divergence from some center point (e.g deviation from mean)
    * qualitative
    * categorical data (no ordering)

.

    library(RColorBrewer)
    cols <- brewer.pal(3, "BuGn")
    pal <- colorRampPalette(cols)
    image(volcano, col=pal(20))



## <a name="plotdevices"></a>Plotting: devices

    pdf(filename)
    dev.off()  <-- remember to use this when plotting to files!
    window()
    dev.cur()
    dev.set
    dev.copy
    dev.copy2pdf


    Further graphing resources: 
    R Graph Gallery
    ggplot2,ggplot2 basic introduction
    lattice package,lattice introduction
    R bloggers
    http://www.biostat.wisc.edu/~kbroman/topten_worstgraphs/
    https://www.facebook.com/notes/facebook-engineering/visualizing-friendships/469716398919
    How to display data badly
    The visual display of quantitative information
    Creating more effective graphs
    R Graphics Cookbook
    ggplot2: Elegant Graphics for Data Analysis
    Flowing Data




## <a name="knitr"></a>Knitr

    library(knitr)
    knit
    knit2html
    browseURL
    
    ```{r codeChunkName,echo=FALSE}
    ```{r codeChunkName,echo=FALSE,results="hide"}
    ```{r codeChunkName,fig.height=4}
    ```{r codeChunkName,results="asis"}
    ```{r codeChunkName,cache=TRUE}
    
    ```{r setopts,echo=FALSE}
    opts_chunk$set(echo=FALSE, results="hide")
    ```
    
    inline:
    ```{r,echo=FALSE}
    time <- sys.time()
    ```
    The current time is `r time`










