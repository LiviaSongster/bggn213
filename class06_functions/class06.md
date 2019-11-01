class06
================
Livia
10/18/2019

# Class 6: R Functions

This is my class 6 work. This will be
**bold**.

### This is R code, exploring read.table()

``` r
# practice opening test text files to explore different types of read.table
test1 <- read.table("test1.txt", sep=",",header=TRUE) 
test1
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
test2 <- read.table("test2.txt", sep="$",header=TRUE) 
test2
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
test3 <- read.table("test3.txt", sep="",header=FALSE) 
test3
```

    ##   V1 V2 V3
    ## 1  1  6  a
    ## 2  2  7  b
    ## 3  3  8  c
    ## 4  4  9  d
    ## 5  5 10  e

## Section 1: Improving analysis code by writing functions

### Part A. Can you improve this analysis code?

``` r
df <- data.frame(a=1:10, b=seq(200,400,length=10),c=11:20,d=NA)
df
```

    ##     a        b  c  d
    ## 1   1 200.0000 11 NA
    ## 2   2 222.2222 12 NA
    ## 3   3 244.4444 13 NA
    ## 4   4 266.6667 14 NA
    ## 5   5 288.8889 15 NA
    ## 6   6 311.1111 16 NA
    ## 7   7 333.3333 17 NA
    ## 8   8 355.5556 18 NA
    ## 9   9 377.7778 19 NA
    ## 10 10 400.0000 20 NA

``` r
df$a <- (df$a - min(df$a)) / (max(df$a) - min(df$a))
df$b <- (df$b - min(df$a)) / (max(df$b) - min(df$b))
df$c <- (df$c - min(df$c)) / (max(df$c) - min(df$c))
df$d <- (df$d - min(df$d)) / (max(df$d) - min(df$d)) 
df
```

    ##            a        b         c  d
    ## 1  0.0000000 1.000000 0.0000000 NA
    ## 2  0.1111111 1.111111 0.1111111 NA
    ## 3  0.2222222 1.222222 0.2222222 NA
    ## 4  0.3333333 1.333333 0.3333333 NA
    ## 5  0.4444444 1.444444 0.4444444 NA
    ## 6  0.5555556 1.555556 0.5555556 NA
    ## 7  0.6666667 1.666667 0.6666667 NA
    ## 8  0.7777778 1.777778 0.7777778 NA
    ## 9  0.8888889 1.888889 0.8888889 NA
    ## 10 1.0000000 2.000000 1.0000000 NA

Here is my function that achieves the same output:

``` r
# this will rescale one column that you input
rescale <- function(x) {
  (x - range(x)[1]) / (range(x)[2] - range(x)[1])
}
df <- data.frame(a=1:10, b=seq(200,400,length=10),c=11:20,d=NA)
df$a <- rescale(df$a)
df
```

    ##            a        b  c  d
    ## 1  0.0000000 200.0000 11 NA
    ## 2  0.1111111 222.2222 12 NA
    ## 3  0.2222222 244.4444 13 NA
    ## 4  0.3333333 266.6667 14 NA
    ## 5  0.4444444 288.8889 15 NA
    ## 6  0.5555556 311.1111 16 NA
    ## 7  0.6666667 333.3333 17 NA
    ## 8  0.7777778 355.5556 18 NA
    ## 9  0.8888889 377.7778 19 NA
    ## 10 1.0000000 400.0000 20 NA
