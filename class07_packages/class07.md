Class 07 Functions and Packages
================
Livia Songster
10/23/2019

## Revisit our functions from last class

Ctrl+Alt+i is the shortcut Here we can “source” an R script from a
previous class (this one is posted online).

``` r
source("http://tinyurl.com/rescale-R")
```

## Write a function both\_na()

We want to write a function, called both\_na(), that counts how many
positions in the two input vectors, x and y, have na values

``` r
x <- c(1,2,NA,3,NA)
y <- c(NA,3,NA,3,4)
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
# will output the positions of the na values
which(is.na(x))
```

    ## [1] 3 5

``` r
which(is.na(y))
```

    ## [1] 1 3

``` r
both_na<-function(x,y) {
  if(length(x) != length(y)) {
    stop("x and y are not the same length")
  }
  sum(is.na(x) & is.na(y))
}

both_na(x,y)
```

    ## [1] 1

Write a grade function

``` r
# student 1
st1<- c(100, 100, 100, 100, 100, 100, 100, 90)
# student 2
st2 <- c(100, NA, 90, 90, 90, 90, 97, 80)

grade <- function(data) {
  data[is.na(data)] <- 0 # replace NAs with zero???
  drop <- which.min(data)
  gradethis <- data[-drop] # exclude the column specified by drop
  print(mean(gradethis,na.rm=TRUE))
}

grade(st1)
```

    ## [1] 100

``` r
grade(st2)
```

    ## [1] 91

``` r
data <- read.csv("student_homework.csv")
#transpose the data because my function will parse by column not row
data2 <- data.frame(t(data[-1]))
colnames(data2) <- data[, 1]

lapply(data2,grade)
```

    ## [1] 91.75
    ## [1] 82.5
    ## [1] 84.25
    ## [1] 84.25
    ## [1] 88.25
    ## [1] 89
    ## [1] 94
    ## [1] 93.75
    ## [1] 87.75
    ## [1] 79
    ## [1] 86
    ## [1] 91.75
    ## [1] 92.25
    ## [1] 87.75
    ## [1] 78.75
    ## [1] 89.5
    ## [1] 88
    ## [1] 94.5
    ## [1] 82.75
    ## [1] 82.75

    ## $`student-1`
    ## [1] 91.75
    ## 
    ## $`student-2`
    ## [1] 82.5
    ## 
    ## $`student-3`
    ## [1] 84.25
    ## 
    ## $`student-4`
    ## [1] 84.25
    ## 
    ## $`student-5`
    ## [1] 88.25
    ## 
    ## $`student-6`
    ## [1] 89
    ## 
    ## $`student-7`
    ## [1] 94
    ## 
    ## $`student-8`
    ## [1] 93.75
    ## 
    ## $`student-9`
    ## [1] 87.75
    ## 
    ## $`student-10`
    ## [1] 79
    ## 
    ## $`student-11`
    ## [1] 86
    ## 
    ## $`student-12`
    ## [1] 91.75
    ## 
    ## $`student-13`
    ## [1] 92.25
    ## 
    ## $`student-14`
    ## [1] 87.75
    ## 
    ## $`student-15`
    ## [1] 78.75
    ## 
    ## $`student-16`
    ## [1] 89.5
    ## 
    ## $`student-17`
    ## [1] 88
    ## 
    ## $`student-18`
    ## [1] 94.5
    ## 
    ## $`student-19`
    ## [1] 82.75
    ## 
    ## $`student-20`
    ## [1] 82.75
