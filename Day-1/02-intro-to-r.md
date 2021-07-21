### Beyond the basic object classes

As we discussed, one of the major advantages of R is its enormous user base that are continuously developing and releasing packages. Implementing the additional functionality of these packages often requires **more complex data object classes to be created**, which are generally related in some way to one or more of the basic data structures in R that we have discussed.

The general approach used to create these novel classes is referred to as *object-orientated programming (OOP)*. Although we will not go into any detail on OOP in this workshop, it is useful to know that several OOP methods are used to create classes in R.

The [S4 OOP approach](http://adv-r.had.co.nz/S4.html) is perhaps the most relevant to this workshop, as it is heavily used by packages from the [Bioconductor project](http://bioconductor.org/), which we will be using on Day 3. S4 OOP provides a flexible approach to creating objects with multiple *slots*, each with defined names and object classes.

An example of a package that has used an S4 OOP approach to create objects of novel classes is the Bioconductor package [*GenomicRanges*](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html), which provides representation and query of genomic region data in R through object classes such as `GRanges`.

To efficiently represent genomic regions data, `GRanges` class objects, at least 3 key slots (each with their own associated class for that vector) will be needed:
- a chromosome slot: with class `character` (e.g. chr1)
- a start coordinate slot: with class `integer` (e.g. 338833)
- an end coordinate slot: with class `integer` (e.g. 338987)

Creating object classes in this way is desirable as it allows *accessor functions* to be created, which allow very simple interaction with the objects of this class. For example, simply using the `chr()` accessor function to easily extract all chromosome identities of the genomic regions in my GRanges object.

An excellent tutorial describing S4 classes and their use in the [Bioconductor project](http://bioconductor.org/) is available [here](https://bioconductor.org/help/course-materials/2017/Zurich/S4-classes-and-methods.html). While this is not a topic you need understand in great detail, it is worth understanding the basic principles.


## Functions

Beyond the functions implemented in base R and packages that you install, R allows you to create user defined functions, which can perform any functionality that you can define.

Defining your own functions can be useful when you want to perform a specific set of tasks repeatedly on some input(s) and return a defined output. Furthermore, once defined functions become part of your global environment and are therefore preserved for future use they minimize the need for repetitive code.

Functions are created using using `function()` with the assignment operator `<-`. The arguments you use in the `function()` command define the variables that those arguments will be assigned to when you call the function. The last line of the function defines what output is returned.

Let's define a basic function as an example.
```r
# define the function that takes a single argument and returns a single argument
myfun <- function(x){
  y <- x + 1
  return(y)
}

# call the function
myfun(x = 10)

# assign the output to a new variable
var1 <- myfun(x = 10)
var1
```

Functions can have as many arguments as you specify. The names of these arguments are only assigned as variables within the function, however it is good practice to avoid using arguments with the same name as variables already existing in your global environment.

For example, if I already have a variable named `x` in my environment, I should avoid using x to define the name of the arguments to my function.
```r
myfun2 <- function(num1, num2){
  num3 <- num1 + num2
  return(num3)
}

# call the function
myfun2(num1 = 10, num2 = 11)
```

## Loops

Loops are used to iterate over a piece of code multiple times, and can therefore be used to achieve specific tasks. The most often used type of loop in R is the `for()` loop, which will evaluate the contents of the loop for all of the values provided to the `for()` function.

For example:
```r
x <- c(1,2,3,4,5,6,7,8,9)

# define the loop, using [i] to define the elements of x used in each iteration
for(i in 1:length(x)){
  print(x[i] * 10)
}
```

We may wish to save the output of each iteration of the loop to a new variable, which can be achieved using the following approach:
```r
x <- c(1,2,3,4,5,6,7,8,9)

# create variable to save results to
y <- c()

# define and run the loop
for(i in 1:length(x)){
  y[i] <- x[i] * 10
}

# print the new vector to console
y
```

## Basic data visualization

R is a very powerful tool for visualization and provides a large amount of user control over the plotting attributes. Basic visualization in R is achieved using the `plot()` function.

```r
# generate a set of random numbers to plot
x <- rnorm(1000, mean = 10, sd = 2)
y <- rnorm(1000, mean = 20, sd = 1)

# plot x against y to produce a scatter plot
plot(x, y)

# add labels, title, and color
plot(x, y,
     main = "X vs. Y",
     xlab = "X values",
     ylab = "Y values",
     col = "red")
```

R can also be easily used to generate histograms:
```r
# generate a simple histogram for x
hist(x, col = "indianred")

# the breaks argument can be used to change how the intervals are broken up
hist(x, col = "indianred", breaks=10)
hist(x, col = "indianred", breaks=50)
hist(x, col = "indianred", breaks=100)
```

> This is an *extremely* brief introduction to data visualization in R, however there are many excellent online resources for learning more about how to perform effective data visualization in R.

### Visualization specific packages

There are a large number of packages designed for specifically for visualization and are very useful in bioinformatic analyses. We won't cover these here since they are covered extensively elsewhere, but some useful visualization packages to be aware of include:
- [ggplot2](https://ggplot2.tidyverse.org/reference/ggplot.html)
- [ggpubr](http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/)
- [ploty](https://plotly.com/python/)

Importantly, visualization implemented in these packages form the basis for some bioinformatics specific data visualization packages that we will explore later in the workshop.

## Import and export tabular data

Tabular data are often stored a text files where the individual fields containing data points are separated by punctuation points. Three functions exist in base R to facilitate reading in tabular data stored as text files.

- `read.table()` - general function for reading in tabular data with various delimiters
- `read.csv()` - used to read in comma separated values files, where commas are the delimiter
- `read.delim()` - used to read in files in which the delimiters are tabs

Use `read.table()` to read in the *all_counts.txt*  that we used on day 1. Since `read.table()` is a general function for loading in tabular data, we need to specify the correct separator/delimiter value using the `sep` argument. Tab delimited data is specified using `\t` in the sep argument.
```r
# using read.table
counts <- read.table(file = "all_counts.txt",
                     sep = "\t", header = TRUE, stringsAsFactors = FALSE)
### Note 1: header accepts logical value indicating if the first row are column names (default FALSE)
### Note 2: we use stringsAsFactors

# check class, dimensions and structure
class(counts); dim(counts); str(counts)
```

Now use `read.delim()`. An important difference between `read.delim()` and `read.table()` are the default setting for the *sep* and *header* arguments. By deault in `read.delim()`, *sep* is set to `\t` and the header argument is set to `TRUE`, so we do not need to explicity call those arguments.
```r
# using read.delim
counts <- read.table(file = "all_counts.txt", stringsAsFactors=FALSE)

# check class, dimensions and structure
class(counts); dim(counts); str(counts)
```

`read.csv()` is used in exactly the same way `read.delim()`, however the `file` specified must contain data separated by commas and have the extension `.csv`.

When datasets get very large, these base R functions can be quite slow. Although we do not have time to cover them here, the [*readr*](https://readr.tidyverse.org/) and [*data.table*](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html) packages contain functions that introduce extremely fast ways of reading data into an R environment.

### Exporting tabular data

The major functions in base R that exist for writing tabular data to file are `write.table()` and `write.csv()`. Similarly to the read functions, `write.table()` provides a more generalized solution to writing data that requires you to specify the separator value.

In both functions, the first argyument specifies the object in your global environment that you wish to write to file. The second argument defines the absolute or relative path to the location you wish to save this file.
```r
# subset counts for first 5 columns and first 2000 genes
counts_sub <- counts[1:2000, 1:5]

# write to tab delimited file using write.table
write.table(counts_sub, file = "all_counts_sub.txt", sep = "\t")
```

In contrast, `read.csv()` does not require you to set the delimitor value, and by default writes data to comma separated value files (.csv).
```r
write.csv(counts_sub, file = "all_counts_sub.csv")
```

## Save R objects to file

It is possible to save R objects to a file that maintains all of their attributes so that they can be easily loaded back into a new R session. This functionality is achieved using the `save()` and `load()` functions.

`save()` accepts the names of the R objects you wish to write to a file (which will have extension `.Rdata`) as the first arguments, and the file path where you wish to write this file under the *file* argument. For example:
```r
# create some R objects
x <- c(1.63, 2.25, 3.83, 4.99)
y <- c(TRUE, FALSE, TRUE, TRUE)
z <- c("a", "b", "c", "d")

# save all 3 objects to one file
save(x, y, z, file = "my_r_objects.rdata")
```

These objects can then be loaded back into your global environment using the `load()` function with the absolute or relative file path. The objects will appear in your new environment with exactly the same names as in your previous environment.
```r
load(file = "my_r_objects.rdata")
```

Single R objects can be saved and restored in a similar way using the `saveRDS()` and `readRDS()` functions. Files saved using RDS must take on the `.RDS` extension.
```r
# save a single object to a specific file path
saveRDS(x, file = "my_r_object.rds")

# use the assignment operator to assign object to any variable you want
x <- readRDS(file = "my_r_object.rds")

# I changed my mind, I want the object to be assigned to variable `a` in my new env
a <- readRDS(file = "my_r_object.rds")
```
