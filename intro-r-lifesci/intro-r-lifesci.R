
## ----, echo=FALSE, message=FALSE, eval=TRUE------------------------------
# Set eval=TRUE to hide all results and figures.
# This sets defaults. Can change this manually in individual chunks.
# Must load knitr so opts_chunk is in search path.
library(knitr)
opts_chunk$set(results="hide", message=FALSE, fig.show="hide", fig.keep="none")


## ------------------------------------------------------------------------
2+2
5*4
2^3


## ------------------------------------------------------------------------
2+3*4/(5+3)*15/2^2+3*4^2


## ------------------------------------------------------------------------
# Notice that this is a comment. 
# Anything behind a # is "commented out" and is not run.
sqrt(144)
log(1000)


## ------------------------------------------------------------------------
log(1000)
log(1000, base=10)
sqrt(log(1000, base=10))


## ----, eval=FALSE--------------------------------------------------------
## # Some simple numeric vectors:
## 1:5
## 6:10
## 1:5 + 6:10
## 1:100
## 
## # Get some help with the seq() function, then create a vector from 2 to 200 by 2s.
## # Notice how the seq() function works -- the `to` argument will never be exceeded.
## help(seq)
## seq(from=2, to=200, by=4)


## ------------------------------------------------------------------------
x <- 5
x

y <- 42
y

y-x
z <- y-x
z

x <- 1:5
y <- 6:10
x
y
x+y
x*y
x^y


## ------------------------------------------------------------------------
ls()
rm(x)
ls()
x # oops! you should get an error because x no longer exists!
rm(y,z)


## ------------------------------------------------------------------------
class(42)
class(log)
name <- "Stephen"
name
class(name)


## ------------------------------------------------------------------------
toupper(name) # name is an object of class character. methods or functions are associated with certain classes.
toupper(log) # can't run a function that expects character on an object of class function


## ------------------------------------------------------------------------
# Get some help with ?c
x <- c(1,3,5)
x
class(x)
length(x)


## ------------------------------------------------------------------------
y <- c("My", "name", "is", "Stephen")
y
class(y)
length(y)
y <- c(y, "Turner")
y
length(y)


## ------------------------------------------------------------------------
sum(x)


## ------------------------------------------------------------------------
y
z <- c(x,y)
z 
class(z)


## ------------------------------------------------------------------------
z
sum(z)


## ------------------------------------------------------------------------
# Create the vector.
x <- 101:150

# Get the first element.
x[1]

# Get the 42nd element.
x[42]

# Get the 20th through the 25th elements. 
x[20:25]

# If you try to access elements that don't exist, you'll return missing values.
# Missing values are represented as NA
x[45:55] #NA is missing value!


## ------------------------------------------------------------------------
data(mtcars)
class(mtcars)
mtcars


## ------------------------------------------------------------------------
head(mtcars)
length(mtcars)
dim(mtcars)
dim(mtcars)[1] # number of rows (individual cars in the survey)
dim(mtcars)[2] # number of columns (number of variables measured)
str(mtcars)


## ------------------------------------------------------------------------
# display the number of cylinders for each car. 
mtcars$cyl
# first display MPG for all vehicles, then calculate the average.
mtcars$mpg
mean(mtcars$mpg)


## ------------------------------------------------------------------------
head(mtcars)
mtcars[1:4, 1:2]


## ------------------------------------------------------------------------
subset(mtcars, cyl==6)
subset(mtcars, cyl>6)
subset(mtcars, mpg>=20 | disp<100)
subset(mtcars, cyl==6, select=c(mpg, disp))
subset(mtcars, cyl>=6 & mpg>=15, select=c(mpg, cyl, qsec))


## ------------------------------------------------------------------------
# Display the number of cylinders.
mtcars$cyl
with(mtcars, cyl)

# Compute the senseless value described above. Both return the same results.
mtcars$mpg * mtcars$cyl / mtcars$disp
with(mtcars, mpg*cyl/disp)


## ------------------------------------------------------------------------
plot(mtcars$mpg)


## ------------------------------------------------------------------------
hist(mtcars$mpg)
hist(mtcars$mpg, breaks=10)
hist(mtcars$mpg, breaks=10, col="black")


## ------------------------------------------------------------------------
# This would also work, but let's use with().  
# plot(mtcars$disp, mtcars$mpg)
with(mtcars, plot(disp, mpg))


## ------------------------------------------------------------------------
with(mtcars, plot(disp, mpg, pch=16))
with(mtcars, plot(disp, mpg, pch=16,  col="red"))
with(mtcars, plot(disp, mpg, pch=16,  col="red", main="MPG vs Displacement"))
with(mtcars, plot(disp, mpg, pch=16,  col="red", main="MPG vs Displacement", 
                  ylab="Fuel Economy (MPG)", xlab="Displacement (cu. in.)"))


## ----, echo=FALSE--------------------------------------------------------
with(subset(mtcars, cyl>4), plot(disp, hp, pch=16, col="blue",
                                 xlab="Displacement (cu. in.)", ylab="Gross Horsepower", 
                                 main="Horsepower vs Displacement for 6 and 8-cylinder vehicles"))


## ------------------------------------------------------------------------
mtcars_8cyl <- subset(mtcars, cyl==8)
mtcars_8cyl


## ----, eval=FALSE--------------------------------------------------------
## getwd()
## help(write.table)
## help(write.csv)


## ----, eval=FALSE--------------------------------------------------------
## getwd()
## setwd("~/Desktop/R")


## ----, eval=FALSE--------------------------------------------------------
## write.csv(mtcars_8cyl, file="cars8.csv")


## ----, eval=FALSE--------------------------------------------------------
## help(read.table)
## help(read.csv)


## ----, eval=FALSE--------------------------------------------------------
## rm(mtcars_8cyl)
## mtcars_8cyl
## cars8 <- read.table(file="cars8.csv", header=TRUE, sep=",", row.names=1)
## cars8
## rm(cars8)
## cars8 <- read.csv(file="cars8.csv", header=TRUE, row.names=1)
## cars8


## ----, eval=FALSE--------------------------------------------------------
## # Install only once.
## install.packages("ggplot2")
## 
## # Load the package every time you want to use it.
## library(ggplot2)


## ----, eval=FALSE--------------------------------------------------------
## # Download the installer script
## source("http://bioconductor.org/biocLite.R")
## 
## # biocLite() is the bioconductor installer function.
## # Run it without any arguments to install the core packages or update any installed packages.
## # This requires internet connectivity and will take some time!
## biocLite()


## ----, eval=FALSE--------------------------------------------------------
## # Do only once
## source("http://bioconductor.org/biocLite.R")
## biocLite("DESeq2")
## 
## # Every time you need to use the DESeq2 package
## library(DESeq2)


## ----savedata, eval=FALSE, echo=FALSE------------------------------------
## # This chunk only run for creating the data for class.
## library(pasilla)
## data(pasillaGenes)
## write.csv(counts(pasillaGenes), file="pasilla_counts.csv")
## write.csv(pData(pasillaGenes)[, c("condition", "type")], file="pasilla_metadata.csv")


## ------------------------------------------------------------------------
# Bioconductor packages.
# Use the installation instructions at http://www.bioconductor.org/install/
# if you haven't already installed these (install once, load every time)
library(Biobase)
library(DESeq2)


## ------------------------------------------------------------------------
# Load the count data
pasillacounts <- read.csv("data/pasilla_counts.csv", header=TRUE, row.names=1)
head(pasillacounts)

# Load the sample metadata
pasillameta <- read.csv("data/pasilla_metadata.csv", header=TRUE, row.names=1)
pasillameta


## ------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData=pasillacounts, colData=pasillameta, design = ~condition)
dds


## ----rundeseq------------------------------------------------------------
# Normalize, estimate dispersion, fit model
dds <- DESeq(dds)


## ----results_simpledesign, fig.keep='all', fig.show='asis'---------------
# Extract results
res <- results(dds)

# View, reorder, and write out results
head(res)
mcols(res)
res <- res[order(res$padj), ]
head(res)
write.csv(res, file="results/pasilla_results.csv")

# Create MA Plot
plotMA(dds, ylim=c(-2,2))


## ----viz, fig.keep='all', fig.show='asis'--------------------------------
# Transform
rld <- rlogTransformation(dds, blind=TRUE)

# Principal components analysis
plotPCA(rld, intgroup=c("condition", "type"))

# Hierarchical clustering analysis
distrl <- dist(t(assay(rld)))
plot(hclust(distrl))


## ----sessioninfo, results='markup'---------------------------------------
sessionInfo()


