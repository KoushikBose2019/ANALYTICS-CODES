if(!require(installr)) {
  install.packages("installr"); require(installr)} #load / install+load installr

# using the package:
updateR() # this will start the updating process of your R installation. 

# working directory
getwd()
# To set working directory 
setwd("D:/R")
getwd()

# Simple math
5+6-3
12+5/6-3/4+2.
(1+3)*9

# mathematical operations
5+5
5-5
5*5
5/5

# value of pi
pi

# other mathematical functions
25^2
sqrt(25)
abs(-5.6) # to get absolue value
factorial(3) #3!
exp(9) # exponent of 9 
# log functions
log(5, base = 3) # user defined base
log10(5) # base 10
log(5, base=10)
log2(5) # base 2

# assignment operator
a <- 1
b <- a
a = 2
a<-45


l <- c(3, 9, 10)
v <- c("A", "B", "C")
v

# available files in the working directory
dir()
# To see workspace
ls()

# remove some object
rm(l)

# help
?rm

rm(list=ls())

# create a sequence
a <- 2:10

# create a sequence of repeated number
r <- rep(1, 10)

# R Objects

# vector - most basic object in R. It contains elements of same class
# vector - character, numeric, integer, complex, logical 
# use of vectorization
a <- c(2,3,4)
b<-c(11, 12, 13)
c <- a*b
c

# creating a logical vector
k= 5
d <- (k< -5)
d
class(d)

# Know the class of an object
class(b)

# example of list - a vector with different classes of objects
li <- list(c(5, "red",3,"blue"))
li
l <- list(c("a", "b", "c"), 1:10, TRUE)
l
lst <- list(c(1,3,5,7), c("A","B","C"))
lst

# Matrix in R
cells <- c(1,26)
rnames <- c("R1", "R2")
cnames <- c("C1", "C2")
mymatrix <- matrix(cells, nrow=2, ncol=2, byrow=FALSE,dimnames=list(rnames, cnames))
mymatrix

# matrix creation is columnwise
m=matrix(1:10,3,4)
m
m[2,1]
#check attributes of a matrix
attributes(m)

# create a matrix from a vector
m3=matrix(1:6)
m3
dim(m3)=c(2,3)
m3

# create a matrix by binding columns or rows
x=1:6
y=11:16
x
y
T=cbind(x,y) # create by column
T
T[2,1]

L=rbind(x,y)
L
# functions

# numeric functions

# absolute value
abs(x)
# square root
s<-sqrt(abs(x))
# round
round(s, digits=1)
# signif
signif(s, digits=6)
# natural logarithm
log(abs(x))
# common logarithm
log10(abs(x))

age <- c(1,3,5,2,11,9,3,9,12,3)
weight <- c(4.4,5.3,7.2,5.2,8.5,7.3,6.0,10.4,10.2,6.1)
mean(weight)

median(age)
min(weight)
max(age)
range(weight)
sd(weight)
var(weight)

plot(age,weight)
cor(age,weight)

#Data Frames
patientID <- c(1, 2, 3, 4)
age <- c(25, 34, 28, 52)
diabetes <- c("Type1", "Type2", "Type1", "Type1")
status <- c("Poor", "Improved", "Excellent", "Poor")
patientdata <- data.frame(patientID, age, diabetes, status)
patientdata

patientdata[1:2]

patientdata[c("diabetes", "status")]

patientdata$age

table(patientdata$diabetes, patientdata$status)