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


# Dataframes are matrices with different classes of objects
# Importing Data

data()
install.packages("Rcpp")
install.packages("readxl")
install.packages("readr")
#Read data from file
car<- read.csv("D:/Data/Car.csv")

# knowing structure of a data
str(car)

# column heading of a data frame
names(car)
colnames(car)

# dimension of a matrix or a data frame
dim(car)
nrow(car)
ncol(car)

# Display few observations of an object
head(car,3)
tail(car)
car[10:15,] # for viewing specific rows from the middle
car[10:15,2:4]
car[,2:4]

# drop some variable from data frame
car$Sound <- NULL

# getting data out of R work space
getwd()
setwd("D:/R")
write.table(names(car),"names_car.txt")
write.csv(tail(car),"tail_car_1.csv")


# Character Functions
a <- "Kolkata City"
# Extract or replace substrings in a character vector
substr(a, 4, 9)
nchar(a)

# Search for pattern
b <- c("DAM", "Deer", "Cat", "daM")
b
grep("daM", b, fixed=T) # looking for exact match

grep("daM", b, ignore.case=T, fixed=F) # ignoring the case

# Split the elements of character vector
KC <- strsplit(a, " ")
frst <- strsplit(a," ")[[1]][1]
scnd <- strsplit(a," ")[[1]][2]
thrd <- "Amazing"
# paste function
paste("x",1:3,sep="_")
paste("Kolkata is an Amazing ",scnd, sep="")
paste(frst,scnd,thrd,sep=" ")

# Uppercase
toupper(a)
# Lowercase
tolower(a)

# date function 
d <- as.Date("2017-12-23") #as: Force an Object to Belong to a Class
class(d)
# day of week
weekdays(d)
# month of year
months(d)
# quarter of year
quarters(d)
# system time
Sys.time()

View(state.x77)

#apply function:
apply (state.x77, 2, median) # 2 means go by column. 1 would have meant "go by row"

# creating a new data frame
Region <- as.data.frame(state.region)
all.states <- as.data.frame(state.x77)

streg <- cbind(all.states,Region)

# knowing structure of a data
str(all.states)

all.states$Name <- row.names(all.states)
rownames(all.states) <- NULL

some.states <- cbind(all.states, Region)
some.states <- some.states[1:10, c("state.region", "Population", "Income")]

# sorting vectors
(sort(some.states$Population)
# or
sort(some.states$Population, decreasing = TRUE)

# sorting data frames
orderIndex.pop <- order(some.states$Population, decreasing = FALSE)
# for vector ordering
some.states$Population[orderIndex.pop]
# dataframe ordering
some.states[orderIndex.pop,]
some.states[order(some.states$Income),]  #display ALL columns
# for decreasing ordering
some.states[order(some.states$Income, decreasing = T),]


# creating subsets
cold.states <- all.states[all.states$Frost > 150, c("Name", "Frost")]
(large.states <- all.states[all.states$Area >= 100000, c("Name", "Area")

# use of merge function
# intersection set
merge(cold.states, large.states)
# full outer join (union)
merge(cold.states, large.states, all = T)
# left outer join
merge(cold.states, large.states, all.x = T)
# right outer join
merge(cold.states, large.states, all.y= T)

# Other example of subset
data1<-read.csv("D:/R/Class 2/data1.csv")
data2<-read.csv("D:/R/Class 2/data2.csv")
#Create subset where Col1 starts with AT
dataAT <- subset(data1, grepl("^AT", Col1),Col1:Col5)
data2AT <- subset(data2, grepl("^AT", Col1),Col1:Col5)