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
a = 1
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

#Read data from file
car<- read.csv("D:/Data/Car.csv")

# knowing structure of a data
str(car)

# column heading of a data frame
names(car)

# dimension of a matrix or a data frame
dim(car)
nrow(car)
ncol(car)

# Display few observations of an object
head(car,3)
tail(car)
car[10:15,] # for viewing specific rows from the middle

# drop some variable from data frame
car$Sound <- NULL

# getting data out of R work space
getwd()
setwd("D:/R")
write.table(names(car),"names_car.txt")
write.csv(names(car),"names_car_1.csv")


# Character Functions
a <- "Kolkata City"
# Extract or replace substrings in a character vector
substr(a, 4, 9)
nchar(a)

# Search for pattern
b <- c("a", "D", "c", "d")
b
grep("d", b, ignore.case=F, fixed=T) # looking for exact match

grep("d", b, ignore.case=T, fixed=F) # ignoring the case

# Split the elements of character vector
strsplit(a, " ")
strsplit(a," ")[[1]][2]

# paste function
paste("x",1:3,sep="_")
paste("Kol",strsplit(a," ")[[1]][2], sep=" ")

# Uppercase
toupper(a)
# Lowercase
tolower(a)

# date function 
d <- as.Date("1970-01-01") #as: Force an Object to Belong to a Class
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

# knowing structure of a data
str(all.states)

all.states$Name <- row.names(all.states)
rownames(all.states) <- NULL

some.states <- cbind(all.states, Region)
some.states <- some.states[1:10], c("state.region", "Population", "Income")]
some.state <- some.states[1:10, c("income")

# sorting vectors
sort(some.states$Population)
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
large.states <- all.states[all.states$Area >= 100000, c("Name", "Area")]

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
#Create subset where Col1 starts with AT
dataAT <- subset(data1, grepl("AT", Col1),Col1:Col5)

# random sampling
set.seed(10)
sample(1:10,5)
set.seed(8)
sample(1:10,5)
set.seed(10)
sample(1:10,5)

data(iris)
names(iris)
dim(iris)
# creation of training (70%) and test (30%) data
TrainIndex <- sample(1:dim(iris)[1], round(dim(iris)[1]*70/100), replace = FALSE)
Training <- iris[TrainIndex,]
Testing <- iris[-TrainIndex,]

# Built-in data
data()
# to see all colours in R
colors()

# bring the data into the workspace
islands <- islands

# data = "islands"
# knowing structure of a data
str(islands)
write.table(islands,"islands.txt")

# specifying the subset
islands_1 <- islands
islands_2 <- islands[c(1, 2, 9, 27)]
islands_3 <- islands[-c(1, 2)]
islands_4 <- islands[islands < 20]
islands_5 <- islands[c("Timor","Sumatra")]
islands_6 <- mean(islands)
islands_7 <- islands[ifelse(islands>mean(islands), T,F)]

isle <- as.data.frame(islands)
isle$Name <- rownames(isle)
rownames(isle) <- NULL


# data = "iris"
# structure of a data
str(iris)
# subsetting data frames
# from 5th to 10th row
iris[5:10,]
# only "Sepal.Length" and "Sepal.Width" column
iris[,c("Sepal.Length", "Sepal.Width")]

# colnames
colnames(car)

# histogram plot
hist(car$Price)
hist(car$Price, col = "green")

hist(car$Price, seq(0, 80000, 10000),
     col = "red", 
     main = "Histogram",
     xlab = "Price")

# barplot of categorical data
barplot(table(car$Make), col="red")
unique(car$Make)
table(car$Make)
colors()

# pie chart
pie((table(car$Make)))
pie(table(car$Type))

# histogram plot
# absolute frequency
hist(car$Price, seq(0, 80000, 5000), # seq(starting point, end point of range/distribution, interval)
     prob = F, #probability = F -not mandatory to include
     ylim = c(0, 250), # y scale (two-element numeric vector with min and max values)
     xlim = c(0, 60000), # x scale (two-element numeric vector with min and max values)
     main = "Histogram",
     xlab = "Price")
# density
hist(car$Price, seq(0, 80000, 10000), 
     prob = T, # probability = T for density plot i.e. relative frequency (probability of having the frequency)
     #ylim = c(0, 0.00006), # scale of y in density 
     #xlim = c(0, 70000),
     main = "Histogram",
     xlab = "Mileage",
     ylab = "Density on Y",col=rainbow(8))
# Smooth normal distribution plot
curve(dnorm(x, 
            mean=mean(car$Price), 
            sd=sd(car$Price)), 
      add=TRUE, 
      col="red", 
      lwd=1)
# smooth density plot
lines(density(car$Price),col="red", lwd=1)

cap <- read.csv("C:/Brio Training/AEC/capcsv.csv")
hist(cap$Torque, seq(0, 40, 2), prob = T, col = "green")

# scatter plot
plot(Price ~ Mileage, 
     data = car, 
     type = "p", # plotting the data with points or dots
     col = "red", 
     main = "Scatter Plot")

plot(Price ~ Mileage, 
     data = car, 
     type = "p", # plotting the data with points or dots
     col = as.factor(car$Make), # if in data variable is not in factor form, use the variable directly 
     main = "Scatter Plot")

# Best fit line plot (lm - linear model, abline - best fit line)
abline(lm(Price ~ Mileage, 
          data = car))
abline(lm(Price ~ Mileage,
          data = car[car$Make == "Chevrolet",]),col="red")
View(car[car$Make =="Chevrolet",1:2]) # picks up price and mileage variable against the previous line where it picks up all variables in the rows

unique(car$Make) # for unique values in a category
# plotting symbols: "pch"
plot(Price ~ Mileage, data = car, type = "p", col="red", pch = 176)

carCluster<- read.csv("C:/Brio Training/AEC/CarCluster.csv")

plot(MPG ~ Horsepower, 
     data = carCluster, 
     type = "p", # plotting the data with points or dots
     col = "blue", 
     main = "Scatter Plot")

View(carCluster[carCluster$Country =="U.S.",3:6]) # picks up column numbers 3 to 6 where Country = U.S.
      
# scatter plot between flash recovery and Volts after in camera battery
batteriesCSV<- read.csv("C:/Brio Training/AEC/batteriesCSV.csv")
plot(FlashRecov ~ VoltsAfter, 
     data = batteriesCSV, 
     type = "p", # plotting the data with points or dots
     main = "Scatter Plot")

# box plot
boxplot(car$Price)
boxplot(car$Price~car$Type)

carpetCSV<- read.csv("C:/Brio Training/AEC/carpetCSV.csv")
boxplot(carpetCSV$Durability~carpetCSV$Carpet)


# Descriptive Statistics and Hypothesis Testing
# statistical functions

# Mean of a Variable
mean(car$Price)
mean(car$Price, trim=0.05) # mean after excluding 5% extreme obs - 2.5 % fm both sides

# Median of a Variable
median(car$Price)

# Mode of a variable (no built-in function)
getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
#match function: return revised record index vector after counting duplicate matched Mileage, once
#tabulate function: returns the same vector in terms of number of counts of each matched Mileage
#max function: returns the maximum count(s) from the vector
#which.max returns the row number of the max recosummard(s)
v <- car$Mileage
# Calculate the mode using the user function.
result <- getmode(v)
print(result)