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
some.states <- some.states[1:10, c("state.region", "Population", "Income")]

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
data1<-read.csv("D:/R/Class2/data1.csv")
#Create subset where Col1 starts with AT
dataAT <- subset(data1, grepl("^AT", Col1),Col1:Col5)

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

# specifying the subset
islands_1 <- islands
islands_2 <- islands[c(1, 2, 9, 27)]
islands_3 <- islands[-c(1, 2)]
islands_4 <- islands[islands < 20]
islands_5 <- islands[c("Timor","Sumatra")]
islands_6 <- mean(islands)
islands_7 <- islands[ifelse(islands>mean(islands), T,F)]


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

# box plot
boxplot(car$Price)
boxplot(car$Price~car$Type)

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
     ylim = c(0, 0.00006), # scale of y in density
     xlim = c(0, 70000),
     main = "Histogram",
     xlab = "Mileage",
     ylab = "Density on Y",col=rainbow(10))
# Smooth normal distribution plot
curve(dnorm(x,
            mean=mean(car$Price),
            sd=sd(car$Price)),
      add=TRUE,
      col="red",
      lwd=1)
# smooth density plot
lines(density(car$Price),col="red", lwd=1)

# scatter plot
plot(Price ~ Mileage,
     data = car,
     type = "p", # plotting the data with points or dots
     col = rainbow(10),
     main = "Scatter Plot")

plot(Price ~ Mileage,
     data = car,
     type = "p", # plotting the data with points or dots
     col = as.factor(car$Make), # if in data variable is not in factor form, use the variable directly
     main = "Scatter Plot")

# regression line plot (lm - linear model, abline - best fit line)
abline(lm(Price ~ Mileage,
          data = car))
abline(lm(Price ~ Mileage,
          data = car[car$Make == "Chevrolet",]))
View(car[car$Make =="Chevrolet",1:2]) # picks up price and mileage variable against the previous line where it picks up all variables in the rows

unique(car$Make) # for unique values in a category
# plotting symbols: "pch"
plot(Price ~ Mileage, data = car, type = "p", pch = 176)


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
#which.max returns the row number of the max record(s)
v <- car$Mileage
# Calculate the mode using the user function.
result <- getmode(v)
print(result)

# Range of a variable
range(car$Price)
fivenum(car$Price) #Tukey's five number summary (minimum, lower-hinge, median, upper-hinge, maximum)
# for summary of descriptive stat
summary(car$Price)

# Standard Deviation of a variable
sd(car$Price)

# mean by different category
tapply(car$Price, car$Make, mean)
tapply(car$Price, car$Model, median)
tapply(car$Price, car$Type, sd)
# to see unique values of character vector or variable
unique(car$Make)

# package install
# install "psych"
install.packages("psych")

# available packages
library()

# load a package
library(psych)
# or
require(psych)

# or import using menu option 'Packages'

# use the functions from "psych"
describe(car$Price) # mad: median absolute deviation (from the median)
kurtosi(car$Price)
describe(car) # for all variables in dataframe but not character variables
# correlation
cor(car$Price, car$Mileage, method = "pearson")

# parametric tests

# test for normality
shapiro.test(car$Price) #Test statistic W is judged by p-value
install.packages("nortest")
require(nortest)
# Anderson-Darling test for normality for large samples
ad.test(car$Price) # #Test statistic A is judged by p-value. price doesnt follow normality

# one-sample t-test - assuming normality for Price
# H0 : mean(Price) = 20000, H1 : Price not equal to 20000
t.test(car$Price, mu=20000, alternative="two.sided") #Upper CL Mean and Lower CL Mean displayed
t.test(car$Price, mu=20000, alternative="greater") #alternative hypothesis is that mu is > mean
t.test(car$Price, mu=20000, alternative="less") #alternative hypothesis is that mu is < mean

# Independent Sample t-test (continuous numeric variable and categorical variable has two levels)
# H0 : mean(Price) remain same for 2 door and 4 door car

#What is class of carmake
carmake = car$Make %in% c("Cadillac","Saturn")

t.test(Price ~ Doors, data = car,
       var.equal = FALSE,
       alternative = "two.sided",
       conf.level = 0.95
) #Difference of lower CL means and upper CL means are displayed

t.test(Price ~ Doors, data = car,
       var.equal = FALSE,
       alternative = "two.sided",
       conf.level = 0.99,
       subset = Make %in% c("Cadillac","Saturn")
)

car$Type <-as.character(car$Type)
car$Make <-as.character(car$Make)
t.test(Price ~ Doors, data = car[car$Make=="Cadillac" & car$Type=="Sedan",],
       var.equal = FALSE,
       alternative = "two.sided",
       conf.level = 0.99
) #Error: Data not valid, only one level found
car[car$Make=="Cadillac" & car$Type=="Sedan",]

# HOV test
# H0V(Price) for 2 door and 4 door car
# Bartlett's test
bartlett.test(Price ~ Doors, data = car) #Significance of test statistics (K^2) is judged by the p-value
# Fligner-Killeen test
fligner.test(Price ~ Doors, data = car)

# one way ANOVA
# H0 : mean(Price) remain same for all types of car
# types of car:
table(car$Type)
anova(lm(Price ~ Type, data = car))

# HoV and Anova for Carpet Durability
bartlett.test(Durability ~ Carpet, data = carpetCSV)
with(carpetCSV, levene.test(Durability, Carpet))
anova(lm(Durability ~ Carpet, data = carpetCSV))

# non-parametric tests
# u-test
# Mann-Whitney U test or Wilcoxon sign rank test (Alterntive to paired t-test when the normality test fails)
shapiro.test(Measurement$Measurment)
tapply(Measurement$Measurment, Measurement$Type, shapiro.test)
wilcox.test(Measurment~Type, data=Measurement)
kruskal.test(Measurment~Type, data=Measurement)

# tests for association (for categorical variables)
# tabulation
MakeVsType <- table(car$Make, car$Type)
# Chi-Test (H0 = attributes are independent)
summary(MakeVsType)

# correlation test (H0 = r = 0 i.e. they are not co-related)
cor.test(car$Price, car$Mileage,
         alternative = "two.sided",
         conf.level = 0.95)
# inference is they are co-related as p<0.05 and for degree of association,you'll have to check the actual value of r.


set.seed(1) # fixing the sample everytime
TrainIndex <- sample(1:dim(car)[1], round(dim(car)[1]*60/100), replace = F)
# replace = F for sampling without replacement
TrainIndex
Training <- car[TrainIndex,] # creating Training data set
Test <- car[-TrainIndex,] # Creating TEst dataset

m.1 <-lm(Price ~ Mileage, data = Training)
summary(m.1)
predict(m.1, Test)
cor(Test$Price,predict(m.1, Test))

#Four graphs that are useful for evaluating the model fit
fit <- lm(weight ~ height, data=women)
par(mfrow=c(2,2))
plot(fit)
#The par(mfrow=c(2,2)) statement is used
#to combine the four plots produced by the plot() function into one large 2x2 graph
#Q-Q: Normality
#Upper left: Residuals vs Predicted should be uncorrelated
#Bottom left: Homoscedasticity or constant variance of residuals
#Bottom right: To see outliers, influencial observation etc.

#Multiple linear regression
states <- as.data.frame(state.x77[,c("Murder", "Population",
                                     "Illiteracy", "Income", "Frost")])
cor(states) #See the relationship two at a time (bivariate correlation)

fit <- lm(Murder ~ Population + Illiteracy + Income + Frost,
          data=states)
summary(fit)

#Multiple linear regression #2
vars <- c("MPG","Weight","Horsepower","Drive_Ratio")
mtcars<-read.csv("D:/R/Class 5/CarCluster.csv")
head(mtcars[vars]) #First 6 observations
describe(mtcars[vars])

fit1 <- lm(MPG ~ Horsepower, data=mtcars)
summary(fit1)
plot(mtcars$Horsepower,studres(fit1))

#Multiple linear regression #2
fit <- lm(MPG ~ Horsepower + Weight + Drive_Ratio, data=mtcars)
summary(fit)

#Non-normality of residuals
library(MASS)
sresid <- studres(fit)
hist(sresid, freq=FALSE,
     main="Distribution of Studentized Residuals")
xfit<-seq(min(sresid),max(sresid),length=40)
yfit<-dnorm(xfit)
lines(xfit, yfit)

#Non-constant error variance test
library(car)
ncvTest(fit)
# plot studentized residuals vs. fitted values
spreadLevelPlot(fit)

#Multi-collinearity
vif(fit) # variance inflation factors
sqrt(vif(fit)) > 2 # problem?


#Aggregate on Mean, SD by Cylinder
aggregate(mtcars[vars], by=list(cyl=mtcars$Cylinders), mean)
aggregate(mtcars[vars], by=list(cyl=mtcars$Cylinders), sd)

#Use of describeBy from Psych
describeBy(mtcars[vars], mtcars$Cylinders)

# Regression analysis with categorical variables
summary(lm(Price ~ Mileage + factor(Make) + factor(Type),data=car))

# logistic regression (generalised linear model)
animaltest <-glm (decision ~ gender, data = Logistic, family="binomial")
summary (animaltest)
# ODDS = Exp(a+bX)
# ODDS ratio = Exp(b)
# Odds ratio
exp(coef(animaltest))

# Deviance is a measure of quality of fit of a generalized linear model.
# Or rather, it's a measure of badness of fit - higher numbers indicate worse fit.
# The null deviance shows how well the response is predicted by the model with nothing but an intercept.
# The residual deviance shows how well the response is predicted by the model when the predictors are included.
# It measures how poorly the model predicts the decisions - the smaller the statistic the better the model.

# The residual deviance is 399.91 on 313 degrees of freedom. We use this to test the overall fit of the model
# by once again treating this as a chi square value.
# p-value = 1 - pchisq(residual deviance, df) comes to 0.0006. The null hypothesis (i.e., the model) is rejected.
# The fitted values are significantly different from the observed values.

# We can also use the measure of reduction in deviance to test how the model performs overall:
# (Ho: No significant reduction in deviance, Ha: Significant reduction in deviance)
# In this example reduction in deviance on adding the predictor variable Gender is 425.57 - 399.91 = 25.66
# with loss of 1 degree of freedom.
# pvalue = 1 - pchisq(reduction of deviance, reduction of df) #if pvalue is >0.05, # then model appears to
# have performed poorly showing no significant reduction in deviance
# In this case, pvalue is < 0.05 showing that there is not enough reason to say that the model performs poorly


# The Akaike writeInformation Criterion (AIC) provides a method for assessing the quality of your model through
# comparison of related models.
# However, the number itself is not meaningful. If you have more than one similar candidate models (where all
# of the variables of the simpler model occur in the more complex models), then you should select the model
# that has the smallest AIC.

# Fisher's scoring algorithm is a derivative of Newton's method for solving maximum likelihood problems
# numerically. The generalised linear model uses the Fisher's scoring technique and we see that the algorithm
# needed 4 iterations to converge.

#Logistic regression - example 2

summary(binary)
# two-way contingency table of categorical outcome and predictors
## we want to make sure there are not 0 cells
xtabs(~ admit + rank, data = binary)

#First, we convert rank to a factor to indicate that rank should be treated as a
##categorical variable.
binary$rank <- factor(binary$rank)
admitlr <-glm (admit ~ gre + gpa + rank, data = binary, family="binomial")
summary (admitlr)
# Odds ratio
exp(coef(admitlr))

#The logistic regression coefficients (Beta) give the change in the log odds of the
#outcome for a one unit change in predictor variable.
#For every one unit change in gre, the log odds of admission (versus nonadmission)
#increases by 0.002 (Beta > 0)
#Having attended an undergraduate institution with a rank of 2. versus an
#institution with a rank of 1, decreases the log odds by 0.675443
#Odds ratio estimates -
#For a one unit increase in gpa, the odds of being admitted to graduate school
#(versus not being admitted) increases by a factor of 2.24


#One-way tables
install.packages("vcd")
library(vcd)
mytable <- with(Arthritis, table(Improved))
mytable
prop.table(mytable)
prop.table(mytable)*100 # you can use round function to trim post decimal places

#Two-way tables
mytable <- xtabs(~ Treatment+Improved, data=Arthritis)

chisq.test(mytable)


# clustering
CarCluster <-read.csv("D:/R/CarCluster.csv")
clus <-kmeans(CarCluster[,c("MPG","Weight", "Drive_Ratio")],3)

# starting with 3 clusters above
clus

# Within cluster sum of square: Sum of squared distances of each data point of the
# cluster to the cluster mean. So three values for cluster 1, 2 and 3
# Total_SS: Sum of squared distances of each data point to the global mean
# (global mean=mean of three final seeds)
# Between_SS: Sum of squared distances of three means to the global mean.
##When computing this, multiply the squared distance of each mean to the global mean
##by the number of data points it represents

CarCluster$c1<-clus$cluster

install.packages("scatterplot3d")
install.packages("rgl")
library("scatterplot3d")
scatterplot3d(CarCluster$MPG, CarCluster$Weight, CarCluster$Drive_Ratio,
              main = "Cluster Analysis",
              xlab = "MPG",
              ylab = "Weight",
              zlab = "Drive Ratio",
              pch = 20, # shape of dots in scatter plot
              type = "h", # for vertical lines to X-Y plane
              highlight.3d=TRUE
)

library(rgl)
open3d()
plot3d(x=CarCluster$MPG,
       y=CarCluster$Weight,
       z=CarCluster$Drive_Ratio,
       col=CarCluster$c1,
       xlab="MPG",
       ylab="Weight",
       zlab="Drive Ratio",
       size=5)

# Factor Analysis
library(foreign)
factor1 <- as.data.frame(read.csv("D:/R/factor1.csv"))
colnames(factor1)

?princomp

m.4 <- princomp(factor1[,c("ITEM13", "ITEM14", "ITEM15", "ITEM16", "ITEM17", "ITEM18", "ITEM19", "ITEM20", "ITEM21", "ITEM22", "ITEM23", "ITEM24")], cor=TRUE)

Or

m.4 <-princomp(factor1[,paste("ITEM",13:24,sep="")],cor=TRUE)

# cor=True is for taking it based on correlation matrix and not covariance matrix etc.

summary(m.4)
# SS loadings contains the eigenvalues associated with the components. The eigenvalues are standardized variance associated with a particular component

loadings(m.4)

# scree plot - showing the Eigenvalues of first 3 factors
windows(10,10)
plot(m.4,type="lines")

# rotate principal components using principal() function
install.packages("psych")
install.packages("GPArotation")
library(psych)
library(GPArotation)
m.5 <- principal(factor1[,c("ITEM13", "ITEM14", "ITEM15", "ITEM16", "ITEM17", "ITEM18", "ITEM19", "ITEM20", "ITEM21", "ITEM22", "ITEM23", "ITEM24")],
                 nfactors=3,
                 rotate="varimax",scores=TRUE)
#h2: component communalities- amount of variance in each variable accounted by the components
#u2: component uniquenesses- amount of variance not accounted by the components (1 - h2)

# print results
m.5
loadings(m.5)


#for loop in R
u1 <- rnorm(30) # create a vector filled with random normal values
print("This loop calculates the square of the first 10 elements of vector u1")

usq<-0
for(i in 1:10)
{
  usq[i]<-u1[i]*u1[i] # i-th element of u1 squared into i-th position of usq
  print(usq[i])
}
print(i)

#Nested for loop
# nested for: multiplication table
mymat = matrix(1:6,3,4) # create a 3 x 4 matrix (of 3 rows and 4 columns)

for(i in 1:dim(mymat)[1])  # for each row
{
  for(j in 1:dim(mymat)[2]) # for each column
  {
    mymat[i,j] = i*j     # assign values based on position: product of two indexes
    print(mymat[i,j])  #changed
  }
  #if (i==2){
    #break()
  #}
}

#if-else condition
mtcars <-read.csv("D:/R/CarCluster.csv")
for(i in 1:nrow(mtcars))

{
  if (mtcars$Horsepower[i]>100)
  {
    mtcars$new[i] <- "High HP"
  }
  else mtcars$new[i] <- "Low HP"
}

mtcars$new <- NULL #changed

#gsub() function replaces all matches of a string, if the parameter is a string vector
#gsub(pattern, replacement, x, ignore.case = FALSE)
# pattern: string to be matched
# replacement: string for replacement
# x: string or string vector
# ignore.case: if TRUE, ignore case
locn<-"Indore U.A."
locn<-gsub("\\b U.A.\\b", "",locn)
locn
#Number replacement
x <- "line 4322: He is now 25 years old, and weights 130.5lbs"
y <- gsub("\\d+","---",x)
y
#Lower case letters replacement
x <- "line 4322: He is now 25 years old, and weights 130lbs"
y <- gsub("[[:lower:]]","-",x)
y
#Vector replacement
x <- c("R Tutorial","PHP Tutorial", "HTML Tutorial")
gsub("Tutorial","Examples",x)


#Duplicated
x <- c(9:20, 1:5, 3:7, 0:8)
## extract unique elements
xu <- x[!duplicated(x)]
## similar, same elements but different order:(just don't print duplicates found from last) #changed
xu2 <- x[!duplicated(x, fromLast = TRUE)]

duplicated(iris)[140:143] /#Row duplicate TRUE or FALSE*/
  anyDuplicated(x)          #Position of First duplicate from start
anyDuplicated(x, fromLast = TRUE) #Position (counted from the first element) of First duplicate from last #changed


#Determinig the Sensitivity and Specificity in LR
training_length=ceiling(0.8*nrow(binary))
training=binary[1:training_length,]
test=binary[(training_length+1):(nrow(binary)),]

library(MASS)

binary_trg=data.frame(training)
binary_trg$rank <- factor(binary_trg$rank)
model1=glm(binary_trg,family=binomial,data=training)

fittedvalues=predict.glm(model1, newdata=test, type="response")
fittedvalueround=rep(0,length(fittedvalues))
for(i in 1:length(fittedvalues))
{
  if(fittedvalues[i]<=0.5)
    fittedvalueround[i]=0
  else
    fittedvalueround[i]=1
}
library(caret)
sensitivity=sensitivity(as.factor(fittedvalueround),as.factor(test$admit))
specificity=specificity(as.factor(fittedvalueround),as.factor(test$admit))
sensitivity <- paste(round(sensitivity, 4)*100, "%", sep="")
specificity <- paste(round(specificity, 4)*100, "%", sep="")


#TIME SERIES
#kings <- scan("http://robjhyndman.com/tsdldata/misc/kings.dat",skip=3) #age of death of 42 successive kings of England
kings <- scan("D:/R/TS/kings.dat",skip=3)
kingstimeseries <- ts(kings)
kingstimeseries

#births <- scan("http://robjhyndman.com/tsdldata/data/nybirths.dat")
births <- scan("D:/R/TS/nybirths.dat")
birthstimeseries <- ts(births, frequency=12, start=c(1946,1)) #1st quarter of 1946

#Plot the Time Series
plot.ts(kingstimeseries)

#Decomposing a time series means separating it into its constituent components,
#which are usually a trend component and an irregular component, and if it is a seasonal
#time series, a seasonal component.

# SIMPLE MOVING AVERAGE OF TIME SERIES
#To estimate the trend component of a non-seasonal time series that can be described
#using an additive model, it is common to use a smoothing method, such as calculating
#the simple moving average of the time series.

install.packages("TTR")
library("TTR")
kingstimeseriesSMA3 <- SMA(kingstimeseries,n=3) #smooth the time series using a simple moving average of order 3
plot.ts(kingstimeseriesSMA3)

#There still appears to be quite a lot of random fluctuations in the time series
#smoothed using a simple moving average of order 3. We can use higher order, say 8
kingstimeseriesSMA8 <- SMA(kingstimeseries,n=8)
plot.ts(kingstimeseriesSMA8)

# SIMPLE EXPONENTIAL SMOOTHING

#If you have a time series that can be described using an additive model with constant
#level and no seasonality, you can use simple exponential smoothing to make short-term
#forecasts.
#rain <- scan("http://robjhyndman.com/tsdldata/hurst/precip1.dat",skip=1) #total annual rainfall in inches for London, from 1813-1912
rain <- scan("D:/R/TS/precip1.dat",skip=1)
rainseries <- ts(rain,start=c(1813))
plot.ts(rainseries)

#To make forecasts using simple exponential smoothing in R, we can fit a simple
#exponential smoothing predictive model using the "HoltWinters()" function in R
#The HoltWinters() function returns a list variable, that contains several named elements.
rainseriesforecasts <- HoltWinters(rainseries, beta=FALSE, gamma=FALSE)
rainseriesforecasts

#In the example above, we have stored the output of the HoltWinters() function in the
#list variable "rainseriesforecasts". The forecasts made by HoltWinters() are stored
#in a named element of this list variable called "fitted"
rainseriesforecasts$fitted
plot(rainseriesforecasts)
#The plot shows the original time series in black, and the forecasts as a red line.

#sum of squared errors for the in-sample forecast errors
rainseriesforecasts$SSE

#We can make forecasts for further time points by using the "forecast.HoltWinters()"
#function in the R "forecast" package.
library("forecast")
rainseriesforecasts2 <- forecast.HoltWinters(rainseriesforecasts, h=8)
rainseriesforecasts2
#the forecasted rainfall for 1920 is about 24.68 inches, with a 95% prediction interval
#of (16.24, 33.11).
plot.forecast(rainseriesforecasts2)
#forecasts for 1913-1920 are plotted as a blue line, the 80% prediction interval as
#an deeper shaded area, and the 95% prediction interval as a lighter shaded area

#If you have a time series that can be described using an additive model with increasing
#or decreasing trend and no seasonality, you can use Holt's exponential smoothing to make
#short-term forecasts.
#skirts <- scan("http://robjhyndman.com/tsdldata/roberts/skirts.dat",skip=5) #Annual diameter of women's skirts at the hem, from 1866 to 1911
skirts <- scan("D:/R/TS/skirts.dat",skip=5)
skirtsseries <- ts(skirts,start=c(1866))
plot.ts(skirtsseries)
#Notice increase in hem diameter from about 600 in 1866 to about 1050 in 1880, and that
#afterwards the hem diameter decreased to about 520 in 1911

#To use HoltWinters() for Holt's exponential smoothing, we need to set the parameter
#gamma=FALSE

skirtsseriesforecasts <- HoltWinters(skirtsseries, gamma=FALSE)
skirtsseriesforecasts
#The estimated value of alpha is 0.84, and of beta is 1.00. These are both high, telling
#us that both the estimate of the current value of the level, and of the slope b of the
#trend component, are based mostly upon very recent observations in the time series.

plot(skirtsseriesforecasts)

#Predictions for 1912 to 1930 (19 more data points)
skirtsseriesforecasts2 <- forecast.HoltWinters(skirtsseriesforecasts, h=19)
plot.forecast(skirtsseriesforecasts2)
#The forecasts are shown as a blue line, with the 80% prediction intervals as deeper
#shaded area, and the 95% prediction intervals as a lighter shaded area.

#time series that can be described using an additive model with increasing
#or decreasing trend and Seasonality as well
retail=read.csv("D:/R/TIME_SERIES.csv",header=TRUE)
rev<-retail$Revenue
ats=ts(rev,start=c(2013),frequency=12) #frequency: No of obs per unit of time
plot(ats,type="l",col="green")
atsfct<-HoltWinters(ats) #Holt Winter's exponential smoothing method

atsfct2<-forecast(atsfct,h=3)
plot.forecast(atsfct2)
lines(atsfct2$fitted,col="red")
