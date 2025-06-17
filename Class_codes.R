#week 2: Introduction
# for assigning a variable use <- 
# for comment use harsh(#)
# for concatenation/ combine use c 
# for modulus %%
# for integer division %/%
# for natural log use log()
#for log to the base 10 use log(10)
# for assignment use ==

#Atomic vectors: one dimensional object that consists of elements of the same 
#                 data types(integer, double, logical, boolean, character, complex)

#how to create vectors
rm(list = ls())   #why?? what is the point of this?
my_num <- c(5, 1, 43, 5)
print(my_num)   #viewing object my_num

my_char <- c("Rice", "Fries", "Burger", "Yam", "Banku",
             "Fufu")                                  #characters/ strings/ texts
my_logi <- c(TRUE, FALSE, TRUE, TRUE, FALSE)          #logical
expenditure <- c(4000, 34232, NA, 34353, 90230, NA)   #use of NA(not applicable/ not available)
seq1 <- seq(0, 20, 4)                                         #seq(from,to, interval/step)
seq2 <- seq(0, 20, length.out = 11)                   #alternate way to create a sequence
rep1 <- rep(5, 4)

rep2 <- rep(c(2, 5, 7), 4)
rep3 <- rep(c(2, 5, 7), each=4)

#subset from vector
ages <- c(45, 34, 21, 54, 23)
names(ages) <- c("Kofi", "Ama", "Lexia", "Anna", "Emma", "Noah")

ages[c(1, 3)]
ages["Anna"]
ages[c("Ama", "Emma")]
ages["Noah"] <- 25
ages

length(ages)
mean(ages)
median(ages)
sd(ages)
var(ages)
max(ages)
min(ages)
range(ages)
sum(ages)
prod(ages)  #product, multiplies everything
quantile(ages, 0.25)

##conversions
#as.numeric -- to coerce to numeric vector
#as.character
#as.logical

Males <- c("Derrick", "Elikem", "Emma", "Andrews")
paste(Males, "likes fufu")

Females <- c("Oye", "Janet", "Aseye", "Zia")
paste(Males, "is always seen with", Females)

nchar(Males)   #number of characters
gsub("rr", "ll", Males)   #substitue any "rr" in Males with "ll"
grep("mm", Males, ignore.case = TRUE)

BP <- c(100, 102, 104, 120, 98, 105, NA)
mean(BP, na.rm = TRUE)
fivenum(BP)   #5 number summary:min, q1, median, q3, max
is.na(BP)     
sum(is.na(BP))    #to count the number of NA's in case we have large data


## MATRICES
Mat_A <- matrix(c(3, -1, 0, 5, 2, 1, 3, 4, 2), nrow = 3, ncol = 3, byrow = TRUE)
Row_1 <- c(1, 2, 3)
Row_2 <- c(4, 5, 6)
Row_3 <- c(7, 8, 9)
Mat_B <- rbind(Row_1, Row_2, Row_3)   #rbind is to bind rows

Row1 <- c(1, 4, 7)
Row2 <- c(2, 5, 8)
Row3 <- c(3, 6, 9)
MatB <- cbind(Row1, Row2, Row3)             #cbind is to bind columns

  
Age<- c(25, 20, 28)
weight <- c(1.67, 1.7, 1.74)
height <- c(64, 58, 70)
Mat_C <- cbind(Age, weight, height)             #cbind is to bind columns
rownames(Mat_C) <- c("Afua", "Kojo", "Atinga") #for naming rows
print(Mat_C)

Mat_D <- diag(c(1, 1, 1))                    #for diagonalising matrices (usually for identity matrix)
Mat_D[1,1]     #[row_number, column_number]
Mat_D[c(1, 2), c(2, 3)]
Mat_D[c(1:3), ]

Mat_C[,"Age"]
Mat_C["Kojo",]

##matrix operation
rowSums(Mat_A) 
dim(Mat_A)   #for dimension of matrix
colSums(Mat_A)

##matrix transpose
t(Mat_A)

#multiplication of matrices
Mat_A %*% Mat_B  #use %*% for matrix multiplication

#for determinant of matrix
det(Mat_A)

##solve systems of equations
Mat_1 <- matrix(c(1, 2, 4, 5), nrow = 2, byrow = TRUE)
Mat_2 <- matrix(c(3, 6), nrow = 2, byrow = TRUE)
solve(Mat_1, Mat_2)


###creating an array
array_1 <- array(data = seq(1:20), dim = c(4, 3, 2))

Gender <- c("M", "F", "M", "M", "F")
Gender <- factor(c("M", "F", "M", "M", "F"))
table(Gender)
prop.table(table(Gender))       ##probability of each gender


##creating a dataframe
Dat_fram1 <- data.frame(Name = c("Yaw", "Afia", "Aku"),
                        Gender = factor(c("Male", "Female", "Female")),
                        Ages = c(23, 21, 19),
                        GPA = c(2.01, 3.26, 1.98))
str(Dat_fram1)
Dat_fram1$Name
mean(Dat_fram1$GPA) 

Dat_fram1$Weight <- c(50.5, 54.0, 45)
head(Dat_fram1, 1)

###creating a list
List_1 <- list(Mat_1, Dat_fram1$Ages, Dat_fram1)
str(List_1)
List_1[1]   #displays the 2by2 matrix, since it is the first
length(List_1)
is.list(List_1)
is.list(Dat_fram1)  #how is this a data frame
is.data.frame(List_1)
is.data.frame(Dat_fram1)



###PACKAGES
install.packages("tidyverse")
install.packages("infer")
install.packages("fitdistrplus")
install.packages("actuar")
install.packages("car")
install.packages("haven")   #helps read data from SASS, 

##use library() to attach package
library(tidyverse)
library(infer)


##importing data in R
count_data <- read.csv(file = file.choose(), 
                       header = TRUE, sep = ",",
                       dec = ".")                   ## after running you get to select a csv file
str(count_data)

count_data1 <- read.table(file = "clipboard", header = TRUE, fill = TRUE)

##alternatively;
count_data2 <- read.csv(file = "file path\\ whatever\\whatever",
                        header = T,
                        sep = ",",
                        dec = ".")

###Simulating vectors
##categorical
set.seed(464)
Gend <- sample(c("Male", "Female"), size = 15, replace = T,
               prob = c(0.65, 0.35))      ##probabilities w/ male before females

Age <- round(rnorm(15, mean = 18, sd = 8), 2)          ##2 is for the dp for the round function
weight <- signif(runif(15, min = 40, max = 62), 4) ##4 is for significant figure for the signi function, also 2 dp

Data_1 <- data.frame(Educ_level = sample(c("None", "Basic", "Secondary", "Tertiary"), size = 40, 
                                         replace = TRUE, prob = c(0.2, 0.3, 0.35, 0.15)),
                                         Income = round(rnorm(40, mean = 1500, sd = 850), 2),
                                         Height = round(runif(40, min = 1.42, max = 1.90), 2))

head(Data_1, 10)  #output first 10 data
str(Data_1)
Data_1$Educ_level <- factor(Data_1$Educ_level, levels = c("None", "Basic", "Secondary", "Tertiary"))
str(Data_1$Educ_level)
write.csv(Data_1, file = "/Users/adelaideyeboah/Desktop/Stat464")   ##can you write a csv, if you have not imported data
                                                                      ## from that folder yet.




####functions in R
# 1. Assign it a name
# 2. argument
# 3. body
#### writing functions
Add_func <- function(x, y)
{
  b <- x + y
  return(b)
}

Add_func(x = 4, y = 6)  
Add_func(4, 6)
Add_func(x = 4, 6)

Derrick_function <- function(m = 2, n)
{
  p <- m^3 + 3*m + n
  return(p)
}
Derrick_function(4, 3)

##in functions if you want to use a default argument call them in the argument.
##you can assign another value to default arguments, it will over ride

greetings <- function(name = "Kumah", message = "Good day")
{
  greet <- paste(message, name)
  print(greet)
}

greetings(name = "Obolo")

Alge_func <- function(x, y)
{
  add <- x + y
  subt <- x - y
  sq_diff <- x^2 - y^2
  return(c(add, subt, sq_diff))
}

Alge_func(x = 5, y = 2)

##alternatively, we can give the rhem names
Alge_func <- function(x, y)
{
  add <- x + y
  subt <- x - y
  sq_diff <- x^2 - y^2
  res <- c(add, subt, sq_diff)
  names(res) <- c("Addition", "Subtraction", "Square Difference")
  return(res)
}

Alge_func(x = 5, y = 2)

##write functions for:
###quadratic
#cubic
#rational function
#mean
#standard deviation

##quadratic
Quad_func <- function(x)
{
  quad_func <- x**2 + x + 7
  print(quad_func)
}
Quad_func(2)

##cubic 
Cubic_func <- function(y)
{
  cube_func <- 2*y^3 - 5*y^2 + y -7
  print(cube_func)
}
Cubic_func(4)

##rational function
Rat_func <- function(m)
{
  numerator <- m^2 + 2*m +1
  denominator <- m -3
  ratnl_func <- numerator/ denominator
  return(ratnl_func)
}
Rat_func(4)


##mean





#######LOOPS
####FOR LOOPS
for (i in 1:10)
{
  print(i^2)
}

for (x in 1:5)
{
  for (y in 5:12)
  {
    print(Alge_func(x, y))
  }
}
##try and add a column for when i is a number then 
##we see the addition, subtraction & square difference
 
##WHILE LOOPS
numb <- 2
while (numb <= 30)
{
  numb <- numb + 2
  print(numb)
} 
##look at code again


##REPEATED loop 
#repeats till forever, so it is advisable to add the if command to stop it.
repeat
{
  print(numb)
  numb <- numb + 2
  if (numb >= 40)
    break
}
 

##NEXT STATEMENT
#USED IN LOOPS/ STATEMNETS TO STOP/ BREAK THE CODE.



#Apply loops, 
##Lapply function <-  used on either lists or vectors
lapply(1:5, function(x) x^2 - 4) 
#retuns a list, hence we have to unlist it
unlist(lapply(1:5, function(x) x^2 - 4))

##sapply <- returns a vector
sapply(1:5, function(x) x^2 - 4)

##apply <- used only on matrices and arrays
apply(Mat_A, 1, FUN = function(x) sum(x^2))   #1 is for row and 2- column
##squared the entries of the matrix & sum the rows 

##mapply <-  multivariate version of lapply
mapply(Alge_func, x = 1:5, y = 5:12)

##tapply <- usually on vectors


##CONDITIONAL STATEMENTS
x <- 4
if (x>=2)
{
  print("number is greater than 2")
}

func_1 <- function(x)
{
  if (x>= 2)
  {
    print("The number is greater than 2.")
  }
  else
  {
    print("The number is less than 2")
  }
}

vec_1 <- c(5, 1, 3, 7, 0)
func_1(vec_1[1])
func_1(vec_1[2])   #does not output anything because 1<2
func_1(vec_1[3])
func_1(vec_1[4])
func_1(vec_1[5])    #does not output anything because 1<2

##alternatively
vec_1 <- c(5, 1, 3, 7, 0)
sapply(vec_1, func_1)

##alternatively (best)
for (i in vec_1)
{
    func_1(i)
}

p_val_int <- function(p_val, alpha = 0.05)
{
  if (p_val >= alpha)
  {
    print("Fail to reject H0")
  }
  else{
    print("Reject H0")
  }
}

p_val_int(0.04)
##generating numbers as p_avlue
pvalue <- runif(10, 0, 1)  
for (i in pvalue)
{
  p_val_int(i)
}

p_val_int1 <- function(p_val, alpha = 0.05)
{
  ifelse(p_val >= alpha, "Fail to Reject H0", "Reject H0")
}
for (i in pvalue)
{
  print(p_val_int1(i))
}


###assignment score
# 80 + <-  excellent
# 70+  very good
# 50≠ <- average
# below <- poor

perf_eval <- function(x)
{
  if (x > 80 )
  {
    Perf <- "Excellent"
  }
  else if (x < 80 & x > 70)
  {
    Perf <- "Very Good"
  }
  else if (x < 70 & x > 50)
  {
    Perf <- "Average"
  }
  else
  {
    Perf <- "Poor performance"
  }
    return(c(x, Perf))
}

perf_eval(85)
sapply(c(72, 45, 56, 84), perf_eval)

###ASSIGNMENT 1
# 1.(a) write an function that will generate grades based on the UG undergraduate grading system
#  (b) apply this on any 25 randomly generated marks.
# 2. (a) Do a similar for class categories based on GPA
#    (b) Generate 20 random GPA & apply the function to them




###DATA VISUALIZATION

#HISTOGRM <-  uses the hist function 
hist(iris$Sepal.Length, col = "beige", xlab = "Sepal Length", main = "")

##BOXPLOT 
boxplot(iris$Sepal.Length, col = "skyblue")












