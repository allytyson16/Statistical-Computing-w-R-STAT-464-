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





 





