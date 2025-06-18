###ASSIGNMENT 1
# 1.(a) write an function that will generate grades based on the UG undergraduate grading system
#  (b) apply this on any 25 randomly generated marks.
# 2. (a) Do a similar for class categories based on GPA
#    (b) Generate 20 random GPA & apply the function to them


## QUESTION 1
# A- 80-100
# B+ - 75-79
# B - 70-74
# C+ = 65-69
# C - 60-64
# D+ - 55-59
# D - 50-54
# E - 40-49
# F - below

set.seed(123)
grad_system <- function(y)
{
  if (y >= 80)
  {
    Grade <- "A"
  }
  else if (y < 80 & y >= 75)
  {
    Grade <- "B+"
  }
  else if (y < 75 & y >= 70)
  {
    Grade <- "B"
  }
  else if (y < 70 & y >= 65)
  {
    Grade <- "C+"
  }
  else if (y < 65 & y >= 60)
  {
    Grade <- "C"
  }
  else if (y < 60 & y >= 55)
  {
    Grade <- "D+"
  }
  else if (y < 55 & y >= 50)
  {
    Grade <- "D"
  }
  else if (y < 50 & y >= 45)
  {
    Grade <- "E"
  }
  else
  {
    Grade <- "F"
  }
}

gen_marks <- sample(30:94, size = 25, replace = TRUE)
mark_grades <- sapply(gen_marks, grad_system)
scores <- data.frame(gen_marks, mark_grades)


## QUESTION 2
# First class - 3.60- 4.0
# Second Upper - 3.0-3.59
# Second Lower - 2.0 - 2.99
# Third  - below

set.seed(123)
gpa_cate <- function(w)
{
  if (w <= 4.0 & w >= 3.6)
  {
    gpa <- "First Class"
  }
  else if (w < 3.6 & w >= 3.0)
  {
    gpa <- "Second Class, Upper Division"
  }
  else if (w < 3.0 & w >= 2.0)
  {
    gpa <- "Second Class, Lower Division"
  }
  else if (w < 2.0 & w >= 1.5)
  {
    gpa <-"Third Class"
  }
  else if (w < 1.5 & 1.0)
  {
    gpa <- "Pass"
  }
  else
  {
    gpa <- "Fail, No Reward"
  }
}

gen_gpas <- signif(runif(20, min = 2, max = 4), 3)
gpas <- sapply(gen_gpas, gpa_cate)
gpa_class <- data.frame(gen_gpas, gpas)
head(gpa_class)
