#exercise 1.1

View(attenu)
#find which lines have station data NA
which(is.na(attenu$station))
#clean attenu file with no NAs
attenu_cleaned <- attenu[!is.na(attenu$station),]
#first 6 lines + dimension of cleaned version
head(attenu_cleaned)
dim(attenu_cleaned)

#exercise 1.2
View(Theoph)
Theoph_2 <- Theoph
str(Theoph_2)
#find median dosage
med_dose <- median(Theoph_2$Dose)
#add new column Dose Class
Theoph_2$Dose_Class <- ifelse(Theoph_2$Dose >= med_dose, "high", "low")
#first 6 lines + dimension of new version
head(Theoph_2)
dim(Theoph_2)

#exercise 1.3
getwd()
setwd("/Users/echotang/qbio490/qbio_data_analysis_echo/week3_homework")
starbucks <- read.csv("starbucks.csv")
is.na(starbucks)
is_row_empty <- rowSums(is.na(starbucks[,2:7])) == 0
nrow(starbucks)
dim(starbucks)
length(is_row_empty)
#clean starbucks data
starbucks_cleaned <- starbucks[is_row_empty == TRUE,]
#compare calories + carbs
plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, xlab = "Carbs per drink in grams", ylab = "Calories per drink")
max(starbucks_cleaned$Calories) #max of 430 calories
which(starbucks_cleaned$Calories == max(starbucks_cleaned$Calories))
plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, xlab = "Carbs per drink in grams", ylab = "Calories per drink", col = ifelse(starbucks_cleaned$Calories == max(starbucks_cleaned$Calories), "red", "black"))
#Starbucks signature hot choclate

max(starbucks_cleaned$Fat)
which(starbucks_cleaned$Fat == max(starbucks_cleaned$Fat))
starbucks_cleaned$is_highest_fat <- ifelse(starbucks_cleaned$Fat == max(starbucks_cleaned$Fat), "TRUE", "FALSE")
plot(starbucks_cleaned$Carb, starbucks_cleaned$Calories, xlab = "Carbs per drink in grams", ylab = "Calories per drink", col = ifelse(starbucks_cleaned$Fat == max(starbucks_cleaned$Fat), "red", "black"))


#exercise 1.4
batting <- read.csv("batting.csv")
#players scoring 3 or more home runs
nrow(batting[batting$HR >= 3,])
plot(batting$yearID, batting$HR, xlab = "Year", ylab = "Number of home runs")

LA_Angels <- batting[batting$teamID == "LAA",]
plot(LA_Angels$yearID, LA_Angels$HR, xlab = "Year", ylab = "Number of home runs")
ATL_PIT <- batting[batting$teamID == "ATL" | batting$teamID == "PIT",]
plot(ATL_PIT$yearID, ATL_PIT$HR, xlab = "Year", ylab = "Number of home runs", col = ifelse(ATL_PIT$teamID == "ATL", "red", "black"))

#exercise 1.5
easy_plot <- function(x, y, color_data){
  med = median(color_data)
  levels = character(92)
  levels <- ifelse(color_data >= med, "high", "low")
  levels = factor(levels)
  plot(x, y, col = levels, pch = 20)
  print(cor.test(x,y))
}
#use easy_plot on starbucks_cleaned and batting data
easy_plot(starbucks_cleaned$Calories, starbucks_cleaned$Fat, starbucks_cleaned$Carb)

batting <- read.csv("batting.csv")
easy_plot(batting$G, batting$R, batting$H)
#All the variables are mostly positively correlated. Those labeled "high" usually had high x and y values and were colored black. Those labeled "low" had low x and y values and were colored red

#exercise 2.1
View(iris)
#The dataset describes the sepal length, sepal width, petal length, and petal width of different species of irises.
dim(iris)
#It contains 150*5 = 750 observations
#It has 5 variables per data point

#exercise 2.2
#There are 4 continuous variables and 1 categorical variable. Sepal length and width, and petal length and width are numeric, whereas species is character.

#exercise 2.3
hist(iris$Sepal.Length, xlab = "Sepal length", main = "Histogram of sepal length")
hist(iris$Sepal.Width, xlab = "Sepal width", main = "Histogram of sepal width")
hist(iris$Petal.Length, xlab = "Petal length", main = "Histogram of petal length")
hist(iris$Petal.Width, xlab = "Petal width", main = "Histogram of petal width")
#Sepal length has the lowest amount of bins. Both petal measurement data seem to be bimodal. 

#exercise 2.4
sepal_width_mean <- mean(iris$Sepal.Width)
iris_copy <- iris
sepal_comp = character(150)
sepal_comp = ifelse(iris_copy$Sepal.Width >= sepal_width_mean, "wide-sepaled", "narrow-sepaled")
#Create new column with that info
iris_copy$Sepal.Category <- sepal_comp
#boxplot with sepal width
boxplot(iris_copy$Sepal.Width ~ iris_copy$Sepal.Category, xlab = "Sepal width per sepal categorization", ylab = "Sepal width")

#exercise 2.5
#Versicolor and virginica look the most similar while setosa looks the most unique.
pairs(iris[1:4], pch = 21, col = c("red", "green3", "blue")[unclass(iris$Species)])



