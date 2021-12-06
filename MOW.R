# Autorzy:
# Katarzyna Muter
# Maciej Rosoł
# MOW - Metody Odkrywania Wiedzy

library(stringr)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(randomForest)
library(caret)
library(e1071)
library(rpart)
library(R.basic)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# rowData <- read.table("spliceATrainKIS.dat", header=TRUE, sep=",")

rowData <- read.table("spliceDTrainKIS.dat", header=TRUE, sep=",")

cname <- colnames(rowData)
rozciecie <- strtoi(substr(cname,2,nchar(cname)))

Data <- data.frame('Label' = rowData[cname][seq(1,dim(rowData)[1],2),1], 'Value' = rowData[cname][seq(2,dim(rowData)[1],2),1])
#Data <- data.frame('Label' = rowData[cname][seq(1,dim(rowData)[1],2),1], 'Value' = rowData[cname][seq(2,dim(rowData)[1],2),1])

litery <- c('A','C','T','G')
cechy = matrix(0, nrow = dim(Data)[1], ncol = 4*nchar(toString(Data$Value[1]))) # wystąpienie danej litery na danym miejcu. Macierz, gdzie każda kolumna oznacza wystąpienie danej litery
                                                    # na danej pozycji we fragmencie RNA. 4 litery * 15 pozycji = 60 kolumn. Liczba wierszy = liczbie próbek

for (i in 1:dim(Data)[1]){ 
  l_l = 0
  for (l in litery) {
    for (m in 1:nchar(toString(Data$Value[1]))){
      if (l == substr(toString(Data$Value[i]),m,m))
        cechy[i, m + l_l*nchar(toString(Data$Value[i]))] = 1
      else 
        cechy[i, m + l_l*nchar(toString(Data$Value[i]))] = 0
    }
    l_l = l_l + 1
  }
}

cechy1 = matrix(0, nrow = dim(Data)[1], ncol = 4) # łączna liczba wystąpień danych liter we fragmencie kodu RNA
for (i in 1:dim(Data)[1]){ 
  for (l in 1:4) {
     cechy1[i, l] = str_count(toString(Data$Value[i]),litery[l])
    }
}

cechy2 = matrix(0, nrow = dim(Data)[1], ncol = 4) # procent danej litery po lewej stronie miejsca rozcięcia
for (i in 1:dim(Data)[1]){
  for (l in 1:4){
    cechy2[i,l] = str_count(substr(toString(Data$Value[i]),1,rozciecie),litery[l])/nchar(substr(toString(Data$Value[i]),1,rozciecie))
  }
}

cechy3 = matrix(0, nrow = dim(Data)[1], ncol = 4) # procent danej litery po prawej stronie miejsca rozcięcia
for (i in 1:dim(Data)[1]){
  for (l in 1:4){
    cechy3[i,l] = str_count(substr(toString(Data$Value[i]),rozciecie,nchar(toString(Data$Value[1]))),litery[l])/nchar(substr(toString(Data$Value[i]),rozciecie,nchar(toString(Data$Value[1]))))
  }
}

cechy4 = matrix(0, nrow = dim(Data)[1], ncol = nchar(toString(rowData[2,1]))*(nchar(toString(rowData[2,1]))-1)*4*4) # jendal litera na jednym miejscu i druga na drugim
cechy5 = matrix(0, nrow = dim(Data)[1], ncol = nchar(toString(rowData[2,1]))*(nchar(toString(rowData[2,1]))-1)*4*4) # jendal litera na jednym miejscu lub druga na drugim

for (p in 1:dim(Data)[1]){
  k = 1
  for (i in 1:nchar(toString(rowData[2,1]))){
    for (j in 1:nchar(toString(rowData[2,1]))){
      if (i!=j){
        
        for(l1 in litery){
          for(l2 in litery){
            
            if(substr(toString(Data$Value[p]),i,i)==l1 & substr(toString(Data$Value[p]),j,j)==l2){
              cechy4[p,k] = 1
            }
            if(substr(toString(Data$Value[p]),i,i)==l1 || substr(toString(Data$Value[p]),j,j)==l2){
              cechy5[p,k] = 1
            }
            k = k+1
          }
        }
      }
    }
  }
}
# connecting features -----------------------------------------------------

features <- cbind(cechy,cechy1,cechy2,cechy3,cechy4,cechy5)
# features <- cbind(cechy,cechy1,cechy2,cechy3)

y <- Data$Label
y<-as.numeric(levels(y))[y]
features <- zscore(features)

# Pearson correlation -----------------------------------------------------

corr <- numeric(dim(features)[2])
for (i in 1:dim(features)[2]){
  corr[i] <- cor(features[,i], y, method = c("pearson"))
  
}
corr <- abs(corr)

sP <- sort(corr, partial = NULL, na.last = NA, decreasing = TRUE,
        method = c("auto", "shell", "quick", "radix"), index.return = TRUE)
selected_pearson <- features[,sP$ix[1:ceiling(dim(features)[2]/1)]]
# selected_pearson <- features[,sP$ix[1:ceiling(dim(features)[2]/10)]]
# selected_pearson <- features[,sP$ix[1:ceiling(dim(features)[2]*0.003)]]

# F-score -----------------------------------------------------------------

mean0 <- colMeans(features[y==0,])
mean1 <- colMeans(features[y==1,])

var0 <- colVars(features[y==0,])
var1 <- colVars(features[y==1,])

F_score <- (mean0-mean1)**2/(var0+var1)

sF <- sort(F_score, partial = NULL, na.last = NA, decreasing = TRUE,
          method = c("auto", "shell", "quick", "radix"), index.return = TRUE)
selected_Fscore <- features[,sF$ix[1:ceiling(dim(features)[2]/1)]]
# selected_Fscore <- features[,sF$ix[1:ceiling(dim(features)[2]/10)]]
# selected_Fscore <- features[,sF$ix[1:ceiling(dim(features)[2]*0.003)]]

# T-test ------------------------------------------------------------------

p_values <- numeric(dim(features)[2])
for (i in 1:dim(features)[2]){
  a<-1
  tryCatch({ a<- t.test(features[y==0,i],features[y==1,i])$p.value},
           error= function(cond) {a <- 1},
           warning = function(war) {a <- 1 },
           finally =  {p_values[i] <- a}
           )
}

sT <- sort(p_values, partial = NULL, na.last = NA, decreasing = FALSE,
           method = c("auto", "shell", "quick", "radix"), index.return = TRUE)
selected_ttest <- features[,sT$ix[1:ceiling(dim(features)[2]/1)]]
# selected_ttest <- features[,sT$ix[1:ceiling(dim(features)[2]/10)]]
# selected_ttest <- features[,sT$ix[1:ceiling(dim(features)[2]*0.003)]]

# Random Forest -----------------------------------------------------------

# ustawienie seed, żeby uzyskać powtarzalność analiz
set.seed(123)
train_index <- sample(1:nrow(features), 0.8 * nrow(features)) # wybranie indeksów danych trenujących
test_index <- setdiff(1:nrow(features), train_index) # wybranie indeksów danych testujących

# PEARSON

# Wybranie danych trenujących i testujących dla cech wybranych według współczynnika korelacji Pearsona
xTrainPearson <- selected_pearson[train_index,] 
xTestPearson <- selected_pearson[test_index,]

yTrain <- y[train_index]
yTest <- y[test_index]

yTrain <- factor(yTrain)
yTest <- factor(yTest)

acc = 0
sens = 0
spec = 0

TP = 0
TN = 0
FP = 0
FN = 0

for (i in 1:5){ 
  RF <- randomForest(data.frame(xTrainPearson), y=data.frame(yTrain)) # utworzenie modelu lasu losowego
  yhatTest = predict(RF, xTestPearson) # wyliczenie etykiet dla wartości testowych
  yhatTrain = predict(RF, xTrainPearson) # wyliczenie etykiet dla wartości trenujących
  
  cmRF_Pearson = confusionMatrix(yTest, yhatTest)
  acc = acc + cmRF_Pearson$overall[1]
  sens = sens + cmRF_Pearson$byClass[1]
  spec = spec + cmRF_Pearson$byClass[2]
  
  TP = TP + cmRF_Pearson$table[4]
  TN = TN + cmRF_Pearson$table[1]
  FP = FP + cmRF_Pearson$table[2]
  FN = FN + cmRF_Pearson$table[3]
#  cmRF_Pearson_test = confusionMatrix(yTrain, yhatTrain)
  print(cmRF_Pearson)  # macierz błędu danych testowych
#  print(cmRF_Pearson_test)  # macierz błędu danych trenujących
}

accPRF = acc/5
sensPRF = sens/5
specPRF = spec/5

TP_P = round(TP/5)
TN_P = round(TN/5)
FN_P = round(FN/5)
FP_P = round(FP/5)

# F-SCORE

xTrainFscore <- selected_Fscore[train_index,]
xTestFscore <- selected_Fscore[test_index,]

acc = 0
sens = 0
spec = 0

TP = 0
TN = 0
FP = 0
FN = 0

for (i in 1:5){
  
  RF <- randomForest(xTrainFscore, y=yTrain)
  yhatTest = predict(RF, xTestFscore)
  yhatTrain = predict(RF, xTrainFscore)

  cmRF_Fscore = confusionMatrix(yTest, yhatTest)
  acc = acc + cmRF_Fscore$overall[1]
  sens = sens + cmRF_Fscore$byClass[1]
  spec = spec + cmRF_Fscore$byClass[2]
  
  TP = TP + cmRF_Fscore$table[4]
  TN = TN + cmRF_Fscore$table[1]
  FP = FP + cmRF_Fscore$table[2]
  FN = FN + cmRF_Fscore$table[3]
# cmRF_Fscore_test = confusionMatrix(yTrain, yhatTrain)
print(cmRF_Fscore)  
# print(cmRF_Fscore_test)  
}
accFRF = acc/5
sensFRF = sens/5
specFRF = spec/5

TP_F = round(TP/5)
TN_F = round(TN/5)
FN_F = round(FN/5)
FP_F = round(FP/5)

# t-test

xTrain_ttest <- selected_ttest[train_index,]
xTest_ttest <- selected_ttest[test_index,]
acc = 0
sens = 0
spec = 0

TP = 0
TN = 0
FP = 0
FN = 0

for (i in 1:5){
  RF <- randomForest(xTrain_ttest, y=yTrain)
  yhatTest = predict(RF, xTest_ttest)
  yhatTrain = predict(RF, xTrain_ttest)

  cmRF_ttest = confusionMatrix(yTest, yhatTest)
  acc = acc + cmRF_ttest$overall[1]
  sens = sens + cmRF_ttest$byClass[1]
  spec = spec + cmRF_ttest$byClass[2]

  TP = TP + cmRF_ttest$table[4]
  TN = TN + cmRF_ttest$table[1]
  FP = FP + cmRF_ttest$table[2]
  FN = FN + cmRF_ttest$table[3]
# cmRF_ttest_test = confusionMatrix(yTrain, yhatTrain)
print(cmRF_ttest)  
# print(cmRF_ttest_test) 
}
accTRF = acc/5
sensTRF = sens/5
specTRF = spec/5

TP_T = round(TP/5)
TN_T = round(TN/5)
FN_T = round(FN/5)
FP_T = round(FP/5)

# SVM ---------------------------------------------------------------------

# PEARSON

svmModel = svm(xTrainPearson, yTrain, , kernel = 'polynomial', degree = 5)
yhatTest = predict(svmModel, xTestPearson)
yhatTrain = predict(svmModel, xTrainPearson)

cmSVM_Pearson = confusionMatrix(yTest, yhatTest)
cmSVM_Pearson_test = confusionMatrix(yTrain, yhatTrain)
print(cmSVM_Pearson)  
print(cmSVM_Pearson_test) 

# F-SCORE

svmModel = svm(xTrainFscore, yTrain, , kernel = 'polynomial', degree = 5)
yhatTest = predict(svmModel, xTestFscore)
yhatTrain = predict(svmModel, xTrainFscore)

cmSVM_Fscore = confusionMatrix(yTest, yhatTest)
cmSVM_Fscore_test = confusionMatrix(yTrain, yhatTrain)
print(cmSVM_Fscore)  
print(cmSVM_Fscore_test)  

# t-test

svmModel = svm(xTrain_ttest, yTrain, , kernel = 'polynomial', degree = 5)
yhatTest = predict(svmModel, xTest_ttest)
yhatTrain = predict(svmModel, xTrain_ttest)

cmSVM_ttest = confusionMatrix(yTest, yhatTest)
cmSVM_ttest_test = confusionMatrix(yTrain, yhatTrain)
print(cmSVM_ttest)  
print(cmSVM_ttest_test) 

# Decision Tree -----------------------------------------------------------

# PEARSON

Tree <- rpart(yTrain~., data = data.frame(xTrainPearson), method = 'class')
yhatTest = predict(Tree, data.frame(xTestPearson),type = 'class')
yhatTrain = predict(Tree, data.frame(xTrainPearson),type = 'class')

cmDT_Pearson = confusionMatrix(yTest, yhatTest)
cmDT_Pearson_test = confusionMatrix(yTrain, yhatTrain)
print(cmDT_Pearson)  
print(cmDT_Pearson_test)

# F-SCORE

Tree <- rpart(yTrain~., data = data.frame(xTrainFscore), method = 'class')
yhatTest = predict(Tree, data.frame(xTestFscore),type = 'class')
yhatTrain = predict(Tree, data.frame(xTrainFscore),type = 'class')

cmDT_Fscore = confusionMatrix(yTest, yhatTest)
cmDT_Fscore_test = confusionMatrix(yTrain, yhatTrain)
print(cmDT_Fscore)  
print(cmDT_Fscore_test) 

# t-test

Tree <- rpart(yTrain~., data = data.frame(xTrain_ttest), method = 'class')
yhatTest = predict(Tree, data.frame(xTest_ttest),type = 'class')
yhatTrain = predict(Tree, data.frame(xTrain_ttest),type = 'class')

cmDT_ttest = confusionMatrix(yTest, yhatTest)
cmDT_ttest_test = confusionMatrix(yTrain, yhatTrain)
print(cmDT_ttest)  
print(cmDT_ttest_test) 
