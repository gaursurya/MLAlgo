###################  Kaggle Prudential Insurance Implementation on R  ###################################

rm(list=ls());
graphics.off();
options(warn=-1);

## Importing Files

df <- read_csv("C:/Users/Surya-PC/Desktop/Data Mining/Prudential.csv")

df <- df[sample(nrow(df)),]


## looking for NaN data
install.packages('Amelia')
library('Amelia')
## For Continous variable
df_Cont<-df[,names(df) %in% c("Product_Info_4","Ins_Age","Ht","Wt","BMI"
,"Employment_Info_1","Employment_Info_4","Employment_Info_6","Family_Hist_2","Family_Hist_3","Family_Hist_4",
 "Family_Hist_5","Medical_History_1")]
missmap(df_Cont, col=c("black", "grey"), legend=TRUE,
y.cex = .8,x.cex = .8,csvar = NULL)

# Data Processning

# Cutoff percentage for "Nan"
Cutoff=.25

#Column Which have "NaN"  More tha Cutoff Percentage
apply(df,2,function(col)sum(is.na(col))/length(col))>Cutoff

# Selecting columns which contain less than NaN as per Cut_Off percentage.
df<-df[which((apply(df, 2,
function(col)sum(is.na(col))/length(col))>Cutoff)==FALSE)]


### Converting int into Factor Level

df$Product_Info_1<-as.factor(df$Product_Info_1)
df$Product_Info_2<-as.factor(df$Product_Info_2)
df$Product_Info_3<-as.factor(df$Product_Info_3) 
df$Product_Info_5<-as.factor(df$Product_Info_5) 
df$Product_Info_6<-as.factor(df$Product_Info_6) 
df$Product_Info_7<-as.factor(df$Product_Info_7)
df$Employment_Info_2<-as.factor(df$Employment_Info_2)
df$Employment_Info_3<-as.factor(df$Employment_Info_3)
df$Employment_Info_5<-as.factor(df$Employment_Info_5) 
df$InsuredInfo_1<-as.factor(df$InsuredInfo_1)
df$InsuredInfo_2<-as.factor(df$InsuredInfo_2)
df$InsuredInfo_3<-as.factor(df$InsuredInfo_3)
df$InsuredInfo_4<-as.factor(df$InsuredInfo_4) 
df$InsuredInfo_5<-as.factor(df$InsuredInfo_5)
df$InsuredInfo_6<-as.factor(df$InsuredInfo_6)
df$InsuredInfo_7<-as.factor(df$InsuredInfo_7) 
df$Insurance_History_1<-as.factor(df$Insurance_History_1)
df$Insurance_History_2<-as.factor(df$Insurance_History_2) 
df$Insurance_History_3<-as.factor(df$Insurance_History_3)
df$Insurance_History_4<-as.factor(df$Insurance_History_4) 
df$Insurance_History_7<-as.factor(df$Insurance_History_7)
df$Insurance_History_8<-as.factor(df$Insurance_History_8)
df$Insurance_History_9<-as.factor(df$Insurance_History_9)
df$Response<-as.factor(df$Response)


#Mean Imputation



# Product_Info_4
df$Product_Info_4 <- ifelse(is.na(df$Product_Info_4), 
mean(df$Product_Info_4, na.rm=TRUE),df$Product_Info_4)

# Ins_Age
df$Ins_Age <- ifelse(is.na(df$Ins_Age), 
mean(df$Ins_Age, na.rm=TRUE),df$Ins_Age)

# Ht
df$Ht <- ifelse(is.na(df$Ht), 
  mean(df$Ht, na.rm=TRUE),df$Ht)

#wt
df$Wt <- ifelse(is.na(df$Wt), 
mean(df$Wt, na.rm=TRUE),df$Wt)

#BMI
df$BMI <- ifelse(is.na(df$BMI), 
mean(df$BMI, na.rm=TRUE),df$BMI)

#Employment_Info_1
df$Employment_Info_1 <- ifelse(is.na(df$Employment_Info_1 ), 
 mean(df$Employment_Info_1 , na.rm=TRUE),df$Employment_Info_1 )

#Employment_Info_4
df$Employment_Info_4 <- ifelse(is.na(df$Employment_Info_4), 
mean(df$Employment_Info_4, na.rm=TRUE),df$Employment_Info_4)

#Employment_Info_6
df$Employment_Info_6 <- ifelse(is.na(df$Employment_Info_6), 
 mean(df$Employment_Info_6, na.rm=TRUE),df$Employment_Info_6)

#Medical_History_1
df$Medical_History_1 <- ifelse(is.na(df$Medical_History_1), 
mean(df$Medical_History_1, na.rm=TRUE),df$Medical_History_1)


#### Working on Continous Data

########################### Continous/Nominal Variable DataFrame ##########################
df_Cont<-df[,names(df) %in% c("Product_Info_4","Ins_Age","Ht","Wt","BMI"
,"Employment_Info_1","Employment_Info_4","Employment_Info_6","Family_Hist_2","Family_Hist_3","Family_Hist_4",
"Family_Hist_5","Medical_History_1")]

#  dim(df_Cont)
## Boxplot For Continous Variables
boxplot(df_Cont$Product_Info_4,df_Cont$Ins_Age,df_Cont$Ht,df_Cont$Wt,df_Cont$BMI
        ,df_Cont$Employment_Info_1,df_Cont$Employment_Info_4,df_Cont$Employment_Info_6,
        df_Cont$Medical_History_1,col = c("#FF000099", "#FF6D0099"))

boxplot(df_Cont$Product_Info_4,df_Cont$Ins_Age,df_Cont$Ht,df_Cont$Wt,df_Cont$BMI
        ,df_Cont$Employment_Info_1,df_Cont$Employment_Info_4,df_Cont$Employment_Info_6,col = c("#FF000099", "#FF6D0099","blue","red"))

### Heat-Map of Continous Data

library(corrplot)
corrplot(cor(df_Cont), method = "circle")


##### wt and BmI having correlation of .85 so dropping one of it.#############


##### Dropping Id as It doen't make sense
df<-df[,!names(df) %in% c("BMI","Id")]


# Bar plot of Response Variable
library('ggplot2')
ggplot(df, aes(x=as.factor(df$Response) )) + geom_bar( fill=rgb(0.1,0.2,0.3,0.4
,.5,.6) )
ggplot(df, aes(x=as.factor(df$Product_Info_2) )) + geom_bar( fill=rgb(0.1,0.2,0.3,0.4
                                                                ,.5,.6) )
ggplot(df, aes(x=as.factor(df$Insurance_History_8) )) + geom_bar( fill=rgb(0.1,0.2,0.3,0.4
                                                                      ,.5,.6) )

##### Using Random Forest to do the Feature Selection with Default Parameter ######

library(randomForest)

rf_prud<-randomForest(Response~.,data=df,ntree=100)

# Variable Importance
varImpPlot(rf_prud,
           sort = T,
           n.var = 15,
           main = "Top 15 - Variable Importance")
importance(rf_prud)

### Selecting top 15 features after Features processning

df<-df[,names(df) %in% c("Wt","Ins_Age","Employment_Info_1","Product_Info_2",
                          "Product_Info_4","Ht","Medical_History_1","Medical_History_2"
                          ,"Employment_Info_6","InsuredInfo_3","Employment_Info_2",
                          "Employment_Info_4","Family_Hist_1","Insurance_History_8",
                          "Insurance_History_4","Response")]

lapply(df, class)

###########################CROSS-VALIDATION(75% Train set and 25% Test set)#######################

set.seed(1234)
index     <- 1:nrow(df)
testindex <- sample(index, trunc(length(index)/4))

dftrain  <- df[-testindex,]
dftest   <- df[testindex,]


########## Multinomial logistic Regression ##########
install.packages('nnet')
library(nnet)

model1<-multinom(formula=Response~.,data =dftrain)

# Predict the model

pred_model1<-predict(model1,dftest)

# Misclassification Error--0.5952

cm<-table(predict(model1,dftest),dftest$Response)
1-sum(diag(cm))/sum(cm)


#2-tailed z test
z <- summary(model1)$coefficients/summary(model1)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2
p


###Support Vector Machine Radial Kernel
library(e1071)


model2<-svm(Response~.,data=dftrain,kernel='radial',cross=5,gamma=.1,cost=10)

#Confusin Matrix and Misclassification Error

predict(model2,dftest)
cm_svm<-table(predicted=predict(model2,dftest),Actual=dftest$Response)


# Missclassification of Radial---0.5784
1-sum(diag(cm_svm))/sum(cm_svm)
## tuning of support vector Machine

#tuning of Radial Kernel
tune_model2<-tune.svm(Response~.,data=dftrain,gamma=seq(.1,1,.2),cost=seq(1,10,2))

plot(tune_model2)


tune_model2<-svm(Response~.,data=dftest,kernel='radial',cross=5,gamma=.1,cost=1)



cc<-table(predicted=predict(tune_model,dftest),Actual=dftest$Response)


1-sum(diag(cc))/sum(cc)


########### Support vector machine using Linear kernel  #######################

model3<-svm(Response~.,data=dftrain,kernel='linear',cross=5,cost=10)

#Confusin Matrix and Misclassification Error

cm_svm_linear<-table(predicted=predict(model3,dftest),Actual=dftest$Response)


# Missclassification of linear---0.5976
1-sum(diag(cm_svm_linear))/sum(cm_svm_linear)

#cost=seq(.1,1,.1)
#0.5729333 


## Tuning of Linear kernel

tune_model2<-tune.svm(Response~.,data=dftrain,cost=seq(2,20,3),cross=5)
plot(tune_model2)

## Model After tuning

model3<-svm(Response~.,data=dftest,kernel='linear',cross=5,cost=5)



#Confusin Matrix Error After Tuning

cm_tune<-table(predicted=predict(model3,dftest),Actual=dftest$Response)


# Missclassification rate After Tuning
1-sum(diag(cm_tune))/sum(cm_tune)



###### Random Forest #######

rf<-randomForest(Response~.,data = dftrain,mtree=500,mtry=3)
plot(rf)

#Confusin Matrix Error

cm<-table(predicted=predict(rf,dftest),Actual=dftest$Response)


# Missclassification rate ---- 0.5532
1-sum(diag(cm))/sum(cm)
seq(300,400,20)

## Tune Random forest

tune_rf <- tuneRF(y=dftrain$Response,x=dftrain[,-16],
            stepFactor = 1,
            plot = TRUE,
            ntreeTry =250,
            mtryStart = seq(19,35,3),
            trace = TRUE,
            improve = 1.5)


#### Random Forest After tuning

rf_tune<-randomForest(Response~.,data = dftrain,ntree=150,mtry=19)
plot(rf_tune)

#Confusin Matrix Error After tuning

cm_tune<-table(predicted=predict(rf_tune,dftest),Actual=dftest$Response)
# Missclassification rate After Tuning
1-sum(diag(cm_tune))/sum(cm_tune)


### Top 8 Features after tuning of model

# Variable Importance
varImpPlot(rf_tune,
           sort = T,
           n.var = 8,
           main = "Top 8 - Variable Importance")
importance(rf_tune)


