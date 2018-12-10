###################################
######### Expression AUC ##########
###################################

rm(list=ls());
graphics.off();
options(warn=-1);

dir_data = "C:/Users/we007904/Desktop/test_Drug/Data";
DataFileName = "RNASeq_log.csv";
CellLineName = "CellLineName.csv";
DrugResponse = "Drug_AUC.csv";
LimitFile = "AUC_Limit.csv";

### Packages for the models

library(hydroGOF); # MSE
library(pROC);
library(pls); # PLSR
######## Parameters
var_cut = 0.5;
mean_cut= 0.4;
Feature_Num = 500;
rr = 1000;

######## Data preparation
#1. Input Data
Data = as.matrix(read.delim(paste(dir_data,"/", DataFileName,sep=""),sep=",",header=F));
Data = t(Data);
colnames(Data) = Data[1,];
Data = Data[-1,];
rownames(Data) = Data[,1];
Data = Data[,-1];
class(Data) <- "numeric"

idx = c();
for (i in 1:dim(Data)[2]){
  ix = which(is.nan(Data[,i]));
  ix1 = which(is.na(Data[,i]));
  ix = union(ix,ix1);
  if(length(ix)>dim(Data)[1]*0.1){
    idx =  rbind(idx,c(i)) ;
  }else{
    Data[ix,i] = mean(Data[-ix,i]);
  }
}
idx = as.vector(idx);
Data = Data[,-idx];

#2. Drug response
Drug = as.matrix(read.delim(paste(dir_data,"/", DrugResponse,sep=""),sep=",",header=F));
colnames(Drug) = Drug[1,];
Drug = Drug[-1,];
rownames(Drug) = Drug[,1];
Drug = Drug[,-1];
class(Drug) <- "numeric"


idx = c();
for (i in 1:dim(Drug)[2]){
  ix = which(is.nan(Drug[,i]));
  ix1 = which(is.na(Drug[,i]));
  ix = union(ix,ix1);
  if(length(ix)>dim(Drug)[1]*0.1){
    idx =  rbind(idx,c(i)) ;
  }else{
    Drug[ix,i] = mean(Drug[-ix,i]);
  }
}
idx = as.vector(idx);
Drug = Drug[,-idx];


CellLine = intersect(rownames(Drug),rownames(Data));
idx1 = match(CellLine,rownames(Drug));
idx2 = match(CellLine,rownames(Data));
Drug = Drug[idx1,];
Data = t(Data[idx2,]);


Limit = as.matrix(read.delim(paste(dir_data,"/", LimitFile,sep=""),sep=",",header=F));
colnames(Limit) = Limit[1,];
Limit = Limit[-1,];
rownames(Limit) = Limit[,1];
Limit = Limit[,-1];
class(Limit) <- "numeric"

DrugNames = intersect(rownames(Limit),colnames(Drug));
Drug = Drug[,match(DrugNames,colnames(Drug))];
Limit = Limit[match(DrugNames,rownames(Limit)),];

DrugResponse = matrix("good",dim(Drug)[1],dim(Drug)[2]);
for (i in 1:dim(Drug)[2]){
  ix = which(Drug[,i]<mean(Drug[,i]));
  DrugResponse[ix,i] = "bad";
}
########## Data preprocessing

V = rep(0,dim(Data)[1]); # variance
for (i in 1:dim(Data)[1]){
  V[i] = var(Data[i,]);
}
M = rowMeans(Data);
idx_v = which(V > var_cut);
idx_m = order(M,decreasing=T); # sort the gene expression by mean expression
idx_m = idx_m[1:round((length(idx_m)*(1-mean_cut)))]; # index of the genes larger than the mean expression cutoff
ix = intersect(idx_v,idx_m);
Data = Data[ix,];

rm(ix, M, V, idx,idx1,idx2, idx_m, idx_v);

########## Prediction
CC_test = matrix(0,dim(Drug)[2],rr);
MSE_test = matrix(0,dim(Drug)[2],rr);
ROC_test = matrix(0,dim(Drug)[2],rr);


Feature_Num = round(Feature_Num/2);

for (mm in 1:dim(Drug)[2]){
  print(mm)
  AUC = Drug[,mm];
  AUCResponse = DrugResponse[,mm];
  for (nn in 1:rr){
    set.seed(nn);
    
    # 1. Prepare sample index for training, validation, and test data
    ix_p = sample.int(dim(Drug)[1], size = dim(Drug)[1], replace = FALSE); # shuffling the index
    n = round(length(ix_p)*0.33); # 33% of the data for testing
    ix_test = ix_p[1:n]; # index of the test cell lines
    ix_train = ix_p[(n+1):length(ix_p)]; # index of the training cell lines
    
    #################
    for (i in 1:1000){
      t=0;
      nl = length(which(DrugResponse[ix_test,mm]=="good"));
      if(nl==length(ix_test)){t = t+1;}
      if(nl==0){t = t+1;}
      nl = length(which(DrugResponse[ix_train,mm]=="good"));
      if(nl==length(ix_train)){t = t+1;}
      if(nl==0){t = t+1;}
      if(t==0){break;}else{
        ix_p = sample.int(dim(Drug)[1], size = dim(Drug)[1], replace = FALSE); # shuffling the index
        n = round(length(ix_p)*0.33); # 33% of the data for testing
        ix_test = ix_p[1:n]; # index of the test cell lines
        ix_train = ix_p[(n+1):length(ix_p)]; # index of the training cell lines
      }
    }
    
    #################
    
    # 2. Feature selection
    CC = cor(t(Data[,ix_train]),AUC[ix_train],method = "pearson");
    CC[is.na(CC)] = 0;
    idx = order(CC,decreasing=T);
    idx = c(idx[1:Feature_Num],tail(idx,Feature_Num));
    
    TrainData = Data[idx,ix_train]; # training data
    TestData = Data[idx,ix_test]; # test data
    rownames(TrainData) = c(1:dim(TrainData)[1]);
    rownames(TestData) = c(1:dim(TestData)[1]);
    
    # cross-validation (3 fold) for selecting the model
    plsr.predict = plsr(AUC[ix_train] ~. , data = data.frame(t(TrainData)), ncomp=2,validation = "CV", segment = 3, segment.type = "random")
    Y = predict(plsr.predict, t(as.matrix(TestData)),ncomp = 1:plsr.predict$ncomp); # predict on the test data
    Y_plsr = Y[,1,1];# choose the first component for prediction
    CC_test[mm,nn] = cor(Y_plsr,AUC[ix_test],method = "pearson"); #correlation between the predicted ones and the ground truth
    MSE_test[mm,nn] = mse(Y_plsr,AUC[ix_test])
    mod2 = try(roc(AUCResponse[ix_test],Y_plsr,levels=c("good","bad"),direction = ">"),TRUE);
    if(isTRUE(class(mod2)=="try-error")) { 
      ROC_test[mm,nn] = 0.5
    }else{r = mod2;
      ROC_test[mm,nn] = r$auc
    }
  }
}

save(CC_test,MSE_test,ROC_test,file="Expression_PLSR_AUC_CV.RData");

