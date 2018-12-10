## WOrking on Drug_ED50

# Intersection of Column From ED_50 and RNASeq_Log
i<-intersect(Drug_ED50$Drug,colnames(RNASeq_log[,-1]))


ED50_df1<-Drug_ED50[which(Drug_ED50$Drug %in% i==TRUE),]

# Cut-off percentage of columns in Drug_ED50
ED50_cut_off<- .2

 
#Calculate the percentage of NaN in All Columns

apply(ED50_df1,2, 
      function(col)sum(is.na(col))/length(col))>ED50_cut_off

# Selecting columns which contain less than NaN as per Cut_Off percentage.

df_Drug_ED50<-ED50_df1[which((apply(ED50_df1, 2,
function(col)sum(is.na(col))/length(col))>ED50_cut_off)==FALSE)]
#dim(df_Drug_ED50)
# 98 278

-------------------------------------------------------------------------------
  
  # Working on RNASeq_log CSV File.
  
  #View(RNASeq_log)
  
  
fil_RNASEq_Log<-cbind(RNASeq_log[,1]
,RNASeq_log[which(colnames(RNASeq_log[,]) %in% i==TRUE)])

# dim(fil_RNASEq_Log)

# 57071    99
  
df1_RNASEq_Log<-fil_RNASEq_Log
#Calculation MEan of All row
rowMeans(df1_RNASEq_Log[,-1],na.rm = TRUE)
#Combine the Data
df1_RNASEq_Log_m<-cbind(df1_RNASEq_Log,rowMeans(df1_RNASEq_Log[,-1]))
#Change Column name
names(df1_RNASEq_Log_m)[length(names(df1_RNASEq_Log_m))]<- 'Mean'





df2_RNASEq_Log<-fil_RNASEq_Log

#Calculating Variance of all Row

apply(df2_RNASEq_Log[,-1],1,var,na.rm =TRUE)

#Combining the data

df2_RNASEq_Log_v<-cbind(df2_RNASEq_Log,apply(df2_RNASEq_Log[,-1],1,var))

#Changing Column Name

names(df2_RNASEq_Log_v)[length(names(df2_RNASEq_Log_v))]<- 'Variance'

# View(df1_RNASEq_Log_m)
# View(df2_RNASEq_Log_v)
# dim(df2_RNASEq_Log_v)
# dim(df1_RNASEq_Log_m)




cut_off_RNASeq_log<-round(.5*nrow(RNASeq_log))


# Sorting of df1_RNASEq_Log_m and df2_RNASEq_Log_v 

s_df1_RNASEq_Log_m<-df1_RNASEq_Log_m[rev(order
(df1_RNASEq_Log_m$Mean,na.last = FALSE)),]


s_df2_RNASEq_Log_v<-df2_RNASEq_Log_v[rev(order
(df2_RNASEq_Log_v$Variance,na.last = FALSE)),]

# View(s_df1_RNASEq_Log_m)
# View(s_df2_RNASEq_Log_v)
# dim(s_df1_RNASEq_Log_m)
# dim(s_df2_RNASEq_Log_v)



#Last 40% of Data in s_df1_RNASEq_Log_m dataset
t<-tail(s_df1_RNASEq_Log_m$RNASeqLog,cut_off_RNASeq_log)

#Last 40% of Data in s_df2_v dataset
t1<-tail(s_df2_RNASEq_Log_v$RNASeqLog,cut_off_RNASeq_log)



# Union of Data set
u<-union(tail(s_df1_RNASEq_Log_m$RNASeqLog,cut_off_RNASeq_log)
         ,tail(s_df2_RNASEq_Log_v$RNASeqLog,cut_off_RNASeq_log))


#Eleminating union from fil_RNASEq_Log


df_RNASeq_log<-fil_RNASEq_Log[!(fil_RNASEq_Log$RNASeqLog %in% u ),]

# dim(df_RNASeq_log)
# 27215    99
# View(df_RNASeq_log)
# View(df_Drug_ED50)

# For FInding 0 Standard Deviation if there is any remove those features
#apply(df_RNASeq_log[,-1],1,sd,na.rm =TRUE)

#xtest<-cbind(df_RNASeq_log,apply(df_RNASeq_log[,-1],1,sd))





# Transposing of df_RNASeq_log

DF_RNA <- df_RNASeq_log 
n<-DF_RNA$RNASeqLog
DF_RNA<-as.data.frame(t(DF_RNA[,-1]))
colnames(DF_RNA)<-n

######################## Cross Validation Task ###########

#CROSS-VALIDATION.
# Split data set into 80% of Trainning,20% of Testing.
#ncol(df_Drug_ED50)
Radial<-''
drug<-''
Mean.drug<-''
for(xx in 2:ncol(df_Drug_ED50)){ 
for (t in 1:3){
    index     <- 1:nrow(DF_RNA)
    testindex <- sample(index, trunc(length(index)/5))
    
    train_DFRNA  <- DF_RNA[-testindex,]
    test_DFRNA   <- DF_RNA[testindex,]
    
    # dim(train_DFRNA)
    # dim(test_DFRNA)
    # dim(DF_RNA)
    
    # Ordering the column of sample_Drug_ED50 with respect--- 
    # --to Train_DFRNA too find correlation coffecient
    
    df_Drug_ED501<-df_Drug_ED50
    
    df_Drug_ED501$Drug<-ordered(df_Drug_ED501$Drug,levels = rownames(train_DFRNA[]))
    xdf_Drug_ED501<-df_Drug_ED501[order(df_Drug_ED501$Drug),]
    
    # Removing NaN from Drugs.
    xdf_Drug_ED501_d<-xdf_Drug_ED501[!is.na(xdf_Drug_ED501$Drug),1]
    xdf_Drug_ED501_dv<-xdf_Drug_ED501[!is.na(xdf_Drug_ED501$Drug),xx]
    bindtrain<-cbind.data.frame(xdf_Drug_ED501_d,xdf_Drug_ED501_dv)
    #cbind.data.frame(xdf_Drug_ED501_d,xdf_Drug_ED501_d1)
    
    #Check column-1 of View(xdf_Drug_ED501) is in similary sorted to.....   
    #...View(rownames(train_DFRNA))
    
    # Sorting and slicing the data according to test_DFRNA same as DF_Train
    df_Drug_ED502<-df_Drug_ED50
    
    df_Drug_ED502$Drug<-ordered(df_Drug_ED502$Drug,levels = rownames(test_DFRNA[]))
    xdf_Drug_ED502<-df_Drug_ED502[order(df_Drug_ED502$Drug),]
    
    xdf_Drug_ED502_d<-xdf_Drug_ED502[!is.na(xdf_Drug_ED502$Drug),1]
    
    xdf_Drug_ED502_dv<-xdf_Drug_ED502[!is.na(xdf_Drug_ED502$Drug),xx]
    bindtest<-cbind.data.frame(xdf_Drug_ED502_d,xdf_Drug_ED502_dv)
    #Check column-1 of View(xdf_Drug_ED502) is in similary sorted to.....   
    #...View(rownames(test_DFRNA))
    
    #   View(xdf_Drug_ED502)
    
    #which(is.na(xdf_Drug_ED501),)
    
    
    
    
    # Removing Coloumns where sd is zero
train_DFRNA <- train_DFRNA[,(!apply(train_DFRNA , 2 , function(x) sd(x)==0 ))==TRUE ]
    
    # Finding the Correlation of Drugs and Features
    df_corr<-tryCatch(as.data.frame(abs(cor(x 
    = train_DFRNA,y =bindtrain[,2],method = 'pearson',use="complete.obs" ))),silent = TRUE)
    
    # Making rowname as a 1st Column
    
    qpoint<-cbind.data.frame(rownames(df_corr),df_corr)
    
    # Changing the Column name 
    names(qpoint)[1]<- 'RNASeq_log'
    
    
    # Top 500 Correlation tail(sort(df_corr$SW001286),5)
    
    xpoint<-qpoint[qpoint$V1 %in% tail(sort(qpoint$V1),500),]
    
    
    
    # Merging training data for svm
    ytrain_DFRNA<-cbind.data.frame(bindtrain,train_DFRNA[,colnames(train_DFRNA) 
                                                         %in% xpoint$RNASeq_log])
    
    # Removing Nan From Training Data Sets
    ytrain_DFRNA<- ytrain_DFRNA[!(is.na(ytrain_DFRNA[2])),]
    
    # dim(ytrain_DFRNA)
    # View(ytrain_DFRNA)
    
    # Merging testing data for svm prediction
    
    ytest_DFRNA<-cbind.data.frame(bindtest,test_DFRNA
                                  [,colnames(test_DFRNA) %in% xpoint$RNASeq_log])
    
    # Removing Nan From Test Data Sets
    ytest_DFRNA<- ytest_DFRNA[!(is.na(ytest_DFRNA[2])),]
    
    # dim(ytest_DFRNA)
    # View(ytest_DFRNA)
    
    # Removing
    ################################################################################################################################
    # Initializing SVM Library..
    library(e1071)
    # Applying SVM Radial with tuning eps-Regression
svm.model.radial <- svm(ytrain_DFRNA[,-1][[1]]~ ., data = ytrain_DFRNA[,-1],type="eps-regression",cross=5,kernel='radial',
cost=c(2^-6,2^-4,2^-2,2^0,2^2,2^4,2^6,2^8),
gamma=c(2^-13,2^-11,2^-9,2^-7,2^-5,2^-3,2^-1,2^1))
    
    svm.pred.radial <-predict(svm.model.radial,ytest_DFRNA[,-1])
    
    pred.table<-cbind.data.frame(ytest_DFRNA[,1:2],svm.pred.radial)
    
    #View( pred.table)
  
    Radial[t]<-cor(y = pred.table[2],
                   x = pred.table[3],method = 'pearson')
  }
drug[xx]<-(colnames(df_Drug_ED501[xx]))
Mean.drug[xx]<-mean(as.numeric(Radial),na.rm = TRUE)


  }


  






