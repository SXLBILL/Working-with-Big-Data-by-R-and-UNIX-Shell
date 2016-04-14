# Appendix

###############################################
### Check if all pairs of files are matched ###
###############################################

unique(unlist(compareParallel(12,4)))

compare = function(fileindex, block = 1000, iter = 10000){
  # fileindex is the index of file, block is the number of data read 
  ## each time, iter is the number of blocks read
  # the function is to return a number 0 or 1, which represent no error  
  ## or error exists
  # the function is to check if trip_fare_.csv file and trip_data_.csv   
  ## file matched
  
  pattern1 = paste("cut -f 1,2 -d, trip_fare_", fileindex, ".csv",  sep="")
  con20 = pipe(pattern1, "r")
  pattern2 = paste("cut -f 1,2 -d, trip_data_", fileindex, ".csv", sep="")
  con21 = pipe(pattern2, "r")
  # initial the error
  error = 0
  for(i in 1:iter){
    # random sample an index from each block
    j = sample(block,1)
    a = readLines(con20, block)[j]
    b = readLines(con21, block)[j]
    if(a!=b){
      print(paste("file:", fileindex, ",", i*block+j))
      error = 1
      break
    }
  }
  return(error)
}


compareParallel = function(filenum, con){
  # filenum is the filesâ€? total number, con is the number of cluster 
  ## created
  # the function is to use parallel processing to check if all pairs of 
  ## files are matched  
  
  cl = makeCluster(con)
  clusterExport(cl, "compare")
  output = vector("list", filenum/con)
  for(i in 1:(filenum/con)){
    index1 = 1 + con*(i-1)
    index2 = con*i
    output[[i]] = parLapply(cl, c(index1:index2), compare)
  }
  stopCluster(cl)
  output
}


#############################
Version 1: R and shell
#############################

### compute the quantiles ###
# get a big frequency table of Y including all files
ytable_all = bigtableY(12)
cumfreq = cumsum(ytable_all$Freq)
size = sum(ytable_all$Freq)
quantiles = seq(0.1, 0.9, by = 0.1)
# compute the quantiles
sapply(quantiles, function(i) getQuantile(i, ytable_all, size, cumfreq))

### fit a regression and estimate beta0 and beta1 ###
# save all triptime data
sapply(1:12, function(i) getTime(i))
# estimate beta0 and beta1
parameters = predictbeta(12)

##################### Version 1 functions ##################
### Compute the deciles (version 1: shell and R version) ###

getDiff = function(fileindex){
  # fileindex is a trip_fare_.csv fileâ€™s index
  # the function returns a list which contains a vector of difference of   
  ## total_amount less total_tolls and a frequency table of the 
  ### difference
  # the function is to save the frequency table and the vector of
  ## total_amount less total_tolls
  
  # use shell command to read total_amount and total_tolls
  pattern11 = paste("cut -f 10,11 -d, trip_fare_", fileindex, ".csv", sep="")
  y = system(pattern11, intern = TRUE)[-1]
  # split characters and transform character to numeric
  y = strsplit(y, ",")
  y = lapply(y, as.numeric)
  diff = unlist(lapply(y, function(y) y[2]-y[1]))
  # get the frequency table
  fretable = data.frame(table(diff))
  Ys = list(table = fretable, diff = diff)
  save(Ys, file = paste("difflist_", fileindex, ".rda", sep=""))
  Ys
}


bigtableY = function(j){
  # j is an integer which represents the number of files
  # the function returns a big frequency table of Y including all files
  # initial a frequency table
  tabley = data.frame()
  for(i in 1:j){
    tabley = rbind(tabley, getDiff(i)$table)
  }
  # sort the Y
  tabley = tabley[order(tabley$diff),]
  # transform Y to integer
  tabley$diff = as.character(tabley$diff)
  tabley$diff = as.numeric(tabley$diff)
  tabley
}


getQuantile = function(singleQ, dftable, size, cumfreq){
  # singleQ is a percent, dftable is a frequency table, size is the   
  ## total number of the frequency table, cunfreq is cumulated sum of   
  ### the frequency
  # the function returns a quantile
  
  if(size*singleQ == ceiling(size*singleQ)){
    pos1 = size*singleQ
    pos2 = size*singleQ + 1
    Q1 = dftable[which(cumfreq>=pos1)[1],]$diff
    Q2 = dftable[which(cumfreq>=pos2)[1],]$diff
    Q = (Q1+Q2)/2
  }else{
    pos = ceiling(size*singleQ)
    Q = dftable[which(cumfreq>=pos)[1],]$diff
  }
  return(Q)
}


### fit the linear regression ###
getTime = function(fileindex){
  # fileindex is a trip_data_.csv fileâ€™s index
  # the function returns a list which contains a vector of triptime and  
  ## a frequency table of the triptime
  # the function is to save the frequency table and vector of triptime
  
  # use shell command to read triptime
  pattern9 = paste("cut -f 9 -d, trip_data_", fileindex, ".csv", sep="")
  triptime = system(pattern9, intern = TRUE)[-1]
  triptime = as.numeric(triptime)
  # save the column as vector and frequency table
  fretable = data.frame(table(triptime))
  Xs = list(table = fretable, time = triptime)
  save(Xs, file = paste("time_", fileindex, ".rda", sep=""))
}


predictbeta = function(filenum){
  # filenum is the filesâ€? total number to estimate the beta
  # the function is to estimate the beta0 and beta1 of simple regression
  
  # initial the parameter
  df = data.frame()
  sumx = 0
  sumy = 0
  sumxy = 0
  sumx2 = 0
  n= 0
  # calculate parameter of each file and sum them
  for(i in 1:filenum){
    pattern1 = paste("diff_", i, ".rda", sep="")
    load(pattern1)
    pattern2 = paste("time_", i, ".rda", sep="")
    load(pattern2)
    diff = Ys$diff
    triptime = Xs$time
    xy = diff*triptime
    sumxy = sumxy + sum(xy)
    sumx = sumx + sum(triptime)
    sumy = sumy + sum(diff)
    sumx2 = sumx2 + sum(triptime^2)
    n = n + length(diff)
  }
  meanx = sumx/n
  meany = sumy/n
  Sxx = sumx2 - n*(meanx^2)
  beta1.hat = (sumxy - n*meanx*meany)/Sxx
  beta0.hat = meany - beta1.hat * meanx
  list(beta0.hat = beta0.hat, beta1.hat = beta1.hat, 
                                Sxx = Sxx, meanx = meanx, n = n)
}


###########################################
Version 2: R and shell with Parallel
###########################################

# compute the quantiles
# get a big frequency table of Y including all files by parallel
ytable_all = getDifftableParallel(filenum = 12, con = 4)
cumfreq = cumsum(ytable_all$Freq)
size = sum(ytable_all$Freq)
quantiles = seq(0.1, 0.9, by = 0.1)
sapply(quantiles, function(i) getQuantile(i, ytable_all, size, cumfreq))

### fit a regression and estimate beta0 and beta1 ###
# save all triptime data by parallel processing
getTimeParallel(12, 4)
# estimate beta0 and beta1
parameters = predictbeta(12)


##################### Version 2 functions ##########################
# Compute the deciles (version 2: shell and R version using parallel) 
library(snow)
getDiff_table = function(i) getDiff(i)$table

getDifftableParallel = function(filenum, con){
  # filenum is the filesâ€? total number, con is the number of cluster 
  ## created
  # the function is to use parallel computation to get all response 
  ## variable and change them to a big frequency table
  
  # initial the first table
  tabley = data.frame()
  cl = makeCluster(con)
  clusterExport(cl, "getDiff_table")
  for(i in 1:(filenum/con)){
    index1 = 1 + con*(i-1)
    index2 = con*i
    result = parLapply(cl, c(index1:index2), getDiff_table)
    # combine the frequency table of response variable
    tabley = rbind(tabley, do.call(rbind, result))
  }
  stopCluster(cl)
  tabley
}

### fit the linear regression ###
getTimeParallel = function(filenum, con){
  # filenum is the filesâ€? total number, con is the number of cluster 
  ## created
  # the function is to use parallel computation to get all triptime data
  
  cl = makeCluster(con)
  clusterExport(cl, "getTime")
  for(i in 1:(filenum/con)){
    index1 = 1 + con*(i-1)
    index2 = con*i
    parLapply(cl, c(index1:index2), getTime)
  }
  stopCluster(cl)
}

# Other functions are the same as version 1


###########################################
### estimate the standard errors by BLB ###
###########################################

# use BLB to estimate the standard error of estimates
parameters = predictbeta(12)
n = parameters$n
meanx = parameters$meanx
Sxx = parameters$Sxx

Ysigma2.hat = BLBsigma(ytable_all, bag = 10, times = 10)
var_beta0.hat = sqrt(Ysigma2.hat*(1/n + meanx^2/Sxx))
var_beta1.hat = sqrt(Ysigma2.hat/Sxx)

BLBsigma = function(population, bag, times){
  # population is represented by a table, bag is the number of subsets 
  ## of size m =   to sample, times is the times to resample N 
  ### observations from each subsets
  # the function is to estimate the sigma square by BLB
  
  n0 = sum(population$Freq)
  vsigma = integer(bag)
  for(j in 1:bag){
    # From Y1,...,Yn to sample Y1,...,Yb
    ymtable = table(sample(population$diff, n0^(0.6), replace = TRUE, 												prob = population$Freq/n0))
    ymtable = df_table(ymtable)
    n1 = sum(ymtable$Freq)
    sigma2 = 0
    # From Y1,...,Yb to resample Y1*,â€?,Yn* and estimate sigma square
    for(i in 1:times) {
      sigma2 = sigma2 + getSigma2(ymtable, limit = 10, n0)
    }
    vsigma[j] = sigma2/times
  }
  mean(vsigma)
}


getSigma2 = function(ymtable, limit, n0){
  # ymtable is a frequency table of y. limit is to control the sampling 
  ## size each time, n0 is total number of population N
  # the function is to resample Y1*,...,Yn* from Y1,...,Yb and estimate  
  ## the sigma square each time
  
  n1 = sum(ymtable$Freq)
  temp = data.frame()
  for(i in 1:limit){
    yntable = table(sample(ymtable$diff, n0/limit, replace = TRUE, 
                           prob = ymtable$Freq/n1))
    yntable = df_table(yntable)
    yntable = rbind(temp, yntable)
    temp = yntable
  }
  
  sumy2 = sum(yntable$diff^2*yntable$Freq)
  meany = sum(yntable$diff*yntable$Freq)/n0
  sigma2 = (sumy2 - n0*meany^2)/(n0-1)
  sigma2
}


df_table = function(table){
  # the function is to transform a frequency table to a dataframe
  table = data.frame(table)
  names(table) = c("diff", "Freq")
  table$diff = as.character(table$diff)
  table$diff = as.numeric(table$diff)
  table
}


#######################################################
### Task 3: estimate multiple regression estimators ###
#######################################################

filenum = 12
con = 4
data = getDataParallel(filenum, con)
df =data.frame()
for(i in 1:(filenum/con)){
  df = rbind(df, do.call(rbind, data[[i]]))
}

model2 = lm(diff~surcharge + triptime, data = df)
summary(model2)


getData = function(fileindex, block = 1000, iter = 10000){
  # fileindex is the index of file, block is the number of data read 
  ## each time, iter is the number of blocks read
  # the function is to randomly sample observation and return its
  ## x1,x2,y
  
  pattern1 = paste("cut -f 7,10,11 -d, trip_fare_", fileindex,
                   ".csv", sep="")
  con20 = pipe(pattern1, "r")
  pattern2 = paste("cut -f 9 -d, trip_data_", fileindex,".csv", sep="")
  con21 = pipe(pattern2, "r")
  v = integer(iter)
  for(i in 1:iter){
    # random sample a index from each block
    j = sample(block,1)
    a = readLines(con20, block)[j]
    b = readLines(con21, block)[j]
    v[i] = paste(a, ",", b, sep ="")
  }
  # transform character to numeric
  y = strsplit(v,",")
  y = lapply(y, as.numeric)
  # output a dataframe with y,x1,x2
  diff = unlist(lapply(y, function(y) y[3]-y[2]))
  surcharge = unlist(lapply(y, function(y) y[1]))
  triptime = unlist(lapply(y, function(y) y[4]))
  df =data.frame(diff,surcharge,triptime)
}

getDataParallel = function(filenum, con){
  # filenum is the filesâ€? total number, con is the number of cluster 
  ## created
  # the function is to use parallel processing to do function "getData"
  
  cl = makeCluster(con)
  clusterExport(cl, "getData")
  output = vector("list", filenum/con)
  for(i in 1:(filenum/con)){
    index1 = 1 + con*(i-1)
    index2 = con*i
    output[[i]] = parLapply(cl, c(index1:index2), getData)
  }
  stopCluster(cl)
  output
}
