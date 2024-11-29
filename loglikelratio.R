library(igraph)

stationary_ts <- function(ts_data,diff=TRUE,normal=TRUE){
  if (diff & normal){
    new_ts <- apply(ts_data,1,diff)
    new_ts <- t(new_ts)
    sd_ts <- apply(new_ts,1,sd)
    new_ts <- apply(new_ts,1,function(x){x/sd(x)})
    return(list(new_ts,sd_ts))
  }
  if (normal & !diff){
    sd_ts <- apply(ts_data,1,sd)
    new_ts <- apply(ts_data,1,function(x){x/sd(x)})
    return(list(new_ts,sd_ts))
  }
  if (!normal & diff){
    new_ts <- apply(ts_data,1,diff)
    new_ts <- t(new_ts)
    return(new_ts)
  }
}

####### ####### #######
# LOAD DATA 
####### ####### #######
# Non evolving network structure
####### ####### ####### ####### 

# define possible change points Quartertly 2018-2020 (for all possible change points, leave a year before and year after in extremes to allow enough data points but also for lag=8)
change_points <- seq(39,72,3)

# check for possible change points, which time series are all 0s (before or after) to remove corresponding edges
# remove those edges for all change point cases for comparable results
ts_rem_bef <- ts_rem_aft <- c()
for (c in change_points){ 
  data_bef <- data[,3:c]
  data_aft <- data[,(c+1):93]
  ts_rem_bef <- c(ts_rem_bef,which(apply(data_bef,1,sum)==0))
  ts_rem_aft <- c(ts_rem_aft,which(apply(data_aft,1,sum)==0))
}

edge_ind_rem <- unique(c(ts_rem_bef,ts_rem_aft))

data <- data[-edge_ind_rem,]
data_edges=data[,c(1,2)]
data=data[,-c(1,2)]
data <- data.matrix(data)

graph_2dig <- graph_from_edgelist(as.matrix(data_edges),directed = TRUE)
globalalpha <- TRUE
alphaOrder <- 8
betaOrder <- rep(1,8)

nomin_ll_mle <- nomin_ll_c <- c()
est_var_res_bef <- est_var_res_aft <- c()
est_sigma_Ho <- c()
whole <- FALSE
concat <- TRUE
for (c in seq(39-2,72-2,3)){
  print(c("iteration ",c))
  if (whole==FALSE){ # difference and normalise after split
    data_bef <- data[,1:c]
    data_aft <- data[,(c+1):91]
    
    dtn_bef <- stationary_ts(data_bef)
    data_norm_bef <- t(dtn_bef[[1]])
    data_sd_bef <- dtn_bef[[2]] 
    
    dtn_aft <- stationary_ts(data_aft)
    data_norm_aft <- t(dtn_aft[[1]])
    data_sd_aft <- dtn_aft[[2]] 
    
    if(concat==TRUE){
      gnar_alpha8_stage1 <- gnar_edge_fit(cbind(data_norm_bef,data_norm_aft),data_edges,alphaOrder = alphaOrder, betaOrder = betaOrder,net = graph_2dig,lead_lag_mat = lead_lag_2sic,
                                          globalalpha = TRUE,lead_lag_weights = FALSE)
      #ll_all <- logLik(gnar_alpha8_stage1$mod)
      #ll_all <- sum(dnorm(x = gnar_alpha8_stage1$y, mean = predict(gnar_alpha8_stage1$mod), sd = summary(gnar_alpha8_stage1$mod)$sigma,log = TRUE))
      est_sigma_Ho <- c(est_sigma_Ho,(summary(gnar_alpha8_stage1[[1]])$sigma)^2)
    }
  }else{ # difference and normalise before split
    if (c==37){ 
      dtn <- stationary_ts(data)
      data_norm <- t(dtn[[1]])
      data_sd <- dtn[[2]] # keep sd of data_train
    }
    data_norm_bef <- data_norm[,1:(c-1)] # indices for time stamps change due to differencing (lost one value)
    data_norm_aft <- data_norm[,c:90]
  }
  # run GNAR-edge
  gnar_bef_alpha8_stage1 <- gnar_edge_fit(data_norm_bef,data_edges,alphaOrder = alphaOrder, betaOrder = betaOrder,net = graph_2dig,lead_lag_mat = lead_lag_2sic,
                                          globalalpha = TRUE,lead_lag_weights = FALSE)
  
  gnar_aft_alpha8_stage1 <- gnar_edge_fit(data_norm_aft,data_edges,alphaOrder = alphaOrder, betaOrder = betaOrder,net = graph_2dig,lead_lag_mat = lead_lag_2sic,
                                          globalalpha = TRUE,lead_lag_weights = FALSE)
  
  # log-likel for change point with logLik() function (for a model fitted by maximum likelihood)
  ll_bef <- logLik(gnar_bef_alpha8_stage1$mod)
  ll_aft <- logLik(gnar_aft_alpha8_stage1$mod)
  ll_mle <- ll_bef+ll_aft
  
  # log-likel for change point
  ll_bef_c <- sum(dnorm(x = gnar_bef_alpha8_stage1$y, mean = predict(gnar_bef_alpha8_stage1$mod), sd = summary(gnar_bef_alpha8_stage1$mod)$sigma,log = TRUE))
  ll_aft_c <- sum(dnorm(x = gnar_aft_alpha8_stage1$y, mean = predict(gnar_aft_alpha8_stage1$mod), sd = summary(gnar_aft_alpha8_stage1$mod)$sigma,log = TRUE))
  ll_c <- ll_bef_c+ll_aft_c
  
  nomin_ll_mle <- c(nomin_ll_mle,ll_mle)
  nomin_ll_c <- c(nomin_ll_c,ll_c)
  
  # estimated variance of residuals
  est_var_res_bef <- c(est_var_res_bef,(summary(gnar_bef_alpha8_stage1[[1]])$sigma)^2)
  est_var_res_aft <- c(est_var_res_aft,(summary(gnar_aft_alpha8_stage1[[1]])$sigma)^2)
}

# WHOLE=TRUE H0
dtn <- stationary_ts(data)
data_norm <- t(dtn[[1]])
data_sd <- dtn[[2]] 

gnar_alpha8_stage1 <- gnar_edge_fit(data_norm,data_edges,alphaOrder = alphaOrder, betaOrder = betaOrder,net = graph_2dig,lead_lag_mat = lead_lag_2sic,
                                    globalalpha = TRUE,lead_lag_weights = FALSE)
ll_all <- logLik(gnar_alpha8_stage1$mod)
ll_all <- sum(dnorm(x = gnar_alpha8_stage1$y, mean = predict(gnar_alpha8_stage1$mod), sd = summary(gnar_alpha8_stage1$mod)$sigma,log = TRUE))
est_sigma_Ho <- (summary(gnar_alpha8_stage1[[1]])$sigma)^2

# WHOLE=FALSE H0 it's obtained within for loop



####### ####### ####### ####### ####### ####### ####### ####### ####### 
####### RESULTS FROM SPLITTING TS AND THEN MAKING STATIONARY #######
####### ####### ####### ####### ####### ####### ####### ####### ####### 
# > ll_all
# 'log Lik.' -645043.9 (df=17)
ll_all <- -645043.9
# > nomin_ll_mle
nomin_ll_mle_2 <-  c(-577472.6 ,-576328.5, -574869.3 ,-573599.8 ,-575353.4, -575568.3 ,-572600.6, -569940.8)
# > nomin_ll_c
nomin_ll_c_2 <- c(-577472.6, -576328.5, -574869.3 ,-573599.8, -575353.4, -575568.3 ,-572600.6, -569940.8)

# > nomin_ll_c
nomin_ll_c <- c( -571863.5 ,-575547.3 ,-577852.4 ,-577388.4)
# > nomin_ll_mle
nomin_ll_mle <- c(-571863.5 ,-575547.3 ,-577852.4 ,-577388.4)


# Visaulise results
#df <- data.frame(l=nomin_ll_c-ll_all,t=colnames(data)[c(25,31,seq(43,81-2,12))])
df <- data.frame(l=c(nomin_ll_c_2,nomin_ll_c)-ll_all,t=colnames(data)[seq(39-2,72-2,3)])
lab_dates <- c("jan 18" ,"apr 18", "jul 18", "oct 18", "jan 19", "apr 19", "jul 19", "oct 19", "jan 20",
               "apr 20", "jul 20", "oct 20")
ggplot(data=df,aes(x=as.factor(seq(1,12)),y=l,group=1))+geom_line()+scale_x_discrete(name ="Timestamp", 
                                                                                     labels=lab_dates)+geom_point()+
  theme(axis.text.x = element_text(size=14, angle=90))+ylab(expression(lambda))+
  theme(text = element_text(size = 14))



####### ####### ####### ####### ####### ####### ####### ####### ####### 
####### RESULTS FROM SPLITTING TS AND THEN MAKING STATIONARY #######
####### ####### ####### ####### ####### ####### ####### ####### ####### 
#est_var_res_aft
est_var_res_aft2 <- c( 0.6695651, 0.6656061, 0.6605854, 0.6597552, 0.6710336, 0.6727952, 0.6560570, 0.6388466, 0.6576214, 0.6888073, 0.6838660, 0.6906907)

#est_var_res_bef
est_var_res_bef2 <- c( 0.6506199, 0.6503862, 0.6488356, 0.6439039, 0.6446036, 0.6461918, 0.6454378, 0.6438015, 0.6434984, 0.6480242, 0.6593508, 0.6577287)

#(summary(gnar_alpha8_stage1[[1]])$sigma)^2
est_sigma_Ho2 <-  0.6587443


#############################################################################
########### Look into parameters for oct19 possible change point ###########
#############################################################################

### SAMPLE 1: before τ time stamp possible change point

data_bef <- data[,1:(60-2)] # time stamp 60 is october 2019
data_aft <- data[,((60-2)+1):91]

# difference and normalise
dtn_bef <- stationary_ts(data_bef)
data_norm_bef <- t(dtn_bef[[1]])
data_sd_bef <- dtn_bef[[2]] 

dtn_aft <- stationary_ts(data_aft)
data_norm_aft <- t(dtn_aft[[1]])
data_sd_aft <- dtn_aft[[2]] 

# run GNAR-edge
gnar_bef_alpha8_stage1 <- gnar_edge_fit(data_norm_bef,data_edges,alphaOrder = alphaOrder, betaOrder = betaOrder,net = graph_2dig,lead_lag_mat = lead_lag_2sic,
                                        globalalpha = TRUE,lead_lag_weights = FALSE)

gnar_aft_alpha8_stage1 <- gnar_edge_fit(data_norm_aft,data_edges,alphaOrder = alphaOrder, betaOrder = betaOrder,net = graph_2dig,lead_lag_mat = lead_lag_2sic,
                                        globalalpha = TRUE,lead_lag_weights = FALSE)


# plot overlapping intervals
change_df <- data.frame(param=rep(names(gnar_aft_alpha8_stage1$mod$coefficients),2),lower=c(confint(gnar_bef_alpha8_stage1$mod)[,1],confint(gnar_aft_alpha8_stage1$mod)[,1]),
                        upper=c(confint(gnar_bef_alpha8_stage1$mod)[,2],confint(gnar_aft_alpha8_stage1$mod)[,2]),time=c(rep("before",16),rep("after",16)),
                        estimated=c(gnar_bef_alpha8_stage1$mod$coefficients,gnar_aft_alpha8_stage1$mod$coefficients))

library(latex2exp)
lab_param <- c(TeX(r'($\alpha_1$)'),TeX(r'($\alpha_2$)'),TeX(r'($\alpha_3$)'),TeX(r'($\alpha_4$)'),TeX(r'($\alpha_5$)'),TeX(r'($\alpha_6$)'),TeX(r'($\alpha_7$)'),TeX(r'($\alpha_8$)'),
               TeX(r'($\beta_{11}$)'),TeX(r'($\beta_{21}$)'),TeX(r'($\beta_{31}$)'),TeX(r'($\beta_{41}$)'),TeX(r'($\beta_{51}$)'),TeX(r'($\beta_{61}$)'),TeX(r'($\beta_{71}$)'),TeX(r'($\beta_{81}$)'))
ggplot(change_df,aes(param,estimated))+geom_errorbar(aes(ymin=lower,ymax=upper,color=time),position = position_dodge(0.6), width = 0.7)+
  geom_point(aes(color = time), position = position_dodge(0.6)) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"),name=expression(tau)) +theme(axis.text.x = element_text(angle = 90))+
  scale_x_discrete(name ="parameters", labels=lab_param)+ theme(text = element_text(size = 16)) + theme(legend.position="top")


