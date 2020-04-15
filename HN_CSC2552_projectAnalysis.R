setwd("C:/Users/larun/Desktop/Data Science Core/Projects/Urology/Image-analysis/post-hoc-nn-eval/p-val-calcs/")

###############################################
#
#         PACKAGES
#
###############################################

library("lubridate")
library("pROC")
library("PRROC")
library("ggplot2")
library("reshape")
library("RColorBrewer")
library("viridis")

###############################################
#
#         FUNCTIONS
#
###############################################

load_data = function(analysis_name,cv=FALSE){
  train = read.csv(paste0(analysis_name,"_train.csv"),header=TRUE,as.is=TRUE)
  if(cv){
    val = read.csv(paste0(analysis_name,"_val.csv"),header=TRUE,as.is=TRUE)
  }
  test = read.csv(paste0(analysis_name,"_test.csv"),header=TRUE,as.is=TRUE)
  
  if(cv){
    data_triad = list("train" = train,"val" = val,"test" = test)
  } else{
    data_triad = list("train" = train,"test" = test)
  }
  
  data_triad = lapply(data_triad,function(x){fix_gender(x)})
  data_triad = lapply(data_triad,function(x){x$manu = get_manu(x[,"full_ID"]) ; return(x)})
  
  full_dat = Reduce(rbind,data_triad)
  if(cv){
    full_dat$set = c(rep("train",nrow(data_triad[["train"]])),
                     rep("val",nrow(data_triad[["val"]])),
                     rep("test",nrow(data_triad[["test"]])))
    full_dat$Data_Split = factor(full_dat$set,levels = c("train","val","test"),labels = c("Training","Validation","Test"))
  } else{
    full_dat$set = c(rep("train",nrow(data_triad[["train"]])),
                     rep("test",nrow(data_triad[["test"]])))
    full_dat$Data_Split = factor(full_dat$set,levels = c("train","test"),labels = c("Training","Test"))
  }
  
  full_dat$us_num = as.numeric(full_dat$us_num)
  
  out_list = list("data_triad" = data_triad,
                  "full_dat" = full_dat)
  
  return(out_list)
}

make_hydro_df = function(full_df,phn_df,new_cols = c("SFU_grade","APD","ERP","Ureter.Dilation",
                                                     "Etiology","UTI","renal_scan1","renal_scan2",
                                                     "renal_scan3","surgery_type","surgery_type2",
                                                     "anomalies")){
  phn.raw.nodup = phn_df[!duplicated(phn_df$Study.ID),]
  phn.raw.nodup = phn.raw.nodup[!is.na(phn.raw$Study.ID),]
  hydro_only_full_dat = data.frame(matrix(nrow=0,ncol=ncol(full_df)+length(new_cols)))
  names(hydro_only_full_dat) <- c(names(full_df),new_cols)
  
  
  ## Renal scan  = nuclear scan 
  ##   -- can get multiple nuclear scans
  
  row = 1
  k = 1
  for(row in 1:nrow(full_df)){
    if(paste0(full_df$study_id[row],":",full_df$kidney_side[row]) %in% paste0(phn.raw.nodup$Study.ID,":",phn.raw.nodup$kidney_side)){
      X = k
      Pred_val = full_df$Pred_val[row]
      Target = full_df$Target[row]
      age_at_baseline = full_df$age_at_baseline[row]
      date_of_current_us = full_df$date_of_current_us[row]
      date_of_us1 = full_df$date_of_us1[row]
      full_ID = full_df$full_ID[row]
      gender = full_df$gender[row]
      kidney_side = full_df$kidney_side[row]
      study_id = full_df$study_id[row]
      cat("study id:\n")
      cat(paste0(study_id,"\n"))
      us_num = full_df$us_num[row]
      cat("us number: \n")
      cat(paste0(us_num,"\n"))
      manu = full_df$manu[row]
      set = full_df$set[row]
      Data_Split = full_df$Data_Split[row]
      anomalies = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,"Anomalies"]
      if(us_num == 1){
        SFU_grade = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,"SFU.Classification"]
        APD = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,"APD"]
        ERP = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,"ERP.diamater"]
      } else if(us_num == 2){
        SFU_grade = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,"SFU.Grade"]
        APD = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,paste0("APD.",us_num-1)]
        ERP = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,paste0("ERP.diamater.",us_num-1)]
      } else{
        SFU_grade = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,paste0("SFU.Grade.",us_num-2)]
        APD = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,paste0("APD.",us_num-1)]
        ERP = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,paste0("ERP.diamater.",us_num-1)]
      }
      cat("SFU grade:\n")
      cat(paste0(SFU_grade,"\n"))
      Ureter.Dilation = ifelse(is.null(phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,paste0("Ureter.Dilation.",us_num)]),
                               yes = NA,
                               no = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,paste0("Ureter.Dilation.",us_num)])
      Etiology = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,"Etiology"]
      UTI = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,"UTI"]
      renal_scan1 = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,"Renal.Scan.1"]
      renal_scan2 = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,"Renal.Scan.2"]
      renal_scan3 = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,"Renal.Scan.3"]
      surgery_type = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,"Type"]
      surgery_type2 = phn.raw.nodup[phn.raw.nodup$Study.ID == study_id,"Type.1"]
      
      in_row = c(SFU_grade,APD,ERP,
                 Ureter.Dilation,Etiology,UTI,renal_scan1,renal_scan2,
                 renal_scan3,surgery_type,surgery_type2,anomalies)
      
      
      hydro_only_full_dat[k,] <- c(full_df[row,],in_row)
      k = k+1
    }
  }
  
  return(hydro_only_full_dat)  
}
get_manu = function(full_id_vec){
  manu = unlist(lapply(strsplit(full_id_vec,split = "_"),function(x){x[length(x)]}))
  return(manu)
}

fix_gender = function(df){
  df_out = df
  df_out$date_of_current_us[!(df_out$gender %in% c("F","M"))] <- df_out$date_of_us1[!(df_out$gender %in% c("F","M"))]
  df_out$date_of_us1[!(df_out$gender %in% c("F","M"))] <- df_out$kidney_side[!(df_out$gender %in% c("F","M"))]
  df_out$kidney_side[!(df_out$gender %in% c("F","M"))] <- df_out$us_num[!(df_out$gender %in% c("F","M"))]
  df_out$us_num[!(df_out$gender %in% c("F","M"))] <- df_out$gender[!(df_out$gender %in% c("F","M"))]
  df_out$gender[!(df_out$gender %in% c("F","M"))] <- "F"
  
  return(df_out)
}

elapsed_months <- function(end_date, start_date) { ## from: https://stackoverflow.com/questions/1995933/number-of-months-between-two-dates/1996404
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}

get_cutpoint <- function(pos_class_vec,sensitivity = 0.95){
  sorted_vec = sort(pos_class_vec)
  cutpoint_thresh = quantile(sorted_vec,c(1-sensitivity))
  return(cutpoint_thresh)
}

get_pred_class <- function(pred_vals,threshold = 0.5){
  return(ifelse(pred_vals > threshold,yes = 1,no = 0))
}

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

revise_full_data = function(full_df,data_triad,cv,scale_model=NULL){
  data_triad_sub = lapply(data_triad,function(x){x[paste0(x[,"study_id"],":",x[,"kidney_side"]) %in% paste0(full_df$study_id,":",full_df$kidney_side),]})
  if(cv){
    full_df$set = c(rep("train",nrow(data_triad_sub[["train"]])),
                    rep("val",nrow(data_triad_sub[["val"]])),
                    rep("test",nrow(data_triad_sub[["test"]])))
    full_df$Data_Split = factor(full_df$set,
                                levels = c("train","val","test"),
                                labels = c("Training","Validation","Test"))
  } else{
    full_df$set = c(rep("train",nrow(data_triad_sub[["train"]])),
                    rep("test",nrow(data_triad_sub[["test"]])))
    full_df$Data_Split = factor(full_df$set,
                                levels = c("train","test"),
                                labels = c("Training","Test"))
  }
  full_df$Target.f = factor(full_df$Target,
                            levels = c(0,1),
                            labels = c("No Surgery","Surgery"))
  # full_df$date_of_us1.date = as.Date(full_df$date_of_us1)
  # full_df$date_of_current_us.date = as.Date(full_df$date_of_current_us)
  # full_df$date_of_current_us.date[full_df$date_of_current_us.date > "2020-01-01"] = NA
  # full_df$us_1_yr = year(full_df$date_of_us1.date) ##  can extract year from lubridate function
  # full_df$us_yr = year(full_df$date_of_current_us.date) ##  can extract year from lubridate function
  # full_df$us_yr[full_df$us_yr > 2020] = NA
  # 
  # full_df$age_at_us = full_df$age_at_baseline + elapsed_months(end_date = full_df$date_of_current_us.date,start_date = full_df$date_of_us1.date)
  # full_df$us_num.f = factor(full_df$us_num,levels = 1:10)
  
  if(!is.null(scale_model)){
    full_df$Scaled_Pred = predict(object = scale_model,newdata=full_df,type="response")
  }
  
  return(full_df)  
}

get_platt_model = function(pth_file_root,is_cv = FALSE){
  my_analysis_name = pth_file_root
  out_full = load_data(analysis_name = my_analysis_name,cv = is_cv)
  
  full_data_triad = out_full$data_triad
  full_dat = out_full$full_dat
  full_cv = is_cv
  
  full_1 = revise_full_data(full_df = full_dat,
                            data_triad = full_data_triad,
                            cv = full_cv)
  
  platt_model = glm(Target ~ Pred_val,
                    data = full_1[full_1$Data_Split == "Training",],
                    family = binomial(link=logit))
  
  scaled_df = revise_full_data(full_df = full_dat,
                               data_triad = full_data_triad,
                               cv = full_cv,
                               scale_model = platt_model)
  
  out_list = list(scaled_df,platt_model)
  names(out_list) = c("scaled_df","platt_model")
  return(out_list)  
}

get_scaled_df = function(pth_file_root,is_cv = FALSE){
  my_analysis_name = pth_file_root
  out_full = load_data(analysis_name = my_analysis_name,cv = is_cv)
  
  full_data_triad = out_full$data_triad
  full_dat = out_full$full_dat
  full_cv = is_cv
  
  full_1 = revise_full_data(full_df = full_dat,
                            data_triad = full_data_triad,
                            cv = full_cv)
  
  platt_model = glm(Target ~ Pred_val,
                    data = full_1[full_1$Data_Split == "Training",],
                    family = binomial(link=logit))
  
  scaled_df = revise_full_data(full_df = full_dat,
                               data_triad = full_data_triad,
                               cv = full_cv,
                               scale_model = platt_model)
  
  return(scaled_df)  
}

get_auc = function(df,sfu_grade = NULL,data_split = "Test"){
  df = df[!is.na(df$Target),]
  
  if(is.na(var(df$Target[df$Data_Split == data_split]))){
    my_var = 0
  } else{
    my_var = var(df$Target[df$Data_Split == data_split])
  }
  
  if(my_var > 0){
    if(is.null(sfu_grade)){
      my_auc = auc(df$Target[df$Data_Split == data_split],
                   df$Scaled_Pred[df$Data_Split == data_split])  
    } else{
      my_auc = auc(df$Target[df$SFU_grade == sfu_grade & df$Data_Split == data_split],
                   df$Scaled_Pred[df$SFU_grade == sfu_grade & df$Data_Split == data_split])  
    }
  } else{
    my_auc = NA
  }
  
  return(my_auc)
}

get_f1 = function(df,threshold = 0.05,data_split = "Training",scaled = TRUE){
  df = df[!is.na(df$Target),]
  if(is.na(var(df$Target[df$Data_Split == data_split]))){
    my_var = 0
  } else{
    my_var = var(df$Target[df$Data_Split == data_split])
  }
  
  if(my_var > 0){
    if(scaled){
      no.surg.scores = na.omit(df$Scaled_Pred[df$Target == 0 & df$Data_Split == data_split])
      surg.scores = na.omit(df$Scaled_Pred[df$Target == 1 & df$Data_Split == data_split])
    } else{
      no.surg.scores = na.omit(df$Pred_val[df$Target == 0 & df$Data_Split == data_split])
      surg.scores = na.omit(df$Pred_val[df$Target == 1 & df$Data_Split == data_split])
    } 
    
    TP = sum(surg.scores >= 0.05)
    FN = sum(surg.scores < 0.05)
    TN = sum(no.surg.scores < 0.05)
    FP = sum(no.surg.scores >= 0.05)
    
    FNR = sum(surg.scores < 0.05)/length(surg.scores)
    FPR = sum(no.surg.scores >= 0.05)/length(no.surg.scores)
    
    precision = TP/(TP+FP)
    recall = TP/(TP+FN) ## same as sensitivity
    specificity = TN/(TN+FP)
    
    f1 = 2*(precision*recall)/(precision+recall)
    
    out_list = list(TP,FNR,TN,FPR,precision,recall,specificity,f1)
    names(out_list) = c("TP","FNR","TN","FPR","precision","recall","specificity","f1")
    return(out_list)
  } else{
    out_list = list(NA,NA,NA,NA,NA,NA,NA,NA)
    names(out_list) = c("TP","FNR","TN","FPR","precision","recall","specificity","f1")
    return(out_list)
  }
    
}

get_auprc = function(df,data_split = "Training",scaled = TRUE){
  df = df[!is.na(df$Target),]
  if(is.na(var(df$Target[df$Data_Split == data_split]))){
    my_var = 0
  } else{
    my_var = var(df$Target[df$Data_Split == data_split])
  }
  
  if(my_var > 0){
    if(scaled){
      no.surg.scores = na.omit(df$Scaled_Pred[df$Target == 0 & df$Data_Split == data_split])
      surg.scores = na.omit(df$Scaled_Pred[df$Target == 1 & df$Data_Split == data_split])
    } else{
      no.surg.scores = na.omit(df$Pred_val[df$Target == 0 & df$Data_Split == data_split])
      surg.scores = na.omit(df$Pred_val[df$Target == 1 & df$Data_Split == data_split])
    }
    
    my_pr = pr.curve(scores.class0 = no.surg.scores,scores.class1 = surg.scores)
    return(my_pr$auc.integral)    
  } else{
    return(NA)
  }

}

get_roc = function(df,data_split = "Test"){
  if(var(df$Target[df$Data_Split == data_split]) > 0){
    df_sub = df[df$Data_Split == data_split,]  
    myroc = roc(df$Target[order(df$Scaled_Pred)],
                df$Scaled_Pred[order(df$Scaled_Pred)])
    return(myroc)    
  } else{
    return(NA)
  }

}

model_comparison = function(mod_a_pth,mod_b_pth){
  ## TRUE MOD COMPAIRSON
  ## AUROC
  true_auroc_a = get_auc(df = mod_a_pth,data_split = "Test")
  true_auroc_b = get_auc(df = mod_b_pth,data_split = "Test") 
  
  true_auroc_diff = true_auroc_a - true_auroc_b
  
  ## AUPRC
  true_auprc_a = get_auprc(df = mod_a_pth,data_split = "Test")
  true_auprc_b = get_auprc(df = mod_b_pth,data_split = "Test") 
  
  true_auprc_diff = true_auprc_a - true_auprc_b
  
  ## PERM MOD COMPARISON  
  test_mod_a = mod_a_pth[mod_a_pth$Data_Split == "Test",]
  test_mod_b = mod_b_pth[mod_b_pth$Data_Split == "Test",]
  perm_auroc_diffs = rep(NA,999)
  perm_auprc_diffs = rep(NA,999)
  for(seed_val in 1:999){
    set.seed(seed_val)
    idx = sample(x = 1:nrow(test_mod_a),size = nrow(test_mod_a),replace = TRUE)
    perm_auroc_a = get_auc(df = test_mod_a[idx,], data_split = "Test")
    perm_auroc_b = get_auc(df = test_mod_b[idx,], data_split = "Test")
    
    perm_auprc_a = get_auprc(df = test_mod_a[idx,],sfu_grade = NULL, data_split = "Test")
    perm_auprc_b = get_auprc(df = test_mod_b[idx,],sfu_grade = NULL, data_split = "Test")
    
    perm_auroc_diffs[seed_val] = perm_auroc_a - perm_auroc_b
    perm_auprc_diffs[seed_val] = perm_auprc_a - perm_auprc_b
  }
  
  par(mfrow=c(1,2))
  
  hist(perm_auroc_diffs,main = "AUROC")
  abline(v = 0,col = "black")
  abline(v = true_auroc_diff,col = "red")
  
  hist(perm_auprc_diffs,main = "AUPRC")
  abline(v = 0,col = "black")
  abline(v = true_auprc_diff,col = "red")
  
}

get_ci = function(mod_a_pth,split = "Test",thresh = 0.05,scaled_pred = TRUE){
  ## TRUE MOD COMPAIRSON
  ## AUROC
  true_auroc_a = get_auc(df = mod_a_pth,data_split = split)
  
  ## AUPRC
  true_auprc_a = get_auprc(df = mod_a_pth,data_split = split,scaled = scaled_pred)
  
  ## F1
  true_f1_score = get_f1(df = mod_a_pth,data_split = split,threshold = thresh)$f1

  ## FN
  true_fn_score = get_f1(df = mod_a_pth,data_split = split,threshold = thresh)$FNR
  
  ## FP
  true_fp_score = get_f1(df = mod_a_pth,data_split = split,threshold = thresh)$FPR
  
  
  ## PERM MOD COMPARISON  
  test_mod_a = mod_a_pth[mod_a_pth$Data_Split == split,]
  perm_aurocs = rep(NA,999)
  perm_auprcs = rep(NA,999)
  perm_f1s = rep(NA,999)
  perm_fns = rep(NA,999)
  perm_fps = rep(NA,999)
  for(seed_val in 1:999){
    set.seed(seed_val)
    idx = sample(x = 1:nrow(test_mod_a),size = nrow(test_mod_a),replace = TRUE)
    perm_aurocs[seed_val] = get_auc(df = test_mod_a[idx,], data_split = split)
    
    perm_auprcs[seed_val] = get_auprc(df = test_mod_a[idx,],data_split = split,scaled = scaled_pred)

    f1_out = get_f1(df = test_mod_a[idx,],data_split = split,scaled = scaled_pred,threshold = thresh)
        
    perm_f1s[seed_val] = f1_out$f1
    
    perm_fns[seed_val] = f1_out$FNR
    
    perm_fps[seed_val] = f1_out$FPR
  }
  
  par(mfrow=c(2,2))

  if(!all(is.na(perm_aurocs))){
    hist(perm_aurocs,main = "AUROC")
    abline(v = 0,col = "black")
    abline(v = true_auroc_a,col = "red")
  }

  if(!all(is.na(perm_auprcs))){
    hist(perm_auprcs,main = "AUPRC")
    abline(v = 0,col = "black")
    abline(v = true_auprc_a,col = "red")
  }
  
  # if(!all(is.na(perm_f1s))){
  #   hist(perm_f1s,main = "F1")
  #   abline(v = 0,col = "black")
  #   abline(v = true_f1_score,col = "red")
  # }
  
  if(!all(is.na(perm_fns))){
    hist(perm_fns,main = "FN")
    abline(v = 0,col = "black")
    abline(v = true_fn_score,col = "red")
  }

  if(!all(is.na(perm_fps))){
    hist(perm_fps,main = "FP")
    abline(v = 0,col = "black")
    abline(v = true_fp_score,col = "red")
  }
  
  par(mfrow=c(1,1))
  
  sorted_auroc = sort(na.omit(perm_aurocs))
  auroc_lcl = quantile(sorted_auroc,0.025)
  auroc_ucl = quantile(sorted_auroc,0.975)
  
  sorted_auprc = sort(na.omit(perm_auprcs))
  auprc_lcl = quantile(sorted_auprc,0.025)
  auprc_ucl = quantile(sorted_auprc,0.975)

  sorted_f1 = sort(na.omit(perm_f1s))
  f1_lcl = quantile(sorted_f1,0.025)
  f1_ucl = quantile(sorted_f1,0.975)
  
  sorted_fn = sort(na.omit(perm_fns))
  fn_lcl = quantile(sorted_fn,0.025)
  fn_ucl = quantile(sorted_fn,0.975)

  sorted_fp = sort(na.omit(perm_fps))
  fp_lcl = quantile(sorted_fp,0.025)
  fp_ucl = quantile(sorted_fp,0.975)
  
  outlist = list("auroc" = true_auroc_a,
                 "auroc_lcl" = auroc_lcl,
                 "auroc_ucl" = auroc_ucl,
                 "auprc" = true_auprc_a,
                 "auprc_lcl" = auprc_lcl,
                 "auprc_ucl" = auprc_ucl,
                 "f1" = true_f1_score,
                 "f1_lcl" = f1_lcl,
                 "f1_ucl" = f1_ucl,
                 "fn" = true_fn_score,
                 "fn_lcl" = fn_lcl,
                 "fn_ucl" = fn_ucl,                 
                 "fp" = true_fp_score,
                 "fp_lcl" = fp_lcl,
                 "fp_ucl" = fp_ucl)
  return(outlist)
}

get_trans_sag_ens = function(trans_df,sag_df,cv=TRUE){
  if(cv){
    vars_to_extract = c("study_id","us_num","Fold","kidney_side","Pred_val","Target","Scaled_Pred","Data_Split")
  } else{
    vars_to_extract = c("study_id","us_num","kidney_side","Pred_val","Target","Scaled_Pred","Data_Split")
  }
  
  trans_sub = trans_df[,vars_to_extract]
  sag_sub = sag_df[,vars_to_extract]
  
  names(trans_sub) = paste0(names(trans_sub),"_trans")
  names(sag_sub) = paste0(names(sag_sub),"_sag")
  
  if(cv){
    trans_sub$study_us_num = paste0(trans_sub$study_id_trans,"_",
                                    trans_sub$us_num_trans,"_",
                                    trans_sub$kidney_side_trans,"_",
                                    trans_sub$Data_Split_trans,"_",
                                    trans_sub$Fold_trans)
    sag_sub$study_us_num = paste0(sag_sub$study_id_sag,"_",
                                  sag_sub$us_num_sag,"_",
                                  sag_sub$kidney_side_sag,"_",
                                  sag_sub$Data_Split_sag,"_",
                                  sag_sub$Fold_sag)
  } else{
    trans_sub$study_us_num = paste0(trans_sub$study_id_trans,"_",trans_sub$us_num_trans,"_",trans_sub$kidney_side_trans)
    sag_sub$study_us_num = paste0(sag_sub$study_id_sag,"_",sag_sub$us_num_sag,"_",sag_sub$kidney_side_sag)
  }
  
  merged_dfs = merge(trans_sub,sag_sub,by = "study_us_num",all = TRUE)
  merged_dfs$Target = merged_dfs$Target_trans
  merged_dfs$Target[is.na(merged_dfs$Target_trans)] <- merged_dfs$Target_sag[is.na(merged_dfs$Target_trans)]
  merged_dfs$Data_Split = merged_dfs$Data_Split_sag
  merged_dfs$Data_Split[is.na(merged_dfs$Data_Split_sag)] = merged_dfs$Data_Split_trans[is.na(merged_dfs$Data_Split_sag)]
  
  pred_mod = glm(Target ~ Pred_val_trans + Pred_val_sag, data = merged_dfs[merged_dfs$Data_Split == "Training",],family=binomial)
  merged_dfs$Scaled_Pred = predict(pred_mod,newdata = merged_dfs,type = "response")
  merged_dfs$Scaled_Pred[is.na(merged_dfs$Scaled_Pred) & is.na(merged_dfs$Scaled_Pred_trans)] = merged_dfs$Scaled_Pred_sag[is.na(merged_dfs$Scaled_Pred) & is.na(merged_dfs$Scaled_Pred_trans)]
  merged_dfs$Scaled_Pred[is.na(merged_dfs$Scaled_Pred) & is.na(merged_dfs$Scaled_Pred_sag)] = merged_dfs$Scaled_Pred_trans[is.na(merged_dfs$Scaled_Pred) & is.na(merged_dfs$Scaled_Pred_sag)]
  
  
  # merged_dfs$Scaled_Pred = sapply(1:nrow(merged_dfs),function(i){mean(na.omit(c(merged_dfs$Scaled_Pred_sag[i],
  #                                                      merged_dfs$Scaled_Pred_sag[i])))})
  # merged_dfs$Pred_val = merged_dfs$Pred_val_sag
  # merged_dfs$Data_Split = merged_dfs$Data_Split_trans
  
  return(merged_dfs)    
}

rename_mods = function(names_vec,suffix){
  return(c(names_vec[1],paste0(names_vec[2:length(names_vec)],suffix)))
}

get_big_ens_df = function(our_mod_ens_df,densenet_mod_ens_df,
                          resnet_mod_ens_df,vgg_mod_ens_df){
  names(our_mod_ens_df) = rename_mods(names(our_mod_ens_df),"_ourmod")
  names(densenet_mod_ens_df) = rename_mods(names(densenet_mod_ens_df),"_densenet")
  names(resnet_mod_ens_df) = rename_mods(names(resnet_mod_ens_df),"_resnet")
  names(vgg_mod_ens_df) = rename_mods(names(vgg_mod_ens_df),"_vgg")
  
  mod_list = list(our_mod_ens_df,densenet_mod_ens_df,
                  resnet_mod_ens_df,vgg_mod_ens_df)
  
  merged_ens_mods_df = Reduce(function(x,y){merge(x,y,by = "study_us_num")},
                              mod_list)
  
  merged_ens_mods_df$Target = merged_ens_mods_df$Target_densenet
  merged_ens_mods_df$Data_Split = merged_ens_mods_df$Data_Split_densenet
  
  ens_model = glm(Target ~ Scaled_Pred_ourmod + 
                    Scaled_Pred_resnet + Scaled_Pred_vgg + 
                    Scaled_Pred_densenet,
                  data = merged_ens_mods_df[merged_ens_mods_df$Data_Split == "Training",],
                  family = binomial)
  
  merged_ens_mods_df$Scaled_Pred = predict(ens_model,newdata = merged_ens_mods_df,type = "response")  
  
  return(merged_ens_mods_df)
}

get_big_st_ens_df = function(our_mod_ens_df,densenet_mod_ens_df,
                             resnet_mod_ens_df,vgg_mod_ens_df){
  names(our_mod_ens_df) = rename_mods(names(our_mod_ens_df),"_ourmod")
  names(densenet_mod_ens_df) = rename_mods(names(densenet_mod_ens_df),"_densenet")
  names(resnet_mod_ens_df) = rename_mods(names(resnet_mod_ens_df),"_resnet")
  names(vgg_mod_ens_df) = rename_mods(names(vgg_mod_ens_df),"_vgg")
  
  mod_list = list(our_mod_ens_df,densenet_mod_ens_df,
                  resnet_mod_ens_df,vgg_mod_ens_df)
  
  merged_ens_mods_df = Reduce(function(x,y){merge(x,y,by = "study_us_num")},
                              mod_list)
  
  merged_ens_mods_df$Target = merged_ens_mods_df$Target_densenet
  merged_ens_mods_df$Data_Split = merged_ens_mods_df$Data_Split_densenet
  
  ens_model = glm(Target ~ Scaled_Pred_sag_ourmod + Scaled_Pred_trans_ourmod + 
                    Scaled_Pred_sag_resnet + Scaled_Pred_trans_resnet + 
                    Scaled_Pred_sag_vgg + Scaled_Pred_trans_vgg + 
                    Scaled_Pred_sag_densenet + Scaled_Pred_trans_densenet,
                  data = merged_ens_mods_df[merged_ens_mods_df$Data_Split == "Training",],
                  family = binomial)
  
  merged_ens_mods_df$Scaled_Pred = predict(ens_model,newdata = merged_ens_mods_df,type = "response")  
  
  return(merged_ens_mods_df)
}

get_auc_df = function(score_df, data_split, my_thresh = 0.05, var_of_int = NULL){
  if(is.null(var_of_int)){
    cat("Please specify variable of interest")
  } else{
    uniq_var_vals = unique(score_df[,var_of_int])
    
    df_cols = c(var_of_int,"auroc","auroc_lcl","auroc_ucl",
                "auprc","auprc_lcl","auprc_ucl",
                "f1","f1_lcl","f1_ucl",
                "fn","fn_lcl","fn_ucl",
                "fp","fp_lcl","fp_ucl")
    auc_df = data.frame(matrix(nrow = 0,ncol = length(df_cols)))
    names(auc_df) <- df_cols
    
    row = 1
    for(my_val in uniq_var_vals){
      auc_df[row,] = rep(NA,ncol(auc_df))
      auc_df[row,var_of_int] = my_val
      
      if(is.na(var(na.omit(score_df[score_df[,var_of_int] == my_val,"Target"])))){
        my_var = 0
      } else{
        my_var = var(na.omit(score_df[score_df[,var_of_int] == my_val,"Target"]))
      }
      
      if(my_var > 0){
        cat(na.omit(score_df[score_df[,var_of_int] == my_val,"Target"]))
        aucs = get_ci(score_df[score_df[,var_of_int] == my_val,],split = data_split,thresh = my_thresh)
        auc_df[row,"auroc"] = aucs$auroc
        auc_df[row,"auroc_lcl"] = aucs$auroc_lcl
        auc_df[row,"auroc_ucl"] = aucs$auroc_ucl
        auc_df[row,"auprc"] = aucs$auprc
        auc_df[row,"auprc_lcl"] = aucs$auprc_lcl
        auc_df[row,"auprc_ucl"] = aucs$auprc_ucl
        auc_df[row,"f1"] = aucs$f1
        auc_df[row,"f1_lcl"] = aucs$f1_lcl
        auc_df[row,"f1_ucl"] = aucs$f1_ucl
        auc_df[row,"fn"] = aucs$fn
        auc_df[row,"fn_lcl"] = aucs$fn_lcl
        auc_df[row,"fn_ucl"] = aucs$fn_ucl
        auc_df[row,"fp"] = aucs$fp
        auc_df[row,"fp_lcl"] = aucs$fp_lcl
        auc_df[row,"fp_ucl"] = aucs$fp_ucl
      } else{
        auc_df[row,"auroc"] = NA
        auc_df[row,"auroc_lcl"] = NA
        auc_df[row,"auroc_ucl"] = NA
        auc_df[row,"auprc"] = NA
        auc_df[row,"auprc_lcl"] = NA
        auc_df[row,"auprc_ucl"] = NA
        auc_df[row,"f1"] = NA
        auc_df[row,"f1_lcl"] = NA
        auc_df[row,"f1_ucl"] = NA
        auc_df[row,"fn"] = NA
        auc_df[row,"fn_lcl"] = NA
        auc_df[row,"fn_ucl"] = NA
        auc_df[row,"fp"] = NA
        auc_df[row,"fp_lcl"] = NA
        auc_df[row,"fp_ucl"] = NA
      }
      row = row + 1
    }
  }
  return(auc_df)
}

###
###     Clinical patient data 
###

phn.raw = read.csv("C:/Users/larun/Desktop/Data Science Core/Projects/Urology/Image-analysis/post-hoc-nn-eval/20190524_combinedObstRef.csv", #sep = "\t",
                   header=TRUE,as.is=TRUE)

## relevant dataframe = phn.raw ; want vcug1 (or all VCUG variables really)
phn.raw$kidney_side = ifelse(test = (phn.raw$Laterality == "Bilateral"),
                             yes = phn.raw$If.bilateral..which.is.the.most.severe.kidney.,
                             no = phn.raw$Laterality)
str(phn.raw$kidney_side)

###
###     Scores from predictive model 
###

predict_mod_scaled = get_scaled_df(pth_file_root = "prehdict_20190802_vanilla_siamese_dim256_c1_full")
str(predict_mod_scaled)

  ## confirming Test set AUROC/AUPRC/FN/FP/etc
get_ci(predict_mod_scaled,split = "Training")
get_ci(predict_mod_scaled,split = "Test")

###
###     Provider and postal code info
###

ppc = read.csv("C:/Users/larun/Desktop/Data Science Core/Projects/Urology/Image organization Nov 2019/dist_provider_labels_20200414.csv",
               header=TRUE,as.is=TRUE)
head(ppc)

predict_mod_scaled$provider = NA
predict_mod_scaled$provider = ppc$Who.indicated.Surgery[match(predict_mod_scaled$study_id,ppc$REDCap.Number)]
table(predict_mod_scaled$provider[predict_mod_scaled$Data_Split == "Training"])
table(predict_mod_scaled$provider[predict_mod_scaled$Data_Split == "Test"])

predict_mod_scaled$pc = NA
predict_mod_scaled$pc = ppc$Postal.Code[match(predict_mod_scaled$study_id,ppc$REDCap.Number)]
table(predict_mod_scaled$pc[predict_mod_scaled$Data_Split == "Training"])
table(predict_mod_scaled$pc[predict_mod_scaled$Data_Split == "Test"])

###
###     Preparing variables in dataframe for analysis
###

  # US year
predict_mod_scaled$us_year = as.numeric(substr(x = predict_mod_scaled$date_of_current_us,start = 1,stop = 4))
predict_mod_scaled$us_year[predict_mod_scaled$us_year < 1990 | predict_mod_scaled$us_year > 2020] = NA

  # Patient age 
predict_mod_scaled$pt_age = elapsed_months(as.Date(predict_mod_scaled$date_of_current_us), as.Date(predict_mod_scaled$date_of_us1)) + predict_mod_scaled$age_at_baseline
predict_mod_scaled$pt_age_bin = NA
predict_mod_scaled$pt_age_bin[predict_mod_scaled$pt_age < 12] = "<1"
predict_mod_scaled$pt_age_bin[predict_mod_scaled$pt_age > 12 & predict_mod_scaled$pt_age < 60] = "1-5"
predict_mod_scaled$pt_age_bin[predict_mod_scaled$pt_age > 60] = ">5"

table(predict_mod_scaled$pt_age_bin[predict_mod_scaled$Data_Split == "Training"])
table(predict_mod_scaled$pt_age_bin[predict_mod_scaled$Data_Split == "Test"])

  # Distance from hospital 
predict_mod_scaled$postal_code = substr(x = predict_mod_scaled$pc,1,3)
table(predict_mod_scaled$postal_code)
predict_mod_scaled$pc_area = NA
predict_mod_scaled$pc_area[substr(predict_mod_scaled$postal_code,1,1) == "M"] = "Toronto"
predict_mod_scaled$pc_area[substr(predict_mod_scaled$postal_code,1,1) == "L"] = "GTA suburbs"
predict_mod_scaled$pc_area[substr(predict_mod_scaled$postal_code,1,1) == "K"] = "Eastern ON"
predict_mod_scaled$pc_area[substr(predict_mod_scaled$postal_code,1,1) == "P"] = "Northern ON"
predict_mod_scaled$pc_area[substr(predict_mod_scaled$postal_code,1,1) == "S"] = "SK"

table(predict_mod_scaled$pc_area[predict_mod_scaled$Data_Split == "Training"])
table(predict_mod_scaled$pc_area[predict_mod_scaled$Data_Split == "Test"])

  # Manufacturer labels
predict_mod_scaled$manu.f = factor(predict_mod_scaled$manu,
                                   levels = c("philips-medical-systems",
                                              "toshiba-mec-us","samsung-medison-co-ltd", 
                                              "atl","acuson",          
                                              "ge-medical-systems","siemens",                
                                              "ge-healthcare","toshiba-mec"),
                                   labels = c("Philips",
                                              "Toshiba1","Samsung", 
                                              "ATL","Acuson",          
                                              "GE1","Siemens",                
                                              "GE2","Toshiba2"))
  # Raw error
predict_mod_scaled$error = abs(predict_mod_scaled$Scaled_Pred - predict_mod_scaled$Target)

  # 0.05 thresholded predicted values
predict_mod_scaled$pred_target = NA
predict_mod_scaled$pred_target[predict_mod_scaled$Scaled_Pred < 0.05] = 0
predict_mod_scaled$pred_target[predict_mod_scaled$Scaled_Pred > 0.05] = 1

  # Provider factor
predict_mod_scaled$prov.f = factor(predict_mod_scaled$provider, levels = c(1,2,3,4,5),
                                   labels = c("surgeon","fellow","residents","NP","pediatrician"))


  ##
  ##      OVERALL SUMMARY STATISTICS
  ##

dim(predict_mod_scaled[predict_mod_scaled$Data_Split == "Training",])
length(unique(predict_mod_scaled$study_id[predict_mod_scaled$Data_Split == "Training"]))

table(predict_mod_scaled$Target[predict_mod_scaled$Data_Split == "Training"])/sum(table(predict_mod_scaled$Target[predict_mod_scaled$Data_Split == "Training"]))
table(predict_mod_scaled$gender[predict_mod_scaled$Data_Split == "Training"])/sum(table(predict_mod_scaled$gender[predict_mod_scaled$Data_Split == "Training"]))
range(na.omit(predict_mod_scaled$us_year[predict_mod_scaled$Data_Split == "Training"]))

dim(predict_mod_scaled[predict_mod_scaled$Data_Split == "Test",])
length(unique(predict_mod_scaled$study_id[predict_mod_scaled$Data_Split == "Test"]))
table(predict_mod_scaled$Target[predict_mod_scaled$Data_Split == "Test"])/sum(table(predict_mod_scaled$Target[predict_mod_scaled$Data_Split == "Test"]))
table(predict_mod_scaled$gender[predict_mod_scaled$Data_Split == "Test"])/sum(table(predict_mod_scaled$gender[predict_mod_scaled$Data_Split == "Test"]))
range(na.omit(predict_mod_scaled$us_year[predict_mod_scaled$Data_Split == "Test"]))

  # sex split
table(predict_mod_scaled$gender[predict_mod_scaled$Data_Split == "Training"])/sum(table(predict_mod_scaled$gender[predict_mod_scaled$Data_Split == "Training"]))
table(predict_mod_scaled$gender[predict_mod_scaled$Data_Split == "Test"])/sum(table(predict_mod_scaled$gender[predict_mod_scaled$Data_Split == "Test"]))

  # Target summary
theme_set(
  theme_classic(base_size = 20)
)
ggplot(predict_mod_scaled,aes(x = Target.f,y = Scaled_Pred,fill = Target.f)) + 
  geom_violin() + geom_hline(aes(yintercept = 0.05)) + 
  facet_grid(~Data_Split) + xlab("Target") + ylab("Predicted p(surgery)") + scale_fill_viridis(discrete = TRUE,name = "Target")

##
##      ANALYSIS OF BIAS IN TEST SET SCORES
##

## VARIABLES TO ASSESS BIAS IN: 
#   - patient age 
#   - patient sex
#   - which provider ordered surgery*
#   - distance patient lives from the hospital*
#   - ultrasound machine
#   - year

# age -- update to age at US, not just age at baseline
age_auc_test_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Test",var_of_int = "pt_age_bin")
predict_mod_scaled[predict_mod_scaled$pt_age_bin == "1-5" & predict_mod_scaled$Data_Split == "Test",]

age_auc_train_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Training",var_of_int = "pt_age_bin")
age_auc_train_df

predict_mod_scaled[predict_mod_scaled$pt_age_bin == ">5" & predict_mod_scaled$Data_Split == "Training",]

# sex
sex_auc_test_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Test",var_of_int = "gender")
sex_auc_test_df

sex_auc_train_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Training",var_of_int = "gender")
sex_auc_train_df


# machine
manu_auc_train_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Training",var_of_int = "manu")
manu_auc_test_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Test",var_of_int = "manu")

theme_set(
  theme_classic(base_size = 20)
)

ggplot(predict_mod_scaled,aes(x = us_year,fill = manu.f)) + 
  geom_bar(position = "fill") + xlab("Year") + ylab("Proportion") + scale_fill_viridis(discrete = TRUE,name = "Ultrasound\nMachine")# + theme_bw()

# US year
us_yr_auc_test_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Test",var_of_int = "us_year")
us_yr_auc_test_df

us_yr_auc_train_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Training",var_of_int = "us_year")
us_yr_auc_train_df


# US number
boxplot(error ~ us_num,predict_mod_scaled,las = 2,xlab="")

predict_mod_scaled[predict_mod_scaled$Data_Split == "Test" & predict_mod_scaled$us_num == 4,]

us_num_auc_train_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Training",var_of_int = "us_num")
us_num_auc_train_df
us_num_auc_test_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Test",var_of_int = "us_num")
us_num_auc_test_df

theme_set(
  theme_classic(base_size = 20)
)
ggplot(predict_mod_scaled[predict_mod_scaled$Data_Split == "Test" & predict_mod_scaled$us_num < 5,],aes(x = Target.f,y = Scaled_Pred, fill = Target.f)) + 
  geom_violin() + geom_hline(aes(yintercept = 0.05)) + facet_grid(~us_num) + xlab("Target") + ylab("Predicted p(surgery)") + 
  scale_fill_viridis(discrete = TRUE,name = "Target")

# Postal code
us_pc_auc_test_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Test",var_of_int = "pc_area")
us_pc_auc_test_df

us_pc_auc_train_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Training",var_of_int = "pc_area")
us_pc_auc_train_df

temp = predict_mod_scaled[predict_mod_scaled$Data_Split == "Training" & predict_mod_scaled$pc_area == "Northern ON",]
temp[!is.na(temp$pc_area),]

table(predict_mod_scaled$pc_area,
      predict_mod_scaled$Target.f,
      predict_mod_scaled$Data_Split)

table(predict_mod_scaled$pc_area,
      predict_mod_scaled$Data_Split)


# Provider
table(predict_mod_scaled$provider)
prov_auc_test_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Test",var_of_int = "prov.f")
prov_auc_test_df

prov_auc_train_df = get_auc_df(score_df = predict_mod_scaled,data_split = "Training",var_of_int = "prov.f")
prov_auc_train_df

table(predict_mod_scaled$prov.f,
      predict_mod_scaled$Target.f,
      predict_mod_scaled$Data_Split)

table(predict_mod_scaled$prov.f,
      predict_mod_scaled$Data_Split)

##***************** EXTRA UN-USED CODE *****************## 

##
##   GRAPHS OF DATA OVER TIME
##
#   - ultrasound number
ggplot(predict_mod_scaled[predict_mod_scaled$us_year > 2005,],aes(x = us_year,fill = Target.f)) + 
  geom_bar(position = "fill")

ggplot(predict_mod_scaled[predict_mod_scaled$us_year > 2005,],aes(x = us_year,fill = gender)) + 
  geom_bar(position = "fill")


# ggplot(predict_mod_scaled[predict_mod_scaled$us_year > 2005,],aes(x = us_year,fill = manu.f)) + 
#   geom_bar(position = "fill") + xlab("Year") + ylab("Proportion") + theme_bw() + scale_fill_discrete(name = "Manufacturer")
# 
# ggplot(predict_mod_scaled,aes(x = us_year,fill = manu)) + 
#   geom_bar(position = "fill") + xlab("Year") + ylab("Proportion") + theme_bw()
# 
# predict_mod_scaled$us_num.f = factor(predict_mod_scaled$us_num,levels = 1:10)
# 
# ggplot(predict_mod_scaled[!is.na(predict_mod_scaled$us_num.f) & predict_mod_scaled$us_year > 2005,],aes(x = us_year,fill = us_num.f)) + 
#   geom_bar(position = "fill") + xlab("Year") + ylab("Proportion") + theme_bw()


##
##  LOGISTIC REGRESSION (QUASI-BINOMIAL) TO PREDICT ERROR LEVEL

error_mod = glm(error ~ us_num + us_year + factor(manu) + gender + age_at_baseline + Target,
                data = predict_mod_scaled, family = "quasi")
summary(error_mod)

## graph showing variable value proportion change by year 
vars_to_graph = c("us_num.f","manu","gender")

grid_df = data.frame(matrix(ncol = 3,nrow = 0))
# names(grid_df) = c("year","variable","proportion")

my_var = vars_to_graph[1]
for(my_var in vars_to_graph){
  my_tbl = table(predict_mod_scaled[predict_mod_scaled$us_year > 2005,my_var],
                 predict_mod_scaled[predict_mod_scaled$us_year > 2005,"us_year"])
  mlt_tbl = melt(my_tbl/sum(my_tbl))  
  grid_df = rbind(grid_df,mlt_tbl)  
}
names(grid_df) = c("Variable","Year","value")

ggplot(grid_df[grid_df$Year > 2015,],aes(x = Year,y = Variable, fill = value)) + geom_tile() + scale_fill_gradient(low = "dodgerblue1",high = "Red")
