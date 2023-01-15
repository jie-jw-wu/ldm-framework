trunc_weights <- function(weights) {
    x <- cbind(x1 = trunc(weights[1]*100)/100, x2=trunc(weights[2]*100)/100, x3=1.000-trunc(weights[1]*100)/100-trunc(weights[2]*100)/100)

    if (!identical(rowSums(x), 1)) {
      x[1,3]=x[1,3] - (rowSums(x)-1)
    }
    
    return(c(x[1,1],x[1,2],x[1,3]))
}

# get weights for clicks, impressions, and click rate
# @p pairwise comparision matrix used for calculating the weights
# @weight_method the criteria weighting method used
get_weights <- function(p, weight_method){
    if (weight_method == 1) {
        # Mean Weight (MW)
        return(c(1/3,1/3,1/3))
    } else if (weight_method == 2) {
        # Standard Deviation Method
        tao = apply(p, 2, sd)
        weights = tao/sum(tao)
        if (sum(weights)!=1) {
            return(trunc_weights(weights))
        }
        
        return (weights)
    } else if (weight_method == 3) {
        # Statistical Deviation Method
        avg = mean(p)
        n = nrow(p)
        weights = c(0,0,0)
        for (j in 1:3) {
            for (i in 1:n) {
                weights[j] = weights[j] + (p[i,j]-avg)*(p[i,j]-avg)/n
            }
        }
        sum_w = sum(weights)
        weights = weights / sum_w
        if (sum(weights)!=1) {
            return(trunc_weights(weights))
        }
        return (trunc_weights(weights))#(weights)
    } else if (weight_method == 4) {
        # Entropy Method
        m = nrow(p)
        n = 3
        new_p = p+matrix(1,m,n)
        p_normalized = (new_p) / sum(new_p)
        e = rep(0, n)
        for (j in 1:n) {
            for (i in 1:m) {
                if (abs(p_normalized[i,j])>= 0.00001) {
                e[j] = e[j] - p_normalized[i,j]*log(p_normalized[i,j])/log(m)
                }
            }
        }
        w = (rep(1, n) - e) / (n - sum(e))
        if (sum(w)!=1) {
            return(trunc_weights(w))
        }
        return (w) 
    } else if (weight_method == 5) {
        # SMART
        return (c(55/165, 10/165, 100/165))
    } else if (weight_method == 6) {
        # Ranking sum
        return (c(2/6, 1/6, 3/6))
    } else if (weight_method == 7) {
        # Ranking Recipocal
        weights = c(0.5/(1+0.5+0.33), 0.33/(1+0.5+0.33), 1/(1+0.5+0.33))
        
        if (sum(weights)!=1) {
            return(trunc_weights(weights))
        }
        return (weights)
    } else if (weight_method == 8) {
        # Ranking exponent
        return (c(4/14, 1/14, 9/14)) 
    } else if (weight_method == 9) {
        # AHP
        return (c(0.221, 0.05, 0.729))
    } else {
        return (c(1/3,1/3,1/3))
    }
}

get_change_pct <- function(a,b){
    if (is.null(b)) {
        return ((a-0.001) / 0.001)
    } else if (b == 0) {
        return ((a-0.001) / 0.001)
    } else {
        return ((a-b)/b)
    }
}

get_pairwise_comparison <- function(test_group, control){
    p = matrix(nrow=length(test_group), ncol=3)
    for (i in 1:length(test_group)) {
        t = test_group[[i]]
        c = test_group[[control]]
        p[i,1] = get_change_pct(t$clicks, c$clicks)
        p[i,2] = get_change_pct(t$impressions, c$impressions)
        p[i,3] = get_change_pct(t$click_rate, c$click_rate)
    }
    return(p)
}

# do criteria weighting and analysis of alternative (MCDM)
# return a array where a[i] means whether (1 or 0) there's a correct predicted decision in top i results.
apply_ddm <- function(test_group, control, num_result, weight_method, mcdm_method) {
    p = get_pairwise_comparison(test_group, control)
    n = length(test_group)
    a = data.frame("ai"=c(),"first_place"=c())
    weights = get_weights(p, weight_method)
    
    if (mcdm_method == 1) { 
        # WSM
        for (i in 1:n){
            ai = 0
            for (j in 1:3){
                ai = ai + p[i,j] * weights[i]
            }
            a=rbind(a, list("ai"=ai, "first_place"=test_group[[i]]$first_place))
        }
    } else if (mcdm_method == 2) {
        # TOPSISLinear
        #r=TOPSISLinear(p,c(x[1,1],x[1,2],x[1,3]),rep('max', 3))
        r=TOPSISLinear(p,weights,rep('max', 3)) # changed on 6/5/2022. it may change the results, TODO: verify if the result changes
      
        for (i in 1:n){
            if (i!=r$Alternatives[i]) {
                print("error i!=r$Alternatives[i]")
            } else {
                a=rbind(a, list("ai"=n-r$Ranking[i], "first_place"=test_group[[i]]$first_place))
            }
        }
    } else if (mcdm_method == 3) {
        # TOPSISVector
        r=TOPSISVector(p,weights,rep('max', 3))
        for (i in 1:n){
            if (i!=r$Alternatives[i]) {
                print("error i!=r$Alternatives[i]")
            } else {
                a=rbind(a, list("ai"=n-r$Ranking[i], "first_place"=test_group[[i]]$first_place))
            }
        }
    } else if (mcdm_method == 4) {
        # MMOORA
        r=MMOORA(p,weights,rep('max', 3))
        for (i in 1:n){
            if (i!=r$Alternatives[i]) {
                print("error i!=r$Alternatives[i]")
            } else {
                a=rbind(a, list("ai"=n-r$Ranking[i], "first_place"=test_group[[i]]$first_place))
            }
        }
    } else if (mcdm_method == 5) {
        # VIKOR
        v=0.5
        r=VIKOR(p,weights,rep('max', 3),v)
        for (i in 1:n){
            if (i!=r$Alternatives[i]) {
                print("error i!=r$Alternatives[i]")
            } else {
                a=rbind(a, list("ai"=n-r$Ranking[i], "first_place"=test_group[[i]]$first_place))
            }
        }
    } else if (mcdm_method == 6) {
        # WPM
        r=WASPAS(p,weights,rep('max', 3),0)
        for (i in 1:n){
            if (i!=r$Alternatives[i]) {
                print("error i!=r$Alternatives[i]")
            } else {
                a=rbind(a, list("ai"=n-r$Ranking[i], "first_place"=test_group[[i]]$first_place))
            }
        }
    } else if (mcdm_method == 7) {
        # WSM
        r=WASPAS(p,weights,rep('max', 3),1)
        for (i in 1:n){
            if (i!=r$Alternatives[i]) {
                print("error i!=r$Alternatives[i]")
            } else {
                a=rbind(a, list("ai"=n-r$Ranking[i], "first_place"=test_group[[i]]$first_place))
            }
        }
        # paper:  "Optimization of Weighted Aggregated Sum Product Assessment"
    }
    
    ordered_a = a[order(-a$ai),]
    
    res = rep(1, num_result)
    for (j in 1:min(num_result, n)) {
        if (ordered_a$first_place[j]==FALSE) {
            res[j] = 0
        } else {
            break
        }
    }
    return(res)
}



get_results <- function(weight_methods, mcdm_methods, min_idx = 0, max_idx = 30000) {
  
  knitr::opts_chunk$set(echo = TRUE)
  library(tidyverse)
  library(broom)
  library(MCDM)
  # Loading in the data
  rawdata <- read_csv("C:\\Users\\wuj\\Documents\\WJ\\Datasets\\upworthy-archive-exploratory-packages-03.12.2020.csv")
  
  reduced_data <- 
      rawdata %>% 
      select(clickability_test_id,
             clicks,
             impressions,
             #X1,
             winner,
             first_place)
  
  reduced_data$click_rate <- reduced_data$clicks/reduced_data$impressions
  grouped_data <- reduced_data %>% arrange(desc(clickability_test_id))

  id = grouped_data[1,"clickability_test_id"]
  num_result = 10
  accuracy = matrix(data=NA, nrow=length(weight_methods)*length(mcdm_methods), ncol=num_result)
  current = 1
  cat("min_idx, max_idx=",min_idx, max_idx, " nrow(grouped_data)=", nrow(grouped_data))
  
  weights_example = rep(0, num_result)
  for (weight_method in weight_methods) {
      for (mcdm_method in mcdm_methods) {
          test_group = c()
          num_group = 0
          num_correct = rep(0, num_result)
          for(i in 1:nrow(grouped_data)) {
            row <- grouped_data[i,]
            if ((id!=row$clickability_test_id) || (i==nrow(grouped_data))) {
                # ddm framework
                if (i >= min_idx && i < max_idx) {
                  control = 1
                  res = apply_ddm(test_group, control, num_result, weight_method, mcdm_method)
                  num_correct = num_correct + res
                  num_group = num_group + 1
                }
              
                test_group = c()
                id = row$clickability_test_id
            }
            test_group = c(test_group, list(row))
          }
          num_correct = num_correct/num_group
          for (k in 1:num_result) {
            accuracy[current,k] = num_correct[k] 
          }
          current=current+1
      }
  }

return(accuracy)
}


# 8 criteria weighting with fixed MCDM
plot_results_based_on_mcdm <- function(mcdm_method, title) {
  r1=get_results(c(1), c(mcdm_method))
  r2=get_results(c(2), c(mcdm_method))
  r3=get_results(c(3), c(mcdm_method))
  r4=get_results(c(4), c(mcdm_method))
  r5=get_results(c(5), c(mcdm_method))
  r6=get_results(c(6), c(mcdm_method))
  r7=get_results(c(7), c(mcdm_method))
  r8=get_results(c(8), c(mcdm_method))
  r9=get_results(c(9), c(mcdm_method))
  
  colors = c("black", "red", "green1", "blue", "cyan", "magenta", "gray", "orange", "green4")
  pchs = c(3,15,16,17,18,4,19,5,11)
  names = c("Mean Weight", "Standard Deviation","Statistical Deviation","Entropy","SMART","Ranking sum","Ranking Recipocal","Ranking exponent", "AHP")
  
  
  plot(1:10, r1, type='b', ylab='Precision @ K', xlab='Top K Deployment Decision Candidates', col=colors[1],pch=pchs[1], main=title, ylim=c(0.5,1))
  lines(1:10, r2, pch=pchs[2], col=colors[2], type='b')
  lines(1:10, r3, pch=pchs[3], col=colors[3], type='b')
  lines(1:10, r4, pch=pchs[4], col=colors[4], type='b')
  lines(1:10, r5, pch=pchs[5], col=colors[5], type='b')
  lines(1:10, r6, pch=pchs[6], col=colors[6], type='b')
  lines(1:10, r7, pch=pchs[7], col=colors[7], type='b')
  lines(1:10, r8, pch=pchs[8], col=colors[8], type='b')
  lines(1:10, r9, pch=pchs[9], col=colors[9], type='b')
  legend("bottomright", inset=.05, names, fill=colors)
}



# 6 MCDM with fixed Mean Weight
plot_results_based_on_criteria_weighting <- function(weight_method, title) {
  r1=get_results(c(weight_method), c(2))
  r2=get_results(c(weight_method), c(3))
  r3=get_results(c(weight_method), c(4))
  r4=get_results(c(weight_method), c(5))
  r5=get_results(c(weight_method), c(6))
  r6=get_results(c(weight_method), c(7))
  
  colors = c("black", "red", "green3", "blue", "cyan", "magenta", "gray")
  pchs = c(3,15,16,17,18,4,19)
  names = c("TOPSISLinear","TOPSISVector","MMOORA","VIKOR","WPM","WSM")
  
  ymin = 0.5 #min(r1[1],r2[1],r3[1],r4[1],r5[1],r6[1])-0.05
  plot(1:10, r1, type='b', ylab='Precision @ K', xlab='Top K Deployment Decision Candidates', col=colors[1],pch=pchs[1], main=title, ylim=c(ymin,1))
  lines(1:10, r2, pch=pchs[2], col=colors[2], type='b')
  lines(1:10, r3, pch=pchs[3], col=colors[3], type='b')
  lines(1:10, r4, pch=pchs[4], col=colors[4], type='b')
  lines(1:10, r5, pch=pchs[5], col=colors[5], type='b')
  lines(1:10, r6, pch=pchs[6], col=colors[6], type='b')
  legend("bottomright", names, fill=colors)
}


plot_results_based_on_cross_validation <- function(weight_method, title) {
  cv_result <- vector(mode="list", length=6)
  names(cv_result) <- c("TOPSISLinear","TOPSISVector","MMOORA","VIKOR","WPM","WSM")
  
  cv_acc_result <- vector(mode="list", length=6)
  names(cv_acc_result) <- c("TOPSISLinear","TOPSISVector","MMOORA","VIKOR","WPM","WSM")
  
  cv_acc_stats <- vector(mode="list", length=6)
  names(cv_acc_stats) <- c("TOPSISLinear","TOPSISVector","MMOORA","VIKOR","WPM","WSM")
  
  size = 1000
  fold_num = 21
  
  for (method in 1:6) {
    cv_result[[method]] <- rep(0,6)
    cv_acc_result[[method]] <- rep(0,0)
  }
  for (fold in 0:fold_num) {
    min_idx = fold * size
    max_idx = (fold + 1) * size
    if (fold == fold_num) {
      max_idx = (fold_num + 2) * size # merge the rest ~500 items in the fold 21
    }
    r1=get_results(c(weight_method), c(2), min_idx, max_idx)
    r2=get_results(c(weight_method), c(3), min_idx, max_idx)
    r3=get_results(c(weight_method), c(4), min_idx, max_idx)
    r4=get_results(c(weight_method), c(5), min_idx, max_idx)
    r5=get_results(c(weight_method), c(6), min_idx, max_idx)
    r6=get_results(c(weight_method), c(7), min_idx, max_idx)
    # sort 6 values and accumulate the num for the <mcdm method, position> pair
    
    df <- data.frame(
      method = 1:6,
      accuracy = c(r1[1,1], r2[1,1], r3[1,1], r4[1,1], r5[1,1], r6[1,1])
    )
    
    # add the accuracy for this fold 
    cv_acc_result[[1]] <- c(cv_acc_result[[1]], r1[1,1])
    cv_acc_result[[2]] <- c(cv_acc_result[[2]], r2[1,1])
    cv_acc_result[[3]] <- c(cv_acc_result[[3]], r3[1,1])
    cv_acc_result[[4]] <- c(cv_acc_result[[4]], r4[1,1])
    cv_acc_result[[5]] <- c(cv_acc_result[[5]], r5[1,1])
    cv_acc_result[[6]] <- c(cv_acc_result[[6]], r6[1,1])
    
    # sort the df by accuracy
    sorted_df <- df[
      with(df, order(accuracy, decreasing=TRUE)),
      ]
    
    rank = 1
    for (pos in 1:6) {
      method = sorted_df[pos,1]
      if (pos > 1 && sorted_df[pos-1,2] > sorted_df[pos,2]+0.01) {
        rank = rank + 1
      }
      cv_result[[method]][rank] = cv_result[[method]][rank] + 1 
    }
  }
  
  for (method in 1:6) {
    res=cv_acc_result[[method]]
    cv_acc_stats[[method]]=c(summary(res), sd(res))
  }
  return(cv_acc_stats)
}

############################################################################################################
# A
res_stats <- vector(mode="list", length=9)

for (weighting_method in 1:9) {
  res = plot_results_based_on_cross_validation(weighting_method, "CV of MDCM with Fixed Criteria Weighting")
  res_stats[[weighting_method]] = res
}
print(res_stats)

############################################################################################################
# B
# "Mean Weight", "Standard Deviation","Statistical Deviation","Entropy","SMART","Ranking Sum","Ranking Recipocal","Ranking Exponent", "AHP"
# plot_results_based_on_criteria_weighting(1,"Deployment Decision Results based on Mean Weight")
# plot_results_based_on_criteria_weighting(2,"Deployment Decision Results based on Standard Deviation")
# plot_results_based_on_criteria_weighting(3,"Deployment Decision Results based on Statistical Deviation")
# plot_results_based_on_criteria_weighting(4,"Deployment Decision Results based on Entropy")
# plot_results_based_on_criteria_weighting(5,"Deployment Decision Results based on SMART")
# plot_results_based_on_criteria_weighting(6,"Deployment Decision Results based on Ranking Sum")
# plot_results_based_on_criteria_weighting(7,"Deployment Decision Results based on Ranking Recipocal")
# plot_results_based_on_criteria_weighting(8,"Deployment Decision Results based on Ranking Exponent")
# plot_results_based_on_criteria_weighting(9,"Deployment Decision Results based on AHP")

############################################################################################################
# C
#"WSM bug", "TOPSISLinear","TOPSISVector","MMOORA","VIKOR","WPM","WSM"
# plot_results_based_on_mcdm(mcdm_method=2, title='Deployment Decision Results based on TOPSISLinear')
# plot_results_based_on_mcdm(mcdm_method=3, title='Deployment Decision Results based on TOPSISVector')
# plot_results_based_on_mcdm(mcdm_method=4, title="Deployment Decision Results based on MMOORA")
# plot_results_based_on_mcdm(mcdm_method=5, title="Deployment Decision Results based on VIKOR")
# plot_results_based_on_mcdm(mcdm_method=6, title="Deployment Decision Results based on WPM")
# plot_results_based_on_mcdm(mcdm_method=7, title="Deployment Decision Results based on WSM")

############################################################################################################
# instruction:
# 1. enable/uncomment only 1 of A or B or C
# 2. select All (ctrl+A), click Run button. (select all is needed to include most updated code in this file)
############################################################################################################
