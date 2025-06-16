#t-HDP _ conditional algorithm - NN model - univariate
library(progress)#to draw the progress bar

#markovian dependence for L layers, univariate layers
#normal kernel and normal base,
#var of the kernel is fixed to 1 and the mean of the base to 0
#hyper are the two concentration parameters and the variance of the base
#if alpharandom=TRUE the two first entries of hyper are only the initialization
#data is the dataset nXL 
#H0 and H are the number of mixture components used in the approximation 
#of the infinite mixtures
telescopic_HDP_NN_uni = function(data, hyper=c(0.1, 0.1, 0.1), 
                                            alpharandom = FALSE,
                                            H0 = 5, H = 5, totiter = 1000){
  #has to work for many X and one param
  kernel_eval_joint = function(x, hyper){
    return(sum(dnorm(x, hyper, 1, log = TRUE)))
  }
  
  #has to work for one X and many param
  kernel_eval = function(x, hyper){
    return(dnorm(x, hyper, 1, log = TRUE))
  }
  
  posterior_sampler = function(x, s2 = 0.1){
    n_temp = length(x)
    if(n_temp>0){
      xbar = mean(x)
      return(rnorm(1, mean = n_temp*s2 / (n_temp*s2+1) * xbar,
                   sd = sqrt(s2 / (n_temp*s2 + 1))) )
    }else{
      return(rnorm(1,0, sqrt(s2)))
    }
  }
  
  s2 = hyper[3] #variance of the base measure
  L = dim(data)[2] #number of layers
  n_tot = dim(data)[1] #number of items
  
  #the data in a matrix X, which is nxL
  X = data 
  
  #RANDOM PROB.S
  pim = array(NA, dim = c(H0, H, L)) #weights (in log scale) 
  b = array(NA, dim = c(H0, H, L)) #sticks for pi
  
  pi0 = matrix(NA, nrow = H0, ncol = L) #weights of the common overall in HDP
  b0 = matrix(NA, nrow = H0, ncol = L) #sticks for pi0 (in log scale)
  theta0 = matrix(NA, nrow = H0, ncol = L) #atoms
  
  #PARTITION 
  m = matrix(NA, nrow = L, ncol = n_tot) #clustering configuration
  c = array(NA, dim = c(L, n_tot)) #auxiliary: tables per item
  #c[l,i] is the table where item i sits, 
  #the table labels are distinct across restaurants at different layers
  #note that there are at most H tables in in each restaurant
  #and at most H0 restaurants, 
  #thus the maximum number of tables is H*H0 at each layer
  k = array(0, dim = c(L, H0*H)) #auxiliary: dishes per table 
  #k[l,t] is the dish served at table t at layer l 
  q = matrix(NA, nrow = L, ncol = H0) #dish freq
  #q[l,h] how many tables at layer l serve dish h
  n = matrix(NA, nrow = L, ncol = H*H0) #table freq
  #n[l,t] how many items sits at table t at layer l
  nbar = array(NA, dim = c(L, H0, H))#table freq remapped
  #nbar[l,r,t] how many items sits at the t-th table in restaurant r at layer l
  #those are remapped in the sense the here t is repeated across restaurants
  
  m_saved = array(NA,c(totiter, L, n_tot)) #cluster allocation
  
  #initialize
  for(l in 1:L){
    m[l, ] = stats::kmeans(as.matrix(X[,l]),min(10,H0))$cluster
  }
  c[1, ] = m[1, ] #at first layer set the tables to the cluster allocation, 
  #cause there is a single restaurant
  for(l in 2:L){ #at subsequent layer set tables to the finer cluster allocation,
    #such that at different restaurants we do not have same label for two tables, 
    #so that if two subjects eat the same dish but in two different restaurant 
    #they are sitted to different tables
    c[l, ] = (m[l-1,]-1)*H + m[l, ]
  }
  #finally set the dishes of each table so to obtain the cluster allocation
  for(l in 1:L){
    for (cc in 1:(H*H0)){
      if (cc %in% c[l,]){
        k[l,cc] = m[l,c[l,]==cc][1]
      }else{
        k[l,cc] = sample(1:H0, 1)
      }
    }
  }
  
  
  q[,] = 0
  for(l in 1:L){
    
    temp = rep(0, H*H0)
    temp[unique(c[l,])] = 1
    for(h in 1:H0){
      q[l,h] = sum((k[l,]==h)*temp)
    }
    
    for(cc in 1:(H*H0)){
      n[l,cc] = sum(c[l,]==cc)
    }
    
    for(h1 in 1:H0){
      for(h2 in 1:H){
        nbar[l, h1, h2] = n[l, (h1-1)*H + h2]
      }
    }
  }
  
  
  
  if(!alpharandom){
    alpha0 = hyper[1]; alpha = hyper[2]
    pb <- progress_bar$new(
      format = " MCMC [:bar] :percent Estimated completion time: :eta",
      total = totiter, clear = FALSE, width= 100)
    
    #mcmc
    for (iter in 1:totiter){
      pb$tick()
      for (l in 1:L){
        ql_sum = c(rev(cumsum(rev(q[l,])))[2:H0],0)
        for (h1 in 1:H0){
          b0[h1, l] = rbeta(1, 1 + q[l,h1], alpha0 + ql_sum[h1])
          if(h1>1){
            pi0[h1, l] = log(b0[h1, l]) + sum(log(1 - b0[1:(h1-1), l]))
          }else{
            pi0[h1, l] = log(b0[h1, l])
          }
          
          theta0[h1, l] = posterior_sampler(X[m[l,]==h1,l], s2)
          
          nbar_sum = c(rev(cumsum(rev(nbar[l,h1,])))[2:H],0)
          for(h2 in 1:H){
            b[h1, h2, l] = rbeta(1, 1 + nbar[l,h1,h2], alpha + nbar_sum[h2])
            if(h2>1){
              pim[h1, h2, l] = log(b[h1, h2, l]) + sum(log(1 - b[h1, 1:(h2-1),l]))
            }else{
              pim[h1, h2, l] = log(b[h1, h2, l])
            }
          }
        }
      }
      for (l in 1:L){
        for (i in 1:n_tot){
          if(l>1){past = m[l-1, i]}else{past = 1}
          fut = rep(0,H)
          if(l<L){
            for (cc in 1:H){
              for (d in 1:H0){
                fut[cc] = fut[cc] + 
                  (k[l+1,(k[l,cc]-1)*H + d] == m[l+1,i])*exp(pim[k[l,cc],d,l+1])
              } 
            }
            fut = log(fut)
          }
          prob = (pim[past,,l]) + kernel_eval(X[i,l], theta0[ ,l]) + fut
          if(max(prob)==-Inf){prob[]=1}
          prob = prob - max(prob)
          prob = exp(prob) 
          if(sum(exp(fut))!=0){ #if fut is zero it cannot move
            c[l, i] = sample(((past-1)*H+1):(past*H), 1, prob = prob)
          }
        }
        for (cc in 1:(H*H0)){
          prob = NULL
          for(kk in 1:H0){
            prob[kk] = pi0[kk,l] + kernel_eval_joint(X[c[l,]==cc,l], theta0[kk,l])   #check what happens if there is no one in cc
          }
          prob = prob - max(prob)
          prob = exp(prob) 
          k[l, cc] = sample(1:H0, 1, prob = prob)
        }
        for (i in 1:n_tot){
          m[l,i] = k[l,c[l,i]]
        }
      } 
      q[,] = 0
      for(l in 1:L){
        
        temp = rep(0, H*H0)
        temp[unique(c[l,])] = 1
        for(h in 1:H0){
          q[l,h] = sum((k[l,]==h)*temp)
        }
        
        for(cc in 1:(H*H0)){
          n[l,cc] = sum(c[l,]==cc)
        }
        
        for(h1 in 1:H0){
          for(h2 in 1:H){
            nbar[l, h1, h2] = n[l, (h1-1)*H + h2]
          }
        }
      }
      #print(c(iter, rand.index(m[1,],true_layer1),rand.index(m[2,],true_layer2) ) )
      m_saved[iter,,] = m
    }
    return(m_saved)
  }else{
    alpha0 = rep(hyper[1],L); alpha = rep(hyper[2],L)
    pb <- progress_bar$new(
      format = " MCMC [:bar] :percent Estimated completion time: :eta",
      total = totiter, clear = FALSE, width= 100)
    
    #mcmc
    for (iter in 1:totiter){
      pb$tick()
      for (l in 1:L){
        ql_sum = c(rev(cumsum(rev(q[l,])))[2:H0],0)
        for (h1 in 1:H0){
          b0[h1, l] = rbeta(1, 1 + q[l,h1], alpha0[l] + ql_sum[h1])
          if(h1>1){
            pi0[h1, l] = log(b0[h1, l]) + sum(log(1 - b0[1:(h1-1), l]))
          }else{
            pi0[h1, l] = log(b0[h1, l])
          }
          
          theta0[h1, l] = posterior_sampler(X[m[l,]==h1,l], s2)
          
          nbar_sum = c(rev(cumsum(rev(nbar[l,h1,])))[2:H],0)
          for(h2 in 1:H){
            b[h1, h2, l] = rbeta(1, 1 + nbar[l,h1,h2], alpha[l] + nbar_sum[h2])
            if(h2>1){
              pim[h1, h2, l] = log(b[h1, h2, l]) + sum(log(1 - b[h1, 1:(h2-1),l]))
            }else{
              pim[h1, h2, l] = log(b[h1, h2, l])
            }
          }
        }
      }
      for (l in 1:L){
        for (i in 1:n_tot){
          if(l>1){past = m[l-1, i]}else{past = 1}
          fut = rep(0,H)
          if(l<L){
            for (cc in 1:H){
              for (d in 1:H0){
                fut[cc] = fut[cc] + 
                  (k[l+1,(k[l,cc]-1)*H + d] == m[l+1,i])*exp(pim[k[l,cc],d,l+1])
              } 
            }
            fut = log(fut)
          }
          prob = (pim[past,,l]) + kernel_eval(X[i,l], theta0[ ,l]) + fut
          if(max(prob)==-Inf){prob[]=1}
          prob = prob - max(prob)
          prob = exp(prob) 
          if(sum(exp(fut))!=0){ #if fut is zero it cannot move
            c[l, i] = sample(((past-1)*H+1):(past*H), 1, prob = prob)
          }
        }
        for (cc in 1:(H*H0)){
          prob = NULL
          for(kk in 1:H0){
            prob[kk] = pi0[kk,l] + kernel_eval_joint(X[c[l,]==cc,l], theta0[kk,l])   #check what happens if there is no one in cc
          }
          prob = prob - max(prob)
          prob = exp(prob) 
          k[l, cc] = sample(1:H0, 1, prob = prob)
        }
        for (i in 1:n_tot){
          m[l,i] = k[l,c[l,i]]
        }
      } 
      q[,] = 0
      for(l in 1:L){
        
        temp = rep(0, H*H0)
        temp[unique(c[l,])] = 1
        for(h in 1:H0){
          q[l,h] = sum((k[l,]==h)*temp)
        }
        
        for(cc in 1:(H*H0)){
          n[l,cc] = sum(c[l,]==cc)
        }
        
        for(h1 in 1:H0){
          for(h2 in 1:H){
            nbar[l, h1, h2] = n[l, (h1-1)*H + h2]
          }
        }
        alpha[l] = rgamma(1, 3 + H*H0, rate = 3 - sum(log(1-b[,,l])))
        if(alpha[l]<0.1){alpha[l]=0.1}
        alpha0[l] = rgamma(1, 3 + H0, rate = 3 - sum(log(1-b0[,l])))
        if(alpha0[l]<0.1){alpha0[l]=0.1}
      }
      
      #print(c(iter, rand.index(m[1,],true_layer1),rand.index(m[2,],true_layer2) ) )
      m_saved[iter,,] = m
    }
    return(m_saved) 
  }
}