#telescopic_HDP_NNIG_uni_VariantwithDP
#VARIANT WITH DP AT FIRST LAYER 
library(progress)#to draw the progress bar

telescopic_HDP_NNIG_uni_VariantwithDP = function(data, hyper=c(0.1, 0.1, 0.1, 0.1, 0.1),
                                   alpharandom = FALSE,
                                   H0 = 5, H = 5, totiter = 1000){
  #has to work for many X and one param
  kernel_eval_joint = function(x, hyper1, hyper2){
    return(sum(dnorm(x, hyper1, sqrt(hyper2), log = TRUE)))
  }
  
  #has to work for one X and many param
  kernel_eval = function(x, hyper1, hyper2){
    return(dnorm(x, hyper1, sqrt(hyper2), log = TRUE))
  }
  
  posterior_sampler = function(x, k0 = 0.1, a = 0.1, b = 0.1){
    n_temp = length(x)
    if(n_temp>0){
      xbar = mean(x)
      sum_squared = sum((x - xbar)^2)
      sigma = 1 / rgamma(1, a + n_temp / 2, b + sum_squared/2 +
                           xbar^2/2 * n_temp*k0 / (k0+n_temp))
      mu = rnorm(1, mean = n_temp / (n_temp + k0) * xbar,
                 sd = sqrt(sigma / (n_temp+k0)))
      return( list(mu = mu, sigma = sigma) )
    }else{
      sigma = 1 / rgamma(1, a , b )
      mu = rnorm(1, mean = 0,
                 sd = sqrt(sigma / k0))
      return( list(mu = mu, sigma = sigma) )
    }
  }
  
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
  sigma0 = matrix(NA, nrow = H0, ncol = L) #atoms #var
  
  #PARTITION 
  m = matrix(NA, nrow = L, ncol = n_tot) #clustering configuration
  c = array(NA, dim = c(L, n_tot)) #auxiliary: tables per item
  k = array(0, dim = c(L, H0*H)) #auxiliary: dishes per table
  q = matrix(NA, nrow = L, ncol = H0) #dish freq
  n = matrix(NA, nrow = L, ncol = H*H0)#table freq
  nbar = array(NA, dim = c(L, H0, H))#table freq remapped
  
  m_saved = array(NA,c(totiter, L, n_tot))
  
  #initialize
  for(l in 1:L){
    m[l, ] = stats::kmeans(as.matrix(X[,l]),min(10,H0))$cluster
  }
  k1 = m[1, ] #tables and cluster are the same at layer 1
  for(l in 2:L){
    c[l, ] = (m[l-1,]-1)*H + m[l, ]
  }
  for(l in 2:L){
    for (cc in 1:(H*H0)){
      if (cc %in% c[l,]){
        k[l,cc] = m[l,c[l,]==cc][1]
      }else{
        k[l,cc] = sample(1:H0, 1)
      }
    }
  }
  
  
  q[,] = 0
  #layer 1 = DP
  
  for(h in 1:H0){
    q[1,h] = sum((k1==h)) #num of items in each cluster at layer 1
  }
  
  for(l in 2:L){
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
    k0 = hyper[3]; a_sigma = hyper[4]; b_sigma = hyper[5]
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
          out = posterior_sampler(X[m[l,]==h1,l], k0, a_sigma, b_sigma)
          theta0[h1, l] = out$mu
          sigma0[h1, l] = out$sigma
          
          nbar_sum = c(rev(cumsum(rev(nbar[l,h1,])))[2:H],0)
          if(l>1){
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
      }
      #for (l in 1:L){
      for (l in 2:L){ #sample table for all layers except 1
        for (i in 1:n_tot){
          past = m[l-1, i]
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
          prob = (pim[past,,l]) + kernel_eval(X[i,l], theta0[ ,l], sigma0[, l]) + fut
          if(max(prob)==-Inf){prob[]=1}
          prob = prob - max(prob)
          prob = exp(prob) 
          if(sum(exp(fut))!=0){ #if fut is zero it cannot move
            c[l, i] = sample(((past-1)*H+1):(past*H), 1, prob = prob)
          }
        }
      }
      
      #first layer
      for (i in 1:n_tot){
        fut = rep(0,H0)
        for (kk in 1:H0){
          for (d in 1:H0){
            fut[kk] = fut[kk] + (k[2, (kk-1)*H + d] == m[2,i])*exp(pim[kk,d,2])
            #print (fut[kk])
          } 
        }
        fut = log(fut)
        prob = pi0[,1] + kernel_eval(X[i,1], theta0[,1], sigma0[, 1]) + fut  
        if(max(prob)==-Inf){prob[]=1}
        prob = prob - max(prob)
        prob = exp(prob) 
        if(sum(exp(fut))!=0){ #if fut is zero it cannot move
          k1[i] = sample(1:H0, 1, prob = prob)
        }
      }
      for (i in 1:n_tot){
        m[1,i] = k1[i]
      }
      
      for (l in 2:L){  
        for (cc in 1:(H*H0)){
          prob = NULL
          for(kk in 1:H0){
            prob[kk] = pi0[kk,l] + kernel_eval_joint(X[c[l,]==cc,l], theta0[kk,l], sigma0[kk,l])   #check what happens if there is no one in cc
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
      #layer 1 = DP
      
      for(h in 1:H0){
        q[1,h] = sum((k1==h)) #num of items in each cluster at layer 1
      }
      
      for(l in 2:L){
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
      #print(c[1,])
      m_saved[iter,,] = m
    }
    return(m_saved)
  }else{
    alpha0 = rep(hyper[1],L); alpha = rep(hyper[2],L)
    k0 = hyper[3]; a_sigma = hyper[4]; b_sigma = hyper[5]
    pb <- progress_bar$new(
      format = " MCMC [:bar] :percent Estimated completion time: :eta",
      total = totiter, clear = FALSE, width= 100)
    
    #mcmc
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
        out = posterior_sampler(X[m[l,]==h1,l], k0, a_sigma, b_sigma)
        theta0[h1, l] = out$mu
        sigma0[h1, l] = out$sigma
        
        nbar_sum = c(rev(cumsum(rev(nbar[l,h1,])))[2:H],0)
        if(l>1){
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
    }
    #for (l in 1:L){
    for (l in 2:L){ #sample table for all layers except 1
      for (i in 1:n_tot){
        past = m[l-1, i]
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
        prob = (pim[past,,l]) + kernel_eval(X[i,l], theta0[ ,l], sigma0[, l]) + fut
        if(max(prob)==-Inf){prob[]=1}
        prob = prob - max(prob)
        prob = exp(prob) 
        if(sum(exp(fut))!=0){ #if fut is zero it cannot move
          c[l, i] = sample(((past-1)*H+1):(past*H), 1, prob = prob)
        }
      }
    }
    
    #first layer
    for (i in 1:n_tot){
      fut = rep(0,H0)
      for (kk in 1:H0){
        for (d in 1:H0){
          fut[kk] = fut[kk] + (k[2, (kk-1)*H + d] == m[2,i])*exp(pim[kk,d,2])
          #print (fut[kk])
        } 
      }
      fut = log(fut)
      prob = pi0[,1] + kernel_eval(X[i,1], theta0[,1], sigma0[, 1]) + fut  
      if(max(prob)==-Inf){prob[]=1}
      prob = prob - max(prob)
      prob = exp(prob) 
      if(sum(exp(fut))!=0){ #if fut is zero it cannot move
        k1[i] = sample(1:H0, 1, prob = prob)
      }
    }
    for (i in 1:n_tot){
      m[1,i] = k1[i]
    }
    
    for (l in 2:L){  
      for (cc in 1:(H*H0)){
        prob = NULL
        for(kk in 1:H0){
          prob[kk] = pi0[kk,l] + kernel_eval_joint(X[c[l,]==cc,l], theta0[kk,l], sigma0[kk,l])   #check what happens if there is no one in cc
        }
        if(max(prob)==-Inf){prob[]=1}
        prob = prob - max(prob)
        prob = exp(prob) 
        k[l, cc] = sample(1:H0, 1, prob = prob)
      }
      for (i in 1:n_tot){
        m[l,i] = k[l,c[l,i]]
      }
    }
    
    q[,] = 0
    #layer 1 = DP
    
    for(h in 1:H0){
      q[1,h] = sum((k1==h)) #num of items in each cluster at layer 1
    }
    
    for(l in 2:L){
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
