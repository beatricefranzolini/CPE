#t-HDP _ conditional algorithm - NN model - univariate
library(progress)#to draw the progress bar

#triangular dependence for 3 layers, univariate layers
#normal kernel and normal-normal base
#H0 and H are the number of mixture components used in the approximation 
#of the infinite mixtures

telescopic_HDP_NN_uni_3L1P <- function(data, colayer, hyper=c(0.1, 0.1, 0.1),
                                       H0 = 10, H = 10, totiter = 1000){
  kernel_eval_joint = function(x, hyper){
    return(sum(dnorm(x, hyper, 1, log = TRUE)))
  }
  
  #has to work for one X and many param
  kernel_eval = function(x, hyper){
    return(rowSums(dnorm(x, hyper, 1, log = TRUE)))
  }
  
  posterior_sampler = function(x, s2 = 0.1){
    n_temp = length(x)
    if(n_temp>0){
      xbar = mean(x)
      return(rnorm(1, mean = n_temp*s2 / (n_temp*s2+1) * xbar,
                   sd = sqrt(1 / (n_temp*s2 + 1))) )
    }else{
      return(rnorm(1,0, sqrt(s2)))
    }
  }
  
  pb <- progress_bar$new(
    format = " MCMC [:bar] :percent Estimated completion time: :eta",
    total = totiter, clear = FALSE, width= 100)
  
  k0 = hyper[3]
  L = length(unique(colayer)) #number of layers
  Layerdim = table(colayer) #number of column for each layer
  n_tot = dim(data)[1] #number of items
  
  #the data in a matrix X, which is nxL
  X = data 
  
  #RANDOM PROB.S
  pim = array(NA, dim = c(H0, H, L)) #weights (in log scale)
  b = array(NA, dim = c(H0, H, L)) #sticks for pi
  
  pi0 = matrix(NA, nrow = H0, ncol = L) #weights of the common overall in HDP
  b0 = matrix(NA, nrow = H0, ncol = L) #sticks for pi0 (in log scale)
  theta0 = array(NA, dim = c(H0, dim(data)[2])) #atoms
  
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
    col_ends = cumsum(Layerdim)[l]
    col_init = col_ends - Layerdim[l] + 1
    m[l, ] = kmeans(as.matrix(X[,col_init:col_ends]), 2)$cluster
    m[l, ] = rep(1, n_tot)
    if(max(m[l,])>H0){print("ERROR: initialization has too many clusters, increase H0")}
  }
  c[1, ] = m[1, ]
  for(l in 2:L){
    c[l, ] = (m[l-1,]-1)*H + m[l, ]
  }
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
  
  alpha0 = hyper[1]; alpha = hyper[2]
  
  #mcmc
  for (iter in 1:totiter){
    pb$tick()
    for (l in 1:L){
      col_ends = cumsum(Layerdim)[l]
      col_init = col_ends - Layerdim[l] + 1
      ql_sum = c(rev(cumsum(rev(q[l,])))[2:H0],0)
      for (h1 in 1:H0){
        b0[h1, l] = rbeta(1, 1 + q[l,h1], alpha0 + ql_sum[h1])
        if(h1>1){
          pi0[h1, l] = log(b0[h1, l]) + sum(log(1 - b0[1:(h1-1), l]))
        }else{
          pi0[h1, l] = log(b0[h1, l])
        }
        theta0[h1, col_init:col_ends ] = posterior_sampler(X[m[l,]==h1, col_init:col_ends])
        
        
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
    
    #prior is gamma(3,3)
    alpha0 = rgamma(1, 3, 3 - sum(log(1-b0)))
    alpha = rgamma(1, 3, 3 - sum(log(1-b)))
    print(alpha0, alpha)
    
    
    for (l in 1:L){
      col_ends = cumsum(Layerdim)[l]
      col_init = col_ends - Layerdim[l] + 1
      for (i in 1:n_tot){
        if(l>1){past = m[1, i]}else{past = 1}
        fut = rep(0,H)
        if(l==1){
          for (cc in 1:H){
            for (d in 1:H0){
              fut[cc] = fut[cc] + 
                (k[l+1,(k[l,cc]-1)*H + d] == m[l+1,i])*exp(pim[k[l,cc],d,l+1]) +
                (k[l+2,(k[l,cc]-1)*H + d] == m[l+2,i])*exp(pim[k[l,cc],d,l+2])
            } 
          }
          fut = log(fut)
        }
        prob = (pim[past,,l]) + kernel_eval(X[i,col_init:col_ends ], 
                                            theta0[ ,col_init:col_ends]) + fut
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
          prob[kk] = pi0[kk,l] + kernel_eval_joint(X[c[l,]==cc,col_init:col_ends], 
                                                   theta0[kk,col_init:col_ends])
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
    m_saved[iter,,] = m
    print(rand.index(m[1,],true_layer[,1]))
  }
  
  return(m_saved)
}

