#t-HDP _ conditional algorithm - NNIW model - multivariate
library(progress)#to draw the progress bar

#triangular dependence for 3 layers, multivariate layers
#normal kernel and normal-inverse-chisquared base (correlation within cluster = 0)
#H0 and H are the number of mixture components used in the approximation 
#of the infinite mixtures (It works only for H0 = H)
telescopic_HDP_NNIX_multi_3L1P <- function(data, colayer, hyper=c(1, 1, 1),
                                               H0 = 10, H = 10, totiter = 1000){
  
  
  #has to work for many X and one param
  #dmvnorm(x, mean = rep(0, p), sigma = diag(p), log = FALSE, checkSymmetry = TRUE)
  kernel_eval_joint = function(x, mean, sigma){
    if(length(x)==0){return(0)}else{
      temp = 0 
      for (i in 1:dim(x)[1]){
        temp = temp + dnorm(x[i,], mean, sqrt(diag(sigma)), log = TRUE)
      }
      return(temp)}
  }
  
  #has to work for one X and many param
  kernel_eval = function(x, mean, sigma){
    temp = NULL; p = length(x)
    for (jj in (1:dim(mean)[1])){
      temp = c(temp, sum(dnorm(x, mean[jj,], sqrt(diag(sigma[jj,,])), log = TRUE)))
    }
    return(temp)
  }
  
  
  posterior_sampler_old = function(x, mu0_prior = NULL, k0 = 1, Lambda0 = NULL, v0 = NULL){
    if(!is.matrix(x)){x = matrix(x, ncol=length(x))}
    n_temp = dim(x)[1]; p = dim(x)[2] 
    if(is.null(mu0_prior)){mu0_prior = rep(0,p)}
    if(is.null(Lambda0)){Lambda0 = 0.5*diag(p)}
    if(is.null(v0)){v0 = p}
    vn = v0 + n_temp; kn = k0 + n_temp
    if(n_temp == 1){
      S = matrix(0, nrow = p, ncol = p)
      xbar = colMeans(x)
    }else if(n_temp == 0){
      S = matrix(0, nrow = p, ncol = p)
      xbar = mu0_prior
    }else{
      S = cov(x) * (n_temp - 1)
      xbar = colMeans(x)
    }
    Lambdan = Lambda0 + S + 
      k0*n_temp / kn * (xbar - mu0_prior)%*%t((xbar - mu0_prior))
    mun = (k0 * mu0_prior + n_temp * xbar) / kn
    #Sigma_post = diag(diag(rinvwishart(vn, (Lambdan))))
    Sigma_post = diag(diag(MCMCpack::riwish(vn, (Lambdan))))
    mu_post = rnorm(p, mun, sqrt(diag(kn^(-1)*Sigma_post)) )
    return(list("mu" = mu_post, 
                "Sigma" = Sigma_post))
  }
  
  posterior_sampler = function(x, mu0_prior = NULL, k0 = 1, Lambda0 = NULL, v0 = NULL){
    if(!is.matrix(x)){x = matrix(x, ncol=length(x))}
    n_temp = dim(x)[1]; p = dim(x)[2] 
    if(is.null(mu0_prior)){mu0_prior = rep(0,p)}
    if(is.null(Lambda0)){Lambda0 = 0.5*diag(p)}
    if(is.null(v0)){v0 = p}
    vn = v0 + n_temp; kn = k0 + n_temp
    if(n_temp == 1){
      S = matrix(0, nrow = p, ncol = p)
      xbar = colMeans(x)
    }else if(n_temp == 0){
      S = matrix(0, nrow = p, ncol = p)
      xbar = mu0_prior
    }else{
      S = cov(x) * (n_temp - 1)
      xbar = colMeans(x)
    }
    Lambdan = diag( (Lambda0*v0 + S + 
                       k0*n_temp / kn * (xbar - mu0_prior)%*%t((xbar - mu0_prior))))
    mun = (k0 * mu0_prior + n_temp * xbar) / kn
    Sigma_post = diag(diag(rinvchisq(p, vn))*Lambdan)
    mu_post = rnorm(p, mun, sqrt(diag(kn^(-1)*Sigma_post)) )
    return(list("mu" = mu_post, 
                "Sigma" = Sigma_post))
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
  mu0 = array(NA, dim = c(H0, dim(data)[2])) #atoms
  Sigma0 = array(NA, dim = c(H0, dim(data)[2], dim(data)[2]))
  
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
    m[l, ] = kmeans(as.matrix(X[,col_init:col_ends]),2)$cluster
    #m[l, ] = rep(1, n_tot)
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
  
  alpha0 = rep(hyper[1],3); alpha = rep(hyper[1],3)
  temp12 = NULL; temp13 = NULL
  adj_temp12 = NULL; adj_temp13 = NULL
  #mcmc
  for (iter in 1:totiter){
    pb$tick()
    for (l in 1:L){
      col_ends = cumsum(Layerdim)[l]
      col_init = col_ends - Layerdim[l] + 1
      ql_sum = c(rev(cumsum(rev(q[l,])))[2:H0],0)
      for (h1 in 1:H0){
        b0[h1, l] = rbeta(1, 1 + q[l,h1], alpha0[l] + ql_sum[h1])
        if(h1>1){
          pi0[h1, l] = log(b0[h1, l]) + sum(log(1 - b0[1:(h1-1), l]))
        }else{
          pi0[h1, l] = log(b0[h1, l])
        }
        temp = posterior_sampler(X[m[l,]==h1, col_init:col_ends])
        mu0[h1, col_init:col_ends ] = temp$mu
        Sigma0[h1, col_init:col_ends, col_init:col_ends] = temp$Sigma
        #print(Sigma0)
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
      alpha[l] = rgamma(1, 1 + H*H0, rate = 1 - sum(log(1-b[,,l])))
      if(alpha[l]<0.1){alpha[l]=0.1}
      alpha0[l] = rgamma(1, 1 + H0, rate = 1 - sum(log(1-b0[,l])))
      if(alpha0[l]<0.1){alpha0[l]=0.1}
      
    }
    
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
                                            mu0[ ,col_init:col_ends], 
                                            Sigma0[,col_init:col_ends, col_init:col_ends]) + fut
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
          prob[kk] = pi0[kk,l] + kernel_eval_joint(X[c[l,]==cc,col_init:col_ends, drop =FALSE], 
                                                   mu0[kk,col_init:col_ends],
                                                   Sigma0[kk,col_init:col_ends, col_init:col_ends])
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
    for (layer in (1:L)){
      write.table(t(m[layer,]), file = paste("layer",layer,".csv",sep=""), append = TRUE, row.names = FALSE,
                  col.names= FALSE)
    }
    write.table(t(c(alpha0, alpha)), file = paste("concentrations.csv"), append = TRUE, row.names = FALSE,
                col.names= FALSE)
    #print(rand.index(m[1,],true_layer[,1]))
    temp12 = c(temp12, rand.index(m[1,],m[2,]))
    temp13 = c(temp13, rand.index(m[1,],m[3,]))
    adj_temp12 = c(adj_temp12, adj.rand.index(m[1,],m[2,]))
    adj_temp13 = c(adj_temp13, adj.rand.index(m[1,],m[3,]))
    print(c( mean(adj_temp12), mean(adj_temp13)))
    print(table(m[1,]))
    print(table(m[2,]))
    print(table(m[3,]))
  }
  write.table(t(c(temp12, temp13)), file = paste("rand_indexes_chains.csv"), row.names = FALSE,
              col.names= FALSE)
  return(m_saved)
}