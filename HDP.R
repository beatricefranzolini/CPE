library(progress)#to draw the progress bar

#compute the marginal likelihood of a norm-norm model (marg lik of a cluster) univariate obs
margnn <- function(y, m = 0, s2 = 1, sigma2 = 1, log = TRUE){
  sigma = sqrt(sigma2); s = sqrt(s2)
  n = length(y); Y2 = sum(y**2); Y = sum(y)
  p = - n / 2 * log(2 * pi) - (n - 1) * log(sigma) - 1/2 * log(n * s2 + sigma2) - 
    Y2 / (2 * sigma2) - m**2 / (2 * s2) + (s * Y / sigma + sigma * m / s ) ** 2 / (2 * (n * s2 + sigma2))
  return(p)
} 

# #compute the marginal likelihood of a norm-norm model (marg lik of a cluster) multivariate obs indep. kernel indep. base
# margnn_multi <- function(y, m = 0, s2 = 1, sigma2 = 1, log = TRUE){
#   J = dim(y)[2]; n = dim(y)[1]
#   ybarn = colSums(y)
#   y2barn = sum(y**2)
#   Sigmanew =  solve(n + solve(diag(J)*s2)) 
#   
#   p = - n / 2 * log(2 * pi) - 1/2 * log(2 * pi * J * s2) +
#     1 / 2 * log(2 * pi * det(Sigmanew) ) +
#     1 / 2 * (t(ybarn) %*% Sigmanew %*% ybarn - y2barn)
#   return(p)
# } 

margnn_multi <- function(y, m = 0, s2 = 1, sigma2 = 1, log = TRUE) {
  n <- nrow(y)
  J <- ncol(y)
  
  ybar <- colMeans(y)
  diff_y <- sweep(y, 2, ybar)
  ssq <- sum(diff_y^2)
  
  # quadratic term for posterior mean
  diff_mu <- ybar - m
  quad_mu <- sum(diff_mu^2)
  
  # log marginal likelihood
  logp <- -n * J / 2 * log(2 * pi * sigma2) +
    - J / 2 * log(s2 / (s2 + sigma2 / n)) +
    - 0.5 * ssq / sigma2 +
    - 0.5 * n * quad_mu / (sigma2 + n * s2)
  
  if (log) return(logp) else return(exp(logp))
}


#MCMC for the (constant) HDP estimating the posterior distribution 
#with dim(data)[2]-variate ind. normal kernel and normal ind base measure
#to get the HDP of Camerlenghi 2018 with no groups, set group = rep(1,n)
#to get the HDP of Teh 2010 set group to the known groups.
HDP_normal_normal <- function(data, group, hyper=c(0.1, 0.1, 0, 0.1, 1), c = FALSE, totiter = 1000){
  pb <- progress_bar$new(
    format = " MCMC [:bar] :percent Estimated completion time: :eta",
    total = totiter, clear = FALSE, width= 100)
  
  #data = scale(data)
  alpha = hyper[1]
  alpha0 = hyper[2]
  m = hyper[3]
  s2 = hyper[4]
  sigma2 = hyper[5]
  n = length(group)
  #initialization
  if(!c[1]){
    c = kmeans(data, centers = 2, nstart = 20)$cluster
  }#c are tables not clusters 
  k = c #each table serves is own dish 
  k_saved = NULL
  for (iter in (1:totiter)){
    pb$tick()
    #sample t###################################################################
    for (i in (1:n)){
      j = group[i]
      restjc = c[group == j & seq(1,n)!=i]
      y_minus_i = data[group == j & seq(1,n)!=i,,drop=FALSE]
      nc = table(restjc)
      c_unique = unique(restjc)
      prob = NULL; h = 1
      for (cc in c_unique){
        prob[h] = margnn_multi(rbind(y_minus_i[t(restjc == cc),,drop=FALSE], data[i,]), m, s2, sigma2, log = TRUE ) - 
          margnn_multi(y_minus_i[t(restjc == cc),,drop = FALSE], m, s2, sigma2, log = TRUE ) +
          log(nc[as.character(cc)])
        h = h + 1
      }
      restjk = k[-i]
      data_minus_i = data[-i,,drop=FALSE]
      nk = rowSums(table(restjk, c[-i]) != 0)
      k_unique = unique(restjk)
      prob_inner = NULL; l = 1
      for (kk in k_unique){
        prob_inner[l] = margnn_multi(rbind(data_minus_i[t(restjk == kk),,drop=FALSE], data[i,]), m, s2, sigma2, log = TRUE ) - 
          margnn_multi(data_minus_i[t(restjk == kk),,drop = FALSE], m, s2, sigma2, log = TRUE ) +
          log(nk[as.character(kk)])
        l = l + 1
      }
      prob_inner[l] = margnn_multi(data[i,,drop = FALSE], m, s2, sigma2, log = TRUE) + log(alpha0)
      prob[h] = log( sum( exp(prob_inner) ) / ( sum(nk) + alpha0 )) + log(alpha)
      prob = prob - max(prob)
      prob = exp(prob) / sum( exp(prob) )
      c[i] = sample(c(c_unique, setdiff(1:n,c[-i])[1]), 1, prob = prob)
      if (sum(c[-i]==c[i])>0){
        k[i] = restjk[c[-i]==c[i]][1]
      }else{
        prob_inner = prob_inner - max(prob_inner)
        prob_inner = exp(prob_inner) / sum( exp(prob_inner) )
        k[i] = sample(c(k_unique, setdiff(1:n,k_unique)[1]), 1, prob = prob_inner)
      }
    }
    #sample k#################################################################
    for (tt in unique(c)){
      restjk = k[c!=tt]
      data_minus_tt = data[c!=tt,,drop=FALSE]
      nk = rowSums(table(restjk, c[c!=tt]) != 0)
      k_unique = unique(restjk)
      prob_inner = NULL; l = 1
      for (kk in k_unique){
        prob_inner[l] = margnn_multi(rbind(data_minus_tt[t(restjk == kk),,drop=FALSE], data[c==tt,,drop=FALSE]), m, s2, sigma2, log = TRUE ) - 
          margnn_multi(data_minus_tt[t(restjk == kk),,drop = FALSE], m, s2, sigma2, log = TRUE ) +
          log(nk[as.character(kk)])
        l = l + 1
      }
      prob_inner[l] = margnn_multi(data[c==tt,,drop = FALSE], m, s2, sigma2, log = TRUE) + log(alpha0)
      prob_inner = prob_inner - max(prob_inner)
      prob_inner = exp(prob_inner) / sum( exp(prob_inner) )
      k[c==tt] = sample(c(k_unique, setdiff(1:n,k_unique)[1]), 1, prob = prob_inner)
    }
    k_saved = rbind(k_saved,k)
    #print(k)
    #print(table(k))
  }
  return(k_saved)
}