#functions from mcclust.ext 
#source: https://github.com/sarawade/mcclust.ext.git 
### Calculate Credible ball

credibleball=function(c.star,cls.draw,c.dist=c("VI","Binder"),alpha=0.05){
  
  n=length(c.star)
  c.dist <- match.arg(c.dist, choices=c.dist)
  
  #Distance functions
  dist.binder=function(c1,c2){
    f=0
    for(i in 1:n){
      f=f+sum(abs((c1==c1[i])-(c2==c2[i])))
    }
    f=f/(n^2)
    return(f)
  }
  
  dist.vi=function(c1,c2){
    f=0
    for(i in 1:n){
      ind1=(c1==c1[i])
      ind2=(c2==c2[i])
      f=f+(log2(sum(ind1))+log2(sum(ind2))-2*log2(sum(ind1*ind2)))/n
    }
    return(f)
  }
  
  #Compute distance between optimal and samples
  M=nrow(cls.draw)
  d=rep(0,M)
  if(c.dist=="Binder") d=apply(cls.draw,1,dist.binder,c2=c.star)
  if(c.dist=="VI") d=apply(cls.draw,1,dist.vi,c2=c.star)
  sd=sort(d,decreasing=F,index.return=T)
  ind.star=ceiling((1-alpha)*M)
  
  cb=cls.draw[sd$ix[1:ind.star],]
  cb.dist=sd$x[1:ind.star]
  
  # Extremes of credible ball
  c.horiz=matrix(cb[(cb.dist==cb.dist[ind.star]),],ncol=n)
  k.cb=apply(cb,1,max)
  min.ind=which(k.cb==min(k.cb))
  c.uppervert=matrix(cb[min.ind[cb.dist[min.ind]==cb.dist[min.ind[length(min.ind)]]],],ncol=n)
  max.ind=which(k.cb==max(k.cb))
  c.lowervert=matrix(cb[max.ind[cb.dist[max.ind]==cb.dist[max.ind[length(max.ind)]]],],ncol=n)
  
  output=list(c.star=c.star,c.horiz=c.horiz,c.uppervert=c.uppervert,c.lowervert=c.lowervert,dist.horiz=cb.dist[ind.star],dist.uppervert=cb.dist[min.ind[length(min.ind)]],dist.lowervert=cb.dist[max.ind[length(max.ind)]])
  class(output)="credibleball"
  return(output)
}

### Compute optimal partition that minimizes the posterior expected loss,
### by performing greedy search: at every iteration, consider the L-closest 
### ancestors and the L-closest descendents

greedy=function(psm,cls.draw=NULL,loss=NULL,start.cl=NULL,maxiter=NULL,L=NULL,suppress.comment=TRUE){
  
  if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
    stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}
  
  n=nrow(psm)
  if(is.null(loss)) loss="VI.lb"
  if(is.null(start.cl)) start.cl=1:n
  if(is.null(maxiter)) maxiter=2*n
  if(is.null(L)) L=2*n
  
  if(loss=="VI" & is.null(cls.draw)) stop("cls.draw must be provided if loss=''VI''")
  
  EVI_lb_local=function(c){
    f=0
    for(i in 1:n){
      ind=(c==c[i])
      f=f+(log2(sum(ind))+log2(sum(psm[i,]))-2*log2(sum(ind*psm[i,])))/n
    }
    return(f)
  }
  EVI_local=function(c){
    M=nrow(cls.draw)
    f=0
    for(i in 1:n){
      ind=(c==c[i])
      f=f+log2(sum(ind))
      for(m in 1:M){
        indm=(cls.draw[m,]==cls.draw[m,i])
        f=f+(log2(sum(indm))-2*log2(sum(ind*indm)))/M
      }
    }
    f=f/n
    return(f)
  }
  EBL_local=function(c){
    f=0
    for(i in 1:n){
      f=f+sum(abs((c[i]==c)-psm[i,]))
    }
    f=f/(n^2)
    return(f)
  }  
  
  #Extra functions
  
  c_combine=function(c,i,j){
    c[c==i|c==j]=min(i,j)
    c[c>max(i,j)]=c[c>max(i,j)]-1
    return(c)
  }
  
  dist_merge_ij=function(ni,nj){
    d=0
    if(loss=="VI.lb"||loss=="VI"){d=((ni+nj)/n)*log2((ni+nj)/n)-(ni/n)*log2(ni/n)-(nj/n)*log2(nj/n)}
    if(loss=="Binder"){d=((ni+nj)^2-(ni)^2-(nj)^2)/(n^2)}
    return(d)
  }
  dist_split_i=function(x,ni){
    d=0
    if(loss=="VI.lb"||loss=="VI"){d=(ni/n)*log2(ni/n)-(x/n)*log2(x/n)-((ni-x)/n)*log2((ni-x)/n)}
    if(loss=="Binder"){d=((ni)^2-(x)^2-(ni-x)^2)/(n^2)}
    return(d)
  }
  
  #Function which given a configuration, finds the L closests configurations and
  # selects the one with the smallest EBL
  local_explore=function(c_star,val_star){
    k=max(c_star)
    nj=rep(0,k)
    for(j in 1:k){
      nj[j]=sum(c_star==j)
    }
    snj_ind=list()
    unj=unique(nj)
    unj=sort(unj)
    U=length(unj)
    lnj=rep(0,U)
    for(i in 1:U){
      snj_ind[[i]]=which(nj==unj[i])
      lnj[i]=length(snj_ind[[i]])
    }
    c_opt=c_star
    val_opt=val_star
    #Merge two clusters
    #Compute distance of merge any two clusters
    if(k>1){
      m_ind=1:U
      if(lnj[1]==1){m_ind=m_ind[-1]}
      d_1=apply(matrix(unj[m_ind],length(m_ind),1),1,dist_merge_ij,nj=unj[1])
      d_mat=rbind(d_1,rep(1,length(m_ind)),m_ind)
      if((U-(lnj[U]==1))>1){
        for(i in 2:(U-(lnj[U]==1))){
          m_ind=i:U
          if(lnj[i]==1){m_ind=m_ind[-1]}
          d_i=apply(matrix(unj[m_ind],length(m_ind),1),1,dist_merge_ij,nj=unj[i])
          d_mat=cbind(d_mat,rbind(d_i,rep(i,length(m_ind)),m_ind))
        }
      }
      sd=sort(d_mat[1,],index.return=T)
      d_mat=matrix(d_mat[,sd$ix],nrow=3)
      colind=1
      l=0
      ub=min(L,choose(k,2))
      while(l<ub){
        i=d_mat[2,colind]
        h=d_mat[3,colind]
        reps=0
        if(i!=h){reps=lnj[i]*lnj[h]}
        if(i==h){reps=choose(lnj[i],2)}
        nj_indi=snj_ind[[i]]
        if(i==h){nj_indi=nj_indi[-lnj[i]]}
        for(i_ind in 1:length(nj_indi)){
          n1_ind=nj_indi[i_ind]
          if(i!=h){nj_indh=snj_ind[[h]]}
          if(i==h){nj_indh=snj_ind[[i]][(i_ind+1):lnj[i]]}
          for(h_ind in 1:length(nj_indh)){
            n2_ind=nj_indh[h_ind]
            #proposed partition
            c_p=c_combine(c_star,n1_ind,n2_ind)
            #compute loss
            if(loss=="VI.lb"){val_p=EVI_lb_local(c_p)}
            if(loss=="VI"){val_p=EVI_local(c_p)}
            if(loss=="Binder"){val_p=EBL_local(c_p)}
            if(val_p<val_opt){
              c_opt=c_p
              val_opt=val_p
            }
          }
        }
        #Update l and colind
        colind=colind+1
        l=l+reps
      }
    }
    #Spliting two clusters
    #Compute distance of splitting any clusters
    if(k<n){
      sind=1+(unj[1]==1)
      m_ind=1:floor(unj[sind]/2)
      d_1=apply(matrix(m_ind,length(m_ind),1),1,dist_split_i,ni=unj[sind])
      d_mat=rbind(d_1,rep(sind,length(m_ind)),m_ind)
      numsp=apply(matrix(m_ind,length(m_ind),1),1,choose,n=unj[sind])
      if((unj[sind]%%2)==0){numsp[length(numsp)]=numsp[length(numsp)]/2}
      numsp=sum(numsp)*lnj[sind]
      if(sind<U){
        for(i in (sind+1):U){
          m_ind=1:floor(unj[i]/2)
          d_i=apply(matrix(m_ind,length(m_ind),1),1,dist_split_i,ni=unj[i])
          d_mat=cbind(d_mat,rbind(d_i,rep(i,length(m_ind)),m_ind))
          numsp=c(numsp,apply(matrix(m_ind,length(m_ind),1),1,choose,n=unj[i]))
          if((unj[i]%%2)==0){numsp[length(numsp)]=numsp[length(numsp)]/2}
          numsp=numsp[1]+sum(numsp[-1])*lnj[i]
        }
      }
      sd=sort(d_mat[1,],index.return=T)
      d_mat=matrix(d_mat[,sd$ix],nrow=3)
      colind=1
      l=0
      ub=min(L,numsp)
      while(l<ub){
        i=d_mat[2,colind]
        nj_new=d_mat[3,colind]
        reps=choose(unj[i],nj_new)
        if(nj_new==(unj[i]/2)){reps=reps/2}
        for(j in 1:lnj[i]){
          ind_set=c(1:nj_new)
          for(h in 1:reps){
            c_p=c_star
            c_p[c_star==snj_ind[[i]][j]][ind_set]=k+1
            #Compute expected loss
            val_p=0
            if(loss=="VI.lb"){val_p=EVI_lb_local(c_p)}
            if(loss=="VI"){val_p=EVI_local(c_p)}
            if(loss=="Binder"){val_p=EBL_local(c_p)}
            if(val_p<val_opt){
              c_opt=c_p
              val_opt=val_p
            }
            if(h<reps){
              #Update set
              ind_set[nj_new]=ind_set[nj_new]+1
              if(ind_set[nj_new]>unj[i]){
                updateind=which(ind_set>=c((unj[i]-nj_new+1):unj[i]))[1]
                ind_set[c((updateind-1):nj_new)]=ind_set[(updateind-1)]+c(1:(nj_new-updateind+2))
              }
            }
          }
        }
        colind=colind+1
        l=l+reps*lnj[i]	
        
      }
    }
    return(list(c_star=c_opt,val_star=val_opt))
  }
  
  #Start at last
  c_star=start.cl
  val_star=0
  if(loss=="VI.lb") val_star=EVI_lb_local(c_star)
  if(loss=="VI") val_star=EVI_local(c_star)
  if(loss=="Binder") val_star=EBL_local(c_star)
  
  it=1
  stop_ind=F
  while((it<maxiter)&(!stop_ind)){
    opt=local_explore(c_star,val_star)
    if(opt$val_star==val_star){
      stop_ind=T
    }
    else{
      val_star=opt$val_star
      c_star=opt$c_star
      it=it+1
      if(!suppress.comment){cat(paste("Iteration=",it," k=",max(c_star), " Loss=",round(val_star,4),"\n"))}
    }
  }
  
  output=list(cl=c_star,value=val_star,iter.greedy=it)
  return(output)
}

### Finds optimal partition minimizes the posterior expected Binder's loss

minbinder=function(psm, cls.draw=NULL, method=c("avg","comp","draws","laugreen","greedy","all"), 
                   max.k=NULL, include.lg=FALSE, include.greedy=FALSE, start.cl.lg=NULL,start.cl.greedy=NULL,tol=0.001, 
                   maxiter=NULL,l=NULL, suppress.comment=TRUE){
  
  if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
    stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}   	
  
  method <- match.arg(method, choices=method)
  if(method %in% c("draws","all") & is.null(cls.draw)) stop("cls.draw must be provided if method=''draws''")
  
  if(method != "greedy"){
    res=minbinder(psm, cls.draw,method=method,max.k=max.k,include.lg=include.lg,start.cl=start.cl.lg, tol=tol)
    n=nrow(psm)
    res$value=res$value/(n^2)*2
    if(method!="all"| (method=="all"&!include.greedy)) {
      if(method=="all") res=c(res,list(method="all"))
      class(res)="c.estimate"
      return(res)
    }
  } 
  
  if(method=="all" & is.null(start.cl.greedy)) {
    ind=which.min(res$value)
    start.cl.greedy=res$cl[ind,]
  }
  res.greedy <- greedy(psm,loss="Binder",start.cl=start.cl.greedy, maxiter=maxiter,L=l,suppress.comment=suppress.comment)
  if(method=="greedy") {
    res.greedy=c(res.greedy,list(method="greedy"))
    class(res.greedy)="c.estimate"
    return(res.greedy)  
  }
  
  res$value <- c(res$value, res.greedy$value)
  res$cl <- rbind(res$cl,res.greedy$cl)
  res$cl[1,] <-res$cl[which.min(res$value),]
  res$value[1] <- min(res$value)
  res=c(res,list(method="all"))
  rownames(res$cl)[nrow(res$cl)] <- names(res$value)[length(res$value)] <- "greedy"    
  if(include.greedy) res=c(res,list(iter.greedy=res.greedy$iter.greedy))
  class(res)="c.estimate"
  return(res)     
}

### Finds optimal partition minimizes the lower bound to the Variation of 
### Information obtain from Jensen's inequality where the expectation and 
### log are reversed.


minVI=function(psm, cls.draw=NULL, method=c("avg","comp","draws","greedy","all"), 
               max.k=NULL, include.greedy=FALSE, start.cl=NULL, maxiter=NULL,l=NULL, suppress.comment=TRUE){
  
  if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
    stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}
  
  method <- match.arg(method, choices=method)
  if(method %in% c("draws","all") & is.null(cls.draw)) stop("cls.draw must be provided if method=''draws''")
  
  if(method == "avg" | method == "all"){
    if(is.null(max.k)) max.k <- ceiling(dim(psm)[1]/8)
    hclust.avg=hclust(as.dist(1-psm), method="average")
    cls.avg= t(apply(matrix(1:max.k),1,function(x) cutree(hclust.avg,k=x)))
    VI.avg= VI.lb(cls.avg,psm)
    val.avg <- min(VI.avg)
    cl.avg <- cls.avg[which.min(VI.avg),] 
    if(method== "avg")  {
      output=list(cl=cl.avg, value=val.avg, method="avg")
      class(output)="c.estimate"
      return(output)
    }
  }
  
  if(method == "comp" | method == "all"){
    if(is.null(max.k)) max.k <- ceiling(dim(psm)[1]/8)
    hclust.comp <- hclust(as.dist(1-psm), method="complete")
    cls.comp <-  t(apply(matrix(1:max.k),1,function(x) cutree(hclust.comp,k=x)))
    VI.comp <- VI.lb(cls.comp,psm)
    val.comp <- min(VI.comp)
    cl.comp <- cls.comp[which.min(VI.comp),] 
    if(method== "comp")  {
      output=list(cl=cl.comp, value=val.comp, method="comp")
      class(output)="c.estimate"
      return(output)
    }
  }
  
  if(method == "draws" | method == "all"){
    n=ncol(psm)
    EVI_lb_local=function(c){
      f=0
      for(i in 1:n){
        ind=(c==c[i])
        f=f+(log2(sum(ind))+log2(sum(psm[i,]))-2*log2(sum(ind*psm[i,])))/n
      }
      return(f)
    }
    VI.draws=apply(cls.draw,1,EVI_lb_local)
    val.draws <- min(VI.draws)
    cl.draw <- cls.draw[which.min(VI.draws),] 
    names(cl.draw) <- NULL
    if(method== "draws") {
      output=list(cl=cl.draw, value=val.draws, method="draws")
      class(output)="c.estimate"
      return(output)
    }
  }
  
  if(method == "greedy" | (method == "all" & include.greedy)){
    if(method=="all" & is.null(start.cl)) {
      ind=which.min(c(val.avg, val.comp, val.draws))
      start.cl=rbind(cl.avg,cl.comp,cl.draw)[ind,]
    }
    res.greedy <- greedy(psm,loss="VI.lb",start.cl=start.cl, maxiter=maxiter,L=l,suppress.comment=suppress.comment)
    if(method=="greedy"){
      res.greedy=c(res.greedy,list(method="greedy"))
      class(res.greedy)="c.estimate"
      return(res.greedy)
    }  
  }
  
  vals <- c(val.avg, val.comp, val.draws)
  cls <- rbind(cl.avg,cl.comp,cl.draw)
  if(include.greedy){
    vals <- c(vals,res.greedy$value)
    cls <- rbind(cls,res.greedy$cl)
  }
  cls <- rbind(cls[which.min(vals),], cls)
  vals <- c(min(vals), vals)
  if(include.greedy){ rownames(cls) <- names(vals) <- c("best","avg","comp","draws","greedy")
  } else rownames(cls) <- names(vals) <- c("best","avg","comp","draws")
  colnames(cls) <- NULL    
  res <- list(cl=cls, value=vals, method="all")
  
  if(include.greedy) res=c(res,list(iter.greedy=res.greedy$iter.greedy))
  class(res)="c.estimate"
  return(res)    
}

plot.credibleball=function(x,data=NULL,dx=NULL,xgrid=NULL,dxgrid=NULL,...){
  if(is.null(data)){
    stop("data must be supplied")}
  if(!is.data.frame(data)){
    stop("data must be a data.frame")}
  p=ncol(data)
  n=nrow(data)
  k.u=apply(x$c.uppervert,1,max)
  n.u=nrow(x$c.uppervert)
  k.l=apply(x$c.lowervert,1,max)
  n.l=nrow(x$c.lowervert)
  k.h=apply(x$c.horiz,1,max)
  n.h=nrow(x$c.horiz)
  par(ask = TRUE)
  if(p==1){
    x1=data[,1]
    if(is.null(dxgrid)){
      if(is.null(xgrid)){
        dout=density(x1)
      }
      else{
        dout=density(data,n=length(xgrid),from=xgrid[1],to=xgrid[length(xgrid)])
      }
      xgrid=dout$x
      dxgrid=dout$y
    }
    #Compute density estimate at x1
    x2=dx
    if(is.null(x2)){
      aout=approx(xgrid,dxgrid,xout=x1)
      x2=aout$y}
    
    #Upper
    cl=rainbow(k.u)
    par(mfrow=c(1,n.u),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
    for(i in 1:n.u){	
      plot(xgrid,dxgrid,type="l",xlab=names(data),ylab="density",main="Credible Ball: Upper Vertical Bound",...)
      for(j in 1:k.u){
        points(x1[x$c.uppervert[i,]==j], x2[x$c.uppervert[i,]==j],col=cl[j],...)
      }}
    #Lower
    cl=rainbow(k.l)
    par(mfrow=c(1,n.l),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
    for(i in 1:n.l){	
      plot(xgrid,dxgrid,type="l",xlab=names(data),ylab="density",main="Credible Ball: Lower Vertical Bound",...)
      for(j in 1:k.l){
        points(x1[x$c.lowervert[i,]==j], x2[x$c.lowervert[i,]==j],col=cl[j],...)
      }}
    #Horiziontal
    par(mfrow=c(1,n.h),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
    for(i in 1:n.h){
      cl=rainbow(k.h[i])	
      plot(xgrid,dxgrid,type="l",xlab=names(data),ylab="density",main="Credible Ball: Horizontal Bound",...)
      for(j in 1:k.h[i]){
        points(x1[x$c.horiz[i,]==j], x2[x$c.horiz[i,]==j],col=cl[j],...)
      }}
  }
  if(p==2){
    x1=data[,1]
    x2=data[,2]
    #Upper
    cl=rainbow(k.u)
    par(mfrow=c(1,n.u),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
    for(i in 1:n.u){
      plot(x1,x2,xlab=names(data)[1],ylab=names(data)[2],main="Credible Ball: Upper Vertical Bound",...)
      for(j in 1:k.u){
        points(x1[x$c.uppervert[i,]==j], x2[x$c.uppervert[i,]==j],col=cl[j],...)
      }}
    #Lower
    cl=rainbow(k.l)
    par(mfrow=c(1,n.l),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
    for(i in 1:n.l){	
      plot(x1,x2,xlab=names(data)[1],ylab=names(data)[2],main="Credible Ball: Lower Vertical Bound",...)
      for(j in 1:k.l){
        points(x1[x$c.lowervert[i,]==j], x2[x$c.lowervert[i,]==j],col=cl[j],...)
      }}
    #Horiziontal
    par(mfrow=c(1,n.h),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
    for(i in 1:n.h){
      cl=rainbow(k.h[i])	
      plot(x1,x2,xlab=names(data)[1],ylab=names(data)[2],main="Credible Ball: Horizontal Bound",...)
      for(j in 1:k.h[i]){
        points(x1[x$c.horiz[i,]==j], x2[x$c.horiz[i,]==j],col=cl[j],...)
      }}
  }
  if(p>2){
    x.pca=princomp(data,scores=T)
    x1=x.pca$scores[,1]
    x2=x.pca$scores[,2]
    #Upper
    cl=rainbow(k.u)
    par(mfrow=c(1,n.u),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
    for(i in 1:n.u){
      plot(x1,x2,xlab="PC 1",ylab="PC 2",main="Credible Ball: Upper Vertical Bound",...)
      for(j in 1:k.u){
        points(x1[x$c.uppervert[i,]==j], x2[x$c.uppervert[i,]==j],col=cl[j],...)
      }}
    #Lower
    cl=rainbow(k.l)
    par(mfrow=c(1,n.l),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
    for(i in 1:n.l){	
      plot(x1,x2,xlab="PC 1",ylab="PC 2",main="Credible Ball: Lower Vertical Bound",...)
      for(j in 1:k.l){
        points(x1[x$c.lowervert[i,]==j], x2[x$c.lowervert[i,]==j],col=cl[j],...)
      }}
    #Horiziontal
    par(mfrow=c(1,n.h),mar=c(2.5,2.5,1.5,.5)+.1, mgp=c(1.5, .5, 0))
    for(i in 1:n.h){
      cl=rainbow(k.h[i])	
      plot(x1,x2,xlab="PC 1",ylab="PC 2",main="Credible Ball: Horizontal Bound",...)
      for(j in 1:k.h[i]){
        points(x1[x$c.horiz[i,]==j], x2[x$c.horiz[i,]==j],col=cl[j],...)
      }}
  }	
}

plot.c.estimate=function(x,data=NULL,dx=NULL,xgrid=NULL,dxgrid=NULL,...){
  if(is.null(data)){
    stop("data must be supplied")}
  if(!is.data.frame(data)){
    stop("data must be a data.frame")}
  p=ncol(data)
  n=nrow(data)
  if(!is.matrix(x$cl)) x$cl=matrix(x$cl,nrow=1)
  k=apply(x$cl,1,max)
  n.c=nrow(x$cl)
  par(ask = TRUE)
  method.names=x$method
  if(method.names=="all"){method.names=names(x$value)}
  if(p==1){
    x1=data[,1]
    if(is.null(dxgrid)){
      if(is.null(xgrid)){
        dout=density(x1)
      }
      else{
        dout=density(data,n=length(xgrid),from=xgrid[1],to=xgrid[length(xgrid)])
      }
      xgrid=dout$x
      dxgrid=dout$y
    }
    #Compute density estimate at x1
    x2=dx
    if(is.null(x2)){
      aout=approx(xgrid,dxgrid,xout=x1)
      x2=aout$y}
    for(i in 1:n.c){
      cl=rainbow(k[i])
      plot(xgrid,dxgrid,type="l",xlab=names(data),ylab="density",main=paste("Method:",method.names[i]),...)
      for(j in 1:k[i]){
        points(x1[x$cl[i,]==j], x2[x$cl[i,]==j],col=cl[j],...)
      }}
    
  }
  if(p==2){
    x1=data[,1]
    x2=data[,2]
    for(i in 1:n.c){
      cl=rainbow(k[i])
      plot(x1,x2,xlab=names(data)[1],ylab=names(data)[2],main=paste("Method:",method.names[i]),...)
      for(j in 1:k[i]){
        points(x1[x$cl[i,]==j], x2[x$cl[i,]==j],col=cl[j],...)
      }}
  }
  if(p>2){
    x.pca=princomp(data,scores=T)
    x1=x.pca$scores[,1]
    x2=x.pca$scores[,2]
    for(i in 1:n.c){
      cl=rainbow(k[i])
      plot(x1,x2,xlab="PC 1",ylab="PC 2",main=paste("Method:",method.names[i]),...)
      for(j in 1:k[i]){
        points(x1[x$cl[i,]==j], x2[x$cl[i,]==j],col=cl[j],...)
      }}
  }	
}

plotpsm=function(psm,method="complete",...){
  
  if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
    stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}
  
  #sort by heirarchical clustering
  hc=hclust(as.dist(1-psm), method = method, members = NULL)
  psm_hc=psm
  n=nrow(psm)
  psm_hc[1:n,]=psm_hc[hc$order,]
  psm_hc[,1:n]=psm_hc[,hc$order]
  
  image(1:n,1:n,1-psm_hc,col=heat.colors(20),...)
}

summary.credibleball=function(object,...){
  cb=object
  k.u=apply(cb$c.uppervert,1,max)
  n.u=nrow(cb$c.uppervert)
  k.l=apply(cb$c.lowervert,1,max)
  n.l=nrow(cb$c.lowervert)
  k.h=apply(cb$c.horiz,1,max)
  n.h=nrow(cb$c.horiz)
  
  t.u=list(n.u)
  for(i in 1:n.u){
    t=table(cb$c.star,cb$c.uppervert[i,],dnn=c("Estimate","Upper vertical bound"))
    t.u[[i]]=addmargins(t)
  }
  t.l=list(n.l)
  for(i in 1:n.l){
    t=table(cb$c.star,cb$c.lowervert[i,],dnn=c("Estimate","Lower vertical bound"))
    t.l[[i]]=addmargins(t)
  }
  t.h=list(n.h)
  for(i in 1:n.h){
    t=table(cb$c.star,cb$c.horiz[i,],dnn=c("Estimate","Horizontal bound"))
    t.h[[i]]=addmargins(t)
  }
  
  output=list(n.u=n.u,k.u=k.u,t.u=t.u,n.l=n.l,k.l=k.l,t.l=t.l,n.h=n.h,k.h=k.h,t.h=t.h,dist.uppervert=cb$dist.uppervert,dist.lowervert=cb$dist.lowervert,dist.horiz=cb$dist.horiz)
  class(output)="summary.credibleball"
  return(output)
}

print.summary.credibleball=function(x,...){
  cat("The credible ball characterizes the uncertainty in the clustering esitmate")
  cat(".\n")
  cat("It can be summarized with:")
  cat("\n")
  cat("1. upper vertical bound: partitions in the ball with the fewest clusters that are most distant,")
  cat("\n")
  cat("2. lower vertical bound: partitions in the ball with the most clusters that are most distant,")
  cat("\n")
  cat("3. horizontal bound: partitions in the ball with the greatest distance.")
  cat("\n")	
  d.u=round(x$dist.uppervert,2)
  d.l=round(x$dist.lowervert,2)
  d.h=round(x$dist.horiz,2)
  if(x$n.u==1){
    if(x$k.u!=1) {
      cat(paste("The upper vertical bound has",x$k.u,"clusters with a distance of", d.u))
      cat(".\n")
    }
    else {
      cat(paste("The upper vertical bound has 1 cluster with a distance of", d.u))
      cat(".\n")
    }
  }
  else{
    cat(paste("The", x$n.u, "upper vertical bounds have", x$k.u[1], "clusters with a distance of", d.u))
    cat(".\n")
  }
  
  if(x$n.l==1){
    cat(paste("The lower vertical bound has", x$k.l, "clusters with a distance of", d.l))
    cat(".\n")
  }
  else{
    cat(paste("The", x$n.l, "lower vertical bounds have", x$k.l[1], "clusters with a distance of", d.l))
    cat(".\n")
  }
  
  if(x$n.h==1){
    cat(paste("The horizontal bound has", x$k.h, "clusters with a distance of", d.h))
    cat(".\n")
  }
  else{
    cat(paste("The", x$n.h, "horizontal bounds have "))
    for(i in 1:(x$n.h-1)){
      cat(x$k.h[i])
      if(i<(x$n.h-1)){ cat(", ")}
    }
    cat(paste(" and",x$k.h[x$n.h], "clusters with a distance of", d.h))
    cat(".\n")
  }
  
  cat("\nA cross tabulation with the point estimate of the partition and the bounds.\n")
  cat("\n")
  for (i in 1:x$n.u){
    print(x$t.u[[i]])
  }
  cat("\n")
  for (i in 1:x$n.l){
    print(x$t.l[[i]])
  }
  cat("\n")
  for (i in 1:x$n.h){
    print(x$t.h[[i]])
  }
  
}

summary.c.estimate=function(object,...){
  x=object
  if(!is.matrix(x$cl)) x$cl=matrix(x$cl,nrow=1)
  k=apply(x$cl,1,max)
  n.c=nrow(x$cl)
  
  t=list(n.c)
  onevec=rep("size",ncol(x$cl))
  for(i in 1:n.c){
    t[[i]]=table(onevec,x$cl[i,],dnn=c("","cluster"))
  }
  
  output=list(method=x$method,k=k,n.c=n.c,t=t,value=x$value)
  class(output)="summary.c.estimate"
  return (output)
}

print.summary.c.estimate=function(x,...){
  
  if(x$method!="all"){
    cat("The partition estimate found with the",x$method, "method has a posterior expected loss of\n", round(x$value,2),"and contains",x$k,"clusters of sizes:\n")
    print(x$t[[1]])
  }
  if(x$method=="all"){
    cat("The best partition estimate has a posterior expected loss of\n", round(x$value[1],2)," and contains",x$k[1],"clusters of sizes:\n")
    print(x$t[[1]])
    cat("\n")
    method.names=names(x$value)
    for(i in 2:x$n.c){
      cat("The partition estimate found with the",method.names[i], "method has a posterior expected loss of\n", round(x$value[i],2)," and contains",x$k[i],"clusters of sizes:\n")
      print(x$t[[i]])
      cat("\n")
    }
  }
}

## Computes the lower bound to the posterior expected Variation of Information

VI.lb=function(cls,psm){
  
  if(any(psm !=t(psm)) | any(psm >1) | any(psm < 0) | sum(diag(psm)) != nrow(psm) ){
    stop("psm must be a symmetric matrix with entries between 0 and 1 and 1's on the diagonals")}
  
  if(is.vector(cls)) cls <- t(cls)
  
  n=dim(psm)[1]
  
  VI.lb.compute=function(c){
    f=0
    for(i in 1:n){
      ind=(c==c[i])
      f=f+(log2(sum(ind))+log2(sum(psm[i,]))-2*log2(sum(ind*psm[i,])))/n
    }
    return(f)
  }
  output=apply(cls,1,VI.lb.compute)
  return(output)
}

## Computes the posterior expected Variation of Information

VI=function(cls,cls.draw){
  
  if(is.vector(cls)) cls <- t(cls)
  
  n=dim(cls.draw)[2]
  M=dim(cls.draw)[1]
  
  VI.compute=function(c){
    f=0
    for(i in 1:n){
      ind=(c==c[i])
      f=f+log2(sum(ind))
      for(m in 1:M){
        indm=(cls.draw[m,]==cls.draw[m,i])
        f=f+(log2(sum(indm))-2*log2(sum(ind*indm)))/M
      }
    }
    f=f/n
    return(f)
  }
  output=apply(cls,1,VI.compute)
  return(output)
}


