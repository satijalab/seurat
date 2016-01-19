
project_map=function(z,x_old,sum_X_old,x_old_tsne,P_tol=5e-6,perplexity=30) {
  sum_z=sum(z^2)
  #x_old=test2@pca.rot[,1:5]
  #sum_X_old=rowSums((x_old^2))
  D_org=sum_z+(-2*as.matrix(x_old)%*%t(z)+sum_X_old)
  P=d2p_cell(D_org,perplexity)

  nn_points = which(P> P_tol);    #Use only a small subset of points to comupute embedding. This keeps all the points that are proximal to the new point
  X_nn_set = x_old[nn_points,]  #Original points                                       
  y_nn_set = x_old_tsne[nn_points,];  #Computed embeddings                     
  P_nn_set = P[nn_points,]; #Probabilities
  y_new0 = (t(as.matrix(y_nn_set))%*%t(as.matrix(rbind(P_nn_set,P_nn_set))))[,1] #Initial guess of point as a weighted average
  sink("/dev/null")
  y_new=optim(y_new0,KullbackFun,gr=NULL,y_nn_set,P_nn_set,method = "Nelder-Mead")
  sink()
  #plot(test2@tsne.rot)
  #points(y_new$par[1],y_new$par[2],col="red",cex=1,pch=16)
  #points(test2@tsne.rot[cell.num,1],test2@tsne.rot[cell.num,2],col="blue",cex=1,pch=16)
  #return(dist(as.matrix(rbind(y_new$par,test2@tsne.rot[cell.num,]))))
  return(y_new$par)
}

d2p_cell=function(D,u=15,tol=1e-4) {
  betamin=-Inf; betamax=Inf;
  tries = 0; tol=1e-4; beta=1
  beta.list=Hbeta(D,beta); h=beta.list[[1]]; thisP=beta.list[[2]]; flagP=beta.list[[3]]
  hdiff = h - log(u);
  while (abs(hdiff) > tol && tries < 50) {
    if (hdiff > 0) {
      betamin = beta;
      if (betamax==Inf) {
        beta = beta * 2;
      } else
      { 
        beta = (beta + betamax) / 2;
      }
    } else {
      betamax = beta;
      if (betamin==-Inf) {
        beta = beta / 2;
      } else {
        beta = (beta + betamin) / 2;
      }
    }
    beta.list=Hbeta(D,beta); h=beta.list[[1]]; thisP=beta.list[[2]]; flagP=beta.list[[3]]
    hdiff = h - log(u);
    tries = tries + 1
  }
  # set the final row of p
  P=thisP
  #Check if there are at least 10 points that are highly similar to the projected point
  return(P)
}


KullbackFun=function(z,y,P) {

  #Computes the Kullback-Leibler divergence cost function for the embedding x in a tSNE map
  #%P = params{1};              %Transition probabilities in the original space. Nx1 vector
  #%y = params{2};              %tSNE embeddings of the training set. Nx2 vector

  print(z); print(dim(y))
  Cost0 = sum(P * log(P));          #Constant part of the cost function
  #Compute pairwise distances in embedding space
  sum_z=sum(z^2)
  sum_y=rowSums((y^2))
  D_yz=sum_z+(-2*as.matrix(y)%*%t(matrix(z,nrow=1))+sum_y)
  Q = 1/(1+D_yz);
  Q = Q/sum(Q)               #Transition probabilities in the embedded space
  Cost = Cost0 - sum(P*log(Q)); #% - 100 * sum(Q .* log(Q));
  return(Cost)
}


Hbeta=function(D, beta) {
    flagP=1;
    P = exp(-D * beta);
    sumP = sum(P);
    if (sumP<1e-8) { #In this case it means that no point is proximal.
      P = ones(length(P),1) / length(P);
      sumP = sum(P);
      flagP=0;
    }
    H = log(sumP) + beta * sum(D * P) / sumP;
    P = P / sumP;
    return(list(H, P, flagP))
}