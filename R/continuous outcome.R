
#'@title sGIM for the continuous outcome
#'@description This function can integrate summary information from external source
#'@usage information.borrowing_lgim_normal(beta_initial,gamma_initial,sigma2_initial,theta_ext_ini,theta_ext_est,v_ext,y_int,x_int,used_index)
#'@param beta_initial the initial values of the main parameter vector
#'@param gamma_initial the initial values of the shift parameter
#'@param sigma2_initial the initial value of the variance parameter
#'@param theta_ext_ini the initial values of parameters in the reduced model
#'@param theta_ext_est the summary information (parameter estimates) from the external source
#'@param v_ext the variance-covariance matrix from the external source
#'@param y_int the outcome vector from the internal data
#'@param x_int the covariate matrix from the internal data
#'@param used_index index for parameters in the reduced model in relation to the full model
#'@return A list of parameter estimates
#'@export
#'@import mvtnorm rootSolve MASS Matrix sandwich

information.borrowing_lgim_normal<-function(beta_initial,gamma_initial,sigma2_initial,
                                     theta_ext_ini,theta_ext_est,v_ext,
                                     y_int,x_int,used_index)
{
  n<-nrow(x_int)

  #estimating function phi
  ee_ib<-function(x_int,theta_ext,used_index,beta,gamma,sigma2)
  {
    p_i<-x_int%*%beta
    p_red_i<-x_int[,used_index]%*%theta_ext
    #estimating function
    ee<-x_int[,used_index]*c(gamma*sigma2+x_int%*%(beta)-p_red_i)*c(exp(x_int%*%(beta)*gamma)) #some extra terms of X,Z left, should add!!
    phi_ee<-t(ee)
    #derivative of estimating function
    dphi_ee<-lapply(1:n,function(x)
    {
      dphi_ee1<-x_int[x,used_index]%*%t(x_int[x,])*c(exp(x_int[x,]%*%(beta)*gamma))+
        x_int[x,used_index]%*%t(x_int[x,])*c(gamma*sigma2+x_int[x,]%*%(beta)-p_red_i[x])*c(exp(x_int[x,]%*%(beta)*gamma))*c(gamma) #derivative with respect to beta #some extra terms of X,Z left, should add!!
      dphi_ee2<-x_int[x,used_index]*sigma2*c(exp(x_int[x,]%*%(beta)*gamma))+
        x_int[x,used_index]*c(gamma*sigma2+x_int[x,]%*%(beta)-p_red_i[x])*c(exp(x_int[x,]%*%(beta)*gamma))*c(x_int[x,]%*%(beta))#derivative with respect to gamma1 #some extra terms of X,Z left, should add!!
      dphi_ee22<-x_int[x,used_index]*gamma*c(exp(x_int[x,]%*%(beta)*gamma))#derivative with respect to sigma2 #some extra terms of X,Z left, should add!!
      dphi_ee3<--x_int[x,used_index]%*%t(x_int[x,used_index])*c(exp(x_int[x,]%*%(beta)*gamma)) #derivative with respect to theta #some extra terms of X,Z left, should add!!

      cbind(dphi_ee1,dphi_ee2,dphi_ee22,dphi_ee3)
    }
    )
    dphi_ee<-do.call(cbind,dphi_ee)
    return(cbind(phi_ee,dphi_ee))
  }

  #lambda_find
  lambda_find<-function(beta,gamma,theta_ext,sigma2)
  {
    ZZ<-ee_ib(x_int,theta_ext,used_index,beta,gamma,sigma2)[,1:(n)]
    dim(ZZ)
    apply(ZZ,1,mean)

    gamma_e<-1
    c<-0
    lambda<-rep(0,nrow(ZZ))
    tol<-10e-7
    Delta_old<-0

    repeat{
      rl<-R1der(lambda,ZZ)
      rll<-R2der(lambda,ZZ)
      Delta<--ginv(rll)%*%rl
      if(mean(abs(Delta))<tol | mean(Delta-Delta_old)==0 | c>100)
      {break}else{
        repeat{
          mm<-0
          repeat{
            delta<-gamma_e*Delta
            index_1<-apply(ZZ,2,function (xx)
            {ifelse(1+t(lambda+delta)%*%as.matrix(xx,ncol=1)<=1/n,1,0)}
            )
            if (sum(index_1)>0)
            {gamma_e<-gamma_e/2
            mm<-mm+1}else{break}}
          index_2<-ifelse(R0der(lambda+delta,ZZ)-R0der(lambda,ZZ)<0,1,0)
          if (index_2==1)
          {gamma_e<-gamma_e/2}else{break}
        }
        Delta_old<-Delta
      }
      lambda<-lambda+delta
      c<-c+1
      gamma_e<-(c)^(-0.5)
    }
    lambda
  }

  ##first derivative of -log EL
  R1der<-function(lambda,ZZ)
  {
    apply(ZZ,2,function(xx)
    {as.matrix(xx,ncol=1)/as.vector((1+t(lambda)%*%as.matrix(xx,ncol=1)))})%*%rep(1,ncol(ZZ))
  }

  #second derivative of -log EL
  R2der<-function(lambda,ZZ)
  {
    r2der<-0
    for(i in 1:ncol(ZZ))
    {
      r2der_i<--as.matrix(ZZ[,i],ncol=1)%*%t(as.matrix(ZZ[,i],ncol=1))/as.vector(1+t(lambda)%*%as.matrix(ZZ[,i],ncol=1))^2
      r2der<-r2der+r2der_i
    }
    r2der
  }

  #-log EL
  R0der<-function(lambda,ZZ)
  {
    apply(ZZ,2, function (xx) {log(as.vector(1+t(lambda)%*%as.matrix(xx,ncol=1)))})%*%rep(1,ncol(ZZ))
  }

  #estimate beta finally
  beta_ee<-function(beta,gamma,sigma2,theta_ext,theta_ext_est,v_ext,used_index,x_int)
  {

    #parameter_all<-c(beta,gamma)

    lambda<-lambda_find(beta,gamma,theta_ext,sigma2)


    M_lh<-function(beta,gamma,sigma2,theta_ext){
      m_sum<-0
      lambda<-lambda_find(beta,gamma,theta_ext,sigma2)
      total<-ee_ib(x_int=x_int,theta_ext=theta_ext,
                   used_index=used_index,beta=beta,gamma=gamma,sigma2=sigma2)
      ZZ<-total[,1:(n)]

      #try
      L<-0
      f<-0
      for (i in 1:n)
      {
        L_i<--log(1+t(matrix(lambda,ncol=1))%*%ZZ[,i])
        f_i<--0.5*log(sigma2)- (y_int[i]-x_int[i,]%*%beta)^2/(2*sigma2) #done
        m_sum<-m_sum+L_i+f_i

        #try
        L<-L+L_i
        f<-f+f_i
      }
      m_sum-t(theta_ext-theta_ext_est)%*%ginv(v_ext)%*%(theta_ext-theta_ext_est)/2
    }

    l<-1
    repeat
    {
      total<-ee_ib(x_int=x_int,theta_ext=theta_ext,
                   used_index=used_index,beta=beta,gamma=gamma,sigma2=sigma2)

      #total<-g_fct(beta=beta,beta_external=beta_external,y_main=y_main,x_main=x_main, cali.matrix=cali.matrix, ind.ex=ind.ex,cal.score=cal.score)

      #total<-g_fct(beta=beta,beta_external=beta_external,y_main=y_main,x_main=x_main)
      #total<-wgeef(beta=beta,adata,r=r,id=id,dist=dist,time=time,n=n)
      ZZ<-total[,1:(n)]
      ZZ_d<-total[,(n+1):((ncol(x_int)+2+length(used_index)+1)*(n-1)+ncol(x_int)+2+length(used_index)+1)]
      tau<-1
      s_ee_sum<-0
      ds_ee_sum<-0
      for1_sum<-0
      for2_sum<-0
      for3_sum<-0
      for (i in 1:(n))
      {
        scaler<-(1/(1+t(matrix(lambda,ncol=1))%*%ZZ[,i]))

        s_ee_i<-matrix(x_int[i,],ncol=1)%*%(y_int[i]-x_int[i,]%*%beta)/(sigma2) #done
        s_ee_sigma2_i<-(y_int[i]-x_int[i,]%*%beta)^2/(2*(sigma2^2))-1/(2*sigma2)
        s_ee_i<-c(s_ee_i,0,s_ee_sigma2_i,rep(0,length(used_index)))

        ds_ee1_i<--matrix(x_int[i,],ncol=1)%*%matrix(x_int[i,],nrow=1)/sigma2 #done
        ds_ee1_sigma2_i<--(y_int[i]-x_int[i,]%*%beta)^2/(sigma2^3)+1/(2*(sigma2^2))
        ds_ee1_betasigma2_i<--matrix(x_int[i,],ncol=1)%*%(y_int[i]-x_int[i,]%*%beta)/(sigma2^2)
        ds_ee1_i<-cbind(ds_ee1_i,0,c(ds_ee1_betasigma2_i),matrix(0,nrow=nrow(ds_ee1_i),
                                                                 ncol=length(used_index)))
        ds_ee_i<-rbind(ds_ee1_i,
                       matrix(0,nrow=length(used_index)+1+1,ncol=ncol(ds_ee1_i)))
        ds_ee_i[-c(1:ncol(x_int)),c(1:ncol(x_int))]<-t(ds_ee_i[c(1:ncol(x_int)),-c(1:ncol(x_int))]) #change here, there is take derivative of sigma with respect to beta
        ds_ee_i[-c(1:ncol(x_int)),-c(1:ncol(x_int))][2,2]<-ds_ee1_sigma2_i

        for1_i<-matrix(t(ZZ_d[,((i-1)*(ncol(x_int)+2+length(used_index))+1):
                                (i*(ncol(x_int)+2+length(used_index)))]),
                       nrow=ncol(x_int)+2+length(used_index))*as.vector(scaler)
        for2_i<-ZZ[,i]%*%t(ZZ[,i])*as.vector(scaler^2)
        for3_i<-t(for1_i)

        s_ee_sum<-s_ee_sum+s_ee_i
        ds_ee_sum<-ds_ee_sum+ds_ee_i
        for1_sum<-for1_sum+for1_i
        for2_sum<-for2_sum+for2_i
        for3_sum<-for3_sum+for3_i

      }
      M_b<-s_ee_sum-for1_sum%*%matrix(lambda,ncol=1)-c(rep(0,ncol(x_int)+2),ginv(v_ext)%*%(theta_ext-theta_ext_est)) #done
      adjust_2ed<-diag(0,ncol(x_int)+2+length(used_index))
      adjust_2ed[-c(1:(ncol(x_int)+2)),-c(1:(ncol(x_int)+2))]<-ginv(v_ext)
      m_bb<-ds_ee_sum-for1_sum%*%ginv(for2_sum)%*%for3_sum-adjust_2ed #it has been checked, should be negative sign!  #recheck the last term, positive or negative
      delta<-ginv(m_bb)%*%M_b
      beta_update<-beta-tau*delta[1:ncol(x_int)]
      theta_ext_update<-theta_ext-tau*delta[-c(1:(ncol(x_int)+2))]

      ##keep gamma small
      gamma_update<-gamma-tau*delta[(ncol(x_int)+1)]
      sigma2_update<-sigma2-tau*delta[(ncol(x_int)+2)]
      ll<-1
      repeat{
        if(mean(abs(gamma_update))<1 & sigma2_update>0)
        {break}else{
          ll<-ll+1
          gamma_update<-gamma-tau/ll*delta[(ncol(x_int)+1)]
          sigma2_update<-sigma2-tau/ll*delta[(ncol(x_int)+2)]
          beta_update<-beta-tau/ll*delta[1:ncol(x_int)]
          theta_ext_update<-theta_ext-tau/ll*delta[-c(1:(ncol(x_int)+2))]

        }
      }



      jj<-1
      repeat{
        criterion_update<-M_lh(beta_update,gamma_update,sigma2_update,theta_ext_update)
        criterion_old<-M_lh(beta,gamma,sigma2,theta_ext)

        if(!is.na(criterion_update) & !is.na(criterion_old) & criterion_update >criterion_old | jj>6)
        {
          break
        }else{
          tau<-tau/2
          beta_update<-beta-tau*delta[1:ncol(x_int)]
          ##keep gamma small
          gamma_update<-gamma-tau/ll*delta[(ncol(x_int)+1)]
          sigma2_update<-sigma2-tau/ll*delta[(ncol(x_int)+2)]
          theta_ext_update<-theta_ext-tau*delta[-c(1:(ncol(x_int)+2))]

          ll<-1
          repeat{
            if(mean(abs(gamma_update))<1 & sigma2_update>0){break}else{
              ll<-ll+1
              gamma_update<-gamma-tau/ll*delta[(ncol(x_int)+1)]
              sigma2_update<-sigma2-tau/ll*delta[(ncol(x_int)+2)]
              beta_update<-beta-tau/ll*delta[1:ncol(x_int)]
              theta_ext_update<-theta_ext-tau/ll*delta[-c(1:(ncol(x_int)+2))]

            }
          }
          jj<-jj+1
        }
      }

      if(mean(abs(beta_update-beta))<1e-3 & abs(sigma2_update-sigma2)<1e-4 | l>30)
      {break}else{
        beta<-beta_update
        gamma<-gamma_update
        sigma2<-sigma2_update
        theta_ext<-theta_ext_update
        lambda<-lambda_find(beta,gamma,theta_ext,sigma2)
        l<-l+1
      }

    }

    total<-ee_ib(x_int=x_int,theta_ext=theta_ext,
                 used_index=used_index,beta=beta,gamma=gamma,sigma2=sigma2)
    ZZ<-total[,1:(n)]
    Pi<-1/(1+t(matrix(lambda,ncol=1))%*%ZZ)/n

    return(list(beta_final=beta_update,gamma_final=gamma_update,sigma2_final=sigma2_update,theta_final=theta_ext_update,l=l,lambda=lambda,Pi=Pi))
  }

  final_results<-beta_ee(beta_initial,gamma_initial,sigma2_initial,theta_ext_ini,
                         theta_ext_est,v_ext,used_index,x_int)
  lambda_final<-final_results$lambda
  beta_est<-final_results$beta_final
  gamma_est<-final_results$gamma_final
  sigma2_est<-final_results$sigma2_final
  theta_est<-final_results$theta_final
  Pi<-final_results$Pi
  list(beta_est=beta_est,
       gamma_est=gamma_est,
       sigma2_est=sigma2_est,
       theta_est=theta_est,
       lambda=lambda_final,
       Pi=Pi)
}







#'@title Asymptotic variance based on sGIM for the continuous outcome
#'@description This function calculate asymptotic variance of sGIM
#'@usage asymptotic_var_general_normal(beta_est,gamma_est,sigma2_est,theta_est,v_ext,y_int,x_int,used_index)
#'@param beta_est the estimates of the main parameter vector
#'@param gamma_est the estimates of the shift parameter
#'@param sigma2_est the initial value of the variance parameter
#'@param theta_est the estimates of parameters in the reduced model
#'@param v_ext the variance-covariance matrix from the external source
#'@param y_int the outcome vector from the internal data
#'@param x_int the covariate matrix from the internal data
#'@param used_index index for parameters in the reduced model in relation to the full model
#'@return The asymptotic variance-covariance matrix
#'@export
#'@import mvtnorm rootSolve MASS Matrix sandwich
asymptotic_var_general_normal<-function(beta_est=beta_est,gamma_est=gamma_est,
                                 sigma2_est=sigma2_est,
                                 theta_est=theta_est,v_ext=v_ext,
                                 y_int=y_int,x_int=x_int,used_index=used_index)
{
  n<-nrow(x_int)
  ds_ee1_mean<-0
  Jtt_mean<-0
  Jtu_mean<-0
  for(x in 1:n)
  {
    # s_ee_i<-matrix(x_int[x,],ncol=1)%*%(y_int[x]-(1+exp(-x_int[x,]%*%beta_est))^(-1))
    # s_ee_i<-c(s_ee_i,0,rep(0,length(used_index)))

    p_i<-x_int[x,]%*%beta_est
    p_red_i<-x_int[x,used_index]%*%theta_est
    #estimating function
    ee<-x_int[x,used_index]*c(gamma_est*sigma2_est+x_int[x,]%*%(beta_est)-p_red_i)*c(exp(x_int[x,]%*%(beta_est)*gamma_est))
    phi_ee_i<-(ee)

    ###-E(gg)
    Egg_i<--phi_ee_i%*%t(phi_ee_i)/n
    Jtt_mean<-Jtt_mean+Egg_i

    #derivative of estimating function
    # dphi_ee1<-x_int[x,used_index]%*%t(x_int[x,]) #derivative with respect to beta
    # dphi_ee2<-x_int[x,used_index]*sigma2_est #derivative with respect to gamma1
    # dphi_ee22<-x_int[x,used_index]*gamma_est #derivative with respect to sigma2
    # dphi_ee3<--x_int[x,used_index]%*%t(x_int[x,used_index]) #derivative with respect to theta


    dphi_ee1<-x_int[x,used_index]%*%t(x_int[x,])*c(exp(x_int[x,]%*%(beta_est)*gamma_est))+
      x_int[x,used_index]%*%t(x_int[x,])*c(gamma_est*sigma2_est+x_int[x,]%*%(beta_est)-p_red_i)*c(exp(x_int[x,]%*%(beta_est)*gamma_est))*c(gamma_est)#derivative with respect to beta #some extra terms of X,Z left, should add!!
    dphi_ee2<-x_int[x,used_index]*sigma2_est*c(exp(x_int[x,]%*%(beta_est)*gamma_est))+
      x_int[x,used_index]*c(gamma_est*sigma2_est+x_int[x,]%*%(beta_est)-p_red_i)*c(exp(x_int[x,]%*%(beta_est)*gamma_est))*c(x_int[x,]%*%(beta_est))#derivative with respect to gamma1 #some extra terms of X,Z left, should add!!
    dphi_ee22<-x_int[x,used_index]*gamma_est*c(exp(x_int[x,]%*%(beta_est)*gamma_est))#derivative with respect to sigma2 #some extra terms of X,Z left, should add!!
    dphi_ee3<--x_int[x,used_index]%*%t(x_int[x,used_index])*c(exp(x_int[x,]%*%(beta_est)*gamma_est)) #derivative with respect to theta #some extra terms of X,Z left, should add!!


    Jtu_i<-cbind(dphi_ee1,dphi_ee2,dphi_ee22,dphi_ee3)/n

    # dphi_ee1<-x_int[x,used_index]%*%t(x_int[x,]) #with respect to beta
    # dphi_ee2<-x_int[x,used_index]*c(1-p_red_i)*exp(gamma_est)*c(p_i) #with respect to gamma
    # dphi_ee3<--x_int[x,used_index]%*%t(x_int[x,used_index]) #with respective theta

    # Jtu_i<-cbind(dphi_ee1,dphi_ee2,dphi_ee3)/n

    Jtu_mean<-Jtu_mean+Jtu_i

    ###E(ss)
    # ds_ee1_i<-matrix(x_int[x,],ncol=1)%*%matrix(x_int[x,],nrow=1)/n
    # ds_ee1_mean<-ds_ee1_mean+ds_ee1_i

    ds_ee1_i<-matrix(x_int[x,],ncol=1)%*%matrix(x_int[x,],nrow=1)/sigma2_est #done
    ds_ee1_sigma2_i<-(y_int[x]-x_int[x,]%*%beta_est)^2/(sigma2_est^3)-1/(2*(sigma2_est^2))
    ds_ee1_betasigma2_i<-matrix(x_int[x,],ncol=1)%*%(y_int[x]-x_int[x,]%*%beta_est)/(sigma2_est^2)
    ds_ee1_i<-cbind(ds_ee1_i,0,c(ds_ee1_betasigma2_i))
    ds_ee_i<-rbind(ds_ee1_i,
                   matrix(0,nrow=1+1,ncol=ncol(ds_ee1_i)))
    ds_ee_i[-c(1:ncol(x_int)),c(1:ncol(x_int))]<-t(ds_ee_i[c(1:ncol(x_int)),-c(1:ncol(x_int))]) #change here, there is take derivative of sigma with respect to beta
    ds_ee_i[-c(1:ncol(x_int)),-c(1:ncol(x_int))][2,2]<-ds_ee1_sigma2_i
    ds_ee1_mean<-ds_ee1_mean+ds_ee_i/n
  }

  ##calculate variance
  ###Juu<-ds_ee1_mean
  Juu<-matrix(0,nrow=ncol(x_int)+2+nrow(v_ext),
              ncol=ncol(x_int)+2+nrow(v_ext))
  Juu[1:(ncol(x_int)+2),1:(ncol(x_int)+2)]<-ds_ee1_mean
  Juu[-c(1:(ncol(x_int)+2)),-c(1:(ncol(x_int)+2))]<-1/n*ginv(v_ext)

  J_V<-matrix(0,nrow=ncol(x_int)+2+2*nrow(v_ext),
              ncol=ncol(x_int)+2+2*nrow(v_ext))
  J_V[1:nrow(v_ext),]<-cbind(Jtt_mean,Jtu_mean)
  J_V[,1:nrow(v_ext)]<-t(cbind(Jtt_mean,Jtu_mean))
  J_V[-(1:nrow(v_ext)),-(1:nrow(v_ext))]<-Juu

  I_V<-bdiag(-Jtt_mean,Juu) #the last block in I_V should be changed if only partial V is available

  V_asymptotic<-ginv(J_V)%*%I_V%*%ginv(J_V)/n

  return(V_asymptotic=V_asymptotic)
}







#'@title GIM for the continuous outcome
#'@description This function can integrate summary information from external source
#'@usage information.borrowing_gim_normal(beta_initial,theta_ext_ini,theta_ext_est,v_ext,y_int,x_int,used_index)
#'@param beta_initial the initial values of the main parameter vector
#'@param theta_ext_ini the initial values of parameters in the reduced model
#'@param theta_ext_est the summary information (parameter estimates) from the external source
#'@param v_ext the variance-covariance matrix from the external source
#'@param y_int the outcome vector from the internal data
#'@param x_int the covariate matrix from the internal data
#'@param used_index index for parameters in the reduced model in relation to the full model
#'@return A list of parameter estimates
#'@export
#'@import mvtnorm rootSolve MASS Matrix sandwich
information.borrowing_gim_normal<-function(beta_initial,
                                    theta_ext_ini,theta_ext_est,v_ext,
                                    y_int,x_int,used_index)
{
  n<-nrow(x_int)

  #estimating function phi
  ee_ib<-function(x_int,theta_ext,used_index,beta)
  {
    p_i<-x_int%*%beta
    p_red_i<-x_int[,used_index]%*%theta_ext
    #estimating function
    ee<-x_int[,used_index]*c(p_i-p_red_i)
    phi_ee<-t(ee)
    #derivative of estimating function
    dphi_ee<-lapply(1:n,function(x)
    {
      dphi_ee1<-x_int[x,used_index]%*%t(x_int[x,]) #derivative with respect to beta
      #dphi_ee2<-x_int[x,used_index]*c(1-p_red_i[x])*exp(gamma)*c(p_i[x]) #derivative with respect to gamma
      dphi_ee3<--x_int[x,used_index]%*%t(x_int[x,used_index]) #derivative with respect to theta

      cbind(dphi_ee1,dphi_ee3)
    }
    )
    dphi_ee<-do.call(cbind,dphi_ee)
    return(cbind(phi_ee,dphi_ee))
  }

  #lambda_find
  lambda_find<-function(beta,theta_ext)
  {
    ZZ<-ee_ib(x_int,theta_ext,used_index,beta)[,1:(n)]
    dim(ZZ)
    apply(ZZ,1,mean)

    gamma_e<-1
    c<-0
    lambda<-rep(0,nrow(ZZ))
    tol<-10e-7
    Delta_old<-0

    repeat{
      rl<-R1der(lambda,ZZ)
      rll<-R2der(lambda,ZZ)
      Delta<--ginv(rll)%*%rl
      if(mean(abs(Delta))<tol | mean(Delta-Delta_old)==0 | c>100)
      {break}else{
        repeat{
          mm<-0
          repeat{
            delta<-gamma_e*Delta
            index_1<-apply(ZZ,2,function (xx)
            {ifelse(1+t(lambda+delta)%*%as.matrix(xx,ncol=1)<=1/n,1,0)}
            )
            if (sum(index_1)>0)
            {gamma_e<-gamma_e/2
            mm<-mm+1}else{break}}
          index_2<-ifelse(R0der(lambda+delta,ZZ)-R0der(lambda,ZZ)<0,1,0)
          if (index_2==1)
          {gamma_e<-gamma_e/2}else{break}
        }
        Delta_old<-Delta
      }
      lambda<-lambda+delta
      c<-c+1
      gamma_e<-(c)^(-0.5)
    }
    lambda
  }

  ##first derivative of -log EL
  R1der<-function(lambda,ZZ)
  {
    apply(ZZ,2,function(xx)
    {as.matrix(xx,ncol=1)/as.vector((1+t(lambda)%*%as.matrix(xx,ncol=1)))})%*%rep(1,ncol(ZZ))
  }

  #second derivative of -log EL
  R2der<-function(lambda,ZZ)
  {
    r2der<-0
    for(i in 1:ncol(ZZ))
    {
      r2der_i<--as.matrix(ZZ[,i],ncol=1)%*%t(as.matrix(ZZ[,i],ncol=1))/as.vector(1+t(lambda)%*%as.matrix(ZZ[,i],ncol=1))^2
      r2der<-r2der+r2der_i
    }
    r2der
  }

  #-log EL
  R0der<-function(lambda,ZZ)
  {
    apply(ZZ,2, function (xx) {log(as.vector(1+t(lambda)%*%as.matrix(xx,ncol=1)))})%*%rep(1,ncol(ZZ))
  }

  #estimate beta finally
  beta_ee<-function(beta,theta_ext,theta_ext_est,v_ext,used_index,x_int)
  {

    #parameter_all<-c(beta,gamma)

    lambda<-lambda_find(beta,theta_ext)


    M_lh<-function(beta,theta_ext){
      m_sum<-0
      lambda<-lambda_find(beta,theta_ext)
      total<-ee_ib(x_int=x_int,theta_ext=theta_ext,
                   used_index=used_index,beta=beta)
      ZZ<-total[,1:(n)]

      #try
      L<-0
      f<-0
      for (i in 1:n)
      {
        L_i<--log(1+t(matrix(lambda,ncol=1))%*%ZZ[,i])
        f_i<--0.5*log(mean((y_int-x_int%*%beta)^2))
        m_sum<-m_sum+L_i+f_i

        #try
        L<-L+L_i
        f<-f+f_i
      }
      m_sum-t(theta_ext-theta_ext_est)%*%ginv(v_ext)%*%(theta_ext-theta_ext_est)/2
    }



    l<-1
    repeat
    {
      total<-ee_ib(x_int=x_int,theta_ext=theta_ext,
                   used_index=used_index,beta=beta)

      #total<-g_fct(beta=beta,beta_external=beta_external,y_main=y_main,x_main=x_main, cali.matrix=cali.matrix, ind.ex=ind.ex,cal.score=cal.score)

      #total<-g_fct(beta=beta,beta_external=beta_external,y_main=y_main,x_main=x_main)
      #total<-wgeef(beta=beta,adata,r=r,id=id,dist=dist,time=time,n=n)
      ZZ<-total[,1:(n)]
      ZZ_d<-total[,(n+1):((ncol(x_int)+length(used_index)+1)*(n-1)+ncol(x_int)+length(used_index)+1)]
      tau<-1
      s_ee_sum<-0
      ds_ee_sum<-0
      for1_sum<-0
      for2_sum<-0
      for3_sum<-0
      for (i in 1:(n))
      {
        scaler<-(1/(1+t(matrix(lambda,ncol=1))%*%ZZ[,i]))

        s_ee_i<-matrix(x_int[i,],ncol=1)%*%(y_int[i]-x_int[i,]%*%beta)
        s_ee_i<-c(s_ee_i,rep(0,length(used_index)))

        ds_ee1_i<--matrix(x_int[i,],ncol=1)%*%matrix(x_int[i,],nrow=1)
        ds_ee1_i<-cbind(ds_ee1_i,matrix(0,nrow=nrow(ds_ee1_i),
                                        ncol=length(used_index)))
        ds_ee_i<-rbind(ds_ee1_i,
                       matrix(0,nrow=length(used_index),ncol=ncol(ds_ee1_i)))

        for1_i<-matrix(t(ZZ_d[,((i-1)*(ncol(x_int)+length(used_index))+1):
                                (i*(ncol(x_int)+length(used_index)))]),
                       nrow=ncol(x_int)+length(used_index))*as.vector(scaler)
        for2_i<-ZZ[,i]%*%t(ZZ[,i])*as.vector(scaler^2)
        for3_i<-t(for1_i)

        s_ee_sum<-s_ee_sum+s_ee_i
        ds_ee_sum<-ds_ee_sum+ds_ee_i
        for1_sum<-for1_sum+for1_i
        for2_sum<-for2_sum+for2_i
        for3_sum<-for3_sum+for3_i

      }
      M_b<-s_ee_sum-for1_sum%*%matrix(lambda,ncol=1)-c(rep(0,ncol(x_int)),ginv(v_ext)%*%(theta_ext-theta_ext_est))
      adjust_2ed<-diag(0,ncol(x_int)+length(used_index))
      adjust_2ed[-c(1:(ncol(x_int))),-c(1:(ncol(x_int)))]<-ginv(v_ext)
      m_bb<-ds_ee_sum-for1_sum%*%ginv(for2_sum)%*%for3_sum-adjust_2ed #it has been checked, should be negative sign!  #recheck the last term, positive or negative
      delta<-ginv(m_bb)%*%M_b
      beta_update<-beta-tau*delta[1:ncol(x_int)]
      theta_ext_update<-theta_ext-tau*delta[-c(1:(ncol(x_int)))]

      ##keep gamma small
      # gamma_update<-gamma-tau*delta[ncol(x_int)+1]
      # ll<-1
      # repeat{
      #   if(abs(gamma_update)<2)
      #   {break}else{
      #     ll<-ll+1
      #     gamma_update<-gamma-tau/ll*delta[ncol(x_int)+1]
      #     beta_update<-beta-tau/ll*delta[1:ncol(x_int)]
      #     theta_ext_update<-theta_ext-tau/ll*delta[-c(1:(ncol(x_int)+1))]
      #
      #   }
      # }



      jj<-1
      repeat{
        criterion_update<-M_lh(beta_update,theta_ext_update)
        criterion_old<-M_lh(beta,theta_ext)

        if(!is.na(criterion_update) & !is.na(criterion_old) & criterion_update >criterion_old | jj>6)
        {
          break
        }else{
          tau<-tau/2
          beta_update<-beta-tau*delta[1:ncol(x_int)]
          ##keep gamma small
          #gamma_update<-gamma-tau*delta[ncol(x_int)+1]
          theta_ext_update<-theta_ext-tau*delta[-c(1:(ncol(x_int)))]

          # ll<-1
          # repeat{
          #   if(abs(gamma_update)<2){break}else{
          #     ll<-ll+1
          #     gamma_update<-gamma-tau/ll*delta[ncol(x_int)+1]
          #     beta_update<-beta-tau/ll*delta[1:ncol(x_int)]
          #     theta_ext_update<-theta_ext-tau/ll*delta[-c(1:(ncol(x_int)+1))]
          #
          #   }
          # }
          jj<-jj+1
        }
      }

      if(mean(abs(beta_update-beta))<1e-3 | l>30)
      {break}else{
        beta<-beta_update
        #gamma<-gamma_update
        theta_ext<-theta_ext_update
        lambda<-lambda_find(beta,theta_ext)
        l<-l+1
      }

    }
    return(list(beta_final=beta_update,theta_final=theta_ext_update,l=l,lambda=lambda))
  }

  final_results<-beta_ee(beta_initial,theta_ext_ini,
                         theta_ext_est,v_ext,used_index,x_int)
  lambda_final<-final_results$lambda
  beta_est<-final_results$beta_final
  #gamma_est<-final_results$gamma_final
  theta_est<-final_results$theta_final
  list(beta_est=beta_est,
       #gamma_est=gamma_est,
       theta_est=theta_est,
       lambda=lambda_final)
}
