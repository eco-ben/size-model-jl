#===========================
#Run model
#===========================

run_model<-function(params, initial.run = FALSE, fish.herbs = TRUE){
  
  with(params, {
    
    
    #---------------------------------------------------------------------------------------
    # Functions called in integration
    #---------------------------------------------------------------------------------------
    
    ui0 = 10^pp
    Fref=((min.fishing.size-xmin)/dx)+1
    
    # Function to construct the spectrum at the start of integration
    u.init.f = function(x, ui0, init.slope) { return (ui0*10^(init.slope*x)) }
    
    # Function to build a lookup table for diet preference of all combinations of predator 
    # and prey body size: diet preference (in the predator spectrum only)
    
    phi.f=function(q){
      q=q
      q0=q0
      sd.q=sd.q
      qmax=qmax
      qmin=qmin
      
      phi=ifelse(q>0,exp(-(q-q0)*(q-q0)/(2*sd.q*sd.q))/(sd.q*sqrt(2.0*pi)),0)   # normalise feeding kernel to sum to 1
      #phi=ifelse(q<=qmax,phi,0)                                                # remove truncation
      #phi=ifelse(q>=qmin,phi,0)
      return(phi)
    }
    
    # Functions to build lookup tables for components of integration which remain constant
    gphi.f = function(q)        { return( 10^(-q)       * phi.f(q )) }  #growth
    mphi.f = function(q)        { return( 10^(alpha*q2) * phi.f(q2)) }	#mortality
    
    # Function to build lookup table for components of 10^(alpha*x)
    expax.f    = function(x, alpha.u) { return( 10^(alpha*x )) }
    
    # Function to compute convolution products: can be for growth or mortality depending on phi
    convolution.f = function(phi, u)
    {   res = matrix(NA,length(x),1)
        for (i in 1:length(res))
        {   res[i] = 0.0
            res[i] = sum(phi[,i] * u * dx)
        }
        return(res)
    }
    
      
    #---------------------
    #Availability function
    #---------------------
  
  
   vulnerability <- function(refuge_density, competitor_density){
     res = 1
     
     for(i in 1:length(x)){
       
       if (competitor_density[i] == 0) {res[i] = 1}
       
       else {res[i] = 1 - (refuge_density[i] / competitor_density[i])}
       
       if (res[i] < min.A) {res[i] = min.A}
       else {res[i] = res[i]}
     }
     return(res)
   }
        
    #--------------------------------
    # Growth and Mortality Equations
    #--------------------------------
    
    ## predators
    death.u<-function(u,A,expax,mphi,a.u,pref.pel){return((a.u*pref.pel*A*expax)*(u*dx)%*%(mphi))}    # use faster matrix method instead of convolution loop function
     
    ## herbivores
    death.h<-function(u,A,expax,mphi,a.h,pref.herb){return((a.h*pref.herb*A*expax)*(u*dx)%*%(mphi))}
    
    ## detritivores
    death.v<-function(u,A,expax,mphi,a.v,pref.ben){return((a.v*pref.ben*A*expax)*(u*dx)%*%(mphi))}

    
    ## detritus output (g.m-3.yr-1)
    out.w<-function(A.v,A.h,v,w,h,alpha.h,pref.det){return(sum((A.v*(10^(x*0.75))*w*v)*dx)+sum((A.h*pref.det*(10^(x*alpha.h))*w*h)*dx))} # biomass density removed by detritivores per year    
    ##Added a pref.det and pref.alg function so that herbivores are both feeding on both types of resource
    
#----------------------------------
#Algal grazing dynamics

    ## algae output (g.m-2.yr-1)
    out.a <- function(A.h, alpha.h, h, a, pref.alg){return(sum(pref.alg*A.h*(10^(x*alpha.h)*a*h)*dx))} #biomass density removed by herbivores ?per day if using daily data or per hyear if using yearly data? 
            
    #---------------------------------------------------------------------------------------
    # Initialising matrices
    #---------------------------------------------------------------------------------------
    
      
    # q1 is a square matrix holding the log(predatorsize/preysize) for all combinations of sizes
    q1  = matrix(NA, length(x), length(y))
    for (i in 1:length(y)) { q1[,i] = y[i] - x}
    
    # q2 is the reverse matrix holding the log(preysize/predatorsize) for all combinations of sizes
    q2 = matrix(-q1, length(x), length(y))	
    
    # matrix for recording the two size spectra over time
    V = U = H   = array(0, c(length(x), N))
    
    # vector to hold detrtitus biomass density (g.m-3)
    W = rep(0,N)
    
    #Vector to hold algal biomass density (g.m-3)
    A = rep(0,N)
    
    # matrix for keeping track of growth and ingested food:
    GG.v = GG.u = GG.h = array(0, c(length(x), N))  

    # matrix for keeping track of vulnerability through time
    a.u = a.h = array(0, c(length(x), N))
   
    
##*************
    # Optional matrices for keeping track of reproduction from ingested food
    #R.v = R.u = R.h = array(0, c(length(x), N))
    
    # matrix for keeping track of predation mortality
    PM.v = PM.u = PM.h  = array(0, c(length(x), N))   
    
    # matrix for keeping track of  total mortality (Z)
    Z.v = Z.u = Z.h = array(0, c(length(x), N))
    
    # matrix for keeping track of senescence mortality and other (intrinsic) mortality
    SM.v = SM.u  = SM.h  =OM.v = OM.u  = OM.h = array(0, length(x))
    
    # empty vector to hold fishing mortality rates at each size class
    Fvec_pred = rep(0,length(x))
    Fvec_herb = rep(0,length(x))
    
    #Empty vector to hold fisheries yield
    Yp  = array(0, c(length(x), N))  
    Yh = array(0, c(length(x),N))
    
    # short hand for matrix indexing
    idx=2:end
    
    # lookup tables for terms in the integrals which remain constant over time
    gphi  = gphi.f(q1)
    mphi  = mphi.f(q2)
    
    # lookup table for components of 10^(alpha*x)
    expax = expax.f(x, alpha)
    
    # vectors for storing slopes of regression
    slope.v.lm = slope.u.lm = rep(0,(N-1))
    slope.v.nls = slope.u.nls = rep(0,(N-1))
    
    #---------------------------------------------------------------------------------------
    # Numerical integration
    #---------------------------------------------------------------------------------------
    
    ##INITIAL VALUES
    
    if (initial.run==T) {
      
      # if running to equilibrium use given initial values
      
      U[1:(ref-1),1]<-u.init.f(x,ui0,r.plank)[1:(ref-1)]*flow    # (phyto+zoo)plankton size spectrum  
      
      U[ref:end,1]<-u.init.f(x,ui0,init.slope = pred.slope)[ref:end] 
      
      H[ref.herb:end,1]<-u.init.f(x,ui0,init.slope = herb.slope)[ref.herb:end]*herb.prod
       
      V[refinv:(ref.det-1),1]<-u.init.f(x,invert.prod,invert.slope)[refinv:(ref.det-1)]
      V[ref.det:inv_end,1]<-u.init.f(x,invert.prod,init.slope = invert.slope)[ref.det:inv_end]
      
      W[1]<-W.init                                                          # initial detritus biomass density (g.m-3) 
      A[1]<-A.init                                                          # initial algal biomass density (g.m-3)
      
    } else {
      
      # use values from non-complex initial run
      
      U[1:(ref-1),1]<-u.init.f(x,ui0,r.plank)[1:(ref-1)]*flow                    # (phyto+zoo)plankton size spectrum  
      U[ref:end,1]<-initial.res$Preds[ref:end]                                      # initial consumer size spectrum
      H[ref.herb:end,1]<-initial.res$Herbs[ref.herb:end]
      V[ref.det:inv_end,1]<-initial.res$Invs[ref.det:inv_end]   # initial detritivore spectrum
      W[1]<-initial.res$W
      A[1]<-initial.res$A
      
    }
    
    # intrinsic natural mortality
    OM.u<-mu0*10^(-0.25*x)      
    OM.v<-mu0*10^(-0.25*x)      
    OM.h<-mu0*10^(-0.25*x)
    
    # senescence mortality rate to limit large fish from building up in the system
    # same function as in Law et al 2008, with chosen parameters gives similar M2 values as in Hall et al. 2006
    SM.u=k.sm*10^(p.s*(x-xs))
    SM.v=k.sm*10^(p.s*(x-xs))
    SM.h=k.sm*10^(p.s*(x-xs))
    
    
    #Fishing mortality at each size
    Fvec_pred[Fref:end] = rep(Fmort_pred,length(Fvec_pred[Fref:end]))  
    Fvec_herb[Fref:end] = rep(Fmort_herb,length(Fvec_herb[Fref:end]))
    
    # iteration over time
    for (i in 1:(N-1)) {
      
      
##----------------------------------------
##Density dependent vulnerability function
a.u[,i] <- vulnerability(refuge, U[,i]+H[,i])      
a.h[,i] <- vulnerability(refuge, H[,i]+U[,i])
a.v <- 1

#------------------------
##Feeding rate functions
      
      #Feeding on predators
      f.pel <- as.vector((a.u[,i]*pref.pel*A.u*expax)*(U[,i]*dx)%*%(gphi))
 
      #Feeding on herbivores
      f.herb <- as.vector((a.h[,i]*pref.herb*A.u*expax)*(H[,i]*dx)%*%(gphi))      
      
      #Feeding on invertebrates
      f.ben <- as.vector((pref.ben*A.u*expax)*(V[,i]*dx)%*%(gphi))
      
      #Feeding on algae
      f.alg<-((1/10^x)*(pref.alg*A.h*10^(x*alpha.h)*A[i]))

      #Feeding on detritus
      f.det<-((1/10^x)*(A.v*10^(x*0.75)*W[i]))
      
      #Feeding on detritus by herbivores 
      f.det_H<-((1/10^x)*(pref.det*A.h*10^(x*alpha.h)*W[i]))


#------------------------
##New growth integrals

      #Growth of predators
      GG.u[,i] <- K.u*(f.pel) + K.h*(f.herb) + K.v*(f.ben)

      #Growth of herbivores
      GG.h[,i] <- K.a*(f.alg) + K.d*(f.det_H) #This now has a lower conversion efficiency for detritus

      #Growth of invertebrates 
      GG.v[,i] <- K.d*(f.det) #This now has a lower conversion efficiency for detritus

#------------------------

      # Predator death integrals 
      # predation mortality
      PM.u[,i]<-as.vector(death.u(U[,i],A.u,expax,mphi,a.u[,i],pref.pel))
      Z.u[,i]<- PM.u[,i]+ OM.u + SM.u + Fvec_pred
      
      # Herbivore death integral
      PM.h[,i]<-as.vector(death.h(U[,i],A.u,expax,mphi,a.h[,i],pref.herb))
      Z.h[,i]<-PM.h[,i]+ OM.h + SM.h + Fvec_herb
     
      # Benthos death integral 
      PM.v[,i]<-as.vector(death.v(U[,i],A.u,expax,mphi,a.v,pref.ben))
      Z.v[,i]<-PM.v[,i]+ OM.v + SM.v
      
#------------------------
##New reproduction values - optional

#R.u[,i] <- r.u*(f.pel) + r.h*(f.herb) + r.v*(f.ben)

#R.h[,i] <- r.a*(f.alg) 

#R.v[,i] <- r.d*(f.det) 

#------------------------

      # biomass density eaten by pred (g.m-3.yr-1)
      eatenbypred<-10^x*(a.u[,i]*pref.pel*A.u*expax*(U[,i]%*%gphi)*U[,i] + a.v*pref.ben*A.u*expax*(V[,i]%*%gphi)*U[,i] + a.h[,i]*pref.herb*A.u*expax*(H[,i]%*%gphi)*U[,i])
      
                        (a.u[,i]*pref.pel*A.u*expax)*(U[,i]*dx)%*%(gphi)


      # biomass density defecated by pred (g.m-3.yr-1)
      defbypred<-10^x*(def.high*a.u[,i]*pref.pel*A.u*expax*(U[,i]%*%gphi)*U[,i] + def.low*a.v*pref.ben*A.u*expax*(V[,i]%*%gphi)*U[,i] + 
      def.low*a.h[,i]*pref.herb*A.u*expax*(H[,i]%*%gphi)*U[,i])
      
      Yp[Fref:end,i]<-Fvec_pred[Fref:end]*U[Fref:end,i]*10^x[Fref:end]
      Yh[Fref:end,i]<-Fvec_herb[Fref:end]*H[Fref:end,i]*10^x[Fref:end]
            
      ##Detritus Biomass Density Pool - fluxes in and out (g.m-3.yr-1) of detritus pool and solve for detritus biomass density in next time step 
      
      #considering pelagic faeces as input as well as dead bodies from both pelagic and benthic communities 
      #and phytodetritus (dying sinking phytoplankton)
      
      if (det.coupling==1.0) input.w<-sinking.rate*(sum(defbypred[ref:end]*dx)+ sum(OM.u[1:end]*U[1:end,i]*10^(x[1:end])*dx) + sum(SM.u[1:end]*U[1:end,i]*10^(x[1:end])*dx)
                                                    + sum(OM.h[1:end]*H[1:end,i]*10^(x[1:end])*dx) + sum(SM.h[1:end]*H[1:end,i]*10^(x[1:end])*dx)) +
                                                      sum(OM.v[1:end]*V[1:end,i]*10^(x[1:end])*dx) + sum(SM.v[1:end]*V[1:end,i]*10^(x[1:end])*dx) 
                                                    #******error fixed 18.11.14 - was missing an *dx - end of second line***********************      

      
      if (det.coupling==0.0) input.w<-sum(OM.v[1:end]*V[1:end,i]*10^(x[1:end])*dx) + sum(SM.v[1:end]*V[1:end,i]*10^(x[1:end])*dx) 
      
      
      output.w<-out.w(v=V[,i],A.v=A.v, w=W[i], h=H[,i], A.h=A.h, pref.det=pref.det, alpha.h=alpha.h) 
      
      # change in detritus biomass density (g.m-3.yr-1)
      dW<-input.w - (output.w) 
      #print(dW)
      
      W[i+1]=W[i] + dW*dt        #biomass density of detritus g.m-2
      
#---------------
#Algal dynamics
      
      output.a <- out.a(a = A[i], A.h = A.h, alpha.h = alpha.h, h = H[,i], pref.alg=pref.alg)

      input.a = alr
      da <- input.a - (output.a)
      
      A[i+1]= A[i]+ da*dt

#---------------

      # Pelagic Predator Density (nos.m-3)- solve for time + dt using implicit time Euler upwind finite difference (help from Ken Andersen and Richard Law)
      
      # Matrix setup for implict differencing
      Ai.u<-Bi.u<-Si.u<-array(0, c(length(x), 1))   
      
      idx=(ref+1):end  #shorthand for matrix referencing
      #This is just the locations in the matrix of predatory fish bigger than plankton, and bigger than eggs (hence the +1)
      
      Ai.u[idx]<-(1/log(10))*-GG.u[idx-1,i]*dt/dx
      Bi.u[idx]<-1+(1/log(10))*GG.u[idx,i]*dt/dx +Z.u[idx,i]*dt
      Si.u[idx]<-U[idx,i]
      
      # Boundary condition at upstream end
      Ai.u[ref]<-0
      Bi.u[ref]<-1
      Si.u[ref]<-U[ref,1]
      
      # Invert matrix
      # hold plankton constant thru time on left hand boundary
      U[1:(ref-1),i+1]<-U[1:(ref-1),1]
      
      #hold small inverts constant through time
      V[refinv:(ref.det-1),i+1]<-V[refinv:(ref.det-1),1]
      
       
#-------------------------
##New recruitment function

#if (pred.stock.recruitment==T) {

      #U[ref,i+1]<-U[ref,i] + (sum(R.u[(ref+1):end,i]*10^x[(ref+1):end]*U[(ref+1):end,i]*dx)*dt)/(dx*10^x[ref]) - (dt/dx)*(1/log(10))*(GG.u[ref,i])*U[ref,i]-dt*Z.u[ref,i]*U[ref,i]
      #energy for recruitment for each size (bigger than egg, hence ref+1) multiplied by its size, by its abundance
  ##**********Need to work through the last part of this - don't really get all the dt / dx stuff

#} else {
  
      ##Original code for no recruitment from stock
      #recruitment at smallest consumer mass  - hold constant       
      U[ref,i+1]<-U[ref,1]
      
#}
#-------------------------      
      
      
      #loop calculation
      for (j in (ref+1):(end)){
        
        #U[j,i+1]<-U[j,i]-(dt/dx)*(1/log(10))*(GG.u[j,i])*U[j,i] + (dt/dx)*(1/log(10))*GG.u[j-1,i]*U[j-1,i] -dt*Z.u[j,i]*U[j,i]
        U[j,i+1]<-(Si.u[j]-Ai.u[j]*U[j-1,i+1])/Bi.u[j]
      }
      
      
      ## Herbivore Density (nos.m-3) - same algorithm as above
      
      Ai.h<-Bi.h<-Si.h<-array(0, c(length(x), 1))   
      idx=(ref.herb+1):end  #shorthand for matrix referencing
      
      Ai.h[idx]<-(1/log(10))*-GG.h[idx-1,i]*dt/dx
      Bi.h[idx]<-1+(1/log(10))*GG.h[idx,i]*dt/dx +Z.h[idx,i]*dt
      Si.h[idx]<-H[idx,i]
      
      # boundary condition at upstream end
      
      Ai.h[ref.herb]<-0
      Bi.h[ref.herb]<-1
      Si.h[ref.herb]<-H[ref.herb,1]  
      
      # invert matrix

#-------------------------
##New recruitment function

      #H[ref.herb,i+1]<-H[ref.herb,i]+ sum(R.h[(ref.herb+1):end,i]*10^x[(ref.herb+1):end]*H[(ref.herb+1):end,i]*dx)*dt/(dx*10^x[ref.herb]) - (dt/dx)*(1/log(10))*(GG.h[ref.herb,i])*H[ref.herb,i]-dt*Z.h[ref.herb,i]*H[ref.herb,i]
      
      ##Original code for no recruitment from stock
      # recruitment at smallest detritivore mass  - hold constant 
      H[ref.herb,i+1]<-H[ref.herb,1]
#-------------------------

      #loop calculation
      for (j in (ref.herb+1):(end)){
        
        #V[j,i+1]<-V[j,i]-(dt/dx)*(1/log(10))*(GG.v[j,i])*V[j,i] + (dt/dx)*(1/log(10))*GG.v[j-1,i]*V[j-1,i]-dt*Z.v[j,i]*V[j,i]
        H[j,i+1]<-(Si.h[j]-Ai.h[j]*H[j-1,i+1])/Bi.h[j]
      }		
      
      
      ## Benthic Detritivore Density (nos.m-3) - same algorithm as above
      
      Ai.v<-Bi.v<-Si.v<-array(0, c(length(x), 1))   
      idx=(ref.det+1):inv_end  #shorthand for matrix referencing
      
      Ai.v[idx]<-(1/log(10))*-GG.v[idx-1,i]*dt/dx
      Bi.v[idx]<-1+(1/log(10))*GG.v[idx,i]*dt/dx +Z.v[idx,i]*dt
      Si.v[idx]<-V[idx,i]
      
      # boundary condition at upstream end
      
      Ai.v[ref.det]<-0
      Bi.v[ref.det]<-1
      Si.v[ref.det]<-V[ref.det,1]  
      
      # invert matrix
      
#-------------------------
##New recruitment function

      #V[ref.det,i+1]<-V[ref.det,i]+ sum(R.v[(ref.det+1):end,i]*10^x[(ref.det+1):end]*V[(ref.det+1):end,i]*dx)*dt/(dx*10^x[ref.det]) - (dt/dx)*(1/log(10))*(GG.v[ref.det,i])*V[ref.det,i]-dt*Z.v[ref.det,i]*V[ref.det,i]
  
      ##Original code for no recruitment from stock
      # recruitment at smallest detritivore mass  - hold constant 
      V[ref.det,i+1]<-V[ref.det,1]      
#-------------------------
    
      
      #loop calculation
      for (j in (ref.det+1):(inv_end)){
        
        #V[j,i+1]<-V[j,i]-(dt/dx)*(1/log(10))*(GG.v[j,i])*V[j,i] + (dt/dx)*(1/log(10))*GG.v[j-1,i]*V[j-1,i]-dt*Z.v[j,i]*V[j,i]
        V[j,i+1]<-(Si.v[j]-Ai.v[j]*V[j-1,i+1])/Bi.v[j]
      }  	
      
    } 
      
    
    #Biomasses through time
    Pred_bio<-colSums(U[ref:end,]*dx*10^x[ref:end]) 
    Herb_bio<-colSums(H[ref.herb:end,]*dx*10^x[ref.herb:end]) 
    Inv_bio<-colSums(V[refinv:inv_end,]*dx*10^x[refinv:inv_end]) 
    
    #To match data of fish > 5cm
    Pred_dat_bio_5<-colSums(U[ref.dst_5:ref.den,]*dx*10^x[ref.dst_5:ref.den])
    Herb_dat_bio_5<-colSums(H[ref.dst_5:ref.den,]*dx*10^x[ref.dst_5:ref.den])
    
    #To match data of fish > 5cm
    Pred_dat_bio_10<-colSums(U[ref.dst_10:ref.den,]*dx*10^x[ref.dst_10:ref.den])
    Herb_dat_bio_10<-colSums(H[ref.dst_10:ref.den,]*dx*10^x[ref.dst_10:ref.den])
    
    #Productivity through time
    P_prod_gr<-colSums(10^x[ref:end]*GG.u[ref:end,]*U[ref:end,]*dx)
    H_prod_gr<-colSums(10^x*GG.h*H*dx)
    I_prod_gr<-colSums(10^x*GG.v*V*dx)
            
    #Fisheries productivity through time
    FP_prod_gr<-colSums(10^x[Fref:end]*GG.u[Fref:end,]*U[Fref:end,]*dx)
    FH_prod_gr<-colSums(10^x[Fref:end]*GG.h[Fref:end,]*H[Fref:end,]*dx)
    
    #Algal dynamics through time
    Al_bio <- A
    Det_bio <- W
    
    #Consumption rates through time
    totaleatenbypred<-sum(eatenbypred[(ref+1):end]*dx)
    eatenbysize <- eatenbypred[(ref+1):end]#*dx
    totaleatenbyben<-sum((1/K.d)*10^x[(ref.det+1):end]*GG.v[(ref.det+1):end,N-1]*V[(ref.det+1):end,N-1]*dx)   #consumption rate g.m-3.yr-1

    #Proportional refuges

     return(list(Body.size = x, Preds=U[,N-1], P_grth=GG.u[, N-1], P_mrt=PM.u[, N-1],
                  Herbs=H[,N-1], H_grth=GG.h[,N-1], H_mrt=PM.h[,N-1],
                  Invs=V[,N-1], I_grth=GG.v[,N-1], I_mrt=PM.v[,N-1], 
                  Pred_gm = Pred_bio[length(Pred_bio)], Herb_gm = Herb_bio[length(Herb_bio)], Inv_gm = Inv_bio[length(Inv_bio)],
                  Pred_dat_gm_5 = Pred_dat_bio_5[length(Pred_dat_bio_5)], Herb_dat_gm_5 = Herb_dat_bio_5[length(Herb_dat_bio_5)],
                  Pred_prod = P_prod_gr[length(P_prod_gr)-1], Herb_prod = H_prod_gr[length(H_prod_gr)-1], Inv_prod = I_prod_gr[length(I_prod_gr)-1],
                  Fpred_prod = FP_prod_gr[length(FP_prod_gr)-1], Fherb_prod = FH_prod_gr[length(FH_prod_gr)-1],  
                  Yield_p = Yp[,N-1], Yield_h = Yh[,N-1], W = W[N-1], A = A[N-1], p_bio_time = Pred_bio, h_bio_time = Herb_bio, i_bio_time = Inv_bio, 
                  A_time = Al_bio, Det_time = Det_bio,
                  Pred_consump = totaleatenbypred,
                  pred_by_size = eatenbysize,
                  Invert_consump = totaleatenbyben,
                  Pred_ref = a.u[,N-1],
                  Herb_ref = a.h[,N-1],
                  params=params))
      
      
    })  
  # end with(params)
  
  
}#end size-based model function