# function to run size spectrum model from https://github.com/alicerogers/Fisheries-and-degrading-reefs

function run_model(params, initial_run = false, fish_herbs = true)
    ui0 = 10^pp
    Int64(Fref)=((min_fishing_size-xmin)/dx)+1
    
    # Function to construct the spectrum at the start of integration
    function u_init_f(x, ui0, init_slope) return ui0*10 .^ collect(init_slope*x) end
    
    # Function to build a lookup table for diet preference of all combinations of predator 
    # and prey body size: diet preference (in the predator spectrum only)
    
    function phi_f(q; q0=q0, sd_q=sd_q, qmax=qmax, qmin=qmin)
      
      phi = ifelse.(q .> 0, exp.(-(q .- q0).^2 ./ (2 .* sd_q.^2)) ./ (sd_q .* sqrt(2π)), 0)
      return phi
    end
    
    # Functions to build lookup tables for components of integration which remain constant
    function gphi_f(q) return 10^(-q)       * phi_f(q ) end  #growth
    function mphi_f(q2; alpha=alpha) return 10^(alpha*q2) * phi_f(q2) end	#mortality
    
    # Function to build lookup table for components of 10^(alpha*x)
    function expax_f(x, alpha) return 10 .^ collect(alpha*x ) end
    
    # Function to compute convolution products: can be for growth or mortality depending on phi
    function convolution_f(phi, u)
    res = matrix(NA,length(x),1)
        for i in 1:length(res)
            res[i] = 0.0
            res[i] = sum(phi[:,i] * u * dx)
            
        end
        return res
    end
    
      
    #---------------------
    #Availability function
    #---------------------
  
  
    function vulnerability(refuge_density, competitor_density; x = x, min_A=min_A)
        res = zeros(Float64, length(x))
     
        for i in 1:length(x)
        
            if competitor_density[i] == 0 
                res[i] = 1
            else 
                res[i] = 1 - (refuge_density[i] / competitor_density[i])
            end

            if res[i] < min_A
                res[i] = min_A
            else 
                res[i] = res[i]
            end
        end
        return res
    end
        
    #--------------------------------
    # Growth and Mortality Equations
    #--------------------------------
    
    ## predators
    function death_u(u,A,expax,mphi,a_u,pref_pel) return (a_u .* pref_pel .* A .* expax) .* ((u .* dx)' * mphi |> vec) end    # use faster matrix method instead of convolution loop function
     
    ## herbivores
    function death_h(u,A,expax,mphi,a_h,pref_herb) return (a_h .* pref_herb .* A .* expax) .* ((u .* dx)' * mphi |> vec) end
    
    ## detritivores
    function death_v(u,A,expax,mphi,a_v,pref_ben) return (a_v .* pref_ben .* A .* expax) .* ((u .* dx)' * mphi |> vec) end

    
    ## detritus output (g_m-3.yr-1)
    function out_w(A_v,A_h,v,w,h,alpha_h,pref_det) return (sum((A_v .* (10 .^ (x*0.75)) .* w .* v) .* dx) + sum((A_h .* pref_det .* (10 .^ (x*alpha_h)) .* w .* h) .* dx)) end # biomass density removed by detritivores per year    
    ##Added a pref_det and pref_alg function so that herbivores are both feeding on both types of resource
    
#----------------------------------
#Algal grazing dynamics

    ## algae output (g_m-2.yr-1)
    function out_a(A_h, alpha_h, h, a, pref_alg) return (sum(pref_alg .* A_h .* (10 .^ (x*alpha_h) .* a .* h) .* dx)) end #biomass density removed by herbivores ?per day if using daily data or per hyear if using yearly data? 
            
    #---------------------------------------------------------------------------------------
    # Initialising matrices
    #---------------------------------------------------------------------------------------
    
      
    # q1 is a square matrix holding the log(predatorsize/preysize) for all combinations of sizes
    q1  = Matrix{Union{Float64, Missing}}(missing, length(x), length(y))
    for i in 1:length(y) 
        q1[:,i] = y[i] .- x 
    end
    
    # q2 is the reverse matrix holding the log(preysize/predatorsize) for all combinations of sizes
    q2 = -q1	
    
    # matrix for recording the two size spectra over time
    V = U = H   = zeros(Float64, length(x), N)
    
    # vector to hold detrtitus biomass density (g_m-3)
    W = zeros(N)
    
    #Vector to hold algal biomass density (g_m-3)
    A = zeros(N)
    
    # matrix for keeping track of growth and ingested food:
    GG_v = GG_u = GG_h = zeros(length(x), N)  

    # matrix for keeping track of vulnerability through time
    a_u = a_h = zeros(length(x), N)
   
    
##*****
    # Optional matrices for keeping track of reproduction from ingested food
    #R_v = R_u = R_h = array(0, c(length(x), N))
    
    # matrix for keeping track of predation mortality
    PM_v = PM_u = PM_h  = zeros(length(x), N)   
    
    # matrix for keeping track of  total mortality (Z)
    Z_v = Z_u = Z_h = zeros(length(x), N)
    
    # matrix for keeping track of senescence mortality and other (intrinsic) mortality
    SM_v = SM_u  = SM_h  =OM_v = OM_u  = OM_h = zeros(length(x))
    
    # empty vector to hold fishing mortality rates at each size class
    Fvec_pred = zeros(length(x))
    Fvec_herb = zeros(length(x))
    
    #Empty vector to hold fisheries yield
    Yp  = zeros(length(x), N)  
    Yh = zeros(length(x),N)
    
    # short hand for matrix indexing
    idx = 2:size_end
    
    # lookup tables for terms in the integrals which remain constant over time
    gphi  = gphi_f(q1)
    mphi  = mphi_f(q2)
    
    # lookup table for components of 10^(alpha*x)
    expax = expax_f(x, alpha)
    
    # vectors for storing slopes of regression
    slope_v_lm = slope_u_lm = zeros(N-1)
    slope_v_nls = slope_u_nls = zeros(N-1)
    
    #---------------------------------------------------------------------------------------
    # Numerical integration
    #---------------------------------------------------------------------------------------
    
    ##INITIAL VALUES
    
    if initial_run
      
      # if running to equilibrium use given initial values
      
      U[1:(Int64(ref)-1),1] = u_init_f(x,ui0,r_plank)[1:(Int64(Int64(ref))-1)] .* flow    # (phyto+zoo)plankton size spectrum  
      
      U[Int64(ref):end,1] = u_init_f(x,ui0, pred_slope)[Int64(ref):end] 
      
      H[Int64(ref_herb):end,1] = u_init_f(x,ui0, herb_slope)[Int64(ref_herb):end] .* herb_prod
       
      V[Int64(refinv):(Int64(ref_det)-1),1] = u_init_f(x,invert_prod,invert_slope)[Int64(refinv):(Int64(ref_det)-1)]
      V[Int64(ref_det):inv_end,1] = u_init_f(x,invert_prod,invert_slope)[Int64(ref_det):inv_end]
      
      W[1] = W_init                                                          # initial detritus biomass density (g_m-3) 
      A[1] = A_init                                                          # initial algal biomass density (g_m-3)
      
    else
      
      # use values from non-complex initial run
      
      U[1:(Int64(ref)-1),1] = u_init_f(x,ui0,r_plank)[1:(Int64(ref)-1)]*flow                    # (phyto+zoo)plankton size spectrum  
      U[Int64(ref):end,1] = initial_res.Preds[Int64(ref):end]                                      # initial consumer size spectrum
      H[Int64(ref_herb):end,1] = initial_res.Herbs[Int64(ref_herb):end]
      V[Int64(ref_det):inv_end,1] = initial_res.Invs[Int64(ref_det):inv_end]   # initial detritivore spectrum
      W[1] = initial_res.W
      A[1] = initial_res.A
      
    end
    
    # intrinsic natural mortality
    OM_u = mu0*10 .^ (-0.25*x)      
    OM_v = mu0*10 .^ (-0.25*x)      
    OM_h = mu0*10 .^ (-0.25*x)
    
    # senescence mortality rate to limit large fish from building up in the system
    # same function as in Law et al 2008, with chosen parameters gives similar M2 values as in Hall et al. 2006
    SM_u = k_sm*10 .^ (p_s*(x .- xs))
    SM_v = k_sm*10 .^ (p_s*(x .- xs))
    SM_h = k_sm*10 .^ (p_s*(x .- xs))
    
    
    #Fishing mortality at each size
    Fvec_pred[Int64(Fref):end] = fill(Fmort_pred,length(Fvec_pred[Int64(Fref):end]))  
    Fvec_herb[Int64(Fref):end] = fill(Fmort_herb,length(Fvec_herb[Int64(Fref):end]))
    
    # iteration over time
    @showprogress for i in 1:(N-1)
        ##----------------------------------------
        ##Density dependent vulnerability function
        a_u[:,i]  =  vulnerability(refuge, U[:,i]+H[:,i])      
        a_h[:,i]  =  vulnerability(refuge, H[:,i]+U[:,i])
        a_v  =  1

        #------------------------
        ##Feeding rate functions
      
        #Feeding on predators
        f_pel  =  (a_u[:, i] .* pref_pel .* A_u .* expax) .* ((U[:, i] .* dx)' * gphi |> vec)
    
        #Feeding on herbivores
        f_herb  =  (a_h[:,i] .* pref_herb .* A_u .* expax) .* ((H[:,i] .* dx)' * gphi |> vec)      
        
        #Feeding on invertebrates
        f_ben  =  (pref_ben .* A_u .* expax) .* ((V[:,i] .* dx)' * (gphi) |> vec)
        
        #Feeding on algae
        f_alg = (((1/10) .^ x) .* (pref_alg*A_h*10 .^ (x*alpha_h) .* A[i]))

        #Feeding on detritus
        f_det = (((1/10) .^ x) .* (A_v*10 .^ (x*0.75) .* W[i]))
        
        #Feeding on detritus by herbivores 
        f_det_H = (((1/10) .^ x) .* (pref_det*A_h*10 .^ (x*alpha_h) .* W[i]))


        #------------------------
        ##New growth integrals

        #Growth of predators
        GG_u[:,i]  =  K_u .* (f_pel) .+ K_h .* (f_herb) .+ K_v .* (f_ben)

        #Growth of herbivores
        GG_h[:,i]  =  K_a .* (f_alg) .+ K_d .* (f_det_H) #This now has a lower conversion efficiency for detritus

        #Growth of invertebrates 
        GG_v[:,i]  =  K_d .* (f_det) #This now has a lower conversion efficiency for detritus

        #------------------------

      # Predator death integrals 
      # predation mortality
      PM_u[:,i] = death_u(U[:,i],A_u,expax,mphi,a_u[:,i],pref_pel)
      Z_u[:,i] =  PM_u[:,i] .+ OM_u .+ SM_u .+ Fvec_pred
      
      # Herbivore death integral
      PM_h[:,i] = death_h(U[:,i],A_u,expax,mphi,a_h[:,i],pref_herb)
      Z_h[:,i] = PM_h[:,i] .+ OM_h .+ SM_h .+ Fvec_herb
     
      # Benthos death integral 
      PM_v[:,i] = death_v(U[:,i],A_u,expax,mphi,a_v,pref_ben)
      Z_v[:,i] = PM_v[:,i] .+ OM_v .+ SM_v
      
        #------------------------
        ##New reproduction values - optional

        #R_u[:,i]  =  r_u*(f_pel) + r_h*(f_herb) + r_v*(f_ben)

        #R_h[:,i]  =  r_a*(f_alg) 

        #R_v[:,i]  =  r_d*(f_det) 

        #------------------------

      # biomass density eaten by pred (g_m-3.yr-1)
      eatenbypred = 10 .^ x .* (a_u[:,i] .* pref_pel .* A_u .* expax .* (U[:,i]' * gphi |> vec) .* U[:,i] .+ a_v .* pref_ben .* A_u .* expax .* (V[:,i]' * gphi |> vec) .* U[:,i] .+ a_h[:,i] .* pref_herb .* A_u .* expax .* (H[:,i]' * gphi |> vec) .* U[:,i])
      
    #                     (a_u[:,i]*pref_pel*A_u*expax)*(U[:,i]*dx)*(gphi)


      # biomass density defecated by pred (g_m-3.yr-1)
      defbypred = 10 .^ x .* (def_high .* a_u[:,i] .* pref_pel .* A_u .* expax .* (U[:,i]' * gphi |> vec) .* U[:,i] .+ def_low .* a_v .* pref_ben .* A_u .* expax .* (V[:,i]' * gphi |> vec) .* U[:,i] .+ 
        def_low .* a_h[:,i] .* pref_herb .* A_u .* expax .* (H[:,i]' * gphi |> vec) .* U[:,i])
      
      Yp[Int64(Fref):end,i] = Fvec_pred[Int64(Fref):end] .* U[Int64(Fref):end,i] .* 10 .^ x[Int64(Fref):end]
      Yh[Int64(Fref):end,i] = Fvec_herb[Int64(Fref):end] .* H[Int64(Fref):end,i] .* 10 .^ x[Int64(Fref):end]
            
      ##Detritus Biomass Density Pool - fluxes in and out (g_m-3.yr-1) of detritus pool and solve for detritus biomass density in next time step 
      
      #considering pelagic faeces as input as well as dead bodies from both pelagic and benthic communities 
      #and phytodetritus (dying sinking phytoplankton)
      
        if det_coupling==1.0
            input_w = sinking_rate * (
                sum(defbypred[Int64(ref):end]*dx) + sum(OM_u[1:end] .* U[1:end,i] .* 10 .^ (x[1:end]) .* dx) + sum(SM_u[1:end] .* U[1:end,i] .* 10 .^ (x[1:end]) .* dx)
                                                    + sum(OM_h[1:end] .* H[1:end,i] .* 10 .^ (x[1:end]) .* dx) + sum(SM_h[1:end] .* H[1:end,i] .* 10 .^ (x[1:end]) .* dx)) +
                                                      sum(OM_v[1:end] .* V[1:end,i] .* 10 .^ (x[1:end]) .* dx) + sum(SM_v[1:end] .* V[1:end,i] .* 10 .^ (x[1:end]) .* dx
            ) 
                                                    #**error fixed 18.11.14 - was missing an *dx - end of second line*********    
        end  

      
        if det_coupling==0.0 
            input_w = sum(OM_v[1:end] .* V[1:end,i] .* 10 .^ (x[1:end]) .* dx) + sum(SM_v[1:end] .* V[1:end,i] .* 10 .^ (x[1:end]) .* dx) 
        end
      
      
      output_w = out_w(A_v, A_h, V[:, i], W[i], H[:, i], alpha_h, pref_det) 
      
      # change in detritus biomass density (g_m-3.yr-1)
      dW = input_w - (output_w) 
      #print(dW)
      
      W[i+1]=W[i] + dW*dt        #biomass density of detritus g_m-2
      
        #---------------
        #Algal dynamics
      
      output_a  =  out_a(A_h, alpha_h, H[:, i], A[i], pref_alg)

      input_a = alr
      da  =  input_a - (output_a)
      
      A[i+1]= A[i]+ da*dt

        #---------------

      # Pelagic Predator Density (nos_m-3)- solve for time + dt using implicit time Euler upwind finite difference (help from Ken Andersen and Richard Law)
      
      # Matrix setup for implict differencing
      Ai_u = Bi_u = Si_u = zeros(Float64, length(x), 1)   
      
      idx=(Int64(ref)+1):size_end  #shorthand for matrix referencing
      #This is just the locations in the matrix of predatory fish bigger than plankton, and bigger than eggs (hence the +1)
      
      Ai_u[idx] = (1/log(10)) .* -GG_u[idx .- 1,i] .* (dt/dx)
      Bi_u[idx] = 1 .+ (1/log(10)) .* GG_u[idx,i] .* (dt/dx) .+ Z_u[idx,i] .* dt
      Si_u[idx] = U[idx,i]
      
      # Boundary condition at upstream end
      Ai_u[Int64(ref)] = 0
      Bi_u[Int64(ref)] = 1
      Si_u[Int64(ref)] = U[Int64(ref),1]
      
      # Invert matrix
      # hold plankton constant thru time on left hand boundary
      U[1:(Int64(ref)-1),i+1] = U[1:(Int64(ref)-1),1]
      
      #hold small inverts constant through time
      V[Int64(refinv):(Int64(ref_det)-1),i+1] = V[Int64(refinv):(Int64(ref_det)-1),1]
      
       
        #-------------------------
        ##New recruitment function

        #if (pred_stock_recruitment==T) {

            #U[Int64(ref),i+1] = U[Int64(ref),i] + (sum(R_u[(Int64(ref)+1):end,i]*10^x[(Int64(ref)+1):end]*U[(Int64(ref)+1):end,i]*dx)*dt)/(dx*10^x[Int64(ref)]) - (dt/dx)*(1/log(10))*(GG_u[Int64(ref),i])*U[Int64(ref),i]-dt*Z_u[Int64(ref),i]*U[Int64(ref),i]
            #energy for recruitment for each size (bigger than egg, hence Int64(ref)+1) multiplied by its size, by its abundance
        ##****Need to work through the last part of this - don't really get all the dt / dx stuff

        #} else {
  
      ##Original code for no recruitment from stock
      #recruitment at smallest consumer mass  - hold constant       
      U[Int64(ref),i+1] = U[Int64(ref),1]
      
        #}
        #-------------------------      
      
      
      #loop calculation
      for j in (Int64(ref)+1):size_end
        
        #U[j,i+1] = U[j,i]-(dt/dx)*(1/log(10))*(GG_u[j,i])*U[j,i] + (dt/dx)*(1/log(10))*GG_u[j-1,i]*U[j-1,i] -dt*Z_u[j,i]*U[j,i]
        U[j,i+1] = (Si_u[j]-Ai_u[j]*U[max(1, j-1),i+1]) / Bi_u[j]
      end
      
      
      ## Herbivore Density (nos_m-3) - same algorithm as above
      
      Ai_h = Bi_h = Si_h = zeros(length(x), 1)   
      idx=(Int64(ref_herb)+1):size_end  #shorthand for matrix referencing
      
      Ai_h[idx] = (1/log(10)) .* -GG_h[idx .- 1,i] .* (dt/dx)
      Bi_h[idx] = 1 .+ (1/log(10)) .* GG_h[idx,i] .* dt/dx .+ Z_h[idx,i] .* dt
      Si_h[idx] = H[idx,i]
      
      # boundary condition at upstream end
      
      Ai_h[Int64(ref_herb)] = 0
      Bi_h[Int64(ref_herb)] = 1
      Si_h[Int64(ref_herb)] = H[Int64(ref_herb),1]  
      
      # invert matrix

        #-------------------------
        ##New recruitment function

            #H[Int64(ref_herb),i+1] = H[Int64(ref_herb),i]+ sum(R_h[(Int64(ref_herb)+1):end,i]*10^x[(Int64(ref_herb)+1):end]*H[(Int64(ref_herb)+1):end,i]*dx)*dt/(dx*10^x[Int64(ref_herb)]) - (dt/dx)*(1/log(10))*(GG_h[Int64(ref_herb),i])*H[Int64(ref_herb),i]-dt*Z_h[Int64(ref_herb),i]*H[Int64(ref_herb),i]
            
            ##Original code for no recruitment from stock
            # recruitment at smallest detritivore mass  - hold constant 
            H[Int64(ref_herb),i+1] = H[Int64(ref_herb),1]
        #-------------------------

        #loop calculation
        for j in (Int64(ref_herb)+1):size_end
            
            #V[j,i+1] = V[j,i]-(dt/dx)*(1/log(10))*(GG_v[j,i])*V[j,i] + (dt/dx)*(1/log(10))*GG_v[j-1,i]*V[j-1,i]-dt*Z_v[j,i]*V[j,i]
            H[j,i+1] = (Si_h[j]-Ai_h[j]*H[max(1, j-1),i+1])/Bi_h[j]
        end
      
      
      ## Benthic Detritivore Density (nos_m-3) - same algorithm as above
      
      Ai_v = Bi_v = Si_v = zeros(length(x), 1)   
      idx=(Int64(ref_det)+1):inv_end  #shorthand for matrix referencing
      
      Ai_v[idx] = (1/log(10)) .* -GG_v[idx .- 1,i] .* dt/dx
      Bi_v[idx] = 1 .+ (1/log(10)) .* GG_v[idx,i] .* dt/dx .+ Z_v[idx,i] .* dt
      Si_v[idx] = V[idx,i]
      
      # boundary condition at upstream end
      
      Ai_v[Int64(ref_det)] = 0
      Bi_v[Int64(ref_det)] = 1
      Si_v[Int64(ref_det)] = V[Int64(ref_det),1]  
      
      # invert matrix
      
        #-------------------------
        ##New recruitment function

            #V[Int64(ref_det),i+1] = V[Int64(ref_det),i]+ sum(R_v[(Int64(ref_det)+1):end,i]*10^x[(Int64(ref_det)+1):end]*V[(Int64(ref_det)+1):end,i]*dx)*dt/(dx*10^x[Int64(ref_det)]) - (dt/dx)*(1/log(10))*(GG_v[Int64(ref_det),i])*V[Int64(ref_det),i]-dt*Z_v[Int64(ref_det),i]*V[Int64(ref_det),i]
        
            ##Original code for no recruitment from stock
            # recruitment at smallest detritivore mass  - hold constant 
            V[Int64(ref_det),i+1] = V[Int64(ref_det),1]      
        #-------------------------
    
      
      #loop calculation
      for j in (Int64(ref_det)+1):(inv_end)
        
        #V[j,i+1] = V[j,i]-(dt/dx)*(1/log(10))*(GG_v[j,i])*V[j,i] + (dt/dx)*(1/log(10))*GG_v[j-1,i]*V[j-1,i]-dt*Z_v[j,i]*V[j,i]
        V[j,i+1] = (Si_v[j]-Ai_v[j]*V[j-1,i+1])/Bi_v[j]
      end
      
    end
      
    
    #Biomasses through time
    Pred_bio = vec(sum(U[Int64(ref):end,:] .* dx .* 10 .^ x[Int64(ref):end], dims = 1)) 
    Herb_bio = vec(sum(H[Int64(ref_herb):end,:] .* dx .* 10 .^ x[Int64(ref_herb):end], dims = 1)) 
    Inv_bio = vec(sum(V[Int64(refinv):inv_end,:] .* dx .* 10 .^ x[Int64(refinv):inv_end], dims = 1)) 
    
    #To match data of fish > 5cm
    Pred_dat_bio_5 = vec(sum(U[Int64(ref_dst_5):Int64(ref_den),:] .* dx .* 10 .^ x[Int64(ref_dst_5):Int64(ref_den)], dims = 1))
    Herb_dat_bio_5 = vec(sum(H[Int64(ref_dst_5):Int64(ref_den),:] .* dx .* 10 .^ x[Int64(ref_dst_5):Int64(ref_den)], dims = 1))
    
    #To match data of fish > 5cm
    Pred_dat_bio_10 = vec(sum(U[Int64(ref_dst_10):Int64(ref_den),:] .* dx .* 10 .^ x[Int64(ref_dst_10):Int64(ref_den)], dims = 1))
    Herb_dat_bio_10 = vec(sum(H[Int64(ref_dst_10):Int64(ref_den),:] .* dx .* 10 .^ x[Int64(ref_dst_10):Int64(ref_den)], dims = 1))
    
    #Productivity through time
    P_prod_gr = vec(sum(10 .^ x[Int64(ref):end] .* GG_u[Int64(ref):end,:] .* U[Int64(ref):end,:] .* dx, dims = 1))
    H_prod_gr = vec(sum(10 .^ x .* GG_h .* H .* dx, dims = 1))
    I_prod_gr = vec(sum(10 .^ x .* GG_v .* V .* dx, dims = 1))
            
    #Fisheries productivity through time
    FP_prod_gr = vec(sum(10 .^ x[Int64(Fref):end] .* GG_u[Int64(Fref):end,:] .* U[Int64(Fref):end,:] .* dx, dims = 1))
    FH_prod_gr = vec(sum(10 .^ x[Int64(Fref):end] .* GG_h[Int64(Fref):end,:] .* H[Int64(Fref):end,:] .* dx, dims = 1))
    
    #Algal dynamics through time
    Al_bio  =  A
    Det_bio  =  W
    
    #Consumption rates through time
    totaleatenbypred = sum(eatenbypred[(Int64(ref)+1):end]*dx)
    eatenbysize  =  eatenbypred[(Int64(ref)+1):end]#*dx
    totaleatenbyben = sum((1/K_d)*10^x[(Int64(ref_det)+1):end]*GG_v[(Int64(ref_det)+1):end,N-1]*V[(Int64(ref_det)+1):end,N-1]*dx)   #consumption rate g_m-3.yr-1

    #Proportional refuges

     return list(Body_size = x, Preds=U[:,N-1], P_grth=GG_u[:, N-1], P_mrt=PM_u[:, N-1],
                  Herbs=H[:,N-1], H_grth=GG_h[:,N-1], H_mrt=PM_h[:,N-1],
                  Invs=V[:,N-1], I_grth=GG_v[:,N-1], I_mrt=PM_v[:,N-1], 
                  Pred_gm = Pred_bio[length(Pred_bio)], Herb_gm = Herb_bio[length(Herb_bio)], Inv_gm = Inv_bio[length(Inv_bio)],
                  Pred_dat_gm_5 = Pred_dat_bio_5[length(Pred_dat_bio_5)], Herb_dat_gm_5 = Herb_dat_bio_5[length(Herb_dat_bio_5)],
                  Pred_prod = P_prod_gr[length(P_prod_gr)-1], Herb_prod = H_prod_gr[length(H_prod_gr)-1], Inv_prod = I_prod_gr[length(I_prod_gr)-1],
                  Fpred_prod = FP_prod_gr[length(FP_prod_gr)-1], Fherb_prod = FH_prod_gr[length(FH_prod_gr)-1],  
                  Yield_p = Yp[:,N-1], Yield_h = Yh[:,N-1], W = W[N-1], A = A[N-1], p_bio_time = Pred_bio, h_bio_time = Herb_bio, i_bio_time = Inv_bio, 
                  A_time = Al_bio, Det_time = Det_bio,
                  Pred_consump = totaleatenbypred,
                  pred_by_size = eatenbysize,
                  Invert_consump = totaleatenbyben,
                  Pred_ref = a_u[:,N-1],
                  Herb_ref = a_h[:,N-1],
                  params=params)  
end