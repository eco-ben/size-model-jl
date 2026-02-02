
#---------------------------------------------------------------------------------------
# Parameters that can vary
#---------------------------------------------------------------------------------------

refuges  =  CSV.read("test_refuges.csv", DataFrame)
refuge  =  refuges[:,8]

#Algal dynamics
alr = 109#0.2 #Maximum algal growth rate per year
AK = 6 #Carrying capacity of algae
min_A=0#0.1               #proportion of population with reduced vulnerability (range: 1-0.1, 1 = no complexity)

#Primary productivity and search rates
pp=-0.5 #1.07            #Primary productivity 
#A_u = 10
A_u_reef= 6.4 #10            #Yearly rate volume searched by predators pre settlement
A_u_pel = 10            #Yearly rate volume searched by predators post settlement
A_u = 6.4
A_v= 0.1 #A_u_reef / 10               #Yearly rate volume searched by invertebrates (should be 10x less than predators?)
A_h= 0.2 #1.25                   #Yearly rate volume searched by herbivores
flow = 10

#Theoretical complexity
setsp=1000        #slope of decrease (settlement onto reef structure)
emsp=100          #increasing slope (emergence from reef structure)
settle=0.1        #weight when settle to reef structure
emerge=1000       #maximum refuge size in body mass (range explored = 1-1500) 
min_A=0

#Feeding preferences
pref_pel=0.33           #Preference for "pelagic" carnivorous prey 
pref_ben=0.5            #Preference for "benthic" invertebrate prey
pref_herb=0.17          #Preference for herbivorous prey

pref_det = 0
pref_alg = 1

#Fishing
Fmort_pred=0            #Level of fishing mortality rate yr^(-1)
Fmort_herb=0
min_fishing_size=1      #Minimum log10 body size fished
max_fishing_size=3.5    #Maximum log10 body size fished

#---------------------------------------------------------------------------------------
# Fixed parameters
#---------------------------------------------------------------------------------------

q0=2.0                  #Mean log10 predator prey mass ratio  100:1.(beta in paper)
sd_q=1.0                #0.6 to 1.0 for lognormal prey preference function. 5=B&R ref value (sigma, standard deviation)
qmax=q0 + 2*sd_q        #Truncation of the feeding kernel 95% distribution between +/- 2*sd_q
qmin=q0 - 2*sd_q

det_coupling=1.0
sinking_rate=0.8        # 35% of PP material reaches the seafloor davies & Payne 1983 Marine biology
alpha=0.75#0.82              # exponent for metabolic requirements plus swimming for predators(Ware et al 1978)
                            # NOTE:this exponent =0.75 for whole organism basal (sedentary) metabolic rate (see growth_v) from Peters (1983) and Brown et al. (2004) for detritivores
alpha_h=0.75#0.6046

K_u=0.15                 #Gross growth conversion efficiency for organisms in the "predator" spectrum Ware (1978)
K_v=0.15                 #Gross growth conversion efficiency for organisms in the "detritivore" spectrum
K_h=0.15                 #Gross growth conversion efficiency for organisms in the "herbivore" spectrum
K_d=0.1                 #Gross growth conversion efficiency for detritus
def_high=0.4            #Fraction of ingested food that is defecated (Peters,1983)
def_low=0.4             #Low = low quality (K) food, high = high quality (K) food

K_a=0.1

mu0=0.2	                #Residual natural mortality
k_sm =0.1               #Constant for senescence mortality 
xs=3                    #Size at sensenscence e
p_s=0.3  			          #Exponent of senescence

#Initial size spectra properties
#ui0=10^(pp)      #Initial intercept of plankton size spectrum, from P_McCloghrie ERSEM output.
#vi0=sinking_rate*ui0   #Initial intercept of detritivore size spectrum  - same fraction as sinking rate
r_plank=-1.0            #Initial (fixed) slope of phyto-zooplankton spectrum
pred_slope=-1#-0.75           #Initial slope of fish spectrum, same as that used in B&R
herb_slope= -1#-0.64          #Initial slope of detrtivore spectrum
invert_slope=-0.36
                                
                            #Proportion of initial plankton intercept that depicts the abundance of larvae / eggs of 
pred_prod = 1#(10^0.56)
herb_prod= 1#(10^-0.04) #1.65           #Herbivorous fish
invert_prod = 10^2.1#(10^-0.49) # 5.4
#sinking_rate         #Invertebrates

W_init=5                #Initial detritus density per scenario
A_init=25

#-------------------------------------------
##New parameters for reproductive investment

r_u=0.1
r_h=0.1
r_v=0.1
r_a=0.1
r_d=0.1


#---------------------------------------------------------------------------------------
# Parameters for numerical integration
#---------------------------------------------------------------------------------------
dx=0.1                  #Size increment after discretization for integration (log body weight) 
xmin=-12                #Minimum log10 body size of plankton
x1=-1.5#-3.5                   #Minimum log10 body size in predators
x1_det=-4               #Minimum log10 body size in dynamics benthic detritivores
x1_herb= -1.5 #-3.5           #Minimum log10 body size of herbivores
xmax=3.5                #Maximum log10 body size of predators
xmax2=3                 #Maximum log10 body size before senescence kicks in (departure form linearity)

##******NEW
inv_max = 3.5
invmin = -6


## Vector with size bins 
x    = xmin:dx:xmax
y    = x
size_end  = length(x)

refinv  =  first(findall(collect(x) .== invmin))

#******NEW

inv_end  =  first(findall(collect(x) .== inv_max))

tstepdays=1.0                             #Timestep in days
tmaxyears=10                           #Maximum number of years for model run

N= Int64(365.0*tmaxyears/tstepdays)	  #Total number of steps in a period
dt=tmaxyears/N                            #Time step (rates per year)   

ref=((x1-xmin)/dx)+1                      #Location of smallest size consumer
#phyto_ref=((-6-xmin)/dx)+1                #Position in x of upper size of phytoplankton
ref_det=((x1_det-xmin)/dx)+1              #Position in x of x1.det
ref_herb=((x1_herb-xmin)/dx)+1            #Position in x of x1.herb

ref2=((xmax2-xmin)/dx)+1                  #Position in x of xmax2 (largest benthic detritivores)

#Setting data ranges for field measurement

#Considering fish bigger than 5cm - or log10 weight > 0
dat_start_5=0
#Considering fish bigger than 10cm - or log10 weight > 1
dat_start_10=1
dat_end = 3.5#4.5

ref_dst_5 = ((dat_start_5-xmin)/dx)+1      # position in x of fish bigger than ~5cm (surveyable)
ref_dst_10 = ((dat_start_5-xmin)/dx)+1      # position in x of fish bigger than ~10cm (surveyable)

ref_den = ((dat_end-xmin)/dx)+1                          # position in x of largest surveyed fish

Fref=((min_fishing_size-xmin)/dx)+1   # position in F vector corresponding to smallest size fished
Fref2=((max_fishing_size-xmin)/dx)+1