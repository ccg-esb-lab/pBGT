using Distributed


@everywhere module ModelParameters
using Distributed





@everywhere cE=1e6;   #Efficiency (2.5e5)
@everywhere Vmax=(2.5e-8);  #Maximal uptake rate  (1e-5)  ######The smaller the slower that populations grows
@everywhere Km=.25;     #Half-saturation constant (0.05) #

@everywhere cATP=cE*Vmax*40 #*3.5e-5  #Critical ATP concentration for division (3.5e-6) for min

@everywhere mu=19.    #mean plasmid copy number
@everywhere cv=0.3
@everywhere p_Ndist=Normal(mu,mu*cv)
@everywhere starting_variation=true


@everywhere rep_delay=3  #chances of plasmid replication in time step 5 super ok, 3 ok


#@everywhere knoise=0.5
@everywhere kNoise_dist=Normal(0.,.1)
@everywhere antibiotic_action=0.75
@everywhere dNormATP=Normal(0.5,0.15)
@everywhere pCost=0.003   #Plasmid Cost (.003 -> 19 plasmids) (1-.943)/19

@everywhere mA=429.03790087463557  #Slope of the linear regression of the antbiotic MIC
@everywhere bA=42.41399416909739 #Intercept of the linerar regression of the antibiotic MIC

@everywhere function getMIC(mui)
mic=round((mA*mui)+bA)
mui==0. ? mic=4. : 0
return mic
end

@everywhere pMIC=getMIC(mu)



@everywhere antibiotic_deg=true

@everywhere degAA=1e-10     #Degradation rate of antibiotic 
#@everywhere degAA=0     #Degradation rate of antibiotic 
#

@everywhere iniCells=1000

@everywhere verbose=false;

################### This parameters for faster runs and less bacteria 
#@everywhere degAA=1e-8     #Degradation rate of antibiotic 
#@everywhere cE=2.5e5;   #Efficiency (2.5e5)
#@everywhere Vmax=(1e-7);  #Maximal uptake rate  (1e-5)  ######The smaller the slower that populations grows
#@everywhere Km=.25;     #Half-saturation constant (0.05)
#@everywhere gparams=[cE,Vmax,Km];
#@everywhere cATP=cE*Vmax*45 #*3.5e-5  #Critical ATP concentration for division (3.5e-6) for min





end