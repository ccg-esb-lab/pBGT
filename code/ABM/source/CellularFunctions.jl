@everywhere module CellularFunctions
using Distributed


@everywhere function plasmids_dynamics(b::Bacteria) 
    
    pR=b.plasmids[1]    #Bacteria expected plasmid copy number (mu)
    pR<0 ? pR=0 : 0
    
    if(b.phenotype==-1) # WildType
        pR=0;
        pA=0;
    elseif(b.phenotype==0)  #Chromosomal Resistance
        pR=0;
        pA=b.plasmids[2]
    else  #Plasmid-bearing
        pR=b.plasmids[1]    #This bacteria max number of plasmids (mu_i)
        pA=b.plasmids[2]   #This bacteria current number of plasmids
        
        pR<0 ? pR=0 : 0
        pA==0 ? pR=0 : 0

        (pA>pR && pR>0) ?  pR=pA : 0   ##check if needed
        if(pA>0)
            for td=1:rep_delay
                
                p_rep=1-(pA/pR)  #Probability of plasmid replication

                rand() < p_rep ? pA=pA+1 : 0
            end  
        end
    end
    pA<0 ? pA=0 : 0
    
    b.plasmids=[pR,pA]
    
end

@everywhere function antibiotic_dynamics(b::Bacteria,aA::Float64)

    pC=b.plasmids[1]
    pA=b.plasmids[2]
    
    rA=0.;
    rA=mA*pA+bA     # Resistance to antibiotic 
    
    if(pC==0&&(pA==0))
        ###sad thing for a man to do
        rA=4 #baseline resistance
        
    else
        #pA==0 ? rA-=bA : 0
        pA==0 ? rA=4 : 0
        
    end
    
    sA=aA/rA        #Succeptibility to antibiotic 
    
    ### Probability of being killed 
    kn=rand(kNoise_dist)
    aliveA= (sA < (1+ (kn)) ? true : false)
    dead=!(aliveA)
    
    return dead
end

@everywhere function grow(b::Bacteria, R,aA)     
    
    pA=b.plasmids[2]
   
    modC=1-(pCost*pA)  #Plasmid Cost
    modC<0 ? modC=0 : 0
    
    c=cE*modC  #Efficiency
    
    
    Vmaxi=Vmax  #Max growth rate
    Kmi=Km #Half-saturation constant
            
    uR=(Vmaxi*R)/(Kmi+R)  #Resource Uptake Function
    gATP=c*uR           #ATP gained
    b.ATP+=gATP         
    
    daA=(pA*degAA)*aA  
        
    return uR,gATP,daA
end



@everywhere function divide(b::Bacteria,lid::Int)
    
    #dNormATP=Normal(0.5,.1)   #Distribution so the resources splits unnequally 
    r=0.
    while (r<=0)
     r=rand(dNormATP)   
    end
    #r=.5
    b.ATP=b.ATP*(.95)  ##Cost of Division
    an=b.ATP
    b.ATP=b.ATP*r
    
    if(b.phenotype==1)
        dpA=0  # number of plasmid received by daughter 

        ######segregate plasmids
        pA=b.plasmids[2]

        for i=1:pA
            rand()<0.5 ? dpA+=1 : 0
        end

        b.plasmids[2]=pA-dpA

        dpA==0 ? echo("Bacteria $lid received 0 plasmids through segregation") : 1
        dpA==pA ? echo("Bacteria $(b.id) was left with 0 plasmids through segregation") : 1

        ######  b.phenotype
        b_d0=Bacteria(lid,b.phenotype,b.level+1,0,0,b.id,an-b.ATP,false,[0,dpA],[],[]) 
    elseif(b.phenotype==-1)
        b_d0=Bacteria(lid,b.phenotype,b.level+1,0,0,b.id,an-b.ATP,false,[0,0],[],[])        
    elseif(b.phenotype==0) 
         b_d0=Bacteria(lid,b.phenotype,b.level+1,0,0,b.id,an-b.ATP,false,[0,b.plasmids[3]],[],[])           
    end
    b_d=deepcopy(b_d0)
    b_d0=Bacteria
    return b_d
end
    








end #end of module
