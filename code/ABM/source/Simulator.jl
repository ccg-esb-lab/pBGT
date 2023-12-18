@everywhere module Simulator
using Distributed

@everywhere function sampleBs(bk::Array{Bacteria},N::Int,numNews::Int)
    
    bs=Bacteria[]      
    alive_indexs=[]  
    #bi=Bacteria
    for (idx,b) in enumerate(bk)
        b.isDead==false ? push!(alive_indexs,idx) : 0
    end
    numAlive=length(alive_indexs)
    i=1
    while((i<=numAlive)&(i<=numNews))
    #for i=1:numNews
        ni=round(Int,rand()*numAlive)
        ni==0 ? ni=1 : 0
        ni>numAlive ? ni=numAlive : 0
        this_index=alive_indexs[ni]
        bi=deepcopy(bk[this_index])
        bi.level=0
        bi.born=0
        bi.ATP=bi.ATP*.25
        bi.dead=Inf
        bi.mother=bi.id  #Keeep track of previous day
        bi.id=i
        bi.divisions=[]
        bi.plasmid_dist=[]
        push!(bi.plasmid_dist,[0;bi.plasmids])
        push!(bs,bi)
        i+=1
        end
    
    return bs  
end





@everywhere function getInitialPopulation(iniCells,mu,phenotype,variation=true)
    b0s= Bacteria[]
    if(mu==1.)
        phenotype=0 
    end
    if(mu==0.)
        phenotype=-1
    end
    
    for i=1:iniCells
        plasmids=Int[]
        if(phenotype==1)
            pR=Int(round(rand(p_Ndist)))
            pR<0 ? pR=0 : 0
            
            plasmids=[pR,pR]
        else
            phenotype==-1 ? plasmids=[0,0] :  0
            phenotype==0 ? (plasmids=[1,1]) :  0 
        end
        
        if(variation)
            cATPi=rand(dNormATP)
        else
             cATPi=cATP*0.5
        end
         b0 = Bacteria(i,phenotype,0,0,Inf,0,cATPi,false,plasmids,[],[]) 
        #plasmids_dynamics(b0) 
        push!(b0.plasmid_dist,[0;b0.plasmids]) 
        push!(b0s,b0)
    end
    return b0s
end;

@everywhere function simulate(t::Int,bs::Array{Bacteria}, R::Float64, aA::Float64,antibiotic_action,lid::Int)
   verbose==true ? tic() : 0
    R<1e-7 ? (return R,aA,lid) : 0
    gATP=0.  # ATP gained
    
    shuffle!(bs)
    
    for i=1:length(bs)   

        if(bs[i].isDead==true)
            continue
        end

        uR,gATP,daA=grow(bs[i], R,aA)
            
        R-=uR
        aA-=daA
        
        bs[i].plasmids[2]<bs[i].plasmids[1] ? plasmids_dynamics(bs[i])  : 0
        
        bs[i].plasmids==bs[i].plasmid_dist[end][2:end]  ? 0 : (push!(bs[i].plasmid_dist,[t;bs[i].plasmids])  ) 
        
            
        if(bs[i].ATP>cATP*antibiotic_action)
            bs[i].isDead = antibiotic_dynamics(bs[i],aA)
            (bs[i].isDead==true) ? (bs[i].dead=t; echo("Bacteria $(i) dead at $(t) "); continue;) : 0
        end
            
                
        if( bs[i].ATP>cATP)
            
           # bs[i].isDead = antibiotic_dynamics(bs[i],aA)
           #(bs[i].isDead==true) ? (bs[i].dead=t; echo("Bacteria $(i) dead at $(t) "); continue;) : 0
            lid+=1  
            
            b0=divide(bs[i],lid)
            b0.born=t
            b0.dead=Inf
        
            pR=round(rand(p_Ndist))
            pR<0 ? (pR=0;echo("pR 0!")) : 0
            bs[i].plasmids[1]=deepcopy(pR)

            pR0=round(rand(p_Ndist))
            pR0<0 ? (pR0=0; echo("pR0 0!") ) : 0
            b0.plasmids[1]=deepcopy(pR0)

            push!(b0.plasmid_dist,[t;b0.plasmids]) 
            push!(bs[i].divisions,t) 
            push!(bs,b0)
        end #divide
        
    end #bs

    verbose==true ? toc() : 0
    return R,aA,lid;
end

@everywhere function simulate_nodaughters(t::Int,bs::Array{Bacteria}, R::Float64, aA::Float64,antibiotic_action,lid::Int)
   verbose==true ? tic() : 0
   
    gATP=0.  # ATP gained
    
    shuffle!(bs)
    
    for i=1:length(bs)   

        if(bs[i].isDead==true)
            continue
        end

        uR,gATP,daA=grow(bs[i], R,aA)
            
        R-=uR
        aA-=daA
        
            
        
        bs[i].plasmids[2]<bs[i].plasmids[1] ? plasmids_dynamics(bs[i])  : 0
        
        
        bs[i].plasmids==bs[i].plasmid_dist[end][2:end]  ? 0 : (push!(bs[i].plasmid_dist,[t;bs[i].plasmids])  ) 
        #println(maxgATP,gATP)    
            
        if(bs[i].ATP>cATP*antibiotic_action)
            bs[i].isDead = antibiotic_dynamics(bs[i],aA)
            (bs[i].isDead==true) ? (bs[i].dead=t; echo("Bacteria $(i) dead at $(t) "); continue;) : 0
        end
            
                
        if( bs[i].ATP>cATP)
            
           # bs[i].isDead = antibiotic_dynamics(bs[i],aA)
           #(bs[i].isDead==true) ? (bs[i].dead=t; echo("Bacteria $(i) dead at $(t) "); continue;) : 0
            lid+=1  
            
            b0=divide(bs[i],lid)
            b0.born=t
            b0.dead=Inf
        
            pR=round(rand(p_Ndist))
            pR<0 ? (pR=0;echo("pR 0!")) : 0
            bs[i].plasmids[1]=deepcopy(pR)

            pR0=round(rand(p_Ndist))
            pR0<0 ? (pR0=0; echo("pR0 0!") ) : 0
            b0.plasmids[1]=deepcopy(pR0)

            push!(b0.plasmid_dist,[t;b0.plasmids]) 
            push!(bs[i].divisions,t) 
            #push!(bs,b0)
        end #divide
        
    end #bs

    verbose==true ? toc() : 0
    return R,aA,lid;
end;






end # end of module
