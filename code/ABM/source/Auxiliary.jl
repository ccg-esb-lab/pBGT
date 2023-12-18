@everywhere module Auxiliary
using Distributed


@everywhere function divide_equal(b::Bacteria,lid::Int)
    
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

        dpA=Int(round(pA/2))
        

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
@everywhere function simulate_equal(t::Int,bs::Array{Bacteria}, R::Float64, aA::Float64,antibiotic_action,lid::Int)
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
            
            b0=divide_equal(bs[i],lid)
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



@everywhere function linearize_counts(bst,time_minsx,bin_edges,norm)
    Ntime_Pcounts_bins=binarize_timePcounts(bst,time_minsx,bin_edges,norm)
    thisd=Ntime_Pcounts_bins[:,end]
    all=[]
    for (i,x) in enumerate(bin_edges)
        #print("$i-$x-$(thisd[i]), ")
       if(thisd[i]!=0) 
           xx=repeat([x],Int(thisd[i]))
            all=[all;xx]
        end
    end
    
    return all
end

@everywhere function test_dist(this_dist_test,aA)
    global p_Ndist=this_dist_test
    global mu

    #antibiotic_action=0.65

    println("Time\tTotal\tResource         .\tN alive Antibiotic")

    R0=1.0

    this_bst=getInitialPopulation(iniCells, mu, 1,starting_variation);
    pdist=get_timePlamids([1],this_bst)
    pcounts0=linearize_counts(pdist,[0],bin_edges,norm)
    
    lid=length(this_bst)

    for ti in time_mins
        if (ti%30==0)|(ti==1)
            tparam=0
            for b in this_bst
                b.isDead==false ? tparam+=1 : 0
            end
            print(ti/60,"\t",length(this_bst),"\t",R0,"\t",tparam,"\t",aA,"\r")
            flush(stdout)
        end
        R0,aA,lid=simulate(ti, this_bst, R0, aA,antibiotic_action,lid)
    end
    println()
    pdist=get_timePlamids(time_mins,this_bst)
    pcountsA=linearize_counts(pdist,time_mins,bin_edges,norm)
    

   return pcounts0,pcountsA
end    

@everywhere function get_mean_from_bincounts(pcounts,bins)
    meanV=[]
    
    for it in collect(1:size(pcounts)[2])
        time_pcoutns=pcounts[:,it]
        tmean=0
        cn=0
        for (ib,b) in enumerate(bins)
            tmean+= time_pcoutns[ib]*b
            time_pcoutns[ib]>0 ?  cn+=time_pcoutns[ib] : 0
        end
        tmean=round(tmean/cn,digits=2)
        push!(meanV,tmean)
    end
    
    return meanV
    
end

@everywhere function plotStreambinned(Ntime_Pcounts_bins0,Ndays_time_mins0,time_mins0,Ndays0,bin_edges0)
    ncolors=length(Ntime_Pcounts_bins0[:,1])

    tc0=palette(:BuGn,ncolors+1 )# )
    graycolor=palette(:grays,3)[2]
    tcp=[c for c in tc0[2:end]]
    tcp[1]=graycolor


    bin_edgesl=["0"]
    for (i,x) in enumerate(bin_edges0[2:end])
        nbs="$(string(bin_edges0[i]+1))-$(string(bin_edges0[i+1]))"
        push!(bin_edgesl,nbs)
    end
    bin_edgesl=reshape(bin_edgesl,1,length(bin_edgesl))


    newtickslabels=collect(0:Ndays0)
    newticks=newtickslabels*time_mins0[end]
    totals=[sum(Ntime_Pcounts_bins0[:,t]) for t in Ndays_time_mins0]
    maxBac=maximum(totals)
    p=plot(size=(1200,400))

    p=areaplot!(p,Ndays_time_mins0,transpose(Ntime_Pcounts_bins0),label=bin_edgesl,color_palette=tcp)

    plot!(p,xticks=(newticks,newtickslabels),xlabel="Time (season)",xlim=(0,Ndays_time_mins0[end]))
    plot!(p,ylim=(100,maxBac),ylabel="Number of bacteria (log)",
    #yaxis=(formatter=y->string(round( y / 10^6,digits=2))),
    #formatter = :scientific
    )
    #annotate!([(0, maxBac * 1.05, Plots.text(L"\times10^{6}", 11, :black, :center))])
    plot!(p,legendtitle="PCN",legend=:outerright,legendfontsize=12,legendtitlefontsize=16)
    return p
end


@everywhere function plotHeatbinned(Ntime_Pcounts_bins0,Ndays_time_mins0,time_mins0,Ndays0,bin_edges0,loged,wcbar)

    #tcp=cgrad([colorant"white",colorant"blue"])
    tcp=cgrad([       colorant"cyan",
            colorant"blue", 
            colorant"darkslateblue",
            colorant"indigo",
            colorant"purple",
            colorant"orangered",
            colorant"red",
            ]) 
    
    #tcp=palette(:cool,1000 )
    
    pmean=get_mean_from_bincounts(Ntime_Pcounts_bins0,bin_edges0) 
    Ntime_Pcounts_bins00=replace(Ntime_Pcounts_bins0, 0=>NaN)
    
    newticks=collect(0:30:time_mins0[end])
    
    p=plot(size=(1200,300))
    if(loged)
        p=heatmap(Ndays_time_mins0,bin_edges0,log10.(Ntime_Pcounts_bins00),colorbar=wcbar, 
        color=tcp,alpha=1,clim=(log10(0),log10(1)))
        plot!(p,cbartitle="Log density",legend=:bottomright)
    else
        p=heatmap(Ndays_time_mins0,bin_edges0,(Ntime_Pcounts_bins00),colorbar=wcbar, color=tcp,alpha=1,
        clim=(0.0,1)
        )
    plot!(p,cbartitle="Density",legend=:bottomright)
    end
    plot!(p,Ndays_time_mins0,pmean,color="black",linewidth=3,label="mean PCN")
    plot!(p,xticks=newticks,xlabel="Time (minutes)",xlim=(0,Ndays_time_mins0[end]))
    plot!(p,ylim=(-1,bin_edges0[end]),ylabel="PCN")

    return p
end



@everywhere function getInitialPopulation2(iniCells,tmu,phenotype,variation=true, sigma=5.0)
    
    b0s= Bacteria[]

    if(mu==1.)
        phenotype=0 
    end
    if(mu==0.)
        phenotype=-1
    end
    p_Ndist = Normal(tmu, sigma)
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

@everywhere function test_survival2(tmu, sigma, ic)
    
    
    muMIC = getMIC(mu)
    pop_mu = getInitialPopulation2(iniCells, tmu, 1, true, sigma)
    this_antibiotic = muMIC * ic
    n_alive = 0

    for b in pop_mu
        b.isDead = antibiotic_dynamics(b, this_antibiotic)
        if !b.isDead
            n_alive += 1
        end
    end

    return n_alive
end



end  # end module