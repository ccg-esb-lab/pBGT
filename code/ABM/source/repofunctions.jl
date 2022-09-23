function plotTimeSeriesAll(ts, PCNseries)
    plot(0, 0, alpha=0.)
    for j in collect(1:length(PCNseries)) #
        ps=PCNseries[j,:][1]
        if ps[1]>0
           lcolor="black"
            acolor=1.
            lwidth=4
        else
           lcolor="grey"
            acolor=0.5
            lwidth=1
        end
        
        xjiggle=rand()/10.
        yjiggle=rand()/5.
        
        
        scatter!(ts[ps.>-1][1:1]+xjiggle, ps[ps.>-1][1:1]+yjiggle, s=:solid, 
            markershape=:circle, color=lcolor, 
            alpha=acolor, legend=:none, markersize=2, linewidth=lwidth)
        plot!(ts[ps.>-1]+xjiggle, ps[ps.>-1]+yjiggle, alpha=acolor, linetype=:steppre, 
            color=lcolor, legend=:none,
            xlabel = "Time",ylabel="PCN",
            ylims = (0,40), grid=false, linewidth=lwidth)
    end
        savefig("$(pathRun)sim_mu$(mu)_T$(T)/timeseries_all_sigma$(floor(sigma*1000))e-3.pdf") 
end


function fitTimeSeriesDistribution(ti, PCNseries)

    bins=collect(0.:40.)

        PCN_t=Float64[]
        for j in collect(1:length(PCNseries))
            #print(PCNseries[j])
            push!(PCN_t,float(PCNseries[j][ti]))
        end

        x=PCN_t[PCN_t.>-1]
        nfit=fit(Normal, (x))
        #print(nfit)
    return nfit
end

function plotTimeSeriesDistribution(time_mins, PCNseries)

    bins=collect(0.:40.)
    #t_dist=[1,60,120]
    t_dist=time_mins

    for ti in t_dist
        PCN_t=Float64[]
        for j in collect(1:length(PCNseries))
            #print(PCNseries[j])
            push!(PCN_t,float(PCNseries[j][ti]))
        end

        #x=PCN_t[PCN_t.>-1]
        #nfit=fit(Normal, (x))
        #print(nfit)
        
        histogram(PCN_t[PCN_t.>-1], bins=bins, ylabel="PCN", xlabel="Number of cells", color="black", alpha=0.2, legend=:none, grid=false)
        
        #plot(nfit, fill=(0, .5,:orange))
        
        #savefig("$(pathRun)sim_mu$(mu)_T$(T)/dist_sigma$(floor(sigma*1000))e-3_t$(ti).pdf") # save the most recent fig as fn
        println("Exporting $(pathRun)sim_mu$(mu)_T$(T)/dist_sigma$(floor(sigma*1000))e-3_t$(ti).pdf") # save the most recent fig as fn
        
    end
    
end

#plotTimeSeriesDistribution([1, T], this_PCNseries)

@everywhere function getSurvivors(this_PCNseries, pcn_survs)
    

    this_PCNfracsurv=Float64[]
    for pcn_surv in pcn_survs

        num_survivors=0
        for j in collect(1:length(this_PCNseries))
            if this_PCNseries[j][end]>=pcn_surv
                num_survivors=num_survivors+1
            end
        end
        frac_survivors=num_survivors/length(this_PCNseries)
        #println("$(pcn_surv) Survivors: $(num_survivors)/$(length(this_PCNseries)) -> $(frac_survivors)")
        push!(this_PCNfracsurv, frac_survivors)
    end
    return this_PCNfracsurv 
end

