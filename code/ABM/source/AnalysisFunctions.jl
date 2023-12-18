@everywhere module AnalysisFunctions
using Distributed
using Printf


@everywhere function alivebact_indexs(this_bst0)
    survivors0=[]
    for (ib,b) in enumerate(this_bst0)
        (b.dead==Inf) ? push!(survivors0,ib) : 0
    end
    return survivors0
end

@everywhere function get_index_Id_dict(this_bst0)
    tis=Dict()
    for (bi,b) in enumerate(this_bst0)
        tis[b.id]=bi
    end
    return tis
end

@everywhere function get_timePlamidsAliveBact(time_mins0,this_bst0)
    maxplasmids=get_max_plasmidsEver(this_bst0)
    pdict=Dict{Int,Int}()
    
    for t=0:maxplasmids  pdict[t]=0 end
    #pdict=Dict(sort(collect(pdict))) #does nothing to the result
    time_Pcounts0=[deepcopy(pdict) for x in time_mins0]
    ididx_dict=get_index_Id_dict(this_bst0)
    aliveI=alivebact_indexs(this_bst0)
    for bi in aliveI
        b=this_bst0[bi]
        tborn=b.born
        while(tborn>0)
            t0=b.born
            b.dead==Inf ? t1=Int(time_mins0[end]) : t1=Int(prev_born)
            this_pdist=b.plasmid_dist
            t_changes=[time_plasmids[1] for time_plasmids in this_pdist]
            p_changes=[time_plasmids[3] for time_plasmids in this_pdist]
            p=this_pdist[1][3]
            t0==0 ? t0=1 : 0
            
            for t=t0:t1
                ti=findall(x->x==t, t_changes)
                ti!=[] ? (p=p_changes[ti][1]) : 0
                time_Pcounts0[t][p]+=1
            end
            mother=b.mother
            tborn=b.born
            tborn>0 ? (newi=ididx_dict[mother];b=this_bst0[newi]) : 0
            prev_born=b.born
        end
    end
    return time_Pcounts0
end

@everywhere function get_Pseries_survivours_indexs(this_bst0)
    survivors=[]
    for (ib,b) in enumerate(this_bst0)
        (b.born==0) & (b.dead==Inf) ? push!(survivors,ib) : 0
    end
    return survivors
end

@everywhere function get_Pseries_killed_indexs(this_bst0,tKill1,tKill2)
    killed=[]
    for (ib,b) in enumerate(this_bst0)
        (b.born==0) & (b.dead<=tKill2) & (b.dead>=tKill1) ? push!(killed,ib) : 0
    end
    return killed
end


@everywhere function get_thisBac_Pseries(this_bst0,index)
    ts=[]
    ps=[]
    for tps in this_bst0[index].plasmid_dist
        push!(ts,tps[1])
        push!(ps,tps[3])
    end
    return ts,ps
end




@everywhere function get_max_plasmidsEver(this_bst0)
    maxplasmids0=0
    for b in this_bst0
        for time_plasmids in b.plasmid_dist
            maxplasmids0< time_plasmids[2] ? maxplasmids0=time_plasmids[2] : 0
        end
    end
    return maxplasmids0
end



@everywhere function get_timePlamids(time_mins0,this_bst0)
    maxplasmids=get_max_plasmidsEver(this_bst0)
    pdict=Dict{Int,Int}()
    #sort(collect(pdict), by = x->x[1])
    for t=0:maxplasmids  pdict[t]=0 end
    #pdict=Dict(sort(collect(pdict))) #does nothing to the result
    
    time_Pcounts0=[deepcopy(pdict) for x in time_mins0]

    for b in this_bst0
        t0=b.born
        
        b.dead==Inf ? t1=Int(max(time_mins0[end],b.plasmid_dist[end][1])) : t1=Int(max(b.dead,b.plasmid_dist[end][1]))
        this_pdist=b.plasmid_dist
        t_changes=[time_plasmids[1] for time_plasmids in this_pdist]
        p_changes=[time_plasmids[3] for time_plasmids in this_pdist]
        p=this_pdist[1][3]
        t0==0 ? t0=1 : 0
        for t=t0:t1
            ti=findall(x->x==t, t_changes)
            ti!=[] ? (p=p_changes[ti][1]) : 0
            #print(time_Pcounts0[t][p])
            time_Pcounts0[t][p]+=1
            #println("\t$t-$p--->$(time_Pcounts0[t][p])")
            #xxxccc
        end
    end

    return time_Pcounts0
end


@everywhere function binarize_timePcounts(time_Pcounts0,time_mins0,bin_edges0,normalize=false)
    time_Pcounts_bins0=zeros(length(bin_edges0),length(time_mins0))
    
    #pdict=Dict(sort(collect(pdict)))
    for (t,this_timedict) in enumerate(time_Pcounts0)
        tkeys=sort(collect(keys(this_timedict)))
        for thisK in tkeys 
            total=sum(values(this_timedict))
            prev_edge=bin_edges0[1]
            this_nplasmid_counts=this_timedict[thisK]
            for (ib,bin_number) in enumerate(bin_edges0)
                
                normalize ? tosum=this_nplasmid_counts/total : tosum=this_nplasmid_counts

                (thisK>prev_edge)&(thisK<=bin_number) ? time_Pcounts_bins0[ib,t]+=tosum : 0
                
                #time_Pcounts_bins0[ib,t]==0 ? time_Pcounts_bins0[ib,t]=0 : 0 
                prev_edge=bin_number
            end   
        end
    end
    return time_Pcounts_bins0
end

@everywhere function get_timePlamidsatDivisions(time_mins0,this_bst0)
    maxplasmids=get_max_plasmidsEver(this_bst0)
    pdict=Dict{Int,Int}()
    #sort(collect(pdict), by = x->x[1])
    for t=0:maxplasmids  pdict[t]=0 end
    #pdict=Dict(sort(collect(pdict))) #does nothing to the result
    
    time_Pcounts0=[deepcopy(pdict) for x in time_mins0]

    for b in this_bst0
        t0=b.born
        
        b.dead==Inf ? t1=Int(max(time_mins0[end],b.plasmid_dist[end][1])) : t1=Int(max(b.dead,b.plasmid_dist[end][1]))
        this_pdist=b.plasmid_dist
        this_divisions=b.divisions
        t_changes=[time_plasmids[1] for time_plasmids in this_pdist]
        p_changes=[time_plasmids[3] for time_plasmids in this_pdist]
        p=this_pdist[1][3]
        t0==0 ? t0=1 : 0
        for t=t0:t1
            
            ti=findall(x->x==t, t_changes)
            ti!=[] ? (p=p_changes[ti][1]) : 0
            isDivTime=t in this_divisions
            #print(time_Pcounts0[t][p])
            isDivTime ? time_Pcounts0[t][p]+=1 : 0
            #println("\t$t-$p--->$(time_Pcounts0[t][p])")
            #xxxccc
        end
    end

    return time_Pcounts0
end

@everywhere function get_meanPlasmidsTime(time_Pcounts0)
    time_Pmean=[]
    for x1 in time_Pcounts0
            #x1=sort(collect(x1), by = x->x[1])
            this_keys=keys(x1)
            this_keys=sort(collect(this_keys))
            
            real_counts=[0]
            for this_k in this_keys
                for ntimes=1:x1[this_k]
                    push!(real_counts,this_k)
                end
            end
            
            this_mean=mean(real_counts)
            this_mean >0 ? (real_counts=real_counts[real_counts.>0];this_mean=mean(real_counts)) : 0
            push!(time_Pmean,this_mean)
        end
    time_Pmean
end


@everywhere function describeBacteria(b::Bacteria)
    strLevel=""
    for i=1:b.level
        strLevel=string(strLevel,' ')
    end
    print("$strLevel ID: $(b.id), Phenotype $(b.phenotype) ATP: $(b.ATP) Level: $(b.level) Born: $(b.born) Dead: $(b.dead) Mother: $(b.mother)  isDead:$(b.isDead) plasmids:$(b.plasmids) divisions:$(b.divisions)\n")
    for k=1:length(b.plasmid_dist)
        println("$(b.plasmid_dist[k][1]):$(b.plasmid_dist[k]), ") 
    end
    println()
end


@everywhere function describeCellsOrdered(bs)
    tids=[]
    for b in bs
        push!(tids,b.id)
    end
    tids=sort(tids)
    for this_id in tids
        for b in bs
            b.id==this_id ? (describeBacteria(b); break ) : 0
        end
    end
end


@everywhere function normalize_FullTimeSeries(time_Pcounts,time_mins,bin_edges)
    time_Pcounts_bins0=zeros(length(bin_edges),length(time_mins))

    for (t,this_counts) in enumerate(time_Pcounts)
        for (ib,bin_number) in enumerate(bin_edges)
            nps=length(this_counts)
            tosum=1/nps
            #tosum=1
            for p in this_counts
                p==bin_number ? time_Pcounts_bins0[ib,t]+=tosum : 0
            end
            #time_Pcounts_bins0[ib,t]==0 ? time_Pcounts_bins0[ib,t]=0 : 0 
        end    
    end
    return time_Pcounts_bins0
end    


@everywhere function get_Nbacteria_Time(tbs,timev)
    nbact=zeros(length(timev))
    
    for b in tbs
        thisb=b.born
        thisb==0 ? thisb=1 : 0 
        for ti=thisb:timev[end]
           nbact[ti]+=1 
        end
    end
    return nbact
end




end  #end of module