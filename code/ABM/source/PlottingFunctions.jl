@everywhere module PlottingFunctions
using Distributed




@everywhere function get_Graph_data(this_bs::Array{Bacteria})
    all_p=[]
    rootsC=[]
    max_lvl=0
    for b in this_bs
        max_lvl<b.level ? max_lvl=b.level : 0
        b.level==0 ? push!(rootsC,"$(b.id).$(b.divisions[1])") : 0
        for bsit in b.plasmid_dist
            push!(all_p,bsit[3])
        end
    end
    max_plasmids=maximum(all_p)
    min_plasmids=minimum(all_p)
    n_cells=length(this_bs)
    println("Max PCN: ",max_plasmids)
    println("Min PCN: ",min_plasmids)
    println("Roots: ",rootsC)
    println("N cells: ", n_cells)
    println("Max lvls: ", max_lvl)
    
    return rootsC,max_plasmids,min_plasmids,max_lvl,n_cells
end

@everywhere function get_node_sizes(Gi)
    list_sizesA=[]
    for node in Gi.nodes(data=true)
        if haskey(node[2],"n_size")
    
            this_sz=node[2]["n_size"]
        else
            this_sz=0
        end
        push!(list_sizesA,this_sz)
    end

    return list_sizesA
end

@everywhere function make_Graph_NodesAtDivision(this_bs::Array{Bacteria},time_mins0)

    Gt_C=ntx.DiGraph()   #Initialize directed Graph
    ########This section makes the tracks of each cell
    for b in this_bs
        id=b.id
        born=b.born
        name0="$id.$born"

        this_p=0
        last_name=name0
        this_divs_plus=[b.born;b.divisions;time_mins0[end]]
        #this_divs_plus=[b.divisions;time_mins[end]]
        #println(last_name,this_divs_plus)
        for t=b.born:time_mins0[end]
            for bsit in b.plasmid_dist
                if bsit[1]==t 
                    this_p=bsit[3]
                end
            end
            if t in this_divs_plus # && last_name!=name0
                name="$id.$t"
                #println(name," ",last_name," ",born)
                Gt_C.add_node(name,n_plasmid=this_p)
                Gt_C.add_edge(last_name,name)
                last_name=name
            end
        #println(last_name,name0)   
        end
    end
    
    ##### This section eliminates not division nodes  #########
    to_remove=[]
    to_link=[]
    for b in this_bs
        name="$(b.id).$(b.born)"
        for node in Gt_C.nodes
            if node==name
                push!(to_remove,node)
                for ng in Gt_C.neighbors(node)
                 #   println(node,"*",ng)
                    if node!=ng
                        push!(to_link,ng)
                    end
                end
            end
        end
    end

    #println(to_remove)
    for rnode in to_remove
        Gt_C.remove_node(rnode)
    end

    #println("Single tracks\n",Gt_C.nodes(data=true))
    #println(to_link)
    
    ############### This section
    
    for b in this_bs
        id=b.id
        born=b.born
        mother=b.mother
        bsit=[]
        for bm in this_bs
            if bm.id==mother
               for k=1:length(bm.plasmid_dist)
                    bsit=bm.plasmid_dist[k]
                    this_pl=bsit[3]
                    if(bsit[1]>=born)
                        break;break
                    end
               end
            end
        end

        for link_node in to_link
            link_id=split(link_node,".")[1]
            link_id=parse(Int,link_id)
            #println(link_id)
            if(link_id==id)
                name=link_node
               # println(name)
                mother_name="$mother.$born"
                Gt_C.add_edge(mother_name,name)
            end
        end


        #Gt_C.node[name]["plasmids"]=this_pl
    end

    Gt_C.add_node("0.0",n_plasmid=0)    # This line roots the inital nodes to a fake node to improve visualization
    return Gt_C
end

@everywhere function get_node_colors(G,gmax)
    green_lvl= colormap("greens",gmax+2)
    list_colorsC=[]
    list_plasmidsC=[]
    for node in G.nodes(data=true)
        #println(node)
        this_pl=node[2]["n_plasmid"]
        this_g='#'*Colors.hex(green_lvl[this_pl+1])
        push!(list_plasmidsC,this_pl)
        push!(list_colorsC,this_g)
    end

    return list_colorsC
end

@everywhere function make_Graph_All_and_Divisions(this_bs::Array{Bacteria},time_mins0)

    Gt_allC=ntx.DiGraph()
    ########This section makes the tracks of each cell
    for b in this_bs
        id=b.id
        born=b.born
        name0="$id-$born"

        this_p=0
        last_name=name0

        this_divs_plus=[b.born;b.divisions;time_mins0[end]]

        for t=b.born:time_mins0[end]
            for bsit in b.plasmid_dist
                if bsit[1]==t
                    this_p=bsit[3]
                end
            end
            if t in this_divs_plus
                name="$id-$t"
                Gt_allC.add_node(name,n_plasmid=this_p,n_size=20)
                Gt_allC.add_edge(last_name,name)
                last_name=name
            else
                name="$id-$t"
                Gt_allC.add_node(name,n_plasmid=this_p,n_size=20)
                Gt_allC.add_edge(last_name,name)
                last_name=name
            end
        end
    end

    for b in this_bs
        id=b.id
        born=b.born
        mother=b.mother
        bsit=[]
        for bm in this_bs
            if bm.id==mother
               for k=1:length(bm.plasmid_dist)
                    bsit=bm.plasmid_dist[k]
                    this_pl=bsit[3]
                    if(bsit[1]>=born)
                        break;break
                    end
               end
            end
        end

        mother_name="$mother-$born"
        name="$id-$born"
        Gt_allC.add_edge(mother_name,name)
    end
    Gt_allC.add_node("0-0",n_plasmid=0,n_size=20)


    ######################################################3
    for b in this_bs
        id=b.id
        born=b.born
        name0="$id-$born"

        this_p=0
        last_name=name0
        this_divs_plus=[b.born;b.divisions;time_mins0[end]]
        #this_divs_plus=[b.divisions;time_mins[end]]
        #println(last_name,this_divs_plus)
        for t=b.born:time_mins0[end]

            for bsit in b.plasmid_dist
                if bsit[1]==t 

                    this_p=bsit[3]
                end
            end
            if t in this_divs_plus # && last_name!=name0

                name="$id-$t"
                #println(name," ",last_name," ",born)
                Gt_allC.add_node(name,n_plasmid=this_p,n_size=100)
              #  Gt_allC.add_edge(name,last_name)
                last_name=name

            end
        #println(last_name,name0)   
        end
    end

    to_remove=[]
    to_link=[]
    for b in this_bs
        name="$(b.id)-$(b.born)"
        for node in Gt_allC.nodes
            if node==name
                push!(to_remove,node)
                for ng in Gt_allC.neighbors(node)
                 #   println(node,"*",ng)
                    if node!=ng
                        push!(to_link,ng)
                    end
                end
            end
        end
    end

    for rnode in to_remove
        #Gt_allC.remove_node(rnode)
        for node in Gt_allC.nodes(data=true)
            if node[1]==rnode
               p_d=node[2]
               this_p=p_d["n_plasmid"] 
               Gt_allC.add_node(node[1],n_plasmid=this_p,n_size=20) 
            end
        end
    end

    for b in this_bs
        id=b.id
        born=b.born
        mother=b.mother
        bsit=[]
        for bm in this_bs
            if bm.id==mother
               for k=1:length(bm.plasmid_dist)
                    bsit=bm.plasmid_dist[k]
                    this_pl=bsit[3]
                    if(bsit[1]>=born)
                        break;break
                    end
               end
            end
        end

        for link_node in to_link
            link_id=split(link_node,"-")[1]
            link_id=parse(Int,link_id)
            #println(link_id)
            if(link_id==id)
                name=link_node
               # println(name)
                mother_name="$mother-$born"
      #          Gt_allC.add_edge(name,mother_name)
            end
        end
    end

    Gt_allC.add_node("0-0",n_plasmid=0)

    return Gt_allC
end

@everywhere function get_node_colorsSharp(G,gmin,gmax)
    green_lvl= colormap("greens",gmax-gmin+3)
    list_colorsC=[]
    list_plasmidsC=[]
    for node in G.nodes(data=true)
        #println(node)
        this_pl=node[2]["n_plasmid"]
        #println(node,this_pl)
        mplasmid=this_pl-gmin+2
        mplasmid<1 ? mplasmid=1 : 0
        this_g='#'*Colors.hex(green_lvl[mplasmid])
        push!(list_plasmidsC,this_pl)
        push!(list_colorsC,this_g)
    end

    return list_colorsC
end

@everywhere function make_Graph_All(this_bs::Array{Bacteria},time_mins0)
    Gt_all=ntx.DiGraph()
    ########This section makes the tracks of each cell
    for b in this_bs
        id=b.id
        born=b.born
        name0="$id-$born"

        this_p=0
        last_name=name0

        this_divs_plus=[b.born;b.divisions;time_mins0[end]]
        
        maxt=Int(min(time_mins0[end],b.dead))
        for t=b.born:maxt
            for bsit in b.plasmid_dist
                if bsit[1]<=t
                    this_p=bsit[3]
                end
            end
            if t in this_divs_plus
                name="$id-$t"
                Gt_all.add_node(name,n_plasmid=this_p,n_size=80,isDead=b.isDead)
                Gt_all.add_edge(last_name,name)
                last_name=name

            else
                name="$id-$t"
                Gt_all.add_node(name,n_plasmid=this_p,n_size=20,isDead=b.isDead)
                Gt_all.add_edge(last_name,name)
                last_name=name

            end

        end
    end


    #println("Single tracks\n",Gt_all.nodes(data=true))

    for b in this_bs
        id=b.id
        born=b.born
        mother=b.mother
        bsit=[]
        for bm in this_bs
            if bm.id==mother
               for k=1:length(bm.plasmid_dist)
                    bsit=bm.plasmid_dist[k]
                    this_pl=bsit[3]
                    if(bsit[1]>=born)
                        break;break
                    end
               end
            end
        end


        mother_name="$mother-$born"
        name="$id-$born"
        Gt_all.add_edge(mother_name,name)

        #Gt_all.node[name]["plasmids"]=this_pl
    end

    # println("Joined tracks\n",Gt_all.nodes(data=true))
    Gt_all.add_node("0-0",n_plasmid=0,n_size=20)
    return Gt_all
end



##########################################





end