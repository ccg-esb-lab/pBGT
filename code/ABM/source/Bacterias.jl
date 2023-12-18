#using Distributed
@everywhere module Bacterias
using Distributed


@everywhere mutable struct Bacteria
        id::Int32
        phenotype::Int8
        level::Int16
        born::Int32
        dead::Float32
        mother::Int32
        ATP::Float32
        isDead::Bool    
        plasmids::Array{Int16}     #contains  [real # ,current #]
        plasmid_dist::Array{Array{Int32}}
        divisions::Array{Int32}
    end

@everywhere function echo(s::AbstractString)
        verbose==true ? println(s) : 0
    end



@everywhere function chahca()
    println(rep_delay)
    end



end #end of module