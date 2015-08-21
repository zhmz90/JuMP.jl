#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

using Base.Meta
# multivarate "dictionary" used for collections of variables/constraints

abstract JuMPContainer{T,N}
abstract JuMPArray{T,N} <: JuMPContainer{T,N} # Q is true if all index sets are of the form 1:n

type IndexPair
    idxvar
    idxset
end
#= Generated on the fly
type JuMPArray{T}
    innerArray::Array{T,N}
    name::String
    indexsets
end
=#
type JuMPDict{T,N} <: JuMPContainer{T,N}
    tupledict::Dict{NTuple{N,Any},T}
end

type JuMPContainerData
    name
    indexsets
    indexexprs::Vector{IndexPair}
    condition
end

#JuMPDict{T,N}(name::String) =
#    JuMPDict{T,N}(Dict{NTuple{N},T}(), name)

Base.getindex(d::JuMPDict, t...) = d.tupledict[t]
Base.setindex!(d::JuMPDict, value, t...) = (d.tupledict[t] = value)

function Base.map{T,N}(f::Function, d::JuMPDict{T,N})
    ret = Base.return_types(f, @compat(Tuple{T}))
    R = (length(ret) == 1 ? ret[1] : Any)
    x = JuMPDict(Dict{NTuple{N,Any},R}())
    for (k,v) in d.tupledict
        x.tupledict[k] = f(v)
    end
    return x
end

Base.isempty(d::JuMPContainer) = isempty(_innercontainer(d))

# generate and instantiate a type which is indexed by the given index sets
# the following types of index sets are allowed:
# 0:K -- range with compile-time starting index
# S -- general iterable set
export @gendict
macro gendict(instancename,T,idxpairs,idxsets...)
    N = length(idxsets)
    allranges = all(s -> (isexpr(s,:(:)) && length(s.args) == 2), idxsets)
    if allranges
        truearray = all(idxsets) do rng
            isexpr(rng,:(:)) && # redundant, but might as well
            rng.args[1] == 1 && # starts from 1
            (length(rng.args) == 2 || rng.args[2] == 1) # steps of one
        end
        if truearray
            sizes = Expr(:tuple, [esc(rng.args[end]) for rng in idxsets]...)
            :($(esc(instancename)) = Array($T, $sizes))
        else
            typename = symbol(string("JuMPArray",gensym()))
            dictnames = Array(Symbol,N)
            # JuMPArray
            offset = Array(Int,N)
            for i in 1:N
                if isa(idxsets[i].args[1],Int)
                    offset[i] = 1 - idxsets[i].args[1]
                else
                    error("Currently only ranges with integer compile-time starting values are allowed as index sets. $(idxsets[i].args[1]) is not an integer in range $(idxsets[i]).")
                end
            end
            typecode = :(type $(typename){T} <: JuMPArray{T,$N}; innerArray::Array{T,$N}; end)
            getidxlhs = :(Base.getindex(d::$(typename)))
            setidxlhs = :(Base.setindex!(d::$(typename),val))
            getidxrhs = :(Base.getindex(d.innerArray))
            setidxrhs = :(Base.setindex!(d.innerArray,val))
            maplhs = :(Base.map(f::Function,d::$(typename)))
            maprhs = :($(typename)(map(f,d.innerArray)))
            wraplhs = :(JuMPContainer_from(d::$(typename),inner)) # helper function that wraps array into JuMPArray of similar type
            wraprhs = :($(typename)(inner))
            for i in 1:N
                varname = symbol(string("x",i))

                push!(getidxlhs.args,:($varname))
                push!(setidxlhs.args,:($varname))

                push!(getidxrhs.args,:(isa($varname, Int) ? $varname+$(offset[i]) : $varname ))
                push!(setidxrhs.args,:($varname+$(offset[i])))

            end

            badgetidxlhs = :(Base.getindex(d::$(typename),wrong...))
            badgetidxrhs = :(if isempty(d)
                                error("Wrong number of indices")
                            else
                                data = first(values(d)).m.varData[d]
                                error("Wrong number of indices for ",data.name, ", expected ",length(data.indexsets))
                            end)

            funcs = :($getidxlhs = $getidxrhs; $setidxlhs = $setidxrhs;
                      $maplhs = $maprhs; $badgetidxlhs = $badgetidxrhs;
                      $wraplhs = $wraprhs)
            geninstance = :($(esc(instancename)) = $(typename)(Array($T)))
            for i in 1:N
                push!(geninstance.args[2].args[2].args, :(length($(esc(idxsets[i])))))
            end
            eval(Expr(:toplevel, typecode))
            eval(Expr(:toplevel, funcs))

            return geninstance
        end

        #= TODO: Use this with code in JuMPArray.jl once Julia can make it efficient
        if all([length(s.args) == 3 for s in idxsets])
            # has step
            for i in 1:N
                if length(idxsets[i].args) == 2
                    push!(idxsets[i].args, idxsets[i].args[2])
                    idxsets[i].args[2] = 1
                end
            end
        end
        geninstance = :($(esc(instancename)) = JuMPArray(Array($T),$(quot(instancename)),$(esc(Expr(:tuple,idxsets...)))))
        for i in 1:N
            push!(geninstance.args[2].args[2].args, :(length($(esc(idxsets[i])))))
        end
        return geninstance
        =#
    else
        # JuMPDict
        return :(
            $(esc(instancename)) = JuMPDict{$T,$N}(Dict{NTuple{$N},$T}())
        )
    end
end

# duck typing approach -- if eltype(innerArray) doesn't support accessor, will fail
for accessor in (:getDual, :getLower, :getUpper)
    @eval $accessor(x::Union(JuMPContainer,Array)) = map($accessor,x)
end

_similar(x::Array) = Array(Float64,size(x))
_similar{T}(x::Dict{T}) = Dict{T,Float64}()

_innercontainer(x::JuMPArray) = x.innerArray
_innercontainer(x::JuMPDict)  = x.tupledict

function _getValueInner(x)
    vars = _innercontainer(x)
    vals = _similar(vars)
    if !isempty(vals)
        m = first(values(x)).m
        data = m.varData[x]
        warnedyet = false
        for I in eachindex(vars)
            tmp = _getValue(vars[I])
            if isnan(tmp) && !warnedyet
                warn("Variable value not defined for entry of $(data.name). Check that the model was properly solved.")
                warnedyet = true
            end
            vals[I] = tmp
        end
    end
    vals
end

JuMPContainer_from(x::JuMPDict,inner) = JuMPDict(inner)

function getValue(x::JuMPContainer)
    getvalue_warn(x)
    JuMPContainer_from(x,_getValueInner(x))
end

# delegate zero-argument functions
for f in (:(Base.endof), :(Base.ndims), :(Base.length), :(Base.abs), :(Base.start))
    @eval $f(x::JuMPArray) = $f(x.innerArray)
end

Base.first(x::JuMPDict)  =  first(x.tupledict)
Base.length(x::JuMPDict) = length(x.tupledict)

Base.ndims{T,N}(x::JuMPDict{T,N}) = N
Base.abs(x::JuMPDict) = map(abs, x)
# delegate one-argument functions
Base.size(x::JuMPArray)   = size(x.innerArray)
Base.size(x::JuMPArray,k) = size(x.innerArray,k)
Base.issym(x::JuMPArray) = issym(x.innerArray)

function Base.start(x::JuMPContainer)
    warn("Iteration over JuMP containers is deprecated. Use keys(d) and values(d) instead")
    start(x.tupledict)
end

function Base.next(x::JuMPArray,k)
    var, gidx = next(x.innerArray, k)
    keys = _next_index(x,k)
    tuple(keys..., x[keys...]), gidx
end

function _next_index{T,N}(x::JuMPArray{T,N}, k)
    @assert !isempty(x)
    m = first(x.innerArray).m
    data = m.varData[x]
    ntuple(N) do index
        cprod = 1
        for i in 1:(index-1)
            cprod *= length(data.indexsets[i])
        end
        locidxset = data.indexsets[index]
        idx = Compat.ceil(Int, mod1(k, cprod*length(locidxset)) / cprod)
        locidxset[idx]
    end
end

function Base.next(x::JuMPDict,k)
    ((idx,var),gidx) = next(x.tupledict,k)
    return (tuple(idx..., var), gidx)
end

Base.done(x::JuMPArray,k) = done(x.innerArray,k)
Base.done(x::JuMPDict,k)  = done(x.tupledict,k)

Base.eltype{T}(x::JuMPContainer{T}) = T

Base.full(x::JuMPContainer) = x

# keys/vals iterations for JuMPContainers
Base.keys(d::JuMPDict)    = keys(d.tupledict)
Base.values(d::JuMPDict)  = values(d.tupledict)

Base.keys(d::JuMPArray)   = KeyIterator(d)
Base.values(d::JuMPArray) = ValueIterator(d.innerArray)

# TODO: add keys/values for Array
Base.values(d::Array) = d

# Wrapper type so that you can't access the values directly
type ValueIterator{T,N}
    x::Array{T,N}
end
Base.start(it::ValueIterator)   =  start(it.x)
Base.next(it::ValueIterator, k) =   next(it.x, k)
Base.done(it::ValueIterator, k) =   done(it.x, k)
Base.length(it::ValueIterator)  = length(it.x)

type KeyIterator{JA<:JuMPArray}
    x::JA
end
Base.start(it::KeyIterator)   =  start(it.x.innerArray)
Base.next(it::KeyIterator, k) =  _next_index(it.x, k), next(it.x.innerArray, k)[2]
Base.done(it::KeyIterator, k) =   done(it.x.innerArray, k)
Base.length(it::KeyIterator)  = length(it.x.innerArray)
