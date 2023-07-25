
"""
    DBCM{T,N} <: AbstractMaxEntropyModel

Type definition for the Directed Binary Configuration Model (DBCM) model.
"""
mutable struct DBCM{T,N} <: AbstractMaxEntropyModel where {T<:Union{Graphs.AbstractGraph, Nothing}, N<:Real}
    "Graph type, can be any subtype of AbstractGraph, but will be converted to SimpleDiGraph for the computation" # can also be empty
    const G::T 
    "Vector holding all maximum likelihood parameters for reduced model (α ; β)"
    const θᵣ::Vector{N}
    "Exponentiated maximum likelihood parameters for reduced model ( xᵢ = exp(-αᵢ) ) linked with out-degree"
    const xᵣ::Vector{N}
    "Exponentiated maximum likelihood parameters for reduced model ( yᵢ = exp(-βᵢ) ) linked with in-degree"
    const yᵣ::Vector{N}
    "Outdegree sequence of the graph" # evaluate usefulness of this field later on
    const d_out::Vector{Int}
    "Indegree sequence of the graph" # evaluate usefulness of this field later on
    const d_in::Vector{Int}
    "Reduced outdegree sequence of the graph"
    const dᵣ_out::Vector{Int}
    "Reduced indegree sequence of the graph"
    const dᵣ_in::Vector{Int}
    "Indices of non-zero elements in the reduced outdegree sequence"
    const dᵣ_out_nz::Vector{Int}
    "Indices of non-zero elements in the reduced indegree sequence"
    const dᵣ_in_nz::Vector{Int}
    "Frequency of each (outdegree, indegree) pair in the graph"
    const f::Vector{Int}
    "Indices to reconstruct the degree sequence from the reduced degree sequence"
    const d_ind::Vector{Int}
    "Indices to reconstruct the reduced degree sequence from the degree sequence"
    const dᵣ_ind::Vector{Int}
    "Expected adjacency matrix" # not always computed/required
    Ĝ::Union{Nothing, Matrix{N}}
    "Variance of the expected adjacency matrix" # not always computed/required
    σ::Union{Nothing, Matrix{N}}
    "Status indicators: parameters computed, expected adjacency matrix computed, variance computed, etc."
    const status::Dict{Symbol, Real}
    "Function used to computed the log-likelihood of the (reduced) model"
    fun::Union{Nothing, Function}
end


Base.show(io::IO, m::DBCM{T,N}) where {T,N} = print(io, """UBCM{$(T), $(N)} ($(m.status[:d]) vertices, $(m.status[:d_unique]) unique degree pairs, $(@sprintf("%.2f", m.status[:cᵣ])) compression ratio)""")

"""Return the reduced number of nodes in the UBCM network"""
Base.length(m::DBCM) = length(m.dᵣ)


"""
    DBCM(G::T; precision::N=Float64, kwargs...) where {T<:Graphs.AbstractGraph, N<:Real}
    DBCM(;d_out::Vector{T}, d_in::Vector{T}, precision::Type{<:AbstractFloat}=Float64, kwargs...)

Constructor function for the `DBCM` type. 
    
By default and dependng on the graph type `T`, the definition of in- and outdegree from ``Graphs.jl`` is applied. 
If you want to use a different definition of degrees, you can pass vectors of degrees sequences as keyword arguments (`d_out`, `d_in`).
If you want to generate a model directly from degree sequences without an underlying graph , you can simply pass the degree sequences as arguments (`d_out`, `d_in`).
If you want to work from an adjacency matrix, or edge list, you can use the graph constructors from the ``JuliaGraphs`` ecosystem.

# Examples     
```jldoctest
# generating a model from a graph


# generating a model directly from a degree sequence


# generating a model directly from a degree sequence with a different precision


# generating a model from an adjacency matrix


# generating a model from an edge list


```

See also [`Graphs.outdegree`](@ref), [`Graphs.indegree`](@ref), [`SimpleWeightedGraphs.outdegree`](@ref), [`SimpleWeightedGraphs.indegree`](@ref).
"""
function DBCM(G::T; d_out::Vector=Graphs.outdegree(G), 
                    d_in::Vector=Graphs.indegree(G), 
                    precision::Type{N}=Float64, 
                    kwargs...) where {T,N<:AbstractFloat}
    T <: Union{Graphs.AbstractGraph, Nothing} ? nothing : throw(TypeError("G must be a subtype of AbstractGraph or Nothing"))
    length(d_out) == length(d_in) ? nothing : throw(DimensionMismatch("The outdegree and indegree sequences must have the same length"))
    # coherence checks
    if T <: Graphs.AbstractGraph # Graph specific checks
        if !Graphs.is_directed(G)
            @warn "The graph is undirected, while the DBCM model is directed, the in- and out-degree will be the same"
        end

        if T <: SimpleWeightedGraphs.AbstractSimpleWeightedGraph
            @warn "The graph is weighted, while DBCM model is unweighted, the weight information will be lost"
        end

        Graphs.nv(G) == 0 ? throw(ArgumentError("The graph is empty")) : nothing
        Graphs.nv(G) == 1 ? throw(ArgumentError("The graph has only one vertex")) : nothing
        Graphs.nv(G) != length(d_out) ? throw(DimensionMismatch("The number of vertices in the graph ($(Graphs.nv(G))) and the length of the degree sequence ($(length(d))) do not match")) : nothing
    end
    # coherence checks specific to the degree sequences
    length(d_out) == 0 ? throw(ArgumentError("The degree sequences are empty")) : nothing
    length(d_out) == 1 ? throw(ArgumentError("The degree sequences only contain a single node")) : nothing
    maximum(d_out) >= length(d_out) ? throw(DomainError("The maximum outdegree in the graph is greater or equal to the number of vertices, this is not allowed")) : nothing
    maximum(d_in)  >= length(d_in)  ? throw(DomainError("The maximum indegree in the graph is greater or equal to the number of vertices, this is not allowed")) : nothing

    # field generation
    dᵣ, d_ind , dᵣ_ind, f = np_unique_clone(collect(zip(d_out, d_in)), sorted=true)
    dᵣ_out = [d[1] for d in dᵣ]
    dᵣ_in =  [d[2] for d in dᵣ]
    dᵣ_out_nz = findall(!iszero, dᵣ_out)
    dᵣ_in_nz  = findall(!iszero, dᵣ_in)
    Θᵣ = Vector{precision}(undef, 2*length(dᵣ))
    xᵣ = Vector{precision}(undef, length(dᵣ))
    yᵣ = Vector{precision}(undef, length(dᵣ))
    status = Dict{Symbol, Real}(:params_computed=>false,            # are the parameters computed?
                                :G_computed=>false,                 # is the expected adjacency matrix computed and stored?
                                :σ_computed=>false,                 # is the standard deviation computed and stored?
                                :cᵣ => length(dᵣ)/length(d_out),    # compression ratio of the reduced model
                                :d_unique => length(dᵣ),            # number of unique (outdegree, indegree) pairs in the reduced model
                                :d => length(d_out)                 # number of vertices in the original graph 
                )
    
    return DBCM{T,precision}(G, Θᵣ, xᵣ, yᵣ, d_out, d_in, dᵣ_out, dᵣ_in, dᵣ_out_nz, dᵣ_in_nz, f, d_ind, dᵣ_ind, nothing, nothing, status, nothing)
end

DBCM(; d_out::Vector{T}, d_in::Vector{T}, precision::Type{N}=Float64, kwargs...) where {T<:Signed, N<:AbstractFloat} = DBCM(nothing; d_out=d_out, d_in=d_in, precision=precision, kwargs...)


"""
    L_DBCM_reduced(θ::Vector, k_out::Vector, k_in::Vector, F::Vector, nz_out::Vector, nz_in::Vector, n::Int=length(k_out))

Compute the log-likelihood of the reduced DBCM model using the exponential formulation in order to maintain convexity.

The arguments of the function are:
    - `θ`: the maximum likelihood parameters of the model ([α; β])
    - `k_out`: the reduced outdegree sequence
    - `k_in`: the reduced indegree sequence
    - `F`: the frequency of each pair in the degree sequence
    - `nz_out`: the indices of non-zero elements in the reduced outdegree sequence
    - `nz_in`: the indices of non-zero elements in the reduced indegree sequence
    - `n`: the number of nodes in the reduced model

The function returns the log-likelihood of the reduced model. For the optimisation, this function will be used to
generate an anonymous function associated with a specific model.

# Examples
```jldoctest
# Generic use:
julia> k_out  = [1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5];
julia> k_in   = [2, 3, 4, 1, 3, 5, 2, 4, 1, 2, 4, 0, 4];
julia> F      = [2, 2, 1, 1, 1, 2, 3, 1, 1, 2, 2, 1, 1];
julia> θ      = rand(length(k_out));
julia> nz_out = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13];
julia> nz_in  = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13];
julia> n      = length(k_out);
julia> L_DBCM_reduced(θ, k_out, k_in, F, nz_out, nz_in, n)

# Use with UBCM model:
julia> G = 
julia> model = DBCM(G);
julia> model_fun = θ -> L_DBCM_reduced(θ, model.dᵣ_out, model.dᵣ_in, model.f, model.dᵣ_out_nz, model.dᵣ_in_nz, model.status[:d_unique])
julia> model_fun(model.Θᵣ)
```
"""
function L_DBCM_reduced(θ::Vector, k_out::Vector, k_in::Vector, F::Vector, nz_out::Vector, nz_in::Vector, n::Int=length(k_out))
    α = @view θ[1:n]
    β = @view θ[n+1:end]
    res = zero(eltype(θ))
    for i ∈ nz_out
        @inbounds res -= F[i] * k_out[i] * α[i]
        for j ∈ nz_in
            if i ≠ j 
                @inbounds res -= F[i] * F[j]       * log(1 + exp(-α[i] - β[j]))
            else
                @inbounds res -= F[i] * (F[i] - 1) * log(1 + exp(-α[i] - β[j]))
            end
        end
    end

    for j ∈ nz_in
        @inbounds res -= F[j] * k_in[j]  * β[j]
    end

    return res
end


"""
    L_DBCM_reduced(m::DBCM)

Return the log-likelihood of the DBCM model `m` based on the computed maximum likelihood parameters.

TO DO: include check for parameters computed
"""
L_DBCM_reduced(m::DBCM) = L_DBCM_reduced(m.θᵣ, m.dᵣ_out, m.dᵣ_in, m.f, m.dᵣ_out_nz, m.dᵣ_in_nz, m.status[:d_unique])


"""
    ∇L_DBCM_reduced!(∇L::AbstractVector, θ::AbstractVector, k_out::AbstractVector, k_in::AbstractVector, F::AbstractVector, nz_out::Vector, nz_in::Vector, x::AbstractVector, y::AbstractVector,n::Int)

Compute the gradient of the log-likelihood of the reduced DBCM model using the exponential formulation in order to maintain convexity.

For the optimisation, this function will be used togenerate an anonymous function associated with a specific model. The function 
will update pre-allocated vectors (`∇L`,`x` and `y`) for speed. The gradient is non-allocating.

The arguments of the function are:
    - `∇L`: the gradient of the log-likelihood of the reduced model
    - `θ`: the maximum likelihood parameters of the model ([α; β])
    - `k_out`: the reduced outdegree sequence
    - `k_in`: the reduced indegree sequence
    - `F`: the frequency of each pair in the degree sequence
    - `nz_out`: the indices of non-zero elements in the reduced outdegree sequence
    - `nz_in`: the indices of non-zero elements in the reduced indegree sequence
    - `x`: the exponentiated maximum likelihood parameters of the model ( xᵢ = exp(-αᵢ) )
    - `y`: the exponentiated maximum likelihood parameters of the model ( yᵢ = exp(-βᵢ) )
    - `n`: the number of nodes in the reduced model

# Examples
```jldoctest
# Explicit use with DBCM model:

# Use within optimisation.jl framework:
julia> fun =   (θ, p)
julia> ∇fun! = (∇L, θ, p)
julia> θ₀ =  # initial condition
julia> foo = MaxEntropyGraphs.Optimization.OptimizationProblem(fun, grad=∇fun!)
julia> prob  = MaxEntropyGraphs.Optimization.OptimizationFunction(prob, θ₀)
julia> method = MaxEntropyGraphs.OptimizationOptimJL.NLopt.LD_LBFGS()
julia> solve(prob, method)
...
```
"""
function ∇L_DBCM_reduced!(  ∇L::AbstractVector, θ::AbstractVector, 
                            k_out::AbstractVector, k_in::AbstractVector, 
                            F::AbstractVector, 
                            nz_out::Vector, nz_in::Vector,
                            x::AbstractVector, y::AbstractVector,
                            n::Int)
    # set pre-allocated values
    α = @view θ[1:n]
    β = @view θ[n+1:end]
    @simd for i in eachindex(α) # to obtain a non-allocating function <> x .= exp.(-α), y .= exp.(-β)
        @inbounds x[i] = exp(-α[i])
        @inbounds y[i] = exp(-β[i])
    end
    # reset gradient to zero
    ∇L .= zero(eltype(∇L))
    
    # part related to α
    @simd for i ∈ nz_out
        fx = zero(eltype(∇L))
        for j ∈ nz_in
            if i ≠ j
                @inbounds c = F[i] * F[j]
            else
                @inbounds c = F[i] * (F[j] - 1)
            end
            @inbounds fx += c * y[j] / (1 + x[i] * y[j])
        end
        @inbounds ∇L[i] = x[i] * fx - F[i] * k_out[i]
    end
    # part related to β
    @simd for j ∈ nz_in
        fy = zero(eltype(∇L))
        for i ∈ nz_out
            if i≠j
                @inbounds c = F[i] * F[j]
            else
                @inbounds c = F[i] * (F[j] - 1)
            end
            @inbounds fy += c * x[i] / (1 + x[i] * y[j])
        end
        @inbounds ∇L[n+j] = y[j] * fy - F[j] * k_in[j]
    end

    return ∇L
end


"""
    ∇L_DBCM_reduced_minus!(args...)

Compute minus the gradient of the log-likelihood of the reduced DBCM model using the exponential formulation in order to maintain convexity. Used for optimisation in a non-allocating manner.

See also [`∇L_DBCM_reduced!`](@ref)
"""
function ∇L_DBCM_reduced_minus!(∇L::AbstractVector, θ::AbstractVector,
                                k_out::AbstractVector, k_in::AbstractVector, 
                                F::AbstractVector, 
                                nz_out::Vector, nz_in::Vector,
                                x::AbstractVector, y::AbstractVector,
                                n::Int)
    # set pre-allocated values
    α = @view θ[1:n]
    β = @view θ[n+1:end]
    @simd for i in eachindex(α) # to obtain a non-allocating function <> x .= exp.(-α), y .= exp.(-β)
        @inbounds x[i] = exp(-α[i])
        @inbounds y[i] = exp(-β[i])
    end
    # reset gradient to zero
    ∇L .= zero(eltype(∇L))

    # part related to α
    @simd for i ∈ nz_out
        fx = zero(eltype(∇L))
        for j ∈ nz_in
            if i ≠ j
                @inbounds c = F[i] * F[j]
            else
                @inbounds c = F[i] * (F[j] - 1)
            end
            @inbounds fx -= c * y[j] / (1 + x[i] * y[j])
        end
        @inbounds ∇L[i] = x[i] * fx + F[i] * k_out[i]
    end
    # part related to β
    @simd for j ∈ nz_in
        fy = zero(eltype(∇L))
        for i ∈ nz_out
            if i≠j
                @inbounds c = F[i] * F[j]
            else
                @inbounds c = F[i] * (F[j] - 1)
            end
            @inbounds fy -= c * x[i] / (1 + x[i] * y[j])
        end
        @inbounds ∇L[n+j] = y[j] * fy + F[j] * k_in[j]
    end

    return ∇L
end


"""
    DBCM_reduced_iter!(θ::AbstractVector, k_out::AbstractVector, k_in::AbstractVector, F::AbstractVector, nz_out::Vector, nz_in::Vector,x::AbstractVector, y::AbstractVector, G::AbstractVector, H::AbstractVector, n::Int)

Computer the next fixed-point iteration for the DBCM model using the exponential formulation in order to maintain convexity.
The function is non-allocating and will update pre-allocated vectors (`θ`, `x`, `y`, `G` and `H`) for speed.

The arguments of the function are:
    - `θ`: the maximum likelihood parameters of the model ([α; β])
    - `k_out`: the reduced outdegree sequence
    - `k_in`: the reduced indegree sequence
    - `F`: the frequency of each pair in the degree sequence
    - `nz_out`: the indices of non-zero elements in the reduced outdegree sequence
    - `nz_in`: the indices of non-zero elements in the reduced indegree sequence
    - `x`: the exponentiated maximum likelihood parameters of the model ( xᵢ = exp(-αᵢ) )
    - `y`: the exponentiated maximum likelihood parameters of the model ( yᵢ = exp(-βᵢ) )
    - `G`: buffer for out-degree related computations
    - `H`: buffer for in-degree related computations
    - `n`: the number of nodes in the reduced model


# Examples
```jldoctest
# Use with DBCM model:
julia> G = 
julia> model = DBCM(G);
julia> G = zeros(eltype(model.Θᵣ), length(model.xᵣ);
julia> H = zeros(eltype(model.Θᵣ), length(model.yᵣ);
julia> x = zeros(eltype(model.Θᵣ), length(model.xᵣ);
julia> y = zeros(eltype(model.Θᵣ), length(model.yᵣ);
julia> DBCM_FP! = θ -> DBCM_reduced_iter!(θ, model.dᵣ_out, model.dᵣ_in, model.f, model.dᵣ_out_nz, model.dᵣ_in_nz, x, y, G, H, model.status[:d_unique])
julia> UBCM_FP!(model.Θᵣ)
```
"""
function DBCM_reduced_iter!(θ::AbstractVector, 
                            k_out::AbstractVector, k_in::AbstractVector, 
                            F::AbstractVector, 
                            nz_out::Vector, nz_in::Vector,
                            x::AbstractVector, y::AbstractVector, 
                            G::AbstractVector,  n::Int) # H::AbstractVector,
    α = @view θ[1:n]
    β = @view θ[n+1:end]
    @simd for i in eachindex(α) # to obtain a non-allocating function <> x .= exp.(-α), y .= exp.(-β) (1.8μs, 6 allocs -> 1.2μs, 0 allocs)
        @inbounds x[i] = exp(-α[i])
        @inbounds y[i] = exp(-β[i])
    end
    G .= zero(eltype(G))
    #H .= zero(eltype(H))
    # part related to α
    @simd for i ∈ nz_out
        for j ∈ nz_in
            if i ≠ j
                @inbounds G[i] += F[j]        * y[j] / (1 + x[i] * y[j])
            else
                @inbounds G[i] += (F[j] - 1)  * y[j] / (1 + x[i] * y[j])
            end
        end
        @inbounds G[i] = -log(k_out[i] / G[i])
    end
    # part related to β
    @simd for j ∈ nz_in
        for i ∈ nz_out
            if i ≠ j
                @inbounds G[j+n] += F[i]        * x[i] / (1 + x[i] * y[j])
                #@inbounds H[j] += F[i]        * x[i] / (1 + x[i] * y[j])
            else
                @inbounds G[j+n] += (F[i] - 1)  * x[i] / (1 + x[i] * y[j])
                #@inbounds H[j] += (F[i] - 1)  * x[i] / (1 + x[i] * y[j])
            end
        end
        @inbounds G[n+j] = -log(k_in[j] / G[j+n])
        #@inbounds θ[n+j] = -log(k_in[j] / H[j])
    end

    return G
    #return θ
end



"""
    initial_guess(m::DBCM, method::Symbol=:degrees)

Compute an initial guess for the maximum likelihood parameters of the DBCM model `m` using the method `method`.

The methods available are: `:degrees` (default), `:degrees_minor`, `:random`, `:uniform`, `:chung_lu`.
"""
function initial_guess(m::DBCM{T,N}; method::Symbol=:degrees) where {T,N}
    #N = typeof(m).parameters[2]
    if isequal(method, :degrees)
        return Vector{N}(vcat(-log.(m.dᵣ_out), -log.(m.dᵣ_in)))
    elseif isequal(method, :degrees_minor)
        isnothing(m.G) ? throw(ArgumentError("Cannot compute the number of edges because the model has no underlying graph (m.G == nothing)")) : nothing
        return Vector{N}(vcat(-log.(m.dᵣ_out ./ (sqrt(Graphs.ne(m.G)) + 1)), -log.(m.dᵣ_in ./ (sqrt(Graphs.ne(m.G)) + 1)) ))
    elseif isequal(method, :random)
        return Vector{N}(-log.(rand(N, 2*length(m.dᵣ_out))))
    elseif isequal(method, :uniform)
        return Vector{N}(-log.(0.5 .* ones(N, 2*length(m.dᵣ_out))))
    elseif isequal(method, :chung_lu)
        isnothing(m.G) ? throw(ArgumentError("Cannot compute the number of edges because the model has no underlying graph (m.G == nothing)")) : nothing
        return Vector{N}(vcat(-log.(m.dᵣ_out ./ (2 * Graphs.ne(m.G))), -log.(m.dᵣ_in ./ (2 * Graphs.ne(m.G)))))
    else
        throw(ArgumentError("The initial guess method $(method) is not supported"))
    end
end


"""
    set_xᵣ!(m::DBCM)

Set the value of xᵣ to exp(-αᵣ) for the DBCM model `m`
"""
function set_xᵣ!(m::DBCM)
    if m.status[:params_computed]
        αᵣ = @view m.θᵣ[1:m.status[:d_unique]]
        m.xᵣ .= exp.(-αᵣ)
    else
        throw(ArgumentError("The parameters have not been computed yet"))
    end
end

"""
    set_yᵣ!(m::DBCM)

Set the value of yᵣ to exp(-βᵣ) for the DBCM model `m`
"""
function set_yᵣ!(m::DBCM)
    if m.status[:params_computed]
        βᵣ = @view m.θᵣ[m.status[:d_unique]+1:end]
        m.yᵣ .= exp.(-βᵣ)
    else
        throw(ArgumentError("The parameters have not been computed yet"))
    end
end


"""
    Ĝ(m::DBCM)

Compute the expected adjacency matrix for the DBCM model `m`
"""
function Ĝ(m::DBCM{T,N}) where {T,N}
    # check if possible
    m.status[:params_computed] ? nothing : throw(ArgumentError("The parameters have not been computed yet"))
    
    # get network size => this is the full size
    n = m.status[:d] 
    # initiate G
    G = zeros(N, n, n)
    # initiate x and y
    x = m.xᵣ[m.dᵣ_ind]
    y = m.yᵣ[m.dᵣ_ind]
    # compute G
    for i = 1:n
        @simd for j = 1:n
            if i≠j
                @inbounds xiyj = x[i]*y[j]
                @inbounds G[i,j] = xiyj/(1 + xiyj)
            end
        end
    end

    return G    
end


"""
    set_Ĝ!(m::DBCM)

Set the expected adjacency matrix for the DBCM model `m`
"""
function set_Ĝ!(m::DBCM)
    m.Ĝ = Ĝ(m)
    m.status[:G_computed] = true
    return m.Ĝ
end


"""
    σˣ(m::DBCM{T,N}) where {T,N}

Compute the standard deviation for the elements of the adjacency matrix for the DBCM model `m`.

**Note:** read as "sigma star"
"""
function σˣ(m::DBCM{T,N}) where {T,N}
    # check if possible
    m.status[:params_computed] ? nothing : throw(ArgumentError("The parameters have not been computed yet"))
    # check network size => this is the full size
    n = m.status[:d]
    # initiate G
    σ = zeros(N, n, n)
    # initiate x and y
    x = m.xᵣ[m.dᵣ_ind]
    y = m.yᵣ[m.dᵣ_ind]
    # compute σ
    for i = 1:n
        @simd for j = i+1:n
            @inbounds xiyj =  x[i]*y[j]
            @inbounds xjyi =  x[j]*y[i]
            @inbounds res[i,j] = sqrt(xiyj)/(1 + xiyj)
            @inbounds res[j,i] = sqrt(xjyi)/(1 + xjyi)
        end
    end

    return σ
end

"""
    set_σ!(m::DBCM)

Set the standard deviation for the elements of the adjacency matrix for the DBCM model `m`
"""
function set_σ!(m::DBCM)
    m.σ = σˣ(m)
    m.status[:σ_computed] = true
    return m.σ
end


"""
    rand(m::DBCM; precomputed=false)

Generate a random graph from the DBCM model `m`.

Keyword arguments:
- `precomputed::Bool`: if `true`, the precomputed expected adjacency matrix (`m.Ĝ`) is used to generate the random graph, otherwise the maximum likelihood parameters are used to generate the random graph on the fly. For larger networks, it is 
  recommended to not precompute the expected adjacency matrix to limit memory pressure.

# Examples
```jldoctest
# generate a DBCM model from the karate club network
julia> G = MaxEntropyGraphs.Graphs.SimpleGraphs.smallgraph(:karate);
julia> model = MaxEntropyGraphs.DBCM(G);
# compute the maximum likelihood parameters
using NLsolve
x_buffer = zeros(length(model.dᵣ));G_buffer = zeros(length(model.dᵣ));
FP_model! = (θ::Vector) -> MaxEntropyGraphs.DBCM_reduced_iter!(θ, model.dᵣ, model.f, x_buffer, G_buffer);
sol = fixedpoint(FP_model!, θ₀, method=:anderson, ftol=1e-12, iterations=1000);
model.Θᵣ .= sol.zero;
model.status[:params_computed] = true;
set_xᵣ!(model);
# set the expected adjacency matrix
MaxEntropyGraphs.set_Ĝ!(model);
# sample a random graph
julia> rand(model)
{34, 78} undirected simple Int64 graph
```
"""
function rand(m::DBCM; precomputed::Bool=false)
    if precomputed
        # check if possible to use precomputed Ĝ
        m.status[:G_computed] ? nothing : throw(ArgumentError("The expected adjacency matrix has not been computed yet"))
        # generate random graph
        #G = Graphs.SimpleGraphFromIterator(  Graphs.Edge.([(i,j) for i = 1:m.status[:d] for j in i+1:m.status[:d] if rand()<m.Ĝ[i,j]]))
        G = Graphs.SimpleDiGraphFromIterator( Graphs.Edge.([(i,j) for i = 1:m.status[:d] for j in 1:m.status[:d] if (rand()<m.Ĝ[i,j] && i≠j)  ]))
    else
        # check if possible to use parameters
        m.status[:params_computed] ? nothing : throw(ArgumentError("The parameters have not been computed yet"))
        # initiate x and y
        x = m.xᵣ[m.dᵣ_ind]
        y = m.yᵣ[m.dᵣ_ind]
        # generate random graph
        # G = Graphs.SimpleGraphFromIterator(Graphs.Edge.([(i,j) for i = 1:m.status[:d] for j in i+1:m.status[:d] if rand()< (x[i]*x[j])/(1 + x[i]*x[j]) ]))
        G = Graphs.SimpleDiGraphFromIterator(Graphs.Edge.([(i,j) for i = 1:m.status[:d] for j in   1:m.status[:d] if (rand() < (x[i]*y[j])/(1 + x[i]*y[j]) && i≠j) ]))
    end

    # deal with edge case where no edges are generated for the last node(s) in the graph
    while Graphs.nv(G) < m.status[:d]
        Graphs.add_vertex!(G)
    end

    return G
end


"""
    rand(m::DBCM, n::Int; precomputed=false)

Generate `n` random graphs from the DBCM model `m`. If multithreading is available, the graphs are generated in parallel.

Keyword arguments:
- `precomputed::Bool`: if `true`, the precomputed expected adjacency matrix (`m.Ĝ`) is used to generate the random graph, otherwise the maximum likelihood parameters are used to generate the random graph on the fly. For larger networks, it is 
  recommended to not precompute the expected adjacency matrix to limit memory pressure.

# Examples
```jldoctest
# generate a DBCM model from the karate club network
julia> G = MaxEntropyGraphs.Graphs.SimpleGraphs.smallgraph(:karate);
julia> model = MaxEntropyGraphs.DBCM(G);
# compute the maximum likelihood parameters
using NLsolve
x_buffer = zeros(length(model.dᵣ));G_buffer = zeros(length(model.dᵣ));
FP_model! = (θ::Vector) -> MaxEntropyGraphs.DBCM_reduced_iter!(θ, model.dᵣ, model.f, x_buffer, G_buffer);
sol = fixedpoint(FP_model!, θ₀, method=:anderson, ftol=1e-12, iterations=1000);
model.Θᵣ .= sol.zero;
model.status[:params_computed] = true;
set_xᵣ!(model);
# set the expected adjacency matrix
MaxEntropyGraphs.set_Ĝ!(model);
# sample a random graph
julia> rand(model, 10)
10-element Vector{Graphs.SimpleGraphs.SimpleDiGraph{Int64}}
```
"""
function rand(m::DBCM, n::Int; precomputed::Bool=false)
    # pre-allocate
    res = Vector{Graphs.SimpleDiGraph{Int}}(undef, n)
    # fill vector using threads
    Threads.@threads for i in 1:n
        res[i] = rand(m; precomputed=precomputed)
    end

    return res
end



"""
    solve_model!(m::DBCM)

Compute the likelihood maximising parameters of the DBCM model `m`. 

By default the parameters are computed using the fixed point iteration method with the degree sequence as initial guess.
"""
function solve_model!(m::DBCM{T,N};  # common settings
                                method::Symbol=:fixedpoint, 
                                initial::Symbol=:degrees,
                                maxiters::Int=1000, 
                                verbose::Bool=false,
                                # NLsolve.jl specific settings (fixed point method)
                                ftol::Real=1e-8,
                                # optimisation.jl specific settings (optimisation methods)
                                abstol::Union{Number, Nothing}=nothing,
                                reltol::Union{Number, Nothing}=nothing,
                                AD_method::Symbol=:AutoZygote,
                                analytical_gradient::Bool=false) where {T,N}
    # initial guess
    θ₀ = initial_guess(m, method=initial)
    # find Inf values
    ind_inf = findall(isinf, θ₀)
    if method==:fixedpoint
        # initiate buffers
        x_buffer = zeros(N, length(m.dᵣ_out)); # buffer for x = exp(-α)
        y_buffer = zeros(N, length(m.dᵣ_in));  # buffer for y = exp(-β)
        G_buffer = zeros(N, length(m.θᵣ)); # buffer for G(x)
        # define fixed point function
        FP_model! = (θ::Vector) -> DBCM_reduced_iter!(θ, m.dᵣ_out, m.dᵣ_in, m.f, m.dᵣ_out_nz, m.dᵣ_in_nz, x_buffer, y_buffer, G_buffer, m.status[:d_unique])
        # obtain solution
        θ₀[ind_inf] .= zero(N);
        sol = NLsolve.fixedpoint(FP_model!, θ₀, method=:anderson, ftol=ftol, iterations=maxiters);
        if NLsolve.converged(sol)
            if verbose 
                @info "Fixed point iteration converged after $(sol.iterations) iterations"
            end
            m.θᵣ .= sol.zero;
            m.θᵣ[ind_inf] .= Inf;
            m.status[:params_computed] = true;
            set_xᵣ!(m);
            set_yᵣ!(m);
        else
            throw(ConvergenceError(method, nothing))
        end
    else
        if analytical_gradient
            # initiate buffers
            x_buffer = zeros(N, length(m.dᵣ_out)); # buffer for x = exp(-α)
            y_buffer = zeros(N, length(m.dᵣ_in));  # buffer for y = exp(-β)
            # initialise gradient buffer
            #gx_buffer = similar(θ₀)
            # define gradient function for optimisation.jl
            grad! = (G, θ, p) -> ∇L_DBCM_reduced_minus!(G, θ, m.dᵣ_out, m.dᵣ_in, m.f, m.dᵣ_out_nz, m.dᵣ_in_nz, x_buffer, y_buffer, m.status[:d_unique]);
        end
        # define objective function and its AD method
        f = AD_method ∈ keys(AD_methods)            ? Optimization.OptimizationFunction( (θ, p) ->   -L_DBCM_reduced(θ, m.dᵣ_out, m.dᵣ_in, m.f, m.dᵣ_out_nz, m.dᵣ_in_nz, m.status[:d_unique]),
                                                                                         AD_methods[AD_method],
                                                                                         grad = analytical_gradient ? grad! : nothing)                      : throw(ArgumentError("The AD method $(AD_method) is not supported (yet)"))
        prob = Optimization.OptimizationProblem(f, θ₀);
        # obtain solution
        sol = method ∈ keys(optimization_methods)   ? Optimization.solve(prob, optimization_methods[method], abstol=abstol, reltol=reltol)                                                : throw(ArgumentError("The method $(method) is not supported (yet)"))
        # check convergence
        if Optimization.SciMLBase.successful_retcode(sol.retcode)
            if verbose 
                @info """$(method) optimisation converged after $(@sprintf("%1.2e", sol.solve_time)) seconds (Optimization.jl return code: $("$(sol.retcode)"))"""
            end
            m.Θᵣ .= sol.u;
            m.status[:params_computed] = true;
            set_xᵣ!(m);
            set_yᵣ!(m);
        else
            throw(ConvergenceError(method, sol.retcode))
        end
    end

    return m
end
