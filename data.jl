module Basis

struct BasisSet
    nstates  :: Int64
    veemax   :: Vector{Int64}
    lambda   :: Vector{Float64}
    ess      :: Vector{Float64}
    jaylist  :: Vector{Vector{Vector{Float64}}}
    eledimen:: Vector{Int64}
    vibdimen :: Vector{Int64}
    rotdimen :: Vector{Int64}
end

function BasisSet(
    veemax_s :: Vector{Float64},
    lambda_s :: Vector{Float64},
    ess_s    :: Vector{Float64},
    jaylist  :: Vector{Float64}
    )
    nstates = length(lambda_s)
    jaylist_s = Vector{Vector{Vector{Float64}}}(undef, nstates)
    for state=1:nstates
        ess = ess_s[state]
        lambda = lambda_s[state]
        jaylist_sigma = Vector{Vector{Float64}}(undef, floor(Int64, 2*ess + 1))
        for (sig, sigma) in enumerate(-ess:ess)
            absomega = abs(lambda + sigma)
            jaylist_sigma[sig] = [jay for jay in jaylist if jay ∉ absomega-floor(absomega):absomega-1.0]
        end
        jaylist_s[state] = jaylist_sigma
    end
    eledimen_s = Vector{Float64}(undef, nstates)
    vibdimen_s = Vector{Float64}(undef, nstates)
    rotdimen_s = Vector{Float64}(undef, nstates)
    @. eledimen_s = floor(Int64, 2*ess_s + 1) 
    @. vibdimen_s  = floor(Int64, veemax_s)
    @. rotdimen_s = (2 * ess_s + 1)*length(jaylist_s)
    BasisSet(
        nstates, veemax_s, lambda_s, ess_s, jaylist_s,
        eledimen_s, vibdimen_s, rotdimen_s
    )
end

# function quantum_numbers(basis::BasisSet, i)
#     statei = (i ÷ basis.nstates) + 1
#     i %= basis.nstates + 1
#     elei = (i ÷ basis.eledimen[statei]) + 1
#     i %= basis.eledimen[statei] + 1
#     vibi = i ÷ basis.vibdimen[statei] + 1
#     i %= basis.vibdimen[statei] + 1
#     forbidjays = basis.forbidjays[statei][elei]
#     for jay in basis.jaylist
#         if jay in forbidjays
#             continue
#         else
#             ndimen += floor(Int64, 2*jay + 1)
        
            
# end

end