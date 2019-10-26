# HMC sampler for MCMC
using Distributions

function HMC(U::Function, grad_U::Function, ϵ::Float64, L::Int, current_q)
    q = current_q
    p = rand(Distributions.MvNormal(size(q)[1], 1))
    current_p = p
    
    # make a half step for momentum at the beginning
    p = p - ϵ * grad_U(q) / 2

    # alternate full steps for position and momentum
    for i = 1:L
        # make a full step for position q
        q = q + ϵ * p

        # make a full step for the momentum except at end of trajectory
        if i != L
            p = p - ϵ * grad_U(q)
        end
    end

    # make a half step for momentum at the end
    p = p - ϵ * grad_U(q) / 2

    # negate momentum at the end to make proposal symmetric
    p = -p

    # evaluate potential and kinetic energies at start and end of trajectory
    current_U = U(current_q)
    current_K = sum(current_p.^2) / 2
    proposed_U = U(q)
    proposed_K = sum(p.^2) / 2

    # accept or reject the state at end of trajectory, returning either 
    # the position at the end of the trajectory or the initial position
    
    if (rand() < exp(current_U - proposed_U + current_K - proposed_K))
        return q
    else
        return current_q
    end
end

function HMC_sample(U::Function, grad_U::Function, ϵ::Float64, L::Int, current_q::T, n::Int) where T
    samples  = T[]
    for i = 1:n
        current_q = HMCp(U, grad_U, ϵ, L, current_q)
        push!(samples, current_q)
    end 
    return samples
end

