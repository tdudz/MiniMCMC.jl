module MiniMCMC

using Zygote

include("samplers/HMC.jl")

export HMC, HMC_loop

end # module
