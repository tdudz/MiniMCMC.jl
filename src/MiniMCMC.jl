module MiniMCMC

using Zygote

include("samplers/HMC.jl")

export HMC, HMC_sample

end # module
