include("PROC.jl")
include("FC.jl")
include("HALF.jl")
include("INT.jl")
println()
println()

function FIND_ALPHA(m,s)
         FC_alpha = FC(m,s)
         INT_alpha = INT(m,s)
         HALF_alpha = HALF(m,s)


alphas = [FC_alpha INT_alpha HALF_alpha]
min_alpha = minimum(alphas)


if      FC_alpha <= min_alpha
        println("Using: FC Method")
        println()
        PROC(m,s,min_alpha)

elseif  INT_alpha <= min_alpha
        println("Using: INT Method, finding an upper bound")
        VINT_proof(m,s,min_alpha)
        println()
        PROC(m,s,min_alpha)

elseif  HALF_alpha <= min_alpha
        println("Using: HALF Method, finding a lower bound")
        Half_proof(m,s,min_alpha)
        println()
        PROC(m,s,min_alpha)

end
println()
println("α is equal to ",min_alpha,"")
end
