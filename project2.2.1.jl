include("PROC.jl")
include("FC.jl")
include("HALF.jl")
include("INT.jl")
include("EBM.jl")
include("HBM.jl")
include("helper_functions.jl")
println()

function FIND_ALPHA(m,s)

         FC_alpha = FC(m,s)
         INT_alpha = INT(m,s)
         HALF_alpha = HALF(m,s)
         EBM_alpha = EBM(m,s)
         HBD_alpha = HBM(m,s)

alphas = [FC_alpha INT_alpha HALF_alpha EBM_alpha HBD_alpha]
min_alpha = minimum(alphas)
 #println(min_alpha)
if min_alpha ==1
   print("   ",m,"   |   ",s,"   |  trivial   |  = ",1)

#println("   M   |   S   |   Method   |      α   ")
#println("---------------------------------------")

elseif      FC_alpha <= min_alpha

        if VProc(m,s,FC_alpha)
           print("   ",m,"   |   ",s,"   |     FC     |  = ",FC_alpha,"")
           #println("f(",m,",",s,") = ",FC_alpha, "   Using: Easy Buddy Match Method")
        else
           print("   ",m,"   |   ",s,"   |     FC     |  ≤ ",FC_alpha,"")
           #println("f(",m,",",s,") ≤ ",FC_alpha, "   Using: Easy Buddy Match Method")
        end

elseif  INT_alpha <= min_alpha

        if VProc(m,s,INT_alpha)
           print("   ",m,"  |   ",s,"   |    INT     |  = ",INT_alpha,"")
        else
           print("   ",m,"  |   ",s,"   |    INT     |  ≤ ",INT_alpha,"")
        end

        #VINT_proof(m,s,min_alpha)
        #println()
        #PROC(m,s,min_alpha)

elseif  HALF_alpha <= min_alpha

        if VProc(m,s,HALF_alpha)
           print("   ",m,"   |   ",s,"   |   HALF    |  = ",HALF_alpha,"")
        else
           print("   ",m,"   |   ",s,"   |   HALF    |  ≤ ",HALF_alpha,"")
        end


elseif EBM_alpha <= min_alpha

       if VProc(m,s,EBM_alpha)
          print("   ",m,"  |   ",s,"  |    EBM     |  = ",EBM_alpha,"")
       else
          print("   ",m,"  |   ",s,"  |    EBM     |  ≤ ",EBM_alpha,"")
       end


elseif HBD_alpha <= min_alpha

       if VProc(m,s,HBD_alpha)
          print("   ",m,"  |   ",s,"  |    HBM     |  = ",HBD_alpha,"")
       else
          print("   ",m,"  |   ",s,"  |    HBM     |  ≤ ",HBD_alpha,"")
       end

end
end



println("   M   |   S   |   Method   |      α   ")
println("---------------------------------------")
for s = 5:60
   for m = s+1:70
      FIND_ALPHA(m,s)
      println()
   end
end
