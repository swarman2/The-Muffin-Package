include("PROC.jl")
include("FC.jl")
include("HALF.jl")
include("INT.jl")
include("EBM.jl")
include("HBM.jl")
include("helper_functions.jl")
include("Mid.jl")
include("GAP.jl")
println()

function FIND_ALPHA(m,s)

         FC_alpha = FC(m,s)
         INT_alpha = INT(m,s)
         HALF_alpha = HALF(m,s)
         EBM_alpha = EBM(m,s)
         HBD_alpha = HBM(m,s)
         MID_alpha = MID(m,s)
         println(MID(m,s))
         #GAP_alpha = GAP(m,s)

alphas = [FC_alpha INT_alpha HALF_alpha EBM_alpha HBD_alpha MID_alpha] #GAP_alpha]
min_alpha = minimum(alphas)
 #println(min_alpha)
if min_alpha ==1
   print("   ",m,"   |   ",s,"   |  trivial   |  = ",1)

#println("   M   |   S   |   Method   |      α   ")
#println("---------------------------------------")
else
      print("   ",m,"   |   ",s,"    |")

   if      FC_alpha <= min_alpha

          # if VProc(m,s,FC_alpha)
              print("     FC     ")
             # print("   ",m,"   |   ",s,"   |     FC     |  = ",FC_alpha,"")
              #println("f(",m,",",s,") = ",FC_alpha, "   Using: Easy Buddy Match Method")
           #else
            #  print("   ",m,"   |   ",s,"   |     FC     |  ≤ ",FC_alpha,"")
              #println("f(",m,",",s,") ≤ ",FC_alpha, "   Using: Easy Buddy Match Method")
           end
   if  INT_alpha <= min_alpha
    print("    INT     ")
      #     if VProc(m,s,INT_alpha)
             # print("   ",m,"  |   ",s,"   |    INT     |  = ",INT_alpha,"")
      #     else
      #        print("   ",m,"  |   ",s,"   |    INT     |  ≤ ",INT_alpha,"")
      #     end

           #VINT_proof(m,s,min_alpha)
           #println()
           #PROC(m,s,min_alpha)
   end
   if  HALF_alpha <= min_alpha
      print("   HALF    ")
           #if VProc(m,s,HALF_alpha)
            #  print("   ",m,"   |   ",s,"   |   HALF    |  = ",HALF_alpha,"")
           #else
            #  print("   ",m,"   |   ",s,"   |   HALF    |  ≤ ",HALF_alpha,"")
           end


   if EBM_alpha <= min_alpha
      print("    EBM     ")
      #    if VProc(m,s,EBM_alpha)
      #       print("   ",m,"  |   ",s,"  |    EBM     |  = ",EBM_alpha,"")
      #    else
      #       print("   ",m,"  |   ",s,"  |    EBM     |  ≤ ",EBM_alpha,"")
          end


   if HBD_alpha <= min_alpha
   print("    HBM     ")
      #    if VProc(m,s,HBD_alpha)
      #       print("   ",m,"  |   ",s,"  |    HBM     |  = ",HBD_alpha,"")
      #    else
      #       print("   ",m,"  |   ",s,"  |    HBM     |  ≤ ",HBD_alpha,"")
      #    end

   end

   if MID_alpha <= min_alpha
      print("    MID     ")
   end

   #if GAP_alpha <= min_alpha
   #   print("    GAP     ")
   #end


   if VProc(m,s,min_alpha)
      print("|  = ",min_alpha,"")
   else
      print("|  ≤ ",min_alpha,"")
   end
end
end



println("   M   |   S   |   Method   |      α   ")
println("---------------------------------------")
for s = 3:60
   for m = s+1:70
      if(gcd(m,s)==1)
         FIND_ALPHA(m,s)
         println()
      end
   end
end
