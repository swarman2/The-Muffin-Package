include("FC.jl")
include("HALF.jl")
include("INT.jl")
include("EBM.jl")
include("HBM.jl")
include("helper_functions.jl")
include("MID.jl")
include("GAP.jl")
include("Train.jl")

function UPPER_BOUND(m,s)
       start_time = time()
       V, sᵥ,sᵥ₋₁ = SV(m,s)
       t_multi = 0
       t_solv =0
       FC_alpha = 1
       INT_alpha = 1
       HALF_alpha = 1
       EBM_alpha = 1
       HBD_alpha = 1
       MID_alpha = 1
       GAP_alpha = 1
       TRAIN_alpha = 1
       if s>=3 && s<=9
          if (m==11 && s==5) || (m==7 && s==6)||(m==8 && s==7)||(m==19 && s==7)||(m==10 && s==9)||(m==11 && s==9)||(m==29&&s==9)||(m==38 && s==9)||(m==47 && s==9)
             min_alpha = INT(m,s)
             INT_alpha = min_alpha
          else
             min_alpha = FC(m,s)
             FC_alpha = min_alpha
          end
       else
          FC_alpha = FC(m,s)
          INT_alpha = INT(m,s)
          HALF_alpha = HALF(m,s)
          EBM_alpha = EBM(m,s)
          HBD_alpha = HBM(m,s)
          MID_alpha, MID_endpoints = MID(m,s, 1//2, true)
          alphas = [FC_alpha INT_alpha HALF_alpha EBM_alpha HBD_alpha MID_alpha]
          min_alpha = minimum(alphas)
          GAP_alpha=1 #forward defintion
          TRAIN_alpha = 1
          if(min_alpha != 1//3)
             GAP_alpha, GAP_endpoints = GAP(m,s,min_alpha, true)
             if GAP_alpha < min_alpha
                min_alpha = GAP_alpha
             end
          end
          if min_alpha!=1//3
             TRAIN_alpha = TRAIN(m,s,min_alpha)
             if TRAIN_alpha< min_alpha
                min_alpha = TRAIN_alpha
             end
          end
       end

       total_time = time()-start_time
       (minutes, seconds) = fldmod(total_time, 60)
       (hours, minutes) = fldmod(minutes, 60)
       time_str = @sprintf("%02d:%02d:%0.2f", hours, minutes, seconds)
       return min_alpha, time_str
end
for s=3:60
   for m=s+1:70
      if gcd(m,s)==1
         min_alpha, time_str = UPPER_BOUND(m,s)
         println(m,"  | ",s,"  |  ≤ ",min_alpha," | ",time_str)
      end
   end
end
