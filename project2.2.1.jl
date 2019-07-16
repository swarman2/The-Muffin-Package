include("PROC.jl")
include("FC.jl")
include("HALF.jl")
include("INT.jl")
include("EBM.jl")
include("HBM.jl")
include("helper_functions.jl")
include("Mid.jl")
include("GAP.jl")
include("one_third_formula.jl")
using Printf
println()

function FIND_ALPHA(m,s,A=0, time_limit_multi = Inf, time_limit_solv = Inf, file = 0, perc_array =0)
   start_time = time()
   V, sᵥ,sᵥ₋₁ = SV(m,s)
   X_conj = false
   d = Int64(m - s)
   k =Int64(floor((s/(3d))))
   if s%3d ==0
      k = k-1
   end
   #println(V," ",d," ",k)
   a = Int64(s-3d*k)
   if s>=3 && s<=9
      FC_alpha = 1
      INT_alpha = 1
      HALF_alpha = 1
      EBM_alpha = 1
      HBD_alpha = 1
      MID_alpha = 1
      GAP_alpha = 1
      proc_bool = true
      if (m==11 && s==5) || (m==7 && s==6)||(m==8 && s==7)||(m==19 && s==7)||(m==10 && s==9)||(m==11 && s==9)||(m==29&&s==9)||(m==38 && s==9)||(m==47 && s==9)
         min_alpha = INT(m,s)
         INT_alpha = min_alpha
      else
         min_alpha = FC(m,s)
         FC_alpha = min_alpha
      end
   else
      if k>=1 && A!=0 && haskey(A,a+d) && V==3
         X_conj = true
         X =A[a+d]*a
         min_alpha_X = (d*k+X)//(3*d*k+a)
         #println(min_alpha)
      elseif k>=1 && a<=2 && A!=0 && V==3
         A[a+d]=1//2
         X_conj = true
         X =A[a+d]*a
         min_alpha_X = (d*k+X)//(3*d*k+a)
      end
      FC_alpha = FC(m,s)
      INT_alpha = INT(m,s)
      HALF_alpha = HALF(m,s)
      EBM_alpha = EBM(m,s)
      HBD_alpha = HBM(m,s)
      MID_alpha = MID(m,s)
      alphas = [FC_alpha INT_alpha HALF_alpha EBM_alpha HBD_alpha MID_alpha]
      min_alpha = minimum(alphas)
      if !X_conj
         GAP_alpha=1 #forward defintion
         if(min_alpha != 1//3)
            GAP_alpha = GAP(m,s,min_alpha)
            if GAP_alpha < min_alpha
               min_alpha = GAP_alpha
            end
         end
      else
         #println(min_alpha_X)
         GAP_alpha =1
         if min_alpha != min_alpha_X
          if  VGAPV3(m,s,min_alpha_X) && min_alpha_X < min_alpha
             min_alpha = min_alpha_X
             GAP_alpha = min_alpha
         end
         if min_alpha != min_alpha_X
            GAP_alpha = GAP(m,s)
            if min_alpha > GAP_alpha
               min_alpha = GAP_alpha
            end
         end
       end
      end
      if min_alpha == 1//3
         proc_bool = one_thrd(m,s,0)
      else
         #println(min_alpha)
         proc_start_time=time()
         proc_bool, err_mess = VProc(m,s,min_alpha, time_limit_solv, time_limit_multi)
         proc_end_time=time()
      end
   end
   if proc_bool == true && k==0 && A!=0
      A[a+d] = min_alpha
   end
    #println(min_alpha)
   if min_alpha ==1
      print("   ",m,"   |   ",s,"   |  trivial   |  = ",1)

   #println("   M   |   S   |   Method   |      α   ")
   #println("---------------------------------------")
   else

         print("   ",m,"   |   ",s,"    |")# PROC_TIME =  ")
         if file !=0
            str = @sprintf("   %i   |   %i    |",m,s)
            write(file,str)
         end
         #@printf("%0.2f",proc_end_time-proc_start_time)

      if      FC_alpha <= min_alpha

             # if VProc(m,s,FC_alpha)
                 print("     FC     ")
                 if file!=0
                     write(file,"     FC     ")
                  end
                # print("   ",m,"   |   ",s,"   |     FC     |  = ",FC_alpha,"")
                 #println("f(",m,",",s,") = ",FC_alpha, "   Using: Easy Buddy Match Method")
              #else
               #  print("   ",m,"   |   ",s,"   |     FC     |  ≤ ",FC_alpha,"")
                 #println("f(",m,",",s,") ≤ ",FC_alpha, "   Using: Easy Buddy Match Method")
              end
      if  INT_alpha <= min_alpha
       print("    INT     ")
       if file !=0
          write(file,"    INT     ")
       end
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
         if file !=0
            write(file,"   HALF    ")
         end
              #if VProc(m,s,HALF_alpha)
               #  print("   ",m,"   |   ",s,"   |   HALF    |  = ",HALF_alpha,"")
              #else
               #  print("   ",m,"   |   ",s,"   |   HALF    |  ≤ ",HALF_alpha,"")
              end


      if EBM_alpha <= min_alpha
         print("    EBM     ")
            if file !=0
               write(file,"    EBM     ")
            end
         #    if VProc(m,s,EBM_alpha)
         #       print("   ",m,"  |   ",s,"  |    EBM     |  = ",EBM_alpha,"")
         #    else
         #       print("   ",m,"  |   ",s,"  |    EBM     |  ≤ ",EBM_alpha,"")
             end


      if HBD_alpha <= min_alpha
      print("    HBM     ")
         if file !=0
            str = @sprintf("    HBM     ")
            write(file,str)
         end
         #    if VProc(m,s,HBD_alpha)
         #       print("   ",m,"  |   ",s,"  |    HBM     |  = ",HBD_alpha,"")
         #    else
         #       print("   ",m,"  |   ",s,"  |    HBM     |  ≤ ",HBD_alpha,"")
         #    end

      end

      if MID_alpha <= min_alpha
         print("    MID     ")
         if file !=0
            str = @sprintf("    MID     " )
            write(file,str)
         end
      end

      if GAP_alpha <= min_alpha
         print("    GAP     ")
         if file !=0
            str = @sprintf("    GAP     ")
            write(file,str)
         end
      end

      if proc_bool == true
         if perc_array !=0
            if FC_alpha == min_alpha
               perc_array[1] = perc_array[1]+1
            elseif HALF_alpha == min_alpha
               perc_array[2]=perc_array[2]+1
            elseif INT_alpha == min_alpha
               perc_array[3]=perc_array[3]+1
            elseif MID_alpha == min_alpha
               perc_array[4]=perc_array[4]+1
            elseif EBM_alpha == min_alpha
               perc_array[5]=perc_array[5]+1
            elseif HBD_alpha == min_alpha
               perc_array[6]=perc_array[6]+1
            elseif GAP_alpha == min_alpha
               perc_array[7]=perc_array[7]+1
            end
         end

         print("|  = ",min_alpha,"")
         if file !=0
            str = @sprintf("|  = %i/%i",numerator(min_alpha),denominator(min_alpha))
            write(file,str)
         end
      elseif proc_bool == -1
         print("|  ≤ ",min_alpha,"|   TIME OUT")
         if file !=0
            str = @sprintf("|  ≤ %i/%i  |   TIME OUT",numerator(min_alpha),denominator(min_alpha))
            write(file,str,"\n")
            write(file,err_mess)
         end
         if perc_array !=0
            perc_array[8] = perc_array[8]+1
         end

      else
         print("|  ≤ ",min_alpha)
         if file !=0
            str = @sprintf("|  ≤ %i/%i",numerator(min_alpha),denominator(min_alpha))
            write(file,str)
         end
         if perc_array !=0
            perc_array[9] = perc_array[9]+1
         end
      end
      @printf("     | runtime:  %0.2f",time()-start_time)
      if file !=0
         str = @sprintf("     | runtime:  %0.2f",time()-start_time)
         write(file,str)
      end
   end
end
#FIND_ALPHA(107,13)

function FIND_ALL(max_m, max_s)
   time_limit_multi = 120
   time_limit_solv =  240
   perc_array = zeros(9)
   filename = "D:/Documents/Muffins/"*string(max_m)*"_"*string(max_s)*"_"*string(time_limit_multi)*"s_"*string(time_limit_solv)*"s.txt"
   file = open(filename, "w")
   A=Dict{Int64,Rational}()
   timed_out=Array{Int64}(undef,0)
   start=time()
      println("   M   |   S   |   Method   |      α   ")
      println("---------------------------------------")
      str = @sprintf("   M   |   S   |   Method   |      α   \n---------------------------------------\n")
      write(file, str)
      for s = 3:max_s
         for m = s+1:max_m
            if gcd(m,s)==1
                     t=time()
                     FIND_ALPHA(m,s,A,time_limit_multi, time_limit_solv, file, perc_array)
                     elapsed = time()-t
                     println()
                     println(file,"\n")
                     #if elapsed >1
                     #   println(elapsed)
                     #end
                     #if elapsed >10
                     #   append!(timed_out,m)
                     #   append!(timed_out,s)
                     #end
            end
         end
      end
   println(time()-start)
   println(file,time()-start)
   @printf("\nFC: %0.2f %% ",100*(perc_array[1])/sum(perc_array))
   @printf("\nHALF: %0.2f %% ",100*(perc_array[2])/sum(perc_array))
   @printf("\nINT: %0.2f %% ",100*(perc_array[3])/sum(perc_array))
   @printf("\nMID: %0.2f %% ",100*(perc_array[4])/sum(perc_array))
   @printf("\nEBM: %0.2f %% ",100*(perc_array[5])/sum(perc_array))
   @printf("\nHBM: %0.2f %% ",100*(perc_array[6])/sum(perc_array))
   @printf("\nGAP: %0.2f %% ",100*(perc_array[7])/sum(perc_array))
   @printf("\nTIME OUT: %0.2f %% ",100*(perc_array[8])/sum(perc_array))
   @printf("\n≠: %0.2f %% ",100*(perc_array[9])/sum(perc_array))

   @printf(file, "\nFC: %0.2f %% ",100*(perc_array[1])/sum(perc_array))
   @printf(file, "\nHALF: %0.2f %% ",100*(perc_array[2])/sum(perc_array))
   @printf(file, "\nINT: %0.2f %% ",100*(perc_array[3])/sum(perc_array))
   @printf(file, "\nMID: %0.2f %%",100*(perc_array[4])/sum(perc_array))
   @printf(file, "\nEBM: %0.2f %% ",100*(perc_array[5])/sum(perc_array))
   @printf(file, "\nHBM: %0.2f %% ",100*(perc_array[6])/sum(perc_array))
   @printf(file, "\nGAP: %0.2f %% ",100*(perc_array[7])/sum(perc_array))
   @printf(file, "\nTIME OUT: %0.2f %% ",100*(perc_array[8])/sum(perc_array))
   @printf(file, "\n≠: %0.2f %% ",100*(perc_array[9])/sum(perc_array))
   close(file)
end
FIND_ALL(210,200)
#FIND_ALPHA(33,26)
