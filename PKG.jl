include("PROC.jl")
include("FC.jl")
include("HALF.jl")
include("INT.jl")
include("EBM.jl")
include("HBM.jl")
include("helper_functions.jl")
include("MID.jl")
include("GAP.jl")
include("Train.jl")
include("one_third_formula.jl")
using Printf
using DataFrames
println()
# load Taro - Pkg to read Excel Data
using CSV

using DelimitedFiles
function FIND_ALPHA(m,s, time_limit_multi = Inf, time_limit_solv = Inf, perc_array =0, scott_alpha = 0)

   start_time = time()
   V, sᵥ,sᵥ₋₁ = SV(m,s)
   t_multi = 0
   t_solv =0
   if s>=3 && s<=9
      FC_alpha = 1
      INT_alpha = 1
      HALF_alpha = 1
      EBM_alpha = 1
      HBD_alpha = 1
      MID_alpha = 1
      GAP_alpha = 1
      TRAIN_alpha = 1
      proc_bool = true
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
      if scott_alpha != 0
         proc_bool = (scott_alpha == min_alpha)
         t_multi = 0
         t_solv = 0
      elseif min_alpha == 1//3
         proc_bool = one_thrd(m,s,0)
      elseif min_alpha == MID_alpha && MID_endpoints !=1
         unique!(MID_endpoints)
         proc_bool, err_mess, t_multi, t_solv = VProc(m,s,min_alpha, time_limit_solv, time_limit_multi, MID_endpoints)
      elseif min_alpha == GAP_alpha && GAP_endpoints !=1
         unique!(GAP_endpoints)
         proc_bool, err_mess, t_multi, t_solv = VProc(m,s,min_alpha, time_limit_solv, time_limit_multi, GAP_endpoints)
      else
         proc_bool, err_mess, t_multi, t_solv = VProc(m,s,min_alpha, time_limit_solv, time_limit_multi, 0)
      end

   end
   if min_alpha ==1 && lcm(m,s)==1 || s==1
      print("   ",m,"   |   ",s,"   |  trivial   |  = ",1)
      str = " trivial "
      str_eq = " = "
      min_alpha = 1
   else
      print("   ",m,"   |   ",s,"    |")# PROC_TIME =  ")
      str=""
      if FC_alpha <= min_alpha
         print(" FC ")
         str = str * @sprintf(" FC ")
      end

      if INT_alpha <= min_alpha
        print(" INT ")
        str = str * @sprintf(" INT ")
      end

      if  HALF_alpha <= min_alpha
         print(" HALF ")
         str =str* @sprintf(" HALF ")
      end

      if EBM_alpha <= min_alpha
         print(" EBM ")
         str =str* @sprintf(" EBM ")
      end


      if HBD_alpha <= min_alpha
         print(" HBM ")
         str =str* @sprintf(" HBM ")
      end

      if MID_alpha <= min_alpha
         print(" MID ")
         str =str* @sprintf(" MID " )
      end

      if GAP_alpha <= min_alpha
         print(" GAP ")
         str = str*@sprintf(" GAP ")
      end

      if TRAIN_alpha <= min_alpha
         print(" TRAIN ")
         str = str *@sprintf(" TRAIN ")
      end
      str_eq = ""
      if proc_bool == true
         if perc_array !=0
            if FC_alpha == min_alpha
               perc_array[1] = perc_array[1]+1
            end
            if HALF_alpha == min_alpha
               perc_array[2]=perc_array[2]+1
            end
            if INT_alpha == min_alpha
               perc_array[3]=perc_array[3]+1
            end
            if MID_alpha == min_alpha
               perc_array[4]=perc_array[4]+1
            end
            if EBM_alpha == min_alpha
               perc_array[5]=perc_array[5]+1
            end
            if HBD_alpha == min_alpha
               perc_array[6]=perc_array[6]+1
            end
            if GAP_alpha == min_alpha
               perc_array[7]=perc_array[7]+1
            end
            if TRAIN_alpha == min_alpha
               perc_array[11]=perc_array[11]+1
            end
            perc_array[10]=perc_array[10]+1
         end
         str_eq = " = "
         print("|  = ",min_alpha,"")

      elseif proc_bool == -1
         str_eq = " TIME OUT "
         print("|  ≤ ",min_alpha,"|   TIME OUT")
         if perc_array !=0
            perc_array[8] = perc_array[8]+1
         end
      else
         str_eq = " ≤ "
         print("|  ≤ ",min_alpha)
         if perc_array !=0
            perc_array[9] = perc_array[9]+1
         end
      end
      @printf("     | runtime:  %0.2f",time()-start_time)
   end
   return min_alpha, str, str_eq, round(time()-start_time,digits = 2), t_multi, t_solv
end
#FIND_ALPHA(107,13)

function FIND_ALL(max_m,min_s, max_s, fileKey=0, have_file = false, Scott_data = false)
   time_limit = 600
   println("FILE KEY = ",fileKey)
   time_limit_multi = time_limit
   time_limit_solv =  time_limit
   perc_array = zeros(11)
   if Scott_data
      df = CSV.read("C:/Users/Steph/Downloads/scott_bigrun.csv")
   end

   #fileSolv = open(dirname(@__FILE__)*"\\July30_ScottVerified_Muffins_solv_time_"* string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".txt","w")
   #fileMulti = open(dirname(@__FILE__)*"\\July30_ScottVerified_Muffins_multi_time_"* string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".txt","w")
   #file_backup =open(dirname(@__FILE__)*"\\backup.txt","w")
   if Scott_data
      file = open(dirname(@__FILE__)*"\\July30_ScottVerified_"* string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".txt","w")
   else
      file = open(dirname(@__FILE__)*"\\July30_PROCVerified_"* string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".txt","w")
   end
   #file2 = open(dirname(@__FILE__)*"\\July30_ScottVerified_not_match_SCOTT_"* string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".txt","w")
   row_df = 1
   #println("test")
   #timed_out=Array{Int64}(undef,0)
   start=time()
      println("   M   |   S   |   Method   |      α   |       runtime   ")
      println("---------------------------------------")
      str = @sprintf("   M   |   S   |   Method   |      α   \n---------------------------------------\n")
      #DF = DataFrame(A = Int64[], B = String[])
      #data_frame = DataFrame( m = Int64[], s = Int64[] , Method = String[] ,Equality = String[] , alpha = Any[] , Runtime = Any[] )
      for s = min_s:max_s
         for m = s+1:max_m
            if gcd(m,s)==1
                     #t=time()
                  #   println(file,"test1")
                  #   flush(file)
                     if Scott_data
                        while !(df[:m][row_df]==m && df[:s][row_df]==s)
                           row_df=row_df +1
                        end
                        str = df[Symbol("f(m,s)")][row_df]
                        str_arr=split(str, "/")
                        numerat = parse(Int64, str_arr[1])
                        denominat = parse(Int64, str_arr[2])
                        scott_alpha = numerat//denominat
                        #if scott_alpha != alpha
                        #   @printf(file2, "\nm: %i  s: %i  | Scott alpha: %i/%i  | other alpha: %i/%i", m,s,numerat,denominat, numerator(alpha), denominator(alpha))
                        #   flush(file2)
                        #   p = true
                        #end
                     else
                        scott_alpha =0
                     end
                     alpha,string,str_eq,time, t_multi, t_solv = FIND_ALPHA(m,s,time_limit_multi, time_limit_solv, perc_array, scott_alpha)
                  #   println(file,"test2")
                  #   flush(file)
                     #push!(data_frame,(m,s,string,str_eq,alpha,time))
                     #t_multi = round(t_multi, digits=2)
                     #t_solv = round(t_solv, digits =2)
                     #p = false
                     #if t_multi > 480
                     #   println(fileMulti, m," | ",s," | multi: ",t_multi," | solv: ", t_solv)
                     #   p = true
                     #end
                     #if t_solv > 480
                     #   println(fileSolv, m," | ",s," | multi: ",t_multi," | solv: ", t_solv)
                     #   p = true
                     #end
                     if !Scott_data
                        println(file, m," | ",s," | ",string," |",str_eq," ",alpha," | multi: ",t_multi," | solv: ", t_solv,"  | time: ",time)
                     else
                        println(file, m," | ",s," | ",string," |",str_eq," ",alpha," | time: ",time)
                     end
                     flush(file)
                     #println(fileSolv, m," | ",s," | multi: ",t_multi," | solv: ", t_solv)
                     #elapsed = time()-t
                     println()


            end
         end
      end
      #close(fileSolv)
      #close(fileMulti)
      #display(data_frame)
   #   println("STATS: stats are not auto saved to a file")
   TIME = (time()-start)
   (minutes, seconds) = fldmod(TIME, 60)
   (hours, minutes) = fldmod(minutes, 60)
   println("FILE KEY = ",fileKey)
   @printf("%02d:%02d:%0.2f", hours, minutes, seconds)
   @printf("\nFC: %0.2f %% ",100*(perc_array[1])/perc_array[10])
   @printf("\nHALF: %0.2f %% ",100*(perc_array[2])/perc_array[10])
   @printf("\nINT: %0.2f %% ",100*(perc_array[3])/perc_array[10])
   @printf("\nMID: %0.2f %% ",100*(perc_array[4])/perc_array[10])
   @printf("\nEBM: %0.2f %% ",100*(perc_array[5])/perc_array[10])
   @printf("\nHBM: %0.2f %% ",100*(perc_array[6])/perc_array[10])
   @printf("\nGAP: %0.2f %% ",100*(perc_array[7])/perc_array[10])
   @printf("\nTRAIN: %0.2f %% ",100*(perc_array[11])/perc_array[10])
   @printf("\nTIME OUT: %0.2f %% ",100*(perc_array[8])/perc_array[10])
   @printf("\n≤: %0.2f %% ",100*(perc_array[9])/perc_array[10])

   @printf(file, "\n%02d:%02d:%0.2f", hours, minutes, seconds)
   @printf(file, "\nFC: %0.2f %% ",100*(perc_array[1])/perc_array[10])
   @printf(file, "\nHALF: %0.2f %% ",100*(perc_array[2])/perc_array[10])
   @printf(file, "\nINT: %0.2f %% ",100*(perc_array[3])/perc_array[10])
   @printf(file, "\nMID: %0.2f %% ",100*(perc_array[4])/perc_array[10])
   @printf(file, "\nEBM: %0.2f %% ",100*(perc_array[5])/perc_array[10])
   @printf(file, "\nHBM: %0.2f %% ",100*(perc_array[6])/perc_array[10])
   @printf(file, "\nGAP: %0.2f %% ",100*(perc_array[7])/perc_array[10])
   @printf(file, "\nTRAIN: %0.2f %% ",100*(perc_array[11])/perc_array[10])
   @printf(file, "\nTIME OUT: %0.2f %% ",100*(perc_array[8])/perc_array[10])
   @printf(file, "\n≤: %0.2f %% ",100*(perc_array[9])/perc_array[10])


   close(file)
   #close(file2)
   #println(perc_array[9])
   if have_file
      #str = dirname(@__FILE__)*"\\Muffins_m"*string(max_m)*"_s"*string(max_s)*"_"*string(time_limit)*"_output.csv"
      #CSV.write(dirname(@__FILE__)*"\\Muffins_m"*string(max_m)*"_s"*string(max_s)*"_output.csv", data_frame)
      #println("file created at: ",str)
   end
end
#INT(11,9)
#PROC(11,10,7//20)
#VMID(11,10,7//20,true)
#println(rand(1000:9999))
   FIND_ALL(110,3,100,rand(1000:9999),true, true)
#FIND_ALPHA(58,45)
#PROC(19,15,7//20)
#VProc(49,10,41//90)
#VMID(23,20,7//20,true, false)
#FIND_ALPHA(54,47)
#VGAPV3(58,45,22//63,1)
#PROC(67,52,127//364)
#endpoints = VGAPV3(67,52,127//364, true,true)
#VProc(67,52,127//364,60,60,endpoints)
#VMID(107,13,365//754, true)
#VGAPV3(67,21,118//259,1)
#VGAPV3(79,21,123//266,1)
#VGAPV3(101,24,359//768,1)
#FIND_ALPHA(81,13)
#FIND_ALPHA(17,10)
#VProc(97,13,85//182)
#FIND_ALPHA(96,19)
if false
VProc(67,21,41//90,300,300,0,1)
println("Our answer")
FIND_ALPHA(67,21)

VProc(83,26,77//169,300,300,0,1)
println("Our answer")
FIND_ALPHA(83,26)

VProc(69,32,31//72,300,300,0,1)
println("Our answer")
FIND_ALPHA(69,32)

VProc(91,34,197//442,300,300,0,1)
println("Our answer")
FIND_ALPHA(91,34)

VProc(110,41,311//697,300,300,0,1)
println("Our answer")
FIND_ALPHA(110,41)

VProc(71,44,31//77,300,300,0,1)
println("Our answer")
FIND_ALPHA(71,44)

VProc(101,47,727//1692,300,300,0,1)
println("Our answer")
FIND_ALPHA(101,47)

VProc(85,52,21//52,300,300,0,1)
println("Our answer")
FIND_ALPHA(85,52)

VProc(95,59,427//1062,300,300,0,1)
println("Our answer")
FIND_ALPHA(95,59)

VProc(97,60,29//72,300,300,0,1)
println("Our answer")
FIND_ALPHA(97,60)

VProc(101,62,25//62,300,300,0,1)
println("Our answer")
FIND_ALPHA(101,62)

VProc(103,63,178//441,300,300,0,1)
println("Our answer")
FIND_ALPHA(103,63)

VProc(107,66,133//330,300,300,0,1)
println("Our answer")
FIND_ALPHA(107,66)

println()
println("SOLVER ")
println()
VProc(106,15,106//225,300,300,0,1)

VProc(109,18,109//234,300,300,0,1)
end
