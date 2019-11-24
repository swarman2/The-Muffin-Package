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
using Dates
using Plots
using CSV
#using DelimitedFiles
#Input -> m amnt muffins
#     -> s amt students
#     time_limit_multi, time_limit_solv are time limits used in PROC
#     perc_arry and overlap_dic strop data
#     scott_alpha is scotts answer or 0 if using PROC
#     print_f determines if it prints with or without spaces inbetween the methods
#     stop_at_first determines if it runs all the methods or checks to see if it has right anwer everytime
function FIND_ALPHA(m,s, time_limit_multi = Inf, time_limit_solv = Inf, perc_array =0, scott_alpha = 0, overlap_dic=0,print_f = true,stop_at_first = false)

   start_time = time()
   V, sᵥ,sᵥ₋₁ = SV(m,s)
   t_multi = 0
   t_solv =0
   FC_alpha = 1
   INT_alpha = 1
   HALF_alpha = 1
   EBM_alpha = 1
   HBM_alpha = 1
   MID_alpha = 1
   GAP_alpha = 1
   TRAIN_alpha = 1
   proc_bool = true
   length_B = 0
   proc_bool = false
   min_alpha = 1
   if scott_alpha != 0 && stop_at_first
      FC_alpha = FC(m,s)
       min_alpha = FC_alpha
      proc_bool = (scott_alpha == min_alpha)
      if proc_bool!=true
         INT_alpha = INT(m,s)
         if INT_alpha < min_alpha
            min_alpha = INT_alpha
         end
         proc_bool = (scott_alpha == min_alpha)
      end
      if proc_bool!=true
         HALF_alpha = HALF(m,s)
         if HALF_alpha < min_alpha
            min_alpha = HALF_alpha
         end
         proc_bool = (scott_alpha == min_alpha)
      end
      if proc_bool!=true
         EBM_alpha = EBM(m,s)
         if EBM_alpha < min_alpha
            min_alpha = EBM_alpha
         end
         proc_bool = (scott_alpha == min_alpha)
      end
      if proc_bool!=true
         HBM_alpha = HBM(m,s)
         if HBM_alpha < min_alpha
            min_alpha = HBM_alpha
         end
         proc_bool = (scott_alpha == min_alpha)
      end
      if proc_bool!=true
         MID_alpha, MID_endpoints = MID(m,s, 1//2, true)
         if MID_alpha < min_alpha
            min_alpha = MID_alpha
         end
         proc_bool = (scott_alpha == min_alpha)
      end
      if proc_bool!=true
         GAP_alpha, GAP_endpoints = GAP(m,s,min_alpha,true)
         if GAP_alpha < min_alpha
            min_alpha = GAP_alpha
         end
         proc_bool = (scott_alpha == min_alpha)
      end
      if proc_bool != true
         TRAIN_alpha = TRAIN(m,s,min_alpha)
         if TRAIN_alpha < min_alpha
            min_alpha = TRAIN_alpha
         end
         proc_bool = (scott_alpha == min_alpha)
      end
      t_multi = 0
      t_solv = 0
   elseif s>=3 && s<=9
      #if (m==11 && s==5) || (m==7 && s==6)||(m==8 && s==7)||(m==19 && s==7)||(m==10 && s==9)||(m==11 && s==9)||(m==29&&s==9)||(m==38 && s==9)||(m==47 && s==9)
      #   min_alpha = INT(m,s)
      #   INT_alpha = min_alpha
      #else
      #   min_alpha = FC(m,s)
      #   FC_alpha = min_alpha
      #end
      FC_alpha = FC(m,s)
      INT_alpha = INT(m,s)
      HALF_alpha = HALF(m,s)
      EBM_alpha = EBM(m,s)
      HBM_alpha = HBM(m,s)
      MID_alpha, MID_endpoints = MID(m,s, 1//2, true)
      alphas = [FC_alpha INT_alpha HALF_alpha EBM_alpha HBM_alpha MID_alpha]
      min_alpha = minimum(alphas)
      if !stop_at_first
            GAP_alpha, GAP_endpoints = GAP(m,s,min_alpha, true)
            TRAIN_alpha = TRAIN(m,s,min_alpha)
         end
      proc_bool = true
   else
      FC_alpha = FC(m,s)
      INT_alpha = INT(m,s)
      HALF_alpha = HALF(m,s)
      EBM_alpha = EBM(m,s)
      HBM_alpha = HBM(m,s)
      MID_alpha, MID_endpoints = MID(m,s, 1//2, true)
      alphas = [FC_alpha INT_alpha HALF_alpha EBM_alpha HBM_alpha MID_alpha]
      min_alpha = minimum(alphas)
      if stop_at_first
         if min_alpha == 1//3
            proc_bool = one_thrd(m,s,0)
         elseif min_alpha == MID_alpha && MID_endpoints !=1
            unique!(MID_endpoints)
            proc_bool, err_mess, t_multi, t_solv, length_B = VProc(m,s,min_alpha, time_limit_solv, time_limit_multi, MID_endpoints,0)
         else
            proc_bool, err_mess, t_multi, t_solv, length_B = VProc(m,s,min_alpha, time_limit_solv, time_limit_multi, 0,0)
         end
      end
      if !proc_bool
         #if(min_alpha != 1//3)
            GAP_alpha, GAP_endpoints = GAP(m,s,min_alpha, true)

            if GAP_alpha < min_alpha
               min_alpha = GAP_alpha
            end
         #end
         if stop_at_first
            if min_alpha == 1//3
               proc_bool = one_thrd(m,s,0)
            elseif min_alpha == MID_alpha && MID_endpoints !=1
               unique!(MID_endpoints)
               proc_bool, err_mess, t_multi, t_solv, length_B = VProc(m,s,min_alpha, time_limit_solv, time_limit_multi, MID_endpoints,0)
            elseif min_alpha == GAP_alpha && GAP_endpoints !=1
               unique!(GAP_endpoints)
               println(GAP_endpoints)
               proc_bool, err_mess, t_multi, t_solv, length_B = VProc(m,s,min_alpha, time_limit_solv, time_limit_multi, GAP_endpoints,0)
            else
               proc_bool, err_mess, t_multi, t_solv, length_B = VProc(m,s,min_alpha, time_limit_solv, time_limit_multi, 0,0)
            end
         end
         if !proc_bool
            #if min_alpha!=1//3
               TRAIN_alpha = TRAIN(m,s,min_alpha)
               if TRAIN_alpha< min_alpha
                  min_alpha = TRAIN_alpha
               end
            #end
            if scott_alpha !=0
               proc_bool=(min_alpha == scott_alpha)
            elseif min_alpha == 1//3
               proc_bool = one_thrd(m,s,0)
            elseif min_alpha == MID_alpha && MID_endpoints !=1
               unique!(MID_endpoints)
               proc_bool, err_mess, t_multi, t_solv, length_B = VProc(m,s,min_alpha, time_limit_solv, time_limit_multi, MID_endpoints,0)
            elseif min_alpha == GAP_alpha && GAP_endpoints !=1
               unique!(GAP_endpoints)
               proc_bool, err_mess, t_multi, t_solv, length_B = VProc(m,s,min_alpha, time_limit_solv, time_limit_multi, GAP_endpoints,0)
            else
               proc_bool, err_mess, t_multi, t_solv, length_B = VProc(m,s,min_alpha, time_limit_solv, time_limit_multi, 0,0)
            end
         end
      else
         #limit to one method
         if FC_alpha == min_alpha
            INT_alpha = 1
            HALF_alpha = 1
            EBM_alpha = 1
            HBM_alpha =1
            MID_alpha = 1
         elseif INT_alpha == min_alpha
            FC_alpha = 1
            HALF_alpha = 1
            EBM_alpha = 1
            HBM_alpha =1
            MID_alpha = 1
         elseif HALF_alpha == min_alpha
            INT_alpha = 1
            FC_alpha = 1
            EBM_alpha = 1
            HBM_alpha =1
            MID_alpha = 1
         elseif EBM_alpha == min_alpha
            INT_alpha = 1
            HALF_alpha = 1
            FC_alpha = 1
            HBM_alpha =1
            MID_alpha = 1
         elseif HBM_alpha == min_alpha
            INT_alpha = 1
            HALF_alpha = 1
            EBM_alpha = 1
            FC_alpha =1
            MID_alpha = 1
         elseif MID_alpha == min_alpha
            INT_alpha = 1
            HALF_alpha = 1
            EBM_alpha = 1
            HBM_alpha =1
            FC_alpha = 1
         end
      end

   end
   if min_alpha ==1 && lcm(m,s)==1 || s==1
      print("   ",m,"   |   ",s,"   |  trivial   |  = ",1)
      str = " trivial "
      str_eq = " = "
      min_alpha = 1
   else
      #print("   ",m,"   |   ",s,"    |")# PROC_TIME =  ")
      str=""
      if FC_alpha == min_alpha
      #   print(" FC ")
         str = str * @sprintf(" FC ")
      elseif print_f
         str = str * @sprintf("    ")
      end

      if INT_alpha == min_alpha
       # print(" INT ")
      #  GAP_alpha = min_alpha
        str = str * @sprintf(" INT ")
     elseif print_f
        str = str * @sprintf("     ")
      end

      if  HALF_alpha == min_alpha
      #   print(" HALF ")
         str =str* @sprintf(" HALF ")
      elseif print_f
         str = str*@sprintf("      ")
      end

      if EBM_alpha == min_alpha
      #   print(" EBM ")
         str =str* @sprintf(" EBM ")
      elseif print_f
         str = str*@sprintf("     ")
      end


      if HBM_alpha == min_alpha
      #   print(" HBM ")
         str =str* @sprintf(" HBM ")
      elseif print_f
         str = str*@sprintf("     ")
      end

      if MID_alpha == min_alpha
      #   print(" MID ")
      #   GAP_alpha = min_alpha
         str =str* @sprintf(" MID " )
      elseif print_f
         str = str*@sprintf("     ")
      end

      if GAP_alpha == min_alpha
      #   print(" GAP ")
         str = str*@sprintf(" GAP ")
      elseif print_f
         str = str*@sprintf("     ")
      end

      if TRAIN_alpha == min_alpha
      #   print(" TRAIN ")
         str = str *@sprintf(" TRAIN ")
      elseif print_f
         str = str*@sprintf("       ")
      end
      str_eq = ""
      if proc_bool == true
         if overlap_dic != 0
            if FC_alpha == min_alpha
               if haskey(overlap_dic, "FC")
                  overlap_dic["FC"]=overlap_dic["FC"]+1
               else
                  overlap_dic["FC"]=1
               end
            end
            if FC_alpha == min_alpha || HALF_alpha == min_alpha
               if haskey(overlap_dic, "FC HALF")
                  overlap_dic["FC HALF"]=overlap_dic["FC HALF"]+1
               else
                  overlap_dic["FC HALF"]=1
               end
            end
            if FC_alpha == min_alpha || HALF_alpha == min_alpha || INT_alpha == min_alpha
               if haskey(overlap_dic, "FC HALF INT")
                  overlap_dic["FC HALF INT"]=overlap_dic["FC HALF INT"]+1
               else
                  overlap_dic["FC HALF INT"]=1
               end
            end
            if FC_alpha == min_alpha || HALF_alpha == min_alpha || INT_alpha == min_alpha || MID_alpha == min_alpha
               if haskey(overlap_dic, "FC HALF INT MID")
                  overlap_dic["FC HALF INT MID"]=overlap_dic["FC HALF INT MID"]+1
               else
                  overlap_dic["FC HALF INT MID"]=1
               end
            end
            if FC_alpha == min_alpha || HALF_alpha == min_alpha || INT_alpha == min_alpha || MID_alpha == min_alpha || EBM_alpha == min_alpha
               if haskey(overlap_dic, "FC HALF INT MID EBM")
                  overlap_dic["FC HALF INT MID EBM"]=overlap_dic["FC HALF INT MID EBM"]+1
               else
                  overlap_dic["FC HALF INT MID EBM"]=1
               end
            end
            if FC_alpha == min_alpha || HALF_alpha == min_alpha || INT_alpha == min_alpha || MID_alpha == min_alpha || EBM_alpha == min_alpha || HBM_alpha == min_alpha
               if haskey(overlap_dic, "FC HALF INT MID EBM HBM")
                  overlap_dic["FC HALF INT MID EBM HBM"]=overlap_dic["FC HALF INT MID EBM HBM"]+1
               else
                  overlap_dic["FC HALF INT MID EBM HBM"]=1
               end
            end
            if FC_alpha == min_alpha || HALF_alpha == min_alpha || INT_alpha == min_alpha || MID_alpha == min_alpha || EBM_alpha == min_alpha || HBM_alpha == min_alpha || GAP_alpha == min_alpha
               if haskey(overlap_dic, "FC HALF INT MID EBM HBM GAP")
                  overlap_dic["FC HALF INT MID EBM HBM GAP"]=overlap_dic["FC HALF INT MID EBM HBM GAP"]+1
               else
                  overlap_dic["FC HALF INT MID EBM HBM GAP"]=1
               end
            end
            if FC_alpha == min_alpha || HALF_alpha == min_alpha || INT_alpha == min_alpha || MID_alpha == min_alpha || EBM_alpha == min_alpha || HBM_alpha == min_alpha || GAP_alpha == min_alpha|| TRAIN_alpha == min_alpha
               if haskey(overlap_dic, "FC HALF INT MID EBM HBM GAP TRAIN")
                  overlap_dic["FC HALF INT MID EBM HBM GAP TRAIN"]=overlap_dic["FC HALF INT MID EBM HBM GAP TRAIN"]+1
               else
                  overlap_dic["FC HALF INT MID EBM HBM GAP TRAIN"]=1
               end
            end
         end
      else
         if overlap_dic !=0
            if haskey(overlap_dic, "not solved")
               overlap_dic["not solved"] = overlap_dic["not solved"]+1
            else
               overlap_dic["not solved"] =1
            end
         end
      end
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
            if HBM_alpha == min_alpha
               perc_array[6]=perc_array[6]+1
            end
            if GAP_alpha == min_alpha
               perc_array[7]=perc_array[7]+1
            end
            if TRAIN_alpha == min_alpha
               perc_array[11]=perc_array[11]+1
            end
         end
         str_eq = " = "

         #print("|  = ",numerator(min_alpha),"/",denominator(min_alpha),"")

      elseif proc_bool == -1
         str_eq = " TIME OUT "
      #   print("|  ≤ ",numerator(min_alpha),"/",denominator(min_alpha),"|   TIME OUT")
         if perc_array !=0
            perc_array[8] = perc_array[8]+1
         end
      else
         str_eq = " ≤ "
      #   print("|  ≤ ",numerator(min_alpha),"/",denominator(min_alpha))
         if perc_array !=0
            perc_array[9] = perc_array[9]+1
         end
      end
      @printf("  %-4d  |  %-4d  |  %s  | %s %5d/%-5d | %0.2f %3s",m,s,str,str_eq,numerator(min_alpha),denominator(min_alpha),time()-start_time,"sec")
   end
   if perc_array !=0
      perc_array[10]=perc_array[10]+1
   end
   return min_alpha, str, str_eq, round(time()-start_time,digits = 2), round(t_multi,digits=2),round(t_solv, digits =2), length_B
end
#FIND_ALPHA(107,13)


#outputs alpha from m = s+1 to max_m and s= min_s to max_s
#the file key is just a random number to help prevent files getting written over
# have file is a bool that determines whether a text file is created
#Scott data is a bool that determines wheter we use Scott or PROC to verify
# csv_file is a bool that determines if we output to a csv
#graphs determines if graphs are saved
#time limit is the limit on multi and solv
#m_le_s2 determines if m is upper-bounded by s squared
#stop at first decides whether to check all the methods or check PROC or SCOTT after each method (for PROC still runs all the fast methods first)
function FIND_ALL(max_m,min_s, max_s, fileKey=0, have_file = false, Scott_data = false, csv_file = false, graphs = false, time_limit = 1000, m_le_s2 = true,stop_at_first = false)
   println()
   time_limit_multi = time_limit
   time_limit_solv =  time_limit
   perc_array = zeros(11)
   dict = Dict{String, Int64}()
   overlap_dic = Dict{String,Int64}()
   numPieces_arr = Array{Int64}(undef,0)
   solveTime_arr = Array{Float64}(undef,0)
   multiTime_arr = Array{Float64}(undef,0)
   Time_arr = Array{Float64}(undef,0)
   removed = Array{String}(undef,0)
   temp_holding = Array{Int64}(undef,0)
   if Scott_data && graphs
      #s == _ mod(4)
      s0 = Array{Int64}(undef,0)
      s1 = Array{Int64}(undef,0)
      s2 = Array{Int64}(undef,0)
      s3 = Array{Int64}(undef,0)
      m0 = Array{Int64}(undef,0)
      m1 = Array{Int64}(undef,0)
      m2 = Array{Int64}(undef,0)
      m3 = Array{Int64}(undef,0)
   end
   if Scott_data &&have_file || graphs || csv_file
      c_dir =dirname(@__FILE__)*"/../DATA/"# current directory
   elseif have_file || graphs || csv_file
      c_dir =dirname(@__FILE__)*"/../DATA/"# current directory
   end
   #if have_file || graphs || csv_file
   #   mkdir(c_dir)
   #end

   if Scott_data
      df = open(dirname(@__FILE__)*"/../DATA/scott_all.txt", "r")
      if have_file
         file2 = open(c_dir * "/not_match_SCOTT_.txt","w")
         println(file2,"---------------------------------------------")
         @printf(file2,"|  %3d ≤ s ≤ %-3d    s < m < %-5d           |",min_s, max_s, max_m)
         @printf(file2,"\n|     m and s that did not match SCOTT      |")
         println(file2,"\n---------------------------------------------")
         println(file2)
      end
   end
   if  have_file
      file = open(c_dir *fileKey,"w")
      println(file,"---------------------------------------------")
      @printf(file,"|      %3d ≤ s ≤ %-3d    s < m < %-5d        |",min_s, max_s, max_m)
      if Scott_data
         @printf(file,"\n|              SCOTT verified                |")


      else
         @printf(file,"\n|               PROC verified                |")
      end
      @printf(file,"\n|              time limit = %-10d     |",time_limit)
      if m_le_s2
         @printf(file,"\n|              m ≤ s*s                       |")
      end
      println(file,"\n---------------------------------------------")
      println(file)
   end

   row_df = 1
      start=time()
      println("   m    |   s    |                   Method(s)                  |        α        |       runtime   ")
      println("------------------------------------------------------------------------------------------------------")
         str = @sprintf("   m    |   s    |                   Method(s)                  |        α        |       runtime   \n------------------------------------------------------------------------------------------------------\n")
      if have_file

         print(file,str)
      end

      if csv_file
         data_frame = DataFrame( m = Int64[], s = Int64[] , Method = String[] ,Equality = String[] , alpha = Any[] , Runtime = Any[] )
      end
      if m_le_s2
         bool=false
      else
         bool = true
      end
      for s = min_s:max_s
         for m = s+1:max_m
            if gcd(m,s)==1 && (m <= s*s || bool)
               if Scott_data

                  str = readline(df)
                  str_arr=split(str, " ")
                  while !(parse(Int64,str_arr[1])==m && parse(Int64,str_arr[2])==s)

                     try
                        if parse(Int64,str_arr[2])>s
                           s=s+1
                           m=s+1
                        end
                     catch
                        println("ERROR: NOT ENOUGH SCOTT DATA ")
                        return
                     end
                     str = readline(df)
                     str_arr=split(str, " ")

                  end
                  #str = df[Symbol("f(m,s)")][row_df]
                  str_arr2=split(str_arr[3], "/")
                  numerat = parse(Int64, str_arr2[1])
                  denominat = parse(Int64, str_arr2[2])
                  scott_alpha = numerat//denominat

               else
                  scott_alpha =0
               end
               alpha,strng,str_eq,time, t_multi, t_solv, num_pieces = FIND_ALPHA(m,s,time_limit_multi, time_limit_solv, perc_array, scott_alpha,overlap_dic,true,stop_at_first)

               if Scott_data
                  if scott_alpha != alpha
                      V,sᵥ,sᵥ₋₁=SV(m,s)
                      if V*sᵥ > (V-1)*sᵥ₋₁
                        temp_str = " V "
                     else
                        temp_str = " V-1 "
                     end
                     if have_file
                        @printf(file2, "\nm: %i  s: %i  | Scott alpha: %i/%i  | other alpha: %i/%i | split: %s" , m,s,numerat,denominat, numerator(alpha), denominator(alpha),temp_str)
                        flush(file2)
                     end
                  end
                  methods = split(strng)
                  if Scott_data && graphs
                     if !("FC" in methods) && alpha !=1
                        if s%4 == 0
                           push!(m0,m)
                           push!(s0, s)
                        elseif s%4 == 1
                           push!(m1,m)
                           push!(s1,s)
                        elseif s%4 == 2
                           push!(m2,m)
                           push!(s2,s)
                        else
                           push!(m3,m)
                           push!(s3,s)
                        end
                     end
                  end
               end
               if (time > 300 && num_pieces <40)|| (num_pieces > 100 && time < 15)
                  push!(temp_holding,m)
                  push!(temp_holding,s)
                  push!(temp_holding,num_pieces)
               end
               if graphs
                  if str_eq != " TIME OUT "
                     append!(numPieces_arr, round(num_pieces))
                     append!(Time_arr, time )
                     append!(multiTime_arr, t_multi )
                     append!(solveTime_arr, t_solv )
                  else
                     remove_string = "Removed m = "*string(m)*"  s = "*string(s)*" from graph due to time out"
                     push!(removed, remove_string)
                  end
               end
               #num = string(numerator(alphaa))
               #denom = string(denominator(alphaa))
               #alpha = num*"/"*denom
               if haskey(dict,strng)
                  dict[strng]=dict[strng]+1
               else
                  dict[strng]=1
               end
               if csv_file
                  (minutes, seconds) = fldmod(time, 60)
                  (hours, minutes) = fldmod(minutes, 60)
                  time_str = @sprintf("%02d:%02d:%0.2f", hours, minutes, seconds)
                  push!(data_frame,(m,s,strng,str_eq,alpha,time_str))

               end
               if !Scott_data && have_file
                  @printf(file,"\n  %-4d  |  %-4d  |  %-32s  | %s %5d/%-5d |   %0.2f ",m,s,strng,str_eq,numerator(alpha),denominator(alpha),time)

               elseif have_file
                  @printf(file,"\n  %-4d  |  %-4d  |  %-32s  | %s %5d/%-5d | runtime:  %0.2f",m,s,strng,str_eq,numerator(alpha),denominator(alpha),time)
               end
               if have_file
                  flush(file)
               end
               println()


            end
         end
      end
   TIME = (time()-start)
   (minutes, seconds) = fldmod(TIME, 60)
   (hours, minutes) = fldmod(minutes, 60)
   if Scott_data && have_file
      close(file2)
   end
   if Scott_data && graphs
      default(size = (2000,1000))
      default(legend=false)
      p1 = plot(s0, m0, seriestype =:scatter,xlims = (0,max_s+2), ylims = maximum(m0)+2,xlabel = "s",ylabel = "non FC m", xticks = 0:5:max_s, yticks = 0:5:max_m, title = "s==0 mod 4")
      savefig(p1,c_dir*"/Plot_(s==0 mod 4).png")
      p2 =plot(s1,m1, seriestype =:scatter,xlims = (0,max_s+2), ylims = maximum(m1)+2,xlabel = "s",ylabel = "non FC m",xticks = 0:5:max_s, yticks = 0:5:max_m, title = "s==1 mod 4")
      savefig(p2,c_dir*"/Plot_(s==1 mod 4).png")
      p3 =plot(s2,m2,seriestype =:scatter,xlims = (0,max_s+2), ylims = maximum(m2)+2,xlabel = "s",ylabel = "non FC m",xticks = 0:5:max_s, yticks = 0:5:max_m, title = "s==2 mod 4")
      savefig(p3,c_dir*"/Plot_(s==2 mod 4).png")
      p4 = plot(s3,m3,seriestype =:scatter,xlims = (0,max_s+2), ylims = maximum(m3)+2,xlabel = "s",ylabel = "non FC m",xticks = 0:5:max_s, yticks = 0:5:max_m, title = "s==3 mod 4")
      savefig(p4,c_dir*"/Plot_(s==3 mod 4).png")
      display(plot(p1,p2,p3,p4))
      println("Graphs saved to: ")
      println(c_dir*"/Plot_(s==0 mod 4).png")
      println(c_dir*"/Plot_(s==1 mod 4).png")
      println(c_dir*"/Plot_(s==2 mod 4).png")
      println(c_dir*"/Plot_(s==3 mod 4).png")
   end
   println()
   println("---------------- STATS ----------------")
   @printf("%s %02d:%02d:%0.2f","Total time: ", hours, minutes, seconds)
   println()
   println()
   print("Amount of times each method produced the correct alpha")
   @printf("\nFC: %i/%i ==> %0.2f %% ",perc_array[1],perc_array[10],100*(perc_array[1])/perc_array[10])
   @printf("\nHALF: %i/%i ==> %0.2f %% ",perc_array[2],perc_array[10],100*(perc_array[2])/perc_array[10])
   @printf("\nINT: %i/%i == >%0.2f %% ",perc_array[3],perc_array[10],100*(perc_array[3])/perc_array[10])
   @printf("\nMID: %i/%i ==> %0.2f %% ",perc_array[4],perc_array[10],100*(perc_array[4])/perc_array[10])
   @printf("\nEBM: %i/%i ==> %0.2f %% ",perc_array[5],perc_array[10],100*(perc_array[5])/perc_array[10])
   @printf("\nHBM: %i/%i ==> %0.2f %% ",perc_array[6],perc_array[10],100*(perc_array[6])/perc_array[10])
   @printf("\nGAP:%i/%i ==>  %0.2f %% ",perc_array[7],perc_array[10],100*(perc_array[7])/perc_array[10])
   @printf("\nTRAIN: %i/%i ==> %0.2f %% ",perc_array[11],perc_array[10],100*(perc_array[11])/perc_array[10])
   println()
   println()
   println("Amount of times the correct alpha was not found")
   @printf("TIME OUT [%i sec]: %i/%i ==> %0.2f %% ",time_limit,perc_array[8],perc_array[10],100*(perc_array[8])/perc_array[10])
   @printf("\nIncorrect upper-bound: %i/%i ==> %0.2f %% ",perc_array[9],perc_array[10],100*(perc_array[9])/perc_array[10])
   println()
   println()
   if have_file
      println(file)
      println(file)
      println(file,"---------------- STATS ----------------")
      println(file)
      @printf(file,"%s %02d:%02d:%0.2f","Total time: ", hours, minutes, seconds)
      println(file)
      println(file)
      println(file,"Amount of times each method produced the correct alpha")
      @printf(file,"\nFC: %i/%i ==> %0.2f %% ",perc_array[1],perc_array[10],100*(perc_array[1])/perc_array[10])
      @printf(file,"\nHALF: %i/%i ==> %0.2f %% ",perc_array[2],perc_array[10],100*(perc_array[2])/perc_array[10])
      @printf(file,"\nINT: %i/%i == >%0.2f %% ",perc_array[3],perc_array[10],100*(perc_array[3])/perc_array[10])
      @printf(file,"\nMID: %i/%i ==> %0.2f %% ",perc_array[4],perc_array[10],100*(perc_array[4])/perc_array[10])
      @printf(file,"\nEBM: %i/%i ==> %0.2f %% ",perc_array[5],perc_array[10],100*(perc_array[5])/perc_array[10])
      @printf(file,"\nHBM: %i/%i ==> %0.2f %% ",perc_array[6],perc_array[10],100*(perc_array[6])/perc_array[10])
      @printf(file,"\nGAP:%i/%i ==>  %0.2f %% ",perc_array[7],perc_array[10],100*(perc_array[7])/perc_array[10])
      @printf(file,"\nTRAIN: %i/%i ==> %0.2f %% ",perc_array[11],perc_array[10],100*(perc_array[11])/perc_array[10])
      println(file)
      println(file)
      println(file,"Amount of times the correct alpha was not found")
      @printf(file,"TIME OUT [ %i sec]: %i/%i ==> %0.2f %% ",time_limit, perc_array[8],perc_array[10],100*(perc_array[8])/perc_array[10])
      @printf(file,"\nIncorrect upper-bound: %i/%i ==> %0.2f %% ",perc_array[9],perc_array[10],100*(perc_array[9])/perc_array[10])
      println(file)
      println(file)
      println(file,"Amount of times methods overlapped: ")
   end
   println()
   println("Amount of times methods overlapped: ")

   for k in keys(dict)
      println(k," ==> ",dict[k],"/",perc_array[10]," = ",round(100*dict[k]/perc_array[10],digits=2),"%")
      if have_file
         println(file, k," ==> ",dict[k],"/",perc_array[10]," = ",round(100*dict[k]/perc_array[10],digits = 2),"%")
      end
   end
   if have_file
      println(file)
      println(file,"Percentage of correct alphas found adding in methods one by one:")
   end
   println()
   println("Percentage of correct alphas found adding in methods one by one:")
   for k in sort(collect(keys(overlap_dic)))

      println(k,"    =>    ",overlap_dic[k],"/",perc_array[10]," = ", round(100*overlap_dic[k]/perc_array[10],digits =2),"%")
      if have_file
         println(file,k,"    =>    ",overlap_dic[k],"/",perc_array[10]," = ", round(100*overlap_dic[k]/perc_array[10],digits =2),"%")
         println(file)
      end
      println()
   end
   if have_file
      close(file)
   end

   println()
   println()
   if have_file
      println("txt file created at: Data/", fileKey,".txt")
   end
   if csv_file
      str =c_dir*"/data.csv"
      CSV.write(str, data_frame)
      println()
      println("csv file created at: ",str)
   end
   #println("-------------------------TEMP DATA-----------------------")
   #i=1
   #while i<length(temp_holding)
   #   println("m:  ",temp_holding[i],"  s: ",temp_holding[i+1],"  num pieces: ",temp_holding[i+2])
   #   i = i+3
   #end
   if !Scott_data && graphs
      for i = 1:length(removed)
         println(removed[i])
      end
      p1 = plot(numPieces_arr, Time_arr, seriestype =:scatter, title = "Solver Variables vs Total Time")
      savefig(p1,c_dir*"/TotalTime.png")
      p2 = plot(numPieces_arr, multiTime_arr, seriestype =:scatter, title = "Solver Variables vs Multi Time")
      savefig(p2,c_dir*"/MultiTime.png")
      p3 =plot(numPieces_arr, solveTime_arr,seriestype =:scatter, title = "Solver Variables vs Solve Time")
      savefig(p3,c_dir*"/SolvTime.png")
      display(plot(p1,p2,p3))
      println("Graphs saved to: ")
      println(c_dir*"/TotalTime.png")
      println(c_dir*"/MultiTime.png")
      println(c_dir*"/SolvTime.png")
   end
end
