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
function FIND_ALPHA(m,s, time_limit_multi = Inf, time_limit_solv = Inf, perc_array =0, scott_alpha = 0, overlap_dic=0,print_f = true)

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
   min_alpha = 1
   if scott_alpha != 0
      FC_alpha = FC(m,s)
       min_alpha = FC_alpha
      proc_bool = (scott_alpha == min_alpha)
      if proc_bool!=true
         INT_alpha = INT(m,s)
         min_alpha = INT_alpha
         proc_bool = (scott_alpha == min_alpha)
      end
      if proc_bool!=true
         HALF_alpha = HALF(m,s)
         min_alpha = HALF_alpha
         proc_bool = (scott_alpha == min_alpha)
      end
      if proc_bool!=true
         EBM_alpha = EBM(m,s)
         min_alpha = EBM_alpha
         proc_bool = (scott_alpha == min_alpha)
      end
      if proc_bool!=true
         HBM_alpha = HBM(m,s)
         min_alpha = HBM_alpha
         proc_bool = (scott_alpha == min_alpha)
      end
      if proc_bool!=true
         MID_alpha, MID_endpoints = MID(m,s, 1//2, true)
         min_alpha = MID_alpha
         proc_bool = (scott_alpha == min_alpha)
      end
      if proc_bool!=true
         GAP_alpha, GAP_endpoints = GAP(m,s,min_alpha,true)
         min_alpha = GAP_alpha
         proc_bool = (scott_alpha == min_alpha)
      end
      if proc_bool != true
         #println(min_alpha)
         TRAIN_alpha = TRAIN(m,s,min_alpha)
         min_alpha = TRAIN_alpha
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
   else
      FC_alpha = FC(m,s)
      INT_alpha = INT(m,s)
      HALF_alpha = HALF(m,s)
      EBM_alpha = EBM(m,s)
      HBM_alpha = HBM(m,s)
      MID_alpha, MID_endpoints = MID(m,s, 1//2, true)
      alphas = [FC_alpha INT_alpha HALF_alpha EBM_alpha HBM_alpha MID_alpha]
      min_alpha = minimum(alphas)
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


      if min_alpha == 1//3
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
   if min_alpha ==1 && lcm(m,s)==1 || s==1
      print("   ",m,"   |   ",s,"   |  trivial   |  = ",1)
      str = " trivial "
      str_eq = " = "
      min_alpha = 1
   else
      #print("   ",m,"   |   ",s,"    |")# PROC_TIME =  ")
      str=""
      if FC_alpha <= min_alpha
      #   print(" FC ")
         str = str * @sprintf(" FC ")
      elseif print_f
         str = str * @sprintf("    ")
      end

      if INT_alpha <= min_alpha
       # print(" INT ")
      #  GAP_alpha = min_alpha
        str = str * @sprintf(" INT ")
     elseif print_f
        str = str * @sprintf("     ")
      end

      if  HALF_alpha <= min_alpha
      #   print(" HALF ")
         str =str* @sprintf(" HALF ")
      elseif print_f
         str = str*@sprintf("      ")
      end

      if EBM_alpha <= min_alpha
      #   print(" EBM ")
         str =str* @sprintf(" EBM ")
      elseif print_f
         str = str*@sprintf("     ")
      end


      if HBM_alpha <= min_alpha
      #   print(" HBM ")
         str =str* @sprintf(" HBM ")
      elseif print_f
         str = str*@sprintf("     ")
      end

      if MID_alpha <= min_alpha
      #   print(" MID ")
      #   GAP_alpha = min_alpha
         str =str* @sprintf(" MID " )
      elseif print_f
         str = str*@sprintf("     ")
      end

      if GAP_alpha <= min_alpha
      #   print(" GAP ")
         str = str*@sprintf(" GAP ")
      elseif print_f
         str = str*@sprintf("     ")
      end

      if TRAIN_alpha <= min_alpha
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
      @printf("  %-4d  |  %-4d  |  %-32s  | %s %5d/%-5d | %0.2f %3s",m,s,str,str_eq,numerator(min_alpha),denominator(min_alpha),time()-start_time,"sec")
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
function FIND_ALL(max_m,min_s, max_s, fileKey=0, have_file = false, Scott_data = false, csv_file = false, graphs = false, time_limit = 1000, m_le_s2 = true)
   println()
   day = string(Dates.today())
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

   if Scott_data
      df = open("D:\\Documents\\Muffins\\scott_all.txt", "r")
      file2 = open(dirname(@__FILE__)*"\\..\\DATA\\"*day*"_not_match_SCOTT_"* string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".txt","w")
   end
   if Scott_data && have_file
      file_name = dirname(@__FILE__)*"\\..\\DATA\\"*day*"_ScottVerified_"* string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".txt"
      file = open(dirname(@__FILE__)*"\\..\\DATA\\"*day*"_ScottVerified_"* string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".txt","w")
   elseif have_file
      file_name = dirname(@__FILE__)*"\\..\\DATA\\"*day*"_PROCVerified_"* string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".txt"
      file = open(dirname(@__FILE__)*"\\..\\DATA\\"*day*"_PROCVerified_"* string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".txt","w")
   end
   if !Scott_data
      all =  open(dirname(@__FILE__)*"\\..\\DATA\\BIGRUN_all_"*string(fileKey)*".tex","w")
      non_fc =  open(dirname(@__FILE__)*"\\..\\DATA\\BIGRUN_non_FC_"*string(fileKey)*".tex","w")
   end
   row_df = 1
      start=time()
      println("   m    |   s    |                   Method(s)                  |        α        |       runtime   ")
      println("------------------------------------------------------------------------------------------------------")
         str = @sprintf("   m    |   s    |                   Method(s)                  |        α        |       runtime   \n------------------------------------------------------------------------------------------------------\n")
      if have_file

         print(file,str)
      end
      print(all,"  m    |   s    |                   Method(s)                  |        α        \n------------------------------------------------------------------------------------")
      print(non_fc,"  m    |   s    |                   Method(s)                  |        α         \n-------------------------------------------------------------------------------------")
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
               alpha,strng,str_eq,time, t_multi, t_solv, num_pieces = FIND_ALPHA(m,s,time_limit_multi, time_limit_solv, perc_array, scott_alpha,overlap_dic)
            #   println(num_pieces)
               if !Scott_data
                  @printf(all,"   \n%-4d  |  %-4d  |   %-32s | %i/%i",m,s,strng,numerator(alpha),denominator(alpha))
                  methods = split(strng)
                  if !("FC" in methods)
                     @printf(non_fc,"   \n%-4d  |  %-4d  |   %-32s | %i/%i",m,s,strng,numerator(alpha),denominator(alpha))
                  end
               end
               if Scott_data
                  if scott_alpha != alpha
                      V,sᵥ,sᵥ₋₁=SV(m,s)
                      if V*sᵥ > (V-1)*sᵥ₋₁
                        temp_str = " V "
                     else
                        temp_str = " V-1 "
                     end
                     @printf(file2, "\nm: %i  s: %i  | Scott alpha: %i/%i  | other alpha: %i/%i | split: %s" , m,s,numerat,denominat, numerator(alpha), denominator(alpha),temp_str)
                     flush(file2)
                  end
                  methods = split(strng)
                  if Scott_data && graphs
                     if !("FC" in methods) && alpha !=1
                        if s%4 == 0
                           if s in s0
                              m0[end] = m0[end]+1
                           else
                              push!(m0,1)
                              push!(s0, s)
                           end
                        elseif s%4 == 1
                           if s in s1
                              m1[end]=m1[end]+1
                           else
                              push!(m1,1)
                              push!(s1,s)
                           end
                        elseif s%4 == 2
                           if s in s2
                              m2[end]=m2[end]+1
                           else
                              push!(m2,1)
                              push!(s2,s)
                           end
                        else
                           if s in s3
                              m3[end]=m3[end]+1
                           else
                              push!(m3,1)
                              push!(s3,s)
                           end
                        end
                     end
                  end
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
                  @printf(file,"\n  %-4d  |  %-4d  |  %-32s  | %s %5d/%-5d |  multi: %0.2f  -  solv: %0.2f  - runtime:  %0.2f ",m,s,strng,str_eq,numerator(alpha),denominator(alpha),t_multi, t_solv,time)

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
   if !Scott_data
      close(all)
      close(non_fc)
   end
   if Scott_data && graphs
      p1 = plot(s0, m0, seriestype =:scatter, title = "Number of non FC per s (s==0 mod 4)")
      savefig(p1,dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_(s==0 mod 4)"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".png")
      p2 =plot(s1,m1, seriestype =:scatter, title = "Number of non FC per s (s==1 mod 4)")
      savefig(p2,dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_(s==1 mod 4)"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".png")
      p3 =plot(s2,m2,seriestype =:scatter, title = "Number of non FC per s (s==2 mod 4)")
      savefig(p3,dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_(s==2 mod 4)"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".png")
      p4 = plot(s3,m3,seriestype =:scatter, title = "Number of non FC per s (s==3 mod 4)")
      savefig(p4,dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_(s==3 mod 4)"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".png")
      display(plot(p1,p2,p3,p4))
      println("Graphs saved to: ")
      println(dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_(s==0 mod 4)"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".png")
      println(dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_(s==1 mod 4)"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".png")
      println(dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_(s==2 mod 4)"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".png")
      println(dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_(s==3 mod 4)"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".png")
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
      println(file,"Amount of times each method produced the correct alpha")
      @printf(file,"%s %02d:%02d:%0.2f","Total time: ", hours, minutes, seconds)
      @printf(file,"%02d:%02d:%0.2f", hours, minutes, seconds)
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
      println("txt file created at: ",file_name)
   end
   if csv_file && Scott_data
      str = dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Muffins_ScottVerified_"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".csv"
      CSV.write(str, data_frame)
      println()
      println("csv file created at: ",str)
   elseif csv_file
      str = dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Muffins_PROCVerified_"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*".csv"
      CSV.write(str, data_frame)
      println()
      println("file created at: ",str)
   end
   if !Scott_data && graphs
      for i = 1:length(removed)
         println(removed[i])
      end
      p1 = plot(numPieces_arr, Time_arr, seriestype =:scatter, title = "Pieces vs Total Time")
      savefig(p1,dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*"_TotalTime.png")
      p2 = plot(numPieces_arr, multiTime_arr, seriestype =:scatter, title = "Pieces vs Multi Time")
      savefig(p2,dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*"_MultiTime.png")
      p3 =plot(numPieces_arr, solveTime_arr,seriestype =:scatter, title = "Pieces vs Solve Time")
      savefig(p3,dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*"_SolvTime.png")
      display(plot(p1,p2,p3))
      println("Graphs saved to: ")
      println(dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*"_TotalTime.png")
      println(dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*"_MultiTime.png")
      println(dirname(@__FILE__)*"\\..\\DATA\\"*day*"_Plot_"*string(max_m)*"_"*string(max_s)*"_"*string(fileKey)*"_SolvTime.png")
   end
end
