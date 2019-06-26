include("src/combinations.jl")
include("src/partitions.jl")
include("helper_functions.jl")
using JuMP
using GLPK

function VProc(m,s,alphaa)

    num=numerator(alphaa)
    denom=denominator(alphaa)
    denom=lcm(s,denom)
    lower_bound = alphaa
    upper__bound = (1 - alphaa)

    upper_bound_num = upper__bound * denom
    lower_bound_num = lower_bound * denom


      B = collect(lower_bound_num:1:upper_bound_num)

    A = Array{Int64,1}(undef,0)


    vec_1 = Multiset(B,Int64(denom),2)
    #println(vec_1,"  ",typeof(vec_1))
    vec_2 = Multiset(B,Int64(denom),3)
    #prinlnt("1.5")
    if(vec_1 !=0 && vec_2 !=0)
        for i=1:length(vec_2)
            push!(vec_1,vec_2[i])
        end
    elseif vec_1 == 0
        vec_1 = vec_2
    end
    #println(vec_1)
    for i=1:length(vec_1)
        for low_upp = lower_bound_num:upper_bound_num
          count_elem = count(i -> (i==low_upp),vec_1[i])
          append!(A,count_elem)
        end
    end
    shape_1 = Int64(length(A)/(length(B)))
    mat_1 = (reshape(A,(length(B)),shape_1))
    mat_1=unique(mat_1,dims=2)


    student_muiltiset_sum = ((m//s)*denom)

    AA = Array{Int64,1}(undef,0)

    min_length=floor(((m/s)*denom)/maximum(B))
    max_length=ceil(((m/s)*denom)/minimum(B))
    V=Int64(ceil(2m/s))


    vec_1 = Multiset(B,Int64((m//s)*denom),V)
    #println(vec_1,"  ",typeof(vec_1))
    vec_2 = Multiset(B,Int64((m//s)*denom),V-1)
    #prinlnt("1.5")
    if(vec_1 !=0 && vec_2 !=0)
        for i=1:length(vec_2)
            push!(vec_1,vec_2[i])
        end
    elseif vec_1 == 0
        vec_1 = vec_2
    end
    #println(vec_1)
    for i=1:length(vec_1)
        for low_upp = lower_bound_num:upper_bound_num
          count_elem = count(i -> (i==low_upp),vec_1[i])
          append!(AA,count_elem)
        end
    end
    #println("2")
    shape_2 = Int64(length(AA)/(length(B)))
    mat_2 = (reshape(AA,(length(B)),shape_2))
    mat_2=unique(mat_2,dims=2)
    #println(mat_2)

    final_mat = [mat_1 (-mat_2)]

    row_1, col_1 = size(mat_1)
    row_2, col_2 = size(mat_2)

    a = (zeros(Int64,col_1+col_2))
    b = (ones(Int64,col_1+col_2))

    for i = 1:col_1
        a[i] = 1
        b[i] = 0
    end

    final_mat = [final_mat; a'; b']

    c = zeros(Int64,row_1)
    final_mat_2 = [c;m;s]

    model=Model(with_optimizer(GLPK.Optimizer))
    @variable(model, x[i=1:col_1+col_2],Int)
    @constraint(model,con_1,x.>=0)
    @constraint(model,con_2,final_mat*x .==final_mat_2)
    optimize!(model)

    if(!has_values(model))
        return false
    else
        return true
    end
end


function PROC(m,s,alphaa)
    println("Procedure for f(",m,", ",s,") ≥ ",alphaa)
    if m%s == 0
        println("s divides m, give all students whole muffin(s)")
      return 1
    end
      num=numerator(alphaa)
      denom=denominator(alphaa)
      denom=lcm(s,denom)
      lower_bound = alphaa
      upper__bound = (1 - alphaa)

      upper_bound_num = upper__bound * denom
      lower_bound_num = lower_bound * denom


      B = collect(lower_bound_num:1:upper_bound_num)

      A = Array{Int64,1}(undef,0)

      vec_1 = Multiset(B,Int64(denom),2)
      #println(vec_1,"  ",typeof(vec_1))
      vec_2 = Multiset(B,Int64(denom),3)
      #prinlnt("1.5")
      if(vec_1 !=0 && vec_2 !=0)
          for i=1:length(vec_2)
              push!(vec_1,vec_2[i])
          end
      elseif vec_1 == 0
          vec_1 = vec_2
      end
      #TODO fix if vec_2 is zero
      #println(vec_1)
      for i=1:length(vec_1)
          for low_upp = lower_bound_num:upper_bound_num
            count_elem = count(i -> (i==low_upp),vec_1[i])
            append!(A,count_elem)
          end
      end

      shape_1 = Int64(length(A)/(length(B)))
      mat_1 = (reshape(A,(length(B)),shape_1))
      mat_1=unique(mat_1,dims=2)


      student_muiltiset_sum = ((m//s)*denom)

      AA = Array{Int64,1}(undef,0)

      min_length=floor(((m/s)*denom)/maximum(B))
      max_length=ceil(((m/s)*denom)/minimum(B))
     #println(min_length,"  ",max_length)
     V=Int64(ceil(2m/s))

      vec_1 = Multiset(B,Int64((m//s)*denom),V)
      #println(vec_1,"  ",typeof(vec_1))
      vec_2 = Multiset(B,Int64((m//s)*denom),V-1)
      #prinlnt("1.5")
      if(vec_1 !=0 && vec_2 !=0)
          for i=1:length(vec_2)
              push!(vec_1,vec_2[i])
          end
      elseif vec_1 == 0
          vec_1 = vec_2
      end
      #println(vec_1)
      for i=1:length(vec_1)
          for low_upp = lower_bound_num:upper_bound_num
            count_elem = count(i -> (i==low_upp),vec_1[i])
            append!(AA,count_elem)
          end
      end
    #println("test")
      shape_2 = Int64(length(AA)/(length(B)))
      mat_2 = (reshape(AA,(length(B)),shape_2))
      mat_2=unique(mat_2,dims=2)
      display(mat_2)

      final_mat = [mat_1 (-mat_2)]

      row_1, col_1 = size(mat_1)
      row_2, col_2 = size(mat_2)

      a = (zeros(Int64,col_1+col_2))
      b = (ones(Int64,col_1+col_2))

      for i = 1:col_1
          a[i] = 1
          b[i] = 0
      end

      final_mat = [final_mat; a'; b']

      c = zeros(Int64,row_1)
      final_mat_2 = [c;m;s]
      #display(final_mat)
      #println("*********")
      #display(final_mat_2)
    #println("test2")
      model=Model(with_optimizer(GLPK.Optimizer))
      @variable(model, x[i=1:col_1+col_2],Int)
      @constraint(model,con_1,x.>=0)
      @constraint(model,con_2,final_mat*x .==final_mat_2)
      optimize!(model)
    #println("test3")
      if(!has_values(model))
          println("No procedure found")
          return false
      end
      B=B//denom

      for i=1:col_1
          if(value.(x)[i]!=0)
          print("Cut ",Int64(value.(x)[i])," muffins {  ")

          for j=1:row_1
              while(mat_1[j,i]>0)
              print(B[j],"  ")
              mat_1[j,i]=mat_1[j,i]-1
          end
      end
      println("}")
      end
      end

      for i=1:col_2
          if(value.(x)[i+col_1]!=0)
          print("Give ",Int64(value.(x)[i+col_1])," students {  ")

          for j=1:row_1
              while(mat_2[j,i]>0)
              print(B[j],"  ")
              mat_2[j,i]=mat_2[j,i]-1
          end
      end
      println("}")
      end
      end
     println()
     println()
end

#for k=1:20
#    println("*********  k = ",k,"  ********* ")
#    PROC(3k+3,3k+2,((k+1)//(3k+2)))
#end
#VProc(61,19,313//684)

#VProc(5,3,5//12)
#VProc(13,5,13//30)
#VProc(17,15,7//20)
