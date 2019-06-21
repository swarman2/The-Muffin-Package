include("Combinatorics.jl-master/src/combinations.jl")
include("Combinatorics.jl-master/src/partitions.jl")

using JuMP
using GLPK


function PROC(m,s,alphaa)

println("Procedure for f(",m,", ",s,") ≥ ",alphaa)

  num=numerator(alphaa)
  denom=denominator(alphaa)
  denom=lcm(s,denom)
  lower_bound = alphaa
  upper_bound = 1 - alphaa

  upper_bound_num = upper_bound * denom
  lower_bound_num = lower_bound * denom

  if upper_bound_num - lower_bound_num < 2
    B = collect(lower_bound_num:1:upper_bound_num + 1)
  else
    B = collect(lower_bound_num:1:upper_bound_num)
  end

  A = Array{Int64,1}(undef,0)
  min_length=floor(((denom))/maximum(B))
  max_length=ceil(((denom))/minimum(B))
  comp_muf = (collect(partitions(Int64(denom),Int64(min_length))))

  for i=min_length+1:max_length
      append!(comp_muf,(collect(partitions(Int64(denom),Int64(i)))))
  end

  for i = 1:length(comp_muf)
      for j = 1:length(comp_muf[i])
          vec_1 = comp_muf[i]
          filter_1 = filter(x -> (x >= minimum(B) && x <= maximum(B)), vec_1)

          if length(filter_1) == length(vec_1)

             for low_upp = lower_bound_num:upper_bound_num
                 count_elem = count(i -> (i==low_upp),filter_1)
                 append!(A,count_elem)
             end
         end
      end
  end


  shape_1 = Int64(length(A)/length(B))
  mat_1 = (reshape(A,length(B),shape_1))
  mat_1=unique(mat_1,dims=2)

  student_muiltiset_sum = ((m/s)*denom)

  AA = Array{Int64,1}(undef,0)

  min_length=floor(((m/s)*denom)/maximum(B))
  max_length=ceil(((m/s)*denom)/minimum(B))
  comp_stu = (collect(partitions(Int64((m/s)*denom),Int64(min_length))))

  for i=min_length+1:max_length
      append!(comp_stu,(collect(partitions(Int64((m/s)*denom),Int64(i)))))
  end

  for i = 1:length(comp_stu)
     for j = 1:length(comp_stu[i])
         vec_1 = comp_stu[i]
         filter_1 = filter(x -> (x >= minimum(B) && x <= maximum(B)), vec_1)

         if length(filter_1) == length(vec_1)
            for low_upp = lower_bound_num:upper_bound_num
                count_elem = count(i -> (i==low_upp),filter_1)
                append!(AA,count_elem)

            end
         end
     end
  end

  shape_2 = Int64(length(AA)/length(B))
  mat_2 = (reshape(AA,length(B),shape_2))
  mat_2=unique(mat_2,dims=2)


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

PROC(5,3,5//12)
PROC(13,5,13//30)
PROC(17,15,7//20)
