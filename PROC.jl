#include("src/combinations.jl")
#include("src/partitions.jl")
include("helper_functions.jl")
using Cbc
using GLPK
using JuMP
using Printf
function VProc(m,s,alphaa, time_limit_solv = 60, time_limit_multi =60)
#Based on the book algorithm page 42
    start_time=time()
    #denom is the least common mulitple between s and the
    #denominator of alpha (it is the denominator of the
    #fractions in the set of possible piece sizes [later called B])
    num=numerator(alphaa)
    denom=denominator(alphaa)
    denom=lcm(s,denom)

    #lower bound and upper bound are the min and max piece sizes
    lower_bound = alphaa
    upper__bound = (1 - alphaa)

    #clear the fractions in the piece size array (B)
    upper_bound_num = upper__bound * denom
    lower_bound_num = lower_bound * denom

    #B is the set of possible piece sizes
    B = collect(lower_bound_num:1:upper_bound_num)

    #use Mulitset function (in helper_functions) to find
    #the possible mulitsets of piece sizes that sum to denom (1 whole muffin)
    #of sizes 2 (muffins cut into two pieces) or 3 (muffins cut into 3 pieces)
    vec_1 = Multiset(B,Int64(denom),2, time_limit_multi)
    if vec_1 == "time out"
        println("MULTISET 1 TIME OUT:  length(B) = ",length(B),"  sets sum to: ",Int64(denom),"  sets of size: ",2)
        str = @sprintf("MULTISET 2 TIME OUT:  length(B) = %i   sets sum to: %i  sets of size: %i",length(B),Int64(denom),2)
        return -1,str
    end
    vec_2 = 0
    if alphaa<1//3
        vec_2 = Multiset(B,Int64(denom),3,time_limit_multi)
        if vec_2 == "time out"
            println("MULTISET 2 TIME OUT:  length(B) = ",length(B),"  sets sum to: ",Int64(denom),"  sets of size: ",3)
            str = @sprintf("MULTISET 2 TIME OUT:  length(B) = %i   sets sum to: %i  sets of size: %i",length(B),Int64(denom),3)
            return -1,str
        end
    end
    #Mulitset returns an array of arrays of the possible mulitsets
    #or 0 if no mulitsets found, combine the arrays of size 2 and 3
    #into one array of arrays
    if(vec_1 !=0 && vec_2 !=0)
        for i=1:length(vec_2)
            push!(vec_1,vec_2[i])
        end
    elseif vec_1 == 0
        vec_1 = vec_2
    elseif vec_1==0 && vec_2==0
        return false
    end
    if alphaa == 1//3
        push!(vec_1, [denom/3; denom/3; denom/3])
    end
    #display(vec_1)
    #count how many times an element of B occurs in a vector
    # and convert to different form (ex B = {5 6 7} V={6 6}-> newV {0 2 0})
    A = Array{Int64,1}(undef,0)
    for i=1:length(vec_1)
        for low_upp = lower_bound_num:upper_bound_num
          count_elem = count(i -> (i==low_upp),vec_1[i])
          append!(A,count_elem)
        end
    end
    #reshape inot a matrix
    shape_1 = Int64(length(A)/(length(B)))
    mat_1 = (reshape(A,(length(B)),shape_1))
    mat_1=unique(mat_1,dims=2)

    V=Int64(ceil(2m/s))
    #use Mulitset function (in helper_functions) to find
    #the possible mulitsets of studdent distributions that sum to denom *m/s
#    println("length of B = ",length(B),"  m//s *denom = ",Int64((m//s)*denom))
    vec_1 = Multiset(B,Int64((m//s)*denom),V,time_limit_multi)
    if vec_1 == "time out"
        println("MULTISET 3 TIME OUT:  length(B) = ",length(B),"  sets sum to: ",Int64((m//s)*denom),"  sets of size: ",V)
        str = @sprintf("MULTISET 3 TIME OUT:  length(B) = %i   sets sum to: %i  sets of size: %i",length(B),Int64((m//s)*denom),V)
        return -1,str
    end
    vec_2 = Multiset(B,Int64((m//s)*denom),V-1,time_limit_multi)
    if vec_2 == "time out"
        println("MULTISET 4 TIME OUT:  length(B) = ",length(B),"  sets sum to: ",Int64((m//s)*denom),"  sets of size: ",V-1)
        str = @sprintf("MULTISET 4 TIME OUT:  length(B) = %i   sets sum to: %i  sets of size: %i",length(B),Int64((m//s)*denom),V-1)
        return -1,str
    end
#    println("finished mulitsets")
    if(vec_1 !=0 && vec_2 !=0)
        for i=1:length(vec_2)
            push!(vec_1,vec_2[i])
        end
    elseif vec_1 == 0
        vec_1 = vec_2
    elseif vec_1==0 && vec_2==0
        return false
    end

    #count how many occurances and convert
    AA = Array{Int64,1}(undef,0)
    for i=1:length(vec_1)
        for low_upp = lower_bound_num:upper_bound_num
          count_elem = count(i -> (i==low_upp),vec_1[i])
          append!(AA,count_elem)
        end
    end

    #reshape
    shape_2 = Int64(length(AA)/(length(B)))
    mat_2 = (reshape(AA,(length(B)),shape_2))
    mat_2=unique(mat_2,dims=2)

    #this is the system on page 42 with students moved to the
    #left hand side (so it can be set = to 0)
    final_mat = [mat_1 (-mat_2)]

    #add two rows of 0's and 1's
    #(one with 1's under muffin and one with 1's under students)
    row_1, col_1 = size(mat_1)
    row_2, col_2 = size(mat_2)
    a = (zeros(Int64,col_1+col_2))
    b = (ones(Int64,col_1+col_2))
    for i = 1:col_1
        a[i] = 1
        b[i] = 0
    end
    final_mat = [final_mat; a'; b']

    #set it equal to final_mat_2 (0's then num muff then num students)
    c = zeros(Int64,row_1)
    final_mat_2 = [c;m;s]
    #println("time: ",time()-start_time)
    #solve the systems
    model=Model(with_optimizer(Cbc.Optimizer,logLevel=0, seconds = time_limit_solv))
    @variable(model, x[i=1:col_1+col_2],Int)
    @constraint(model,con_1,x.>=0)
    @constraint(model,con_2,final_mat*x .==final_mat_2)
#    display(final_mat)
#    println("*******************")
#    display(final_mat_2)
    optimize!(model)
    #println("time2: ",time()-start_time)
    #return true if it has a solution, false otherwise
    if (termination_status(model)==MOI.TIME_LIMIT)
        println("SOLVER TIME OUT  dim(A) = ",size(final_mat))
        row,col = size(final_mat)
        str = @sprintf("SOLVER TIME OUT  dim(A) = (%i, %i) ",row,col)
        return -1,str
    end
    if(termination_status(model)==MOI.OPTIMAL)
        return true,0
    else
        return false,0
    end
end

#PROC is not well commented (but VProc is)
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

    vec_1 = Multiset(B,Int64(denom),2)
    if vec_1 == "time out"
        println(vec_1)
        return false
    end
    vec_2 = Multiset(B,Int64(denom),3)
    if vec_2 == "time out"
        println(vec_2)
        return false
    end
    if(vec_1 !=0 && vec_2 !=0)
        for i=1:length(vec_2)
            push!(vec_1,vec_2[i])
        end
    elseif vec_1 == 0
        vec_1 = vec_2
    elseif vec_1==0 && vec_2==0
        println("No procedure")
    end
    A = Array{Int64,1}(undef,0)
    for i=1:length(vec_1)
        for low_upp = lower_bound_num:upper_bound_num
            count_elem = count(i -> (i==low_upp),vec_1[i])
            append!(A,count_elem)
        end
    end

    shape_1 = Int64(length(A)/(length(B)))
    mat_1 = (reshape(A,(length(B)),shape_1))
    mat_1=unique(mat_1,dims=2)

    V=Int64(ceil(2m/s))
    vec_1 = Multiset(B,Int64((m//s)*denom),V)
    if vec_1 == "time out"
        println(vec_1)
        return false
    end
    vec_2 = Multiset(B,Int64((m//s)*denom),V-1)
    if vec_2 == "time out"
        println(vec_2)
        return false
    end
    AA = Array{Int64,1}(undef,0)
    if(vec_1 !=0 && vec_2 !=0)
        for i=1:length(vec_2)
            push!(vec_1,vec_2[i])
        end
    elseif vec_1 == 0
        vec_1 = vec_2
    elseif vec_1==0 && vec_2==0
        println("No procedure")
    end

    for i=1:length(vec_1)
        for low_upp = lower_bound_num:upper_bound_num
            count_elem = count(i -> (i==low_upp),vec_1[i])
            append!(AA,count_elem)
        end
    end
    shape_2 = Int64(length(AA)/(length(B)))
    mat_2 = (reshape(AA,(length(B)),shape_2))
    mat_2=unique(mat_2,dims=2)
    #display(mat_2)

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

    B=B//denom
    denom=denom*2
    #for k = 1:3
        model=Model(with_optimizer(Cbc.Optimizer, logLevel =0))
        @variable(model, x[i=1:col_1+col_2],Int)
        @constraint(model,con_1,x.>=0)
        @constraint(model,con_2,final_mat*x .==final_mat_2)
        @objective(model, Max, x[1])
        optimize!(model)
        if(termination_status(model)!=MOI.OPTIMAL)
            println("No procedure found")
            return false
        end

        row,col=size(mat_1)
        ogmat_1 = Matrix{Int64}(undef,row,col)
        for i = 1:row
            for j=1:col
                ogmat_1[i,j]=mat_1[i,j]
            end
        end
        row,col=size(mat_2)
        ogmat_2 = Matrix{Int64}(undef,row,col)
        for i = 1:row
            for j=1:col
                ogmat_2[i,j]=mat_2[i,j]
            end
        end
    #    while termination_status(model) == MOI.OPTIMAL

            println("All numbers assumed to have denominator: ",denom)
            for i=1:col_1
                if(value.(x)[i]!=0)
                    print("Cut ",Int64(value.(x)[i])," muffins {  ")

                    for j=1:row_1
                        while(mat_1[j,i]>0)
                            print(Int64(B[j]*denom),"  ")
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
                            print(Int64(B[j]*denom),"  ")
                            mat_2[j,i]=mat_2[j,i]-1
                        end
                    end
                    println("}")
                end
            end
            println()
            println()
            row,col=size(mat_1)
            for i = 1:row
                for j=1:col
                    mat_1[i,j]=ogmat_1[i,j]
                end
            end
            row,col=size(mat_2)
            for i = 1:row
                for j=1:col
                    mat_2[i,j]=ogmat_2[i,j]
                end
            end
    #        y=value.(x)
    #        @constraint(model, x[1]<=y[1]-1)
    #        optimize!(model)
        #    println(termination_status(model))
    #    end
end

#for k=1:20
#    println("*********  k = ",k,"  ********* ")
#    PROC(3k+3,3k+2,((k+1)//(3k+2)))
#end
#VProc(61,19,313//684)

#VProc(5,3,5//12)
#VProc(13,5,13//30)
#VProc(17,15,7//20)
#VProc(28,25,1//3)
#PROC(29,23,49//138)
#PROC(69,32,937//2176)
#PROC(68,53,37//106)
#PROC(107,13,365//754)
#PROC(31,27,103//297)
#PROC(86,17,86//187)
#PROC(21,10,21//50)
#

#PROC(21,10,21//50)
#PROC(31,10,31//70)
#PROC(41,10, 41//90)
#s = time()
#VProc(67,21,41//90)
#println(time()-s)

    PROC(189,83,1143//2656)
