#include("src/combinations.jl")
#include("src/partitions.jl")
include("helper_functions.jl")
using Cbc
using GLPK
using JuMP
using Printf
function VProc(m,s,alphaa, time_limit_solv = 60, time_limit_multi =60, Endpts = 0, proof =0)
#Based on the book algorithm page 42
    #proof = 1
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    x,y=FINDEND(m,s,alphaa,V)
    xbuddy =1-x
    if proof >= 1
        println("\nProcedure for f(",m,", ",s,") = ",alphaa)
        println("***********************************************")
        print("V = ",V,"   ")
        println("s_V = ",sᵥ)
    end
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
    B_1 = collect(lower_bound_num:1:x*denom)
    B_2 = collect(y*denom: upper_bound_num)

    if Endpts !=0
        B=Endpts
        i = 1
        while Endpts[i] <x*denom
            i=i+1
        end
        #println(x*denom)
        #println(i)
        B_1= Endpts[1:i]
        if x!=y
            B_2 = Endpts[i+1:length(Endpts)]
        else
            B_2 = Endpts[i:length(Endpts)]
        end
    end

    #use Mulitset function (in helper_functions) to find
    #the possible mulitsets of piece sizes that sum to denom (1 whole muffin)
    #of sizes 2 (muffins cut into two pieces) or 3 (muffins cut into 3 pieces)
    time_multi = time()
    vec_1 = Multiset(B,Int64(denom),2, time_limit_multi)
    #    vec_1_2 =Multiset(B_2,Int64(denom),2, time_limit_multi)
    #println(B)
    if vec_1 == "time out"
        println("MULTISET 1 TIME OUT:  length(B) = ",length(B),"  sets sum to: ",Int64(denom),"  sets of size: ",2)
        str = @sprintf("MULTISET 2 TIME OUT:  length(B) = %i   sets sum to: %i  sets of size: %i",length(B),Int64(denom),2)
        return -1,str, time()-time_multi, 0
    end
    vec_2 = 0
    if alphaa<1//3
        vec_2_1 = Multiset(B,Int64(denom),3,time_limit_multi)
        if vec_2 == "time out"
            println("MULTISET 2 TIME OUT:  length(B) = ",length(B),"  sets sum to: ",Int64(denom),"  sets of size: ",3)
            str = @sprintf("MULTISET 2 TIME OUT:  length(B) = %i   sets sum to: %i  sets of size: %i",length(B),Int64(denom),3)
            return -1,str, time()-time_multi, 0
        end
    end
    #Mulitset returns an array of arrays of the possible mulitsets
    #or 0 if no mulitsets found, combine the arrays of size 2 and 3
    #into one array of arrays
    if(vec_1 !=0 && vec_2 !=0)
        for i=1:length(vec_2)
            push!(vec_1,vec_2[i])
        end
    elseif vec_1 == 0 && vec_2 !=0
        vec_1 = vec_2
    elseif vec_1==0 && vec_2==0
        return false,0,time()-time_multi, 0
    end
    if alphaa == 1//3
        push!(vec_1, [denom/3; denom/3; denom/3])
    end
    #display(vec_1)
    #count how many times an element of B occurs in a vector
    # and convert to different form (ex B = {5 6 7} V={6 6}-> newV {0 2 0})
    A = Array{Int64,1}(undef,0)
    for i=1:length(vec_1)
        for low_upp in B
          count_elem = count(i -> (i==low_upp),vec_1[i])
          append!(A,count_elem)
        end
    end
    vec_1 = nothing
    vec_2 = nothing
    #reshape inot a matrix
    shape_1 = Int64(length(A)/(length(B)))
    mat_1 = (reshape(A,(length(B)),shape_1))
    mat_1=unique(mat_1,dims=2)
    V=Int64(ceil(2m/s))
    #use Mulitset function (in helper_functions) to find
    #the possible mulitsets of studdent distributions that sum to denom *m/s
#    println("length of B = ",length(B),"  m//s *denom = ",Int64((m//s)*denom))

    vec_1_stu=Vector{Vector{Int64}}(undef,0)
    numSm_B2 = 0
    numLg_B2 = 0
    if length(B_1)==2
        i = V
        #println(V)
        while i>=0
            #println(i*B_1[1]+(V-i)*B_1[2]," == ",m*denom//s)
            if i*B_1[1]+(V-i)*B_1[2] == m*denom//s
                break
            end
            i=i-1
        end
        if i!= 0
            push!(vec_1_stu, [B_1[1]])
            for j=1:i-1
                push!(vec_1_stu[1],B_1[1])
            end
            for j=1:V-i
                push!(vec_1_stu[1],B_1[2])
            end
        else
            push!(vec_1_stu, [B_1[2]])
            for i=1:V-1
                push!(vec_1_stu[1],B_1[2])
            end
        end
        #println(vec_1)
        numSm_B2 = i
        numLg_B2 = V-i
    else
        vec_1_stu = Multiset(B_1,Int64((m//s)*denom),V,time_limit_multi)
        if vec_1_stu == "time out"
            println("MULTISET 3 TIME OUT:  length(B) = ",length(B),"  sets sum to: ",Int64((m//s)*denom),"  sets of size: ",V)
            str = @sprintf("MULTISET 3 TIME OUT:  length(B) = %i   sets sum to: %i  sets of size: %i",length(B),Int64((m//s)*denom),V)
            return -1,str, time()-time_multi, 0
        end

    end

    vec_2_stu=Vector{Vector{Int64}}(undef,0)
    if length(B_2)==2
        i = V-1
        while i>=0
            #println(i*B_2[1]+(V-1-i)*B_2[2]," == ",m*denom//s)
            if i*B_2[1]+(V-1-i)*B_2[2] == m*denom//s
                break
            end
            i=i-1
        end
        if i!= 0
            push!(vec_2_stu, [B_2[1]])
            for j=1:i-1
                push!(vec_2_stu[1],B_2[1])
            end
            for j=1:V-1-i
                push!(vec_2_stu[1],B_2[2])
            end
        else
            push!(vec_2_stu, [B_2[2]])
            for i=1:V-2
                push!(vec_2_stu[1],B_2[2])
            end
        end
        numSm_B2 = i
        numLg_B2 = V-1-i
        #display(vec_2)
        #println(vec_2)
    else
        vec_2_stu = Multiset(B_2,Int64((m//s)*denom),V-1,time_limit_multi)
        if vec_2_stu == "time out"
            println("MULTISET 4 TIME OUT:  length(B) = ",length(B),"  sets sum to: ",Int64((m//s)*denom),"  sets of size: ",V-1)
            str = @sprintf("MULTISET 4 TIME OUT:  length(B) = %i   sets sum to: %i  sets of size: %i",length(B),Int64((m//s)*denom),V-1)
            return -1,str, time()-time_multi, 0
        end
    end

    time_multi = time()- time_multi
    #    println("finished mulitsets")
    #vec_1_length = length(vec_1_stu)
    #vec_2_length = length(vec_2_stu)

    AA_1 = Array{Int64,1}(undef,0)
    for i=1:length(vec_1_stu)
        for low_upp in B_1
          count_elem = count(i -> (i==low_upp),vec_1_stu[i])
          append!(AA_1,count_elem)
        end
    end
    #println(AA)

    #reshape
    shape_2 = Int64(length(AA_1)/(length(B_1)))
    col_mat_1 = shape_2
    #mat_2 = (reshape(AA_1,(length(B_1)),shape_2))
    #mat_2_1=unique(mat_2,dims=2)
    #display(mat_2_1)
    #    if vec_1_length ==1
    #        println(mat_2_1)
    #        println(sᵥ," V students")
    #    end

    if(vec_1_stu !=0 && vec_2_stu !=0)
        for i=1:length(vec_2_stu)
            push!(vec_1_stu,vec_2_stu[i])
        end
    elseif vec_1_stu == 0 && vec_2_stu!=0
        vec_1_stu = vec_2_stu
    elseif vec_1_stu==0 && vec_2_stu==0
        return false, 0, time_multi,0
    end
    #display(vec_1)
    #count how many occurances and convert
    #println(B)
    #println(vec_1)

    AA = Array{Int64,1}(undef,0)
    for i=1:length(vec_1_stu)
        for low_upp in B
          count_elem = count(i -> (i==low_upp),vec_1_stu[i])
          append!(AA,count_elem)
        end
    end
    #println(AA)
    #reshape
    vec_1_stu = nothing
    vec_2_stu = nothing
    shape_2 = Int64(length(AA)/(length(B)))
    mat_2 = (reshape(AA,(length(B)),shape_2))
    mat_2=unique(mat_2,dims=2)
    #println()
    #println(sᵥ)
    #println("*******************")
    #println()
    #display(mat_2)
    #    return
    #this is the system on page 42 with students moved to the
    #left hand side (so it can be set = to 0)

    row_1, col_1 = size(mat_1)
    row_2, col_2 = size(mat_2)


    #add two rows of 0's and 1's
    #(one with 1's under muffin and one with 1's under students)


    k = (m-1)/10
    loc_small = 0
    loc_large = 0
    muff_split = zeros(length(B))
    col_of_mat_1 = 0
    piece_split = zeros(length(B))
    col_of_mat_2 = 0
    if s ==10 && k-floor(k)==0 && k >=2

        k = Int64(k)
        piece_size = (10k+1)//(20k+10) * denom
        for i=1 : length(B)
            if B[i]==piece_size
                loc_small = i
            end
            if B[i] == denom - piece_size
                loc_large = i
            end
        end

        muff_split[loc_small]=1
        muff_split[loc_large] = 1

        for i = 1: col_1
            equal = true
            for j=1: row_1
                if mat_1[j,i] != muff_split[j]
                    equal = false
                end
            end
            if equal == true
                col_of_mat_1 = i
                break
            end
        end

        piece_split[loc_small] = 2k+1

        for i = 1: col_2
            equal = true
            for j=1: row_2
                if mat_2[j,i] != piece_split[j]
                    equal = false
                end
            end
            if equal == true
                col_of_mat_2 = i
                break
            end
        end
    end
    final_mat = [mat_1 (-mat_2)]

    #mat_1 = nothing
    #mat_2 = nothing
    #row_mat_1, col_mat_1 = size(mat_2_1)
    a = (zeros(Int64,col_1+col_2))
    b = (ones(Int64,col_1+col_2))
    c = zeros(Int64,col_1+col_2)
    for i = 1:col_1
        a[i] = 1
        b[i] = 0
    end
    row,col=size(final_mat)
    for i=col_1+col_mat_1+1:col
        b[i]=0
        c[i]=1
    end
    final_mat = [final_mat; a'; b';c']
    #println()
    #println("*******************")
    #display(final_mat)

    #set it equal to final_mat_2 (0's then num muff then num students)
    c = zeros(Int64,row_1)
    final_mat_2 = [c;m;sᵥ; sᵥ₋₁]
    #println("time: ",time()-start_time)
    #solve the systems
    #println("test")

    model=Model(with_optimizer(Cbc.Optimizer,logLevel=0, seconds = time_limit_solv))
    @variable(model, x[i=1:col_1+col_2],Int)
    @constraint(model,con_1,x.>=0)
    @constraint(model,con_2,final_mat*x .==final_mat_2)

    if  length(B_1)==2 && 1-xbuddy<y

        #println("B_1:  ",B_1)
        #println("B:  ",B)
        col_ = 1
        while final_mat[1,col_]!=1
            col_=col_+1
            if col_ > col_1
                break
            end
        end
        _col = 1
        while final_mat[2,_col]!=1
            _col=_col+1
            if _col > col_1
                break
            end
        end
        #display(final_mat)
        @constraint(model, x[col_] == numSm_B2*sᵥ//final_mat[1,col_])
        @constraint(model, x[_col] == numLg_B2*sᵥ//final_mat[2,_col])
    elseif length(B_2)==2 && 1-xbuddy<y

        #println("B_1:  ",B_1)
        #println("B:  ",B)
        col_ = 1
        while final_mat[length(B_1)+1,col_]==0
            col_=col_+1
            if col_ > col_1
                break
            end
        end
        _col = 1
        while final_mat[length(B_1)+2,_col]==0
            _col=_col+1
            if _col > col_1
                break
            end
        end
        #display(final_mat)
        #println(col_,"  ",_col)
        @constraint(model, x[col_ ] == numSm_B2*sᵥ₋₁//final_mat[length(B_1)+1,col_])
        @constraint(model, x[ _col] == numLg_B2*sᵥ₋₁//final_mat[length(B_1)+2,_col])

    end
    if s ==10 && k-floor(k)==0 && k >=2
        #display(B//denom)
        #println()
        #display(mat_1)
        #println("col_of_mat_1 = ",col_of_mat_1)
        #display(mat_2)
        #println("col_of_mat_2 = ",col_of_mat_2)
        #display(final_mat)

        @constraint(model, x[col_of_mat_1] == 4k+2)
        @constraint(model, x[col_of_mat_2+ col_1] == 2)
    end
    #display(final_mat)
    #println("*******************")
    #display(final_mat_2)

    time_solver = time()
    optimize!(model)
    time_solver = time()-time_solver

    #display(value.(x))
    #println("time2: ",time()-start_time)
    #return true if it has a solution, false otherwise
    if (termination_status(model)==MOI.TIME_LIMIT)
        println("SOLVER TIME OUT  dim(A) = ",size(final_mat))
        row,col = size(final_mat)
        str = @sprintf("SOLVER TIME OUT  dim(A) = (%i, %i) ",row,col)
        return -1,str, time_multi, time_solver
    end
    if proof ==1 && termination_status(model)==MOI.OPTIMAL
        B=B//denom
        println( "All numbers assumed to have denominator: ",denom)
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
    elseif proof>=1
        println("NO SOLUTION")
    end
    term_status = termination_status(model)
    if proof ==2
        B=B//denom
        denom=denom
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
        poss_sol = Vector{Vector{Int64}}(undef,0)
        for k = 1:col_1
        #    println("********************************************* k = ",k," *****************************")
            #display(final_mat)
            model=Model(with_optimizer(Cbc.Optimizer, logLevel =0))
            @variable(model, x[i=1:col_1+col_2],Int)
            @constraint(model,con_1,x.>=0)

            @constraint(model,con_2,final_mat*x .==final_mat_2)
            #@constraint(model, x[1] <=6)
            #@constraint(model, x[2] ==7)
            @objective(model, Max, x[k])
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
            while termination_status(model) == MOI.OPTIMAL

                if false
                    println("All numbers assumed to have denominator: ",denom)
                    println(value.(x))
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
                end
                y=value.(x)
                push!(poss_sol,y)
                @constraint(model, x[k]<=y[k] - 1)
                optimize!(model)
            #println(termination_status(model))
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
        end
        a=1
        b=1
        unique!(poss_sol)
        println()
        println(length(poss_sol), " SOLUTIONS")
        println("All numbers assumed to have denominator: ",denom)
        #println(value.(x))

        for j = 1:length(poss_sol)
            for i=1:col_1
                if(poss_sol[j][i]!=0)
                    print("Cut ",Int64(poss_sol[j][i])," muffins {  ")

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
                if(poss_sol[j][i+col_1]!=0)
                    print("Give ",Int64(poss_sol[j][i+col_1])," students {  ")
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
        end && termination_status(model)==MOI.OPTIMAL
    end
    if proof >=1
        println("***********************************************")
    end
    if(termination_status(model)==MOI.OPTIMAL)
        return true,0, time_multi, time_solver
    else
        return false,0, time_multi, time_solver
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
    denom=denom
    poss_sol = Vector{Vector{Int64}}(undef,0)
    for k = 1:col_1
    #    println("********************************************* k = ",k," *****************************")
        #display(final_mat)
        model=Model(with_optimizer(Cbc.Optimizer, logLevel =0))
        @variable(model, x[i=1:col_1+col_2],Int)
        @constraint(model,con_1,x.>=0)

        @constraint(model,con_2,final_mat*x .==final_mat_2)
        #@constraint(model, x[1] <=6)
        #@constraint(model, x[2] ==7)
        @objective(model, Max, x[k])
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
        while termination_status(model) == MOI.OPTIMAL

            if false
                println("All numbers assumed to have denominator: ",denom)
                println(value.(x))
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
            end
            y=value.(x)
            push!(poss_sol,y)
            @constraint(model, x[k]<=y[k] - 1)
            optimize!(model)
        #println(termination_status(model))
        end
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
    a=1
    b=1
    unique!(poss_sol)
    #while a <= length(poss_sol)
    #    b=a+1
    #    while b<=length(poss_sol)
            #println(poss_sol[a][1:col_1])
            #println(poss_sol[b][1:col_1])
            #print(union(poss_sol[a][1:col_1],poss_sol[b][1:col_1]))
            #println()
    #        if length(union(poss_sol[a][1:col_1],poss_sol[b][1:col_1])) == length(unique(poss_sol[a][1:col_1]))
    #            deleteat!(poss_sol, b)
    #        end
    #        b=b+1
        #    println(a," ",b)
    #    end
    #    a=a+1
    #end
    println()
    println(length(poss_sol), " SOLUTIONS")
    println("All numbers assumed to have denominator: ",denom)
    #println(value.(x))

    for j = 1:length(poss_sol)
        for i=1:col_1
            if(poss_sol[j][i]!=0)
                print("Cut ",Int64(poss_sol[j][i])," muffins {  ")

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
            if(poss_sol[j][i+col_1]!=0)
                print("Give ",Int64(poss_sol[j][i+col_1])," students {  ")
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
    end


end


#for k=1:20
#    println("*********  k = ",k,"  ********* ")
#    PROC(3k+3,3k+2,((k+1)//(3k+2)))
#end
#VProc(61,19,313//684)
#PROC(67,21,41//90)
#VProc(5,3,5//12)
#PROC(13,5,13//30)
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

    #PROC(189,83,1143//2656)
#VProc(121,10,121//250)
#start = time()
#for k=1:12
#    println("f(",10*k+1,", 10) ≥ ", (10k+1)//(20k+10),"   ",VProc(10k+1,10,(10k+1)//(20k+10)))
#end
#println(time()-start)

#32 sec with 2 var known #60 sec with no var known
#PROC(21,10,21//50)
#PROC(13,5,13//30)
#PROC(67,21,118//259)
#VProc(81,13,137//286)
#VProc(11,10,7//20)
#VProc(107,13,365//754)
#VProc(11,10,7//20)
#VProc(27,10,9//20)
#VProc(94,25,1087//2350,10,10,0,1
#VProc(83,26,107//234,120,120,0,1)
