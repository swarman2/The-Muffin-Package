#for information on this method see chapter 5 of "The Mathematics of Muffins"
include("helper_functions.jl")
using Cbc
using GLPK
using JuMP
using Printf

#params: m,s,alpha, time limit for the solver, time limit for each mulitset call, endpts (gaps) if prev. found with Mid or GAP
# whether or not to have a proof (0= no proof, 1=procedure, 2=procedure and matrix)
#returns wheter or not a procedure was found (true/false), the time out message (or 0), the time the multiset took, the time the solver took
function VProc(m,s,alphaa, time_limit_solv = 60, time_limit_multi =60, Endpts = 0, proof =0)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    x,y=FINDEND(m,s,alphaa,V)
    xbuddy =1-x
    if proof >= 1
        println("\nProcedure for f(",m,", ",s,") = ",numerator(alphaa),"/",denominator(alphaa))
        println("***********************************************")
        #print("V = ",V,"   ")
        #println("s_V = ",sᵥ)
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
    else

        B = collect(lower_bound_num:1:upper_bound_num)
        B_1 = collect(lower_bound_num:1:x*denom)
        B_2 = collect(y*denom: upper_bound_num)

    end
    l_B = length(B)
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
        return -1,str, time()-time_multi, 0,l_B
    end
    vec_2 = 0
    if alphaa<1//3
        vec_2_1 = Multiset(B,Int64(denom),3,time_limit_multi)
        if vec_2 == "time out"
            println("MULTISET 2 TIME OUT:  length(B) = ",length(B),"  sets sum to: ",Int64(denom),"  sets of size: ",3)
            str = @sprintf("MULTISET 2 TIME OUT:  length(B) = %i   sets sum to: %i  sets of size: %i",length(B),Int64(denom),3)
            return -1,str, time()-time_multi, 0,l_B
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
        return false,0,time()-time_multi, 0,l_B
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

    vec_1_stu=Vector{Vector{Int64}}(undef,0)
    numSm_B2 = 0
    numLg_B2 = 0
    if length(B_1)==2
        i = V
        while i>=0
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
        numSm_B2 = i
        numLg_B2 = V-i
    else
        vec_1_stu = Multiset(B_1,Int64((m//s)*denom),V,time_limit_multi)
        if vec_1_stu == "time out"
            println("MULTISET 3 TIME OUT:  length(B) = ",length(B),"  sets sum to: ",Int64((m//s)*denom),"  sets of size: ",V)
            str = @sprintf("MULTISET 3 TIME OUT:  length(B) = %i   sets sum to: %i  sets of size: %i",length(B),Int64((m//s)*denom),V)
            return -1,str, time()-time_multi, 0,l_B
        end

    end
    vec_2_stu=Vector{Vector{Int64}}(undef,0)
    if length(B_2)==2
        i = V-1
        while i>=0
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
    else
        vec_2_stu = Multiset(B_2,Int64((m//s)*denom),V-1,time_limit_multi)
        if vec_2_stu == "time out"
            println("MULTISET 4 TIME OUT:  length(B) = ",length(B),"  sets sum to: ",Int64((m//s)*denom),"  sets of size: ",V-1)
            str = @sprintf("MULTISET 4 TIME OUT:  length(B) = %i   sets sum to: %i  sets of size: %i",length(B),Int64((m//s)*denom),V-1)
            return -1,str, time()-time_multi, 0,l_B
        end
    end

    time_multi = time()- time_multi
    AA_1 = Array{Int64,1}(undef,0)
    for i=1:length(vec_1_stu)
        for low_upp in B_1
          count_elem = count(i -> (i==low_upp),vec_1_stu[i])
          append!(AA_1,count_elem)
        end
    end

    #reshape
    shape_2 = Int64(length(AA_1)/(length(B_1)))
    col_mat_1 = shape_2

    if(vec_1_stu !=0 && vec_2_stu !=0)
        for i=1:length(vec_2_stu)
            push!(vec_1_stu,vec_2_stu[i])
        end
    elseif vec_1_stu == 0 && vec_2_stu!=0
        vec_1_stu = vec_2_stu
    elseif vec_1_stu==0 && vec_2_stu==0
        return false, 0, time_multi,0,l_B
    end


    AA = Array{Int64,1}(undef,0)
    for i=1:length(vec_1_stu)
        for low_upp in B
          count_elem = count(i -> (i==low_upp),vec_1_stu[i])
          append!(AA,count_elem)
        end
    end
    vec_1_stu = nothing
    vec_2_stu = nothing
    shape_2 = Int64(length(AA)/(length(B)))
    mat_2 = (reshape(AA,(length(B)),shape_2))
    mat_2=unique(mat_2,dims=2)

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

    #set it equal to final_mat_2 (0's then num muff then num students)
    c = zeros(Int64,row_1)
    final_mat_2 = [c;m;sᵥ; sᵥ₋₁]
    #solve the systems
    model=Model(with_optimizer(Cbc.Optimizer,logLevel=0, seconds = time_limit_solv))
    @variable(model, x[i=1:col_1+col_2],Int)
    @constraint(model,con_1,x.>=0)
    @constraint(model,con_2,final_mat*x .==final_mat_2)

    if  length(B_1)==2 && 1-xbuddy<y
        col_ = 1
        while final_mat[1,col_]==0
            col_=col_+1
            if col_ > col_1
                break
            end
        end
        _col = 1
        while final_mat[2,_col]==0
            _col=_col+1
            if _col > col_1
                break
            end
        end

        @constraint(model, x[col_] == numSm_B2*sᵥ//final_mat[1,col_])
        @constraint(model, x[_col] == numLg_B2*sᵥ//final_mat[2,_col])

    elseif length(B_2)==2 && 1-xbuddy<y
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

        @constraint(model, x[col_ ] == numSm_B2*sᵥ₋₁//final_mat[length(B_1)+1,col_])
        @constraint(model, x[ _col] == numLg_B2*sᵥ₋₁//final_mat[length(B_1)+2,_col])

    end
    if s ==10 && k-floor(k)==0 && k >=2
        @constraint(model, x[col_of_mat_1] == 4k+2)
        @constraint(model, x[col_of_mat_2+ col_1] == 2)
    end
    #if proof = 0 solve and don't print anything
    #if proof = 1 solve and print solution
    #if proof = 2 solve and keep looking for solutions
    eps = 10e-5
    if proof == 0
        time_solver = time()
        optimize!(model)
        time_solver = time()-time_solver

        #return true if it has a solution, false otherwise
        if (termination_status(model)==MOI.TIME_LIMIT)
            println("SOLVER TIME OUT  dim(A) = ",size(final_mat))
            row,col = size(final_mat)
            str = @sprintf("SOLVER TIME OUT  dim(A) = (%i, %i) ",row,col)
            return -1,str, time_multi, time_solver,l_B
        end


        term_status = termination_status(model)
    end
    if proof ==1
        time_solver = time()
        optimize!(model)
        time_solver = time()-time_solver
        term_status = termination_status(model)
        #return true if it has a solution, false otherwise
        if (termination_status(model)==MOI.TIME_LIMIT)
            println("SOLVER TIME OUT  dim(A) = ",size(final_mat))
            row,col = size(final_mat)
            str = @sprintf("SOLVER TIME OUT  dim(A) = (%i, %i) ",row,col)
            return -1,str, time_multi, time_solver,l_B
        end
        if termination_status(model) != MOI.OPTIMAL
            println("NO SOLUTION")
        else
            B=B//denom
            y_1 = value.(x)
            println( "All numbers assumed to have denominator: ",denom)
            for i=1:col_1
                if(y_1[i]!=0)
                    if abs(y_1[i] - floor(y_1[i])) <= eps || abs(y_1[i]-ceil(y_1[i]))<=eps# in case of floatint point round error (ex. 3.9999999)
                        y_1[i] = round(y_1[i])
                    end
                    print("Cut ",Int64(y_1[i])," muffins {  ")

                    for j=1:row_1
                        while(mat_1[j,i]>0)
                            print(Int64(B[j]*denom),"  ")
                            mat_1[j,i]=mat_1[j,i]-1
                        end
                    end
                    println("}")
                end
            end
            y = deleteat!(y_1,1:col_1)
            y = filter(z -> z>0, y)
            println()
            for i=1:length(y)
                if(y[i]!=0)
                    if abs(y[i]- floor(y[i])) <= eps || abs(y[i]-ceil(y[i]))<=eps# in case of floatint point round error (ex. 3.9999999)
                        y[i]= round(y[i])
                    end
                    print("Give ",Int64(y[i])," students {  ")
                    for j=1:row_1
                        while(mat_2[j,i]>0)
                            print(Int64(B[j]*denom),"  ")
                            mat_2[j,i]=mat_2[j,i]-1
                        end
                    end
                    println("}")
                end
            end
        end
    end
    if proof ==2
        term_status = 0
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
        time_solver = time()
        for k = 1:col_1

            model=Model(with_optimizer(Cbc.Optimizer, logLevel =0))
            @variable(model, x[i=1:col_1+col_2],Int)
            @constraint(model,con_1,x.>=0)

            @constraint(model,con_2,final_mat*x .==final_mat_2)
            @objective(model, Max, x[k])
            optimize!(model)
            if term_status != MOI.OPTIMAL && termination_status(model)==MOI.OPTIMAL
                term_status = termination_status(model)
            end
            if(termination_status(model)!=MOI.OPTIMAL)
                println("No procedure found")
                return false, 0, 0,0,l_B
            end
            while termination_status(model) == MOI.OPTIMAL
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
        time_solver = time()-time_solver
        a=1
        b=1
        unique!(poss_sol)
        println()
        if length(poss_sol)!= 0
            println(length(poss_sol), " SOLUTIONS")
            println("All numbers assumed to have denominator: ",denom)
        end
        #println(value.(x))

        for j = 1:length(poss_sol)
            for i=1:col_1
                if(poss_sol[j][i]!=0)
                    if abs(poss_sol[i][j] - floor(poss_sol[i][j])) <= eps || abs(poss_sol[i][j]-ceil(poss_sol[i][j]))<=eps
                        poss_sol[i][j]= round(poss_sol[i][j])
                    end
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
                    if abs(poss_sol[i][j] - floor(poss_sol[i][j])) <= eps || abs(poss_sol[i][j]-ceil(poss_sol[i][j]))<=eps
                        poss_sol[i][j]= round(poss_sol[i][j])
                    end
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
            row,col=size(ogmat_1)
            for i = 1:row
                for j=1:col
                    mat_1[i,j]=ogmat_1[i,j]
                end
            end
            row,col=size(ogmat_2)
            for i = 1:row
                for j=1:col
                    mat_2[i,j]=ogmat_2[i,j]
                end
            end
        end
    end
    if proof >=1
        println("***********************************************")
    end
    if(term_status==MOI.OPTIMAL)
        return true,0, time_multi, time_solver,l_B
    else
        return false,0, time_multi, time_solver,l_B
    end
end
#VProc(11,10,7//20,Inf,Inf,0,0)
