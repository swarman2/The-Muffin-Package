include("helper_functions.jl")
using JuMP
using Cbc

function VGAPV3_proof(m,s,alpha)
    println()
    println("PROOF OF f(",m," ",s,") = ",alpha)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    x,y=FINDEND(m,s,alpha,V)
    ybuddy=1-y
    xbuddy=1-x
    println("sᵥ: ",sᵥ,"   sᵥ₋₁: ",sᵥ₋₁)

    if x>=y
        println("INVALID INTERVALS")
        return false
    end
    println(x,"  ***************  ",y)
#    println("x: ",x)
#    println("y: ",y)
#    println("1-x: ",1-x)
#    println("1-y: ",1-y)
    if(V₋₁shares<Vshares)
         #___(_______)________(_____)___|___(_____)____
         #alpha   y-buddy   xbuddy   x  |   y    1-alpha
         #smallShare is an array that holds the possiblities for number of small shares
         num_small_shares = V₋₁shares
        num_large_shares = Vshares - V₋₁shares
        S=sᵥ
        shares=Vshares
        VV=V #which is split
        num_split_shares=Int64(num_large_shares/2)
        endpoints=Array{Rational,2}(undef,0,0)
        endpoints =  [alpha x]
        #endpoints = [endpoints; [xbuddy x]]
        I1 = num_small_shares
        I2 = Int64(num_large_shares/2)
        I3 = I2

#        println("TEST")
    else
        num_small_shares = V₋₁shares-Vshares
        num_large_shares =   Vshares
        S=sᵥ₋₁
        I1=Int64(num_small_shares/2)
        I2=I1
        I3=num_large_shares
        shares=V₋₁shares
        VV=V-1 #which is split
        num_split_shares=Int64(num_small_shares/2)
        endpoints =Array{Rational, 2}(undef,0,0)
        endpoints =[y 1-alpha]
        #endpoints = [endpoints; [xbuddy (1-alpha)]]
        #    sharesInIntervals = Array{Int64}(undef,0)
        #    append!(sharesInIntervals, [Int64(num_small_shares/2) Int64(num_small_shares/2) num_large_shares])

        #    for i = 1:numIntervals
        #        print(endpoints[i,:])
        #        println("  ",sharesInIntervals[i])
        #    end
    end #******************************end else
    #print_Intervals(m,s,alpha)
    endpoints=[alpha x]
    endpoints=[endpoints;y (1-alpha); xbuddy ybuddy]
    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    println("ENDPOINTS BEFORE BUDDY MATCH")
    Display(endpoints)
#    endpoints=buddymatch(endpoints,V,y,m,s)

row,col=size(endpoints)

endpoints=[endpoints; 1//2 1//2 ]


endpoints = buddymatch(endpoints,V,y,m,s)
row,col=size(endpoints)

#for i=1:row
#    if endpoints[i,1]+endpoints[i,2] == (m//s)
#endpoints=[endpoints; [m//2s m//2s]]
#endpoints=sort(collect(Iterators.flatten(endpoints)))
#endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
#endpoints = transpose(endpoints)
#endpoints = buddymatch(endpoints,V,y,m,s)

if V==3
    endpoints = [endpoints; [(m//2s) (m//2s)];[(m//2s) (m//2s)]]
    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    row,col=size(endpoints)

    #special buddy match for closed interval
    for i=1:row -1
        #buddy
        row,col=size(endpoints)
        lower = endpoints[i,2]
        upper = endpoints[i+1,1]

        buddyIn_1 = false
        buddyIn_2=false
        for k=1:row
            row,col=size(endpoints)
            for j=1:2
                if endpoints[k,j] == 1-upper
                    buddyIn_1=true
                end
                if endpoints[k,j]==1-lower
                    buddyIn_2=true
                end
            end
        end
        if buddyIn_1 == false || buddyIn_2 == false
            println("[",1-upper,"  ",1-lower,"] by buddying [",lower,"  ",upper,"]")
            endpoints=[endpoints; [(1-upper) (1-lower)];  [(1-upper) (1-lower)]]
            i=1
            row,col=size(endpoints)
        end
        endpoints=sort(collect(Iterators.flatten(endpoints)))
        endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
        endpoints = transpose(endpoints)
        if lower>=y && V==3 #in the two shares
        #    println("[",m//s-(1-lower),"  ",(m//s)-(1-upper),"] by matching [",1-upper,"  ",1-lower,"]")
            matchIn_1 = false
            matchIn_2=false
            for k=1:row
                row,col=size(endpoints)
                for j=1:2
                    if endpoints[k,j] ==m//s- upper
                        matchIn_1=true
                    end
                    if endpoints[k,j]==m//s-lower
                        matchIn_2=true
                    end
                end
            end
            if matchIn_1 == false || matchIn_2 == false
                println("[",m//s-upper,"  ",m//s-lower,"] by matching [",lower,"  ",upper,"]")
                endpoints=[endpoints; [(m//s-upper) (m//s-lower)];  [(m//s-upper) (m//s-lower)]]
                i=1
                row,col=size(endpoints)
            end
        end
        endpoints=sort(collect(Iterators.flatten(endpoints)))
        endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
        endpoints = transpose(endpoints)
    end
end


    row,col=size(endpoints)

    gap_found=true
    gap_found2=true
    while(gap_found == true)# || gap_found2 == true)
        #buddy match it
        endpoints = buddymatch(endpoints,V,y,m,s)

        row, col= size(endpoints)
        numIntervals = row

        numVIntervals = 0
        for i =1: row
            if endpoints[i,2]<=x
                numVIntervals = numVIntervals+1
            end
        end
        permV=perm(V, numVIntervals)
        permV_ = perm(V-1, row - numVIntervals)


        possV = Vector{Vector{Int64}}()
        possV_ = Vector{Vector{Int64}}()
        println("\n",V," intervals")
        for j=1:numVIntervals
            println("J_",j,": ",endpoints[j,:])
        end
        println("\n",V-1," intervals")
        for j=numVIntervals+1: row
            println("I_",j-numVIntervals,": ",endpoints[j,:])
        end
        #Find the possible distribtions of muffins
        for i=1:length(permV)
            A=permV[i]
            sum_1=0
            sum_2=0

            for j =1:numVIntervals

                sum_1=sum_1+A[j]*endpoints[j,1]
                sum_2=sum_2+A[j]*endpoints[j,2]
            end
            if (sum_1 < m//s) && (sum_2 > m//s)
                append!(possV, permV[i,:])
            end
        end

        #Find the possible distribtions of muffins
        for i=1:length(permV_)
            A=permV_[i]
            sum_1=0
            sum_2=0
            for j =numVIntervals+1:row
                sum_1=sum_1+A[j-numVIntervals]*endpoints[j,1]
                sum_2=sum_2+A[j-numVIntervals]*endpoints[j,2]
            end
            if (sum_1 < m//s) && (sum_2 > m//s)
                append!(possV_, permV_[i,:])
            end
        end
        println("\nPossible ",V," students:")
        println(possV)
        println("\nPossible ",V-1," students:")
        println(possV_)

        mat_V=transpose(hcat(possV...))
    #    endpoints, gap_found = findGaps(mat_V, endpoints,m,s)
        mat_V_=transpose(hcat(possV_...))
        endpoints, gap_found = findGaps(mat_V,mat_V_, endpoints,m,s)
    end #end finding gaps


    println("ENDPOINTS AFTER BUDDY MATCH")
    Display(endpoints)
    row, col= size(endpoints)
    numIntervals = row

    numVIntervals = 0
    for i =1: row
        if endpoints[i,2]<=x
            numVIntervals = numVIntervals+1
        end
    end
    permV=perm(V, numVIntervals)
    permV_ = perm(V-1, row - numVIntervals)


    possV = Vector{Vector{Int64}}()
    possV_ = Vector{Vector{Int64}}()

    #Find the possible distribtions of muffins
    for i=1:length(permV)
        A=permV[i]
        sum_1=0
        sum_2=0
        for j =1:numVIntervals
            sum_1=sum_1+A[j]*endpoints[j,1]
            sum_2=sum_2+A[j]*endpoints[j,2]
        end
        if (sum_1 < m//s) && (sum_2 > m//s)
            append!(possV, permV[i,:])
        end
    end

    #Find the possible distribtions of muffins
    for i=1:length(permV_)
        A=permV_[i]
        sum_1=0
        sum_2=0
        for j =numVIntervals+1:row
            sum_1=sum_1+A[j-numVIntervals]*endpoints[j,1]
            sum_2=sum_2+A[j-numVIntervals]*endpoints[j,2]
        end
        if (sum_1 < m//s) && (sum_2 > m//s)
            append!(possV_, permV_[i,:])
        end
    end
    mat_V=transpose(hcat(possV...))
    mat_V_=transpose(hcat(possV_...))

    symmIntervals = Array{Int64,2}(undef,0,0) # row [i j] if interval i is symm to interval j
    row_endpt,col_endpt = size(endpoints)
    for i = 1:row_endpt
        for j=i+1:row_endpt
            if (endpoints[i,1]+endpoints[j,2])==1 && endpoints[i,2]-endpoints[i,1]==endpoints[j,2]-endpoints[j,1] || (endpoints[i,1]+endpoints[j,2])==m//s && endpoints[i,2]-endpoints[i,1]==endpoints[j,2]-endpoints[j,1]
                if length(symmIntervals) == 0
                    symmIntervals = [i j]
                else
                    symmIntervals = [symmIntervals; i j]
                end
            end
        end
    end
    println("\nSYMMETRIC INTERVALS")
    display(symmIntervals)
#    row_symm, col_symm = size(symmIntervals)
#    #display(symmIntervals)
#    for i=1:row_symm
#        for k=1:2
#            val = symmIntervals[i,k]
#            for j=1:row_symm
#                if j!=i && symmIntervals[j,1]==val
#                    symmIntervals = [symmIntervals; symmIntervals[i,k%2+1] symmIntervals[j,2]]
#                end
#                if j!=i && symmIntervals[j,2]==val
#                    symmIntervals = [symmIntervals; symmIntervals[i,k%2+1] symmIntervals[j,1]]
#                end
#            end
#        end
#    end
#    #display(symmIntervals)
#    sort!(symmIntervals,dims=2)
#    symmIntervals = unique(symmIntervals,dims =2)
#    println("\nSYMMETRIC INTERVALS")
#    display(symmIntervals)
        row_symm, col_symm = size(symmIntervals)

    #display(mat_V)
    #display(mat_V_)
    println("\nPossible distributions of students for the ",V," shares:")
    display(mat_V)
    println("\nPossible distributions of students for the ",V-1," shares:")
    display(mat_V_)
    if length(mat_V)==0
        println("\nNo possible distributions of ",V," shares: alpha ≤ ",alpha)
        return true
    elseif length(mat_V_)==0
        println("\nNo possible distributions of ",V-1," shares: alpha ≤ ",alpha)
        return true
    end
    row_mat,col_mat = size(mat_V)
    row_mat_,col_mat_=size(mat_V_)
    row_mat_new=0
    while(row_mat_new < row_mat + row_mat_)
#        println("test")
        Z= zeros(Int64,1,col_mat)
        mat_V = [mat_V; Z]
        row_mat_new,col_mat = size(mat_V)
    end
    row_mat_new=0
    while(row_mat_new < row_mat + row_mat_)
        Z = zeros(Int64,1,col_mat_)
        mat_V_ = [Z; mat_V_]
        row_mat_new,col_mat_ = size(mat_V_)
    end
    poss_Dist = [mat_V mat_V_]

    println("\nConcationate them in order to make one system of equations")
    display(poss_Dist)

    row_Dist, col_Dist = size(poss_Dist)
    row_endpt, col_endpt = size(endpoints)
    while col_Dist < row_endpt
        poss_Dist = [poss_Dist zeros(row_Dist)]
        row_Dist, col_Dist = size(poss_Dist)
    end

    A=Array{Int64}(undef,0) # A will become the matrix in the eqation Ax=b
    for i=1:row_symm
        append!(A, (poss_Dist[:,symmIntervals[i,1]]-poss_Dist[:,symmIntervals[i,2]]))
    end
    A=reshape(A,Int64(length(A)/row_symm),row_symm)
    A=transpose(A)
    #A is now a matrix with rows correspoding to the symmetric identities
    #ex if the intervals are I1, I2, I3, I4, I5 (|I1|=|I4|,|I2|=|I3|)
    #A = [I1-I4 ; I2-I3]


    println("\nIf intervals i and j are symmetric, subtract col i from col j and add as a row to A")
    A=unique(A, dims = 1)
    display(A)
    row,col=size(A)
    #add zeros for each row of I₁ - I₂
    b=zeros(Int64,row)
    #sum up cols in mat_V
    #sum up cols in mat_V
    if length(mat_V)!=0
        Sum_1 = mat_V[:,1]
        r,c=size(mat_V)
        for i = 2:c
            Sum_1 = Sum_1 + mat_V[:,i]
        end
        #transpose to get right dim
        Sum_1=Sum_1'
        b=[b;Vshares]
        #append as row to A
        A=[A; Sum_1]
    end
    if length(mat_V_)!=0
        Sum_2 = mat_V_[:,1]
        r,c=size(mat_V_)
        for i = 2:c
            Sum_2 = Sum_2 + mat_V_[:,i]
        end
        #transpose to get right dim
        Sum_2=Sum_2'
        b=[b;V₋₁shares]
        #append as row to A
        A=[A;Sum_2]
    end


    println("\nAdd a row of ",V,"'s and a row of ",V-1,"'s \n(these when mulitplied by number of each type of distribution of students will add to number of ",V," students and number of ",V-1," students)")
    #append as rows to A
    display(A)
    #Append a row of ones to A (num stud per distribtion adds to num VV students)

    Ones=(ones(Int64,col))'
#    if(V!=3)
#        A=vcat(A, Ones)
#    end
    split_Intervals_1 = Array{Int64}(undef,0)
    split_Intervals_2 = Array{Int64}(undef,0)
    if VV==V
        row,col=size(endpoints)
        for i = 1:row
            if endpoints[i,1]>=xbuddy && endpoints[i,2]<=1//2
                append!(split_Intervals_1, i)
            end
            if endpoints[i,1]>=1//2 && endpoints[i,2]<=x
                append!(split_Intervals_2,i)
            end
        end
    else
        row,col=size(endpoints)
        for i = 1:row
            if endpoints[i,1]>=y && endpoints[i,2]<=1//2
                append!(split_Intervals_1, i)
            end
            if endpoints[i,1]>=1//2 && endpoints[i,2]<=ybuddy
                append!(split_Intervals_2,i)
            end
        end
    end
    #display(split_Intervals_1)
    #display(split_Intervals_2)
    symm_sum_1 = Array{Int64}(undef,0)
    symm_sum_2 = Array{Int64}(undef,0)
    for i in split_Intervals_1
        if length(symm_sum_1) == 0
            symm_sum_1 = poss_Dist[:,i]
        else
            symm_sum_1 =symm_sum_1+ poss_Dist[:,i]
        end

    end
    for i in split_Intervals_2
        if length(symm_sum_2) == 0
            symm_sum_2 = poss_Dist[:,i]
        else
            symm_sum_2 =symm_sum_2 + poss_Dist[:,i]
        end
    end
    symm_sum_1=symm_sum_1'
    symm_sum_2 = symm_sum_2'
    #display(symm_sum_1)
    #display(symm_sum_2)
    if length(symm_sum_1)!=0
        A=[A;symm_sum_1]
        b=[b;num_split_shares]
    end
    if length(symm_sum_2) !=0
        A=[A;symm_sum_2]
        b=[b; num_split_shares]
    end
    println("\nAdd a rows for the sums of the ",VV," shares symmetric around 1//2")
    display(A)
    println("Now solve the system for the number of each type of distribtuion (= ",b,")")

    row,col=size(A)
    #display(A)
    model=Model(with_optimizer(Cbc.Optimizer, logLevel=0))
    @variable(model, x[i=1:col],Int)

    @constraint(model,con,A*x .==b)
    @constraint(model,con_1,x.>=0)

    optimize!(model)
    if(termination_status(model)==MOI.OPTIMAL)
        display(value.(x))
        println("\nThere is a solution on the naturals, alpha > ",alpha)
        return false
    else
        println("\nThere is no solution on the naturals, alpha ≤ ",alpha)
        return true
    end
end

function VGAPV3(m,s,alpha)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    x,y=FINDEND(m,s,alpha,V)
    ybuddy=1-y
    xbuddy=1-x

    if x>=y

        return false
    end

#    println("x: ",x)
#    println("y: ",y)
#    println("1-x: ",1-x)
#    println("1-y: ",1-y)
    if(V₋₁shares<Vshares)
         #___(_______)________(_____)___|___(_____)____
         #alpha   y-buddy   xbuddy   x  |   y    1-alpha
         #smallShare is an array that holds the possiblities for number of small shares
         num_small_shares = V₋₁shares
        num_large_shares = Vshares - V₋₁shares
        S=sᵥ
        shares=Vshares
        VV=V #which is split
        num_split_shares=Int64(num_large_shares/2)
        endpoints=Array{Rational,2}(undef,0,0)
        endpoints =  [alpha x]
        #endpoints = [endpoints; [xbuddy x]]
        I1 = num_small_shares
        I2 = Int64(num_large_shares/2)
        I3 = I2

#        println("TEST")
    else
        num_small_shares = V₋₁shares-Vshares
        num_large_shares =   Vshares
        S=sᵥ₋₁
        I1=Int64(num_small_shares/2)
        I2=I1
        I3=num_large_shares
        shares=V₋₁shares
        VV=V-1 #which is split
        num_split_shares=Int64(num_small_shares/2)
        endpoints =Array{Rational, 2}(undef,0,0)
        endpoints =[y 1-alpha]
        #endpoints = [endpoints; [xbuddy (1-alpha)]]
        #    sharesInIntervals = Array{Int64}(undef,0)
        #    append!(sharesInIntervals, [Int64(num_small_shares/2) Int64(num_small_shares/2) num_large_shares])

        #    for i = 1:numIntervals
        #        print(endpoints[i,:])
        #        println("  ",sharesInIntervals[i])
        #    end
    end #******************************end else
    #print_Intervals(m,s,alpha)
    endpoints=[alpha x]
    endpoints=[endpoints;y (1-alpha); xbuddy ybuddy]
    #endpoints=sort(collect(Iterators.flatten(endpoints)))
    #endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    #endpoints = transpose(endpoints)

#    endpoints=buddymatch(endpoints,V,y,m,s)

row,col=size(endpoints)

        endpoints=[endpoints; 1//2 1//2 ]


endpoints = buddymatch(endpoints,V,y,m,s)
row,col=size(endpoints)

#for i=1:row
#    if endpoints[i,1]+endpoints[i,2] == (m//s)
#endpoints=[endpoints; [m//2s m//2s]]
#endpoints=sort(collect(Iterators.flatten(endpoints)))
#endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
#endpoints = transpose(endpoints)
#endpoints = buddymatch(endpoints,V,y,m,s)

if V==3
    endpoints = [endpoints; [(m//2s) (m//2s)];[(m//2s) (m//2s)]]
    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    row,col=size(endpoints)

    #special buddy match for closed interval
    for i=1:row -1
        #buddy
        row,col=size(endpoints)
        lower = endpoints[i,2]
        upper = endpoints[i+1,1]

        buddyIn_1 = false
        buddyIn_2=false
        for k=1:row
            row,col=size(endpoints)
            for j=1:2
                if endpoints[k,j] == 1-upper
                    buddyIn_1=true
                end
                if endpoints[k,j]==1-lower
                    buddyIn_2=true
                end
            end
        end
        if buddyIn_1 == false || buddyIn_2 == false

            endpoints=[endpoints; [(1-upper) (1-lower)];  [(1-upper) (1-lower)]]
            i=1
            row,col=size(endpoints)
        end
        endpoints=sort(collect(Iterators.flatten(endpoints)))
        endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
        endpoints = transpose(endpoints)
        if lower>=y && V==3 #in the two shares
        #    println("[",m//s-(1-lower),"  ",(m//s)-(1-upper),"] by matching [",1-upper,"  ",1-lower,"]")
            matchIn_1 = false
            matchIn_2=false
            for k=1:row
                row,col=size(endpoints)
                for j=1:2
                    if endpoints[k,j] ==m//s- upper
                        matchIn_1=true
                    end
                    if endpoints[k,j]==m//s-lower
                        matchIn_2=true
                    end
                end
            end
            if matchIn_1 == false || matchIn_2 == false
                endpoints=[endpoints; [(m//s-upper) (m//s-lower)];  [(m//s-upper) (m//s-lower)]]
                i=1
                row,col=size(endpoints)
            end
        end
        endpoints=sort(collect(Iterators.flatten(endpoints)))
        endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
        endpoints = transpose(endpoints)
    end
end


    row,col=size(endpoints)

    gap_found=true
    gap_found2=true
    while(gap_found == true)# || gap_found2 == true)
        #buddy match it
        endpoints = buddymatch(endpoints,V,y,m,s)

        row, col= size(endpoints)
        numIntervals = row

        numVIntervals = 0
        for i =1: row
            if endpoints[i,2]<=x
                numVIntervals = numVIntervals+1
            end
        end
        permV=perm(V, numVIntervals)
        permV_ = perm(V-1, row - numVIntervals)


        possV = Vector{Vector{Int64}}()
        possV_ = Vector{Vector{Int64}}()

        #Find the possible distribtions of muffins
        for i=1:length(permV)
            A=permV[i]
            sum_1=0
            sum_2=0

            for j =1:numVIntervals

                sum_1=sum_1+A[j]*endpoints[j,1]
                sum_2=sum_2+A[j]*endpoints[j,2]
            end
            if (sum_1 < m//s) && (sum_2 > m//s)
                append!(possV, permV[i,:])
            end
        end

        #Find the possible distribtions of muffins
        for i=1:length(permV_)
            A=permV_[i]
            sum_1=0
            sum_2=0
            for j =numVIntervals+1:row
                sum_1=sum_1+A[j-numVIntervals]*endpoints[j,1]
                sum_2=sum_2+A[j-numVIntervals]*endpoints[j,2]
            end
            if (sum_1 < m//s) && (sum_2 > m//s)
                append!(possV_, permV_[i,:])
            end
        end

        mat_V=transpose(hcat(possV...))
    #    endpoints, gap_found = findGaps(mat_V, endpoints,m,s)
        mat_V_=transpose(hcat(possV_...))
        endpoints, gap_found = findGaps(mat_V,mat_V_, endpoints,m,s)
    end #end finding gaps


    row, col= size(endpoints)
    numIntervals = row

    numVIntervals = 0
    for i =1: row
        if endpoints[i,2]<=x
            numVIntervals = numVIntervals+1
        end
    end
    permV=perm(V, numVIntervals)
    permV_ = perm(V-1, row - numVIntervals)


    possV = Vector{Vector{Int64}}()
    possV_ = Vector{Vector{Int64}}()

    #Find the possible distribtions of muffins
    for i=1:length(permV)
        A=permV[i]
        sum_1=0
        sum_2=0
        for j =1:numVIntervals
            sum_1=sum_1+A[j]*endpoints[j,1]
            sum_2=sum_2+A[j]*endpoints[j,2]
        end
        if (sum_1 < m//s) && (sum_2 > m//s)
            append!(possV, permV[i,:])
        end
    end

    #Find the possible distribtions of muffins
    for i=1:length(permV_)
        A=permV_[i]
        sum_1=0
        sum_2=0
        for j =numVIntervals+1:row
            sum_1=sum_1+A[j-numVIntervals]*endpoints[j,1]
            sum_2=sum_2+A[j-numVIntervals]*endpoints[j,2]
        end
        if (sum_1 < m//s) && (sum_2 > m//s)
            append!(possV_, permV_[i,:])
        end
    end
    mat_V=transpose(hcat(possV...))
    mat_V_=transpose(hcat(possV_...))

    symmIntervals = Array{Int64,2}(undef,0,0) # row [i j] if interval i is symm to interval j
    row_endpt,col_endpt = size(endpoints)
    for i = 1:row_endpt
        for j=i+1:row_endpt
            if (endpoints[i,1]+endpoints[j,2])==1 && endpoints[i,2]-endpoints[i,1]==endpoints[j,2]-endpoints[j,1] || (endpoints[i,1]+endpoints[j,2])==m//s && endpoints[i,2]-endpoints[i,1]==endpoints[j,2]-endpoints[j,1]
                if length(symmIntervals) == 0
                    symmIntervals = [i j]
                else
                    symmIntervals = [symmIntervals; i j]
                end
            end
        end
    end

#    row_symm, col_symm = size(symmIntervals)
#    #display(symmIntervals)
#    for i=1:row_symm
#        for k=1:2
#            val = symmIntervals[i,k]
#            for j=1:row_symm
#                if j!=i && symmIntervals[j,1]==val
#                    symmIntervals = [symmIntervals; symmIntervals[i,k%2+1] symmIntervals[j,2]]
#                end
#                if j!=i && symmIntervals[j,2]==val
#                    symmIntervals = [symmIntervals; symmIntervals[i,k%2+1] symmIntervals[j,1]]
#                end
#            end
#        end
#    end
#    #display(symmIntervals)
#    sort!(symmIntervals,dims=2)
#    symmIntervals = unique(symmIntervals,dims =2)
#    println("\nSYMMETRIC INTERVALS")
#    display(symmIntervals)
        row_symm, col_symm = size(symmIntervals)

    #display(mat_V)
    #display(mat_V_)

    if length(mat_V)==0
        return true
    elseif length(mat_V_)==0
        return true
    end
    row_mat,col_mat = size(mat_V)
    row_mat_,col_mat_=size(mat_V_)
    row_mat_new=0
    while(row_mat_new < row_mat + row_mat_)
#        println("test")
        Z= zeros(Int64,1,col_mat)
        mat_V = [mat_V; Z]
        row_mat_new,col_mat = size(mat_V)
    end
    row_mat_new=0
    while(row_mat_new < row_mat + row_mat_)
        Z = zeros(Int64,1,col_mat_)
        mat_V_ = [Z; mat_V_]
        row_mat_new,col_mat_ = size(mat_V_)
    end
    poss_Dist = [mat_V mat_V_]

    row_Dist, col_Dist = size(poss_Dist)
    row_endpt, col_endpt = size(endpoints)
    while col_Dist < row_endpt
        poss_Dist = [poss_Dist zeros(row_Dist)]
        row_Dist, col_Dist = size(poss_Dist)
    end

    A=Array{Int64}(undef,0) # A will become the matrix in the eqation Ax=b
    for i=1:row_symm
        append!(A, (poss_Dist[:,symmIntervals[i,1]]-poss_Dist[:,symmIntervals[i,2]]))
    end
    A=reshape(A,Int64(length(A)/row_symm),row_symm)
    A=transpose(A)
    #A is now a matrix with rows correspoding to the symmetric identities
    #ex if the intervals are I1, I2, I3, I4, I5 (|I1|=|I4|,|I2|=|I3|)
    #A = [I1-I4 ; I2-I3]



    A=unique(A, dims = 1)
    row,col=size(A)
    #add zeros for each row of I₁ - I₂
    b=zeros(Int64,row)
    #sum up cols in mat_V
    #sum up cols in mat_V
    if length(mat_V)!=0
        Sum_1 = mat_V[:,1]
        r,c=size(mat_V)
        for i = 2:c
            Sum_1 = Sum_1 + mat_V[:,i]
        end
        #transpose to get right dim
        Sum_1=Sum_1'
        b=[b;Vshares]
        #append as row to A
        A=[A; Sum_1]
    end
    if length(mat_V_)!=0
        Sum_2 = mat_V_[:,1]
        r,c=size(mat_V_)
        for i = 2:c
            Sum_2 = Sum_2 + mat_V_[:,i]
        end
        #transpose to get right dim
        Sum_2=Sum_2'
        b=[b;V₋₁shares]
        #append as row to A
        A=[A;Sum_2]
    end



    #Append a row of ones to A (num stud per distribtion adds to num VV students)

    Ones=(ones(Int64,col))'
#    if(V!=3)
#        A=vcat(A, Ones)
#    end
    split_Intervals_1 = Array{Int64}(undef,0)
    split_Intervals_2 = Array{Int64}(undef,0)
    if VV==V
        row,col=size(endpoints)
        for i = 1:row
            if endpoints[i,1]>=xbuddy && endpoints[i,2]<=1//2
                append!(split_Intervals_1, i)
            end
            if endpoints[i,1]>=1//2 && endpoints[i,2]<=x
                append!(split_Intervals_2,i)
            end
        end
    else
        row,col=size(endpoints)
        for i = 1:row
            if endpoints[i,1]>=y && endpoints[i,2]<=1//2
                append!(split_Intervals_1, i)
            end
            if endpoints[i,1]>=1//2 && endpoints[i,2]<=ybuddy
                append!(split_Intervals_2,i)
            end
        end
    end
    #display(split_Intervals_1)
    #display(split_Intervals_2)
    symm_sum_1 = Array{Int64}(undef,0)
    symm_sum_2 = Array{Int64}(undef,0)
    for i in split_Intervals_1
        if length(symm_sum_1) == 0
            symm_sum_1 = poss_Dist[:,i]
        else
            symm_sum_1 =symm_sum_1+ poss_Dist[:,i]
        end

    end
    for i in split_Intervals_2
        if length(symm_sum_2) == 0
            symm_sum_2 = poss_Dist[:,i]
        else
            symm_sum_2 =symm_sum_2 + poss_Dist[:,i]
        end
    end
    symm_sum_1=symm_sum_1'
    symm_sum_2 = symm_sum_2'
    #display(symm_sum_1)
    #display(symm_sum_2)
    if length(symm_sum_1)!=0
        A=[A;symm_sum_1]
        b=[b;num_split_shares]
    end
    if length(symm_sum_2) !=0
        A=[A;symm_sum_2]
        b=[b; num_split_shares]
    end

    row,col=size(A)
    #display(A)
    model=Model(with_optimizer(Cbc.Optimizer, logLevel=0))
    @variable(model, x[i=1:col],Int)

    @constraint(model,con,A*x .==b)
    @constraint(model,con_1,x.>=0)

    optimize!(model)
    if(termination_status(model)==MOI.OPTIMAL)
        return false
    else
        return true
    end

end

function findGaps(possDistV,possDistV_,endpoints,m,s)
    #find gaps
    """
    #this loop sorts the possible shares by which intervals are used
    for j=1:numIntervals
        for i in possInd
            if(X[i][j]!=0)
                println(X[i])
            end
        end
        println()
    end
    """
#    println("TEST")
    lowerbound = 0
    upperbound = 0
    _gap=false
    numVDistributions,numVIntervals=size(possDistV)
    numV_Distributions, numV_Intervals=size(possDistV_)
#    println("\nPOSSDISTV")
#    display(possDistV)
#    println("\nPOSSDISTV_")
#    display(possDistV_)
    #display(possDist)
    newendpoints=endpoints
    dist1=Vector{Vector{Any}}(undef,2)
    temp_endpt = Vector{Rational}(undef,0)
    lowEnd = Array{Rational}(undef,0)
    highEnd = Array{Rational}(undef,0)
    for j=1:numVIntervals
        #println("\nChecking interval: ",j)
        lowEnd = Array{Rational}(undef,0)
        highEnd = Array{Rational}(undef,0)
        for i=1:numVDistributions
            if(possDistV[i,j]!=0)
        #        println("\nUsing student: ",possDistV[i,:])
                #all types of students that use interval j
                #X[i] is the student type
                #sharesInIntervals[j] is the number of shares in that interval
                #endpoints[j,:] is the endpoints of that interval

                #lower bound
                sum=0
            #    println("X[i][j]= ",X[i][j])
            #    println("X[i] = ",X[i])
                for k=1:length(possDistV[i,:])
                    if(k!=j)
                    #    println("endpoints[",k,",2]*possDistV[",i,"][",k,"] = ",endpoints[k,2]*possDistV[i,k])
                        sum=sum+endpoints[k,2]*possDistV[i,k]
                    else
                        sum=sum+endpoints[k,2]*(possDistV[i,k]-1)
                    end
                end
                lowerbound_temp=m//s-sum
            #    println("lowerbound_temp = ", lowerbound_temp)

                #upper bound
                sum=0
                for k=1:length(possDistV[i,:])
                    if(k!=j)
                #        println("endpoints[",k,",1]*possDistV[",i,"][",k,"] = ",endpoints[k,1]*possDistV[i,k])
                        sum=sum+endpoints[k,1]*possDistV[i,k]
                    else
                        sum=sum+endpoints[k,1]*(possDistV[i,k]-1)
                    end
                end
                upperbound_temp=m//s-sum
                #push!(temp_endpt,[lowerbound_temp upperbound_temp])
                if lowerbound_temp<endpoints[j,1]
                    lowerbound_temp=endpoints[j,1]
                end
                if upperbound_temp>endpoints[j,2]
                    upperbound_temp=endpoints[j,2]
                end
                append!(lowEnd,lowerbound_temp)
                append!(highEnd,upperbound_temp)
                #append!(temp_endpt,upperbound_temp)
            #    println("upperbound_temp = ",upperbound_temp)

            end #end if(X[i][j]!=0)
        end  #end  for i in possInd
        temp_endpt = Vector{Rational}(undef,0)
     #  println()
      # println(lowEnd*230)
        #println(highEnd*230)
    #    println()

    while unique(lowEnd)!=lowEnd
        for i =1:length(lowEnd)
            for j=i+1:length(lowEnd)
                if lowEnd[i]==lowEnd[j]
                    if highEnd[i]>=highEnd[j]
                        deleteat!(lowEnd,j)
                        deleteat!(highEnd,j)
                    else
                        deleteat!(lowEnd,i)
                        deleteat!(highEnd,i)
                    end
                    break
                end
            end
        end

    #    println()
    #    println(lowEnd*230)
    #    println(highEnd*230)
    #    println()

    end

        while (length(lowEnd)>0)
            lowIndex=1
            #println(lowEnd*230)
            for i = 1:length(lowEnd)
                if lowEnd[i]<lowEnd[lowIndex]
                    lowIndex=i
                end
            end
    #        println("lowIndex = ",lowIndex)
            push!(temp_endpt,lowEnd[lowIndex])
            push!(temp_endpt,highEnd[lowIndex])
            deleteat!(lowEnd,lowIndex)
            deleteat!(highEnd,lowIndex)
        end
    #    println("**********************")
    #    for i=1:Int64(length(temp_endpt)/2)
    #        print(temp_endpt[2i-1]*230,"  ")
        #end
        #println()
        #println("**********************")
        merge=true
    #    println("LOOKING FOR GAPS IN INTERVAL ",j,":  [",endpoints[j,1],",",endpoints[j,2],"]")
    #    println(temp_endpt*230)
    #println("******************************")
        while merge==true
            merge = false
    #        println(temp_endpt*230)
        #    println(length(temp_endpt))
            for i=1:Int64(length(temp_endpt)/2-1)
                if temp_endpt[2*i] >= temp_endpt[2*i+1]

                    merge=true
                    if (temp_endpt[2*i]<temp_endpt[2*(i+1)])
                        deleteat!(temp_endpt,2i+1)
                        deleteat!(temp_endpt,2i)
                    else
                        deleteat!(temp_endpt,2(i+1))
                        deleteat!(temp_endpt,2i+1)
                    end
                    break
                end
            end
        end
    #    println("######################  ", temp_endpt*230)
    #    println("**********")
    #    println(endpoints[j,1]*230,"   ",endpoints[j,2]*230)

        endpt=filter(x-> x > endpoints[j,1],temp_endpt)
        endpt=filter(x-> x < endpoints[j,2],endpt)

    #    println(endpt)
        #println("-----------------------",endpt*246)
        if length(endpt)%2==0
            for i=1:Int64(length(endpt)/2)
                upperbound=endpt[i]
                lowerbound=endpt[2i]
                if lowerbound>upperbound
                        #display(temp_endpt*246)
            #           println("ADDED: [",upperbound*246,", ",lowerbound*246,"]")
                #    if length(temp_endpt) == 2
                #        temp_endpt=sort(temp_endpt)
                    #    lowerbound=temp_endpt[i][2]
                    #    upperbound = temp_endpt[i][1]
                    #check to see if the lowerbound and upper bound are valid for a gap

                #    if(lowerbound>upperbound && lowerbound>endpoints[j,1] && upperbound < endpoints[j,2])
                        _gap=true
                        newendpoints=[newendpoints; upperbound lowerbound]#; (1-lowerbound) (1-upperbound)]
                    #    println("\n[",lowerbound,"  ",upperbound,"] = ∅")
                    #    println("\n[",lowerbound," ",upperbound,"] empty : ",dist1)
                        newendpoints=sort(collect(Iterators.flatten(newendpoints)))
                        newendpoints=reshape(newendpoints,(2,Int64(length(newendpoints)/2)))
                        newendpoints = transpose(newendpoints)
                        #numVIntervals = numVIntervals+1

                end #end if(lowerbound>upperbound)
            end#end for i = 1:Int64(length(temp_endpt)/2)
        elseif length(endpt)==1
            if temp_endpt[1]<endpoints[j,1] #then change upper bound
                endpoints[j,2]=temp_endpt[2]
            else
                endpoints[j,1]=temp_endpt[1]
            end
        else
            println()
            println("ERROR")
            println(temp_endpt*230)
            println(endpt*230)
            println()
        end

    end #end for j=1:numInterval

    temp_endpt = Array{Rational}(undef,0)
    lowEnd = Array{Rational}(undef,0)
    highEnd = Array{Rational}(undef,0)
    for j=numVIntervals+1:numVIntervals + numV_Intervals
        lowEnd = Array{Rational}(undef,0)
        highEnd = Array{Rational}(undef,0)
        for i=1:numV_Distributions
            if(possDistV_[i,j-numVIntervals]!=0)
                #all types of students that use interval j
                #X[i] is the student type
                #sharesInIntervals[j] is the number of shares in that interval
                #endpoints[j,:] is the endpoints of that interval

                #lower bound
                sum=0
            #    println("X[i][j]= ",X[i][j])
            #    println("X[i] = ",X[i])
                for k=numVIntervals+1:numVIntervals+length(possDistV_[i,:])
                #    println("endpoints [",k,",",2,"] = ",endpoints[k,2])
                    if(k!=j)
            #            println("endpoints[",k,",2]*X[",i,"][",k,"] = ",endpoints[k,2]*X[i][k])
                        sum=sum+endpoints[k,2]*possDistV_[i,k-numVIntervals]
                    else
                        sum=sum+endpoints[k,2]*(possDistV_[i,k-numVIntervals]-1)
                    end
                end
                lowerbound_temp=m//s-sum
            #    println("lowerbound_temp = ", lowerbound_temp)
                #upper bound
                sum=0
                for k=numVIntervals+1:numVIntervals+length(possDistV_[i,:])
                #    println("endpoints [",k,",",1,"] = ",endpoints[k,1])
                    if(k!=j)
                        sum=sum+endpoints[k,1]*possDistV_[i,k-numVIntervals]
                    else
                        sum=sum+endpoints[k,1]*(possDistV_[i,k-numVIntervals]-1)
                    end

                end
                upperbound_temp=m//s-sum
                if lowerbound_temp<endpoints[j,1]
                    lowerbound_temp=endpoints[j,1]
                end
                if upperbound_temp>endpoints[j,2]
                    upperbound_temp=endpoints[j,2]
                end
                append!(lowEnd,lowerbound_temp)
                append!(highEnd,upperbound_temp)
            end #end if(X[i][j]!=0)
        end  #end  for i in possInd
        while unique(lowEnd)!=lowEnd
            for i =1:length(lowEnd)
                for j=i+1:length(lowEnd)
                    if lowEnd[i]==lowEnd[j]
                        if highEnd[i]>=highEnd[j]
                            deleteat!(lowEnd,j)
                            deleteat!(highEnd,j)
                        else
                            deleteat!(lowEnd,i)
                            deleteat!(highEnd,i)
                        end
                        break
                    end
                end
            end
        end
        temp_endpt = Vector{Rational}(undef,0)
        while (length(lowEnd)>0)
            lowIndex=1
            for i = 1:length(lowEnd)
                if lowEnd[i]<lowEnd[lowIndex]
                    lowIndex=i
                end
            end
            push!(temp_endpt,lowEnd[lowIndex])
            push!(temp_endpt,highEnd[lowIndex])
            deleteat!(lowEnd,lowIndex)
            deleteat!(highEnd,lowIndex)
        end
    #    println("**********************")
    #    for i=1:Int64(length(temp_endpt)/2)
    #        print(temp_endpt[2i],"  ")
    #    end
    #    println()
    #    println("**********************")
    #    println("***********************")
        merge=true
        while merge==true
            merge = false
        #    println(temp_endpt)
        #    println(length(temp_endpt))
            for i=1:Int64(length(temp_endpt)/2-1)
                if temp_endpt[2*i] >= temp_endpt[2*i+1]
                    merge = true
                    if (temp_endpt[2*i]<temp_endpt[2*(i+1)])
                        deleteat!(temp_endpt,2i+1)
                        deleteat!(temp_endpt,2i)
                    else
                        deleteat!(temp_endpt,2(i+1))
                        deleteat!(temp_endpt,2i+1)
                    end
                    break
                end
            end
        end

        #println(temp_endpt*246)
        #println("**********")
        #println(endpoints[j,1]*246,"   ",endpoints[j,2]*246)
        endpt=filter(x-> x > endpoints[j,1],temp_endpt)
        endpt=filter(x-> x < endpoints[j,2],endpt)
        #println("-----------------------",endpt*246)
        if length(endpt)%2==0
            for i=1:Int64(length(endpt)/2)
                upperbound=endpt[i]
                lowerbound=endpt[2i]
                if lowerbound>upperbound
                        #display(temp_endpt*246)
            #           println("ADDED: [",upperbound*246,", ",lowerbound*246,"]")
                #    if length(temp_endpt) == 2
                #        temp_endpt=sort(temp_endpt)
                    #    lowerbound=temp_endpt[i][2]
                    #    upperbound = temp_endpt[i][1]
                    #check to see if the lowerbound and upper bound are valid for a gap

                #    if(lowerbound>upperbound && lowerbound>endpoints[j,1] && upperbound < endpoints[j,2])
                        _gap=true
                        newendpoints=[newendpoints; upperbound lowerbound]#; (1-lowerbound) (1-upperbound)]
                    #    println("\n[",lowerbound,"  ",upperbound,"] = ∅")
                    #    println("\n[",lowerbound," ",upperbound,"] empty : ",dist1)
                        newendpoints=sort(collect(Iterators.flatten(newendpoints)))
                        newendpoints=reshape(newendpoints,(2,Int64(length(newendpoints)/2)))
                        newendpoints = transpose(newendpoints)
                        #numVIntervals = numVIntervals+1

                end #end if(lowerbound>upperbound)
            end#end for i = 1:Int64(length(temp_endpt)/2)
        elseif length(endpt)==1
            if temp_endpt[1]<endpoints[j,1] #then change upper bound
                endpoints[j,2]=temp_endpt[2]
            else
                endpoints[j,1]=temp_endpt[1]
            end
        else
            println()
            println("ERROr")
            println(temp_endpt)
            println(endpt)
            println()
        end
    end #end for j=1:numInterval
    return newendpoints, _gap
end



function buddymatch(endpoints, V,y,m,s)
    row,col=size(endpoints)
    i=1
    while i!=row
        #buddy
#        println(i)
        row,col=size(endpoints)
        lower = endpoints[i,2]
        #    println("********  ",i+1,"  ",row)
        upper = endpoints[i+1,1]


        #    println("lower = ",lower,"  upper = ",upper)
        buddyIn_1 = false
        buddyIn_2=false
        for k=1:row
            row,col=size(endpoints)
            for j=1:2
                if endpoints[k,j] == 1-upper
                    buddyIn_1=true
                end
                if endpoints[k,j]==1-lower
                    buddyIn_2=true
                end
            end
        end
        if buddyIn_1 == false && buddyIn_2 == false
        #    println("[",1-upper,"  ",1-lower,"] by buddying [",lower,"  ",upper,"]")
            endpoints=[endpoints; [(1-upper) (1-lower)]]
            i=1
            row,col=size(endpoints)
        end
        endpoints=sort(collect(Iterators.flatten(endpoints)))
        endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
        endpoints = transpose(endpoints)
        if lower>=y && V==3 #in the two shares
            matchIn_1 = false
            matchIn_2=false
            for k=1:row
                row,col=size(endpoints)
                for j=1:2
                    if endpoints[k,j] ==m//s- upper
                        matchIn_1=true
                    end
                    if endpoints[k,j]==m//s-lower
                        matchIn_2=true
                    end
                end
            end
            if matchIn_1 == false && matchIn_2 == false
        #        println("[",m//s-upper,"  ",m//s-lower,"] by matching [",lower,"  ",upper,"]")
                endpoints=[endpoints; [(m//s-upper) (m//s-lower)]]
                i=1
                row,col=size(endpoints)
            end
        end
        endpoints=sort(collect(Iterators.flatten(endpoints)))
        endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
        endpoints = transpose(endpoints)
        i=i+1
#        println(i)
    end

    return endpoints
end
function Display(array)
    largest_denom = 0
    row,col=size(array)
    for i = 1:row
        for j=1:col
            if denominator(array[i,j])>largest_denom
                largest_denom = Int64(denominator(array[i,j]))
            end
        end
    end
    for i = 1:row
        for j=1:col
            if gcd(denominator(array[i,j]),Int64(largest_denom))!=denominator(array[i,j])
                largest_denom = largest_denom*(denominator(array[i,j])/(gcd(denominator(array[i,j]),Int64(largest_denom))))
            end
        end
    end
    matrix=Array{Any}(undef,0)
        for j=1:col
            for i=1:row
            num=numerator(array[i,j])*(largest_denom/denominator(array[i,j]))
            str = string(string(Int64(num))*"/"*string(Int64(largest_denom)))
            append!(matrix,[str])
        end
    end
    matrix = reshape(matrix, Int64(length(matrix)/col),col)
    display(matrix)
end

#VGAPV3(31,19,54//133)
#VGAPV3_proof(54,47,16//47)
#VGAPV3_proof(54,47,603//1786)
#println(VGAPV3_proof(54,47,603//1786))
#println(VGAPV3(54,47,839//2444))
#println(VGAPV3(54,47,603//1786))
#VGAPV3_proof(47,41,227//656)
#VGAPV3_proof(31,19,54//133)

#VGAPV3_proof(47,41,85//246)

#VGAPV3_proof(54,47,16//47)
#VGAPV3(54,47,85//246)
#VGAPV3_proof(47,41,85//246)
#VGAPV3_proof(50,39,209//611)
#VGAPV3_proof(29,23,143//414)


#VGAPV3_proof(29,23,39//115)

#VGAPV3_proof(32,25,17//50)
if false
println(VGAPV3(29,23,39//115)) #proposed alpha 49/138
println(VGAPV3(32,25,17//50))#proposed alpha 44/125
println(VGAPV3(31,27,37//108))
println(VGAPV3(52,31,89//217))
println(VGAPV3(38,33,34//99))
println(VGAPV3(39,34,47//136))
println(VGAPV3(62,37,91/222))
println(VGAPV3(43,39,53//156))
println(VGAPV3(50,39,9//26)) #proposed alpha 68/195
println(VGAPV3(67,40,197//480))
println(VGAPV3(47,41,85//246))
end
#VGAPV3_proof(31,24,151//432)
#VGAPV3_proof(32,25,17//50)
#VGAPV3_proof(11,10,37//110)
#VGAPV3_proof(68,53,37//106)
