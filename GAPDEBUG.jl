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

    if x>y
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

        println("TEST")
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
    #println("ENDPOINTS BEFORE BUDDY MATCH")
    #Display(endpoints)
#    endpoints=buddymatch(endpoints,V,y,m,s)
row,col=size(endpoints)
for i=1:row
    if endpoints[i,1]+ endpoints[i,2]==1
        endpoints=[endpoints; 1//2 1//2 ]
    end

end
endpoints=sort(collect(Iterators.flatten(endpoints)))
endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
endpoints = transpose(endpoints)
println("ENDPOINTS BEFORE BUDDY MATCH")
Display(endpoints)

        possV = Vector{Vector{Int64}}()
        possV_ = Vector{Vector{Int64}}()
numVIntervals = 0
for i =1: row
    if endpoints[i,2]<=x
        numVIntervals = numVIntervals+1
    end
end
row,col=size(endpoints)
permV=perm(V, numVIntervals)
permV_ = perm(V-1, row - numVIntervals)
println("\n",V-1," intervals")
for j=numVIntervals+1: row
    println("I_",j-numVIntervals,": ",endpoints[j,:])
end
#Find the possible distribtions of muffins
println(permV_)
for i=1:length(permV_)
    A=permV_[i]
    sum_1=0
    sum_2=0
    for j =numVIntervals+1:row
        sum_1=sum_1+A[j-numVIntervals]*endpoints[j,1]
        sum_2=sum_2+A[j-numVIntervals]*endpoints[j,2]
    end
    println(A)
    println("Sum_1: ",sum_1,"  Sum_2:  ",sum_2)
    if (sum_1 < m//s) && (sum_2 > m//s)
        append!(possV_, permV_[i,:])
    end
end
println("\nPossible students:")
println(possV_)
endpoints, gap_found = findGaps(mat_V_, endpoints,m,s)
Display(endpoints)

#endpoints = buddymatch(endpoints,V,y,m,s)
    row,col=size(endpoints)

    gap_found=true
    gap_found2=true
    while(gap_found == true || gap_found2 == true)
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
        println("\nPossible students:")
        println(possV)
        mat_V=transpose(hcat(possV...))
        endpoints, gap_found = findGaps(mat_V, endpoints,m,s)
        mat_V_=transpose(hcat(possV_...))
        endpoints, gap_found2 = findGaps(mat_V_, endpoints,m,s)
    end #end finding gaps
    row,col=size(endpoints)

    #for i=1:row
    #    if endpoints[i,1]+endpoints[i,2] == (m//s)
            endpoints = [endpoints; [(m//2s) (m//2s)];[(m//2s) (m//2s)]]
    #    end
    #end

    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    row,col=size(endpoints)

    #special buddy match for closed interval
    for i=1:row
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
    row_symm, col_symm = size(symmIntervals)
    #display(symmIntervals)
    for i=1:row_symm
        for k=1:2
            val = symmIntervals[i,k]
            for j=1:row_symm
                if j!=i && symmIntervals[j,1]==val
                    symmIntervals = [symmIntervals; symmIntervals[i,k%2+1] symmIntervals[j,2]]
                end
                if j!=i && symmIntervals[j,2]==val
                    symmIntervals = [symmIntervals; symmIntervals[i,k%2+1] symmIntervals[j,1]]
                end
            end
        end
    end
    #display(symmIntervals)
    sort!(symmIntervals,dims=2)
    symmIntervals = unique(symmIntervals,dims =2)
    println("\nSYMMETRIC INTERVALS")
    display(symmIntervals)
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
    display(split_Intervals_1)
    display(split_Intervals_2)
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
    A=[A;symm_sum_1]
    A=[A;symm_sum_2]
    b=[b;num_split_shares; num_split_shares]
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


function findGaps(possDist,endpoints,m,s)
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
    lowerbound = 0
    upperbound = 0
    _gap=false
    numDistributions,numIntervals=size(possDist)
    #display(possDist)
    for j=1:numIntervals
        lowerbound=endpoints[j,1]
        upperbound=endpoints[j,2]
        for i=1:numDistributions
            if(possDist[i,j]!=0)
                #all types of students that use interval j
                #X[i] is the student type
                #sharesInIntervals[j] is the number of shares in that interval
                #endpoints[j,:] is the endpoints of that interval

                #lower bound
                sum=0
            #    println("X[i][j]= ",X[i][j])
            #    println("X[i] = ",X[i])
                for k=1:length(possDist[i,:])
                    if(k!=j)
            #            println("endpoints[",k,",2]*X[",i,"][",k,"] = ",endpoints[k,2]*X[i][k])
                        sum=sum+endpoints[k,2]*possDist[i,k]
                    end
                end
                lowerbound_temp=m//s-sum
            #    println("lowerbound_temp = ", lowerbound_temp)
                if lowerbound_temp>lowerbound && lowerbound_temp < endpoints[j,2]
                    lowerbound=lowerbound_temp
                end
                #upper bound
                sum=0
                for k=1:length(possDist[i,:])
                    if(k!=j)
                        sum=sum+endpoints[k,1]*possDist[i,k]
                    end
                end
                upperbound_temp=m//s-sum
                if upperbound_temp <upperbound && upperbound_temp > endpoints[j,1]
                    upperbound = upperbound_temp
                end
            end #end if(X[i][j]!=0)
        end  #end  for i in possInd

        #check to see if the lowerbound and upper bound are valid for a gap

        if(lowerbound>upperbound && lowerbound>endpoints[j,1] && upperbound < endpoints[j,2])
            _gap=true
            endpoints=[endpoints; lowerbound upperbound; (1-lowerbound) (1-upperbound)]

            endpoints=sort(collect(Iterators.flatten(endpoints)))
            endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
            endpoints = transpose(endpoints)
            numIntervals=numIntervals+2
        end #end if(lowerbound>upperbound)
    end #end for j=1:numInterval
    return endpoints, _gap
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
            #println("[",m//s-(1-lower),"  ",(m//s)-(1-upper),"] by matching [",1-upper,"  ",1-lower,"]")
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
            #    println("[",m//s-upper,"  ",m//s-lower,"] by matching [",lower,"  ",upper,"]")
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
                largest_denom = denominator(array[i,j])
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
    matrix = reshape(matrix, Int64(length(matrix)/2),2)
    display(matrix)
end

#VGAPV3(31,19,54//133)
#VGAPV3(54,47,16//47)
#VGAPV3_proof(54,47,603//1786)
#println(VGAPV3_proof(54,47,603//1786))
#println(VGAPV3(54,47,839//2444))
#println(VGAPV3(54,47,603//1786))
#VGAPV3_proof(47,41,227//656)

VGAPV3_proof(47,41,85//246)
#VGAPV3_proof(54,47,16//47)
