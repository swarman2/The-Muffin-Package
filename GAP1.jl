include("helper_functions.jl")
using JuMP
using Cbc

function VGAPV3(m,s,alpha, proof=false)

    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    x,y=FINDEND(m,s,alpha,V)
    ybuddy=1-y
    xbuddy=1-x

    denom=denominator(alpha)
    denom=lcm(s,denom)*2
    if proof
        println()
        println("PROOF OF f(",m,", ",s,") = ",alpha)
        println("There are ",V,"-students and ",V-1,"-students")
        println("s",V,": ",sᵥ,"   s",V-1,": ",sᵥ₋₁)
        println("All numbers assumed to have denominator ",denom)
    end

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
        #endpoints = [endpoints; [xbuddy x]]
        I1 = num_small_shares
        I2 = Int64(num_large_shares/2)
        I3 = I2
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

    end #******************************end else
    #print_Intervals(m,s,alpha)
    endpoints=[alpha*denom x*denom]
    endpoints=[endpoints;y*denom (1-alpha)*denom; xbuddy*denom ybuddy*denom]
    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    if proof
        println("ENDPOINTS BEFORE BUDDY MATCH")
        endpoints=Int(endpoints)
        display(endpoints)
    end
    if x>=y
        if proof
            println("INVALID INTERVALS")
        end
        return false
    end
    endpoints=[endpoints; denom//2 denom//2 ]
    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    #display(endpoints)
    if proof
        println("\nSPLIT AT 1/2 (",Int64(denom),") and  m//2s (",Int64(denom*m//2s),")")
        #endpoints = Int(endpoints)
        #display(endpoints)
    end
    endpoints = buddymatch(endpoints,V,y,m,s,denom,proof)
    #display(endpoints)
    row,col=size(endpoints)
    if V==3
        endpoints = [endpoints; [(m*denom//2s) (m*denom//2s)];[(m*denom//2s) (m*denom//2s)]]

        endpoints=sort(collect(Iterators.flatten(endpoints)))
        endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
        endpoints = transpose(endpoints)
        if proof
        #    println("\nSPLIT AT m//2s (",Int64(denom*m//2s),")")
            endpoints = Int(endpoints)
            display(endpoints)
        end
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
                    if endpoints[k,j] == denom-upper
                        buddyIn_1=true
                    end
                    if endpoints[k,j]==denom-lower
                        buddyIn_2=true
                    end
                end
            end
            if buddyIn_1 == false && buddyIn_2 == false && denom-upper <= endpoints[1,1] && denom-lower >= endpoints[row,col]
                if proof
                    println("[",Int64(denom-upper),"  ",Int64(denom-lower),"] by buddying [",Int64(lower),"  ",Int64(upper),"]")
                end
                endpoints=[endpoints; [(denom-upper) (denom-lower)];  [(denom-upper) (denom-lower)]]
                i=1
                row,col=size(endpoints)
            end
            endpoints=sort(collect(Iterators.flatten(endpoints)))
            endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
            endpoints = transpose(endpoints)
            if lower>=y*denom && V==3 #in the two shares
                matchIn_1 = false
                matchIn_2=false
                for k=1:row
                    row,col=size(endpoints)
                    for j=1:2
                        if endpoints[k,j] ==(m//s)*denom- upper
                            matchIn_1=true
                        end
                        if endpoints[k,j]==(m//s)*denom-lower
                            matchIn_2=true
                        end
                    end
                end
                if matchIn_1 == false && matchIn_2 == false && (m//s)*denom-upper <= endpoints[1,1] && (m//s)*denom-lower >= endpoints[row,col]
                    if proof
                        println("[",Int64((m//s)*denom-upper),"  ",Int64((m//s)*denom-lower),"] by matching [",Int64(lower),"  ",Int64(upper),"]")
                    end
                    endpoints=[endpoints; [((m//s)*denom-upper) ((m//s)*denom-lower)];  [((m//s)*denom-upper) ((m//s)*denom-lower)]]
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
        endpoints = buddymatch(endpoints,V,y,m,s,denom,proof)

        row, col= size(endpoints)
        numIntervals = row

        numVIntervals = 0
        for i =1: row
            if endpoints[i,2]<=x*denom
                numVIntervals = numVIntervals+1
            end
        end
        permV=perm(V, numVIntervals)
        permV_ = perm(V-1, row - numVIntervals)


        possV = Vector{Vector{Int64}}()
        possV_ = Vector{Vector{Int64}}()
        if proof
            endpoints = Int(endpoints)
            println("\n",V," intervals")
            for j=1:numVIntervals
                println("J_",j,": ",endpoints[j,:])
            end
            println("\n",V-1," intervals")
            for j=numVIntervals+1: row
                println("I_",j-numVIntervals,": ",endpoints[j,:])
            end
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
            if (sum_1 < (m//s)*denom) && (sum_2 > (m//s)*denom)
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
            if (sum_1 < (m//s)*denom) && (sum_2 > (m//s)*denom)
                append!(possV_, permV_[i,:])
            end
        end
        if proof
            println("\nPossible ",V," students:")
            for i=1:length(possV)
                PRINT(possV[i])
                println()
            end
            println("\nPossible ",V-1," students:")
            for i=1:length(possV_)
                PRINT(possV_[i])
                println()
            end
                println()
        end

        mat_V=transpose(hcat(possV...))
        #    endpoints, gap_found = findGaps(mat_V, endpoints,m,s)
        mat_V_=transpose(hcat(possV_...))
        endpoints, gap_found = findGaps(mat_V,mat_V_, endpoints,m,s,denom,proof)
    end #end finding gaps

    if proof
        println("NO MORE GAPS")
        #println("ENDPOINTS AFTER BUDDY MATCH")
        #display(endpoints)
    end
    row, col= size(endpoints)
    numIntervals = row

    numVIntervals = 0
    for i =1: row
        if endpoints[i,2]<=x*denom
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
        if (sum_1 < (m//s)*denom) && (sum_2 > (m//s)*denom)
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
        if (sum_1 < (m//s)*denom) && (sum_2 > (m//s)*denom)
            append!(possV_, permV_[i,:])
        end
    end
    mat_V=transpose(hcat(possV...))
    mat_V_=transpose(hcat(possV_...))

    symmIntervals = Array{Int64,2}(undef,0,0) # row [i j] if interval i is symm to interval j
    row_endpt,col_endpt = size(endpoints)
    for i = 1:row_endpt
        for j=i+1:row_endpt
            if (endpoints[i,1]+endpoints[j,2])==denom && endpoints[i,2]-endpoints[i,1]==endpoints[j,2]-endpoints[j,1] || (endpoints[i,1]+endpoints[j,2])==(m*denom)//s && endpoints[i,2]-endpoints[i,1]==endpoints[j,2]-endpoints[j,1]
                if length(symmIntervals) == 0
                    symmIntervals = [i j]
                else
                    symmIntervals = [symmIntervals; i j]
                end
            end
        end
    end
    #println("\nSYMMETRIC INTERVALS")
    #display(symmIntervals)
    row_symm, col_symm = size(symmIntervals)
    if proof
        println("\nPossible distributions of students for the ",V," shares:")
        for i=1:length(possV)
            PRINT(possV[i])
            println()
        end
        println("\nPossible distributions of students for the ",V-1," shares:")
        for i=1:length(possV_)
            PRINT(possV_[i])
            println()
        end
    end
    if length(mat_V)==0
        if proof
            println("\nNo possible distributions of ",V," shares: alpha ≤ ",alpha)
        end
        return true
    elseif length(mat_V_)==0
        if proof
            println("\nNo possible distributions of ",V-1," shares: alpha ≤ ",alpha)
        end
        return true
    end
    row_mat,col_mat = size(mat_V)
    row_mat_,col_mat_=size(mat_V_)
    row_mat_new=0
    while(row_mat_new < row_mat + row_mat_)
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

    #println("\nConcationate them in order to make one system of equations")
    #display(poss_Dist)

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


    #println("\nIf intervals i and j are symmetric, subtract col i from col j and add as a row to A")
    A=unique(A, dims = 1)
    #display(A)
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
    #println("\nAdd a row of ",V,"'s and a row of ",V-1,"'s \n(these when mulitplied by number of each type of distribution of students will add to number of ",V," students and number of ",V-1," students)")
    #append as rows to A
    #display(A)
    #Append a row of ones to A (num stud per distribtion adds to num VV students)
    Ones=(ones(Int64,col))'

    split_Intervals_1 = Array{Int64}(undef,0)
    split_Intervals_2 = Array{Int64}(undef,0)
    if VV==V
        row,col=size(endpoints)
        for i = 1:row
            if endpoints[i,1]>=xbuddy*denom && endpoints[i,2]<=denom//2
                append!(split_Intervals_1, i)
            end
            if endpoints[i,1]>=denom//2 && endpoints[i,2]<=x*denom
                append!(split_Intervals_2,i)
            end
        end
    else
        row,col=size(endpoints)
        for i = 1:row
            if endpoints[i,1]>=y*denom && endpoints[i,2]<=denom//2
                append!(split_Intervals_1, i)
            end
            if endpoints[i,1]>=denom//2 && endpoints[i,2]<=ybuddy*denom
                append!(split_Intervals_2,i)
            end
        end
    end
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
    if length(symm_sum_1)!=0
        A=[A;symm_sum_1]
        b=[b;num_split_shares]
    end
    if length(symm_sum_2) !=0
        A=[A;symm_sum_2]
        b=[b; num_split_shares]
    end
    #println("\nAdd a rows for the sums of the ",VV," shares symmetric around 1//2")
    #display(A)
    if proof
        println("Then solve the system for the number of each type of distribtuion (= ",b,")")
    end

    row,col=size(A)
    #display(A)
    model=Model(with_optimizer(Cbc.Optimizer, logLevel=0))
    @variable(model, x[i=1:col],Int)

    @constraint(model,con,A*x .==b)
    @constraint(model,con_1,x.>=0)

    optimize!(model)
    if(termination_status(model)==MOI.OPTIMAL)
        if proof
            display(value.(x))
            println("\nThere is a solution on the naturals, alpha > ",alpha)
        end
        return false
    else
        if proof
            println("\nThere is no solution on the naturals, alpha ≤ ",alpha)
        end
        return true
    end
end

function findGaps(possDistV,possDistV_,endpoints,m,s,denom,proof = false)
    #find gaps
    _gap=false
    numVDistributions,numVIntervals=size(possDistV)
    numV_Distributions, numV_Intervals=size(possDistV_)

    newendpoints=endpoints
    #temp_endpt = Vector{Rational}(undef,0)
    #lowEnd = Array{Rational}(undef,0)
    #highEnd = Array{Rational}(undef,0)
    for j=1:numVIntervals
        #println("\nChecking interval: ",j)
        lowEnd = Array{Rational}(undef,0)
        highEnd = Array{Rational}(undef,0)
        for i=1:numVDistributions
            if(possDistV[i,j]!=0)
                #lower bound
                sum=0
                for k=1:length(possDistV[i,:])
                    if(k!=j)
                        sum=sum+endpoints[k,2]*possDistV[i,k]
                    else
                        sum=sum+endpoints[k,2]*(possDistV[i,k]-1)
                    end
                end
                lowerbound_temp=(m//s)*denom-sum
                #upper bound
                sum=0
                for k=1:length(possDistV[i,:])
                    if(k!=j)
                        sum=sum+endpoints[k,1]*possDistV[i,k]
                    else
                        sum=sum+endpoints[k,1]*(possDistV[i,k]-1)
                    end
                end
                upperbound_temp=(m//s)*denom-sum
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
        temp_endpt = Vector{Rational}(undef,0)
        #only consider unique lower bounds
        lowEnd_temp = Array{Rational}(undef,0)
        for i=1:length(lowEnd)
            append!(lowEnd_temp,lowEnd[i])
        end
        #println(lowEnd)
        highEnd_temp = Array{Rational}(undef,0)
        for i=1:length(highEnd)
            append!(highEnd_temp,highEnd[i])
        end
        while unique(lowEnd_temp)!=lowEnd_temp
            for i =1:length(lowEnd_temp)
                for j=i+1:length(lowEnd_temp)
                    if lowEnd_temp[i]==lowEnd_temp[j]
                        if highEnd_temp[i]>=highEnd_temp[j]
                            deleteat!(lowEnd_temp,j)
                            deleteat!(highEnd_temp,j)
                        else
                            deleteat!(lowEnd_temp,i)
                            deleteat!(highEnd_temp,i)
                        end
                        break
                    end
                end
            end
        end

        while (length(lowEnd_temp)>0)
            lowIndex=1
            #println(lowEnd*230)
            for i = 1:length(lowEnd_temp)
                if lowEnd_temp[i]<lowEnd_temp[lowIndex]
                    lowIndex=i
                end
            end
            push!(temp_endpt,lowEnd_temp[lowIndex])
            push!(temp_endpt,highEnd_temp[lowIndex])
            deleteat!(lowEnd_temp,lowIndex)
            deleteat!(highEnd_temp,lowIndex)
        end
        merge=true
        while merge==true
            merge = false
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
        endpt=filter(x-> x > endpoints[j,1],temp_endpt)
        endpt=filter(x-> x < endpoints[j,2],endpt)
        int_lowEnd = Array{Int64}(undef,length(lowEnd))
        for i=1:length(lowEnd)
            int_lowEnd[i] = Int64(lowEnd[i])
        end
        lowEnd=int_lowEnd
        int_highEnd = Array{Int64}(undef,length(highEnd))
        for i=1:length(highEnd)
            int_highEnd[i] = Int64(highEnd[i])
        end
        highEnd=int_highEnd
        if length(endpt)%2==0
            for i=1:Int64(length(endpt)/2)
                upperbound=Int64(endpt[i])
                lowerbound=Int64(endpt[2i])
                if lowerbound>upperbound
                    if proof
                        println("NEW GAP FOUND: [",upperbound,",",lowerbound,"]")
                        k=1
                        for i=1:numVDistributions
                            if(possDistV[i,j]!=0)
                            #    println(possDistV[i,:])
                                print(lowEnd[k]," ≤ ")
                                PRINT(possDistV[i,:])
                                println(" ≤ ",highEnd[k])
                                k=k+1
                            end
                        end
                    end
                    _gap=true
                    newendpoints=[newendpoints; upperbound lowerbound]#; (1-lowerbound) (1-upperbound)]
                    newendpoints=sort(collect(Iterators.flatten(newendpoints)))
                    newendpoints=reshape(newendpoints,(2,Int64(length(newendpoints)/2)))
                    newendpoints = transpose(newendpoints)
                end #end if(lowerbound>upperbound)
            end#end for i = 1:Int64(length(temp_endpt)/2)
        elseif length(endpt)==1

            if temp_endpt[1]<endpoints[j,1] #then change upper bound
                _gap=true
                if proof
                    println("GAP EXTENDED: ",Int64(newendpoints[j,2])," changed to ",Int64(temp_endpt[2]))
                end
                newendpoints[j,2]=temp_endpt[2]
                if proof
                    k=1
                    for i=1:numVDistributions
                        if(possDistV[i,j]!=0)
                            print(lowEnd[k]," ≤ ")
                            PRINT(possDistV[i,:])
                            println(" ≤ ",highEnd[k])
                            k=k+1
                        end
                    end
                end
            elseif temp_endpt[1]>endpoints[j,1]
                _gap=true
                if proof
                    println("GAP EXTENDED: ",Int64(newendpoints[j,1])," changed to ",Int64(temp_endpt[1]))
                end
                newendpoints[j,1]=temp_endpt[1]
                if proof
                    k=1
                    for i=1:numVDistributions
                        if(possDistV[i,j]!=0)
                            print(lowEnd[k]," ≤ ")
                            PRINT(possDistV[i,:])
                            println(" ≤ ",highEnd[k])
                            k=k+1
                        end
                    end
                end
            end

        else
            println()
            println("ERROR")
            println()
        end
    end #end for j=1:numInterval

   #move on to the V-1 Intervals
    for j=numVIntervals+1:numVIntervals + numV_Intervals
        lowEnd = Array{Rational}(undef,0)
        highEnd = Array{Rational}(undef,0)
        for i=1:numV_Distributions
            if(possDistV_[i,j-numVIntervals]!=0)
                #lower bound
                sum=0
                for k=numVIntervals+1:numVIntervals+length(possDistV_[i,:])
                    if(k!=j)
                        sum=sum+endpoints[k,2]*possDistV_[i,k-numVIntervals]
                    else
                        sum=sum+endpoints[k,2]*(possDistV_[i,k-numVIntervals]-1)
                    end
                end
                lowerbound_temp=(m//s)*denom-sum
                #upper bound
                sum=0
                for k=numVIntervals+1:numVIntervals+length(possDistV_[i,:])
                    if(k!=j)
                        sum=sum+endpoints[k,1]*possDistV_[i,k-numVIntervals]
                    else
                        sum=sum+endpoints[k,1]*(possDistV_[i,k-numVIntervals]-1)
                    end

                end
                upperbound_temp=(m//s)*denom-sum
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
        lowEnd_temp = Array{Rational}(undef,0)
        for i=1:length(lowEnd)
            append!(lowEnd_temp,lowEnd[i])
        end
        #println(lowEnd)
        highEnd_temp = Array{Rational}(undef,0)
        for i=1:length(highEnd)
            append!(highEnd_temp,highEnd[i])
        end
        while unique(lowEnd_temp)!=lowEnd_temp
            for i =1:length(lowEnd_temp)
                for j=i+1:length(lowEnd_temp)
                    if lowEnd_temp[i]==lowEnd_temp[j]
                        if highEnd_temp[i]>=highEnd_temp[j]
                            deleteat!(lowEnd_temp,j)
                            deleteat!(highEnd_temp,j)
                        else
                            deleteat!(lowEnd_temp,i)
                            deleteat!(highEnd_temp,i)
                        end
                        break
                    end
                end
            end
        end
        temp_endpt=Array{Rational}(undef,0)
        while (length(lowEnd_temp)>0)
            lowIndex=1
            for i = 1:length(lowEnd_temp)
                if lowEnd_temp[i]<lowEnd_temp[lowIndex]
                    lowIndex=i
                end
            end
            push!(temp_endpt,lowEnd_temp[lowIndex])
            push!(temp_endpt,highEnd_temp[lowIndex])
            deleteat!(lowEnd_temp,lowIndex)
            deleteat!(highEnd_temp,lowIndex)
        end
        merge=true
        while merge==true
            merge = false
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
        endpt=filter(x-> x > endpoints[j,1],temp_endpt)
        endpt=filter(x-> x < endpoints[j,2],endpt)
        int_lowEnd = Array{Int64}(undef,length(lowEnd))
        for i=1:length(lowEnd)
            int_lowEnd[i] = Int64(lowEnd[i])
        end
        lowEnd=int_lowEnd
        int_highEnd = Array{Int64}(undef,length(highEnd))
        for i=1:length(highEnd)
            int_highEnd[i] = Int64(highEnd[i])
        end
        highEnd=int_highEnd
        #println("-----------------------",endpt*246)
        if length(endpt)%2==0
            for i=1:Int64(length(endpt)/2)
                upperbound=Int64(endpt[i])
                lowerbound=Int64(endpt[2i])
                if lowerbound>upperbound
                    if proof
                        println("NEW GAP FOUND: [",upperbound,",",lowerbound,"]")
                        k=1
                        for i=1:numV_Distributions
                            if(possDistV_[i,j-numVIntervals]!=0)
                                print(lowEnd[k]," ≤ ")
                                PRINT(possDistV_[i,:])
                                println(" ≤ ",highEnd[k])
                                k=k+1
                            end
                        end
                    end
                    _gap=true
                    newendpoints=[newendpoints; upperbound lowerbound]#; (1-lowerbound) (1-upperbound)]
                    newendpoints=sort(collect(Iterators.flatten(newendpoints)))
                    newendpoints=reshape(newendpoints,(2,Int64(length(newendpoints)/2)))
                    newendpoints = transpose(newendpoints)
                end #end if(lowerbound>upperbound)
            end#end for i = 1:Int64(length(temp_endpt)/2)
        elseif length(endpt)==1
            if temp_endpt[1]<endpoints[j,1] #then change upper bound
                _gap=true
                if proof
                    println("GAP EXTENDED: ",newendpoints[j,2]," changed to ",temp_endpt[2])
                end
                endpoints[j,2]=temp_endpt[2]
                k=1
                if proof
                    for i=1:numV_Distributions
                        if(possDistV_[i,j-numVIntervals]!=0)
                            print(lowEnd[k]," ≤ ")
                            PRINT(possDistV_[i,:])
                            println(" ≤ ",highEnd[k])
                            k=k+1
                        end
                    end
                end
            elseif temp_endpt[1]>endpoints[j,1]
                _gap=true
                if proof
                    println("GAP EXTENDED: ",newendpoints[j,1]," changed to ",temp_endpt[1])
                end
                endpoints[j,1]=temp_endpt[1]
                k=1
                if proof
                    for i=1:numV_Distributions
                        if(possDistV_[i,j-numVIntervals]!=0)
                            print(lowEnd[k]," ≤ ")
                            PRINT(possDistV_[i,:])
                            println(" ≤ ",highEnd[k])
                            k=k+1
                        end
                    end
                end
            end

        else
            println()
            println("ERROr")
            println()
        end
    end #end for j=1:numInterval
    return newendpoints, _gap
end

function buddymatch(endpoints, V,y,m,s,denom,proof = false)
    row,col=size(endpoints)
    i=1
    while i!=row
        #buddy
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
                if endpoints[k,j] == denom-upper
                    buddyIn_1=true
                end
                if endpoints[k,j]==denom-lower
                    buddyIn_2=true
                end
            end
        end
        if buddyIn_1 == false && buddyIn_2 == false && denom-upper >=endpoints[1,1] && denom-lower <= endpoints[row,col]
            if proof
                println("[",Int64(denom-upper),"  ",Int64(denom-lower),"] by buddying [",Int64(lower),"  ",Int64(upper),"]")
            end
            endpoints=[endpoints; [(denom-upper) (denom-lower)]]
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
                    if endpoints[k,j] ==(m//s)*denom- upper
                        matchIn_1=true
                    end
                    if endpoints[k,j]==(m//s)*denom-lower
                        matchIn_2=true
                    end
                end
            end
            if matchIn_1 == false && matchIn_2 == false && (m//s)*denom-upper <= endpoints[1,1] && (m//s)*denom-lower >= endpoints[row,col]
                if proof
                    println("[",Int64((m//s)*denom-upper),"  ",Int64((m//s)*denom-lower),"] by matching [",Int64(lower),"  ",Int64(upper),"]")
                end
                endpoints=[endpoints; [((m//s)*denom-upper) ((m//s)*denom-lower)]]
                i=1
                row,col=size(endpoints)
            end
        end
        endpoints=sort(collect(Iterators.flatten(endpoints)))
        endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
        endpoints = transpose(endpoints)
        i=i+1
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
function PRINT(array)
    print("( ")
    for i=1:length(array)
        for j=1:array[i]
            print(i," ")
        end
    end
    print(")")
end
function Int(matrix)
    row,col=size(matrix)
    mat = Matrix{Int64}(undef, row,col)
    for i=1:row
        for j=1:col
            mat[i,j]=Int64(matrix[i,j])
        end
    end
    return mat
end

if false
println(VGAPV3(29,23,39//115)) #proposed alpha 49/138
println(VGAPV3(32,25,17//50))#proposed alpha 44/125
println(VGAPV3(31,27,37//108))
println(VGAPV3(52,31,89//217))
println(VGAPV3(38,33,34//99))
println(VGAPV3(39,34,47//136))
println(VGAPV3(62,37,91//222))
println(VGAPV3(43,39,53//156))
println(VGAPV3(50,39,9//26)) #proposed alpha 68/195
println(VGAPV3(67,40,197//480))
println(VGAPV3(47,41,85//246))
end
#VGAPV3(38,33,34//99,true)
#VGAPV3(29,23,215//644,true)
#VGAPV3_proof(32,25,17//50)
#VGAPV3_proof(11,10,37//110)
#VGAPV3_proof(68,53,37//106)
#VGAPV3(29,23,172//483,true)
