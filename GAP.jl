##############################
# !!ASSUMES MID DID NOT WORK!!
##############################
include("GAP1.jl")
include("PROC.jl")#for binsearch

using JuMP
using GLPK

#binary search, used in GAP uses PROC and VGAP to find the
#alpha is a sorted array
    function linearSearch(m,s,array)
        for i=1:length(array)
    #        print(array[i])
    #        println("   ",VGAP(m,s,array[i]))
            if VGAP(m,s,array[i])==true
            #    GAP_proof(m,s,array[i])
                return array[i]
            end
        end
    end

    function binSearch(m,s,array)
        return binSearchHelp(m,s,array, 1, length(array))

    end
    function binSearchHelp(m,s,array,front,back)
    #    println("front: ",front,"  back:  ",back)
        if front>back
            return -1
        end
        guessIndex = front + Int64(floor((back-front)/2))

        guess = array[guessIndex]

        if bool_Proc == true && bool_Gap ==true
            return guess
        elseif bool_Proc ==true
            return binSearchHelp(m,s,array, guessIndex+1,back)
        else
            return binSearchHelp(m,s,array,front,guessIndex-1)
        end
    end



function GAP_proof(m,s,alpha)

    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    if V==3
        GAPV3(m,s,alpha)
        return
    end

    x,y=FINDEND(m,s,alpha,V)

    xbuddy = 1-x
    ybuddy = 1-y

    if (x==1-alpha && y==alpha) ||(y==1-alpha && x==alpha)
        println("all one interval alpha could be ≥ ",alpha)
        return false
    end
    if x > y
        println("intervals not disjoint alpha could be ≥ ", alpha)
        return false
    end
    endpoints = print_Intervals(m,s,alpha)
    if(V₋₁shares<Vshares)
         num_small_shares = V₋₁shares
        num_large_shares = Vshares - V₋₁shares
        S=sᵥ
        shares=Vshares
        VV=V #which is split
        num_split_shares=Int64(num_large_shares/2)

    else
        num_small_shares = V₋₁shares-Vshares
        num_large_shares =   Vshares
        S=sᵥ₋₁
        shares=V₋₁shares
        VV=V-1 #which is split
        num_split_shares=Int64(num_small_shares/2)

    end #******************************end else

    #println("Endpoints: ")
    #display(endpoints)
    row, col= size(endpoints)
    numIntervals = row
    X = perm(VV, numIntervals)
    possInd = Array{Int64}(undef,0)

    #Find the possible distribtions of muffins
    for i=1:length(X)
        A=X[i]
        sum_1=0
        sum_2=0
        for j =1:numIntervals
            sum_1=sum_1+A[j]*endpoints[j,1]
            sum_2=sum_2+A[j]*endpoints[j,2]
        end
        if (sum_1 < m//s) && (sum_2 > m//s)
            append!(possInd, i)
        end
    end

    if(length(possInd)==0)
        println("No possible distributions of muffins alpha ≤ ",alpha)
        return true
    else
            println("\nPossible distributions of muffins: ")
            for i in possInd
                println(X[i,:])
            end
    end

    while(true)

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
        _gap=false

        lowerbound = 0
        upperbound = 0
        for j=1:numIntervals
            lowerbound=endpoints[j,1]
            upperbound=endpoints[j,2]
            for i in possInd
                if(X[i][j]!=0)
                    #all types of students that use interval j
                    #X[i] is the student type
                    #sharesInIntervals[j] is the number of shares in that interval
                    #endpoints[j,:] is the endpoints of that interval

                    #lower bound
                    sum=0
                #    println("X[i][j]= ",X[i][j])
                #    println("X[i] = ",X[i])
                    for k=1:length(X[i])
                        if(k!=j)
                #            println("endpoints[",k,",2]*X[",i,"][",k,"] = ",endpoints[k,2]*X[i][k])
                            sum=sum+endpoints[k,2]*X[i][k]
                        end
                    end
                    lowerbound_temp=m//s-sum
                #    println("lowerbound_temp = ", lowerbound_temp)
                    if lowerbound_temp>lowerbound && lowerbound_temp < endpoints[j,2]
                        lowerbound=lowerbound_temp
                    end
                    #upper bound
                    sum=0
                    for k=1:length(X[i])
                        if(k!=j)
                            sum=sum+endpoints[k,1]*X[i][k]
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
            println("\nEndpoints after finding gaps: ")
            display(endpoints)
            println()

            if(_gap==false)
                println("No gaps found")
                println("alpha could be greater than ",alpha)
                return
            end

            #
            X = perm(VV, numIntervals)
            possInd = Array{Int64}(undef,0)
        #    display(X)
        #find possible students
            for i=1:length(X)
                A=X[i]
                sum_1=0
                sum_2=0
                for j =1:numIntervals
                    sum_1=sum_1+A[j]*endpoints[j,1]
                    sum_2=sum_2+A[j]*endpoints[j,2]
                end
                if (sum_1 < m//s) && (sum_2 > m//s)
                    append!(possInd, i)
                end
            end

            if length(possInd)==0
                println("No possible muffin distributions")
                println("alpha ≤ ",alpha)
                return
            end

            poss_Dist=reshape([],0,2) #possibe distributions of muffins
            for i in possInd
                if i == possInd[1]
                    poss_Dist=X[i,:]
                else
                    poss_Dist=[poss_Dist; X[i,:]]
                end
            end
            #get it in the shape I want
            poss_Dist = transpose(reshape(hcat(poss_Dist...), (length(poss_Dist[1]), length(poss_Dist))))
            println("Possible muffin distributions")
            display(poss_Dist)

            symmIntervals = Array{Int64,2}(undef,0,0) # row [i j] if interval i is symm to interval j
            row_endpt,col_endpt = size(endpoints)
            for i = 1:row_endpt
                for j=i+1:row_endpt
                    if ((endpoints[i,1]+endpoints[j,2])== 1 && endpoints[i,2] + endpoints[j,1] == 1) || ((endpoints[i,1]+endpoints[j,2])== m//s && endpoints[i,2] + endpoints[j,1] == m//s && V-1==2)
                        if length(symmIntervals) == 0
                            symmIntervals = [i j]
                        else
                            symmIntervals = [symmIntervals; i j]
                        end
                    end
                end
            end

            row_symm,col_symm=size(symmIntervals)
            A=Array{Int64}(undef,0,) # A will become the matrix in the eqation Ax=b
            for i=1:row_symm
                append!(A, (poss_Dist[:,symmIntervals[i,1]]-poss_Dist[:,symmIntervals[i,2]]))
            end
            A=reshape(A,Int64(length(A)/row_symm),row_symm)
            A=transpose(A)
            #A is now a matrix with rows correspoding to the symmetric identities
            #ex if the intervals are I1, I2, I3, I4, I5 (|I1|=|I4|,|I2|=|I3|)
            #A = [I1-I4 ; I2-I3]

            #Sum_1 is the sum of the left symmetric cols of possDist (these sum to num_split_shares)
            #Sum_2 is the sum of the rigth symmetric cols of possDist (these sum to num_split_shares)
            Sum_1=Array{Int64}(undef,0)
            Sum_2=Array{Int64}(undef,0)
            for i=1:row_symm
                if length(Sum_1)==0
                    Sum_1=poss_Dist[:,symmIntervals[i,1]]
                else
                    Sum_1=Sum_1+poss_Dist[:,symmIntervals[i,1]]
                end
                if length(Sum_2)==0
                    Sum_2=poss_Dist[:,symmIntervals[i,2]]
                else
                    Sum_2=Sum_2+poss_Dist[:,symmIntervals[i,2]]
                end
            end
            #transpose to get right dim
            Sum_1=Sum_1'
            Sum_2=Sum_2'

            #append as rows to A
            A=[A; Sum_1]
            A=[A;Sum_2]

            #Append a row of ones to A (num stud per distribtion adds to num VV students)
            row,col=size(A)
            Ones=(ones(Int64,col))'
            A=vcat(A, Ones)

            model=Model(with_optimizer(GLPK.Optimizer))
            @variable(model, x[i=1:length(possInd)],Int)
            row,col=size(symmIntervals)

            #add zeros for each row of I₁ - I₂
            b=zeros(Int64,row)

            b=[b;num_split_shares; num_split_shares;S]
            @constraint(model,con,A*x .==b)
            @constraint(model,con_1,x.>=0)
            println("System of equations (= ",b)
            display(A)
            optimize!(model)
            if(has_values(model))
                println()
                println("There is a solution on the Naturals")
                #println(value.(x))
                println("Looking for more gaps")
                println()
            else
                println("No solution on the Naturals")
                println("alpha ≤ ",alpha)
                return
            end
        end #end while

end
function VGAP(m,s,alpha)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    if V==3
        return VGAPV3(m,s,alpha)
    end
    x,y=FINDEND(m,s,alpha,V)

    xbuddy = 1-x
    ybuddy = 1-y
    if (x==1-alpha && y==alpha) ||(y==1-alpha && x==alpha)
        return false
    end
    if x > y
        return false
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
        endpoints=Array{Rational,2}(undef,0,0)
        endpoints =  [alpha ybuddy]
        endpoints = [endpoints; [xbuddy 1//2]; [1//2 x]]

    else
        num_small_shares = V₋₁shares-Vshares
        num_large_shares =   Vshares
        S=sᵥ₋₁
        shares=V₋₁shares
        VV=V-1 #which is split
        num_split_shares=Int64(num_small_shares/2)
        endpoints =Array{Rational, 2}(undef,0,0)
        endpoints =[y 1//2]
        endpoints = [endpoints;[1//2 ybuddy]; [xbuddy (1-alpha)]]
    #    sharesInIntervals = Array{Int64}(undef,0)
    #    append!(sharesInIntervals, [Int64(num_small_shares/2) Int64(num_small_shares/2) num_large_shares])

    #    for i = 1:numIntervals
    #        print(endpoints[i,:])
    #        println("  ",sharesInIntervals[i])
    #    end

end #******************************end else

    row, col= size(endpoints)
    numIntervals = row
    X = perm(VV, numIntervals)
    possInd = Array{Int64}(undef,0)

    #Find the possible distribtions of muffins
    for i=1:length(X)
        A=X[i]
        sum_1=0
        sum_2=0
        for j =1:numIntervals
            sum_1=sum_1+A[j]*endpoints[j,1]
            sum_2=sum_2+A[j]*endpoints[j,2]
        end
        if (sum_1 < m//s) && (sum_2 > m//s)
            append!(possInd, i)
        end
    end
    if(length(possInd)==0)
        return true
    end
    while(true)

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
        for j=1:numIntervals
            lowerbound=endpoints[j,1]
            upperbound=endpoints[j,2]
            for i in possInd
                if(X[i][j]!=0)
                    #all types of students that use interval j
                    #X[i] is the student type
                    #sharesInIntervals[j] is the number of shares in that interval
                    #endpoints[j,:] is the endpoints of that interval

                    #lower bound
                    sum=0
                #    println("X[i][j]= ",X[i][j])
                #    println("X[i] = ",X[i])
                    for k=1:length(X[i])
                        if(k!=j)
                #            println("endpoints[",k,",2]*X[",i,"][",k,"] = ",endpoints[k,2]*X[i][k])
                            sum=sum+endpoints[k,2]*X[i][k]
                        end
                    end
                    lowerbound_temp=m//s-sum
                #    println("lowerbound_temp = ", lowerbound_temp)
                    if lowerbound_temp>lowerbound && lowerbound_temp < endpoints[j,2]
                        lowerbound=lowerbound_temp
                    end
                    #upper bound
                    sum=0
                    for k=1:length(X[i])
                        if(k!=j)
                            sum=sum+endpoints[k,1]*X[i][k]
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

            if(_gap==false)
                return false
            end

            #
            X = perm(VV, numIntervals)
            possInd = Array{Int64}(undef,0)
        #    display(X)
        #find possible students
            for i=1:length(X)
                A=X[i]
                sum_1=0
                sum_2=0
                for j =1:numIntervals
                    sum_1=sum_1+A[j]*endpoints[j,1]
                    sum_2=sum_2+A[j]*endpoints[j,2]
                end
                if (sum_1 < m//s) && (sum_2 > m//s)
                    append!(possInd, i)
                end
            end

            if length(possInd)==0
                return true
            end

            poss_Dist=reshape([],0,2) #possibe distributions of muffins
            for i in possInd
                if i == possInd[1]
                    poss_Dist=X[i,:]
                else
                    poss_Dist=[poss_Dist; X[i,:]]
                end
            end
            #get it in the shape I want
            poss_Dist = transpose(reshape(hcat(poss_Dist...), (length(poss_Dist[1]), length(poss_Dist))))

            symmIntervals = Array{Int64,2}(undef,0,0) # row [i j] if interval i is symm to interval j
            row_endpt,col_endpt = size(endpoints)
            for i = 1:row_endpt
                for j=i+1:row_endpt
                    if (endpoints[i,1]+endpoints[j,2])==1 && endpoints[i,2]-endpoints[i,1]==endpoints[j,2]-endpoints[j,1]
                        if length(symmIntervals) == 0
                            symmIntervals = [i j]
                        else
                            symmIntervals = [symmIntervals; i j]
                        end
                    end
                end
            end

            row_symm,col_symm=size(symmIntervals)
            A=Array{Int64}(undef,0,) # A will become the matrix in the eqation Ax=b
            for i=1:row_symm
                append!(A, (poss_Dist[:,symmIntervals[i,1]]-poss_Dist[:,symmIntervals[i,2]]))
            end
            A=reshape(A,Int64(length(A)/row_symm),row_symm)
            A=transpose(A)
            #A is now a matrix with rows correspoding to the symmetric identities
            #ex if the intervals are I1, I2, I3, I4, I5 (|I1|=|I4|,|I2|=|I3|)
            #A = [I1-I4 ; I2-I3]

            #Sum_1 is the sum of the left symmetric cols of possDist (these sum to num_split_shares)
            #Sum_2 is the sum of the rigth symmetric cols of possDist (these sum to num_split_shares)
            Sum_1=Array{Int64}(undef,0)
            Sum_2=Array{Int64}(undef,0)
            for i=1:row_symm
                if length(Sum_1)==0
                    Sum_1=poss_Dist[:,symmIntervals[i,1]]
                else
                    Sum_1=Sum_1+poss_Dist[:,symmIntervals[i,1]]
                end
                if length(Sum_2)==0
                    Sum_2=poss_Dist[:,symmIntervals[i,2]]
                else
                    Sum_2=Sum_2+poss_Dist[:,symmIntervals[i,2]]
                end
            end
            #transpose to get right dim
            Sum_1=Sum_1'
            Sum_2=Sum_2'

            #append as rows to A
            A=[A; Sum_1]
            A=[A;Sum_2]

            #Append a row of ones to A (num stud per distribtion adds to num VV students)
            row,col=size(A)
            Ones=(ones(Int64,col))'
            A=vcat(A, Ones)

            model=Model(with_optimizer(GLPK.Optimizer))
            @variable(model, x[i=1:length(possInd)],Int)
            row,col=size(symmIntervals)

            #add zeros for each row of I₁ - I₂
            b=zeros(Int64,row)

            b=[b;num_split_shares; num_split_shares;S]
            @constraint(model,con,A*x .==b)
            @constraint(model,con_1,x.>=0)
            optimize!(model)
            if(!has_values(model))
                return true
            end
        end #end while
end


function GAP(m,s)
    if m%s==0
        return 1
    end
    array=Array{Rational}(undef,0)
    alph = 1//3
    num=1
    denom=3
    while denom<=m*s
        while alph<1//2
            append!(array,alph)
            num=num+1
            alph = num//denom
        #    println(alph)
    #        println(alph)
        end
        num=Int64(ceil(denom/3))
        denom = denom+1
        while (denom%s!=0)
        #    println(denom)
            denom=denom+1
        end
        alph = num//denom
    #    println(alph)
    end
    sort!(array)
    unique!(array)
#    println(length(array))
    #println(array)
    alpha = linearSearch(m,s,array)
#    println(alpha)
#    println(alpha)
    if alpha==-1
        return 1
    else
        return alpha
    end
end
#GAP(31,19)
#GAP(41,19)
GAP(59,22)
#GAP_proof(31,19,54//133)
#GAP_proof(41,19,131//304)
#GAP_proof(59,22,167//374)
#GAP_proof(41,23,149//368)
#GAP_proof(54,25,151//350)
#GAP_proof(67,25,223//500)
#GAP_proof(59,26,191//442)
#GAP_proof(47,29,117//290)
#GAP_proof(55,31,151//372)
#GAP_proof(67,31,187//434)
#GAP_proof(55,34,151//374)
#VGAP(55,34,151//374)
#VGAP(31,19,233//589)
#GAP_proof(31,19,208//513)
#GAP_proof(54,47,16//47)
#GAP_proof(41,19,204,475)
#perm(3,3)
