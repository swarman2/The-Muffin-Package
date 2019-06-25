##############################
# !!ASSUMES MID DID NOT WORK!!
##############################


include("HALF.jl") #for SV and FINDEND
include("src\\permutations.jl") #for multiset_permutations
using JuMP
using GLPK
function perm(n,r)
    A=Array{Int64,1}(undef,r*(n+1))
    for i=1:length(A)
        A[i]=floor((i-1)/r)
    end
#    print(A)

    X=(collect(multiset_permutations(A,r)))
#    display(X)
#    println("********************")
    X=filter(x->sum(x)==n,X)
    return X

end
function GAP_proof(m,s,alpha)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁

    x,y=FINDEND(m,s,alpha,V)

    xbuddy = 1-x
    ybuddy = 1-y
    if(V₋₁shares<Vshares)
         #___(_______)________(_____)___|___(_____)____
         #alpha   y-buddy   xbuddy   x  |   y    1-alpha
         #smallShare is an array that holds the possiblities for number of small shares
         num_small_shares = V₋₁shares
        num_large_shares = Vshares - V₋₁shares
        println("m  = ",m,"  s = ", s)
        println( sᵥ," ",V,"-students \t",sᵥ₋₁," ",V-1,"-students \t",sᵥ*V," ",V,"-shares \t",sᵥ₋₁*(V-1)," ",V-1,"-shares")
        println()
        println("SPLIT THE ",V," SHARES")
        println("   ( ",num_small_shares," small ",V,"-shares)   (",num_large_shares," large ",V,"-shares)   |   ( ",V₋₁shares,"  ",V-1,"-shares )    ")
        println("  ",alpha,"     ", ybuddy,"  ", xbuddy,"          ",x," | ", y,"              ",1-alpha)
        println()
        println("SPLIT THE ",V," SHARES AGAIN")
        println("   ( ",num_small_shares," " ,V,"-shares)           (",Int64(num_large_shares/2)," large ",V,"-shares |  ",Int64(num_large_shares/2)," large ",V,"-shares)    ")
        println("  ",alpha,"     ", ybuddy,"     ", xbuddy,"          1/2      ",x)
        println()
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
        println("m  = ",m,"  s = ", s)
        println( sᵥ," ",V,"-students \t",sᵥ₋₁," ",V-1,"-students \t",sᵥ*V," ",V,"-shares \t",sᵥ₋₁*(V-1)," ",V-1,"-shares")
        println()
        println("SPLIT THE ",V-1," SHARES")
        println("   ( ",Vshares," ",V,"-shares)  | (",num_small_shares," small ",V-1,"-shares)      ( ",num_large_shares," large ",V-1,"-shares )    ")
        println("  ",alpha,"     ", x," | ", y,"          ",ybuddy,"  ", xbuddy,"              ",1-alpha)
        println()
        println("SPLIT THE ",V-1," SHARES AGAIN")
        println("    (",Int64(num_small_shares/2),"  ",V-1,"-shares | ",Int64(num_small_shares/2),"  ",V-1,"-shares)      ( ",num_large_shares," large ",V-1,"-shares )    ")
        println("     ", y,"    1/2 ", ybuddy,"          ",xbuddy, "        ",1-alpha)
        println()
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
    println("Endpoints: ")
    display(endpoints)
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

    println("\nPossible distributions of muffins: ")
    for i in possInd
        println(X[i,:])
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

    x,y=FINDEND(m,s,alpha,V)

    xbuddy = 1-x
    ybuddy = 1-y
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
                return false
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
VGAP(55,34,151//374)
