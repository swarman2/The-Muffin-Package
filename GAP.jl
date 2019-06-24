
function perm(n,r)
    A=Array{Int64,1}(undef,n*(n+1))
    for i=1:length(A)
        A[i]=floor((i-1)/n)
    end
    #print(A)

    X=(collect(multiset_permutations(A,r)))
#    print(X)
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
        endpoints=Array{Rational}(undef,0)
        append!(endpoints, [alpha ybuddy xbuddy 1//2 1//2 x])
        I1=num_small_shares
        I2=Int64(num_large_shares/2)
        I3=I2
        numIntervals=3
        X = collect(partitions(V,numIntervals))

        possInd = Array{Int64}(undef,0)
        for i=1:length(X)
            A=X[i]
            if (A[1]*alpha+A[2]*xbuddy+A[3]*1//2 < m/s) && (A[1]*ybuddy+A[2]*1/2+A[3]*x> m/s)
                append!(possInd, i)
            end
        end
        S=sᵥ
        matrix=(X[possInd[1]])

        for i=2:length(possInd)
            Y=(X[possInd[i]])

            matrix=[matrix Y]
        end

        row,col=size(matrix)

        matrix=[matrix; transpose(ones(Int64,(col)))]

        println("A:")
        display(matrix)



        m=Model(with_optimizer(GLPK.Optimizer))
        @variable(m, x[i=1:length(possInd)],Int)

        b=[I1;I2;I3;S]
        println("b: ",b)
        @constraint(m,con,matrix*x .==b)
        optimize!(m)
        if(!has_values(m))
            println("No solution to this system on the Naturals: 𝛂 ≤ ",alpha)
        else
            potInt = Array{Ratinal}(undef,0)
            for targetInterval=1: numIntervals
                for i in possInd
                    if(X[i][targetInterval]!=0)

                    end
                end
                println()
            end
        end
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

        endpoints =Array{Rational, 2}(undef,0,0)
        endpoints =[y 1//2]
        endpoints = [endpoints;[1//2 ybuddy]; [xbuddy (1-alpha)]]
        sharesInIntervals = Array{Int64}(undef,0)
        append!(sharesInIntervals, [Int64(num_small_shares/2) Int64(num_small_shares/2) num_large_shares])
        numIntervals = length(sharesInIntervals)
        for i = 1:numIntervals
            print(endpoints[i,:])
            println("  ",sharesInIntervals[i])
        end
        X = perm(V-1, numIntervals)
        possInd = Array{Int64}(undef,0)
        for i=1:length(X)
            A=X[i]
            if (A[1]*y+A[2]*1//2+A[3]*xbuddy < m/s) && (A[1]*1//2+A[2]*ybuddy+A[3]*(1-alpha)> m/s)
                append!(possInd, i)
            end
        end
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
                    if lowerbound_temp>lowerbound
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
                    if upperbound_temp <upperbound
                        upperbound = upperbound_temp
                    end
                end #end if(X[i][j]!=0)
            end  #end  for i in possInd
        #    println()
        #    println("For I[",j,"]   lower = ",lowerbound," upper = ",upperbound)
        #    println()
            if(lowerbound>upperbound && lowerbound>endpoints[j,1] && upperbound < endpoints[j,2])
                endpointsArray=Array{Rational,1}(undef,0)
                row,col=size(endpoints)
                for i=1:row
                    for j=1:col
                        append!(endpointsArray, endpoints[i,j])
                    end
                end
                append!(endpointsArray, upperbound)
                append!(endpointsArray, lowerbound)
                endpoints=[endpoints; lowerbound upperbound; (1-lowerbound) (1-upperbound)]
        #        display(endpoints)
                endpoints=sort(collect(Iterators.flatten(endpoints)))
        #        println(endpoints)
                endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
                endpoints = transpose(endpoints)
                numIntervals=numIntervals+2
            #    endpoints=sort(endpoints, dims=1)
            #    endpoints=sort(endpoints, dims=2)
        #        display(endpoints)
                #endpoints = reshape(endpointsArray,2)
            end #end if(lowerbound>upperbound)
        end #end for j=1:numIntervals
        display(endpoints)
        X = perm(V-1, numIntervals)
        possInd = Array{Int64}(undef,0)
        for i=1:length(X)
            A=X[i]
            sum_1=0
            sum_2=0
            for j =1:numIntervals
                sum_1=sum_1+A[j]*endpoints[j,1]
                sum_2=sum_2+A[j]*endpoints[j,2]
            end

            if (sum_1 < m//s) && (sum_2 > m//s)
                println(sum_1,"   ",sum_2,"    m/s: ",m//s)
                append!(possInd, i)
            end
        end
        for i in possInd
            println(X[i,:])
        end


    end
end

GAP_proof(31,19,54//133)
