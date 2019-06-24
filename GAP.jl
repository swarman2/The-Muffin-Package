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
while true
        row, col= size(endpoints)
        numIntervals = row
        X = perm(VV, numIntervals)
        possInd = Array{Int64}(undef,0)
    #    for i=1:length(X)
    #        A=X[i]
    #        if (A[1]*y+A[2]*1//2+A[3]*xbuddy < m/s) && (A[1]*1//2+A[2]*ybuddy+A[3]*(1-alpha)> m/s)
    #            append!(possInd, i)
    #        end
    #    end
    println("Endpoints: ")
    display(endpoints)
        for i=1:length(X)
            A=X[i]
            sum_1=0
            sum_2=0
            for j =1:numIntervals
                sum_1=sum_1+A[j]*endpoints[j,1]
                sum_2=sum_2+A[j]*endpoints[j,2]
            end

    #            println("******************",sum_1,"   ",sum_2,"    m/s: ",m//s)

            if (sum_1 < m//s) && (sum_2 > m//s)
    #             println("############################",sum_1,"   ",sum_2,"    m/s: ",m//s)
            #    println(sum_1,"   ",sum_2,"    m/s: ",m//s)
                append!(possInd, i)
            end
        end
        println("Possible distributions of muffins: ")
        for i in possInd
            println(X[i,:])
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
        #    println()
    #        println("For I[",j,"]   lower = ",lowerbound," upper = ",upperbound)
        #    println()
            if(lowerbound>upperbound && lowerbound>endpoints[j,1] && upperbound < endpoints[j,2])
    #                println("For I[",j,"]   lower = ",lowerbound," upper = ",upperbound)
                _gap=true
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
        println("Endpoints: ")
        display(endpoints)
        println()
        if(_gap==false)
            println("No gaps found")
            println("alpha could be greater than ",alpha)
            return
        end
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

        #       println("******************",sum_1,"   ",sum_2,"    m/s: ",m//s)

            if (sum_1 < m//s) && (sum_2 > m//s)
        #               println("############################",sum_1,"   ",sum_2,"    m/s: ",m//s)
            #    println(sum_1,"   ",sum_2,"    m/s: ",m//s)
                append!(possInd, i)
            end
        end
        mat_1=reshape([],0,2)
        if length(possInd)==0
            println("No possible muffin distributions")
            println("alpha ≤ ",alpha)
            return
        end
        println("Possible muffin distributions")
        for i in possInd
            if i == possInd[1]
                mat_1=X[i,:]
            else
                mat_1=[mat_1; X[i,:]]
            end
            #println(X[i,:])
        end


        mat_1 = transpose(reshape(hcat(mat_1...), (length(mat_1[1]), length(mat_1))))
        display(mat_1)
        mat_2=[(mat_1[:,1]-mat_1[:,4]);(mat_1[:,2]-mat_1[:,3])]
        mat_2=reshape(mat_2,Int64(length(mat_2)/2),2)
        mat_2=transpose(mat_2)
        row,col=size(mat_2)
        Sum=(mat_1[:,1]+mat_1[:,2]+mat_1[:,3]+mat_1[:,4])'
        mat_2=vcat(mat_2, Sum)
        Ones=(ones(Int64,col))'
        mat_2=vcat(mat_2, Ones)

        #display(mat_2)


        model=Model(with_optimizer(GLPK.Optimizer))
        @variable(model, x[i=1:length(possInd)],Int)
    #    @variable(m, x[i=1:3],Int)
    #    print("size of x: ",size(x))
    #    mat_2=[1 -2 -3; 1 2 1]
        println("System of equations (=[0,0;",shares,";",S,"])")
        display(mat_2)
        b=[0;0;shares;S]
        @constraint(model,con,mat_2*x .==b)
        optimize!(model)
        if(has_values(model))
        #    println("There is a solution on the Naturals")
        #    println("alpha could be greater than ",alpha)
        else
            println("No solution on the Naturals")
            println("alpha ≤ ",alpha)
            return
        end
end


#    end
end

#GAP_proof(31,19,54//133)
#GAP_proof(41,19,131//304)
#GAP_proof(59,22,167//374)
#GAP_proof(41,23,149//368)
#GAP_proof(54,25,151//350)
GAP_proof(67,25,223//500)
#GAP_proof(59,26,191//442)
