#include("src\\permutations.jl")
include("helper_functions.jl")
using JuMP
using Cbc
using LinearAlgebra
function binarySearchMID(m,s,array)
    return binSearchHelpMID(m,s,array,1,length(array))
end
function binSearchHelpMID(m,s,array,front,back)
    #println("front: ",front,"  back: ",back)
    if front>back
        return 1
    end
    guessIndex =Int64(floor((front+back)/2))
    #println("m: ",m,"  s: ",s,"  alpha: ",array[guessIndex])
    if guessIndex == front
        if VMID(m,s,array[back],false)==true
            return array[back]
        else
            return 1
        end
    end
    if VMID(m,s,array[guessIndex],false) == true
        return binSearchHelpMID(m,s,array,front,guessIndex)
    else
        return binSearchHelpMID(m,s,array,guessIndex, back)
    end
end

function VMID(m,s,alpha,proof = false)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    x,y=FINDEND(m,s,alpha,V)
    xbuddy = 1-x
    ybuddy = 1-y
    if x > y
        if proof
            println("intervals not disjoint alpha > ", alpha)
        end
        return false
    end
    if proof
        endpoints = print_Intervals(m,s,alpha)
    else
        if(V₋₁shares<Vshares)
            endpoints =  [alpha ybuddy]
            endpoints = [endpoints; [xbuddy 1//2]; [1//2 x]]
        else
            endpoints =[y 1//2]
            endpoints = [endpoints;[1//2 ybuddy]; [xbuddy (1-alpha)]]
        end
    end
    if(V₋₁shares<Vshares)
         num_small_shares = V₋₁shares
        num_large_shares = Vshares - V₋₁shares
        S=sᵥ
        shares=Vshares
        VV=V #which is split
        I1 = num_small_shares
        I2 = Int64(num_large_shares/2)
        I3 = I2
        #num_split_shares=Int64(num_large_shares/2)

    else
        num_small_shares = V₋₁shares-Vshares
        num_large_shares =   Vshares
        S=sᵥ₋₁
        shares=V₋₁shares
        VV=V-1 #which is split
        I1=Int64(num_small_shares/2)
        I2=I1
        I3=num_large_shares
        #num_split_shares=Int64(num_small_shares/2)

    end #******************************end else

    row, col= size(endpoints)

    numIntervals = row
    X = perm(VV, numIntervals)
    possInd = Array{Int64}(undef,0)


        X = perm(VV, numIntervals)
        possInd = Array{Int64}(undef,0)

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
            if proof
                println("No possible muffin distributions")
                println("alpha ≤ ",alpha)
            end
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
        if proof
            println("Possible muffin distributions")
            row,col=size(poss_Dist)
            for i=1:row
                PRINT(poss_Dist[i,:])
                println()
            end
        end
        A = transpose(poss_Dist)

        #Append a row of ones to A (num stud per distribtion adds to num VV students)
        row,col=size(A)
        Ones=(ones(Int64,col))'
        A=vcat(A, Ones)
        row = row+1
        if proof
            display(A)
        end
        model=Model(with_optimizer(Cbc.Optimizer, logLevel=0))
        @variable(model, x[i=1:length(possInd)],Int)
        b=[I1,I2,I3,S]
        @constraint(model,con,A*x .==b)
        @constraint(model,con_1[i=1:length(possInd)],x[i] >=0)
        if proof
            println("System of equations = ",b)
            display(A)
        end


        optimize!(model)
        if (termination_status(model) == MOI.OPTIMAL)
            if proof
                println()
                println(value.(x))
                println("There is a solution on the Naturals")
            #println(value.(x))
                println("alpha > ",alpha)
                println()
            end
            return false
        end
        if proof
            println("No solution on the Naturals")
            println("alpha ≤ ",alpha)
        end
        return true

end
function MID(m,s,min_al=1//2)
    #println("Test: m ",m,"  s ",s)
    if m%s==0
        return 1
    end
    array=Array{Rational}(undef,0)
    alph = 1//3
    num=1
    denom=3
    while denom<=m*s
        while alph<min_al
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
    array = filter( x -> x > 1/3, array)
    #for i = 1:length(array)
    #    if array[i]==9//26
        #    println(i)
    #    end
    #end

    alpha = binarySearchMID(m,s,array)

    if alpha==-1
        return 1
    else
        return alpha
    end
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
#println(MID(33,26))
#println(VMID(33,26,157//377,true))
#MID(71,44)
