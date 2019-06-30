#include("src\\permutations.jl")
include("helper_functions.jl")
using JuMP
using GLPK
using Cbc
using LinearAlgebra
function linearSearchMID(m,s,array)
    for i=1:length(array)
        println(array[i])
    #    println("   ",MID_proof(m,s,array[i]))
#        (MID_proof(m,s,array[i]))
#        println("test")
        if VMID(m,s,array[i])==true
        #    MID_proof(m,s,array[i])
        #    GAP_proof(m,s,array[i])
            return array[i]
        end
    end
    return -1
end

function binarySearchMID(m,s,array)
    return binSearchHelpMID(m,s,array,1,length(array))
end
function binSearchHelpMID(m,s,array,front,back)
#    println("front: ",front,"  back: ",back)

    if front>back
        return 1
    end
    guessIndex =Int64(floor((front+back)/2))
#    println("m: ",m,"  s: ",s,"  alpha: ",array[guessIndex])
    if guessIndex == front
        if VMID(m,s,array[back])==true
            return array[back]
        else
            return 1
        end
    end

    if VMID(m,s,array[guessIndex]) == true
        return binSearchHelpMID(m,s,array,front,guessIndex)
    else
        return binSearchHelpMID(m,s,array,guessIndex, back)
    end
end

function VMID(m,s,alpha)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁

    x,y=FINDEND(m,s,alpha,V)

    xbuddy = 1-x
    ybuddy = 1-y

    if x > y
        return false
    end
    #endpoints = print_Intervals(m,s,alpha)
    if(V₋₁shares<Vshares)
         num_small_shares = V₋₁shares
        num_large_shares = Vshares - V₋₁shares
        S=sᵥ
        shares=Vshares
        VV=V #which is split
        I1 = num_small_shares
        I2 = Int64(num_large_shares/2)
        I3 = I2
        endpoints=Array{Rational,2}(undef,0,0)
        endpoints =  [alpha ybuddy]
        endpoints = [endpoints; [xbuddy 1//2]; [1//2 x]]
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
        endpoints =Array{Rational, 2}(undef,0,0)
        endpoints =[y 1//2]
        endpoints = [endpoints;[1//2 ybuddy]; [xbuddy (1-alpha)]]
        #num_split_shares=Int64(num_small_shares/2)

    end #******************************end else


    row, col= size(endpoints)

    numIntervals = row
    X = perm(VV, numIntervals)
    possInd = Array{Int64}(undef,0)


    X = perm(VV, numIntervals)
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


    A = transpose(poss_Dist)

    #Append a row of ones to A (num stud per distribtion adds to num VV students)
    row,col=size(A)
    Ones=(ones(Int64,col))'
    A=vcat(A, Ones)
    row = row+1


    model=Model(with_optimizer(Cbc.Optimizer, logLevel = 0))
    @variable(model, x[i=1:length(possInd)],Int)
    b=[I1,I2,I3,S]
    @constraint(model,con,A*x .==b)
    @constraint(model,con_1,x.>=0)

    optimize!(model)
    if (termination_status(model) == MOI.OPTIMAL)
        naturals = true
        for i = 1:length(value.(x))
            if (value.(x)[i]<0)
                naturals = false
            end
        end
        if naturals
            return false
        end
    end
        return true
end


function MID_proof(m,s,alpha)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁

    x,y=FINDEND(m,s,alpha,V)

    xbuddy = 1-x
    ybuddy = 1-y


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

        A = transpose(poss_Dist)

        #Append a row of ones to A (num stud per distribtion adds to num VV students)
        row,col=size(A)
        Ones=(ones(Int64,col))'
        A=vcat(A, Ones)
        row = row+1
        display(A)
        model=Model(with_optimizer(Cbc.Optimizer, logLevel=0))
        @variable(model, x[i=1:length(possInd)],Int)
        b=[I1,I2,I3,S]
        @constraint(model,con,A*x .==b)
        @constraint(model,con_1[i=1:length(possInd)],x[i] >=0)
        println("System of equations = ",b)
        display(A)


        optimize!(model)
        if (termination_status(model) == MOI.OPTIMAL)
            naturals = true
            for i = 1:length(value.(x))
                if (value.(x)[i]<0)
                    naturals = false
                end
            end
            if naturals
                println()
                println(value.(x))
                println("There is a solution on the Naturals")
            #println(value.(x))
                println("alpha ≥ ",alpha)
                println()
                return
            end
        end
            println("No solution on the Naturals")
            println("alpha ≤ ",alpha)
            return

end
function MID(m,s)
    #println("Test: m ",m,"  s ",s)
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
    array = filter( x -> x > 1/3, array)


    alpha = binarySearchMID(m,s,array)

    if alpha==-1
        return 1
    else
        return alpha
    end
end
