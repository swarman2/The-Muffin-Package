#for information on this method see chapter 10 of "The Mathematics of Muffins"
include("helper_functions.jl")
using JuMP
using Cbc

#VMID takes m,s,alpha and if proof is true prints a proof
#ret_endpts is used in PKG to pass the MID endpoints to PROC (speeds up Mulitsets)
function VMID(m,s,alpha,proof = false, ret_endpts =false)
    #setup with some information
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    x,y=FINDEND(m,s,alpha,V)
    xbuddy = 1-x
    ybuddy = 1-y
    denom=denominator(alpha)
    denom=lcm(s,denom)

    #if the intervals are not disjoint return false
    if x > y
        if proof
            println("intervals not disjoint alpha > ", alpha)
        end
        if ret_endpts
            return false, 0
        end
        return false
    end

    #detrmine endpoints (either with print intervals function or manually)
    if proof
        #prints out intervals
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
    #Some variables we use later
    #VV = V if the V shares are split, and V-1 if the V-1 shares are split
    #I1, I2, and I3 are the amounts of shares in the first second and third intervals
    # S = number of students in either V or V-1 (whichever is split)
    if(V₋₁shares<Vshares)
        S=sᵥ
        VV=V #which is split
        I1 = V₋₁shares
        I2 = Int64((Vshares - V₋₁shares)/2)
        I3 = I2
    else
        S=sᵥ₋₁
        shares=V₋₁shares
        VV=V-1 #which is split
        I1=Int64((V₋₁shares-Vshares)/2)
        I2=I1
        I3=Vshares
    end #******************************end else

    row, col= size(endpoints)
    numIntervals = row

    #find possible permutations
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

    #if there are no possible distrubtions return false
    if length(possInd)==0
        if proof
            println("No possible muffin distributions")
            println("alpha ≤ ",alpha)
        end
        if ret_endpts
            endpoints = collect(alpha*denom:1:(1-alpha)*denom)
            array = [xbuddy, ybuddy, 1-xbuddy, y]
            sort!(array)
            array=array*denom
            endpoints = filter(x -> x<=array[1]|| x>=array[2], endpoints)
            endpoints = filter(x -> x<=array[3]|| x>=array[4], endpoints)
            sort!(endpoints)
            return true, endpoints
        end
        return true
    end

    #make a matrix (poss_Dist) where each row is a
    #way to distribute the pieces to the students
    poss_Dist=reshape([],0,2) #possible distributions of muffins
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

    #A will be the matrix in our final matrix equation Ax=b
    #the columns of A are the distrubtions of students
    A = transpose(poss_Dist)

    #Append a row of ones to A (num stud per distribtion adds to num VV students)
    row,col=size(A)
    Ones=(ones(Int64,col))'
    A=vcat(A, Ones)
    row = row+1

    #b = the size of each interval and S (number of students)
    b=[I1,I2,I3,S]


    #set up a model and solve for the vector x
    #x represents the the number of students who get each type of muffin piece distrib.
    model=Model(with_optimizer(Cbc.Optimizer, logLevel=0))
    @variable(model, x[i=1:length(possInd)],Int)

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
            println("alpha > ",alpha)
            println()
        end
        if ret_endpts
            endpoints = collect(alpha*denom:1:(1-alpha)*denom)
            array = [xbuddy, ybuddy, 1-xbuddy, y]
            sort!(array)
            array=array*denom
            endpoints = filter(x -> x<=array[1]|| x>=array[2], endpoints)
            endpoints = filter(x -> x<=array[3]|| x>=array[4], endpoints)
            sort!(endpoints)
            return false, endpoints
        end
        return false
    end
    if proof
        println("No solution on the Naturals")
        println("alpha ≤ ",numerator(alpha),"/",denominator(alpha))
    end
    if ret_endpts
        endpoints = collect(alpha*denom:1:(1-alpha)*denom)
        array = [xbuddy, ybuddy, 1-xbuddy, y]
        sort!(array)
        array=array*denom
        endpoints = filter(x -> x<=array[1]|| x>=array[2], endpoints)
        endpoints = filter(x -> x<=array[3]|| x>=array[4], endpoints)
        sort!(endpoints)

        return true, endpoints
    end
    return true
end
function MID(m,s,min_al=1//2, ret_endpts = false)
    if m%s==0
        return 1,1
    end
    array=Array{Rational}(undef,0)
    #Find all numbers between 1/3 and 1/2 with denominator a multiple of s < m*s
    alph = 1//3
    num=1
    denom=3
    while denom<=m*s
        while alph<=min_al
            append!(array,alph)
            num=num+1
            alph = num//denom
        end
        num=Int64(ceil(denom/3))
        denom = denom+1
        while (denom%s!=0)
            denom=denom+1
        end
        alph = num//denom
    end
    sort!(array)
    unique!(array)
    array = filter( x -> x > 1/3, array)
    alpha, endpoints = binarySearchMID(m,s,array)
    if alpha == 1
        return 1,1
    end
    if ret_endpts
        endpoints = filter(x-> x!= 1//2, endpoints)
        return alpha, endpoints
    else
        return alpha
    end
end

#Helper function used to print arrays in form used by book
#(1 0 1) might go to (5/12 7/12)
function PRINT(array)
    print("( ")
    for i=1:length(array)
        for j=1:array[i]
            print(i," ")
        end
    end
    print(")")
end
#binary search using VMID
function binarySearchMID(m,s,array)
    return binSearchHelpMID(m,s,array,1,length(array))
end
function binSearchHelpMID(m,s,array,front,back)
    if front>back
        return 1
    end
    guessIndex =Int64(floor((front+back)/2))
    if guessIndex == front
        valid, endpoints = VMID(m,s,array[back],false, true)
        if valid==true
            return array[back], endpoints
        else
            return 1,1
        end
    end
    if VMID(m,s,array[guessIndex],false) == true
        return binSearchHelpMID(m,s,array,front,guessIndex)
    else
        return binSearchHelpMID(m,s,array,guessIndex, back)
    end
end
