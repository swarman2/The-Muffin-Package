include("helper_functions.jl")
include("GAP.jl")
using JuMP
using GLPK
#if proof == 3 prints to file
function VTRAIN(m,s,alpha, proof = 0)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    x,y=FINDEND(m,s,alpha,V)
    if proof == 3

        file =open(dirname(@__FILE__)*"/../DATA/"*string(m)*"-"*string(s)*".txt","a+")
        println(file)
        println(file,"*******************************************************")
        println(file,"*                 TRAIN PROOF                          *")
        println(file,"*******************************************************")
        println(file)
    end

    if x>y
    #    print("                                                         intervals not disjoint")
        if proof >=1
            println("Intervals not disjoint, alpha could be > ",alpha)
        end
        if proof == 3
            println(file,"Intervals not disjoint, alpha could be > ",alpha)
            close(file)
        end
        return false
    end
    #if x == alpha || y == (1-alpha)
    #    if proof>=1
    #        println("Weird intervals")
    #    end
    #    return false, 1
    #end
    ybuddy=1-y
    xbuddy=1-x
    denom=denominator(alpha)
    denom=lcm(s,denom)
    #println(denom)

    if V₋₁shares > Vshares
        endpoints = GAP_int(m,s,alpha)

        #this if never really happens
        if endpoints == false
            if proof >0
                println("no endpoints")
            end
            if proof == 3
                println(file,"no endpoints")
                close(file)
            end
            return false
        end
        if proof>=1
            print("The following numbers assumed to have denominator: ")
            println(denom)
            println("Intervals: ")
            display(convert_Int(endpoints))
        end
        if proof == 3
            print(file,"The following numbers assumed to have denominator: ")
            println(file,denom)
            println(file,"Intervals: ")
            row,col = size(endpoints)
            endpoints = convert_Int(endpoints)
            for row_end = 1:row
                for col_end = 1: col
                    print(file,endpoints[row_end,col_end],"  ")
                end
                println(file)
            end
        end
        row,col = size(endpoints)

        kL1 = endpoints[1,1]//denom
        kR1 = endpoints[1,2]//denom
        kL2 = endpoints[2,1]// denom
        kR2 = endpoints[row-1,2]//denom
        kL3 = endpoints[row,1]// denom
        kR3 = endpoints[row,2]//denom

        SS = 2*V*V*s - 2*V*s-4*V*m + 2*m #number of small shares
        a = ceil(SS/(V*s - 2m - (.5 * SS)))
        if a%2 != 0
            a=a+1
        end
        Z = m//s - 1 - (V-2)*(1-alpha)

        LastSmall = endpoints[row-1,2]
        permV=perm(V-1, row)
        possV = Vector{Vector{Int64}}()
        #Find the possible distribtions of muffins
        for i=1:length(permV)
            A=permV[i]
            sum_1=0
            sum_2=0
            equal = true
            #equal = false
            for j =1:row
                sum_1=sum_1+A[j]*endpoints[j,1]
                sum_2=sum_2+A[j]*endpoints[j,2]
                if A[j]!=0 && endpoints[j,1]!=endpoints[j,2]
                    equal=false
                end
            end
            if !equal && (sum_1 < (m//s)*denom) && (sum_2 > (m//s)*denom)
                append!(possV, permV[i,:])
            elseif equal && (sum_1 <= (m//s)*denom) && (sum_2 >= (m//s)*denom)
                append!(possV, permV[i,:])
            end
        end
        numSmallPossible = 0
        mat_V=transpose(hcat(possV...))
        #display(mat_V)
        if length(mat_V) == 0
        #    print("                                                         no student distributions")
            if proof >=1
                println("No types of students, alpha ≤ ",alpha)
            end
            if proof == 3
                println(file,"No types of students, alpha ≤ ",alpha)
                close(file)
            end

            return true
        end
        row_V, col_V = size(mat_V)
        for i = 1: row_V
            max = sum(mat_V[i,:]) - mat_V[i,col_V]
            if max >numSmallPossible
                numSmallPossible = max
            end
        end
        Y = numSmallPossible -1
        B = max(2//(kL2-alpha), 1//(Z+1-alpha))
        if proof >0

            #display(convert_Int(endpoints/2))
            println("kL1: ",Float64(kL1*denom),"  kL2: ",Float64(kL2*denom),"  kL3: ",Float64(kL3*denom))
            println("Z: ",Float64(Z*denom)," Y: ",Y," a = ",a, " SS = ",SS)
            println()
            println("Let B = max(2/(kL2-alpha), 1/(Z+1-alpha))")
            println("--note: A and B do NOT have denominator ",denom)
            println("B = ", numerator(B),"/",denominator(B))
            println("=> A = ", numerator(B-B*alpha -1),"/",denominator(B-B*alpha -1))
        end
        if proof == 3
            println(file,"kL1: ",Float64(kL1*denom),"  kL2: ",Float64(kL2*denom),"  kL3: ",Float64(kL3*denom))
            println(file,"Z: ",Float64(Z*denom)," Y: ",Y," a = ",a, " SS = ",SS)
            println(file,)
            println(file,"Let B = max(2/(kL2-alpha), 1/(Z+1-alpha))")
            println(file,"--note: A and B do NOT have denominator ",denom)
            println(file,"B = ", numerator(B),"/",denominator(B))
            println(file,"=> A = ", numerator(B-B*alpha -1),"/",denominator(B-B*alpha -1))
        end
        A = B-B*alpha -1
        eps = 10e-10
        i=0
        if proof >0
            println()
            println("--note: these numbers do NOT have denominator ",denom)
            println("2(1+floor(",numerator(A),"/",denominator(A)," - ",numerator(B),"/",denominator(B),"(",Float64(numerator(kL1)*(denom/denominator(kL1))),"/",denom," + ",eps,"))) < ",a)
            println(2*(1+floor(A-B*(kL1+eps)))," < ",a)
        end
        if proof == 3
            println(file)
            println(file,"--note: these numbers do NOT have denominator ",denom)
            println(file,"2(1+floor(",numerator(A),"/",denominator(A)," - ",numerator(B),"/",denominator(B),"(",Float64(numerator(kL1)*(denom/denominator(kL1))),"/",denom," + ",eps,"))) < ",a)
            println(file,2*(1+floor(A-B*(kL1+eps)))," < ",a)
            close(file)
        end

        if 2(1+floor(A-B*(kL1+eps)))<a
            return true
        else
            return false
        end

    else
        #println("V case")
        return false

    end

end
function GAP_int(m,s,alpha)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    x,y=FINDEND(m,s,alpha,V)
    ybuddy=1-y
    xbuddy=1-x

    denom=denominator(alpha)
    denom=lcm(s,denom)
    #print_Intervals(m,s,alpha)
    endpoints=[alpha*denom x*denom]
    endpoints=[endpoints;y*denom (1-alpha)*denom; xbuddy*denom ybuddy*denom]
    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    if x>=y
        return false
    end
    endpoints=[endpoints; denom//2 denom//2 ]
    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    endpoints = buddymatch(endpoints,V,y,m,s,denom)
    #display(endpoints)
    row,col=size(endpoints)
    if V==3
        endpoints = [endpoints; [(m*denom//2s) (m*denom//2s)];[(m*denom//2s) (m*denom//2s)]]
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
                    if endpoints[k,j] == denom-upper
                        buddyIn_1=true
                    end
                    if endpoints[k,j]==denom-lower
                        buddyIn_2=true
                    end
                end
            end

            if buddyIn_1 == false && buddyIn_2 == false && denom-upper >= endpoints[1,1] && denom-lower <= endpoints[row,col]
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
                if matchIn_1 == false && matchIn_2 == false && (m//s)*denom-upper >= endpoints[1,1] && (m//s)*denom-lower <= endpoints[row,col]
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
        endpoints = buddymatch(endpoints,V,y,m,s,denom,0)

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
            equal = true
            #equal = false
            for j =1:numVIntervals
                sum_1=sum_1+A[j]*endpoints[j,1]
                sum_2=sum_2+A[j]*endpoints[j,2]
                if A[j]!=0 && endpoints[j,1]!=endpoints[j,2]
                    equal=false
                end
            end
            if !equal && (sum_1 < (m//s)*denom) && (sum_2 > (m//s)*denom)
                append!(possV, permV[i,:])
            elseif equal && (sum_1 <= (m//s)*denom) && (sum_2 >= (m//s)*denom)
                append!(possV, permV[i,:])
            end
        end

        #Find the possible distribtions of muffins
        for i=1:length(permV_)
            A=permV_[i]
            sum_1=0
            sum_2=0
            equal = true
        #equal = false
            for j =numVIntervals+1:row
                sum_1=sum_1+A[j-numVIntervals]*endpoints[j,1]
                sum_2=sum_2+A[j-numVIntervals]*endpoints[j,2]
                if A[j-numVIntervals]!=0 && endpoints[j,1]!=endpoints[j,2]
                    equal=false
                end
            end
            if !equal && (sum_1 < (m//s)*denom) && (sum_2 > (m//s)*denom)
                append!(possV_, permV_[i,:])
            elseif equal && (sum_1 <= (m//s)*denom) && (sum_2 >= (m//s)*denom)
                append!(possV_, permV_[i,:])
            end
        end


        mat_V=transpose(hcat(possV...))
        mat_V_=transpose(hcat(possV_...))
        numVDistributions,numVIntervals=size(mat_V)
        numV_Distributions, numV_Intervals=size(mat_V_)

        gap_found1 = false
        gap_found2 =false

        endpoints, gap_found1 = findGaps(mat_V,1,numVIntervals, endpoints,m,s,denom,0)
        if !gap_found1
            endpoints, gap_found2 = findGaps(mat_V_,numVIntervals+1,numVIntervals + numV_Intervals, endpoints,m,s,denom,0)
        end
        if gap_found1||gap_found2
            gap_found = true
        else
            gap_found=false
        end
    end #end finding gaps
    numVIntervals = 0
    for i =1: row
        if endpoints[i,2]<=x*denom
            numVIntervals = numVIntervals+1
        end
    end
    row,col=size(endpoints)
    if V₋₁shares > Vshares
        return endpoints[numVIntervals+1:row,:]
    else
        return endpoints[1:numVIntervals,:]
    end
end

function TRAIN(m,s,min_al = 1//2)
    if m%s==0
        return 1
    end
    array=Array{Rational}(undef,0)
    alph = 1//3
    append!(array,alph)
    num=1
    denom=3
    while denom<=m*s
        while alph<=min_al
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
    for i=1:length(array)
        if array[i]==41//90
    #        println("*************  ",i)
        end
    end

    #println(array)
    alpha = binarySearchTRAIN(m,s,array)
    #    println(alpha)
    #    println(alpha)
    if alpha==-1
        return 1
    else
        return alpha
    end
end
#binary search, used in GAP uses PROC and VGAP to find the
#alpha is a sorted array

function binarySearchTRAIN(m,s,array)
    return binSearchHelpTRAIN(m,s,array,1,length(array))
end
function binSearchHelpTRAIN(m,s,array,front,back)
    #println("m: ",m,"  s: ",s)
#    println()
    if front>back
        return 1
    end
    guessIndex =Int64(floor((front+back)/2))
    #println("front: ",front,"  back: ",back,"  alpha: ",array[guessIndex])

    if guessIndex == front
#        print("    VTRAIN(",m,",",s,",",array[back],")     ")
        if VTRAIN(m,s,array[back])==true
            return array[back]
        else
            return 1
        end
    end
#print("    VTRAIN(",m,",",s,",",array[guessIndex],")     ")
    if VTRAIN(m,s,array[guessIndex]) == true
        return binSearchHelpTRAIN(m,s,array,front,guessIndex)
    else
        return binSearchHelpTRAIN(m,s,array,guessIndex, back)
    end
end
function pTrain(m,s)
    alpha = TRAIN(m,s)
    println()
    println("f(",m,",",s,") = ",alpha)
end
