include("helper_functions.jl")
include("GAP.jl")
using JuMP
using GLPK

function VTRAIN(m,s,alpha, proof = 0)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    x,y=FINDEND(m,s,alpha,V)
    if x>y
    #    print("                                                         intervals not disjoint")
        return false
    end
    ybuddy=1-y
    xbuddy=1-x
    denom=denominator(alpha)
    denom=lcm(s,denom)*2
    #println(denom)
    if V₋₁shares > Vshares
        #println("test")
        #endpoints   = [y*denom (1-alpha)*denom]
        #display(endpoints)
        endpoints = GAP_int(m,s,alpha)
        #display(endpoints)
        #println(alpha,"  ",x,"  ",y,"  ",1-alpha)
        #println(alpha)
        #display(endpoints)
        if endpoints == false
            if proof >0
        #        println("no endpoints")
            end
        #    print("no endpoints")
            return false
        end
        row,col = size(endpoints)
        #display(endpoints/2)
    #    if row <3
    #        return false
    #    end

        kL1 = endpoints[1,1]//denom
        kR1 = endpoints[1,2]//denom
        kL2 = endpoints[2,1]// denom
        kR2 = endpoints[row-1,2]//denom
        kL3 = endpoints[row,1]// denom
        kR3 = endpoints[row,2]//denom


    #    if xbuddy != endpoints[row,1]//denom
    #        return false
    #    end
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
        #if B == 1//(Z+1-alpha)
        #    print("                        used B = 1/(Z+1-alpha)")
        #    B = 2//(kL2-alpha)
        #    print("used B = 2/(kL2-alpha)")
        #else
        #    print("used B = 2/(kL2-alpha)")
        #end
        Bi = kL1
        if false
            println(B)
            println(Bi)
            Betai = Int64(630*Bi)
            while Bi<kR2
                if !(1 <= floor(B-B*alpha -B*Bi)+B*Z+B*Bi)
                    println("beta_i = ",Betai,"/630    :  1 > floor(B-B*alpha -B*Bi)+B*Z+B*Bi = ",floor(B-B*alpha -B*Bi)," + ",B*Z+B*Bi," = ",floor(B-B*alpha -B*Bi)+B*Z+B*Bi)
                end
                Bi = Bi + 6/denom
                Betai = Betai+3
            end
        end
        while floor(B*(1-alpha-kL1)) >= floor(B*(1-2*alpha))
            model = Model(with_optimizer(GLPK.Optimizer))
            @variable(model,Bi)
            @constraint(model, Bi*(1-alpha-kL1) in MOI.Integer())
            @constraint(model, Bi*(1-alpha-kL1) >= 0)
            @objective(model, Min, Bi)
            optimize!(model)
            #println("test")

            B = value(Bi)
        end
        if proof >0
            println(denom)
            display(endpoints/2)
            println("kL1: ",kL1*denom/2,"  kL2: ",kL2*denom/2,"  kL3: ",kL3*denom/2)
            println("Z: ",Z," Y: ",Y," a = ",a, " SS = ",SS)

            println("Let B = max(2//(kL2-alpha), 1//(Z+1-alpha))")
            println("B = ", B)
        end
        #B = 2//(kL2-alpha)
        A = B-B*alpha -1
        #model=Model(with_optimizer(GLPK.Optimizer))
        #@variable(model, A)
        #@variable(model, B)
        #@constraint(model, A, A== B-B*alpha -1)
        #@constraint(model, B >= x/(Z+1-alpha))
        #@constraint(model, B >= 2/(kL2 - alpha))
        #@objective(model, Min, 1+(A-B*kL1))
        if false
            Bi=kL1
            println("kL1 = ",float(kL1))
            model=Model(with_optimizer(GLPK.Optimizer))
            @variable(model,B_i)
            @constraint(model, kL1<=B_i<=kR2)
            @variable(model, eqn)
            @constraint(model, eqn ==  (B-B*alpha-B*B_i))
            @constraint(model,eqn in MOI.Integer())
            @objective(model, Min, B_i)
            optimize!(model)
            #println(1," <= ",floor(B-B*alpha-B*Bi)+B*Z+B*Bi)
            #println(B-B*alpha-B*Bi)
            #println("Initial B: ",B)
            SUM = 0
            initialB=B
            #guess = (1-floor(B-B*alpha-B*Bi)+B*Z+B*Bi)//(Z//denom + Bi//denom)
            #println("Guess: ",B+guess//denom)
            #while !(1<=floor(B-B*alpha-B*Bi)+B*Z+B*Bi)

                #println(Bi)
                while termination_status(model)==MOI.OPTIMAL
                    #    if proof >0
                    #    println("Error no B_i, increasing B by 1/denom")
                    #    println("B = ",B+1//denom)

                    #    end
                    #    if 1-(floor(B-B*alpha-B*kR2)+B*Z+B*kR2)>=1
                            #println(1-(floor(B-B*alpha-B*kR2)+B*Z+B*kR2))
                    #        println("add 1")
                    #        B=B+1
                    #        SUM=SUM+1
                    #    else
                        #    println("add 1/denom")
                            #B= B+1//denom

                            #SUM=SUM+1//denom
                    #    end
                    #println(B)
                    if false
                        model=Model(with_optimizer(GLPK.Optimizer))
                        @variable(model,B_i)
                        @constraint(model, kL1<=B_i<=kR2)
                        @variable(model, eqn)
                        @constraint(model, eqn ==  (B-B*alpha-B*B_i))
                        @constraint(model,eqn in MOI.Integer())
                        @objective(model, Min, B_i)
                        optimize!(model)
                        #B_i =0
                        #return true #TODO this is a guess
                    end
                    Bi = value(B_i)
                    #Bi = round(Bi, digits = 4)
                    println("*******************************")
                    println("Bi = ",Bi)
                    println("B-B*alpha-B*B_i = ",B-B*alpha-B*Bi)
                    println(floor(B-B*alpha-B*Bi)," + ",B*Z+B*Bi," = ", floor(B-B*alpha-B*Bi)+B*Z+B*Bi)
                    println()

                    Bi = Bi+10e-4
                    #println(10e-4)
                    #Bi = round(Bi, digits = 4)
                    println("Bi + eps = ",Bi)
                    println("B-B*alpha-B*B_i = ",B-B*alpha-B*Bi)
                    if 1 > floor(B-B*alpha-B*Bi)+B*Z+B*Bi && Bi > kL1
                        println("floor(B-B*alpha-B*Bi) + B*Z+B*Bi = ",floor(B-B*alpha-B*Bi)," + ",B*Z+B*Bi," = ", floor(B-B*alpha-B*Bi)+B*Z+B*Bi)
                        println("error")
                    end
                    @constraint(model, B_i >= value(B_i)+10e-4)
                    optimize!(model)
                end
            #println("Bi = ",Bi,"  kL1 = ",kL1)
            #    println("final value of B = ",B)
            #println("diference = ",B-initialB)
            #println("had to add: ",SUM)
                A = B-B*alpha -1

                #    B = 0
            B0 =0
            #    A=0

            println(floor(B*(1-alpha-kL1))," < ",floor(B*(1-2*alpha)))
            while(floor(B*(1-alpha-kL1))>= floor(B*(1-2*alpha)))

            B=B+1//denom
            println(floor(B*(1-alpha-kL1))," < ",floor(B*(1-2*alpha)))
            break
            end
        end
        eps = 10e-10
        #println(eps)
        i=0
        #model = Model(with_optimizer(GLPK.Optimizer))
        #@variable(model,A)
        #@variable(model,B)
        #@constraint(model, A == B-B*alpha-1)
        #@constraint(model, B >= 2//(Z+1-alpha))
        #@constraint(model, B >= 2//(kL2 - alpha))
        #@objective(model, Min, 1+(A-B*kL1))
        #optimize!(model)
        #println(value(A),"   ",value(B))
        if false
            B0 = max(2//(Z+1-alpha),2//(kL2-alpha))
            #println(2//(Z+1-alpha),"  ",2//(kL2-alpha))
            A = B0-B0*alpha -1
            #    println("A_",i," = ",A,"  B_",i," = ",B0)
            i=i+1
            if floor(B0*(1-alpha-kL1)) < floor(B0*(1-2*alpha))
                B = B0
            #    break
            else
                model = Model(with_optimizer(GLPK.Optimizer))
                @variable(model, B1 )

                @constraint(model, B1 >= B0+ eps)
                @constraint(model, B1-2*alpha in MOI.Integer())
                @constraint(model, B1-2*alpha >=0)
                @objective(model, Min, B1)
                optimize!(model)
                if termination_status(model) == MOI.OPTIMAL
                    B0 = B1
                else
                    B0 = 0
            #        break
                end
            end
        end
        #if B != B0
        #    println("error")
        #    return
        #end
        if proof >0
            println("A = ",A,"  B = ",B)
            println(2*(1+floor(A-B*(kL1+eps)))," < ",a)
        end

        if 2(1+floor(A-B*(kL1+eps)))<a
            #println("Yes")
            return true
        else
        #    println("A = ",value(A)," B = ",value(B))
            #println("DK")
            #@constraint(model, B >= value(B)+.001)
            #optimize!(model)
        #    println(termination_status(model))
            return false
        end

    else
        #println("V case")
        return false
        endpoints = GAP_int(m,s,alpha)
        row,col = size(endpoints)
        if row <3
            return false
        end
        kR3 = endpoints[row,col]*denom
        kR2 = endpoints[row-1,col]*denom
        LS = -(2*V*V*s - 2*V*s - 4*V*m + 2*m)
        a = ceil(LS//(2*m - (V-1)*s - (.5)*LS))
        if a%2 !=0
            a = a+1
        end
        Z = m//s - 1 - (V-1)*alpha

        permV=perm(V, row)
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
        numLargePossible = 0
        mat_V=transpose(hcat(possV...))
        row_V, col_V = size(mat_V)
        for i = 1: row_V
            max = sum(mat_V[i,:]) - mat_V[i,1]
            if max >numLargePossible
                numLargePossible = max
            end
        end
        Y = numLargePossible -1
        #println("Z: ",Z," Y: ",Y," a = ",a, " SS = ",SS,"  kL1 = ",kL1, " kL2 = ", kL2)

        model=Model(with_optimizer(GLPK.Optimizer))
        @variable(model, A)
        @variable(model, B)
        @objective(model, Min, 1+(B*kR3 -A))

        @constraint(model,A + B*Z +1 <= 0)
        #@constraint(model,A + B >= 0)
        @constraint(model,(A-1 -B*alpha)*Y <= A+B*Z +1)
        @constraint(model, (B-(kR2)*B -1)>=A)

        optimize!(model)
        #println(termination_status(model))
        #println(value(A)," ",value(B))
        #println(value(A) - value(B)*kL1)
        if !(1 + floor(value(B) - value(A)*2)> floor(value(B)*kR3 - value(A)))
            #println("A and B failed")
            return false
        end
        eps = 1//26
        if 2(1+floor(value(B)*(kR3+eps)-value(A)))<a
            #println("Yes")
            return true
        else
            #println("DK")
            return false
        end
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
    denom=lcm(s,denom)*2
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
        #    endpoints, gap_found = findGaps(mat_V, endpoints,m,s)
        mat_V_=transpose(hcat(possV_...))
        endpoints, gap_found = findGaps(mat_V,mat_V_, endpoints,m,s,denom,0)
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
#VTRAIN(11,5,13//30)
#VTRAIN(31,27,103//297)
#println(TRAIN(67,21))
#VTRAIN(67,21,41//90,1)
#println(TRAIN(83,26))
#println(TRAIN(69,32))
#println(TRAIN(91,34))
#VTRAIN(67,21,41//90,1)
#TRAIN(67,21)

#printlnf(TRAIN(67,21))
#println(TRAIN(29,23))
#VTRAIN(69,32,689//1600)
#VTRAIN(83,26,77//169)
#println(TRAIN(83,26))
#println(TRAIN(69,32))
#VTRAIN(69,32,921//2144,1)

#VTRAIN(67,21,41//90,1)
if false
    pTrain(67,21)
#pTrain(94,25)
pTrain(83,26)
pTrain(69,32)
pTrain(91,34)
pTrain(110,41)
pTrain(71,44)
pTrain(101,47)
pTrain(85,52)
pTrain(95,59)
pTrain(97,60)
pTrain(101,62)
pTrain(103,63)
pTrain(107,66)
end
