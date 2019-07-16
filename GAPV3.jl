include("helper_functions.jl")
function VGAP_V3(m,s,alpha)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁

    x,y=FINDEND(m,s,alpha,V)

    xbuddy = 1-x
    ybuddy = 1-y
    if x>y
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
        endpoints =  [alpha x]
        #endpoints = [endpoints; [xbuddy x]]
        I1 = num_small_shares
        I2 = Int64(num_large_shares/2)
        I3 = I2
        println("TEST")

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
        endpoints =Array{Rational, 2}(undef,0,0)
        endpoints =[y 1-alpha]
        #endpoints = [endpoints; [xbuddy (1-alpha)]]
    #    sharesInIntervals = Array{Int64}(undef,0)
    #    append!(sharesInIntervals, [Int64(num_small_shares/2) Int64(num_small_shares/2) num_large_shares])

    #    for i = 1:numIntervals
    #        print(endpoints[i,:])
    #        println("  ",sharesInIntervals[i])
    #    end

    end #******************************end else
    endpoints = [alpha x]
    endpoints = [endpoints; [y (1-alpha)]]

    lower = x
    upper = y
    match = true
    while match == true
        match = false
        #buddy
        #println("[",1-upper,"  ",1-lower,"] by buddying [",lower,"  ",upper,"]")


        #match
        if 1-upper > y
            endpoints=[endpoints; [(1-upper) (1-lower)]]
        #    println("[",m//s-(1-lower),"  ",(m//s)-(1-upper),"] by matching [",1-upper,"  ",1-lower,"]")
            lower = m//s-(1-lower)
            upper = (m//s)-(1-upper)
            endpoints=[endpoints; [lower upper]]
            match = true
        end
    end

    row,col=size(endpoints)
    for i=1:row
        if endpoints[i,1]+endpoints[i,2] == (m//s)
            endpoints = [endpoints; [(m//2s) (m//2s)]]
        end
        if endpoints[i,1]+endpoints[i,2] == 1
            endpoints = [endpoints; [1//2 1//2]] #split it at 1/2
        end

    end

    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    #display(endpoints)

    row, col= size(endpoints)
    numIntervals = row
    X = perm(VV, numIntervals)
    possInd = Array{Int64}(undef,0)

    #Find the possible students
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
    if (length(possInd)==0)
        return true
    end

    _gap=false
    display(endpoints)
    while(_gap==true)

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
                lower = lowerbound
                upper = upperbound
                match = true
                while match == true
                    match = false
                    #buddy
                    #println("[",1-upper,"  ",1-lower,"] by buddying [",lower,"  ",upper,"]")


                    #match
                    if 1-upper > y
                        endpoints=[endpoints; [(1-upper) (1-lower)]]
                    #    println("[",m//s-(1-lower),"  ",(m//s)-(1-upper),"] by matching [",1-upper,"  ",1-lower,"]")
                        lower = m//s-(1-lower)
                        upper = (m//s)-(1-upper)
                        endpoints=[endpoints; [lower upper]]
                        match = true
                    end
                end
                endpoints = [endpoints; [1//2 1//2]] #split it at 1/2
                row,col=size(endpoints)
                for i=1:row
                    if endpoints[i,1]+endpoints[i,2] == (m//2s)
                        endpoints = [endpoints; [(m//2s) (m//2s)]]
                    end
                end


                endpoints=sort(collect(Iterators.flatten(endpoints)))
                endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
                endpoints = transpose(endpoints)
                numIntervals=numIntervals+2
            end #end if(lowerbound>upperbound)
        end #end for j=1:numInterval

    end#end while
    if(V₋₁shares<Vshares)

        println("ENDPOINTS AFTER FINDING GAPS")
        display(endpoints)
        numIntervals,col=size(endpoints)
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
                if ((endpoints[i,1]+endpoints[j,2])==1 && endpoints[i,2]-endpoints[i,1]==endpoints[j,2]-endpoints[j,1]) || ((endpoints[i,1]+endpoints[j,2])==m//s && endpoints[i,2]-endpoints[i,1]==endpoints[j,2]-endpoints[j,1])
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

        model=Model(with_optimizer(Cbc.Optimizer,logLevel=0))
        @variable(model, x[i=1:length(possInd)],Int)
        row,col=size(symmIntervals)

        #add zeros for each row of I₁ - I₂
        b=zeros(Int64,row)

        b=[b;num_split_shares; num_split_shares;S]
        @constraint(model,con,A*x .==b)
        @constraint(model,con_1,x.>=0)
        optimize!(model)
        if (termination_status(model)!=MOI.OPTIMAL)
            return true
        else
            return false
        end


end

VGAP_V3(55,51,139//408)
