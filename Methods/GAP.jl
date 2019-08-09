include("helper_functions.jl")
using JuMP
using Cbc
#input
#if proof = 1 prints proof, if proof =2 prints proof with matrix
#if ret_endpts = true it returns the endpoints
#m,s,alpha (muffins, students, alpha)

#output
# bool , endpts
# true if GAP can show upperbound false if it can't
# endpts = matrix or 0 if no endpts
function VGAP(m,s,alpha, proof=0, ret_endpts = false)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    x,y=FINDEND(m,s,alpha,V)
    ybuddy=1-y
    xbuddy=1-x
    #if x == alpha || y == (1-alpha)
    #    if proof >=1
    #        println("Weird intervals")
    #    end
    #    return true, 0
    #end
    denom=denominator(alpha)
    denom=lcm(s,denom)
    if proof >= 3
        file = open(dirname(@__FILE__)*"/../DATA/"*string(m)*"-"*string(s)*".txt","a+")
        println(file)
        println(file,"*******************************************************")
        println(file,"*                 GAP PROOF                          *")
        println(file,"*******************************************************")
        println(file)

    else
        file =0
    end
    if proof >=1
        println()
        println("GAP PROOF OF f(",m,", ",s,") = ",numerator(alpha),"/",denominator(alpha))
        println("There are ",V,"-students and ",V-1,"-students")
        println(sᵥ," ",V,"-shares  ",sᵥ₋₁," ",V-1,"-shares")
        println("All numbers assumed to have denominator ",denom)
    end
    if proof >=3
        println(file)
        println(file,"GAP PROOF OF f(",m,", ",s,") = ",numerator(alpha),"/",denominator(alpha))
        println(file,"There are ",V,"-students and ",V-1,"-students")
        println(file,sᵥ," ",V,"-shares  ",sᵥ₋₁," ",V-1,"-shares")
        println(file,"All numbers assumed to have denominator ",denom)
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
        #endpoints = [endpoints; [xbuddy x]]
        I1 = num_small_shares
        I2 = Int64(num_large_shares/2)
        I3 = I2
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

    end #******************************end else
    #print_Intervals(m,s,alpha)
    endpoints=[alpha*denom x*denom]
    endpoints=[endpoints;y*denom (1-alpha)*denom; xbuddy*denom ybuddy*denom]
    Endpoints=[alpha*denom x*denom]
    Endpoints=[Endpoints;y*denom (1-alpha)*denom; xbuddy*denom ybuddy*denom]
    #you will see these next three lines throughout
    #the sort the matrix smallest to largest (left to right top to bottom)
    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    Endpoints=sort(collect(Iterators.flatten(Endpoints)))
    Endpoints=reshape(Endpoints,(2,Int64(length(Endpoints)/2)))
    Endpoints = transpose(Endpoints)
    #######
    if proof >=1
        println("Intervals before Buddy-Match" )
        endpoints=convert_Int(endpoints)
        display(endpoints)
    end
    if proof >=3
        println(file,"Intervals before Buddy-Match" )
        endpoints=convert_Int(endpoints)
        row, col = size(endpoints)
        for i = 1:row
            for j=1:col
                print(file,endpoints[i,j]," ")
            end
            println(file)
        end
    end
    if x>y
        if proof >=1
            println("INVALID INTERVALS")
        end
        if proof >=3
            println(file)
            println(file,"INVALID INTERVALS")
            close(file)
        end

        if ret_endpts
            B = collect(alpha*denom:1:(1-alpha)*denom)
            #display(B)
            #display(endpoints)
            row_e, col_e = size(endpoints)
            while row_e -1 > 0
                B=filter(x-> x<=endpoints[row_e-1,2]|| x>=endpoints[row_e,1],B)
                row_e = row_e-1
            end
            #display(B)
            B= filter(x-> x!= 1//2, B)
            return false, B
        end
        return false
    end
    if ! (denom//2 in endpoints)
        endpoints=[endpoints; denom//2 denom//2 ]
        if proof >= 1
            print("\nSPLIT AT 1/2 (",(denom/2),") ")
        end
        if proof >=3
            print(file,"\nSPLIT AT 1/2 (",(denom/2),") ")
        end
    end

    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    # if we are returning endpoints we need to find them seperatly (so we don't include 1/2 or m/2s)
    if ret_endpts
        Endpoints = endpoints
        gap_found=true
        while(gap_found == true)
            #buddy match it
            Endpoints = buddymatch(Endpoints,V,y,m,s,denom,0,0)
            #Endpoints = buddymatch(Endpoints,V,y,m,s,denom,0)
            row, col= size(Endpoints)
            numIntervals = row
            numVIntervals = 0
            for i =1: row
                if Endpoints[i,2]<=x*denom
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
                for j =1:numVIntervals
                    sum_1=sum_1+A[j]*Endpoints[j,1]
                    sum_2=sum_2+A[j]*Endpoints[j,2]
                end
                if (sum_1 <= (m//s)*denom) && (sum_2 >= (m//s)*denom)
                    append!(possV, permV[i,:])
                end
            end

            #Find the possible distribtions of muffins
            for i=1:length(permV_)
                A=permV_[i]
                sum_1=0
                sum_2=0
                for j =numVIntervals+1:row
                    sum_1=sum_1+A[j-numVIntervals]*Endpoints[j,1]
                    sum_2=sum_2+A[j-numVIntervals]*Endpoints[j,2]
                end
                if (sum_1 <= (m//s)*denom) && (sum_2 >= (m//s)*denom)
                    append!(possV_, permV_[i,:])
                end
            end

            #look for more gaps using student distribtions
            mat_V=transpose(hcat(possV...))
            mat_V_=transpose(hcat(possV_...))
            numVDistributions,numVIntervals=size(mat_V)
            numV_Distributions, numV_Intervals=size(mat_V_)

            gap_found1 = false
            gap_found2 =false

            Endpoints, gap_found1 = findGaps(mat_V,1,numVIntervals, Endpoints,m,s,denom,0)
            if !gap_found1
                Endpoints, gap_found2 = findGaps(mat_V_,numVIntervals+1,numVIntervals + numV_Intervals, Endpoints,m,s,denom,0)
            end
            if gap_found1||gap_found2
                gap_found = true
            else
                gap_found=false
            end
            #Endpoints, g = findGaps(mat_V, mat_V_, Endpoints,m,s,denom,0,false)
        end #end finding gaps
    end

    if proof >=1
        if V == 3 && !((denom*m/2s) in endpoints)
            print("  Split m/2s (",(denom*m/(2s)),")")
        end
        println()
    end
    if proof >=3
        if V == 3 && !((denom*m/2s) in endpoints)
            print(file,"  Split m/2s (",(denom*m/(2s)),")")
        end
        println(file)
    end
    endpoints = buddymatch(endpoints,V,y,m,s,denom,proof,file)
    row,col=size(endpoints)

    #if V=3 split at m/2s and buddy match in a special way
    #lets say we split at x, we make x its own interval
    #when we buddy match we have to make m/s - x its own interval
    if V==3
        if  !((denom*m/2s) in endpoints)
            endpoints = [endpoints; [(m*denom//2s) (m*denom//2s)];[(m*denom//2s) (m*denom//2s)]]
            endpoints=sort(collect(Iterators.flatten(endpoints)))
            endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
            endpoints = transpose(endpoints)
        end
        if proof >=1
        #    println("\nSPLIT AT m//2s (",Int64(denom*m//2s),")")
            endpoints = convert_Int(endpoints)
            display(endpoints)
        end
        if proof >=3
            endpoints = convert_Int(endpoints)
            row, col = size(endpoints)
            for i = 1:row
                for j=1:col
                    print(file,endpoints[i,j]," ")
                end
                println(file)
            end
        end
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
                if proof >=1
                    println("[",float(denom-upper),"  ",float(denom-lower),"] by buddying [",float(lower),"  ",float(upper),"]")
                end
                if proof >=3
                    println(file,"[",float(denom-upper),"  ",float(denom-lower),"] by buddying [",float(lower),"  ",float(upper),"]")
                end
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
                    if proof >=1
                        println("[",Float64((m//s)*denom-upper),"  ",Float64((m//s)*denom-lower),"] by matching [",Float64(lower),"  ",Float64(upper),"]")
                    end
                    if proof >=3
                        println(file,"[",Float64((m//s)*denom-upper),"  ",Float64((m//s)*denom-lower),"] by matching [",Float64(lower),"  ",Float64(upper),"]")
                    end
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

    # While iterate over the endpoints looking for more gaps until no more gaps found
    gap_found=true
    while(gap_found == true)
        #buddy match it
        endpoints = buddymatch(endpoints,V,y,m,s,denom,proof,file)


        #find out the distribtion of intervals (V or V-1)
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
        if proof >=1
            endpoints = convert_Int(endpoints)
            println("\n",V," intervals")
            for j=1:numVIntervals
                if endpoints[j,1]==endpoints[j,2]
                    println("I_",j,": [",endpoints[j,1],", ",endpoints[j,2],"]")
                else
                    println("I_",j,": (",endpoints[j,1],", ",endpoints[j,2],")")
                end
            end
            println("\n",V-1," intervals")
            for j=numVIntervals+1: row
                if endpoints[j,1]==endpoints[j,2]
                    println("J_",j-numVIntervals,": [",endpoints[j,1],", ",endpoints[j,2],"]")
                else
                    println("J_",j-numVIntervals,": (",endpoints[j,1],", ",endpoints[j,2],")")
                end
            end
        end
        if proof >=3
            endpoints = convert_Int(endpoints)
            println(file,"\n",V," intervals")
            for j=1:numVIntervals
                if endpoints[j,1]==endpoints[j,2]
                    println(file,"I_",j,": [",endpoints[j,1],", ",endpoints[j,2],"]")
                else
                    println(file,"I_",j,": (",endpoints[j,1],", ",endpoints[j,2],")")
                end
            end
            println(file,"\n",V-1," intervals")
            for j=numVIntervals+1: row
                if endpoints[j,1]==endpoints[j,2]
                    println(file,"J_",j-numVIntervals,": [",endpoints[j,1],", ",endpoints[j,2],"]")
                else
                    println(file,"J_",j-numVIntervals,": (",endpoints[j,1],", ",endpoints[j,2],")")
                end
            end
        end
        #Find the possible distribtions of muffins
        for i=1:length(permV)
            A=permV[i]
            sum_1=0
            sum_2=0
            equal = true
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
        if proof >=1
            println("\nPossible ",V," students:")
            for i=1:length(possV)
                PRINT(possV[i])
                println()
            end
            println("\nPossible ",V-1," students:")
            for i=1:length(possV_)
                PRINT(possV_[i])
                println()
            end
                println()
        end
        if proof >=3
            println(file,"\nPossible ",V," students:")
            for i=1:length(possV)
                PRINT(possV[i],file)
                println(file)
            end
            println(file,"\nPossible ",V-1," students:")
            for i=1:length(possV_)
                PRINT(possV_[i],file)
                println(file)
            end
                println(file)
        end
        #look for more gaps using student distribtions
        mat_V=transpose(hcat(possV...))
        mat_V_=transpose(hcat(possV_...))
        numVDistributions,numVIntervals=size(mat_V)
        numV_Distributions, numV_Intervals=size(mat_V_)

        gap_found1 = false
        gap_found2 =false

        endpoints, gap_found1 = findGaps(mat_V,1,numVIntervals, endpoints,m,s,denom,proof,true,file)
        if !gap_found1
            endpoints, gap_found2 = findGaps(mat_V_,numVIntervals+1,numVIntervals + numV_Intervals, endpoints,m,s,denom,proof,true,file)
        end
        if gap_found1||gap_found2
            gap_found = true
        else
            gap_found=false
        end
    end #end finding gaps

    if proof >=1
        println("NO MORE GAPS")
        println()

    end
    if proof >=3
        println(file,"NO MORE GAPS")
        println(file)

    end
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
    symmIntervals = Array{Int64,2}(undef,0,0) # row [i j] if interval i is symm to interval j
    row_endpt,col_endpt = size(endpoints)
    #intervals are symmetric if they are the same size and the same distance from either 1/2 or m/2s
    for i = 1:row_endpt
        for j=i+1:row_endpt
            if (endpoints[i,1]+endpoints[j,2])==denom && endpoints[i,2]-endpoints[i,1]==endpoints[j,2]-endpoints[j,1] || (endpoints[i,1]+endpoints[j,2])==(m*denom)//s && endpoints[i,2]-endpoints[i,1]==endpoints[j,2]-endpoints[j,1]
                if length(symmIntervals) == 0
                    symmIntervals = [i j]
                else
                    symmIntervals = [symmIntervals; i j]
                end
            end
        end
    end
    if proof>=1
        row, col= size(endpoints)
        numIntervals = row
        numVIntervals = 0
        for i =1: row
            if endpoints[i,2]<=x*denom
                numVIntervals = numVIntervals+1
            end
        end
        println("\nSYMMETRIC INTERVALS")
        row,col = size(symmIntervals)
        for i = 1:row
            a = symmIntervals[i,1]
            b = symmIntervals[i,2]
            if a > numVIntervals && b > numVIntervals
                println(" J_",a-numVIntervals," = J_",b-numVIntervals)
            elseif b > numVIntervals
                println(" I_",a," = J_",b-numVIntervals)
            elseif a > numVIntervals
                println(" J_",a-numVIntervals," = I_",b)
            else
                println(" I_",a," = I_",b)
            end
        end

        #display(symmIntervals)
    end
    if proof>=3
        row, col= size(endpoints)
        numIntervals = row
        numVIntervals = 0
        for i =1: row
            if endpoints[i,2]<=x*denom
                numVIntervals = numVIntervals+1
            end
        end
        println(file,"\nSYMMETRIC INTERVALS")
        row,col = size(symmIntervals)
        for i = 1:row
            a = symmIntervals[i,1]
            b = symmIntervals[i,2]
            if a > numVIntervals && b > numVIntervals
                println(file," J_",a-numVIntervals," = J_",b-numVIntervals)
            elseif b > numVIntervals
                println(file," I_",a," = J_",b-numVIntervals)
            elseif a > numVIntervals
                println(file," J_",a-numVIntervals," = I_",b)
            else
                println(file," I_",a," = I_",b)
            end
        end

        #display(symmIntervals)
    end
    if length(symmIntervals)==0
        if proof>=1
            println("No symmetric intervals")
            println("f(",m,", ",s,")  ≤ ",numerator(alpha),"/",denominator(alpha))
        end
        if proof >=3
            println(file,"No symmetric intervals")
            println(file,"f(",m,", ",s,")  ≤ ",numerator(alpha),"/",denominator(alpha))
            close(file)
        end
        if ret_endpts
            #denom = denom/2
            B = collect(alpha*denom:1:(1-alpha)*denom)


            row_e, col_e = size(Endpoints)
            while row_e -1 > 0
                B=filter(x-> x<=Endpoints[row_e-1,2]|| x>=Endpoints[row_e,1],B)
                row_e = row_e-1
            end
            #display(B)
            return true, B
        end
        return true
    end
    row_symm, col_symm = size(symmIntervals)
    if proof>=1
        println("\nPossible ",V," students:")
        for i=1:length(possV)
            PRINT(possV[i])
            println()
        end
        println("\nPossible ",V-1," students:")
        for i=1:length(possV_)
            PRINT(possV_[i])
            println()
        end
    end
    if proof>=1
        println(file,"\nPossible ",V," students:")
        for i=1:length(possV)
            PRINT(possV[i],file)
            println(file)
        end
        println(file,"\nPossible ",V-1," students:")
        for i=1:length(possV_)
            PRINT(possV_[i],file)
            println(file)
        end
    end
    #test and make sure there is at least one possible distribtion
    #for both the V and the V-1 shares
    if length(mat_V)==0
        if proof>=1
            println("\nNo possible distributions of ",V," shares: f(",m,", ",s,")  ≤ ",alpha)
        end
        if proof>=3
            println(file,"\nNo possible distributions of ",V," shares: f(",m,", ",s,")  ≤ ",alpha)
            close(file)
        end
        if ret_endpts
            B = collect(alpha*denom:1:(1-alpha)*denom)
            row_e, col_e = size(endpoints)
            while row_e -1 > 0
                B=filter(x-> x<=endpoints[row_e-1,2]|| x>=endpoints[row_e,1],B)
                row_e = row_e-1
            end
            B= filter(x-> x!= 1//2, B)
            i=1
            while i<length(B)
                if B[i]%2 !=0
                    deleteat!(B,i)
                end
                i = i+1
            end
            B= B/2
            return true, B
        end
        return true
    elseif length(mat_V_)==0
        if proof >=1
            println("\nNo possible distributions of ",V-1," shares: f(",m,", ",s,")  ≤ ",alpha)
        end
        if proof >=3
            println(file,"\nNo possible distributions of ",V-1," shares: f(",m,", ",s,")  ≤ ",alpha)
            close(file)
        end
        if ret_endpts
            #denom = denom/2
            B = collect(alpha*denom:1:(1-alpha)*denom)


            row_e, col_e = size(Endpoints)
            while row_e -1 > 0
                B=filter(x-> x<=Endpoints[row_e-1,2]|| x>=Endpoints[row_e,1],B)
                row_e = row_e-1
            end
            #display(B)
            return true, B
        end
        return true
    end
    #make the two matricies the same size
    #so we can form them into one system of eqns
    row_mat,col_mat = size(mat_V)
    row_mat_,col_mat_=size(mat_V_)
    row_mat_new=0
    while(row_mat_new < row_mat + row_mat_)
        Z= zeros(Int64,1,col_mat)
        mat_V = [mat_V; Z]
        row_mat_new,col_mat = size(mat_V)
    end
    row_mat_new=0
    while(row_mat_new < row_mat + row_mat_)
        Z = zeros(Int64,1,col_mat_)
        mat_V_ = [Z; mat_V_]
        row_mat_new,col_mat_ = size(mat_V_)
    end
    poss_Dist = [mat_V mat_V_]

    row_Dist, col_Dist = size(poss_Dist)
    row_endpt, col_endpt = size(endpoints)
    while col_Dist < row_endpt
        poss_Dist = [poss_Dist zeros(row_Dist)]
        row_Dist, col_Dist = size(poss_Dist)
    end

    A=Array{Int64}(undef,0) # A will become the matrix in the eqation Ax=b
    for i=1:row_symm
        append!(A, (poss_Dist[:,symmIntervals[i,1]]-poss_Dist[:,symmIntervals[i,2]]))
    end

    A=reshape(A,Int64(length(A)/row_symm),row_symm)


    A=transpose(A)
    #A is now a matrix with rows correspoding to the symmetric identities
    #ex if the intervals are I1, I2, I3, I4, I5 (|I1|=|I4|,|I2|=|I3|)
    #A = [I1-I4 ; I2-I3]

    A=unique(A, dims = 1)
    if proof == 2
        println("\nIf intervals i and j are symmetric, subtract col i from col j and add as a row to A")
        display(A)
    end
    if proof == 4
        println(file,"\nIf intervals i and j are symmetric, subtract col i from col j and add as a row to A")
        row, col = size(A)
        for i = 1:row
            for j=1:col
                @printf(file,"%4d  ",A[i,j])
            end
            println(file)
        end
    end
    row,col=size(A)
    #add zeros for each row of I₁ - I₂
    b=zeros(Int64,row)
    #sum up cols in mat_V
    #sum up cols in mat_V
    if length(mat_V)!=0
        Sum_1 = mat_V[:,1]
        r,c=size(mat_V)
        for i = 2:c
            Sum_1 = Sum_1 + mat_V[:,i]
        end
        #transpose to get right dim
        Sum_1=Sum_1'
        b=[b;Vshares]
        #append as row to A
        A=[A; Sum_1]
    end
    if length(mat_V_)!=0
        Sum_2 = mat_V_[:,1]
        r,c=size(mat_V_)
        for i = 2:c
            Sum_2 = Sum_2 + mat_V_[:,i]
        end
        #transpose to get right dim
        Sum_2=Sum_2'
        b=[b;V₋₁shares]
        #append as row to A
        A=[A;Sum_2]
    end
    if proof == 2
        println("\nAdd a row of ",V,"'s and a row of ",V-1,"'s \n(these when multiplied by number of each type of distribution of students will add to number of ",V," students and number of ",V-1," students)")
        display(A)
    end
    if proof == 4
        println(file,"\nAdd a row of ",V,"'s and a row of ",V-1,"'s \n(these when multiplied by number of each type of distribution of students will add to number of ",V," students and number of ",V-1," students)")
        row, col = size(A)
        for i = 1:row
            for j=1:col
                @printf(file,"%4d  ",A[i,j])
            end
            println(file)
        end
    end
    #Append a row of ones to A (num stud per distribtion adds to num VV students)
    Ones=(ones(Int64,col))'

    #We know one last bit of info
    #We know the amount of shares in the intervals symm about 1/2

    #find the interverals that are symmetric around 1/2
    #split_Intervals_1 are the intervals on the left of 1/2
    #split_Intervals_2 are the intervals on the right of 1/2
    split_Intervals_1 = Array{Int64}(undef,0)
    split_Intervals_2 = Array{Int64}(undef,0)
    if VV==V
        row,col=size(endpoints)
        for i = 1:row
            if endpoints[i,1]>=xbuddy*denom && endpoints[i,2]<=denom//2
                append!(split_Intervals_1, i)
            end
            if endpoints[i,1]>=denom//2 && endpoints[i,2]<=x*denom
                append!(split_Intervals_2,i)
            end
        end
    else
        row,col=size(endpoints)
        for i = 1:row
            if endpoints[i,1]>=y*denom && endpoints[i,2]<=denom//2
                append!(split_Intervals_1, i)
            end
            if endpoints[i,1]>=denom//2 && endpoints[i,2]<=ybuddy*denom
                append!(split_Intervals_2,i)
            end
        end
    end

    symm_sum_1 = Array{Int64}(undef,0)
    symm_sum_2 = Array{Int64}(undef,0)
    for i in split_Intervals_1
        if length(symm_sum_1) == 0
            symm_sum_1 = poss_Dist[:,i]
        else
            symm_sum_1 =symm_sum_1+ poss_Dist[:,i]
        end

    end
    for i in split_Intervals_2
        if length(symm_sum_2) == 0
            symm_sum_2 = poss_Dist[:,i]
        else
            symm_sum_2 =symm_sum_2 + poss_Dist[:,i]
        end
    end
    symm_sum_1=symm_sum_1'
    symm_sum_2 = symm_sum_2'
    if length(symm_sum_1)!=0
        A=[A;symm_sum_1]
        b=[b;num_split_shares]
    end
    if length(symm_sum_2) !=0
        A=[A;symm_sum_2]
        b=[b; num_split_shares]
    end


    if proof ==2
        println("\nAdd a rows for the sums of the ",VV," shares symmetric around 1/2")
        display(A)
        println("Then solve the system for the number of each type of distribution (= ",b,")")
        display(A)
    elseif proof ==1 || proof ==2
        println("Then solve the system for the number of each type of distribution (= ",b,")")
    end
    if proof ==4
        println(file,"\nAdd a rows for the sums of the ",VV," shares symmetric around 1/2")
        row, col = size(A)
        for i = 1:row
            for j=1:col
                @printf(file,"%4d  ",A[i,j])
            end
            println(file)
        end
        println(file,"Then solve the system for the number of each type of distribution (= ",b,")")
        row, col = size(A)
        for i = 1:row
            for j=1:col
                @printf(file,"%4d  ",A[i,j])
            end
            println(file)
        end
    elseif proof == 3
        println(file, "Then solve the system for the number of each type of distribution (= ",b,")")
    end

    row,col=size(A)

    model=Model(with_optimizer(Cbc.Optimizer, logLevel=0))
    @variable(model, X[i=1:col],Int)

    @constraint(model,con,A*X .==b)
    @constraint(model,con_1,X.>=0)

    optimize!(model)
    if(termination_status(model)==MOI.OPTIMAL)
        if proof >=1
            display(value.(X))
            println("\nThere is a solution on the naturals, failed to prove f(",m,", ",s,")  ≤ ",numerator(alpha),"/",denominator(alpha))
        end
        if proof >=3
            println(file," solution = ",value.(X))
            println(file,"\nThere is a solution on the naturals, failed to prove f(",m,", ",s,")  ≤ ",numerator(alpha),"/",denominator(alpha))
            close(file)
        end

        if ret_endpts
            #denom = denom/2
            B = collect(alpha*denom:1:(1-alpha)*denom)

            row_e, col_e = size(Endpoints)
            while row_e -1 > 0
                B=filter(x-> x<=Endpoints[row_e-1,2]|| x>=Endpoints[row_e,1],B)
                row_e = row_e-1
            end
            #display(B)
            return false, B
        end
        return false
    else
        if proof >=1
            println("\nThere is no solution on the naturals, f(",m,", ",s,")  ≤ ",numerator(alpha),"/",denominator(alpha))
        end
        if proof >=3
            println(file,"\nThere is no solution on the naturals, f(",m,", ",s,")  ≤ ",numerator(alpha),"/",denominator(alpha))
            close(file)
        end
        if ret_endpts
            #denom = denom/2
            B = collect(alpha*denom:1:(1-alpha)*denom)
            row_e, col_e = size(Endpoints)
            while row_e -1 > 0
                B=filter(x-> x<=Endpoints[row_e-1,2]|| x>=Endpoints[row_e,1],B)
                row_e = row_e-1
            end
            #display(B)

            return true, B
        end
        return true
    end
end

#Input
#possible distribution matrix
#front and back determine what rows in endpts to refer to
#endpoints is a matrix of known endpoints
#m,s (muffins, students)
#denom is denominator your working with
#proof = 1 or 2 if you want a proof printed
#gap_exteded = false if you don't want to include extending gaps TODO this can probably be removed we always consider it

#output
# newendpoints, _gap
# returns a matrix which includes any new gaps, and a bool which is true if a gap was found false otherwise
function findGaps(possDist,front,back,endpoints,m,s,denom,proof = 0, gap_extended=true, file = 0)
    #find gaps
    _gap=false
    #    println("NEW CALL")
    numDistributions, numIntervals = size(possDist)
    #we will use endpoints throughout and don't want it to change while working with it
    newendpoints=endpoints
    for j=front:back#1:numVIntervals
        #These arrays will hold the blocked off intervals
        #ex. we will know that [lowEnd[i], highEnd[i]] is not empty
        lowEnd = Array{Rational}(undef,0)
        highEnd = Array{Rational}(undef,0)
        for i=1:numDistributions
            #println("i = ",i,"  j = ",j)
            #display(possDist)
            if(possDist[i,j-front+1]!=0)
                #lower bound
                sum=0
                #use idea described at bottom of pg 148 to
                #find more gaps using student distributions

                for k=front:back
                    if(k!=j)
                        sum=sum+endpoints[k,2]*possDist[i,k-front+1]
                    else
                        sum=sum+endpoints[k,2]*(possDist[i,k-front+1]-1)
                    end
                end
                lowerbound_temp=(m//s)*denom-sum
                #upper bound
                sum=0
                for k=front:back
                    if(k!=j)
                        sum=sum+endpoints[k,1]*possDist[i,k-front+1]
                    else
                        sum=sum+endpoints[k,1]*(possDist[i,k-front+1]-1)
                    end
                end
                upperbound_temp=(m//s)*denom-sum
                if lowerbound_temp<endpoints[j,1]
                    lowerbound_temp=endpoints[j,1]
                end
                if upperbound_temp>endpoints[j,2]
                    upperbound_temp=endpoints[j,2]
                end
                append!(lowEnd,lowerbound_temp)
                append!(highEnd,upperbound_temp)
            end #end if(X[i][j]!=0)
        end  #end  for i in possInd

        #make temp endpts so we don't corrupt orignal data
        temp_endpt = Vector{Rational}(undef,0)

        lowEnd_temp = Array{Rational}(undef,0)
        for i=1:length(lowEnd)
            append!(lowEnd_temp,lowEnd[i])
        end
        #println(lowEnd)
        highEnd_temp = Array{Rational}(undef,0)
        for i=1:length(highEnd)
            append!(highEnd_temp,highEnd[i])
        end
        #only consider unique lower bounds
        while unique(lowEnd_temp)!=lowEnd_temp
            for i =1:length(lowEnd_temp)
                for j=i+1:length(lowEnd_temp)
                    if lowEnd_temp[i]==lowEnd_temp[j]
                        if highEnd_temp[i]>=highEnd_temp[j]
                            deleteat!(lowEnd_temp,j)
                            deleteat!(highEnd_temp,j)
                        else
                            deleteat!(lowEnd_temp,i)
                            deleteat!(highEnd_temp,i)
                        end
                        break
                    end
                end
            end
        end
        #now our data is of the form [lowEnd[1],highEnd[1],lowEnd[2],highEnd[2]...]
        while (length(lowEnd_temp)>0)
            lowIndex=1
            #println(lowEnd*230)
            for i = 1:length(lowEnd_temp)
                if lowEnd_temp[i]<lowEnd_temp[lowIndex]
                    lowIndex=i
                end
            end
            push!(temp_endpt,lowEnd_temp[lowIndex])
            push!(temp_endpt,highEnd_temp[lowIndex])
            deleteat!(lowEnd_temp,lowIndex)
            deleteat!(highEnd_temp,lowIndex)
        end

        #if we have overlaping intervals merge them
        merge=true
        while merge==true
            merge = false
            for i=1:Int64(length(temp_endpt)/2-1)
                if temp_endpt[2*i] >= temp_endpt[2*i+1]
                    merge=true
                    if (temp_endpt[2*i]<temp_endpt[2*(i+1)])
                        deleteat!(temp_endpt,2i+1)
                        deleteat!(temp_endpt,2i)
                    else
                        deleteat!(temp_endpt,2(i+1))
                        deleteat!(temp_endpt,2i+1)
                    end
                    break
                end
            end
        end

        #make sure we only have data in the orignal interval
        endpt=filter(x-> x > endpoints[j,1],temp_endpt)
        endpt=filter(x-> x < endpoints[j,2],endpt)

        #convert lowEnd and highEnd to int data
        int_lowEnd = Array{Any}(undef,length(lowEnd))
        for i=1:length(lowEnd)
            if lowEnd[i]-floor(lowEnd[i])!=0
                int_lowEnd[i] = float(lowEnd[i])
            else
                int_lowEnd[i] = Int64(lowEnd[i])
            end
        end
        lowEnd=int_lowEnd
        int_highEnd = Array{Any}(undef,length(highEnd))
        for i=1:length(highEnd)
            if highEnd[i]-floor(highEnd[i])!=0
                int_highEnd[i] = float(highEnd[i])
            else
                int_highEnd[i] = Int64(highEnd[i])
            end
        end
        highEnd=int_highEnd

        #print out the gaps found
        if length(endpt)%2==0
            for i=1:Int64(length(endpt)/2)
                upperbound=Int64(endpt[i])
                lowerbound=Int64(endpt[2i])
                if lowerbound>upperbound
                    if proof >=1
                        println("NEW GAP FOUND: [",upperbound,",",lowerbound,"]")
                        k=1
                        for i=1:numDistributions

                            if(possDist[i,j-front+1]!=0)
                            #    println(possDistV[i,:])
                                print(lowEnd[k]," ≤ ")
                                PRINT(possDist[i,:])
                                println(" ≤ ",highEnd[k])
                                k=k+1
                            end
                        end
                    end
                    if proof >=3
                        println(file,"NEW GAP FOUND: [",upperbound,",",lowerbound,"]")
                        k=1
                        for i=1:numDistributions

                            if(possDist[i,j-front+1]!=0)
                            #    println(possDistV[i,:])
                                print(file,lowEnd[k]," ≤ ")
                                PRINT(possDist[i,:],file)
                                println(file," ≤ ",highEnd[k])
                                k=k+1
                            end
                        end
                    end
                    _gap=true
                    newendpoints=[newendpoints; upperbound lowerbound]#; (1-lowerbound) (1-upperbound)]
                    newendpoints=sort(collect(Iterators.flatten(newendpoints)))
                    newendpoints=reshape(newendpoints,(2,Int64(length(newendpoints)/2)))
                    newendpoints = transpose(newendpoints)
                end #end if(lowerbound>upperbound)
            end#end for i = 1:Int64(length(temp_endpt)/2)
        elseif length(endpt)==1 && gap_extended

            if temp_endpt[1]<endpoints[j,1] #then change upper bound
                _gap=true
                if proof >=1
                    println("GAP EXTENDED: ",Int64(newendpoints[j,2])," changed to ",Int64(temp_endpt[2]))
                end
                if proof >=3
                    println(file,"GAP EXTENDED: ",Int64(newendpoints[j,2])," changed to ",Int64(temp_endpt[2]))
                end
                newendpoints[j,2]=temp_endpt[2]
                if proof >=1
                    k=1
                    for i=1:numDistributions
                        if(possDist[i,j-front+1]!=0)
                            print(lowEnd[k]," ≤ ")
                            PRINT(possDist[i,:])
                            println(" ≤ ",highEnd[k])
                            k=k+1
                        end
                    end
                end
                if proof >=3
                    k=1
                    for i=1:numDistributions
                        if(possDist[i,j-front+1]!=0)
                            print(file,lowEnd[k]," ≤ ")
                            PRINT(possDist[i,:],file)
                            println(file," ≤ ",highEnd[k])
                            k=k+1
                        end
                    end
                end
            elseif temp_endpt[1]>endpoints[j,1] && gap_extended
                _gap=true
                if proof >=1
                    println("GAP EXTENDED: ",Int64(newendpoints[j,1])," changed to ",Int64(temp_endpt[1]))
                end
                if proof >=3
                    println(file,"GAP EXTENDED: ",Int64(newendpoints[j,1])," changed to ",Int64(temp_endpt[1]))
                end
                newendpoints[j,1]=temp_endpt[1]
                if proof>=1
                    k=1
                    for i=1:numDistributions
                        if(possDist[i,j-front+1]!=0)
                            print(lowEnd[k]," ≤ ")
                            PRINT(possDist[i,:])
                            println(" ≤ ",highEnd[k])
                            k=k+1
                        end
                    end
                end

                if proof>=3
                    k=1
                    for i=1:numDistributions
                        if(possDist[i,j-front+1]!=0)
                            print(file,lowEnd[k]," ≤ ")
                            PRINT(possDist[i,:],file)
                            println(file," ≤ ",highEnd[k])
                            k=k+1
                        end
                    end
                end

            end

        end
    end #end for j=1:numInterval

   #move on to the V-1 Intervals
   #it's the same
return newendpoints, _gap

end

#input
#a matrix of endpoints, V (to determine to buddy or not), y, to determine if in 2 shares
#m,s,denom and proof are as before

#output a matrix which includes any new endpoints
function buddymatch(endpoints, V,y,m,s,denom,proof = 0, file = 0)
    row,col=size(endpoints)
    i=1
    while i!=row
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
        if buddyIn_1 == false && buddyIn_2 == false && denom-upper >=endpoints[1,1] && denom-lower <= endpoints[row,col]
            if proof >=1
                println("[",Float64(denom-upper),"  ",Float64(denom-lower),"] by buddying [",Float64(lower),"  ",Float64(upper),"]")
            end
            if proof >=3
                println(file,"[",Float64(denom-upper),"  ",Float64(denom-lower),"] by buddying [",Float64(lower),"  ",Float64(upper),"]")
            end
            endpoints=[endpoints; [(denom-upper) (denom-lower)]]
            i=1
            row,col=size(endpoints)
        end
        endpoints=sort(collect(Iterators.flatten(endpoints)))
        endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
        endpoints = transpose(endpoints)
        if lower>=y && V==3 #in the two shares
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
                if proof >=1
                    println("[",Float64((m//s)*denom-upper),"  ",Float64((m//s)*denom-lower),"] by matching [",Float64(lower),"  ",Float64(upper),"]")
                end
                if proof >=3
                    println(file,"[",Float64((m//s)*denom-upper),"  ",Float64((m//s)*denom-lower),"] by matching [",Float64(lower),"  ",Float64(upper),"]")
                end
                endpoints=[endpoints; [((m//s)*denom-upper) ((m//s)*denom-lower)]]
                i=1
                row,col=size(endpoints)
            end
        end
        endpoints=sort(collect(Iterators.flatten(endpoints)))
        endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
        endpoints = transpose(endpoints)
        i=i+1
    end
    return endpoints
end

#PRINT prints an array in "converted" form
#ex: array = [2 0 0 1] (lets say that means use 2 I_1 shares and 1 I_4 shares)
# this function would print (1 1 4)
function PRINT(array, file = 0)
    if file != 0
        print(file,"( ")
        for i=1:length(array)
            for j=1:array[i]
                print(file,i," ")
            end
        end
        print(file,")")
    else
        print("( ")
        for i=1:length(array)
            for j=1:array[i]
                print(i," ")
            end
        end
        print(")")
    end
end

#This just converts from rational to int (or float)
function convert_Int(matrix)
    row,col=size(matrix)
    mat = Matrix{Any}(undef, row,col)
    for i=1:row
        for j=1:col
            if matrix[i,j]-floor(matrix[i,j])!=0
                mat[i,j]=float(matrix[i,j])
            else
                mat[i,j]=Int64(matrix[i,j])
            end
        end
    end
    return mat
end

#input
#m,s (muffins, students) min_al current best upper bound
#ret_endpts = true if you want endpoints back

#output alpha, endpoints (if ret_endpts  = true)
function GAP(m,s,min_al = 1//2, ret_endpts = false)
    if m%s==0
        return 1,1
    end
    if min_al >1//2
        min_al = 1//2
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
    filter!(x -> 1/3 <= x , array)
    #println(array)
    #filter!(x ->  x <min_al, array)

    for i = 1:length(array)
        if array[i]==49//144
        #    println(i)
        end
    end
#    println(array)
    alpha, endpoints = binarySearchGAP(m,s,array)
    if alpha == 1
        return alpha,1
    end
    if ret_endpts && endpoints!=0
        #endpoints=sort(collect(Iterators.flatten(endpoints)))
        temp=Array{Int64}(undef,0)
        for i = 1:length(endpoints)
            try
                append!(temp,convert(Int64,endpoints[i]))
            catch

            end
        end
        endpoints = temp
        return alpha, endpoints
    elseif ret_endpts
        return alpha,1
    else
        return alpha
    end
end
#binary search, used in GAP uses PROC and VGAP to find the
#alpha is a sorted array

function binarySearchGAP(m,s,array)
    alpha, endpoints =binSearchHelpGAP(m,s,array,1,length(array))
    return alpha,endpoints
end
function binSearchHelpGAP(m,s,array,front,back)
    #println("m: ",m,"  s: ",s)

    if front>back
        return 1
    end
    guessIndex =Int64(floor((front+back)/2))
    #println("front: ",front,"  back: ",back,"  alpha: ",array[guessIndex])
    if guessIndex == front
        valid, endpoints =VGAP(m,s,array[back], 0, true)
        if valid==true
            return array[back], endpoints
        else
            return 1, 1
        end
    end

    if VGAP(m,s,array[guessIndex]) == true
        return binSearchHelpGAP(m,s,array,front,guessIndex)
    else
        return binSearchHelpGAP(m,s,array,guessIndex, back)
    end
end

function pGap(m,s,alpha)
    x=GAP(m,s)
    println("f(",m,",",s,")    ",x==alpha,"   ",x)
end
