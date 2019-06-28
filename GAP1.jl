include("helper_functions.jl")
using JuMP
using GLPK

function VGAPV3(m,s,alpha)

    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁

    x,y=FINDEND(m,s,alpha,V)

    xbuddy = 1-x
    ybuddy = 1-y
    endpoints=[alpha 1-alpha]
    endpoints=[endpoints;x y; xbuddy ybuddy]

    if xbuddy<ybuddy
        lower=xbuddy
        higher=ybuddy
    else
        lower=ybuddy
        higher=xbuddy
    end
    while lower > y #in the 2 shares
        #match
        low_match=m//s-higher
        high_match=m//s-lower
        lower = 1-high_match
        higher = 1-low_match
        endpoints=[endpoints; low_match high_match ;  lower higher]
        _gap=true
        endpoints=sort(collect(Iterators.flatten(endpoints)))
        endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
        endpoints = transpose(endpoints)
    end

    row,col=size(endpoints)
    for i=1:row
        if endpoints[i,1]+ endpoints[i,2]==1
            endpoints=[endpoints; 1//2 1//2 ]
        end
    end
    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)

    i=1
    while endpoints[i,1] <= m//(2s)
        i=i+1
    end
    j=1
    while endpoints[j,1] <= 1-m//(2s)
        j=j+1
    end

    endpoints=[endpoints; m//(2s) m//(2s); m//(2s) m//(2s); 1-m//(2s) 1-m//(2s);1-m//(2s) 1-m//(2s)]
    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    symmIntervals = Array{Int64,2}(undef,0,0) # row [i j] if interval i is symm to interval j
    row_endpt,col_endpt = size(endpoints)
    for i = 1:row_endpt
        for j=i+1:row_endpt
            if ((endpoints[i,1]+endpoints[j,2])== 1 && endpoints[i,2] + endpoints[j,1] == 1) || ((endpoints[i,1]+endpoints[j,2])== m//s && endpoints[i,2] + endpoints[j,1] == m//s && V-1==2)
                if length(symmIntervals) == 0
                    symmIntervals = [i j]
                else
                    symmIntervals = [symmIntervals; i j]
                end
            end
        end
    end
#    display(symmIntervals)
    numVIntervals = 0
    for i =1: row_endpt
        if endpoints[i,2]<=x
            numVIntervals = numVIntervals+1
        end
    end
    permV=perm(V, numVIntervals)
    if(row_endpt == numVIntervals)
        return 1
    end
    permV_ = perm(V-1, row_endpt - numVIntervals)
    #display(permV)
    #println("*********")
    #display(permV_)

    possV = Vector{Vector{Int64}}()
    possV_ = Vector{Vector{Int64}}()

    #Find the possible distribtions of muffins
    for i=1:length(permV)
        A=permV[i]
        sum_1=0
        sum_2=0
        for j =1:numVIntervals
            sum_1=sum_1+A[j]*endpoints[j,1]
            sum_2=sum_2+A[j]*endpoints[j,2]
        end
        if (sum_1 < m//s) && (sum_2 > m//s)
            append!(possV, permV[i,:])
        end
    end

    #Find the possible distribtions of muffins
    for i=1:length(permV_)
        A=permV_[i]
        sum_1=0
        sum_2=0
        for j =numVIntervals+1:row_endpt
            sum_1=sum_1+A[j-numVIntervals]*endpoints[j,1]
            sum_2=sum_2+A[j-numVIntervals]*endpoints[j,2]
        end
        if (sum_1 < m//s) && (sum_2 > m//s)
            append!(possV_, permV_[i,:])
        end
    end

#    println(possV_)
    mat_V=transpose(hcat(possV...))
    mat_V_=transpose(hcat(possV_...))

    #display(mat_V)
    #display(mat_V_)
    if(length(permV)!=0 && length(permV_)!=0)
        row_V, col_V = size(mat_V)
        row_V_, col_V_=size(mat_V_)
        while row_V < row_V_
            Z = transpose(zeros(Int64, col_V))
            mat_V=[mat_V; Z]
            row_V, col_V = size(mat_V)
        end
        while row_V_ < row_V
            Z = transpose(zeros(Int64, col_V_))
            mat_V_=[mat_V_; Z]
            row_V_, col_V_ = size(mat_V_)
        end
        poss_Dist = [mat_V mat_V_]
    elseif length(permV)!=0
        poss_Dist = mat_V
    else
        poss_Dist = mat_V_
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
    row,col=size(symmIntervals)
    #add zeros for each row of I₁ - I₂
    b=zeros(Int64,row)

    #sum up cols in mat_V
    Sum_1 = mat_V[:,1]
    r,c=size(mat_V)
    for i = 2:c
        Sum_1 = Sum_1 + mat_V[:,i]
    end

    Sum_2 = mat_V_[:,1]
    r,c=size(mat_V_)
    for i = 2:c
        Sum_2 = Sum_2 + mat_V_[:,i]
    end
    b=[b;sᵥ; sᵥ₋₁]
        #transpose to get right dim
    Sum_1=Sum_1'
    Sum_2=Sum_2'

    #append as rows to A
    A=[A; Sum_1]
    A=[A;Sum_2]
    #Append a row of ones to A (num stud per distribtion adds to num VV students)
    row,col=size(A)
    Ones=(ones(Int64,col))'
    if(V!=3)
        A=vcat(A, Ones)
    end


    #display(A)
    model=Model(with_optimizer(GLPK.Optimizer))
    @variable(model, x[i=1:col],Int)

    @constraint(model,con,A*x .==b)
    @constraint(model,con_1,x.>=0)

    optimize!(model)
    if(has_values(model))
        return false
    else
        return true
    end

end

function GAPV3(m,s,alpha)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁

    x,y=FINDEND(m,s,alpha,V)

    xbuddy = 1-x
    ybuddy = 1-y
    endpoints=[alpha 1-alpha]
    endpoints=[endpoints;x y; xbuddy ybuddy]
    println("ENDPOINTS: ")
    display(endpoints)
    if xbuddy<ybuddy
        lower=xbuddy
        higher=ybuddy
    else
        lower=ybuddy
        higher=xbuddy
    end
    while lower > y #in the 2 shares
        #match
        low_match=m//s-higher
        high_match=m//s-lower
        lower = 1-high_match
        higher = 1-low_match
        endpoints=[endpoints; low_match high_match ;  lower higher]
        _gap=true
        endpoints=sort(collect(Iterators.flatten(endpoints)))
        endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
        endpoints = transpose(endpoints)
    end

    row,col=size(endpoints)
    for i=1:row
        if endpoints[i,1]+ endpoints[i,2]==1
            endpoints=[endpoints; 1//2 1//2 ]
        end
    end
    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)

    i=1
    while endpoints[i,1] <= m//(2s)
        i=i+1
    end
    j=1
    while endpoints[j,1] <= 1-m//(2s)
        j=j+1
    end

    endpoints=[endpoints; m//(2s) m//(2s); m//(2s) m//(2s); 1-m//(2s) 1-m//(2s);1-m//(2s) 1-m//(2s)]
    endpoints=sort(collect(Iterators.flatten(endpoints)))
    endpoints=reshape(endpoints,(2,Int64(length(endpoints)/2)))
    endpoints = transpose(endpoints)
    println("ENDPOINTS AFTER BUDDY-MATCH: ")

    display(endpoints)

    symmIntervals = Array{Int64,2}(undef,0,0) # row [i j] if interval i is symm to interval j
    row_endpt,col_endpt = size(endpoints)
    for i = 1:row_endpt
        for j=i+1:row_endpt
            if ((endpoints[i,1]+endpoints[j,2])== 1 && endpoints[i,2] + endpoints[j,1] == 1) || ((endpoints[i,1]+endpoints[j,2])== m//s && endpoints[i,2] + endpoints[j,1] == m//s && V-1==2)
                if length(symmIntervals) == 0
                    symmIntervals = [i j]
                else
                    symmIntervals = [symmIntervals; i j]
                end
            end
        end
    end
#    display(symmIntervals)
    numVIntervals = 0
    for i =1: row_endpt
        if endpoints[i,2]<=x
            numVIntervals = numVIntervals+1
        end
    end
    permV=perm(V, numVIntervals)
    permV_ = perm(V-1, row_endpt - numVIntervals)
    #display(permV)
    #println("*********")
    #display(permV_)

    possV = Vector{Vector{Int64}}()
    possV_ = Vector{Vector{Int64}}()

    #Find the possible distribtions of muffins
    for i=1:length(permV)
        A=permV[i]
        sum_1=0
        sum_2=0
        for j =1:numVIntervals
            sum_1=sum_1+A[j]*endpoints[j,1]
            sum_2=sum_2+A[j]*endpoints[j,2]
        end
        if (sum_1 < m//s) && (sum_2 > m//s)
            append!(possV, permV[i,:])
        end
    end

    #Find the possible distribtions of muffins
    for i=1:length(permV_)
        A=permV_[i]
        sum_1=0
        sum_2=0
        for j =numVIntervals+1:row_endpt
            sum_1=sum_1+A[j-numVIntervals]*endpoints[j,1]
            sum_2=sum_2+A[j-numVIntervals]*endpoints[j,2]
        end
        if (sum_1 < m//s) && (sum_2 > m//s)
            append!(possV_, permV_[i,:])
        end
    end

#    println(possV_)
    mat_V=transpose(hcat(possV...))
    mat_V_=transpose(hcat(possV_...))

    #display(mat_V)
    #display(mat_V_)
    if(length(permV)!=0 && length(permV_)!=0)
        row_V, col_V = size(mat_V)
        row_V_, col_V_=size(mat_V_)
        while row_V < row_V_
            Z = transpose(zeros(Int64, col_V))
            mat_V=[mat_V; Z]
            row_V, col_V = size(mat_V)
        end
        while row_V_ < row_V
            Z = transpose(zeros(Int64, col_V_))
            mat_V_=[mat_V_; Z]
            row_V_, col_V_ = size(mat_V_)
        end
        poss_Dist = [mat_V mat_V_]
    elseif length(permV)!=0
        poss_Dist = mat_V
    else
        poss_Dist = mat_V_
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
    row,col=size(symmIntervals)
    #add zeros for each row of I₁ - I₂
    b=zeros(Int64,row)

    #sum up cols in mat_V
    Sum_1 = mat_V[:,1]
    r,c=size(mat_V)
    for i = 2:c
        Sum_1 = Sum_1 + mat_V[:,i]
    end

    Sum_2 = mat_V_[:,1]
    r,c=size(mat_V_)
    for i = 2:c
        Sum_2 = Sum_2 + mat_V_[:,i]
    end
    b=[b;sᵥ; sᵥ₋₁]
        #transpose to get right dim
    Sum_1=Sum_1'
    Sum_2=Sum_2'

    #append as rows to A
    A=[A; Sum_1]
    A=[A;Sum_2]
    #Append a row of ones to A (num stud per distribtion adds to num VV students)
    row,col=size(A)
    Ones=(ones(Int64,col))'
    if(V!=3)
        A=vcat(A, Ones)
    end


    #display(A)
    model=Model(with_optimizer(GLPK.Optimizer))
    @variable(model, x[i=1:col],Int)

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

end
