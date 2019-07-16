include("helper_functions.jl")
include("GAP.jl")
function VTRAIN(m,s,alpha)
    V,sᵥ,sᵥ₋₁=SV(m,s)
    Vshares=V*sᵥ
    V₋₁shares=(V-1)*sᵥ₋₁
    x,y=FINDEND(m,s,alpha,V)
    ybuddy=1-y
    xbuddy=1-x
    denom=denominator(alpha)
    denom=lcm(s,denom)*2
    if V₋₁shares > Vshares
        endpoints   = [y*denom (1-alpha)*denom]
        #display(endpoints)
        endpoints = GAP_int(m,s,alpha)
        #display(endpoints)
        #Find the possible distribtions of muffins
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



#VTRAIN(11,5,13//30)
VTRAIN(31,27,103//297)
