
#mulitsets of B that sum to T of size k
function Multiset(B,T,k, time_limit=Inf)
    A=Dict{String,Vector{Vector{Int64}}}()

    i=length(B)
    A,sol= Multiset_helper(B,A,i,T,k,time(),time_limit)
    if A==-1
        return "time out"
    end

    return sol
end
function Multiset_helper(B,A,i::Int64,T::Int64,k::Int64, ogTime,time_limit)
    currTime = time()
    if currTime-ogTime>time_limit
        return -1,-1
    end
    if T <=0 || i == 0 || k==0
        return A,0
    elseif B[1]*k ==T
        array=Array{Int64}(undef, k)
        for j=1:length(array)
            array[j]=B[1]
        end
        return A,[array]
    elseif i==1 && T%k != B[1]
        return A,0
    elseif T%B[1]==0 && i==1
        return A,0
    elseif k==1
        inFirst = false
        for j=1:i
            if B[j]==T
                inFirst = true
            end
        end
        if inFirst
            return A,[[T]]
        else
            return A,0
        end
    end
    str=string(i)*string(T) * string(k)
    if(haskey(A,str))
        return A, A[str]
    else
        sol=Vector{Vector{Int64}}()
        A, Temp1=Multiset_helper(B,A, i-1,T,k,ogTime,time_limit)
        if A == -1
            return -1,-1
        end
        A, Temp2=Multiset_helper(B,A, i,Int64(T-B[i]),k-1,ogTime,time_limit)
        if A == -1
            return -1,-1
        end

        if Temp1 !=0 && Temp2 !=0
           sol=Vector{Vector{Int64}}()
           row = length(Temp2)

           for j = 1:row
               push!(sol, [B[i]])
               for k =1: length(Temp2[j])
                   push!(sol[j], Temp2[j][k])
               end
           end

            row = length(Temp1)
            for j = 1:row
               push!(sol,Temp1[j])
            end
        elseif Temp1 !=0
           sol=Vector{Vector{Int64}}()
           row = length(Temp1)
            for j = 1:row
               push!(sol,Temp1[j])
            end
        elseif Temp2 !=0
            sol=Vector{Vector{Int64}}()
            row = length(Temp2)
            for j = 1:row
                push!(sol, [B[i]])
                for k=1:length(Temp2[j])
                    push!(sol[j], Temp2[j][k])
                end
            end
        else
            return A,0
        end
        str=string(i)*string(T) * string(k)
        A[str] =vec(sol)
        return A, sol
    end
end

function print_Intervals(m,s,alpha)
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
        num_split_shares=Int64(num_large_shares/2)
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
        num_split_shares=Int64(num_small_shares/2)
        endpoints =Array{Rational, 2}(undef,0,0)
        endpoints =[y 1//2]
        endpoints = [endpoints;[1//2 ybuddy]; [xbuddy (1-alpha)]]
    end
    return endpoints
end

#returns all permutations of numbers 1 through n of size
# r that sum to n
#(ex 2,3:
#   [[0,0,2],[0,1,1],[0,2,0],[1,0,1],[1,1,0],[2,0,0]]
#)
function perm(n,r)
    B=Vector{Int64}(undef,0)
    if r==0
        return B
    end
    X=Vector{Vector{Int64}}()
    for i=1:r
        push!(B,0)
    end
    l=1
    while (true)
        B[1]=B[1]+1
        i=1

        while i <r
            if sum(B)==n+1#B[i]==n+1
                B[i]=0
                B[i+1]=B[i+1]+1
            end
            i = i +1
        end
        if sum(B) == n
        #    println(B)
            push!(X,[B[1]])
            for k=2:length(B)
                append!(X[l],B[k])
            end
            l=l+1
        #    display(X)
        end
        if sum(B) > n
            break;
        end
    end
    return X

end

#SV takes m,s as input and returns
#V, sᵥ, and sᵥ₋₁
#using the V-Conjecture (pg 66) and
#a linear function from pg 68
#gotten by solving the system:
#   Vsᵥ + (V-1)(sᵥ₋₁) = pieces = 2m
#   sᵥ + sᵥ₋₁ = s

function SV(m,s)
  V=Int64(ceil(2*m/s))
  sᵥ₋₁=V*s-2*m
  sᵥ = 2*m - s*(V-1)
  return V, sᵥ,sᵥ₋₁
end

#FINDEND input m,s alpha and V
#this returns x and y where all V shares are in
#the interval (𝛂, x) and all V-1 shares are in the
#interval (y, 1-𝛂)
#This is adapted from psudocode on page 69
function FINDEND(m,s,alpha,V)
  y=m//s -(1-alpha)*(V-2)
  if y>= (1-alpha)
    y=1-alpha
  end
  if y <= alpha
    y=alpha
  end
  x=m//s - alpha*(V-1)
  if x<=alpha
    x=alpha
  end
  if x>= (1-alpha)
    x=1-alpha
  end
  return x,y
end
