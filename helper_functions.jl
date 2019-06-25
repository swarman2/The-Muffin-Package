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
function perm(n,r)
    A=Array{Int64,1}(undef,r*(n+1))
    for i=1:length(A)
        A[i]=floor((i-1)/r)
    end
#    print(A)

    X=(collect(multiset_permutations(A,r)))
#    display(X)
#    println("********************")
    X=filter(x->sum(x)==n,X)
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
