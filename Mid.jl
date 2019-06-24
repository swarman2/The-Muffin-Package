include("src\\permutations.jl")
include("Half.jl")
using JuMP
using GLPK
function MID_proof(m,s,alpha)
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
        I1=num_small_shares
        I2=Int64(num_large_shares/2)
        I3=I2
        X = perm(V,3)
    #    println(X)
        possInd = Array{Int64}(undef,0)
        for i=1:length(X)
            A=X[i]
            if (A[1]*alpha+A[2]*xbuddy+A[3]*1//2 < m/s) && (A[1]*ybuddy+A[2]*1/2+A[3]*x> m/s)
                append!(possInd, i)
            end
        end
    S=sᵥ
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

        I1=Int64(num_small_shares/2)
        I2=I1
        I3=num_large_shares
        X = perm(V-1,3)
    #    display(X)
    #    println(X)
        possInd = Array{Int64}(undef,0)
        for i=1:length(X)
            A=X[i]
    #        println("Test ",A,": ",A[1]*y+A[2]*1/2+A[3]*xbuddy, "<", m/s, " && ",A[1]*1/2+A[2]*ybuddy+A[3]*(1-alpha),">", m/s)
            if (A[1]*y+A[2]*1/2+A[3]*xbuddy < m/s) && (A[1]*1/2+A[2]*ybuddy+A[3]*(1-alpha)> m/s)
                append!(possInd, i)
            end
        end
        S=sᵥ₋₁
    end
    matrix=(X[possInd[1]])
#    println(matrix)
    for i=2:length(possInd)
        Y=(X[possInd[i]])
#        println(Y)
        matrix=[matrix Y]
    end
#    display(matrix)
    row,col=size(matrix)
#print(size((ones(Int64,(col)))))
    matrix=[matrix; transpose(ones(Int64,(col)))]
# print(size(matrix))
    println("A:")
    display(matrix)



    m=Model(with_optimizer(GLPK.Optimizer))
    @variable(m, x[i=1:length(possInd)],Int)
    #print("size of x: ",size(x))
    b=[I1;I2;I3;S]
    println("b: ",b)
    @constraint(m,con,matrix*x .==b)
    optimize!(m)
    if(!has_values(m))
        println("No solution to this system on the Naturals: 𝛂 ≤ ",alpha)
    else
        println("This solution has a system on the Naturals: INCONCLUSIVE")
    end
end

function MID(m,s,alpha)
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

        I1=num_small_shares
        I2=Int64(num_large_shares/2)
        I3=I2
        X = perm(V,3)
    #    println(X)
        possInd = Array{Int64}(undef,0)
        for i=1:length(X)
            A=X[i]
            if (A[1]*alpha+A[2]*xbuddy+A[3]*1//2 < m/s) && (A[1]*ybuddy+A[2]*1/2+A[3]*x> m/s)
                append!(possInd, i)
            end
        end
    S=sᵥ
    else
        num_small_shares = V₋₁shares-Vshares
       num_large_shares =   Vshares

        I1=Int64(num_small_shares/2)
        I2=I1
        I3=num_large_shares
        X = perm(V-1,3)
    #    println(X)
        possInd = Array{Int64}(undef,0)
        for i=1:length(X)
            A=X[i]
            if (A[1]*y+A[2]*1/2+A[3]*xbuddy < m/s) && (A[1]*1/2+A[2]*ybuddy+A[3]*(1-alpha)> m/s)
                append!(possInd, i)
            end
        end
        S=sᵥ₋₁
    end
    if(length(possInd)==0)
        return Inf
    end
    matrix=(X[possInd[1]])
#    println(matrix)
    for i=2:length(possInd)
        Y=(X[possInd[i]])
#        println(Y)
        matrix=[matrix Y]
    end
#    display(matrix)
    row,col=size(matrix)
#print(size((ones(Int64,(col)))))
    matrix=[matrix; transpose(ones(Int64,(col)))]
# print(size(matrix))


    m=Model(with_optimizer(GLPK.Optimizer))
    @variable(m, x[i=1:length(possInd)],Int)
    #print("size of x: ",size(x))
    b=[I1;I2;I3;S]
    @constraint(m,con,matrix*x .==b)
    optimize!(m)
    if(has_values(m))
        println("alpha could be greater than ",alpha)
    else
        println("alpha ≤ ",alpha)
    end
end

#returns all permutations of numbers 1 through n of size
# r that sum to n
#(ex 2,3:
#   [[0,0,2],[0,1,1],[0,2,0],[1,0,1],[1,1,0],[2,0,0]]
#)

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

MID(23,13,53//130)
MID(23,14,17//42)
MID(33,20,49//120)
MID(37,21,103//252)
MID(59,14,131//280)
