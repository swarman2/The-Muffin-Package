COND(X,a,d)
    if a/3<=X && X<= min(a/2,(a+d)/4)
        return true
    else
        return false
    end
end

HBM(a,d,k)
   possX = Array{Int64}(undef,0)
    if d>=1 && k>=1 && a<=2*d && gcd(a,d)==1
        #10
        X=max((a+2d)//6,(2a-d)//3)
        if COND(X,a,d)
            append!(possX,X)
        end
        #11
        if a!=1 || d!=1
            X=max((a+d)//5,(2a-d)//3,d//2)
            if COND(X,a,d)
                append!(possX,X)
            end
        end
        #12
        if a != (7d//5)
            X=max((3a-2d)//4,(a+2d)//6)
            if COND(X,a,d)
                append!(possX,X)
            end
        end
        #13
        if a<7d//5 || a> 5d//3
            X = max(a-d,(a+2d)//6)
            if COND(X,a,d)
                append!(possX,X)
            end
        end
        #14
        if a!=1 || d!=1
            X=max((2a-d)//3,d//2)
            if COND(X,a,d)
                append!(possX,X)
            end
        end
        #15
        if a<d || a>7d//5
            X=max((3a-2d)//4,(a+d)//5,d//2)
            if COND(X,a,d)
                append!(possX,X)
            end
        end
        #16
        if a!=1 || d!=1
            X= max((2a-d)//3, (a+d)//5)
            if COND(X,a,d)
                append!(possX,X)
            end
        end
        #17
        if a<7d//5
            X=(a+2d)//6
            if COND(X,a,d)
                append!(possX,X)
            end
        end
        #18
        if d>a || a>7d//3
            X=max(a-d,(a+d)//5,d//2)
            if COND(X,a,d)
                append!(possX,X)
            end
        end
        #19
        if a<d//2 || a>7d//5
            X=max((3a-2d)//4,d//2)
            if COND(X,a,d)
                append!(possX,X)
            end
        end
        #20
        if a>d || a<d//3
            X=max((2a-d)//3, (a+2d)//8)
            if COND(X,a,d)
                append!(possX,X)
            end
        end
        #21
        if a>7d//5|| a<d
            X=max((a+d)//5),(3a-2d)//4)
            if COND(X,a,d)
                append!(possX,X)
            end
        end
        if length(possX)==0
            println("No possible X")
            return Inf
        end
        X=minimum(possX)
        return (d*k+X)//(3*d*k+a)
end
