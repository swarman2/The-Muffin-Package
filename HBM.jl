function COND(X,a,d)
    if a/3<=X && X<= min(a/2,(a+d)/4)
        return true
    else
        return false
    end
end

function HBM(m,s)
d = Int64(m-s)
k =Int64( floor(s/(3d)))
a=Int64(s-3d*k)
   possX = Array{Rational}(undef,0)
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
            X=max((a+d)//5,(3a-2d)//4)
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
end
function Phbm(m,s)
    println("f(",m,", ",s,") ≤ ", HBM(m,s))
end

Phbm(25,22)
Phbm(33,29)
Phbm(34,31)
Phbm(38,31)
Phbm(43,35)
Phbm(41,36)
Phbm(43,40)
Phbm(49,40)
Phbm(45,41)
Phbm(49,43)
Phbm(55,48)
Phbm(59,48)
Phbm(52,49)
Phbm(60,49)
Phbm(57,50)
