function Test(A,B)
#    println(1+floor(A-B*(59//126)))
#println()
#println(A + B*(-67//126)," ≥ ",0,"    ", A + B*(-67//126) >= 0)
#println((A+1 -B + B*(41//90))*4," ≤ ", A+B*(-67//126),"    ",(A+1 -B + B*(41//90))*4 <= A+B*(-67//126))
#println(B-(606//1260)*B +1, " ≤ ", A,"    ",B-(606//1260)*B +1 <= A)
#println(floor(A - B*(59//126))," < ", 1+ floor(2*A-B),"    ",floor(A - B*(59//126)) < 1+ floor(2*A-B))
#println()
if (A + B*(-67//126) >= 0 && (A+1 -B + B*(41//90))*4 <= A+B*(-67//126) ) && B-(606//1260)*B +1 <= A &&floor(A - B*(59//126)) < 1+ floor(2*A-B) && floor(A - B*(59//126)) < 1+ floor(2*A-B)
    return A,B
else
    return 0,0
end
end
function foo()
aa = 0//90
bb = 0//90
min = Inf
A = 41*90//90
B = 0//90
aaa = A
if true
while A < 43
    B = 0//90
    if A >= 1+aaa
        aaa=A
        println(A)
    end
    if A == 335//8
        print(A," ")
    end
    while B< 200
#println(A,"  ",B)
        if A == 335//8 && B == 630//8
            println(B)
        end
        if A + B*(-67//126) >= 0
            a,b = Test(A,B)
            if a!=0 || b!=0
                if min > 1+floor(A-B*(59//126))
                    #println("A: ",A," B: ",B)
                    min =  1+floor(A-B*(59//126))
                    println(A,"  ", B)
                    aa = A
                    bb = B
                end
            end
        end
        B = B+1//90

    end
    A= A+1//90
end
end
#aa = 335//8
#aa = 630//8
println("A: ",aa," B: ",bb)
eps = .000001
while eps < 1
    if 2*(1+floor(aa - ((bb) * (59//126 + eps))))<12
        println("yes")
        break
    else
        println(2*(1+floor(aa - ((bb) * (59//126 + eps)))))
        println(2*(1+floor(aa - ((bb) * (59//126)))))
        println(aa - ((bb) * (59//126 + eps)))
        println(aa - ((bb) * (59//126+0.0)))
        println("no")
    end
    eps = eps +.01
end
#Test(8,1)
end
foo()
#Test(335//8,630//8)
