function EBM(m::Int64,s::Int64, proof = false)
    #1
    if m%s == 0
      return 1
    end

    #2
    if ceil(2m//s)>=4
        return 1
    end

    d = Int64(m - s)
    k =Int64(floor((s/(3d))))
    if s%3d ==0
        k = k-1
    end
    a = Int64(s-3d*k)
    if proof
        println("d = ",d," k = ",k," a = ",a)
    end
    #3
    if 2d + 1 <= a && a <= 3d
        if proof
        println("2 * ",d," + 1 ≤ ",a," and ", a," ≤ 3 * ",d)
        println("α ≤ 1/3")
    end
        return 1//3
    end

    #4
    if a == 2d
        if proof
        println(a," = 2 * ",d)
        println("α ≤ ",a,"/2")
    end
        return a//2
    end

    #5
    if 1 <= a <= 2d - 1

        X=minimum([a//2 (a+d)//4])
        if proof
        println("1 ≤ ",a," ≤ 2 * ",d," - 1")
        println("X = min(",a,"/2, ",a+d,"/4)")
        println("X = ",X)
        println("α ≤ ","(",d," * ",k," + ",X,")/(3 * ",d," * ",k," + ",a,") = ",(d*k +X)//(3d*k+a))
    end
        return (d*k +X)//(3d*k+a)
end
end

#Pebm(29,27)
#Pebm(33,31)
#Pebm(35,33)
#Pebm(37,34)
#Pebm(39,37)
#Pebm(41,39)
#Pebm(44,41)
#Pebm(45,43)
#Pebm(46,43)
#Pebm(47,45)
#Pebm(49,45)
#Pebm(49,46)
#Pebm(50,47)
#Pebm(51,49)
#Pebm(53,50)
#Pebm(13,12)
