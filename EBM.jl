function EBM(m::Int,s::Int)

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

    #3
    if 2d + 1 <= a && a <= 3d
        return 1//3
    end

    #4
    if a == 2d
        return a//2
    end

    #5
    if 1 <= a <= 2d - 1
        X=minimum([a//2 (a+d)//4])
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
