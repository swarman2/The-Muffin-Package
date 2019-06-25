function EBM(m::Int,s::Int)
    if ceil(2m/s)>=4
        return 1
    end
    if m%s == 0
      return 1
    end
    d = Int64(m-s)
    k =Int64( floor(s/(3d)))
    a=Int64(s-3d*k)
    if a==d ||a==2d || a==3d
    #    println("Bad Input: a ∈ {d, 2d, 3d}")
        return 1
    end
    if 2d+1<=a && a<=3d
        return 1//3
    end
    if 1<=a && a<= 2d
        X=min(a//2,(a+d)//4)
        return ((d*k+X)//(3d*k+a))
    end
    return 1
end

function Pebm(m,s)#print easy buddy method
    println("f(",m,",",s,") ≤ ",EBM(m,s))
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
