#for information on this method see chapter 4 of "The Mathematics of Muffins"
using Printf
function FC(m,s, proof = false)
    if m%s == 0
      return 1
    end
    min_val =(min(m//Int64(s*ceil(2*m/s)),1-(m//Int64(s*floor(2*m/s)))))
    if proof

      @printf("\nmax{ %i/%i, min{%i/%i * %i/%i , 1- %i/%i * %i/%i}",1,3,m,s,1,ceil(2*m/s),m,s,1,floor(2*m/s))
      @printf("\n= max{%i/%i, %i/%i}",1,3,numerator(min_val),denominator(min_val))
      @printf("\n= %i/%i",numerator(max(1//3, min_val)),denominator(max(1//3, min_val)))
    end
return (max(1//3, min_val))
end
