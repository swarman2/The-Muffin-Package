include("helper_functions.jl")
using JuMP
using GLPK

"""
one_thrd is

INPUTS:
 * m -> number of muffins
 * s -> number of students
 * proof -> boolean:
   (0 - No Proof)
   (1 - Print Proof)

OUTPUTS:
 *
 * Print proof (Depending on input (proof))

"""

function one_thrd(m,s, proof = 0)
    if proof >=1
        println("f(",m,", ",s,") ≥ 1/3 proof")
    end
    if m%s == 0
        if proof >=1
            println("s divides m, give all students whole muffin(s)")
        end
      return 1
    end
    fixed_muffin = Int64(m - s)
    if proof >=2
        println("Divide ",fixed_muffin," muffins into (1/3, 1/3, 1/3)\n")
    end
    denom=lcm(s,3)
    if proof >=1
        println("Numbers assumed to have denominator: ",denom)
    end
    new_muffin = m - fixed_muffin

    pieces = Int64(fixed_muffin * 3)
    floor_pieces = Int64(floor(pieces/s))
    ceil_pieces  = Int64(ceil(pieces/s))

    mat_1 = [floor_pieces ceil_pieces]
    row_1, col_1 = size(mat_1)

    a = (ones(Int64,col_1))'
    final_mat = [mat_1; a]

    mat_3 = [pieces; s]

    model = Model(with_optimizer(GLPK.Optimizer))
      @variable(model, x[i=1:col_1],Int)
      @constraint(model,con_1,x.>=0)
      @constraint(model,con_2,final_mat*x .== mat_3)
      optimize!(model)

      if(!has_values(model))
          if proof >=1
              println("No procedure found")
          end
          return false
      end

    student_share_1 = Array{Any}(undef,Int64(floor_pieces))
    student_share_2 = Array{Any}(undef,Int64(ceil_pieces))

    for k = 1:floor_pieces
        fill!(student_share_1, Rational(1//3))
    end

    for g = 1:ceil_pieces
        fill!(student_share_2, Rational(1//3))
    end

    type1stud=Vector{Vector{Rational}}(undef,0)
    type2stud=Vector{Vector{Rational}}(undef,0)

    for i = 1:1:value.(x)[1]
        if proof >=2
            print("Give student ",Int64(i)," -> (  ")
            for l = 1:length(student_share_1)
                print(Float64(student_share_1[l]*denom),"  ")
            end
            println(")")
        end
        push!(type1stud,student_share_1)
    end

    for j = 1:1:value.(x)[2]
        if proof >=2
            print("Give student ",Int64((value.(x)[1]+j))," -> (  ")
            for l = 1:length(student_share_2)
                print(Float64(student_share_2[l]*denom),"  ")
            end
            println(")")
        end
        push!(type2stud,student_share_2)
    end

    new_muffin = new_muffin - 1
    if proof >=2
        println("\nCut 1 muffin -> (1/2, 1/2)\n")
    end

    x1 = Rational(Rational(m//s) - floor_pieces * 1//3)
    x2 = Rational(Rational(m//s) - ceil_pieces * 1//3)

    x1_half = Rational(x1 - 1//2)
    x1_half_dis = Rational(abs(x1_half - 1//2))

    x2_half = Rational(x2 - 1//2)
    x2_half_dis = Rational(abs(x2_half - 1//2))

    half_to_W = false

    if x1_half_dis < x2_half_dis
        if proof >=2
            println("Give a ",Int64((floor_pieces+2)),"-student the piece of size (1/2)")
        end
   elseif x1_half_dis > x2_half_dis   ## true for 11/5
        half_to_W = true
        if proof >=2
            println("Give a ",Int64((ceil_pieces+2)),"-student the piece of size (1/2)")
        end
    end

    if half_to_W  # 11/5
       if length(type2stud)>=2
          push!(type2stud[1],1//2)
          push!(type2stud[2],1//2)
          missing_piece = m//s-sum(type2stud[1])
       elseif (length(type2stud)==1)
          push!(type2stud[1],1//2)
          push!(type1stud[1],1//2)
          missing_piece = m//s-sum(type2stud[1])
       else
          println("???")
       end
    else
       if length(type1stud)>=2
          push!(type1stud[1],1//2)
          push!(type1stud[2],1//2)
          missing_piece = m//s-sum(type1stud[1])
       elseif (length(type1stud)==1)
          push!(type2stud[1],1//2)
          push!(type1stud[1],1//2)
          missing_piece = m//s-sum(type1stud[1])
       else
          println("???")
       end
    end
    if proof >=2
        for i = 1:length(type1stud)
            print("Give student ",Int64(i)," -> (  ")
            for l = 1:length(type1stud[i])
                print(Float64(type1stud[i][l]*denom),"  ")
            end
            println(")")
        end

        for j = 1:length(type2stud)
            print("Give student ",Int64((j+length(type1stud)))," -> (  ")
            for l = 1:length(type2stud[j])
                print(Float64(type2stud[j][l]*denom),"  ")
            end
            println(")")
        end
    end
    full = false
    muffin_count = m-2
    while !full
        for i=1:length(type1stud)
            if sum(type1stud[i])==m//s - missing_piece
                push!(type1stud[i],missing_piece)
            end
        end
        for i=1:length(type2stud)
            if sum(type2stud[i])==m//s - missing_piece
                push!(type2stud[i],missing_piece)
            end
        end
        x_half = 1
        new_piece = 1-missing_piece
        for i = 1:length(type1stud)
            if sum(type1stud[i])!=m//s
                x_half = abs(1//2 - (m//s - sum(type1stud[i])-new_piece))
                break
            end
        end
        y_half =1
        for i = 1:length(type2stud)
            if sum(type2stud[i])!=m//s
                y_half = abs(1//2 - (m//s - sum(type2stud[i])-new_piece))
                break
            end
        end
        if x_half < y_half
            for i=1:length(type1stud)
                if sum(type1stud[i])!=m//s
                    if length(type1stud)-i >=1
                        push!(type1stud[i],new_piece)
                        push!(type1stud[i+1],new_piece)
                        missing_piece = m//s - sum(type1stud[i])
                        break
                    else
                        push!(type1stud[i],new_piece)
                        for j=1:length(type2stud)
                            if sum(type2stud[j])!=m//s
                                push!(type2stud[j],new_piece)
                            end
                        end
                        missing_piece = m//s - sum(type1stud[i])
                        break
                    end
                end
            end
        else
            for i=1:length(type2stud)
                if sum(type2stud[i])!=m//s
                    if length(type2stud)-i >=1
                        push!(type2stud[i],new_piece)
                        push!(type2stud[i+1],new_piece)
                        missing_piece = m//s - sum(type2stud[i])
                        break
                    else
                        push!(type2stud[i],new_piece)
                        for j=1:length(type1stud)
                            if sum(type1stud[j])!=m//s
                                push!(type1stud[j],new_piece)
                            end
                        end
                        missing_piece = m//s - sum(type2stud[i])
                        break
                    end
                end
            end
        end
        muffin_count = muffin_count -2
        if proof >=2
            if muffin_count >0
                println("\nCut 2 muffins -> (",numerator(new_piece),"/",denominator(new_piece),", ", numerator(1-new_piece),"/",denominator(1-new_piece),")\n")
            else
                println("\nCut 1 muffin -> (",numerator(new_piece),"/",denominator(new_piece),", ", numerator(1-new_piece),"/",denominator(1-new_piece),")\n")
            end
        end

        if proof >=2
            for i = 1:length(type1stud)
                print("Give student ",Int64(i)," -> (  ")
                for l = 1:length(type1stud[i])
                    print(Float64(type1stud[i][l]*denom),"  ")
                end
                println(")")
            end


            for j = 1:length(type2stud)
                print("Give student ",Int64((j+length(type1stud)))," -> (  ")
                for l = 1:length(type2stud[j])
                    print(Float64(type2stud[j][l]*denom),"  ")
                end
                println(")")
            end

            println()
        end
        full = true
        for i=1:length(type2stud)
            if length(type2stud[i])< ceil_pieces +2
                full = false
            end
        end
        for i=1:length(type1stud)
            if length(type1stud[i])< floor_pieces +2
                full = false
            end
        end
    end
    valid = true
    if proof ==1
        distrib = Array{Any}(undef,0)
        while  length(type1stud)>0
            k = type1stud[1]
            num = 0
            for i = 1:length(type1stud)
                if type1stud[i] == k
                    num = num+1
                end
            end
        #    num = count(i -> type1stud[i] == k)
            append!(distrib, num)
            append!(distrib,[k])
            type1stud = filter(i -> i!= k, type1stud)
        end
        i=1
        while i < length(distrib)
            print("Give ", distrib[i]," student(s)  -> (  ")
            for l = 1:length(distrib[i+1])
                print(Float64(distrib[i+1][l]*denom),"  ")
            end
            println(")")
            i = i+2
        end
        println()
        distrib = Array{Any}(undef,0)
        while  length(type2stud)>0
            k = type2stud[1]
            num = 0
            for i = 1:length(type2stud)
                if type2stud[i] == k
                    num = num+1
                end
            end
        #    num = count(i -> type1stud[i] == k)
            append!(distrib, num)
            append!(distrib,[k])
            type2stud = filter(i -> i!= k, type2stud)
        end
        i=1
        while i < length(distrib)
            print("Give ", distrib[i]," student(s)  -> (  ")
            for l = 1:length(distrib[i+1])
                print(Float64(distrib[i+1][l]*denom),"  ")
            end
            println(")")
            i = i+2
        end
        println()
    end
    for i=1:length(type1stud)
        if sum(type1stud[i]) != m//s
            valid = false
            if proof >=1
                println(type1stud[i])
                println("Sum of student ",i," = ",sum(type1stud[i])," ≠ ",m//s)
            end
        end
    end
    for i=1:length(type2stud)
        if sum(type2stud[i]) != m//s
            valid = false
            if proof >=1
                println(type2stud[i])
                println("Sum of student ",i+length(type1stud)," = ",sum(type2stud[i])," ≠ ",m//s)
            end
        end
    end
    if valid
        if proof >=1
            println("All students have an equal amount of muffins")
        end
    end
    return valid
end
