include("PKG.jl")
using Printf
function menu()
    @printf("\t MENU: enter 1, 2 or 3")
    @printf("\n1. Giving m and s what is the largest possible smallest piece")
    @printf("\n2. Giving a range of m and s find all largest possible smallest pieces")
    @printf("\n3. Quit\n")
    answer = readline()
    while !(answer == "1" || answer == "2" || answer =="3" )
        println("Enter a number between 1 and 3")
        answer = readline()
    end
    return answer
end
function main()
    answer = menu()
    while answer != "3"
        if answer == "1"

            @printf("Enter m: ")
            m = readline()
            while true
                try
                    m = parse(Int64,m)
                    break;
                catch
                    @printf("Enter an integer for m: ")
                    m=readline()
                end
            end
            #m=parse(Int64, m)
            @printf("Enter s: ")
            s = readline()
            while true
                 try
                    s = parse(Int64,s)
                    break;
                catch
                    @printf("Enter an integer for s: ")
                    s=readline()
                end
            end

            alpha,string,str_eq,time = FIND_ALPHA(m,s)
            if str_eq == " TIME OUT "
                str_eq = " ≤ "
            end
            #@printf("\nf(%i,%i) %s %i/%i | Methods: %s  ", m,s,str_eq, numerator(alpha), denominator(alpha), string)
            @printf("\nEnter which method for proof, or 'none': ")
            method = readline()
            if method != "none"
                print_proof(method,m,s,alpha)
                println()
                VProc(m,s,alpha, Inf, Inf, 0, 1)
            end


        end
        if answer == "2"
            println("Would you like to store the results in a csv file? (yes, no)")
            file = readline()
            @printf("Enter max m: ")
            m = readline()
            while true
                try
                    m = parse(Int64,m)
                    break;
                catch
                    @printf("Enter an integer for m: ")
                    m=readline()
                end
            end
            #m=parse(Int64, m)
            @printf("Enter min s: ")
            min_s = readline()
            while true
                 try
                    min_s = parse(Int64,min_s)
                    break;
                catch
                    @printf("Enter an integer for min s: ")
                    min_s=readline()
                end
            end
            while min_s<=1
                @printf("Enter min s ≥ 1: ")
                min_s = readline()
                while true
                     try
                        min_s = parse(Int64,min_s)
                        break;
                    catch
                        @printf("Enter an integer for min s: ")
                        min_s=readline()
                    end
                end
            end
            @printf("Enter max s: ")
            max_s = readline()
            while true
                 try
                    max_s = parse(Int64,max_s)
                    break;
                catch
                    @printf("Enter an integer for max s: ")
                    max_s=readline()
                end
            end
            while max_s<min_s
                @printf("Enter max s > min s: ")
                min_s = readline()
                while true
                     try
                        min_s = parse(Int64,min_s)
                        break;
                    catch
                        @printf("Enter an integer for min s: ")
                        min_s=readline()
                    end
                end
            end
            if file == "yes"
                println("The file will be saved as :",dirname(@__FILE__)*"/Muffins_m"*string(m)*"_s"*string(max_s)*"_output.csv")
                FIND_ALL(m,min_s, max_s, true)
            else
                FIND_ALL(m,min_s, max_s, false)
            end
        end
        answer=menu()
    end
end

function print_proof(method, m, s, alpha)
    if method == "FC"
        println("\nFC proof of upperbound")
        print("***********************************************")
        FC(m,s,true)
        println("***********************************************")
    elseif method == "INT"
        println("\nINT proof of upperbound")
        print("***********************************************")
        VINT(m,s,alpha,true)
        println("***********************************************")
    elseif method == "HALF"
        println("\nHALF proof of upperbound")
        println("***********************************************")
        Half_proof(m,s,alpha)
        println("***********************************************")
    elseif method == "MID"
        println("\nMID proof of upperbound")
        println("***********************************************")
        VMID(m,s,alpha,true)
        println("***********************************************")
    elseif method == "EBM"
        println("\nEBM proof of upperbound")
        println("***********************************************")
        EBM(m,s,true)
        println("***********************************************")
    elseif method == "HBM"
        println("\nHBM proof of upperbound")
        println("***********************************************")
        HBM(m,s,true)
        println("***********************************************")
    elseif method == "GAP"
        println("\nGAP proof of upperbound")
        println("***********************************************")
        VGAPV3(m,s,alpha,1)
        println("***********************************************")
    end
end
main()
