
function bitArray2Str(x::Union{BitArray{1},Array{Bool,1}})
    parseString(b::Bool) = b == 1 ? "1" : "0"
    number = ""
    for i = 1:1:length(x)
        number = number * parseString(x[i])
    end
    return number
end

function intArray2Str(x::Array{Int64,1})
    parseString(n::Int64) = n == 1 ? "1" : "0"
    number = ""
    for i = 1:1:length(x)
        number = number * parseString(x[i])
    end
    return number
end


function bin2dec(buf)
    num = 0
    for i in 1:length(buf)
        if buf[length(buf) - i + 1]
            num += 2^(i - 1)
        end
    end
    return num
end

function bin2dec_twoscomp(x::Union{BitArray{1},Array{Bool,1}})
    if x[1] == false
        num = bin2dec(x)
    else
        _xor = Array{Bool,1}(trues(length(x)))
        _sum1 = [falses(length(x) - 1);true]
        x = map(⊻, x, _xor) # invert Bits
        x = map(+, x, _sum1) # sum 1. At this point x is no longer a BitArray but an Array{Int64,1}

        for i = length(x):-1:2
            if x[i] == 2
                x[i] = 0
                x[i - 1] = x[i - 1] + 1
            end
        end
        x[1] == 2 ? x[1] = 0 : x[1] = x[1]

        number = intArray2Str(x)
        num = -parse(Int, number, base = 2)
    end 
    return num
end

