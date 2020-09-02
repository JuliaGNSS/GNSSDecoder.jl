PREAMBLE = BitArray([1,1,0,1,0,0,0,1])

function make_subframe_ninv(id_data)
    return BitArray(vcat(falses(240),falses(8), reverse(id_data), falses(19), falses(22), PREAMBLE))
end
function make_subframe_inv(id_data)
    return BitArray(vcat(falses(240),falses(8), reverse(id_data), falses(19), falses(22), .!PREAMBLE))
end

test_buffer_ninv = Vector{BitArray}(undef, 0)
push!(test_buffer_ninv, make_subframe_ninv([0,0,1]))
push!(test_buffer_ninv, make_subframe_ninv([0,1,0]))
push!(test_buffer_ninv, make_subframe_ninv([0,1,1]))
push!(test_buffer_ninv, make_subframe_ninv([1,0,0]))
push!(test_buffer_ninv, make_subframe_ninv([1,0,1]))

test_buffer_inv = Vector{BitArray}(undef, 0)
push!(test_buffer_inv, make_subframe_inv([0,0,1]))
push!(test_buffer_inv, make_subframe_inv([0,1,0]))
push!(test_buffer_inv, make_subframe_inv([0,1,1]))
push!(test_buffer_inv, make_subframe_inv([1,0,0]))
push!(test_buffer_inv, make_subframe_inv([1,0,1]))

false_frame = zeros(300)

@testset "Test Not Inverted Preamble" begin
    buffer = vcat(test_buffer_ninv[5], test_buffer_ninv[4], test_buffer_ninv[3], test_buffer_ninv[2], test_buffer_ninv[1])
    out = GNSSDecoder.find_preamble(buffer)
    @test out == true


    buffer = vcat(false_frame, test_buffer_ninv[4], test_buffer_ninv[3], test_buffer_ninv[2], test_buffer_ninv[1])
    out = GNSSDecoder.find_preamble(buffer)
    @test out == false

    buffer = vcat(test_buffer_ninv[1], test_buffer_ninv[2], test_buffer_ninv[3], test_buffer_ninv[4], test_buffer_ninv[5])
    out = GNSSDecoder.find_preamble(buffer)
    @test out == false

end

@testset "Test Inverted Preamble" begin
    buffer = vcat(test_buffer_inv[5], test_buffer_inv[4], test_buffer_inv[3], test_buffer_inv[2], test_buffer_inv[1])
    out = GNSSDecoder.find_preamble(buffer)
    @test out == true


    buffer = vcat(false_frame, test_buffer_inv[4], test_buffer_inv[3], test_buffer_inv[2], test_buffer_inv[1])
    out = GNSSDecoder.find_preamble(buffer)
    @test out == false

    buffer = vcat(test_buffer_inv[1], test_buffer_inv[2], test_buffer_inv[3], test_buffer_inv[4], test_buffer_inv[5])
    out = GNSSDecoder.find_preamble(buffer)
    @test out == false
end

@testset "Test Mixed Preamble" begin
    buffer = vcat(test_buffer_inv[5], test_buffer_ninv[4], test_buffer_inv[3], test_buffer_ninv[2], test_buffer_inv[1])
    out = GNSSDecoder.find_preamble(buffer)
    @test out == true


    buffer = vcat(false_frame, test_buffer_ninv[4], test_buffer_ninv[3], test_buffer_inv[2], test_buffer_inv[1])
    out = GNSSDecoder.find_preamble(buffer)
    @test out == false

    buffer = vcat(test_buffer_ninv[1], test_buffer_inv[2], test_buffer_inv[3], test_buffer_inv[4], test_buffer_ninv[5])
    out = GNSSDecoder.find_preamble(buffer)
    @test out == false
end