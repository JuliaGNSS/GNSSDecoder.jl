

"""
    #checks for parity errors. 
    $(SIGNATURES)

    # Details
    # Uses the parity algorithm to compute parity and compare them to given parity bits
"""
function parity_check(word, prev_29, prev_30)
    
        # Parity check to verify the data integrity:
    D_25 = prev_29 ⊻ word[1] ⊻ word[2] ⊻ word[3] ⊻ word[5] ⊻ word[6] ⊻ word[10] ⊻ word[11] ⊻ word[12] ⊻ word[13] ⊻ word[14] ⊻ word[17] ⊻ word[18] ⊻ word[20] ⊻ word[23]
    D_26 = prev_30 ⊻ word[2] ⊻ word[3] ⊻ word[4] ⊻ word[6] ⊻ word[7] ⊻ word[11] ⊻ word[12] ⊻ word[13] ⊻ word[14] ⊻ word[15] ⊻ word[18] ⊻ word[19] ⊻ word[21] ⊻ word[24]
    D_27 = prev_29 ⊻ word[1] ⊻ word[3] ⊻ word[4] ⊻ word[5] ⊻ word[7] ⊻ word[8] ⊻ word[12] ⊻ word[13] ⊻ word[14] ⊻ word[15] ⊻ word[16] ⊻ word[19] ⊻ word[20] ⊻ word[22]
    D_28 = prev_30 ⊻ word[2] ⊻ word[4] ⊻ word[5] ⊻ word[6] ⊻ word[8] ⊻ word[9] ⊻ word[13] ⊻ word[14] ⊻ word[15] ⊻ word[16] ⊻ word[17] ⊻ word[20] ⊻ word[21] ⊻ word[23]
    D_29 = prev_30 ⊻ word[1] ⊻ word[3] ⊻ word[5] ⊻ word[6] ⊻ word[7] ⊻ word[9] ⊻ word[10] ⊻ word[14] ⊻ word[15] ⊻ word[16] ⊻ word[17] ⊻ word[18] ⊻ word[21] ⊻ word[22] ⊻ word[24]
    D_30 = prev_29 ⊻ word[3] ⊻ word[5] ⊻ word[6] ⊻ word[8] ⊻ word[9] ⊻ word[10] ⊻ word[11] ⊻ word[13] ⊻ word[15] ⊻ word[19] ⊻ word[22] ⊻ word[23] ⊻ word[24]

    computed_parity_bits = [D_25, D_26, D_27, D_28, D_29, D_30]
    data_integrity = computed_parity_bits == word[25:30]
    return data_integrity
end