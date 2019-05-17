function paritycheck(word,prev_29,prev_30)

        #Parity check to verify the data integrity:
        D_25 = prev_29 ⊻ word[30] ⊻ word[29] ⊻ word[28] ⊻ word[26] ⊻ word[25] ⊻ word[21] ⊻ word[20] ⊻ word[19] ⊻ word[18] ⊻ word[17] ⊻ word[14] ⊻ word[13] ⊻ word[11] ⊻ word[8]
        D_26 = prev_30 ⊻ word[29] ⊻ word[28] ⊻ word[27] ⊻ word[25] ⊻ word[24] ⊻ word[20] ⊻ word[19] ⊻ word[18] ⊻ word[17] ⊻ word[16] ⊻ word[13] ⊻ word[12] ⊻ word[10] ⊻ word[7]
        D_27 = prev_29 ⊻ word[30] ⊻ word[28] ⊻ word[27] ⊻ word[26] ⊻ word[24] ⊻ word[23] ⊻ word[19] ⊻ word[18] ⊻ word[17] ⊻ word[16] ⊻ word[15] ⊻ word[12] ⊻ word[11] ⊻ word[9]
        D_28 = prev_30 ⊻ word[29] ⊻ word[27] ⊻ word[26] ⊻ word[25] ⊻ word[23] ⊻ word[22] ⊻ word[18] ⊻ word[17] ⊻ word[16] ⊻ word[15] ⊻ word[14] ⊻ word[11] ⊻ word[10] ⊻ word[8]
        D_29 = prev_30 ⊻ word[30] ⊻ word[28] ⊻ word[26] ⊻ word[25] ⊻ word[24] ⊻ word[22] ⊻ word[21] ⊻ word[17] ⊻ word[16] ⊻ word[15] ⊻ word[14] ⊻ word[13] ⊻ word[10] ⊻ word[9] ⊻ word[7]
        D_30 = prev_29 ⊻ word[28] ⊻ word[26] ⊻ word[25] ⊻ word[23] ⊻ word[22] ⊻ word[21] ⊻ word[20] ⊻ word[18] ⊻ word[16] ⊻ word[12] ⊻ word[9] ⊻ word[8] ⊻ word[7]

        computed_parity_bits = [D_30; D_29; D_28; D_27; D_26; D_25]
        data_integrity = computed_parity_bits == word[1:6]
        return data_integrity
end
