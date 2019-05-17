
function find_preamble(buffer, preamble_pos)

   #preamble_mask = [falses(length(buffer) - 308); trues(8); falses(292); trues(8)]
   preamble_mask = [falses(length(buffer) - 308 - preamble_pos); trues(8); falses(292); trues(8); falses(preamble_pos)]

   masked_buffer = map(&, buffer, preamble_mask)

   preamble = [true; true; false; true; false; false; false; true]

   inverted_preamble = .!preamble

#    direct_preamble_mask = [falses(length(buffer) - 308); preamble; falses(292); preamble]
direct_preamble_mask = [falses(length(buffer) - 308 - preamble_pos); preamble; falses(292); preamble; falses(preamble_pos)]
#    inverted_preample_mask = [falses(length(buffer) - 308); inverted_preamble; falses(292); inverted_preamble]
inverted_preample_mask = [falses(length(buffer) - 308 - preamble_pos); inverted_preamble; falses(292); inverted_preamble; falses(preamble_pos)]

   found_preamble = (masked_buffer == direct_preamble_mask) || (masked_buffer == inverted_preample_mask)
   #found_preamble equals true if either the direct or the inverted forms of the preambles are found
   found_inverted_preamble = masked_buffer == inverted_preample_mask
       if !found_preamble
           preamble_pos = preamble_pos + 1
       end
   return found_preamble,found_inverted_preamble, preamble_pos
end
