function pseudo_range(dc::GNSSDecoderState, code_phase)
    code_time = code_phase / GNSSSignals.get_code_frequency(GPSL1)*Hz
    t_tow = dc.num_bits_buffered / GNSSSignals.get_data_frequency(GPSL1)*Hz
    c = 299792458 # Speed of light
    
    p_range = (t_tow + code_time) * c
end
