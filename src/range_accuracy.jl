# The broadcast accuracy of a satellite's ranging signal, in metres: every message
# broadcasts it as an index, and each ICD states the metres it stands for. A
# positioning engine weights its measurements by it, and a navigation file records it.

"""
    get_range_accuracy(state::GNSSDecoderState) -> Union{Nothing,Float64}

The broadcast accuracy of the satellite's ranging signal, in metres: the nominal
value each ICD computes from the broadcast index, from the decoder's validated
data.

| Message                   | Index         | Metres                                                                                          |
|:------------------------- |:------------- |:----------------------------------------------------------------------------------------------- |
| GPS LNAV (L1 C/A)         | `URA_index`   | user range accuracy, `2^(1+N/2)` to index 6, `2^(N-2)` beyond (IS-GPS-200N §20.3.3.3.1.3)       |
| Galileo I/NAV (E1-B, E5b) | `SISA_E1_E5b` | signal-in-space accuracy of the E1-E5b pair (OS SIS ICD Table 91)                               |
| Galileo F/NAV (E5a)       | `SISA_E1_E5a` | signal-in-space accuracy of the E1-E5a pair (OS SIS ICD Table 91)                               |
| BeiDou D1/D2 (B1I, B3I)   | `URAI`        | user range accuracy, `2^(N/2+1)` below index 6, `2^(N-2)` beyond (BDS-SIS-ICD-B1I-3.0 §5.2.4.5) |

GPS rounds the nominal values of indices 1, 3 and 5 to 2.8, 5.7 and 11.3 m, as
IS-GPS-200 asks. Where the satellite broadcasts no accuracy prediction — GPS and
BeiDou index 15, Galileo's NAPA (255) and its spare indices 126-254 — the
accuracy is `Inf`: no prediction, to be weighted as none.

`nothing` until the index is decoded, and where no single figure is broadcast:
on Galileo E6-B, which broadcasts none; on GPS CNAV and CNAV-2, whose nominal URA
is the RSS of an elevation-dependent component (`ura_ed_index`) and a component
that grows with the time since `t_op` (`ura_ned0_index`, `ura_ned1_index`,
`ura_ned2_index`; IS-GPS-705J IRN-004 §20.3.3.1.1 and §20.3.3.2.4), so it depends
on the satellite's elevation and the time it is evaluated at; and on the BeiDou
B-CNAV messages, whose ICDs defer the metres of their SISAI indices to a future
update (BDS-SIS-ICD-B1C-1.0 §7.16). Their indices stay on the decoded data.
"""
get_range_accuracy(state::GNSSDecoderState) = range_accuracy(state.data)

range_accuracy(data::GPSLNAVData) = gps_ura_metres(data.URA_index)
range_accuracy(data::GalileoINAVData) = galileo_sisa_metres(data.SISA_E1_E5b)
range_accuracy(data::GalileoFNAVData) = galileo_sisa_metres(data.SISA_E1_E5a)
range_accuracy(data::BeiDouDNAVData) = beidou_ura_metres(data.URAI)
range_accuracy(::Union{GalileoE6BData,AbstractGPSCNAVData,AbstractBeiDouCNAVData}) = nothing

# The nominal URA of an LNAV index below 15 (IS-GPS-200N §20.3.3.3.1.3), with the
# roundings the ICD asks for.
function gps_nominal_ura(N::Integer)
    N == 1 && return 2.8
    N == 3 && return 5.7
    N == 5 && return 11.3
    N <= 6 ? 2.0^(1 + N / 2) : 2.0^(N - 2)
end

# LNAV's index runs 0-15, 15 being "no accuracy prediction available".
gps_ura_metres(::Nothing) = nothing
gps_ura_metres(N::Integer) = N == 15 ? Inf : gps_nominal_ura(N)

# BDS-SIS-ICD-B1I-3.0 Table 5-4: index 15 is a manoeuvre or no accuracy prediction.
beidou_ura_metres(::Nothing) = nothing
beidou_ura_metres(N::Integer) = N == 15 ? Inf : N < 6 ? 2.0^(N / 2 + 1) : 2.0^(N - 2)

# Galileo OS SIS ICD Table 91: four linear segments to 6 m, then spare indices and
# NAPA (255), none of which is a prediction.
galileo_sisa_metres(::Nothing) = nothing
function galileo_sisa_metres(index::Integer)
    index <= 49 && return 0.01 * index
    index <= 74 && return 0.5 + 0.02 * (index - 50)
    index <= 99 && return 1.0 + 0.04 * (index - 75)
    index <= 125 && return 2.0 + 0.16 * (index - 100)
    return Inf
end
