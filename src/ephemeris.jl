# The broadcast ephemeris and clock, normalised across every message that carries
# one: the per-message decoded fields keep their ICD names (and RINEX writers read
# them), while a propagator reads this one record.

"""
$(TYPEDEF)

The broadcast ephemeris and satellite clock of one decoder, normalised across
the messages that carry them. [`get_broadcast_ephemeris`](@ref) returns it.

The messages state the orbit in two forms. GPS LNAV, Galileo I/NAV and F/NAV and
BeiDou D1/D2 broadcast Keplerian elements: `√A`, a constant mean-motion
correction `Δn` and the absolute `Ω̇`. GPS CNAV and CNAV-2 and BeiDou B-CNAV1-3
broadcast quasi-Keplerian ones: `ΔA` off a reference semi-major axis, its rate
`Ȧ`, `Δn_0` with its rate, and on GPS `ΔΩ̇` off a reference `Ω̇_REF`. This record
holds the quasi-Keplerian form, which contains the other, so one propagator
covers every message: the Keplerian messages report `A = √A²` and zero rates.

At `t_k = t - t_0e` seconds from the reference epoch, folded across a week
crossover, the elements are

```julia
A_k = A + A_dot * t_k
n = sqrt(μ / A^3) + Δn_0 + Δn_0_dot * t_k / 2
M_k = M_0 + n * t_k
```

and the rest follows the user algorithm every ICD here shares (IS-GPS-200
Table 20-IV and §30.3.3.1.3, Galileo OS SIS ICD Table 58, BDS-SIS-ICD-B1I-3.0
§5.2.4.12, BDS-SIS-ICD-B1C-1.0 §7.8.2), with `Ω_dot_e` the Earth rotation rate
of the system's reference frame. The BeiDou GEO satellites
(`is_geostationary`, from [`get_orbit_class`](@ref)) have a branch of their
own: D1/D2 propagates them in an inertial frame tilted 5° off the equator and
then rotates them into BDCS (BDS-SIS-ICD-B1I-3.0 §5.2.4.12).

The clock polynomial is `Δt_SV = a_f0 + a_f1 * (t - t_0c) + a_f2 * (t - t_0c)^2`,
before the relativistic term `F * e * √A * sin(E_k)` and the group delay of the
ranging signal ([`get_group_delay`](@ref)).

Every field is a `Float64` (the flag a `Bool`), so the record is a bits type a
navigation filter can copy, store and compare with `==` to notice a new issue
of data. Angles are in radians, as the decoded fields are.

# Fields

$(TYPEDFIELDS)

# See Also

  - [`get_broadcast_ephemeris`](@ref): the accessor that produces one
"""
struct BroadcastEphemeris
    """
    Reference time of ephemeris (s of week)
    """
    t_0e::Float64
    """
    Semi-major axis at `t_0e` (m)
    """
    A::Float64
    """
    Rate of the semi-major axis (m/s); zero on the Keplerian messages
    """
    A_dot::Float64
    """
    Mean motion difference from the computed value at `t_0e` (rad/s)
    """
    Δn_0::Float64
    """
    Rate of the mean motion difference (rad/s²); zero on the Keplerian messages
    """
    Δn_0_dot::Float64
    """
    Mean anomaly at `t_0e` (rad)
    """
    M_0::Float64
    """
    Eccentricity
    """
    e::Float64
    """
    Argument of perigee (rad)
    """
    ω::Float64
    """
    Longitude of the ascending node at the start of the week (rad)
    """
    Ω_0::Float64
    """
    Rate of right ascension (rad/s), absolute: `Ω̇_REF + ΔΩ̇` on GPS CNAV and CNAV-2
    """
    Ω_dot::Float64
    """
    Inclination at `t_0e` (rad)
    """
    i_0::Float64
    """
    Rate of inclination (rad/s)
    """
    i_dot::Float64
    """
    Sine harmonic correction to the orbit radius (m)
    """
    C_rs::Float64
    """
    Cosine harmonic correction to the orbit radius (m)
    """
    C_rc::Float64
    """
    Sine harmonic correction to the argument of latitude (rad)
    """
    C_us::Float64
    """
    Cosine harmonic correction to the argument of latitude (rad)
    """
    C_uc::Float64
    """
    Sine harmonic correction to the inclination (rad)
    """
    C_is::Float64
    """
    Cosine harmonic correction to the inclination (rad)
    """
    C_ic::Float64
    """
    Reference time of the clock (s of week)
    """
    t_0c::Float64
    """
    Clock bias (s)
    """
    a_f0::Float64
    """
    Clock drift (s/s)
    """
    a_f1::Float64
    """
    Clock drift rate (s/s²)
    """
    a_f2::Float64
    """
    Gravitational parameter of the system's reference frame (m³/s²)
    """
    μ::Float64
    """
    Relativistic correction constant (s/√m)
    """
    F::Float64
    """
    Earth rotation rate of the system's reference frame (rad/s)
    """
    Ω_dot_e::Float64
    """
    Whether the satellite is a BeiDou GEO, which the ICD propagates differently
    """
    is_geostationary::Bool
end

"""
    get_broadcast_ephemeris(state::GNSSDecoderState) -> Union{Nothing,BroadcastEphemeris}

The decoder's validated ephemeris and satellite clock as a
[`BroadcastEphemeris`](@ref), or `nothing` until
[`is_decoding_completed_for_positioning`](@ref). Galileo E6-B, which carries
corrections rather than an ephemeris, answers `nothing` always.

It is built from `data` on each call, so it always reflects the latest
validated issue, and allocates nothing.

# Example

```julia
ephemeris = get_broadcast_ephemeris(state)
if !isnothing(ephemeris) && is_satellite_healthy(state)
    t_k = t - ephemeris.t_0e   # fold across a week crossover
    n = sqrt(ephemeris.μ / ephemeris.A^3) + ephemeris.Δn_0 + ephemeris.Δn_0_dot * t_k / 2
    # …
end
```
"""
function get_broadcast_ephemeris(state::GNSSDecoderState)
    is_decoding_completed_for_positioning(state) || return nothing
    broadcast_ephemeris(
        state.data,
        frame_constants(state.data),
        get_orbit_class(state) === geostationary_orbit,
    )
end

# The constants of the reference frame each constellation states its orbits in: the
# gravitational parameter μ, the relativistic correction F and the Earth rotation
# rate (WGS-84 for GPS, GTRF for Galileo, BDCS for BeiDou).
frame_constants(::AbstractGPSData) = (μ = GPS_μ, F = GPS_F, Ω_dot_e = EARTH_ROTATION_RATE)
frame_constants(::AbstractGalileoData) =
    (μ = GALILEO_μ, F = GALILEO_F, Ω_dot_e = EARTH_ROTATION_RATE)
frame_constants(::AbstractBeiDouData) =
    (μ = BEIDOU_μ, F = BEIDOU_F, Ω_dot_e = BEIDOU_EARTH_ROTATION_RATE)

# The Keplerian messages: `√A`, a constant `Δn` and the absolute `Ω̇`.
function broadcast_ephemeris(
    data::Union{GPSL1CAData,AbstractGalileoEphemerisData,BeiDouDNAVData},
    frame,
    is_geostationary,
)
    sqrt_A = something(data.sqrt_A)
    broadcast_ephemeris(
        data,
        frame,
        is_geostationary,
        sqrt_A^2,
        0.0,
        something(data.Δn),
        0.0,
        something(data.Ω_dot),
    )
end

# GPS CNAV and CNAV-2: `ΔA` and `ΔΩ̇` off fixed reference values (IS-GPS-200N
# §30.3.3.1.3, as IS-GPS-705J and IS-GPS-800J).
const GPS_CNAV_A_REF = 26_559_710.0          # m
const GPS_CNAV_Ω_DOT_REF = -2.6e-9 * GNSS_PI  # rad/s (−2.6e-9 semicircles/s)

function broadcast_ephemeris(data::AbstractGPSCNAVData, frame, is_geostationary)
    broadcast_ephemeris(
        data,
        frame,
        is_geostationary,
        GPS_CNAV_A_REF + something(data.ΔA),
        something(data.A_dot),
        something(data.Δn_0),
        something(data.Δn_0_dot),
        GPS_CNAV_Ω_DOT_REF + something(data.ΔΩ_dot),
    )
end

# BeiDou B-CNAV: `ΔA` off a reference that depends on the orbit class, and `Ω̇`
# outright (BDS-SIS-ICD-B2a-1.0 Table 7-8: 19 bits at 2⁻⁴⁴ π/s), not as a delta.
const BEIDOU_CNAV_A_REF_MEO = 27_906_100.0        # m
const BEIDOU_CNAV_A_REF_IGSO_GEO = 42_162_200.0   # m

# BDS-SIS-ICD-B1C-1.0 Table 7-3 and the B2a/B2b counterparts: 3 MEO, 1 GEO and
# 2 IGSO. The reserved 0 leaves the orbit undefined; the readiness gate keeps it
# out (`is_known_sat_type`), and it throws here rather than pass for an IGSO.
function beidou_cnav_reference_semi_major_axis(sat_type)
    sat_type == 3 && return BEIDOU_CNAV_A_REF_MEO
    (sat_type == 1 || sat_type == 2) && return BEIDOU_CNAV_A_REF_IGSO_GEO
    throw(
        ArgumentError(
            "BeiDou sat_type $sat_type is reserved: no reference semi-major axis",
        ),
    )
end

function broadcast_ephemeris(data::AbstractBeiDouCNAVData, frame, is_geostationary)
    broadcast_ephemeris(
        data,
        frame,
        is_geostationary,
        beidou_cnav_reference_semi_major_axis(something(data.sat_type)) +
        something(data.ΔA),
        something(data.A_dot),
        something(data.Δn_0),
        something(data.Δn_0_dot),
        something(data.Ω_dot),
    )
end

# The elements every message names alike, around the ones each states its own
# way. Positional, so that every argument reaches it concretely typed.
function broadcast_ephemeris(data, frame, is_geostationary, A, A_dot, Δn_0, Δn_0_dot, Ω_dot)
    BroadcastEphemeris(
        something(data.t_0e),
        A,
        A_dot,
        Δn_0,
        Δn_0_dot,
        something(data.M_0),
        something(data.e),
        something(data.ω),
        something(data.Ω_0),
        Ω_dot,
        something(data.i_0),
        something(data.i_dot),
        something(data.C_rs),
        something(data.C_rc),
        something(data.C_us),
        something(data.C_uc),
        something(data.C_is),
        something(data.C_ic),
        something(data.t_0c),
        something(data.a_f0),
        something(data.a_f1),
        something(data.a_f2),
        frame.μ,
        frame.F,
        frame.Ω_dot_e,
        is_geostationary,
    )
end
