# The broadcast ephemeris and clock, normalised across every message that carries one.
# Per-message decoded fields keep their ICD names (RINEX writers read them); a propagator
# reads this one record.

"""
$(TYPEDEF)

The broadcast ephemeris and satellite clock of one decoder, normalised across the
messages that carry them. [`get_broadcast_ephemeris`](@ref) returns it.

The messages state the orbit in two forms. GPS LNAV, Galileo I/NAV and F/NAV and BeiDou
D1/D2 broadcast Keplerian elements: `√A`, a constant mean-motion correction `Δn` and the
absolute `Ω̇`. GPS CNAV and CNAV-2 and BeiDou B-CNAV1-3 broadcast quasi-Keplerian ones:
`ΔA` off a reference semi-major axis, its rate `Ȧ`, `Δn_0` with its rate, and on GPS `ΔΩ̇`
off a reference `Ω̇_REF`. This record holds the quasi-Keplerian form, which contains the
other, so one propagator covers every message: Keplerian messages report `A = (√A)²`
and zero rates.

At `t_k = t - t_0e` seconds from the reference epoch, folded across a week crossover
(GNSSSignals' `fold_week_crossover`), the elements are

```julia
A_k = A + A_dot * t_k
n = sqrt(μ / A^3) + Δn_0 + Δn_0_dot * t_k / 2
M_k = M_0 + n * t_k
```

and the rest follows the user algorithm every ICD here shares (IS-GPS-200 Table 20-IV and
§30.3.3.1.3, Galileo OS SIS ICD Table 68, BDS-SIS-ICD-B1I-3.0 §5.2.4.12,
BDS-SIS-ICD-B1C-1.0 §7.8.2), with `Ω_dot_e` the Earth rotation rate it fixes. BeiDou GEO
satellites (`is_geostationary`, from [`get_orbit_class`](@ref)) branch: D1/D2 propagates
them in an inertial frame tilted 5° off the equator, then rotates them into BDCS
(BDS-SIS-ICD-B1I-3.0 §5.2.4.12).

The clock polynomial is `Δt_SV = a_f0 + a_f1 * (t - t_0c) + a_f2 * (t - t_0c)^2`, before
the relativistic term `F * e * √A * sin(E_k)` and the signal's group delay
([`get_group_delay`](@ref)).

Every field is a `Float64` (the flag a `Bool`), so the record is a bits type that can be
copied and compared with `==` to notice a new issue of data. Angles are in radians.

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
    Gravitational parameter, as the constellation's ICD fixes it (m³/s²)
    """
    μ::Float64
    """
    Relativistic correction constant `F = −2√μ/c²` (s/√m)
    """
    F::Float64
    """
    Earth rotation rate, as the constellation's ICD fixes it (rad/s)
    """
    Ω_dot_e::Float64
    """
    Whether the satellite is a BeiDou GEO, which the ICD propagates differently
    """
    is_geostationary::Bool
end

"""
    get_broadcast_ephemeris(state::GNSSDecoderState) -> Union{Nothing,BroadcastEphemeris}

The decoder's validated ephemeris and satellite clock as a [`BroadcastEphemeris`](@ref),
or `nothing` until the validated data holds the set [`has_positioning_data`](@ref)
requires. Galileo E6-B (corrections, no ephemeris) always answers `nothing`.

Built from `data` on each call, so it reflects the latest validated issue; allocates
nothing.

# Example

```julia
ephemeris = get_broadcast_ephemeris(state)
if !isnothing(ephemeris) && is_sat_healthy(state)
    t_k = fold_week_crossover(t - ephemeris.t_0e)
    n = sqrt(ephemeris.μ / ephemeris.A^3) + ephemeris.Δn_0 + ephemeris.Δn_0_dot * t_k / 2
    # …
end
```
"""
function get_broadcast_ephemeris(state::GNSSDecoderState)
    has_positioning_messages(state.data) || return nothing
    broadcast_ephemeris(
        state.data,
        ephemeris_constants(state.data),
        get_orbit_class(state) === geostationary_orbit,
    )
end

# Constants each ICD fixes for its user algorithm: gravitational parameter μ, Earth
# rotation rate of its geodetic system (WGS-84 for GPS, BDCS for BeiDou; Galileo states
# its own) and the relativistic correction F, which follows from μ.
ephemeris_constants(::AbstractGPSData) =
    (μ = GPS_μ, F = GPS_F, Ω_dot_e = EARTH_ROTATION_RATE)
ephemeris_constants(::AbstractGalileoData) =
    (μ = GALILEO_μ, F = GALILEO_F, Ω_dot_e = EARTH_ROTATION_RATE)
ephemeris_constants(::AbstractBeiDouData) =
    (μ = BEIDOU_μ, F = BEIDOU_F, Ω_dot_e = BEIDOU_EARTH_ROTATION_RATE)

# The Keplerian messages: `√A`, a constant `Δn` and the absolute `Ω̇`.
function broadcast_ephemeris(
    data::Union{GPSLNAVData,AbstractGalileoEphemerisData,BeiDouDNAVData},
    ephemeris_constants,
    is_geostationary,
)
    sqrt_A = something(data.sqrt_A)
    broadcast_ephemeris(
        data,
        ephemeris_constants,
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

function broadcast_ephemeris(
    data::AbstractGPSCNAVData,
    ephemeris_constants,
    is_geostationary,
)
    broadcast_ephemeris(
        data,
        ephemeris_constants,
        is_geostationary,
        GPS_CNAV_A_REF + something(data.ΔA),
        something(data.A_dot),
        something(data.Δn_0),
        something(data.Δn_0_dot),
        GPS_CNAV_Ω_DOT_REF + something(data.ΔΩ_dot),
    )
end

# BeiDou B-CNAV: `ΔA` off an orbit-class reference, and `Ω̇` absolute, not a delta
# (BDS-SIS-ICD-B2a-1.0 Table 7-8: 19 bits at 2⁻⁴⁴ π/s).
const BEIDOU_CNAV_A_REF_MEO = 27_906_100.0        # m
const BEIDOU_CNAV_A_REF_IGSO_GEO = 42_162_200.0   # m

# `sat_type` (B1C-1.0 and B2a-1.0 Table 7-8, B2b-1.0 Table 7-6): 3 MEO, 1 GEO, 2 IGSO.
# The reserved 0 is kept out by the readiness gate (`is_known_sat_type`); it throws here
# rather than pass for an IGSO.
function beidou_cnav_reference_semi_major_axis(sat_type)
    sat_type == 3 && return BEIDOU_CNAV_A_REF_MEO
    (sat_type == 1 || sat_type == 2) && return BEIDOU_CNAV_A_REF_IGSO_GEO
    throw(
        ArgumentError(
            "BeiDou sat_type $sat_type is reserved: no reference semi-major axis",
        ),
    )
end

function broadcast_ephemeris(
    data::AbstractBeiDouCNAVData,
    ephemeris_constants,
    is_geostationary,
)
    broadcast_ephemeris(
        data,
        ephemeris_constants,
        is_geostationary,
        beidou_cnav_reference_semi_major_axis(something(data.sat_type)) +
        something(data.ΔA),
        something(data.A_dot),
        something(data.Δn_0),
        something(data.Δn_0_dot),
        something(data.Ω_dot),
    )
end

# Elements every message names alike, around the ones each states its own way. Positional,
# so every argument arrives concretely typed.
function broadcast_ephemeris(
    data,
    ephemeris_constants,
    is_geostationary,
    A,
    A_dot,
    Δn_0,
    Δn_0_dot,
    Ω_dot,
)
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
        ephemeris_constants.μ,
        ephemeris_constants.F,
        ephemeris_constants.Ω_dot_e,
        is_geostationary,
    )
end
