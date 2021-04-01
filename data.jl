struct ElectronicBasis
    electronic_state            :: Int64
    orbital_momentum_projection :: Int64
    spin_momentum               :: Float64
    spin_momentum_projection    :: Float64
end

struct VibrationalBasis
    electronic_state    :: Int64
    vibrational_number  :: Int64
end

struct RotationalBasis
    total_momentum :: Float64
    total_momentum_projection :: Float64
    total_momentum_molecular  :: Float64
end

struct RovibronicBasis
    electronic  :: ElectronicBasis
    vibrational :: VibrationalBasis
    rotational  :: RotationalBasis
end