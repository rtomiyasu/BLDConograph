sampleName = 'M402-pf-119673'

# (a, b, c, alpha, beta, gamma)
latticeParameters = [1.0654e+001, 2.2895e+001, 2.3003e+001, 6.0600e+001, 8.9377e+001, 8.9119e+001]

# resolution
resolution = 0.1

# DoesPrudentSymmetrySearch 0:quick search, 1:prudent search
doesPrudentSymmetrySearch = 1

# Bravais types to execute Bravais lattice determination (0:No, 1:Yes)
flgsDetermination = {
    'Triclinic' : 1, 'MonoclinicP' : 1,
    'MonoclinicB' : 1, 'OrthorhombicP' : 1,
    'OrthorhombicB' : 1, 'OrthorhombicI' : 1,
    'OrthorhombicF' : 1, 'TetragonalP' : 1,
    'TetragonalI' : 1, 'Rhombohedral' : 1,
    'Hexagonal' : 1, 'CubicP' : 1,
    'CubicI' : 1, 'CubicF' : 1
                                                                }

# AxisForRhombohedralSymmetry <"Rhombohedral" or "Hexagonal"  -->
axisForRhombohedralSymmetry = 'Hexagonal'

# AxisForBaseCenteredSymmetry <"A", "B", or "C">
axisForBaseCenteredSymmetry = 'B'

# epsilon
epsilon = resolution * 2