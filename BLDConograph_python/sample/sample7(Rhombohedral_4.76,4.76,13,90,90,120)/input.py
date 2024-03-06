sampleName = 'Alumina676'

# (a, b, c, alpha, beta, gamma)
latticeParameters = [4.7428e+000, 4.7802e+000, 5.1093e+000, 1.1689e+002, 9.0170e+001, 1.1969e+002]

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