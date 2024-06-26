sampleName = 'Alumina676'

# Parameters a,b,c,alpha,beta,gamma of primitive unit-cell (not necessary to be Niggli-reduced)
latticeParameters = [4.7428e+000, 4.7802e+000, 5.1093e+000, 1.1689e+002, 9.0170e+001, 1.1969e+002]

# This parameter is used to judge if two unit-cell parameters are equivalent or not, by checking if 
# their corresponding metric tensors (sij) and (tij) have a relative difference that satisfies |sij-tij| <= epsilon*max{sij, tij}.
epsilon = 0.1

# 0: quick search, 1: exhaustive search
doesPrudentSymmetrySearch = 1

# Axis for rhombohedral cells: <"Rhombohedral" or "Hexagonal"
axisForRhombohedralSymmetry = 'Hexagonal'

# Unique axis for base-centered monoclinic cells: "A", "B", or "C"
axisForBaseCenteredSymmetry = 'B'
