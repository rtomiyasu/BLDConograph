sampleName = 'sample3'

# Parameters a,b,c,alpha,beta,gamma of primitive unit-cell (not necessary to be Niggli-reduced)
latticeParameters = [1.6332e+001, 1.6348e+001, 1.6357e+001, 1.0941e+002, 1.0952e+002, 1.0934e+002]

# This parameter is used to judge if two unit-cell parameters are equivalent or not, by checking if 
# their corresponding metric tensors (sij) and (tij) have a relative difference that satisfies |sij-tij| <= epsilon*max{sij, tij}.
epsilon = 0.2

# 0: quick search, 1: exhaustive search
doesPrudentSymmetrySearch = 1

# Axis for rhombohedral cells: <"Rhombohedral" or "Hexagonal"
axisForRhombohedralSymmetry = 'Hexagonal'

# Unique axis for base-centered monoclinic cells: "A", "B", or "C"
axisForBaseCenteredSymmetry = 'B'
