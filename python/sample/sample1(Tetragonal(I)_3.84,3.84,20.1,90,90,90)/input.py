sampleName = 'HERMES_Sr327_250K'

# Parameters a,b,c,alpha,beta,gamma of primitive unit-cell (not necessary to be Niggli-reduced)
latticeParameters = [3.8417e+000, 3.8524e+000,1.0389e+001, 1.0041e+002, 1.0053e+002, 9.0046e+001]

# This parameter is used to judge if two unit-cell parameters are equivalent or not, by checking if 
# their corresponding metric tensors (sij) and (tij) have a relative difference that satisfies |sij-tij| <= epsilon*max{sij, tij}.
epsilon = 0.2

# 0: quick search, 1: exhaustive search
doesPrudentSymmetrySearch = 1

# Axis for rhombohedral cells: <"Rhombohedral" or "Hexagonal"
axisForRhombohedralSymmetry = 'Hexagonal'

# Unique axis for base-centered monoclinic cells: "A", "B", or "C"
axisForBaseCenteredSymmetry = 'B'
