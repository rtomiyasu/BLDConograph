sampleName = 'CeO2_matsushita'

# Parameters a,b,c,alpha,beta,gamma of primitive unit-cell (not necessary to be Niggli-reduced)
latticeParameters = [3.7250e+000, 3.7269e+000, 3.7286e+000, 9.0100e+001, 1.1996e+002, 1.1991e+002]

# This parameter is used to judge if two unit-cell parameters are equivalent or not, by checking if 
# their corresponding metric tensors (sij) and (tij) have a relative difference that satisfies |sij-tij| <= epsilon*max{sij, tij}.
epsilon = 0.2

# 0: quick search, 1: exhaustive search
doesPrudentSymmetrySearch = 1

# Axis for rhombohedral cells: <"Rhombohedral" or "Hexagonal"
axisForRhombohedralSymmetry = 'Hexagonal'

# Unique axis for base-centered monoclinic cells: "A", "B", or "C"
axisForBaseCenteredSymmetry = 'B'
