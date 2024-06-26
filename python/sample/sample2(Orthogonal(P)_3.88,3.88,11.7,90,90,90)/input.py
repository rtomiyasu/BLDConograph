sampleName = '1339prbasrcuo_bs1'

# Parameters a,b,c,alpha,beta,gamma of primitive unit-cell (not necessary to be Niggli-reduced)
latticeParameters = [3.8792e+000, 3.8830e+000, 1.1654e+001, 9.0018e+001, 9.0039e+001, 9.0057e+001]

# This parameter is used to judge if two unit-cell parameters are equivalent or not, by checking if 
# their corresponding metric tensors (sij) and (tij) have a relative difference that satisfies |sij-tij| <= epsilon*max{sij, tij}.
epsilon = 0.1

# 0: quick search, 1: exhaustive search
doesPrudentSymmetrySearch = 1

# Axis for rhombohedral cells: <"Rhombohedral" or "Hexagonal"
axisForRhombohedralSymmetry = 'Hexagonal'

# Unique axis for base-centered monoclinic cells: "A", "B", or "C"
axisForBaseCenteredSymmetry = 'B'
