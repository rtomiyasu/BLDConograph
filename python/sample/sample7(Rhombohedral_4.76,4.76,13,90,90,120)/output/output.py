''' << Unit-cell Parameters with the Minimal Distance >>
Bravais_type, Number_of_candidates, a b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
      Triclinic  1  4.7428e+00  4.7802e+00  5.1093e+00  1.1689e+02  9.0170e+01  1.1969e+02  0.0000e+00
  Monoclinic(P)  0
  Monoclinic(C)  9  7.1242e+00  4.7428e+00  5.1093e+00  9.0000e+01  9.6215e+01  9.0000e+01  1.1402e-01
Orthorhombic(P)  0
Orthorhombic(C)  7  4.7839e+00  6.9610e+00  5.1829e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.3652e+00
Orthorhombic(I)  2  4.7839e+00  5.1829e+00  6.9610e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.3652e+00
Orthorhombic(F)  1  6.9610e+00  6.9816e+00  7.1242e+00  9.0000e+01  9.0000e+01  9.0000e+01  9.3912e+00
  Tetragonal(P)  0
  Tetragonal(I)  2  4.9874e+00  4.9874e+00  6.9610e+00  9.0000e+01  9.0000e+01  9.0000e+01  6.0576e+00
   Rhombohedral  2  4.7690e+00  4.7690e+00  1.3074e+01  9.0000e+01  9.0000e+01  1.2000e+02  1.9315e+00
      Hexagonal  0
       Cubic(P)  0
       Cubic(I)  0
       Cubic(F)  1  7.0226e+00  7.0226e+00  7.0226e+00  9.0000e+01  9.0000e+01  9.0000e+01  9.5568e+00
'''
AllCandidates = {
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Triclinic' : [
		'4.7428e+00  4.7802e+00  5.1093e+00  1.1689e+02  9.0170e+01  1.1969e+02  0.0000e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Monoclinic(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Monoclinic(C)' : [
		'7.1242e+00  4.7428e+00  5.1093e+00  9.0000e+01  9.6215e+01  9.0000e+01  1.1402e-01',
		'6.9610e+00  4.7839e+00  5.1829e+00  9.0000e+01  9.5978e+01  9.0000e+01  7.4072e-01',
		'6.9816e+00  4.7802e+00  5.1725e+00  9.0000e+01  9.6031e+01  9.0000e+01  8.5304e-01',
		'4.7839e+00  6.9610e+00  5.1829e+00  9.0000e+01  9.1162e+01  9.0000e+01  5.3179e+00',
		'4.7839e+00  5.1829e+00  6.9610e+00  9.0000e+01  9.0252e+01  9.0000e+01  5.3612e+00',
		'5.1725e+00  6.9816e+00  4.7802e+00  9.0000e+01  9.1332e+01  9.0000e+01  5.3717e+00',
		'6.9816e+00  5.1725e+00  4.7802e+00  9.0000e+01  9.0314e+01  9.0000e+01  5.4266e+00',
		'5.1093e+00  7.1242e+00  4.7428e+00  9.0000e+01  9.0170e+01  9.0000e+01  5.5731e+00',
		'7.1242e+00  5.1093e+00  4.7428e+00  9.0000e+01  9.0062e+01  9.0000e+01  5.5738e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(C)' : [
		'4.7839e+00  6.9610e+00  5.1829e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.3652e+00',
		'4.7839e+00  6.9610e+00  5.1829e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.3652e+00',
		'4.7839e+00  5.1829e+00  6.9610e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.3652e+00',
		'5.1725e+00  6.9816e+00  4.7802e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.4328e+00',
		'5.1725e+00  6.9816e+00  4.7802e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.4328e+00',
		'5.1093e+00  7.1242e+00  4.7428e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.5740e+00',
		'5.1093e+00  7.1242e+00  4.7428e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.5740e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(I)' : [
		'4.7839e+00  5.1829e+00  6.9610e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.3652e+00',
		'4.7802e+00  5.1725e+00  6.9816e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.4328e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(F)' : [
		'6.9610e+00  6.9816e+00  7.1242e+00  9.0000e+01  9.0000e+01  9.0000e+01  9.3912e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Tetragonal(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Tetragonal(I)' : [
		'4.9874e+00  4.9874e+00  6.9610e+00  9.0000e+01  9.0000e+01  9.0000e+01  6.0576e+00',
		'4.9802e+00  4.9802e+00  6.9816e+00  9.0000e+01  9.0000e+01  9.0000e+01  6.0940e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Rhombohedral' : [
		'4.7690e+00  4.7690e+00  1.3074e+01  9.0000e+01  9.0000e+01  1.2000e+02  1.9315e+00',
		'5.0369e+00  5.0369e+00  1.1807e+01  9.0000e+01  9.0000e+01  1.2000e+02  5.3528e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Hexagonal' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Cubic(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Cubic(I)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Cubic(F)' : [
		'7.0226e+00  7.0226e+00  7.0226e+00  9.0000e+01  9.0000e+01  9.0000e+01  9.5568e+00'
	]
}
