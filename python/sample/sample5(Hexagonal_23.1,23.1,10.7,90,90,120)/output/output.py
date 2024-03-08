''' << Unit-cell Parameters with the Minimal Distance >>
Bravais_type, Number_of_candidates, a b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
      Triclinic  1  1.0654e+01  2.2895e+01  2.3003e+01  6.0600e+01  8.9377e+01  8.9119e+01  0.0000e+00
  Monoclinic(P)  1  2.2895e+01  1.0654e+01  2.3003e+01  9.0000e+01  1.1940e+02  9.0000e+01  6.5065e+00
  Monoclinic(C)  8  3.9628e+01  2.3157e+01  1.0654e+01  9.0000e+01  9.0871e+01  9.0000e+01  7.1764e+00
Orthorhombic(P)  0
Orthorhombic(C)  8  2.3157e+01  3.9628e+01  1.0654e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.1568e+01
Orthorhombic(I)  0
Orthorhombic(F)  0
  Tetragonal(P)  0
  Tetragonal(I)  0
   Rhombohedral  0
      Hexagonal  1  2.3157e+01  2.3157e+01  1.0654e+01  9.0000e+01  9.0000e+01  1.2000e+02  2.0549e+01
       Cubic(P)  0
       Cubic(I)  0
       Cubic(F)  0
'''
AllCandidates = {
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Triclinic' : [
		'1.0654e+01  2.2895e+01  2.3003e+01  6.0600e+01  8.9377e+01  8.9119e+01  0.0000e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Monoclinic(P)' : [
		'2.2895e+01  1.0654e+01  2.3003e+01  9.0000e+01  1.1940e+02  9.0000e+01  6.5065e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Monoclinic(C)' : [
		'3.9628e+01  2.3157e+01  1.0654e+01  9.0000e+01  9.0871e+01  9.0000e+01  7.1764e+00',
		'2.2895e+01  4.0082e+01  1.0654e+01  9.0000e+01  9.0881e+01  9.0000e+01  1.0298e+01',
		'4.0082e+01  2.2895e+01  1.0654e+01  9.0000e+01  9.0212e+01  9.0000e+01  1.1366e+01',
		'2.3157e+01  3.9628e+01  1.0654e+01  9.0000e+01  9.0252e+01  9.0000e+01  1.1465e+01',
		'3.9896e+01  2.3003e+01  1.0654e+01  9.0000e+01  9.0652e+01  9.0000e+01  1.7474e+01',
		'2.3003e+01  3.9896e+01  1.0654e+01  9.0000e+01  9.0623e+01  9.0000e+01  1.8383e+01',
		'4.0879e+01  1.0654e+01  2.3157e+01  9.0000e+01  9.0366e+01  9.0000e+01  1.5146e+02',
		'4.1177e+01  1.0654e+01  2.3003e+01  9.0000e+01  9.0891e+01  9.0000e+01  1.5373e+02'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(C)' : [
		'2.3157e+01  3.9628e+01  1.0654e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.1568e+01',
		'2.3157e+01  3.9628e+01  1.0654e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.1568e+01',
		'2.2895e+01  4.0082e+01  1.0654e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.1584e+01',
		'2.2895e+01  4.0082e+01  1.0654e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.1584e+01',
		'2.3003e+01  3.9896e+01  1.0654e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.8765e+01',
		'2.3003e+01  3.9896e+01  1.0654e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.8765e+01',
		'1.0654e+01  4.0879e+01  2.3157e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.5170e+02',
		'1.0654e+01  4.1177e+01  2.3003e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.5514e+02'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(I)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(F)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Tetragonal(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Tetragonal(I)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Rhombohedral' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Hexagonal' : [
		'2.3157e+01  2.3157e+01  1.0654e+01  9.0000e+01  9.0000e+01  1.2000e+02  2.0549e+01'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Cubic(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Cubic(I)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Cubic(F)' : [
	]
}
