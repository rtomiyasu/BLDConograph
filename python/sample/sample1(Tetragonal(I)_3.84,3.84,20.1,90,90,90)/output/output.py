''' << Unit-cell Parameters with the Minimal Distance >>
Bravais_type, Number_of_candidates, a b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
      Triclinic  1  3.8417e+00  3.8524e+00  1.0389e+01  1.0041e+02  1.0053e+02  9.0046e+01  0.0000e+00
  Monoclinic(P)  0
  Monoclinic(C)  9  5.4384e+00  5.4427e+00  1.0389e+01  9.0000e+01  1.0490e+02  9.0000e+01  1.4585e-01
Orthorhombic(P)  0
Orthorhombic(C)  7  3.8417e+00  2.0080e+01  3.8524e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.6424e-01
Orthorhombic(I)  1  3.8417e+00  3.8524e+00  2.0080e+01  9.0000e+01  9.0000e+01  9.0000e+01  5.6424e-01
Orthorhombic(F)  1  5.4384e+00  5.4427e+00  2.0080e+01  9.0000e+01  9.0000e+01  9.0000e+01  8.0605e-01
  Tetragonal(P)  0
  Tetragonal(I)  1  3.8471e+00  3.8471e+00  2.0080e+01  9.0000e+01  9.0000e+01  9.0000e+01  5.6723e-01
   Rhombohedral  2  4.4413e+00  4.4413e+00  3.0237e+01  9.0000e+01  9.0000e+01  1.2000e+02  2.1294e+01
      Hexagonal  0
       Cubic(P)  0
       Cubic(I)  0
       Cubic(F)  0
'''
AllCandidates = {
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Triclinic' : [
		'3.8417e+00  3.8524e+00  1.0389e+01  1.0041e+02  1.0053e+02  9.0046e+01  0.0000e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Monoclinic(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Monoclinic(C)' : [
		'5.4384e+00  5.4427e+00  1.0389e+01  9.0000e+01  1.0490e+02  9.0000e+01  1.4585e-01',
		'2.0080e+01  3.8417e+00  3.8524e+00  9.0000e+01  9.0271e+01  9.0000e+01  2.2566e-01',
		'5.4427e+00  5.4384e+00  1.0407e+01  9.0000e+01  1.0527e+02  9.0000e+01  4.4481e-01',
		'2.0080e+01  3.8524e+00  3.8417e+00  9.0000e+01  9.0118e+01  9.0000e+01  5.1742e-01',
		'3.8524e+00  2.0080e+01  3.8417e+00  9.0000e+01  9.0046e+01  9.0000e+01  5.6399e-01',
		'3.8417e+00  2.0428e+01  3.8524e+00  9.0000e+01  9.0046e+01  9.0000e+01  2.0473e+01',
		'3.8524e+00  2.0436e+01  3.8417e+00  9.0000e+01  9.0046e+01  9.0000e+01  2.0654e+01',
		'3.8524e+00  2.0452e+01  3.8417e+00  9.0000e+01  9.0046e+01  9.0000e+01  2.1103e+01',
		'3.8417e+00  2.0464e+01  3.8524e+00  9.0000e+01  9.0046e+01  9.0000e+01  2.1506e+01'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(C)' : [
		'3.8417e+00  2.0080e+01  3.8524e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.6424e-01',
		'3.8524e+00  2.0080e+01  3.8417e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.6424e-01',
		'3.8524e+00  2.0080e+01  3.8417e+00  9.0000e+01  9.0000e+01  9.0000e+01  5.6424e-01',
		'3.8417e+00  2.0428e+01  3.8524e+00  9.0000e+01  9.0000e+01  9.0000e+01  2.0473e+01',
		'3.8524e+00  2.0436e+01  3.8417e+00  9.0000e+01  9.0000e+01  9.0000e+01  2.0654e+01',
		'3.8524e+00  2.0452e+01  3.8417e+00  9.0000e+01  9.0000e+01  9.0000e+01  2.1103e+01',
		'3.8417e+00  2.0464e+01  3.8524e+00  9.0000e+01  9.0000e+01  9.0000e+01  2.1506e+01'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(I)' : [
		'3.8417e+00  3.8524e+00  2.0080e+01  9.0000e+01  9.0000e+01  9.0000e+01  5.6424e-01'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(F)' : [
		'5.4384e+00  5.4427e+00  2.0080e+01  9.0000e+01  9.0000e+01  9.0000e+01  8.0605e-01'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Tetragonal(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Tetragonal(I)' : [
		'3.8471e+00  3.8471e+00  2.0080e+01  9.0000e+01  9.0000e+01  9.0000e+01  5.6723e-01'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Rhombohedral' : [
		'4.4413e+00  4.4413e+00  3.0237e+01  9.0000e+01  9.0000e+01  1.2000e+02  2.1294e+01',
		'4.4413e+00  4.4413e+00  3.0248e+01  9.0000e+01  9.0000e+01  1.2000e+02  2.1724e+01'
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
	]
}
