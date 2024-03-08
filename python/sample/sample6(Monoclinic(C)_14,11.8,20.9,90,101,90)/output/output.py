''' << Unit-cell Parameters with the Minimal Distance >>
Bravais_type, Number_of_candidates, a b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
      Triclinic  1  9.1568e+00  9.1712e+00  2.0946e+01  8.1634e+01  8.1751e+01  8.0086e+01  0.0000e+00
  Monoclinic(P)  3  9.1712e+00  9.1568e+00  2.0946e+01  9.0000e+01  9.8366e+01  9.0000e+01  4.3962e+01
  Monoclinic(C)  10  1.4031e+01  1.1792e+01  2.0946e+01  9.0000e+01  1.0088e+02  9.0000e+01  7.1525e-01
Orthorhombic(P)  3  9.1568e+00  9.1712e+00  2.0946e+01  9.0000e+01  9.0000e+01  9.0000e+01  5.9119e+01
Orthorhombic(C)  9  9.1712e+00  4.1560e+01  9.1568e+00  9.0000e+01  9.0000e+01  9.0000e+01  7.2822e+01
Orthorhombic(I)  0
Orthorhombic(F)  1  1.1792e+01  1.4031e+01  4.1593e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.2154e+02
  Tetragonal(P)  3  9.1640e+00  9.1640e+00  2.0946e+01  9.0000e+01  9.0000e+01  9.0000e+01  5.9119e+01
  Tetragonal(I)  0
   Rhombohedral  2  1.0116e+01  1.0116e+01  6.1747e+01  9.0000e+01  9.0000e+01  1.2000e+02  6.5561e+01
      Hexagonal  0
       Cubic(P)  0
       Cubic(I)  0
       Cubic(F)  0
'''
AllCandidates = {
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Triclinic' : [
		'9.1568e+00  9.1712e+00  2.0946e+01  8.1634e+01  8.1751e+01  8.0086e+01  0.0000e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Monoclinic(P)' : [
		'9.1712e+00  9.1568e+00  2.0946e+01  9.0000e+01  9.8366e+01  9.0000e+01  4.3962e+01',
		'9.1568e+00  9.1712e+00  2.0946e+01  9.0000e+01  9.8249e+01  9.0000e+01  4.4503e+01',
		'9.1568e+00  2.0946e+01  9.1712e+00  9.0000e+01  9.9914e+01  9.0000e+01  5.5470e+01'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Monoclinic(C)' : [
		'1.4031e+01  1.1792e+01  2.0946e+01  9.0000e+01  1.0088e+02  9.0000e+01  7.1525e-01',
		'4.1560e+01  9.1712e+00  9.1568e+00  9.0000e+01  9.6121e+01  9.0000e+01  4.4832e+01',
		'4.1578e+01  9.1568e+00  9.1712e+00  9.0000e+01  9.6239e+01  9.0000e+01  4.5587e+01',
		'1.1792e+01  1.4031e+01  2.1609e+01  9.0000e+01  1.0576e+02  9.0000e+01  6.0956e+01',
		'9.1712e+00  4.1560e+01  9.1568e+00  9.0000e+01  9.9914e+01  9.0000e+01  6.9892e+01',
		'9.1568e+00  4.1578e+01  9.1712e+00  9.0000e+01  9.9914e+01  9.0000e+01  7.1378e+01',
		'1.1792e+01  1.4031e+01  2.0946e+01  9.0000e+01  9.0100e+01  9.0000e+01  7.8445e+01',
		'9.1712e+00  4.1593e+01  9.1568e+00  9.0000e+01  9.9914e+01  9.0000e+01  8.5941e+01',
		'9.1568e+00  4.3582e+01  9.1712e+00  9.0000e+01  9.9914e+01  9.0000e+01  1.9701e+02',
		'9.1712e+00  4.3593e+01  9.1568e+00  9.0000e+01  9.9914e+01  9.0000e+01  1.9709e+02'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(P)' : [
		'9.1568e+00  9.1712e+00  2.0946e+01  9.0000e+01  9.0000e+01  9.0000e+01  5.9119e+01',
		'9.1568e+00  9.1712e+00  2.0946e+01  9.0000e+01  9.0000e+01  9.0000e+01  5.9119e+01',
		'9.1568e+00  9.1712e+00  2.0946e+01  9.0000e+01  9.0000e+01  9.0000e+01  5.9119e+01'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(C)' : [
		'9.1712e+00  4.1560e+01  9.1568e+00  9.0000e+01  9.0000e+01  9.0000e+01  7.2822e+01',
		'9.1712e+00  4.1560e+01  9.1568e+00  9.0000e+01  9.0000e+01  9.0000e+01  7.2822e+01',
		'9.1568e+00  4.1578e+01  9.1712e+00  9.0000e+01  9.0000e+01  9.0000e+01  7.4249e+01',
		'9.1568e+00  4.1578e+01  9.1712e+00  9.0000e+01  9.0000e+01  9.0000e+01  7.4249e+01',
		'1.1792e+01  1.4031e+01  2.0946e+01  9.0000e+01  9.0000e+01  9.0000e+01  7.8447e+01',
		'1.1792e+01  1.4031e+01  2.0946e+01  9.0000e+01  9.0000e+01  9.0000e+01  7.8447e+01',
		'9.1712e+00  4.1593e+01  9.1568e+00  9.0000e+01  9.0000e+01  9.0000e+01  8.8340e+01',
		'9.1568e+00  4.3582e+01  9.1712e+00  9.0000e+01  9.0000e+01  9.0000e+01  1.9807e+02',
		'9.1712e+00  4.3593e+01  9.1568e+00  9.0000e+01  9.0000e+01  9.0000e+01  1.9815e+02'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(I)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(F)' : [
		'1.1792e+01  1.4031e+01  4.1593e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.2154e+02'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Tetragonal(P)' : [
		'9.1640e+00  9.1640e+00  2.0946e+01  9.0000e+01  9.0000e+01  9.0000e+01  5.9119e+01',
		'9.1640e+00  9.1640e+00  2.0946e+01  9.0000e+01  9.0000e+01  9.0000e+01  5.9119e+01',
		'9.1640e+00  9.1640e+00  2.0946e+01  9.0000e+01  9.0000e+01  9.0000e+01  5.9119e+01'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Tetragonal(I)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Rhombohedral' : [
		'1.0116e+01  1.0116e+01  6.1747e+01  9.0000e+01  9.0000e+01  1.2000e+02  6.5561e+01',
		'1.0116e+01  1.0116e+01  6.3800e+01  9.0000e+01  9.0000e+01  1.2000e+02  2.3464e+02'
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
