''' << Unit-cell Parameters with the Minimal Distance >>
Bravais_type, Number_of_candidates, a b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
       Cubic(F)  0
       Cubic(I)  1  1.8886e+01  1.8886e+01  1.8886e+01  9.0000e+01  9.0000e+01  9.0000e+01  2.0738e+00
       Cubic(P)  0
      Hexagonal  0
   Rhombohedral  4  2.6684e+01  2.6684e+01  1.6384e+01  9.0000e+01  9.0000e+01  1.2000e+02  2.1613e+00
  Tetragonal(I)  3  1.8897e+01  1.8897e+01  1.8862e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.7561e+00
  Tetragonal(P)  0
Orthorhombic(F)  3  1.8862e+01  2.6698e+01  2.6751e+01  9.0000e+01  9.0000e+01  9.0000e+01  2.0413e+00
Orthorhombic(I)  1  1.8862e+01  1.8897e+01  1.8898e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.7554e+00
Orthorhombic(C)  3  1.8862e+01  1.8898e+01  1.8897e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.7554e+00
Orthorhombic(P)  0
  Monoclinic(C)  9  1.6384e+01  2.6698e+01  1.6348e+01  9.0000e+01  1.0963e+02  9.0000e+01  8.2021e-01
  Monoclinic(P)  0
      Triclinic  1  1.6332e+01  1.6348e+01  1.6357e+01  1.0941e+02  1.0952e+02  1.0934e+02  0.0000e+00
'''
AllCandidates = {
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Cubic(F)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Cubic(I)' : [
		'1.8886e+01  1.8886e+01  1.8886e+01  9.0000e+01  9.0000e+01  9.0000e+01  2.0738e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Cubic(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Hexagonal' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Rhombohedral' : [
		'2.6684e+01  2.6684e+01  1.6384e+01  9.0000e+01  9.0000e+01  1.2000e+02  2.1613e+00',
		'2.6727e+01  2.6727e+01  1.6332e+01  9.0000e+01  9.0000e+01  1.2000e+02  3.0533e+00',
		'2.6714e+01  2.6714e+01  1.6348e+01  9.0000e+01  9.0000e+01  1.2000e+02  3.1942e+00',
		'2.6707e+01  2.6707e+01  1.6357e+01  9.0000e+01  9.0000e+01  1.2000e+02  4.2517e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Tetragonal(I)' : [
		'1.8897e+01  1.8897e+01  1.8862e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.7561e+00',
		'1.8879e+01  1.8879e+01  1.8898e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.9868e+00',
		'1.8880e+01  1.8880e+01  1.8897e+01  9.0000e+01  9.0000e+01  9.0000e+01  2.0104e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Tetragonal(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(F)' : [
		'1.8862e+01  2.6698e+01  2.6751e+01  9.0000e+01  9.0000e+01  9.0000e+01  2.0413e+00',
		'1.8898e+01  2.6661e+01  2.6737e+01  9.0000e+01  9.0000e+01  9.0000e+01  2.3679e+00',
		'1.8897e+01  2.6693e+01  2.6707e+01  9.0000e+01  9.0000e+01  9.0000e+01  3.1412e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(I)' : [
		'1.8862e+01  1.8897e+01  1.8898e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.7554e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(C)' : [
		'1.8862e+01  1.8898e+01  1.8897e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.7554e+00',
		'1.8862e+01  1.8898e+01  1.8897e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.7554e+00',
		'1.8897e+01  1.8898e+01  1.8862e+01  9.0000e+01  9.0000e+01  9.0000e+01  1.7554e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Monoclinic(C)' : [
		'1.6384e+01  2.6698e+01  1.6348e+01  9.0000e+01  1.0963e+02  9.0000e+01  8.2021e-01',
		'1.8862e+01  1.8898e+01  1.8897e+01  9.0000e+01  9.0161e+01  9.0000e+01  1.0351e+00',
		'1.6357e+01  2.6751e+01  1.6332e+01  9.0000e+01  1.0952e+02  9.0000e+01  1.1898e+00',
		'1.6384e+01  2.6661e+01  1.6357e+01  9.0000e+01  1.0949e+02  9.0000e+01  1.4161e+00',
		'1.6384e+01  2.6693e+01  1.6332e+01  9.0000e+01  1.0944e+02  9.0000e+01  1.4168e+00',
		'1.8898e+01  1.8862e+01  1.8897e+01  9.0000e+01  9.0114e+01  9.0000e+01  1.4417e+00',
		'1.6348e+01  2.6737e+01  1.6332e+01  9.0000e+01  1.0934e+02  9.0000e+01  1.5906e+00',
		'1.8898e+01  1.8897e+01  1.8862e+01  9.0000e+01  9.0030e+01  9.0000e+01  1.7357e+00',
		'1.6357e+01  2.6707e+01  1.6348e+01  9.0000e+01  1.0941e+02  9.0000e+01  2.2016e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Monoclinic(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Triclinic' : [
		'1.6332e+01  1.6348e+01  1.6357e+01  1.0941e+02  1.0952e+02  1.0934e+02  0.0000e+00'
	]
}
