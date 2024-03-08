''' << Unit-cell Parameters with the Minimal Distance >>
Bravais_type, Number_of_candidates, a b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
      Triclinic  1  3.7250e+00  3.7269e+00  3.7286e+00  1.1998e+02  9.0070e+01  1.1991e+02  0.0000e+00
  Monoclinic(P)  0
  Monoclinic(C)  9  5.2672e+00  3.7310e+00  3.7291e+00  9.0000e+01  9.0057e+01  9.0000e+01  1.4161e-02
Orthorhombic(P)  0
Orthorhombic(C)  9  3.7250e+00  5.2764e+00  3.7286e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.0928e-02
Orthorhombic(I)  3  3.7250e+00  3.7286e+00  5.2764e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.0928e-02
Orthorhombic(F)  1  5.2672e+00  5.2737e+00  5.2764e+00  9.0000e+01  9.0000e+01  9.0000e+01  4.6516e-02
  Tetragonal(P)  0
  Tetragonal(I)  3  3.7300e+00  3.7300e+00  5.2672e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.2839e-02
   Rhombohedral  4  3.7276e+00  3.7276e+00  9.1349e+00  9.0000e+01  9.0000e+01  1.2000e+02  6.3860e-02
      Hexagonal  0
       Cubic(P)  0
       Cubic(I)  0
       Cubic(F)  1  5.2725e+00  5.2725e+00  5.2725e+00  9.0000e+01  9.0000e+01  9.0000e+01  8.4427e-02
'''
AllCandidates = {
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Triclinic' : [
		'3.7250e+00  3.7269e+00  3.7286e+00  1.1998e+02  9.0070e+01  1.1991e+02  0.0000e+00'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Monoclinic(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Monoclinic(C)' : [
		'5.2672e+00  3.7310e+00  3.7291e+00  9.0000e+01  9.0057e+01  9.0000e+01  1.4161e-02',
		'3.7286e+00  5.2764e+00  3.7250e+00  9.0000e+01  9.0070e+01  9.0000e+01  1.9370e-02',
		'5.2764e+00  3.7286e+00  3.7250e+00  9.0000e+01  9.0040e+01  9.0000e+01  2.4151e-02',
		'3.7310e+00  5.2672e+00  3.7291e+00  9.0000e+01  9.0030e+01  9.0000e+01  2.9450e-02',
		'3.7310e+00  3.7291e+00  5.2672e+00  9.0000e+01  9.0020e+01  9.0000e+01  2.9570e-02',
		'3.7286e+00  5.2737e+00  3.7269e+00  9.0000e+01  9.0100e+01  9.0000e+01  3.0352e-02',
		'5.2764e+00  3.7250e+00  3.7286e+00  9.0000e+01  9.0003e+01  9.0000e+01  3.0897e-02',
		'5.2737e+00  3.7286e+00  3.7269e+00  9.0000e+01  9.0060e+01  9.0000e+01  3.5322e-02',
		'5.2737e+00  3.7269e+00  3.7286e+00  9.0000e+01  9.0017e+01  9.0000e+01  4.5016e-02'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(C)' : [
		'3.7250e+00  5.2764e+00  3.7286e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.0928e-02',
		'3.7286e+00  5.2764e+00  3.7250e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.0928e-02',
		'3.7286e+00  5.2764e+00  3.7250e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.0928e-02',
		'3.7310e+00  5.2672e+00  3.7291e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.1163e-02',
		'3.7310e+00  5.2672e+00  3.7291e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.1163e-02',
		'3.7291e+00  3.7310e+00  5.2672e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.1163e-02',
		'3.7286e+00  5.2737e+00  3.7269e+00  9.0000e+01  9.0000e+01  9.0000e+01  4.5800e-02',
		'3.7269e+00  5.2737e+00  3.7286e+00  9.0000e+01  9.0000e+01  9.0000e+01  4.5800e-02',
		'3.7286e+00  5.2737e+00  3.7269e+00  9.0000e+01  9.0000e+01  9.0000e+01  4.5800e-02'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(I)' : [
		'3.7250e+00  3.7286e+00  5.2764e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.0928e-02',
		'3.7291e+00  3.7310e+00  5.2672e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.1163e-02',
		'3.7269e+00  3.7286e+00  5.2737e+00  9.0000e+01  9.0000e+01  9.0000e+01  4.5800e-02'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Orthorhombic(F)' : [
		'5.2672e+00  5.2737e+00  5.2764e+00  9.0000e+01  9.0000e+01  9.0000e+01  4.6516e-02'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Tetragonal(P)' : [
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Tetragonal(I)' : [
		'3.7300e+00  3.7300e+00  5.2672e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.2839e-02',
		'3.7268e+00  3.7268e+00  5.2764e+00  9.0000e+01  9.0000e+01  9.0000e+01  3.6192e-02',
		'3.7278e+00  3.7278e+00  5.2737e+00  9.0000e+01  9.0000e+01  9.0000e+01  4.6669e-02'
	],
	# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell
	'Rhombohedral' : [
		'3.7276e+00  3.7276e+00  9.1349e+00  9.0000e+01  9.0000e+01  1.2000e+02  6.3860e-02',
		'3.7294e+00  3.7294e+00  9.1263e+00  9.0000e+01  9.0000e+01  1.2000e+02  9.2164e-02',
		'3.7276e+00  3.7276e+00  9.1353e+00  9.0000e+01  9.0000e+01  1.2000e+02  9.4100e-02',
		'3.7282e+00  3.7282e+00  9.1322e+00  9.0000e+01  9.0000e+01  1.2000e+02  1.1047e-01'
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
		'5.2725e+00  5.2725e+00  5.2725e+00  9.0000e+01  9.0000e+01  9.0000e+01  8.4427e-02'
	]
}
