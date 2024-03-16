import numpy as np
import os
import sys
sys.path.append ('for_3D')
import importlib
from for_3D.error_stable_bravais_3D import BravaisLatticeDetermination

def make_dir_list (root = 'sample', fname = 'input.py'):
    """ Output all the folders in root containing input.py. """
    dirList = []
    for dir, _, files in os.walk (root): # All the subdirectories of root are checked.
        # dir : folders in root = 'sample'
        # files : files in dir
        if fname in files:
            # The directroy name dir is added if dir contains a file with fname. 
            dirList.append (dir)
    return dirList

def lattice2sym (latticeParameters):
    """ Convert unit-cell parameters to the symmetric matrix:
    [a^2, a*b*cos(gamma), a*c*cos(beta)],
    [a*b*cos(gamma), b^2, b*c*cos(alpha)],
    [a*c*cos(beta), b*c*cos(alpha), c^2]"""
    # Unit cell parameters: 
    a, b, c, alpha, beta, gamma = latticeParameters
    # degree --> radian
    alpha = np.deg2rad (alpha)
    beta = np.deg2rad (beta)
    gamma = np.deg2rad (gamma)
    arr = np.array ([
        [a ** 2, a * b * np.cos (gamma), a * c * np.cos (beta)],
        [a * b * np.cos (gamma), b ** 2, b * c * np.cos (alpha)],
        [a * c * np.cos (beta), b * c * np.cos (alpha), c ** 2]
    ])
    return arr

def sym2lattice (S):
    """ Convert a symmetry matrix to unit-cell parameters a,b,c,α,β,γ"""
    a, b, c = np.sqrt (np.diag (S))
    alpha = np.rad2deg (np.arccos (S[1,2] / (b * c)))
    beta = np.rad2deg (np.arccos (S[0,2] / (a * c)))
    gamma = np.rad2deg (np.arccos (S[0,1] / (a * b)))
    return [a, b, c, alpha, beta, gamma]

def read_input (dir):
    """ Read dir/input.py. """
    sys.path.append (dir) # Add the path to dir.
    # Read input.py.
    import input
    importlib.reload (input)
    
    # Read the input parameters.
    bc_input = BravaisLatticeDetermination.InputType()
    bc_input.Sobs = lattice2sym(input.latticeParameters) # Unit-cell parameters --> symmetric matrix Sobs
    bc_input.doesPrudentSearch   = input.doesPrudentSymmetrySearch
    bc_input.axisForRhombohedralSymmetry = input.axisForRhombohedralSymmetry
    bc_input.axisForBaseCenteredSymmetry = input.axisForBaseCenteredSymmetry
    bc_input.epsilon = input.epsilon
    
    sys.path.remove (dir) # Remove the path to dir.
    return (input.latticeParameters, bc_input)


#-------------------------------------------------------
#      Post process & make output file
#-------------------------------------------------------
def output_file (bDict, dir):
    # Output the results in bDict.
    savePath = os.path.join (dir, 'output', 'output.py')
    with open (savePath, 'w', encoding='utf-8') as f:
        f.write ("''' << Unit-cell Parameters with the Minimal Distance >>\n")
        f.write ('Bravais_type, Number_of_candidates, a b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell\n')
        
        ''' Unit-cell parameters with the minimal distances from the input unit cell.'''
        maxLen = max ([len (name) for name in bDict])
        for name in bDict.keys():
            n = len(bDict[name])
            text = [' ' * (maxLen - len (name)) + name, str (n)]
            if n > 0:
                lattice = sym2lattice( bDict[name][0]['GramMatrix'] )
                d = bDict[name][0]['DistanceFromInput']
                text += ['{:.4e}'.format (lt) for lt in lattice] # S
                text += ['{:.4e}'.format (d)] # distance    
            text = '  '.join (text) + '\n'
            f. write (text)
        f.write ("'''\n")
        
        ''' All the unit-cell parameters'''
        text = 'AllCandidates = {\n'
        for i, name in enumerate (bDict.keys()):
            text += '\t' + '# a, b, c, alpha, beta, gamma, Distance_from_the_input_unit_cell\n'
            text += '\t' + "'" + name + "'" + ' : [\n'
            for j in range(len(bDict[name])):
                lattice = sym2lattice( bDict[name][j]['GramMatrix'] )
                d = bDict[name][j]['DistanceFromInput']
                text += '\t' * 2
                text += "'" + '  '.join (['{:.4e}'.format(lt) for lt in (lattice + [d])]) + "'"
                if j + 1 < len (bDict[name]):
                    text += ','
                text += '\n'
            text += '\t' + ']'
            if i + 1 < len (bDict):
                text += ','
            text += '\n'
        text += '}\n'
        f. write (text)

#------------------------------------------------------
#    main
#------------------------------------------------------
def test_samples (dirList = None):
    # If dirList = None, all the 'input.py' files in 'sample' are searched for and processed.
    if dirList is None:
        dirList = make_dir_list ('sample')
    
    # Loop for all the input.py.
    for dir in dirList:

        if not os.path.exists (os.path.join (dir, 'input.py')):
            print ('There is not input.py in {}...'.format(dir))
            continue
        
        # Read input parameters.
        latticeParams, bc_input = read_input (dir) 
        if not isinstance (bc_input.doesPrudentSearch, bool):
            bc_input.doesPrudentSearch = bool (bc_input.doesPrudentSearch)

        print ('-----folder name: {} -----'.format (dir))
        print ('<< Input parameters of primitive unit-cells: a, b, c, alpha, beta, gamma >>')
        print ('  '.join (['{:.4e}'.format(v) for v in list (latticeParams)]))

        # Make a folder 'output' in dir. Nothing is done if it exists.
        os.makedirs (os.path.join (dir, 'output'), exist_ok = True)

        try:
            # Bravais Lattice Determination
            bc = BravaisLatticeDetermination ()
            bc.set_bravais_class (bc_input)
        except AssertionError as e:
            print (e)
            continue

        # Save the result in output.py.
        output_file (bc.result.toDictionary(bc_input), dir)
        print ('Output file: {}\n'.format('output.py'))

if __name__ == '__main__':
    test_samples ()
