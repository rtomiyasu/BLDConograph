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
def putBravaisCandidateMinDistance (bDict):
    """put lattice parameters with the minimal distance from the input parameters """
    bcDict = {}
    for name in bDict.keys():
        bcDict[name] = {}
        # if no candidate, set NumberOfSolution 0 and continue
        if len (bDict[name]) == 0:
            bcDict[name]['NumberOfCandidates'] = 0
            continue
        bcDict[name]['NumberOfCandidates'] = len(bDict[name])
        bcDict[name]['DistanceFromInputUnitCell'] = bDict[name][0][2]
        bcDict[name]['UnitCellParameter'] = sym2lattice (bDict[name][0][1])
    return bcDict

def putBravaisCandidates (bDict):
    """ put lattice parameters  projection matices and values of distance to dict data format
    output : 
    {lattice pattern name : 
                    'UnitCellParameters' : [
                            [a,b,c,α,β,γ], [a,b,c,α,β,γ],...],
                    'DIstanceFromUnitInputUnitCell' : [
                            d0, d1, d2...]
    ・・・・}"""
    bcDict = {}
    for name in bDict.keys():
        bcDict[name] = {}
        if len (bDict[name]) == 0:
            continue
        bcDict[name]['UnitCellParameters'] = [
            sym2lattice (S) for _, S, _ in bDict[name]]
        bcDict[name]['DistanceFromInputUnitCell'] = [
            d for _, _, d in bDict[name]]
    return bcDict        

def write_lattice_min_dist (f, bDict):
    '''write bravas candidates pair (lattice parameter & distance)
    of minnimum distance of each bravais patterns'''
    # select pairs (lattice param & distance) of min distance
    bcDict = putBravaisCandidateMinDistance (bDict)

    maxLen = max ([len (name) for name in bcDict])
    for name, valuesDict in bcDict.items():
        n = valuesDict['NumberOfCandidates']
        text = [' ' * (maxLen - len (name)) + name, str (n)]
        if n > 0:
            d = valuesDict ['DistanceFromInputUnitCell']
            lattice = valuesDict ['UnitCellParameter']
            text += ['{:.4e}'.format (lt) for lt in lattice]
            text += ['{:.4e}'.format (d)]    
        text = '  '.join (text) + '\n'
        f.write (text)

def write_lattice_candidates (f, bDict):
    '''write all bravais candidates'''
    #prepare the projection matrixces & distance pair of all bravais candidates
    #   for output to file.
    bcDict = putBravaisCandidates (bDict)
    text = 'AllCandidates = {\n'
    for i, (name, valuesDict) in enumerate (bcDict.items()):
        text += '\t' + '# a, b, c, alpha, beta, gamma, Distance from the input unit cell\n'
        text += '\t' + "'" + name + "'" + ' : [\n'

        if 'UnitCellParameters' in valuesDict:
            for j, (lattice, d) in enumerate (zip (
                                valuesDict['UnitCellParameters'],
                                valuesDict['DistanceFromInputUnitCell'])):
                text += '\t' * 2
                text += "'" + '  '.join (['{:.4e}'.format(lt) for lt in (lattice + [d])]) + "'"
                
                if j < len (valuesDict['UnitCellParameters']) - 1:
                    text += ',\n'
                else:
                    text += '\n'
        if i < len (bcDict) - 1:
            text += '\t' + '],\n'
        else:
            text += '\t' + ']\n'
    text += '}'
    f. write (text)

def renewalBravaisTypeName (bDict, axis):
    """rename lattice pattern name of output dict
       ex base-centered monoclinic --> Monoclinic (C)
       axis : axisBaseCenter = A, B or C"""
    bDict = {
        BravaisLatticeDetermination.key2str (k, axis) : v for k, v in bDict.items()}
    return bDict

def output_to_py_file (bDict, dir, axisBaseCenter):
    # make summary document file of BLD Conograph
    # & output to py file
    bDict = renewalBravaisTypeName (bDict, axisBaseCenter)
    savePath = os.path.join (dir, 'output', 'output.py')
    with open (savePath, 'w', encoding='utf-8') as f:
        f.write ("''' << Unit-cell Parameters with the Minimal Distance >>\n")
        f.write ('Bravais type, Number of candidates, a b, c, alpha, beta, gamma, Distance from the input unit cell\n')

        # lattice parameters min distance
        write_lattice_min_dist (f, bDict)        
        f.write ("'''\n")
        # lattice candidates
        write_lattice_candidates (f, bDict)

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

        # save the result to py file with renewal bravais type name
        output_to_py_file (bc.result.BravaisClasses, dir, bc_input.axisForBaseCenteredSymmetry)
        print ('completed all: {} ....\n'.format(dir))

if __name__ == '__main__':
    test_samples ()
