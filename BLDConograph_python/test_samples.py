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

def putSobs (latticeParameters):
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
    out =  (
            input.sampleName,
            input.latticeParameters,
            input.doesPrudentSymmetrySearch,
            input.axisForRhombohedralSymmetry,
            input.axisForBaseCenteredSymmetry,
            input.epsilon
                                            )
    sys.path.remove (dir) # Remove the path to dir.
    return out

def ProcessBLD (Sobs, doesPrudent ,axisRhombohedral, axisBaseCenter ,eps):
    """execute prcess to produce projection candidates of
    each 14 bravais lattice pattern from lattice parameter.
    input :
        Sobs : symmetry matrix converted from lattice param
                [a,b,c,alpha,beta,gamma]
        doesPrudent : flg of prudent search
        axisHRHombohedral : axis for lattice pattern of Rhombohedral
                            ('Rhombohedral or 'Hexagonal)
        eps : threshold to select candidate
    output:
        dict data {lattice pattern name : [(g, S), (g, S),....], ....}
        g : conversion matrix, S : projection"""
    bc = BravaisLatticeDetermination ()
    bc_input = BravaisLatticeDetermination.InputType()
    bc_input.Sobs = Sobs
    bc_input.DoesPrudentSearch = doesPrudent
    bc_input.eps  = eps
    bc_input.AxisForRhombohedralSymmetry = axisRhombohedral
    bc_input.AxisForBaseCenteredSymmetry = axisBaseCenter
    bc.set_bravais_class (bc_input)
    return bc.result.BravaisClasses

#-------------------------------------------------------
#      Post process & make output file
#-------------------------------------------------------
def selectMatrixMinDist (name, bDict):
    """select projection by minimum dist value between projection & gSobstg"""
    dists = np.array ([d for _, _, d in bDict[name]])
    mats = [S for _, S, _ in bDict[name]]
    minIdx = np.argmin (dists)
    mat = mats [minIdx]
    return mat, dists[minIdx], len (mats)
    
def putBravaisCandidates (bDict):
    """
    set projection matices and values of distance to dict data format
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

def putBravaisCandidateMinDistance (bDict):
    """put lattice parameters of min dist form candidates"""
    bcDict = {}
    for name in bDict.keys():
        bcDict[name] = {}
        # if no candidate, set NumberOfSolution 0 and continue
        if len (bDict[name]) == 0:
            bcDict[name]['NumberOfSolutions'] = 0
            continue

        # mat : projection, d : dist, length : num of candidates
        mat, d, length = selectMatrixMinDist (name, bDict)
        
        bcDict[name]['NumberOfSolutions'] = length
        bcDict[name]['DistanceFromUnitCell'] = d
        bcDict[name]['UnitCellParameter'] = sym2lattice (mat)
    return bcDict

def make_ouput_dict (sampleName, bDict):
    """make output of dict format"""
    candidateMinDist = putBravaisCandidateMinDistance (bDict)
    candidates = putBravaisCandidates (bDict)

    outputDict ={
        'SampleName' : sampleName,        
        'BravaisLatticeMinDistance' : candidateMinDist,
        'BravaisCandidates' : candidates
                                                        }
    
    return outputDict    

def write_lattice_min_dist (f, outputDict):
    maxLen = max ([len (name) for name in outputDict['BravaisLatticeMinDistance']])
    for name, valuesDict in outputDict['BravaisLatticeMinDistance'].items():
        n = valuesDict['NumberOfSolutions']
        text = [' ' * (maxLen - len (name)) + name, str (n)]
        if n > 0:
            d = valuesDict ['DistanceFromUnitCell']
            lattice = valuesDict ['UnitCellParameter']
            text += ['{:.4e}'.format (lt) for lt in lattice]
            text += ['{:.4e}'.format (d)]    
        text = '  '.join (text) + '\n'
        f.write (text)

def write_lattice_candidates (f, outputDict):
    text = 'BravaisCandidates = {\n'
    for name, valuesDict in outputDict['BravaisCandidates'].items():
        text += '\t' + '# a, b, c, alpha, beta, gamma, dist\n'
        text += '\t' + "'" + name + "'" + ' : [\n'

        if 'UnitCellParameters' in valuesDict:
            for i, (lattice, d) in enumerate (zip (
                                valuesDict['UnitCellParameters'],
                                valuesDict['DistanceFromInputUnitCell'])):
                text += '\t' * 2
                text += "'" + '  '.join (['{:.4e}'.format(lt) for lt in (lattice + [d])]) + "'"
                
                if i < len (valuesDict['UnitCellParameters']) - 1:
                    text += ',\n'
                else:
                    text += '\n'
        text += '\t' + '],\n'
    text += '}'
    f. write (text)

def renewalBravaisTypeName (outputDict, axis):
    """rename lattice pattern name of output dict
       ex base-centered monoclinic --> Monoclinic (C)
       axis : axisBaseCenter = A, B or C"""
    outputDict['BravaisLatticeMinDistance'] = {
        BravaisLatticeDetermination.key2str (k, axis) : v for k, v in outputDict['BravaisLatticeMinDistance'].items()}
    outputDict['BravaisCandidates'] = {
        BravaisLatticeDetermination.key2str (k, axis) : v for k , v in outputDict['BravaisCandidates'].items()}
    return outputDict

def output_to_py_file (outputDict, dir, axisBaseCenter):
    # make summary document file of BLD Conograph
    # & output to py file
    outputDict = renewalBravaisTypeName (outputDict, axisBaseCenter)
    savePath = os.path.join (dir, 'output', outputDict['SampleName'] + '.out.py')
    comment = '# Result is divided 2 parts of BravaisLatticeMinDistance & BravaisCandidates.\n'
    comment += '# BravaisLatticeMinDistance : Selected Bravais lattice parameters of Minimum distance\n'
    comment +='# BravaisCandidates : Lists of all candidates of lattice parameters\n'
    with open (savePath, 'w', encoding='utf-8') as f:
        f.write ('sampleName = ' + "'" + outputDict['SampleName'] + "'" + '\n')
        f.write (comment)
        f.write ("'''\n<<Lattice Parameters Minimum Distance>>\n")
        f.write ('PatternName numCandidate, a b, c, alpha, beta, gamma, dist\n')

        # lattice parameters min distance
        write_lattice_min_dist (f, outputDict)        
        f.write ("'''\n")
        # lattice candidates
        write_lattice_candidates (f, outputDict)

#------------------------------------------------------
#    main
#------------------------------------------------------
def process_BLDConograph (dirList = None):
    # If dirList = None, all the 'input.py' files in 'sample' are searched for and processed.
    if dirList is None:
        dirList = make_dir_list ('sample')
    
    # Loop for all the input.py.
    for dir in dirList:

        if not os.path.exists (os.path.join (dir, 'input.py')):
            print ('There is not input.py in {}...'.format(dir))
            continue
        
        sampleName, latticeParams, \
        doesPrudent, axisRhom, axisBaseCenter, eps = read_input (dir)
        if not isinstance (doesPrudent, bool):
            doesPrudent = bool (doesPrudent)

        # print logs
        print ('-----folder name {} sample name : {}-----'.format (dir, sampleName))
        print ('Parameters doesPrudent : {}, axis rhombohedral sym : {}, axis base center : {}'.format(
            int (doesPrudent), axisRhom, axisBaseCenter))
        print ('<<lattice params a, b, c, alpha, beta, gamma>>')
        print ('  '.join (['{:.4e}'.format(v) for v in list (latticeParams)]))

        # Make a folder 'output' in dir. Nothing is done if it exists.
        os.makedirs (os.path.join (dir, 'output'), 
                    exist_ok = True)

        # Unit-cell parameters --> symmetric matrix Sobs
        Sobs = putSobs (latticeParams)
        try:
            # Bravais Lattice Determination
            bDict = ProcessBLD (Sobs, doesPrudent,
                        axisRhom, axisBaseCenter, eps)
        except AssertionError as e:
            print (e)
            continue

        # build the dict data for output file
        outputDict = make_ouput_dict (sampleName, bDict)

        # save the result to py file with renewal bravais type name
        output_to_py_file (outputDict, dir, axisBaseCenter)
        print ('completed the BLD process : {} ....\n'.format(dir))

if __name__ == '__main__':
    process_BLDConograph ()
