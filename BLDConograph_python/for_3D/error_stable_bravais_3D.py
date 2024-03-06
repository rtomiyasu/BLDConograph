import numpy as np
import sys
sys.path.append('./for_2D')
import mlist
from for_2D.error_stable_bravais_2D import gauss_algorithm
from Buerger_reduction import Buerger_reduction, dist, check_equiv
from Selling_Delaunay_reduction import Delaunay_reduction


def Delaunay_reduction_of_inverse (S):
    """ Selling reduction of S^{-1}.
    input : 3x3 symmetric positive-definite S. 
    output: 3x3 basis transform matrix g 
             such that (g S g^T)^{-1} is Selling reduced (i.e., all non-diagonal entries are <= 0). """
    Sinv = np.linalg.inv (S)
    h, Sinv_red = Delaunay_reduction (Sinv)
    assert abs( np.linalg.det (h) ) == 1
    # g = (h.T)^-1
    g = np.linalg.inv (h.T).astype (int)
    return g, g.dot (S).dot (g.T)

def permute_matrix (ps = [0,1,2]):
    # Returns a 3x3 permutation matrix (all 6 matrices are available).
    # E.g., if ps = [2,0,1], the returned matrix is [[0,0,1],[1,0,0],[0,1,0]].
    assert (min (ps) >= 0) & (max (ps) < 3)
    mat = np.zeros ((3,3), dtype="int")
    for row, p in enumerate (ps):
        mat [row, p] = 1
    return mat


#-------------------------------------------------
#  Class : BravaisLatticeDetermination
#-------------------------------------------------
class BravaisLatticeDetermination:
    __listNames = ['triclinic', 'primitive monoclinic', 'base-centered monoclinic',
                   'primitive orthorhombic', 'base-centered orthorhombic', 
                   'body-centered orthorhombic', 'face-centered orthorhombic',
                   'primitive tetragonal', 'body-centered tetragonal',
                   'rhombohedral', 'hexagonal',
                   'primitive cubic', 'body-centered cubic', 'face-centered cubic']
    
    @staticmethod
    def key2str(key, AxisForBaseCenteredSymmetry):
        __listStrings = { 'triclinic'                : 'Triclinic', 
                          'primitive monoclinic'     : 'Monoclinic(P)',
                          'base-centered monoclinic' : 'Monoclinic',
                          'primitive orthorhombic'   : 'Orthorhombic(P)',
                          'base-centered orthorhombic'   : 'Orthorhombic(C)', 
                          'body-centered orthorhombic'   : 'Orthorhombic(I)',
                          'face-centered orthorhombic'   : 'Orthorhombic(F)',
                          'primitive tetragonal'     : 'Tetragonal(P)',
                          'body-centered tetragonal' : 'Tetragonal(I)',
                          'rhombohedral'             : 'Rhombohedral',
                          'hexagonal'                : 'Hexagonal',
                          'primitive cubic'          : 'Cubic(P)', 
                          'body-centered cubic'      : 'Cubic(I)', 
                          'face-centered cubic'      : 'Cubic(F)' }
        if key == 'base-centered monoclinic':
            if AxisForBaseCenteredSymmetry == 'A':
                return __listStrings[key] + '(B)'
            elif AxisForBaseCenteredSymmetry == 'B':
                return __listStrings[key] + '(C)'
            elif AxisForBaseCenteredSymmetry == 'C':
                return __listStrings[key] + '(A)'
        elif key in __listStrings.keys():
            return __listStrings[key]
        return ""

    class ConstantType:
        # Matrices for centering.
        h_F = np.array ([[1,1, 0],[1,-1,0],[1,1,2]])
        h_I = np.array ([[1,1,-1],[1,-1,0],[0,0,1]]) # = 2*(h_F^{-1})^T
        h_R = np.array ([[-1,0,1],[0,1,-1],[-1,-1,-1]])
        h_B = np.array ([[ 1,1,0],[1,-1,0],[ 0, 0, 1]])
        # Matrices for error-stable Bravais lattice determination.
        CmP = [np.identity(3, dtype="int"), permute_matrix([0,2,1]), permute_matrix([1,0,2])]
        CmC = { True : mlist.CmC,                   # 69 matrices for an exhaustive search
                False: mlist.CmC[:mlist.len_CmC_1]  # 21 matrices for a quick search
                                        }
        CoF = [h_F, h_F.dot (permute_matrix([0,2,1])), h_F.dot (permute_matrix([1,2,0]))]
        CoI = [h_I, h_I.dot (permute_matrix([0,2,1])), h_I.dot (permute_matrix([1,2,0]))]
        ChR = { True : mlist.ChR,                   # 64 matrices for an exhaustive search
                False: mlist.ChR[: mlist.len_ChR_1] # 16 matrices for a quick search
                          }

    class InputType:
        def __init__ (self,):
            self.Sobs = np.identity(3, dtype="float")
            self.eps = 0.2
            self.DoesPrudentSearch = True
            self.AxisForBaseCenteredSymmetry = 'B'
            self.AxisForRhombohedralSymmetry = 'Hexagonal'

    class ResultType:
        """ 
            BravaisClasses: 14 lists of solutions with the Bravais-type symmetry 
            listNames     : Names of 14 Bravais types 
        """
        def __init__ (self, listNames):
            self.BravaisClasses = {}
            for name in listNames:
                self.BravaisClasses[name] = list()
    
        #----------------------------------------------
        #   Print out functions
        #----------------------------------------------
        def toText(self, Input):
            text = ''
            for name in self.BravaisClasses.keys():
                text += '<< {} : length {} >>'.format(BravaisLatticeDetermination.key2str(name, Input.AxisForBaseCenteredSymmetry), len (self.BravaisClasses[name])) + '\n'
                for arr1, arr2, distance in self.BravaisClasses[name]:
                    arr3 = arr1.dot(Input.Sobs).dot(arr1.T)
                    text += '{},   {},   {},   {}\n'.format(
                                    'g', 'Projection of g.S.g^T', 'g.S.g^T',
                                    'dist=' + str (distance))
                    for i in range (3):
                        text += '( ' + ' '.join ([str (i) for i in arr1[i]]) + ' )　　'
                        text += '( ' +' '.join (['{:.4f}'.format(v) for v in arr2[i]]) + ' )　　'
                        text += '( ' +' '.join (['{:.4f}'.format(v) for v in arr3[i]]) + ' )'
                        text += '\n'
                    text += '\n'
            return text
    
    
    @classmethod
    def primitive_monoclinic (cls, S_obs, gb, eps):
        #--------------------------------------------
        # Primitive monoclinic
        #--------------------------------------------
        ans = []
        for g in cls.ConstantType.CmP:
            g2 = g.dot(gb)
            S_new = g2.dot (S_obs).dot (g2.T)
            S = np.array ([[S_new[0,0],          0, S_new[0,2]],
                           [         0, S_new[1,1],          0],
                           [S_new[0,2],          0, S_new[2,2]]])
            if check_equiv (S_new, S, eps):
                # run the Gauss algorithm for S2.
                S2 = np.array ([[S[0,0], S[0,2]],
                                [S[2,0], S[2,2]]])
                h2, S2 = gauss_algorithm (S2)
                h = np.array ([[h2[0,0], 0, h2[0,1]],
                               [      0, 1,       0],
                               [h2[1,0], 0, h2[1,1]]])
                ans.append ([h.dot (g2), h.dot (S).dot (h.T)])
        return ans
    
    
    @classmethod
    def base_centered_monoclinic (cls, S_obs, gd, eps,
                                doesPrudentSearch = True, AxisForBaseCenteredSymmetry = 'B'):
        #--------------------------------------------
        # base-centered monoclinic
        #--------------------------------------------
        assert AxisForBaseCenteredSymmetry in ['A', 'B', 'C']
        ans = []
        for g in cls.ConstantType.CmC[doesPrudentSearch]:
            g2 = g.dot(gd)
            S_new = g2.dot (S_obs).dot (g2.T)
            S = np.array ([[S_new[0,0],          0, S_new[0,2]],
                           [         0, S_new[1,1],          0],
                           [S_new[0,2],          0, S_new[2,2]]])
            if check_equiv (S_new, S, eps):
                # run the Gauss algorithm for S2.
                S2 = np.array ([[S[0,0], S[0,2]],
                                [S[2,0], S[2,2]]])
                h2, S2 = gauss_algorithm (S2)
                if h2[1,0] % 2 != 0: # if the (2,1)-entry is odd
                    if h2[0,0] % 2 != 0: # if the (1,1)-entry is odd
                        h2 = np.array ([[1,0],[-1,-1]]).dot (h2)
                    else:
                        h2 = np.array ([[0,1],[1,0]]).dot (h2)
                h = np.array ([[h2[0,0], 0, h2[0,1]],
                               [      0, 1,       0],
                               [h2[1,0], 0, h2[1,1]]])
                if AxisForBaseCenteredSymmetry == 'A':
                    h = permute_matrix([1,2,0]).dot (h)
                elif AxisForBaseCenteredSymmetry == 'C':
                    h = permute_matrix([2,0,1]).dot (h)
                ans.append ([h.dot (g2), h.dot (S).dot (h.T)])
        return ans
    
    
    @classmethod
    def face_centered_orthorhombic (cls, S_obs, gd, eps):
        #--------------------------------------------
        # Face centered Orthorhombic (面心直方格子)
        #--------------------------------------------
        ans = []
        for g in cls.ConstantType.CoF:
            g2 = g.dot (gd)
            S_new = g2.dot (S_obs).dot (g2.T)
            # S:= [[s11,0,0],[0,s22,0],[0,0,s33]] using Snew = (sij)
            s11, s22, s33 = np.diag (S_new)
            S = np.array ([[s11,   0,   0],
                           [  0, s22,   0],
                           [  0,   0, s33]])
            # if dist (S_new, S) < eps:
            if check_equiv (S_new, S, eps):
                # Sort the diagonal entries of S, replacing the rows of g
                indices = np.argsort (np.diag (S))
                gp = permute_matrix(indices)
                ans.append ([g2 [indices, :], gp.dot(S).dot(gp.T)])
        return ans
    
    
    @classmethod
    def body_centered_orthorhombic (cls, S_obs, gd2, eps, AxisForBaseCenteredSymmetry = 'B'):
        #--------------------------------------------
        # Body-centered Orthorhombic (体心直方格子)
        #--------------------------------------------
        ans = []
        for g in cls.ConstantType.CoI:
            g2 = g.dot (gd2)
            S_new = g2.dot (S_obs).dot (g2.T)
            # S:= [[s11,0,0],[0,s22,0],[0,0,s33]] using Snew = (sij)
            s11, s22, s33 = np.diag (S_new)
            S = np.array ([[s11,   0,   0],
                           [  0, s22,   0],
                           [  0,   0, s33]])
            if check_equiv (S_new, S, eps):
                # Sort the diagonal entries of S, replacing the rows of g
                indices = np.argsort (np.diag (S))
                gp = permute_matrix(indices)
                ans.append ([g2 [indices, :], gp.dot(S).dot(gp.T)])
        return ans
    
    @classmethod
    def rhombohedral_centering (cls, S_obs, gd, eps,
                                doesPrudentSearch = True,
                                AxisForRhombohedralSymmetry = 'Hexagonal'):
        #--------------------------------------------
        # Output: Ans_hR for Rhombohedral (菱面体格子)
        #--------------------------------------------
        assert AxisForRhombohedralSymmetry in ['Hexagonal', 'Rhombohedral']
        ans_hR = []
        for g in cls.ConstantType.ChR[doesPrudentSearch]:
            g2 = g.dot (gd)
            S_new = g2.dot (S_obs).dot (g2.T)
            # a := (s11 + s22 + s33)/3,
            # d := (s12 + s13 + s23)/3,
            # using the (i; j)-entry sij of S_new.
            s11, s22, s33 = np.diag (S_new)
            s12, s13, s23 = S_new[0,1], S_new[0,2], S_new[1,2]
            a = (s11 + s22 + s33) / 3
            d = (s12 + s13 + s23) / 3
            # S = [[a,d,d],[d,a,d],[d,d,a]]
            S = np.array ([[a, d, d],[d, a, d],[d, d, a]])
            if check_equiv (S_new, S, eps):
                if AxisForRhombohedralSymmetry == 'Hexagonal':
                    h_R = cls.ConstantType.h_R
                    ans_hR.append ([h_R.dot (g2), h_R.dot (S).dot (h_R.T)])
                else:
                    ans_hR.append ([g2, S])
        return ans_hR
    
    
    #----------(3) Projection to higher symmetries:-----------
    @staticmethod
    def primitive_monoclinic_to_orthorhombic (mlist, eps):
        #--------------------------------------------
        # Monoclinic --> Orthorhombic 
        #--------------------------------------------
        ans = []
        for g, Sm in mlist:
            s11, s22, s33 = np.diag (Sm)
            So = np.array ([[s11,   0,   0],
                            [  0, s22,   0],
                            [  0,   0, s33]])
            if check_equiv (Sm, So, eps):
                #Sort the diagonal entries of So
                indices = np.argsort (np.diag (So))
                gp = permute_matrix(indices)
                g_new = g [indices, :] # g_new = g.gp
                ans.append ([g_new, gp.dot (So).dot (gp.T)])
        return ans
    
    @staticmethod
    def base_monoclinic_to_orthorhombic (mlist, eps, AxisForRhombohedralSymmetry):
        #--------------------------------------------
        # Monoclinic --> Orthorhombic 
        #--------------------------------------------
        ans = []
        for g, Sm in mlist:
            s11, s22, s33 = np.diag (Sm)
            So = np.array ([[s11,   0,   0],
                            [  0, s22,   0],
                            [  0,   0, s33]])
            if check_equiv (Sm, So, eps):
                # Transeform the base axis to c.
                if AxisForRhombohedralSymmetry == 'A' or AxisForRhombohedralSymmetry == 'C':
                    if AxisForRhombohedralSymmetry == 'A': # base-axis = b -> c 
                        gp = permute_matrix([2,0,1])
                    else: # base-axis = a -> c
                        gp = permute_matrix([1,2,0])
                    g_new = gp.dot (g)
                    So = gp.dot (So).dot (gp.T)
                else: 
                    g_new = g
                # So(0,0) <= So(1,1)
                if So[0,0] > So[1,1]:
                    gp = permute_matrix([1,0,2])
                    ans.append ([gp.dot (g_new), gp.dot (So).dot (gp.T)])
                else:
                    ans.append ([g_new, So])
        return ans

   
    @staticmethod
    def orthorhombic_to_tetragonal (mlist, eps):
        #--------------------------------------------
        # Orthorhombic --> Tetragonal
        #--------------------------------------------
        ans = []
        for g, So in mlist:
            for i, j, k in [(0,1,2), (0,2,1), (1,2,0)]:
                s = (So[i, i] + So[j, j]) / 2
                St = np.array ([[s, 0,       0],
                                [0, s,       0],
                                [0, 0, So[k,k]]])
                # gnew := t(vi vj vk), using the rows of g = t(v1 v2 v3).
                gp = permute_matrix([i,j,k])
                g_new = g[[i, j, k], :] # g_new = gp.g
                if check_equiv (gp.dot (So).dot (gp.T), St, eps):
                    ans.append ([g_new, St])
        return ans
    
    
    @staticmethod
    def primitive_monoclinic_to_hexagonal (ans_mP, eps):
        #--------------------------------------------
        # Primitive Monoclinic --> Hexagonal 
        #--------------------------------------------
        ans_hP = []
        for g, SmP in ans_mP:
            s11, s22, s33 = np.diag (SmP)
            s13 = SmP[0,2]
            s = s11 + s33 + 2 * s13
            ShP = np.array ([[   s, -s/2,   0],
                             [-s/2,    s,   0],
                             [   0,    0, s22]])
            # gp = ([[1,0,0],[0,0,1],[0,1,0]])
            gp = permute_matrix([0,2,1])
            g_new = gp.dot (g)
            if check_equiv (gp.dot(SmP).dot (gp.T), ShP, eps):
                ans_hP.append ([g_new, ShP])
        return ans_hP
    
    
    @staticmethod
    def orthorhombic_to_cubic (mlist, eps):
        #--------------------------------------------
        # Orthorhombic --> Cubic
        #--------------------------------------------
        ans = []
        for g, So in mlist:
            s = np.mean (np.diag (So))
            Sc = np.array ([[s,0,0],[0,s,0],[0,0,s]])
            if check_equiv (So, Sc, eps):
                ans.append ([g, Sc])
        return ans
    
    
    def __init__ (self,):
        """ Initialize self.result as 14 empty lists. """
        self.result = BravaisLatticeDetermination.ResultType (BravaisLatticeDetermination.__listNames)
        assert BravaisLatticeDetermination.__listNames == list (self.result.BravaisClasses.keys())

    def set_bravais_class (self, Input):
        """ All the results of this method are stored in self.result. """
        """ self.resultに14個の行列リストをセットする. """
        #--------------------------------------------
        # matrix dict {name : matrix list} of:
        # names :
        #    ans_tP (Triclinic)
        #    ans_mP (Primitive Monoclinic)
        #    ans_mC (Base-centered Monoclinic)
        #    ans_oP (Primitive Orthorhombic)
        #    ans_oC (Base-centered Orthorhombic)
        #    ans_oI (Body-centered Orthorhombic)
        #    ans_oF (Face-centered Orthorhombic)
        #    ans_tP (Primitive Tetragonal)
        #    ans_tI (Body-centered Tetragonal)
        #    ans_hR (Rhombohedral)
        #    ans_hP (Hexagonal)
        #    ans_cP (Primitive Cubic)
        #    ans_cI (Body-centered Cubic)
        #    ans_cF (Face-centered Cubic)
        #--------------------------------------------
        
        Sobs = Input.Sobs
        eps  = Input.eps
        doesPrudentSearch = Input.DoesPrudentSearch
        AxisForBaseCenteredSymmetry = Input.AxisForBaseCenteredSymmetry
        AxisForRhombohedralSymmetry = Input.AxisForRhombohedralSymmetry
        # S_buer = gb.Sobs.gb^T is Buerger-reduced.
        # S_del  = gd.Sobs.gd^T is Delaunay-reduced.
        # The inverse of S_del_inv = gd2.Sobs.gd2^T is Delaunay-reduced.
        gb, S_buer = Buerger_reduction (Sobs) 
        gd, S_del  = Delaunay_reduction (S_buer); gd = gd.dot(gb)
        gd2,S_del_inv = Delaunay_reduction_of_inverse (S_buer); gd2 = gd2.dot(gb)

        ## Triclinic
        self.result.BravaisClasses['triclinic'].append([gb, S_buer])

        # --- (primitive / base-centered) monoclinic ---
        ## Primitive monoclinic
        self.result.BravaisClasses['primitive monoclinic'] = self.primitive_monoclinic (Sobs, gb, eps)
        ## Base-centered monoclinic 
        self.result.BravaisClasses['base-centered monoclinic'] = self.base_centered_monoclinic (Sobs, gd, eps,doesPrudentSearch, AxisForBaseCenteredSymmetry)
        
        # --- (face / body-centered) orthorhombic ---
        ## Face-centered orthorhombic
        self.result.BravaisClasses['face-centered orthorhombic'] = self.face_centered_orthorhombic (Sobs, gd, eps)
        ## Body-centered orthorhimbic
        self.result.BravaisClasses['body-centered orthorhombic'] = self.body_centered_orthorhombic (Sobs, gd2, eps)
        
        ## Rhombohedral
        self.result.BravaisClasses['rhombohedral'] = self.rhombohedral_centering (Sobs, gd, eps, doesPrudentSearch, AxisForRhombohedralSymmetry)

        # --- (Monoclinic --> Orthorhombic) ---
        ## Primitive Monoclinic --> Primitive orthorhombic
        self.result.BravaisClasses['primitive orthorhombic'] = self.primitive_monoclinic_to_orthorhombic (self.result.BravaisClasses['primitive monoclinic'], eps)
        # Base-centered Monoclinic --> Base-centered orthorhombic
        self.result.BravaisClasses['base-centered orthorhombic'] = self.base_monoclinic_to_orthorhombic (self.result.BravaisClasses['base-centered monoclinic'], eps, AxisForRhombohedralSymmetry)
        
        # --- (Orthorhombic --> Tetragonal) ---
        ## Primitive Orthorhombic --> Primitive tetragonal
        self.result.BravaisClasses['primitive tetragonal'] = self.orthorhombic_to_tetragonal (self.result.BravaisClasses['primitive orthorhombic'], eps)
        ## Body-centered orthorombic --> Body-centered tetragonal
        self.result.BravaisClasses['body-centered tetragonal'] = self.orthorhombic_to_tetragonal (self.result.BravaisClasses['body-centered orthorhombic'], eps)
        
        # --- (Primitive Monoclinic --> Hexagonal) ---
        ## Primitive Monoclinic --> Hexagonal
        self.result.BravaisClasses['hexagonal'] = self.primitive_monoclinic_to_hexagonal (self.result.BravaisClasses['primitive monoclinic'], eps)
        
        # --- (Orthorhombic --> Cubic) ---
        ## Primtive orthorhombic --> Primitive cubic
        self.result.BravaisClasses['primitive cubic'] = self.orthorhombic_to_cubic (self.result.BravaisClasses['primitive orthorhombic'], eps)
        ## Body-centered orthorhombic --> Body-centered cubic
        self.result.BravaisClasses['body-centered cubic'] = self.orthorhombic_to_cubic (self.result.BravaisClasses['body-centered orthorhombic'], eps)
        ## Face-centered orthorhombic -->  Face-centered cubic
        self.result.BravaisClasses['face-centered cubic'] = self.orthorhombic_to_cubic (self.result.BravaisClasses['face-centered orthorhombic'], eps)
        assert BravaisLatticeDetermination.__listNames == list (self.result.BravaisClasses.keys())
        
        # Set distances in the 3rd entries and sort.
        for _, list_btype in self.result.BravaisClasses.items():               
            for entries in list_btype:
                g = entries[0]
                S = entries[1]
                entries.append(dist(S, g.dot (Sobs).dot (g.T)))
            list_btype.sort(key=lambda x: x[2])


if __name__ == '__main__':
    ndim = 3
    # Make a positive-definite matrix S from a basis matrix B.
    B = [np.random.random() for _ in range (ndim*ndim)]
    B = np.array (B).reshape(ndim,ndim)
    S_input = B.dot (B.T)
    print ("* Input S:")
    print (S_input)

    bc = BravaisLatticeDetermination ()
    bc_input = BravaisLatticeDetermination.InputType()
    
    bc_input.Sobs = S_input
    bc_input.eps  = 0.2
    bc_input.DoesPrudentSearch = False
    bc_input.AxisForBaseCenteredSymmetry = 'A'
    bc_input.AxisForRhombohedralSymmetry = 'Hexagonal'
    
    bc.set_bravais_class (bc_input)
    print( bc.result.toText(bc_input) )
