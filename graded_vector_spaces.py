#!/usr/bin/python

import numpy as np
from vector_spaces import *
import rings as ri
import rationals as QQ

class gvs():
    '''The class of graded vector spaces
    
    Usage:
        V=gvs(k)
        k [ri.field]: field of definition of V
    
    Attributes:
        V.base [ri.field]: field of definition of V
        V.gps [dictionary]: graded pieces of V
            { n [int] : W [vs] - the nth graded piece of V }
    
    Methods:
        V^n:
            Args: n [int]
            Returns: [vs] The nth graded piece of V.
        V*W, V.otimes(W):
            Args: V, W [gvs]
            Returns: [gvs] Tensor product of V and W.
        V+W, V.oplus(W):
            Args: V, W [gvs]
            Returns: [gvs] Direct sum of V and W.
        V.shift(n=1):
            Args: n [int] default=1
            Returns: [gvs] Shift of V by n in the grading.
        print(V):
            Returns: A string showing the dictionary
                { n : (V^n).dim }
        V.setgps(D):
            Args: D [dictionary]
                { n [int] : either dim(V^n) [int] or V^n [vs] }
            Modifies V.gps to incorporate D.
    '''
    def __init__(self,k):
        '''Instantiates a zero graded vector space over k'''
        self.base=k
        self.gps={}

    def __xor__(self,n):
        ''' Accesses the graded piece V^n '''
        if n in self.gps:
            return self.gps[n]
        else:
            return vs(self.base,0)

    def otimes(self,other):
        '''Tensor product of graded vector spaces:

        Note that if e^m_1,...,e^m_{k_m} is a basis of V^m and
        f^n_1,...,f^n_{l_n} is a basis of W^n
        then this tensor product will output the basis in the order
        e^0_1*f^{n+m}_1, e^0_1*f^{n+m}_2,...e^0_1*f^{n+m}_{l_{n+m}},
        e^0_2*f^{n+m}_1, ...
        ...
        e^0_{k_0}*f^{n+m}_1,...
        e^1_1*f^{n+m-1}_1,...
        '''
        k=self.base
        new_gvs=gvs(k)
        V=self.gps
        W=other.gps
        for m in V:
            for n in W:
                if m+n not in new_gvs.gps:
                    new_gvs.gps[m+n]=V[m]*W[n]
                else:
                    new_gvs.gps[m+n]+=V[m]*W[n]
        return new_gvs

    def __mul__(self,other):
        '''Tensor product of graded vector spaces.'''
        return self.otimes(other)
    
    def oplus(self,other):
        '''Direct sum of graded vector spaces'''
        new_gvs=self
        for n in other.gps:
            if n in new_gvs.gps:
                new_gvs.gps[n]+=other.gps[n]
            else:
                new_gvs.gps[n]=other.gps[n]
        return new_gvs

    def __add__(self,other):
        '''Direct sum of graded vector spaces'''
        return self.oplus(other)
    
    def __truediv__(self,other):
        'Quotient graded vector space'
        k=self.base
        V=self.gps
        W=other.gps
        new_gvs=gvs(k)
        for n in V:
            new_gvs.gps[n]=(V[n])/(W[n])

    def shift(self,n=1):
        'Shift a graded vector space by n; default shift is 1'
        k=self.base
        new_gvs=gvs(k)
        for m in self.gps:
            new_gvs[n+m]=self.gps[n]
        return new_gvs

    def __str__(self):
        'Printing a graded vector space'
        graded_dimensions={}
        for n in self.gps:
            graded_dimensions[n]=self.gps[n].dim
        return "%s" % str(graded_dimensions)

    def setgps(self,D):
        '''V.setgps(D) modifies V.gps to incorporate D.'''
        for n in D:
            if type(D[n]) is int:
                self.gps[n]=vs(k,D[n])
            elif type(D[n]) is vs:
                if D[n].base==k:
                    self.gps[n]=D[n]
                else:
                    print("Graded pieces over different fields")
            else:
                print("Tried to create a graded vector space out of bananas.")

            
        V1,W1=A.source,A.target
        V2,W2=D.source,D.target
        V1_test,W2_test=C.source,C.target
        V2_test,W1_test=B.source,B.target
        space_test=(V1==V1_test and V2==V2_test and \
                     W1==W1_test and W2==W2_test)

        if space_test and degree_test:
            block_ghomo=ghomo(d,V1+V2,W1+W2)
            M=np.asarray()

    def Id(self):
        '''Returns the identity homomorphism for this
        graded vector space.'''
        identity_ghomo=ghomo(0,self,self)
        for n in self.gps:
            identity_homomorphism.gps[n]=(self^n).Id()
        return identity_ghomo
    
class ghomo():
    '''The class of graded maps between graded vector spaces.

    F=ghomo(d,V,W) creates an empty graded map of degree d
    between V and W. Its graded pieces can be set with
    F.gps(M), where M is a dictionary
    { n [int] : F^n:V^n-->W^{n+d} [homo] or [np.ndarray] }

    Attributes:
        F.source [gvs]: domain of F
        F.target [gvs]: range of F
        F.gps [dictionary]:
            { n [int] : F.matrix[n] [homo] }
        F.base [ri.field]: field of definition of F

    Methods:
        F*G:
            Args: F, G [homo]
            Returns: [homo] Composition of F with G
        F+G:
            Args: F, G [homo] with same .source and .target
            Returns: [homo] Sum of F with G in hom(source,target)
        F^n:
            Args: n [int] Returns F.gps[n]
        F.oplus(G):
            Args: G [homo]
            Returns: [homo] Block matrix F\oplus G
        F.otimes(G):
            Args: G [homo]
            Returns: [homo] Tensor/Kronecker product F\otimes G
        F.ker():
            Returns: [vs] Kernel of F
        F.image():
            Returns: [vs] Image of F
    '''

    def __init__(self,d,V,W):
        self.base=V.base
        if V.base!=W.base:
            print("Warning: created a homomorphism, ",F,\
                  "  between vector spaces over different fields.")
        self.source,self.target=V,W
        self.degree=d
        self.gps={}

    def __xor__(self,n):
        ''' Accesses the graded piece F^n '''
        k=self.source.base
        if n in self.gps:
            return self.gps[n]
        else:
            M=k.mod(np.zeros(shape=(self.target.dim,self.source.dim)))
            return homo(0,self.source,self.target)

    def __mul__(self,other):
        'Compose graded linear maps'
        d=self.degree
        e=other.degree
        V=self.source
        W=other.target
        new_ghomo=ghomo(d+e,V,W)
        M=self.gps
        N=other.gps
        for n in M:
            if n+d in N:
                # We only bother to keep track of maps
                # between nonzero vector spaces
                new_ghomo.gps[n]=(self^n)*(other^(n+d))
        return new_ghomo
    
    def __add__(self,other):
        'Add linear maps'
        src1,src2=self.source,other.source
        tar1,tar2=self.target,other.target
        d,e=self.degree,other.degree
        if d!=e:
            print("Cannot add graded linear maps with different degrees.")
        else:
            new_ghomo=ghomo(d,src_1,tar1)
            if src1==src2 and tar1==tar2:
                new_ghomo.gps=self.gps
                for n in other.gps:
                    if n in new_ghomo.gps:
                        new_ghomo.gps[n]+=other^n
                    else:
                        new_ghomo.gps[n]=other^n
                return new_ghomo
                else:
                    print("Cannot add graded linear maps with ",\
                          "different sources/targets.")

    def oplus(self,other):
        '''Form a block diagonal matrix of graded linear maps'''
        d=self.degree
        B=ghomo(d,other.source,self.target)
        C=ghomo(d,self.source,other.target)
        return ghomo.block(self,B,C,other)
      
    def otimes(self,other):
        'Tensor/Kronecker/outer product f\otimes g of graded linear maps'
        d=self.degree
        e=other.degree
        A=self.gps
        B=other.gps
        new_source=self.source*other.source
        new_target=self.target*other.target
        new_ghomo=ghomo(d,new_source,new_target)
        for m in A:
            for n in B:
                addendum=A[m].otimes(B[n])
                if m+n not in new_ghomo.gps:
                    new_ghomo.gps[m+n]=addendum
                else:
                    new_ghomo.gps[m+n]=new_ghomo^(m+n).oplus(addendum)
        return new_ghomo

    def ker(self):
        'Kernel of graded linear maps'
        new_gvs=gvs(self.base)
        for n in self.gps:
            new_gvs.gps[n]=(self^n).ker()
        return new_gvs

    def image(self):
        'Image of graded linear maps'
        new_gvs=gvs(self.base)
        d=self.degree
        for n in self.gps:
            new_gvs.gps[n+d]=(self^n).image()
        return new_gvs
            
    def cohomology(self):
        'Cohomology of a graded linear map, that is kernel/image'
        return (self.ker())/(self.image())

    def setgps(self,M):
        for n in M:
            if type(M[n]) is np.ndarray:
                self.gps[n]=homo(V^n,W^(n+d),M[n])
            elif type(M[n]) is homo:
                self.gps[n]=M[n]
            else:
                print("ghomo requires final argument ",\
                      "to be of type np.array or homomomorphism.")

    def rejig_1(self,n=1):
        '''This implements the canonical isomorphism
        hom(M,N)[n] --> hom(M,N[n])
        '''
        new_ghomo=ghomo(self.degree-n,self.source,self.target.shift(n))
        new_ghomo.gps=self.gps
        return new_ghomo
    
    def rejig_2(self,n=1):
        '''This implements the canonical isomorphism
        hom(M,N)[d]-->hom(M[n],N[n])[d]
        f--------->((-1)**(n(d-1))) f
        '''
        new_ghomo=ghomo(self.degree,self.source.shift(n),self.target.shift(n))
        for n in self.gps:
            new_ghomo.gps[n]=self.base.mod((-1)**(n*(self.degree-1)))*self.gps[n]
            # again, we need to make sure ndarrays can cope with this
        return new_ghomo
    
    def rejig_3(self,n=1):
        '''This implements the canonical isomorphism
        hom(M,N)-->hom(M[k],N)[k]
        Remember there are two different ways to implement this
        which differ by a sign - the other is rejig_2(rejig_1(self,-n),n)
        '''
        return rejig_1(rejig_2(self,n),-n)

    def koszulify(self,m=1):
        new_ghomo=ghomo(self.degree,self.source,self.target)
        for n in self.gps:
            new_ghomo.gps[n]=self.base.mod((-1)**(m*(n-1)))*self.gps[n]
            # Hopefully np.arrays can cope with this
        return new_ghomo
    
    @classmethod
    def block(cls,A,B,C,D):
        '''Returns a block matrix from four compatible
        graded homomorphisms.'''
        d,d_other=A.degree,(B.degree,C.degree,D.degree)
        if (x==d for x in d_other):
            block_ghomo=ghomo(d,A.source+D.source,A.target+D.target)
            for n in (A.gps or B.gps or C.gps or D.gps):
                block_ghomo.gps[n]=homo.block(A^n,B^n,C^n,D^n)
            return block_ghomo
        else:
            print("Can't add graded homomorphisms of different degrees.")


                
'''Some random tests
Q=QQ.QQ()
V=gvs(Q)
V.setgps({0:2,3:2})
print(((V.oplus(V.otimes(V)))^3).dim)
N=np.array([[1,1],[0,0]])
M={3:N}
F=ghomo(0,V,V)
F.gps=M
print((F.otimes(F)).ker())
'''
