#!/usr/bin/python

import numpy as np
import vector_spaces as vs
import rings as ri
import rationals as QQ

class graded_vector_space():
    '''The class of graded vector spaces
    
    Usage:
        V=graded_vector_space(k,D)
        k [ri.field]: field of definition of V
        D [dictionary]:
            { n [int] : m [int] or V [vs.vector_space] }
        D[n] is either the nth graded piece of V or else
        its dimension.
    
    Attributes:
        V.base [ri.field]: field of definition of V
        V.graded_pieces [dictionary]:
            { n [int] : W [vs.vector_space] }
    
    Methods:
        V^n:
            Args: n [int]
            Returns: [vs.vector_space] The nth graded piece of V.
        V*W, V.otimes(W):
            Args: V, W [graded_vector_space]
            Returns: [graded_vector_space] Tensor product of V and W.
        V+W, V.oplus(W):
            Args: V, W [graded_vector_space]
            Returns: [graded_vector_space] Direct sum of V and W.
        V.shift(n=1):
            Args: n [int] default=1
            Returns: [graded_vector_space] Shift of V by n in the grading.
        print(V):
            Returns: A string showing the dictionary
                { n : (V^n).dim }
    '''
    def __init__(self,k,D):
        self.base=k
        self.graded_pieces={}
        for n in D:
            if type(D[n]) is int:
                self.graded_pieces[n]=vs.vector_space(k,D[n])
            elif type(D[n]) is vs.vector_space:
                if D[n].base==k:
                    self.graded_pieces[n]=D[n]
                else:
                    print("Graded pieces over different fields")
            else:
                print("Tried to create a graded vector space out of bananas.")

    def __xor__(self,n):
        ''' Accesses the graded piece V^n '''
        if n in self.graded_pieces:
            return self.graded_pieces[n]
        else:
            return vs.vector_space(self.base,0)

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
        new_graded_vector_space=graded_vector_space(k,{})
        V=self.graded_pieces
        W=other.graded_pieces
        for m in V:
            for n in W:
                if m+n not in new_graded_vector_space.graded_pieces:
                    new_graded_vector_space.graded_pieces[m+n]=V[m]*W[n]
                else:
                    new_graded_vector_space.graded_pieces[m+n]+=V[m]*W[n]
        return new_graded_vector_space

    def __mul__(self,other):
        '''Tensor product of graded vector spaces.'''
        return self.otimes(other)
    
    def oplus(self,other):
        '''Direct sum of graded vector spaces'''
        new_graded_vector_space=self
        for n in other.graded_pieces:
            if n in new_graded_vector_space.graded_pieces:
                new_graded_vector_space.graded_pieces[n]+=other.graded_pieces[n]
            else:
                new_graded_vector_space.graded_pieces[n]=other.graded_pieces[n]
        return new_graded_vector_space
    
    def __add__(self,other):
        '''Direct sum of graded vector spaces'''
        return self.oplus(other)
    
    def __truediv__(self,other):
        'Quotient graded vector space'
        k=self.base
        V=self.graded_pieces
        W=other.graded_pieces
        new_graded_vector_space=graded_vector_space(k,{})
        for n in V:
            new_graded_vector_space.graded_pieces[n]=(V[n])/(W[n])

    def shift(self,n=1):
        'Shift a graded vector space by n; default shift is 1'
        k=self.base
        new_graded_vector_space=graded_vector_space(k,{})
        for m in self.graded_pieces:
            new_graded_vector_space[n+m]=self.graded_pieces[n]
        return new_graded_vector_space

    def __str__(self):
        'Printing a graded vector space'
        graded_dimensions={}
        for n in self.graded_pieces:
            graded_dimensions[n]=self.graded_pieces[n].dim
        return "%s" % str(graded_dimensions)

    @classmethod
    def zero(cls,k):
        return graded_vector_space(k,{})
    
class graded_homomorphism():
    '''The class of graded homomorphisms between graded vector spaces.

    F=homomorphism(d,V,W,M) creates a dictionary of linear maps
    F^n:V^n-->W^{n+d} with matrix M[n] [np.ndarray].
    Here V, W are of type [graded_vector_space].

    Attributes:
        F.source [graded_vector_space]: domain of F
        F.target [graded_vector_space]: range of F
        F.graded_maps [dictionary]:
            { n [int] : F.matrix[n] [vs.homomorphism] }
        F.base [ri.field]: field of definition of F

    Methods:
        F*G:
            Args: F, G [homomorphism]
            Returns: [homomorphism] Composition of F with G
        F+G:
            Args: F, G [homomorphism] with same .source and .target
            Returns: [homomorphism] Sum of F with G in hom(source,target)
        F^n:
            Args: n [int] Returns F.graded_maps[n]
        F.oplus(G):
            Args: G [homomorphism]
            Returns: [homomorphism] Block matrix F\oplus G
        F.otimes(G):
            Args: G [homomorphism]
            Returns: [homomorphism] Tensor/Kronecker product F\otimes G
        F.ker():
            Returns: [vector_space] Kernel of F
        F.image():
            Returns: [vector_space] Image of F
    '''

    def __init__(self,d,V,W,M):
        self.base=V.base
        if V.base!=W.base:
            print("Warning: created a homomorphism, ",F,\
                  "  between vector spaces over different fields.")
        self.source,self.target=V,W
        self.degree=d
        self.graded_maps={}
        for n in M:
            if type(M[n]) is np.ndarray:
                self.graded_maps[n]=vs.homomorphism(V^n,W^(n+d),M[n])
            elif type(M[n]) is vs.homomorphism:
                self.graded_maps[n]=M[n]
            else:
                print("graded_homomorphism requires final argument ",\
                      "to be of type np.array or vs.homomomorphism.")

    def __mul__(self,other):
        'Compose graded linear maps'
        d=self.degree
        e=other.degree
        V=self.source
        W=other.target
        new_graded_homomorphism=graded_homomorphism(d+e,V,W,{})
        M=self.graded_maps
        N=other.graded_maps
        for n in M:
            if n+d in N:
                # We only bother to keep track of maps
                # between nonzero vector spaces
                new_graded_homomorphism.graded_maps[n]=M[n]*N[n+d]
        return new_graded_homomorphism
    
    def __add__(self,other):
        'Add linear maps'
        src1,src2=self.source,other.source
        tar1,tar2=self.target,other.target
        d,e=self.degree,other.degree
        if d!=e:
            print("Cannot add graded linear maps with different degrees.")
        else:
            new_graded_homomorphism=graded_homomorphism(d,src_1,tar1,{})
            if src1==src2 and tar1==tar2:
                new_graded_homomorphism.graded_maps=self.graded_maps
                for n in other.graded_maps:
                    if n in new_graded_homomorphism.graded_maps:
                        new_graded_homomorphism.graded_maps[n]+=other.graded_maps[n]
                    else:
                        new_graded_homomorphism.graded_maps[n]=other.graded_maps[n]
                        return new_graded_homomorphism
                else:
                    print("Cannot add graded linear maps with ",\
                          "different sources/targets.")

    def oplus(self,other):
        '''Form a block matrix of graded linear maps'''
        d=self.degree
        if d!=other.degree:
            print("Can't direct sum maps of different degrees!")
        else:
            new_source=self.source+other.source
            new_target=self.target+other.target
            new_graded_homomorphism=graded_homomorphism(d,new_source,\
                                                        new_target,self.graded_maps)
            for n in other.matrix:
                if n in new_graded_homomorphism:
                    new_graded_homomorphism[n]+=other.graded_maps[n]
                else:
                    new_graded_homomorphism[n]=other.graded_maps[n]
            return new_graded_homomorphism

    def __xor__(self,n):
        ''' Accesses the graded piece F^n '''
        if n in self.graded_maps:
            return self.graded_maps[n]
        else:
            return vs.homomorphism(0,self.source,self.target,{})
        
    def otimes(self,other):
        'Tensor/Kronecker/outer product f\otimes g of graded linear maps'
        #Still needs further testing
        d=self.degree
        e=other.degree
        A=self.graded_maps
        B=other.graded_maps
        new_source=self.source*other.source
        new_target=self.target*other.target
        new_graded_homomorphism=graded_homomorphism(d,new_source,new_target,{})
        for m in A:
            for n in B:
                addendum=A[m].otimes(B[n])
                if m+n not in new_graded_homomorphism.graded_maps:
                    new_graded_homomorphism.graded_maps[m+n]=addendum
                else:
                    new_graded_homomorphism.graded_maps[m+n]+=addendum
        return new_graded_homomorphism

    def ker(self):
        'Kernel of graded linear maps'
        new_graded_vector_space=graded_vector_space(self.base,{})
        for n in self.graded_maps:
            new_graded_vector_space.graded_pieces[n]=self.graded_maps[n].ker()
        return new_graded_vector_space

    def image(self):
        'Image of graded linear maps'
        new_graded_vector_space=graded_vector_space(self.base,{})
        d=self.degree
        for n in self.graded_maps:
            new_graded_vector_space.graded_pieces[n+d]=self.graded_maps[n].image()
        return new_graded_vector_space

            
    def cohomology(self):
        'Cohomology of a graded linear map, that is kernel/image'
        return (self.ker())/(self.image())

    @classmethod
    def zero(cls,d,V,W):
        graded_homomorphism(d,V,W,{})

'''Some random tests
Q=QQ.QQ()
V=graded_vector_space(Q,{0:2,3:2})
print(((V.oplus(V.otimes(V)))^3).dim)
N=np.array([[1,1],[0,0]])
M={3:N}
F=graded_homomorphism(0,V,V,M)
print((F.otimes(F)).ker())
'''
