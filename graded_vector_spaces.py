#!/usr/bin/python

import numpy as np
import vector_spaces as vs

class gr_vs():
    '''The class of graded vector spaces
    
    Usage:
        V=gr_vs(D)
        D [dictionary]:
            { n [int] : m [int] or V [vs.vector_space] }
        D[n] is either the nth graded piece of V or else
        its dimension.
    
    Attributes:
        V.graded_pieces [dictionary]:
            { n [int] : W [vs.vector_space] }
    
    Methods:
        V^n:
            Args: n [int]
            Returns: [vs.vector_space] The nth graded piece of V.
        V*W, V.otimes(W):
            Args: V, W [gr_vs]
            Returns: [gr_vs] Tensor product of V and W.
        V+W, V.oplus(W):
            Args: V, W [gr_vs]
            Returns: [gr_vs] Direct sum of V and W.
        V.shift(n=1):
            Args: n [int] default=1
            Returns: [gr_vs] Shift of V by n in the grading.
        print(V):
            Returns: A string showing the dictionary
                { n : (V^n).dim }
    '''
    def __init__(self,D):
        self.graded_pieces={}
        for n in D:
            if type(D[n]) is int:
                self.graded_pieces[n]=vs.vector_space(D[n])
            elif type(D[n]) is vs.vector_space:
                self.graded_pieces[n]=D[n]
            else:
                print("Tried to create a graded vector space out of bananas.")
    def __xor__(self,n):
        ''' Accesses the graded piece V^n '''
        if n in self.graded_pieces:
            return self.graded_pieces[n]
        else:
            return vs.vector_space(0)
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
        new_gr_vs=gr_vs({})
        for m in self.graded_pieces:
            for n in other.graded_pieces:
                if m+n not in new_gr_vs.graded_pieces:
                    new_gr_vs.graded_pieces[m+n]=self.graded_pieces[m]*other.graded_pieces[n]
                else:
                    new_gr_vs.graded_pieces[m+n]+=self.graded_pieces[m]*other.graded_pieces[n]
        return new_gr_vs
    def __mul__(self,other):
        '''Tensor product of graded vector spaces.'''
        return self.otimes(other)
    def oplus(self,other):
        '''Direct sum of graded vector spaces'''
        new_gr_vs=self
        for n in other.graded_pieces:
            if n in new_gr_vs.graded_pieces:
                new_gr_vs.graded_pieces[n]+=other.graded_pieces[n]
            else:
                new_gr_vs.graded_pieces[n]=other.graded_pieces[n]
        return new_gr_vs
    def __add__(self,other):
        '''Direct sum of graded vector spaces'''
        return self.oplus(other)
    def __truediv__(self,other):
        'Quotient graded vector space'
        new_gr_vs={}
        for n in self.graded_pieces:
            new_gr_vs[n]=(self.graded_pieces[n])/(other.graded_pieces[n])
    def shift(self,n=1):
        'Shift a graded vector space by n; default shift is 1'
        new_gr_vs=gr_vs({})
        for m in self.graded_pieces:
            new_gr_vs[n+m]=self.graded_pieces[n]
        return new_gr_vs
    def __str__(self):
        'Printing a graded vector space'
        gr_dim={}
        for n in self.graded_pieces:
            gr_dim[n]=self.graded_pieces[n].dim
        return "%s" % str(gr_dim)
    
class zero_gr_vs(gr_vs):
    'A class of zero graded vector spaces'
    def __init__(self):
        self.graded_pieces={}
    
class gr_hom():
    'Graded linear maps'
    def __init__(self,d,V,W,M):
        '''
        Define a graded linear map V-->W of degree d.
        The matrix of the linear map V^n-->W^{n+d} is accessed as
        self.matrix[n].
        
        Input:
        d is int
        V, W are gr_vs
        M is dictionary of (np.matrix or hom)
        '''
        self.source,self.target=V,W
        self.degree=d
        self.matrix={}
        for n in M:
            if type(M[n]) is np.matrix:
                self.matrix[n]=vs.homomorphism(V[n],W[n+d],M[n])
            elif type(M[n]) is vs.homomorphism:
                self.matrix[n]=M[n]
            else:
                print("gr_hom requires M to be of type np.matrix or hom.")
    def __mul__(self,other):
        'Compose graded linear maps'
        new_gr_hom=gr_hom(self.degree+other.degree,self.source,other.target,{})
        d=self.degree
        for n in self.matrix:
            if n+d in other.matrix:
                # We only bother to keep track of maps between nonzero vector spaces
                new_gr_hom[n]=self.matrix[n]*other.matrix[n+d]
        return vs.homomorphism(self.source,other.target,(self.matrix[n])*(other.matrix[n]))
    def __add__(self,other):
        'Add linear maps'
        src1,src2=self.source,other.source
        tar1,tar2=self.target,other.target
        new_gr_hom=gr_hom(src_1,tar1,{})
        if src1==src2 and tar1==tar2:
            new_gr_hom.matrix=self.matrix
            for n in other.matrix:
                if n in new_gr_hom.matrix:
                    new_gr_hom.matrix[n]+=other.matrix[n]
                else:
                    new_gr_hom.matrix[n]=other.matrix[n]
            return new_gr_hom
        else:
            print("Cannot add graded linear maps with different sources/targets.")
    def oplus(self,other):
        '''
        Form a block matrix of graded linear maps
        '''
        d=self.degree
        if d!=other.degree:
            print("Can't direct sum maps of different degrees!")
        else:
            new_source=self.source+other.source
            new_target=self.target+other.target
            new_gr_hom=gr_hom(d,new_source,new_target,self.matrix)
            for n in other.matrix:
                if n in new_gr_hom:
                    new_gr_hom[n]+=other.matrix[n]
                else:
                    new_gr_hom[n]=other.matrix[n]
            return new_gr_hom
    def otimes(self,other):
        'Tensor/Kronecker/outer product f\otimes g of graded linear maps'
        ### Still needs to be done
        d=self.degree
        if d==other.degree:
            A=self.matrix
            B=other.matrix
            new_source=self.source*other.source
            new_target=self.target*other.target
            return vs.homomorphism(new_source,new_target,np.kron(A,B))# This is just wrong
        else:
            print("Tried to compute the tensor product of maps of different degree.")
    def ker(self):
        'Kernel of graded linear maps'
        new_gr_vs=gr_vs({})
        for n in self.matrix:
            new_gr_vs[n]=self.matrix[n].ker()
    def image(self):
        'Image of graded linear maps'
        new_gr_vs=gr_vs({})
        d=self.degree
        for n in self.matrix:
            new_gr_vs[n+d]=self.matrix[n].image()
    def cohomology(self):
        'Cohomology of a graded linear map, that is kernel/image'
        return (self.ker())/(self.image())

class zero_gr_hom(gr_hom):
    '''
    A class of zero graded linear maps
    (where specifying the source/target/degree is still meaningful).
    '''
    def __init__(self,d,V,W):
        self.source,self.target=V,W
        self.degree=d
        self.matrix={}
