#!/usr/bin/python

import array_manipulations
import numpy as np
import arithmetic as ar
import rings as ri
import rationals as QQ

def rank_nullity(matrix):
    '''Computes the rank and nullity of a 2-d np.array
    via Gaussian elimination over the field of definition'''
    A=matrix
    rows,cols=A.shape
    # We start with a guess of the rank and nullity,
    # setting them to be the biggest they could possibly be.
    rank,nullity=rows,cols
    for n in range(0,rows):
        # We look for the first nonzero element in row n
        found,s,t=array_manipulations.find_nonzero(A[n])
        if found:
            # If there is a nonzero element we know the
            # nullity is at least 1 less than our previous guess.
            nullity-=1
            # We now clear column number s by subtracting a suitable
            # multiple of row n from all subsequent rows.
            for x in range(n+1,rows):
                c=A[x,s]/t
                for y in range(0,cols):
                    A[x,y]-=c*A[n,y]
        else:
            # If we find an empty row, we know the rank is
            # at least 1 less than our previous guess.
            rank-=1
    return rank,nullity

class vs():
    '''The class of vector spaces

    V=vs(k,n) creates a vector space of dimension n [int] over k [ri.field].

    Attributes:
        V.dim [int]: dimension of V
        V.base [ri.field]: field of definition for V

    Methods:
        V*W, V.otimes(W):
            Args: V, W [vs]
            Returns: [vs] Tensor product of V and W
        V+W, V.oplus(W):
            Args: V, W [vs]
            Returns: [vs] Direct sum of V and W
        V/W:
            Args: V, W [vs]
            Returns: [vs] Quotient of V by W
    '''
    def __init__(self,k,d):
        self.dim=d
        self.base=k

    def otimes(self,other):
        '''Tensor product of two vector spaces'''
        if self.base==other.base:
            return vs(self.base,(self.dim)*(other.dim))
        else:
            print("Cannot tensor vector spaces over different fields!")
            
    def __mul__(self,other):
        '''Tensor product of two vector spaces'''
        return self.otimes(other)
    
    def oplus(self,other):
        '''Direct sum of two vector spaces'''
        if self.base==other.base:
            new_vs=vs(self.base,self.dim+other.dim)
            return new_vs
        else:
            print("Cannot add vector spaces over different fields!")
            
    def __add__(self,other):
        '''Direct sum of two vector spaces'''
        return self.oplus(other)
    
    def __truediv__(self,other):
        '''Quotient vector space

        Currently this is just naively returning
        a vector space of dimension equal to the difference provided this
        is positive.
        '''
        if self.dim >= other.dim:
            return vs(self.base,self.dim-other.dim)
        else:
            print("Quotient vector space has negative dimension...")

class homo():
    '''The class of linear maps between vector spaces

    F=homo(V,W,M) creates a linear map
    F:V-->W with matrix M [np.ndarray].
    Here V, W are of type [vs].

    Attributes:
        F.source [vs]: domain of F
        F.target [vs]: range of F
        F.matrix [np.ndarray]: matrix of linear map F
        F.base [ri.field]: field of definition of F

    Methods:
        F*G:
            Args: F, G [homo]
            Returns: [homo] Composition of F with G
        F+G:
            Args: F, G [homo] with same .source and .target
            Returns: [homo] Sum of F with G in hom(source,target)
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
    def __init__(self,V,W,M):
        '''Args: V, W [vs]; M [np.ndarray]'''
        self.source,self.target=V,W
        self.base=self.source.base
        self.matrix=self.base.mod(M)
        
    def __mul__(self,other):
        '''Compose linear maps. Arg: other [homo].'''
        return homo(self.source,other.target,(self.matrix).dot(other.matrix))
    
    def __add__(self,other):
        '''Add linear maps. Arg: other [homo].'''
        src1,src2=self.source,other.source
        tar1,tar2=self.target,other.target
        if src1==src2 and tar1==tar2:
            return homo(src1,tar1,self.matrix+other.matrix)
        else:
            print("Cannot add linear maps with different sources/targets.")
            
    def oplus(self,other):
        '''Form a block matrix. Arg: other [homo].'''
        A=self.matrix
        B=other.matrix
        V1=self.source
        V2=other.source
        W1=self.target
        W2=other.target
        k=V1.base
        new_source=V1+V2
        new_target=W1+W2
        y1,x1=A.shape
        y2,x2=B.shape
        C=k.mod(np.zeros(shape=(y1,x2)))
        D=k.mod(np.zeros(shape=(y2,x1)))
        E=np.asarray(np.bmat([[A,C],[D,B]])) # bmat outputs np.matrix
        return homo(new_source,new_target,E)
    
    def otimes(self,other):
        '''Tensor product of linear map f\otimes g
        (also known as Kronecker or outer product).
        Arg: other [homo].
        '''
        A=self.matrix
        B=other.matrix
        new_source=self.source*other.source
        new_target=self.target*other.target
        return homo(new_source,new_target,\
                    np.asarray(np.kron(A,B))) #np.kron outputs np.matrix
    
    def ker(self):
        '''Find the kernel of linear map.'''
        rank,nullity=rank_nullity(self.matrix)
        return vs(self.base,nullity)
    
    def image(self):
        '''Find the image of linear map'''
        rank,nullity=rank_nullity(self.matrix)
        return vs(self.base,rank)

'''Some simple examples
Q=QQ.QQ()
V=vs(3,Q)
W=vs(3,Q)
u=ri.number(Q,(1,2))
M=np.array([[0,0,0],[0,0,0],[0,0,1]])
N=np.array([[u,0,0],[0,0,0],[1,0,1]])
F=homo(V,W,M)
G=homo(V,W,N)
print((G.otimes(F)).ker().dim)
'''
