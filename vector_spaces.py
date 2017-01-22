#!/usr/bin/python

import numpy as np

def rank_nullity(matrix):
    '''Performs Gaussian elimination to compute the rank and nullity of a np.matrix A'''
    A=matrix
    rows,cols=A.shape
    # We start with a guess of the rank and nullity,
    # setting them to be the biggest they could possibly be.
    rank,nullity=rows,cols
    for n in range(0,rows):
        # We find the coordinates of the nonzero entries in row n
        S=np.nonzero(A[n])
        # This is stored as an array: S[0][x],S[1][x] gives the
        # row,column (respectively) of the xth nonzero entry.
        # As we are working one row at a time, S[0][x] will always be 0.
        # If there are no nonzero entries then S[0] will be empty.
        if S[0].size==0:
            # If we find an empty row, we know the rank is
            # at least 1 less than our previous guess.
            rank-=1
        else:
            # Otherwise, we know the nullity is at least 1 less
            # than our previous guess.
            nullity-=1
            # We now clear column number S[1][0] by subtracting a suitable
            # multiple of row n from all subsequent rows.
            c=S[1][0]
            for x in range(n+1,rows):
                    t=A[x,c]/A[n,c]
                    for y in range(0,cols):
                        A[x,y]-=t*A[n,y]
    return rank,nullity

class vector_space():
    '''The class of vector spaces

    V=vector_space(n) creates a vector space of dimension n [int].

    Attributes:
        V.dim [int]: dimension of V

    Methods:
        V*W, V.otimes(W):
            Args: V, W [vector_space]
            Returns: [vector_space] Tensor product of V and W
        V+W, V.oplus(W):
            Args: V, W [vector_space]
            Returns: [vector_space] Direct sum of V and W
        V/W:
            Args: V, W [vector_space]
            Returns: [vector_space] Quotient of V by W
    '''
    def __init__(self,d):
        self.dim=d
    def otimes(self,other):
        '''Tensor product of two vector spaces'''
        return vector_space((self.dim)*(other.dim))
    def __mul__(self,other):
        '''Tensor product of two vector spaces'''
        return V.otimes(W)
    def oplus(self,other):
        '''Direct sum of two vector spaces'''
        new_vs=vector_space(self.dim+other.dim)
        return new_vs
    def __add__(self,other):
        '''Direct sum of two vector spaces'''
        return oplus(self,other)
    def __truediv__(self,other):
        '''Quotient vector space

        Currently this is just naively returning
        a vector space of dimension equal to the difference provided this
        is positive.
        '''
        if self.dim >= other.dim:
            return vector_space(self.dim-other.dim)
        else:
            print("Quotient vector space has negative dimension...")

class homomorphism():
    '''The class of linear maps between vector spaces

    F=homomorphism(V,W,M) creates a linear map
    F:V-->W with matrix M [np.matrix].
    Here V, W are of type [vector_space].

    Attributes:
        F.source [vector_space]: domain of F
        F.target [vector_space]: range of F
        F.matrix [np.matrix]: matrix of linear map F

    Methods:
        F*G:
            Args: F, G [homomorphism]
            Returns: [homomorphism] Composition of F with G
        F+G:
            Args: F, G [homomorphism] with same .source and .target
            Returns: [homomorphism] Sum of F with G in hom(source,target)
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
    def __init__(self,V,W,M):
        '''Args: V, W [vector_space]; M [np.matrix]'''
        self.source,self.target=V,W
        self.matrix=M
    def __mul__(self,other):
        '''Compose linear maps. Arg: other [homomorphism].'''
        return homomorphism(self.source,other.target,(self.matrix)*(other.matrix))
    def __add__(self,other):
        '''Add linear maps. Arg: other [homomorphism].'''
        src1,src2=self.source,other.source
        tar1,tar2=self.target,other.target
        if src1==src2 and tar1==tar2:
            return homomorphism(src1,tar1,self.matrix+other.matrix)
        else:
            print("Cannot add linear maps with different sources/targets.")
    def oplus(self,other):
        '''Form a block matrix. Arg: other [homomorphism].'''
        A=self.matrix
        B=other.matrix
        new_source=self.source+other.source
        new_target=self.target+other.target
        return homomorphism(new_source,new_target,np.bmat[[A,0],[0,B]])
    def otimes(self,other):
        '''Tensor product of linear maps f\otimes g
        (also known as Kronecker or outer product).
        Arg: other [homomorphism].
        '''
        A=self.matrix
        B=other.matrix
        new_source=self.source*other.source
        new_target=self.target*other.target
        return homomorphism(new_source,new_target,np.kron(A,B))
    def ker(self):
        '''Find the kernel of linear map.'''
        rank,nullity=rank_nullity(self.matrix)
        return vector_space(nullity)
    def image(self):
        '''Find the image of linear map'''
        rank,nullity=rank_nullity(self.matrix)
        return vector_space(rank)
