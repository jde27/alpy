#!/usr/bin/python

import numpy as np
import rings as ri
import rationals as QQ
import finite_fields as fi
import graded_linear_algebra as gla
import a8_categories as acat

############################# REWRITE FROM HERE

class dynkin_graph():
    def __init__(self,V,A):
        '''
        Input a directed graph in the form:
        V - a set
        A - a dictionary of sets
        
        Example:
        V={'A','B','C'}
        A={'A':{'B','C'},'B':{'C'}}
        encodes the graph with three vertices A,B,C and arrows
        from A to B and C and from B to C.
        '''
        self.vertices,self.arrows=G,A
        
    def categorify(self,N):
        'Produce the associated dg-category of dimension N'
        obj=self.vertices
        mor={}
        diff={}
        prod={}
        sph=gvs(self.base)
        sph.setgps({0:1,N:1})  # The cohomology of S^N
        ###################THIS NEEDS MASSIVE EDITING:
        ########## NO LONGER USE MATRICES
        ############# AND have to use setgps
        m_sph={0:np.matrix('1'),N:np.matrix('1 (-1)**N')}
        #The sign comes from Seidel's convention [Seidel, p.8 Eq (1.3)]
        tween1=gvs(self.base)
        tween1.setgps({1:1})   # Rank 1, degree 1
        tween2=gvs(self.base)
        tween2.setgps({N-1:1}) # Rank 1, degree N-1
        i=np.matrix('1')
        for X in obj:
            mor[(X,X)]=sph
            prod[(X,X,X)]=ghomo(0,mor[(X,X)]*mor[(X,X)],mor[(X,X)])
            prod[(X,X,X)].setgps(m_sph)
        for X in obj:
            for Y in self.arrows[X]:
                mor[(X,Y)]=tween1
                mor[(Y,X)]=tween2
                # Need to get signs right! Have not checked them.
                prod[(X,X,Y)]=ghomo(0,mor[(X,Y)]*mor[(X,X)],mor[(X,Y)],{1:i})
                prod[(X,Y,X)]=ghomo(0,mor[(Y,X)]*mor[(X,Y)],mor[(X,X)],{N:i})
                prod[(X,X,Y)]=ghomo(0,mor[(X,Y)]*mor[(X,X)],mor[(X,Y)],{N-1:i})
                prod[(X,Y,Y)]=ghomo(0,mor[(Y,Y)]*mor[(X,Y)],mor[(X,Y)],{1:i})
                prod[(Y,X,Y)]=ghomo(0,mor[(X,Y)]*mor[(Y,X)],mor[(Y,Y)],{N:i})
                prod[(Y,Y,X)]=ghomo(0,mor[(Y,X)]*mor[(Y,Y)],mor[(Y,X)],{N-1:i})
        for X in obj:
            for Y in self.arrows[X]:
                for Z in (self.arrows[X].keys() & self.arrows[Y].keys()):
                    # Confused about the m_2 products from triangles - have asked Ailsa
                    prod[(X,Y,Z)]=ghomo(0,mor[(Y,Z)]*mor[(X,Y)],mor[(X,Z)],{})
        return dg_cat(obj,mor,diff,prod)
