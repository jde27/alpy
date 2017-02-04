#!/usr/bin/python

import rings as ri
import numpy as np
import graded_linear_algebra as gla
import copy

### There is some noticeable code reuse in some of the functions,
### for forming certain tensor products (denoted TV below).

class A8Category:
    '''The class of A_\infty-categories.

    Usage: A=A8Category(base,objects,morphisms,operations)

    Attributes:

        A.base [ri.Field]: Base field over which A is defined.

        A.objects [list]: Objects in A.

        A.morphisms [dict]:
            { (X,Y) [tuple of A.objects] : hom(X,Y) [gla.GradedVectorSpace] }

        A.operations [dict]:
            { (X_0,...,X_d) [tuple of A.objects] :
                \mu^d_A:hom(X_{d-1},X_d) (x) ... (x) hom(X_0,X_1)
                    ----> hom(X_0,X_d) [gla.GradedLinearMap] }

    Methods:

        A.hom(X,Y):
            Args: X, Y in A.objects
            Returns: [gla.GradedVectorSpace] A.morphisms[(X,Y)]
                if defined (and zero otherwise -
                we only store nonzero morphism spaces).

        A.mu((X_0,...,X_d)):
            Args: (X_0,...,X_d) tuple of A.objects
            Returns: [gla.GradedLinearMap] A.operations[(X_0,...,X_d)]
                if defined (and zero otherwise -
                we only store nonzero A_\infty operations).

        A.yoneda(X):
            Args: X in A.objects
            Returns: [A8Module] Yoneda module over A associated to X.
    '''
    def __init__(self,base,objects,morphisms,operations):
        self.base=base
        self.objects=objects
        self.morphisms=morphisms
        self.operations=operations

    def hom(self,X,Y):
        '''Provides safe access to morphism spaces.'''
        if (X,Y) in self.morphisms:
            # The morphism spaces are instances of GradedVectorSpace,
            # and hence complicated nested objects (including
            # a dictionary). We therefore use deepcopy to avoid
            # mutability backwash.
            return copy.deepcopy(self.morphisms[(X,Y)])
        else:
            # If we haven't bothered to define morphisms[(X,Y)]
            # then it's the zero vector space, which can be created
            # like this:
            return gla.GradedVectorSpace(self.base)

    def mu(self,word):
        '''Provides safe access to A_\infty-operations.

        word is a tuple of objects in the category.'''
        if word in self.operations:
            # The A_\infty-operations are instances of GradedLinearMap,
            # and hence complicated nested objects (including
            # a dictionary). We therefore use deepcopy to avoid
            # mutability backwash.
            return copy.deepcopy(self.operations[word])
        else:
            # If we haven't bothered to define mu(word)
            # then it's the zero map, which we construct as follows:
            d=len(word)-1
            # Forming the tensor product
            # TV = hom(X_{d-1},X_d) (x) ... (x) hom(X_0,X_1)
            TV=self.hom(word[d-1],word[d])
            for i in range(1,d):
                TV=TV.otimes(self.hom(word[d-i-1],word[d-i]))
            # Return the zero graded linear map
            # hom(X_{d-1},X_d) (x) ... (x) hom(X_0,X_1) --> hom(X_0,X_d)
            return gla.GradedLinearMap(2-d,TV,self.hom(word[0],word[d]))
        
    def yoneda(self,Q):
        '''Returns the Yoneda module over A associated to Q in A.objects.'''
        A=self
        if Q in A.objects:
            '''If M=Yoneda(Q) then M(X)=hom(X,Q)'''
            modules={X: A.hom(X,Q)
                     for X in A.objects if (X,Q) in A.morphisms}
            '''The operations   M(X) (x) hom(X',X) (x) ...
            are precisely the operations \mu_A on
            hom(X,Q) (x) hom(X',X) (x) ...
            However, we now index them by (...,X',X)
            instead of (...,X',X,Q)'''
            operations={word[:-1]: A.mu(word)
                        for word in A.operations
                        if word[-1]==Q}
            return A8Module(A,modules,operations)
        else:
            print("Yoneda module does not exist: ",Q," is not in ",A,".")
            return None

class A8Module:
    '''The class of A_\infty-modules over an A_\infty-category A.

    Usage: M=A8Module(cat,modules,operations)

    Attributes:

        M.cat [A8Category]: The A_\infty-category A.

        M.modules [dict]:
            { X in [A.objects] : M(X) [gla.GradedVectorSpace] }
    
        M.operations [dict]:
            { (X_0,...,X_{d-1}) in [A.objects] :
                \mu^d_M: M(X_{d-1}) (x) hom(X_{d-2},X_{d-1})
                        (x)...(x) hom(X_0,X_1)----> M(X_0)
                            [gla.GradedLinearMap] }
    
    Methods:

        M.mod(X):
            Args: X in A.objects
            Returns: [gla.GradedVectorSpace] M.module[X]
                if defined (and zero otherwise -
                we only store nonzero vector spaces).

        M.mu((X_0,...,X_{d-1})):
            Args: (X_0,...,X_{d-1}) tuple of A.objects
            Returns: [gla.GradedLinearMap]
                M.operations[(X_0,...,X_{d-1})] if defined
                (and zero otherwise - we only store nonzero
                A_\infty-module operations).
    
        M.twist(X):
            Args: X in A.objects
            Returns: [A8Module] Twist of M around X.
    '''
    def __init__(self,cat,modules,operations):
        self.cat=cat
        self.modules=modules
        self.operations=operations

    def mod(self,X):
        # The module spaces are instances of GradedVectorSpace,
        # and hence complicated nested objects (including
        # a dictionary). We therefore use deepcopy to avoid
        # mutability backwash.
        if X in self.modules:
            return copy.deepcopy(self.modules[X])
        else:
            # If we haven't bothered to define module[X]
            # then it's the zero vector space, which can be created
            # like this:
            return gla.GradedVectorSpace(self.cat.base)

    def mu(self,word):
        # The A_\infty-module operations are instances of
        # GradedLinearMap, and hence complicated nested objects
        # (including a dictionary). We therefore use deepcopy to avoid
        # mutability backwash.
        if word in self.operations:
            return copy.deepcopy(self.operations[word])
        else:
            # If we haven't bothered to define mu(word)
            # then it's the zero map, which we construct as follows:
            d=len(word)
            # Forming the tensor product
            # TV = M(X_{d-1}) (x) hom(X_{d-2},X_{d-1})
            #                 (x) ... (x) hom(X_0,X_1)
            TV=self.mod(word[d-1])
            for i in range(1,d):
                TV=TV.otimes(self.cat.hom(word[d-i-1],word[d-i]))
            # Return the zero graded linear map
            # M(X_{d-1}) (x) hom(X_{d-2},X_{d-1}) (x) ...
            #               ...(x) hom(X_0,X_1) ----> M(X_0)[2-d]
            return gla.GradedLinearMap(2-d,TV,self.mod(word[0]))

    def ltimes(self,cochain_complex):
        '''M.ltimes(Z) returns the tensor product Z(x)M
        of an A_\infty-module M with a chain complex (Z,d).
        '''
        A=self.cat
        Z=cochain_complex.cochains
        d=cochain_complex.differential
        modules={X: Z.otimes(self.mod(X)) for X in self.modules}
        # We first form the differential \mu^1_{Z(x)M}:
        # We have \mu^1_{Z(x)M}(z (x) b) =
        #          = (-1)^{|b|-1} dz (x) b + z (x) db
        # The signs are obtained by Koszulifying the identity on M(X).
        operations={(X,): d.otimes(self.mod(X).eye().koszulify())
                    +Z.eye().otimes(self.mu((X,)))
                    for X in self.modules}
        # Then we add in \mu^d, d\geq 2:
        operations.update({word: Z.eye().otimes(self.mu(word))
                           for word in self.operations
                           if len(word)!=1})
        return A8Module(A,modules,operations)

    def shift(self,m=1):
        '''Returns the A_\infty module shifted in degree by m.

        This is achieved by tensoring with K[m].
        '''
        K=self.cat.base
        cochains=gla.GradedVectorSpace(K)
        cochains.graded_dim={-m:1}
        differential=gla.GradedLinearMap(1,cochains,cochains)
        Z=gla.CochainComplex(cochains,differential)
        return self.ltimes(Z)
        
    
    def cpx(self,Q):
        '''Returns the cochain complex M(Q), \mu^1.'''
        cochains=self.mod(Q)
        differential=self.mu((Q,))
        return gla.CochainComplex(cochains,differential)
        
    def twist(self,Q):
        A=self.cat
        # We first form the tensor product M(Q) (x) yoneda(Q)
        Y=A.yoneda(Q)
        cochain_complex=self.cpx(Q)
        T=Y.ltimes(cochain_complex)
        # We then form the canonical evaluation morphism
        # ev^d(c(x)b,a_{d-1},...,a_1) = \mu^{d+1}_M(c,b,a_{d-1},...,a_1)
        components={word[:-1]: self.mu(word)
                    for word in self.operations
                    if word[-1]==Q and len(word)!=1}
        ev=A8ModuleMap(0,T,self,components)
        # Finally, we return the cone on ev
        return ev.cone()


    def __str__(self):
        return "%s" % str(self.modules)

    def verify(self):
        '''Check d^2=0'''
        A=self.cat
        M=self
        verify_dict={}
        for X in A.objects:
            zero_map1=gla.GradedLinearMap(2,M.mod(X),M.mod(X))
            verify_dict[X]=(M.mu((X,))*M.mu((X,))==zero_map1)
        '''Check Leibniz rule'''
        for X in A.objects:
            for Y in A.objects:
                zero_map2=gla.GradedLinearMap(1,M.mod(Y).otimes(A.hom(X,Y)),M.mod(X))
                verify_dict[(X,Y)]=(M.mu((X,))*M.mu((X,Y))
                                    +M.mu((X,Y))*(M.mu((Y,)).otimes(A.hom(X,Y).eye()))
                                    +M.mu((X,Y))*(M.mod(Y).eye().otimes(A.mu((X,Y))))
                                    ==zero_map2)
        '''Check module-associativity'''
        for X in A.objects:
            for Y in A.objects:
                for Z in A.objects:
                    zero_map3=gla.GradedLinearMap(0,(M.mod(Z).otimes(A.hom(Y,Z)).otimes(A.hom(X,Y))),M.mod(X))
                    verify_dict[(X,Y,Z)]=(M.mu((X,Y))*(M.mu((Y,Z)).otimes(A.hom(X,Y).eye()))
                                          +M.mu((X,Z))*(M.mod(Z).eye().otimes(A.mu((X,Y,Z))))
                                          ==zero_map3)
        if not all(verify_dict):
            print("Not an A_\infty module, sorry.")
        else:
            print("Your luck is in.")
    
class A8ModuleMap:
    '''The class of pre-module homomorphisms t: M-->N between
    A_\infty-modules.

    Usage: t=A8ModuleMap(degree,source,target)

    Attributes:
    
    t.degree [int]: The degree of the morphism t.

    t.source, t.target [A8Module]: The source and target of t.

    t.components [dict]:
        { (X_0,...,X_{d-1}) in A.objects :
                t^d: M(X_{d-1}) (x) hom(X_{d-2},X_{d-1})
                        (x)...(x) hom(X_0,X_1)----> N(X_0)
                            [gla.GradedLinearMap] }

    Methods:

        t.cpt((X_0,...,X_{d-1})):
            Args: X_0,...,X_{d-1} tuple of A.objects
            Returns: [gla.GradedLinearMap] t.components[(X_0,...,X_{d-1})]
                if defined (and zero otherwise -
                we only store nonzero parts of the morphism).

        t.cone() [A8Module]: Returns the cone on the morphism t.
    '''
    def __init__(self,degree,source,target,components):
        self.degree=degree
        self.source=source
        self.target=target
        self.components=components

    def cpt(self,word):
        if word in self.components:
            # The components of a morphism are instances of
            # GradedLinearMap, and hence complicated nested objects
            # (including a dictionary). We therefore use deepcopy to avoid
            # mutability backwash.
            return copy.deepcopy(self.components[word])
        else:
            # If we haven't bothered to define cpt[word]
            # then it's the zero map, which can be created
            # like this:
            d=len(word)
            # Forming the tensor product
            # TV = M(X_{d-1}) (x) hom(X_{d-2},X_{d-1})
            #                 (x) ... (x) hom(X_0,X_1)
            TV=self.source.mod(word[d-1])
            for i in range(1,d):
                TV=TV.otimes(self.source.cat.hom(word[d-i-1],word[d-i]))
            # Return the zero graded linear map
            # M(X_{d-1}) (x) hom(X_{d-2},X_{d-1}) (x) ...
            #               ...(x) hom(X_0,X_1) ----> N(X_0)[1+|t|-d]
            return gla.GradedLinearMap(1+self.degree-d,TV,
                                       self.target.mod(word[0]))

    def cone(self):
        '''Returns the cone on an A_\infty pre-module morphism.

        Would be nice if this checked that the morphism is closed.
        '''
        (M,N)=(self.source,self.target)
        A=M.cat
        # The module for the cone is C(X) = M(X)[1] (+) N(X)
        new_module={X: M.mod(X).shift().oplus(N.mod(X))
                    for X in (M.modules.keys()|N.modules.keys())}
        # The operations are the block matrices
        # (\mu_M  0    )
        # (F      \mu_N)
        # but F and \mu_M need to be suitably shifted
        # (using rejig_*) to make sense.
        L=M.shift()
        zero=A8ModuleMap(1,N,L,{})
        all_keys=(self.components.keys()
                  |M.operations.keys()
                  |N.operations.keys())
        new_operations={word:
                        gla.GradedLinearMap.block(
                            M.mu(word).rejig_2(),zero.cpt(word),
                            self.cpt(word).rejig_3(),N.mu(word))
                        for word in all_keys}
        return A8Module(A,new_module,new_operations)

class DynkinGraph():
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
        self.vertices,self.arrows=V,A

    def categorify(self,K,N):
        'Produce the associated dg-category of dimension N'

        def sph(K,N):
            '''Returns the graded vector space H^*(S^N;K).'''
            V=gla.GradedVectorSpace(K)
            V.graded_dim={0:1,N:1}
            return V
        
        def sph_op(V):
            '''Returns the ring operations on V=H^*(S^N;K).'''
            F=gla.GradedLinearMap(0,V.otimes(V),V)
            K=V.base
            N=max(V.graded_dim.keys())
            i=K.num(np.array([[1]]))
            if N%2==1:
                sigma=-1
            else:
                sigma=1
            j=K.num(np.array([[1,sigma]]))
            F.graded_map={0:i,2:j}
            return F
            
        def pt(K,N):
            '''Returns the graded vector space K in degree N.'''
            V=gla.GradedVectorSpace(K)
            V.graded_dim={N:1}
            return V
            
        def oper(A,B,C):
            # Sometimes in A(x)B --> C there is only one degree where
            # both domain and target have a nonzero summand
            # and the map is just the identity.
            deg=max(((A.otimes(B)).graded_dim.keys())&(C.graded_dim.keys()))
            F=gla.GradedLinearMap(0,A.otimes(B),C)
            i=K.num(np.array([[1]]))
            F.graded_map={deg:i}
            return F

        objects=self.vertices
        arrow=self.arrows
        morphisms={(X,X): sph(K,N) for X in objects}
        morphisms.update({(X,Y): pt(K,0)
                          for X in objects
                          for Y in arrow[X]})
        morphisms.update({(Y,X): pt(K,N)
                          for X in objects
                          for Y in arrow[X]})
        operations={(X,X,X): sph_op(morphisms[(X,X)])
                    for X in objects}

        '''We add in operations
        hom(a,b) (x) hom(a,a) --> hom(a,b)
        hom(b,b) (x) hom(a,b) --> hom(a,b)
        hom(b,a) (x) hom(a,b) --> hom(a,a)
        when a-->b in the directed Dynkin graph.
        '''
        operations.update({(a,b,c):
                           oper(morphisms[(b,c)],morphisms[(a,b)],
                                morphisms[(a,c)])
                           for a in objects
                           for b in objects
                           for c in objects
                           if ((a==b and c in arrow[a]) or
                               (a==c and b in arrow[a]) or
                               (b==c and a in arrow[b]))})
        '''We add in operations
        hom(a,b) (x) hom(b,a) --> hom(b,b)
        hom(a,a) (x) hom(b,a) --> hom(b,a)
        hom(b,a) (x) hom(b,b) --> hom(b,a)
        when a-->b in the directed Dynkin graph.
        '''
        operations.update({(a,b,c):
                           oper(morphisms[(b,c)],morphisms[(a,b)],
                                morphisms[(a,c)])
                           for a in objects
                           for b in objects
                           for c in objects
                           if ((a==b and a in arrow[c]) or
                               (a==c and a in arrow[b]) or
                               (b==c and b in arrow[a]))})
        '''We add in operations
        hom(Y,Z) (x) hom(X,Y) --> hom(X,Z)
        for every triangle X-->Y-->Z, X-->Z in the
        directed Dynkin graph.
        '''
        operations.update({(X,Y,Z):
                           oper(morphisms[(Y,Z)],morphisms[(X,Y)],
                                morphisms[(X,Z)])
                           for X in objects
                           for Y in arrow[X]
                           for Z in arrow[Y]
                           if Z in arrow[X]})
        return A8Category(K,objects,morphisms,operations)
        
    @staticmethod
    def BP(p,q,K,N):
        '''Generates the Dynkin diagram for a Brieskorn-Pham singularity
        of the form x^p+y^q=0.'''
        M=(p-1)*(q-1)
        vertices={m for m in range(1,M+1)}
        def nbhd(m):
            return {m+1,m+p,m+p+1} & {x for x in range(1,M+1)}
        arrows={m : nbhd(m) for m in range(1,M+1)}
        G=DynkinGraph(vertices,arrows)
        return G.categorify(K,N)
