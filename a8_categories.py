#!/usr/bin/python

import rings as ri
import graded_linear_algebra as gla

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

        A.mu(X_0,...,X_d):
            Args: X_0,...,X_d in A.objects
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
        if (X,Y) in self.morphisms:
            return self.morphisms[(X,Y)]
        else:
            # If we haven't bothered to define morphisms[(X,Y)]
            # then it's the zero vector space, which can be created
            # like this:
            return gla.GradedVectorSpace(self.base)

    def mu(self,*args):
        if args in self.operations:
            return self.operations[args]
        else:
            d=len(args)-1
            # Forming the tensor product
            # TV = hom(X_{d-1},X_d) (x) ... (x) hom(X_0,X_1)
            TV=copy.deepcopy(self.hom(args[d-1],args[d]))###
            for i in range(1,d):
                TV=TV.otimes(self.hom(args[d-i-1],arg[d-i]))
            # Return the zero graded linear map
            # hom(X_{d-1},X_d) (x) ... (x) hom(X_0,X_1) --> hom(X_0,X_d)
            return gla.GradedLinearMap(2-d,TV,self.hom(args[0],args[d]))
        
    def yoneda(self,Q):
        '''Returns the Yoneda module over A associated to Q in A.objects.'''
        A=self
        if Q in A.objects:
            modules={X: A.hom(X,Q)
                     for X in A.objects if (X,Q) in A.morphisms}
            operations={word: A.operations(word+(Q,))
                        for word in A.operations
                        if word+(Q,) in A.operations}
            # A.operations is indexed by tuples of objects
            # so we need to turn Q into a 1-tuple.
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

        M.mu(X_0,...,X_{d-1}):
            Args: X_0,...,X_{d-1} in A.objects
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
        if X in self.modules:
            return self.modules[X]
        else:
            # If we haven't bothered to define module[X]
            # then it's the zero vector space, which can be created
            # like this:
            return gla.GradedVectorSpace(self.cat.base)

    def mu(self,*args):
        if args in self.operations:
            return self.operations[args]
        else:
            d=len(args)
            # Forming the tensor product
            # TV = M(X_{d-1}) (x) hom(X_{d-2},X_{d-1})
            #                 (x) ... (x) hom(X_0,X_1)
            TV=self.module(args[d-1])
            for i in range(1,d):
                TV=TV.otimes(self.cat.hom(args[d-i-1],arg[d-i]))
            # Return the zero graded linear map
            # M(X_{d-1}) (x) hom(X_{d-2},X_{d-1}) (x) ...
            #               ...(x) hom(X_0,X_1) ----> M(X_0)[2-d]
            return gla.GradedLinearMap(2-d,TV,self.module(args[0]))

    def .otimes(self,cochain_complex):
        '''Form the tensor product of an A_\infty-module with
        a chain complex (Z,d).
        '''
        A=self.cat
        Z=cochain_complex.cochains
        d=cochain_complex.differential
        modules={X: Z.otimes(self.mod(X))) for X in self.modules}
        # We have \mu^1_{Z(x)M}(z (x) b) =
        #          = (-1)^{|b|-1} dz (x) b + z (x) db
        # The signs are obtained by Koszulifying the identity on M(X).
        operations={X: d.otimes(self.mod(X).Id().koszulify())\
                    for X in self.modules}
        # We have \mu^d_{Z(x)M}(z(x)b,a...) = z(x)\mu^d_M(b,a...)
        # for d\geq 2:
        operations.update({word: Z.Id().otimes(self.mu(word))\
                     for word in (self.operations and not A.objects)})
        
    def cpx(self,Q):
        '''Returns the cochain complex M(Q), \mu^1.'''
        cochains=self.mod(Q)
        differential=self.mu(Q)
        return CochainComplex(cochains,differential)
        
    def twist(self,Q):
        A=self.cat
        # We first form the tensor product M(Q) (x) yoneda(Q)
        Y=A.yoneda(Q)
        cochain_complex=self.cpx(Q)
        T=Y.otimes(cochain_complex)
        # We then form the canonical evaluation morphism
        # ev^d(c(x)b,a_{d-1},...,a_1) = \mu^{d+1}_M(c,b,a_{d-1},...,a_1)
        components={word: self.m(word+Q)\
                    for word in self.operations}
        ev=A8ModuleMap(T,self,components)
        # Finally, we return the cone on ev
        return cone(ev)

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

        t.cpt(X_0,...,X_{d-1}):
            Args: X_0,...,X_{d-1} in A.objects
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

    def cpt(self,*args):
        if args in self.componens:
            return self.operations[args]
        else:
            d=len(args)
            # Forming the tensor product
            # TV = M(X_{d-1}) (x) hom(X_{d-2},X_{d-1})
            #                 (x) ... (x) hom(X_0,X_1)
            TV=self.source.module(args[d-1])
            for i in range(1,d):
                TV=TV.otimes(self.source.cat.hom(args[d-i-1],arg[d-i]))
            # Return the zero graded linear map
            # M(X_{d-1}) (x) hom(X_{d-2},X_{d-1}) (x) ...
            #               ...(x) hom(X_0,X_1) ----> N(X_0)[1+|t|-d]
            return gla.GradedLinearMap(1+self.degree-d,TV,
                                   self.target.module(args[0]))

    def cone(self):
        '''Returns the cone on an A_\infty pre-module morphism.

        Would be nice if this checked that the morphism is closed.
        '''
        M,N=self.source,self.target
        A=M.cat
        # The module for the cone is C(X) = M(X)[1] (+) N(X)
        new_module={X: M.mod(X).shift().oplus(N.mod(X))
                    for X in (M.modules or N.modules)}
        # The operations are the block matrices
        # (\mu_M  0    )
        # (F      \mu_N)
        # but F and \mu_M need to be suitably shifted
        # (using rejig_*) to make sense.
        zero=A8ModuleMap(1,N,M)
        new_operations={word:
                        gla.GradedLinearMap.block(
                            M.mu(word).rejig_2(),zero.cpt(word),
                            ev.cpt(word).rejig_3(),N.mu(word))}
        return A8Module(A,new_module,new_operations)
