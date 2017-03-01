#!/usr/bin/python

import fields as fi
from linear_algebra import *

class A8Category(AlgebraicStructure):
    '''The class of A_\infty-categories.

    USAGE:

        A=A8Category(K,objects,morphisms,operations)
    
    creates an A_\infty-category over the field K, with given objects,
    morphisms and operations.

    ATTRIBUTES:

        A.field [fi.Field]

            The field over which A is defined.

        A.objects [set]

            A set of objects (usually integers, but could be
            anything).

        A.morphisms [dict] {(X,Y): VectorSpace}

            A dictionary indexed by 2-tuples from the set
            A.objects. The entry A.morphisms[(X,Y)] is the space of
            morphisms in A from X to Y. This can also be accessed
            using:

                A[(X,Y)]

            which will return the zero vector space
            if A.morphisms[(X,Y)] has not been set.

        A.operations [dict] {(X_0,...,X_d): LinearMap}

            A dictionary indexed by (d+1)-tuples from the set
            A.objects. The entry A.operations[(X_0,...,X_d)] encodes
            the A_\infty-operation

              \mu_A^d: A[(X_{d-1},X_d)] (x) ... (x) A[(X_0,X_1)]
                                ----> A[(X_0,X_d)]

            This can also be accessed using

                A.mu(X_0,...,X_d)

            which will return 0 if the operation has not been set.

    METHODS:

        A[(X,Y)]

            Returns A.morphisms[(X,Y)] if defined or else zero.

        A.hom(X_0,...,X_d)

            Returns the tensor product:

                A[(X_{d-1},X_d)] (x) ... (x) A[(X_0,X_1)]

            This is indexed by d-tuples of basis elements.

        A.mu(X_0,...,X_d)

            Returns A.operations[(X_0,...,X_d)] if defined or else
            zero.

        A.verify()

            Returns true if the A_\infty-equations are satisfied by
            A. Otherwise raises an exception.

        A.yoneda(X)

            Returns the A_\infty Yoneda module associated to the
            object X.

    '''

    _required_fields=['field','objects','morphisms','operations']

    def __getitem__(self,i):
        '''Returns A.morphisms[i] if defined and zero otherwise.'''
        if i in self.morphisms:
            return self.morphisms[i]
        else:
            return VectorSpace(self.field)
        
    def hom(self,*word):
        '''Returns the tensor product:
        A[(X_{d-1},X_d)] (x) ... (x) A[(X_0,X_1)]
        '''
        d=len(word)-1
        A,K=self,self.field
        tensor_list=()
        for i in range(0,d):
            tensor_list=(A[(word[i],word[i+1])],)+tensor_list
        return VectorSpace.tensor(*tensor_list)

    def mu(self,*word):
        A=self
        # Either we have defined the operation associated to 'word'
        # or else we haven't and we want to return the zero map.
        if word in A.operations:
            return A.operations[word]
        else:
            d=len(word)-1
            return LinearMap(A.hom(*word),A[(word[0],word[d])],2-d)

    def verify(self):
        '''Returns true if the A_\infty-equations hold for A and raises an
        exception otherwise.'''
        A=self
        super_words=set()
        for word in A.operations:
            if A.mu(*word).deg!=3-len(word):
                raise ValueError('Not an A_\infty-module: '+
                                 'degrees of operations are wrong')

        for inner_op in A.operations:
            for outer_op in A.operations:
                pair=(inner_op[0],inner_op[-1])
                for k, val in enumerate(outer_op):
                    if ((k<len(outer_op)-1) and ((val,outer_op[k+1])==pair)):
                        if k+2<len(outer_op):
                            super_words.update({outer_op[:k]+inner_op+outer_op[k+2:]})
                        else:
                            super_words.update({outer_op[:k]+inner_op})

        def check_word(*word):
            d=len(word)-1

            def term(n,m):
                tuple_of_maps=()
                for i in range(1,d-n-m+1):
                    tuple_of_maps=tuple_of_maps+(A[(word[d-i],word[d-i+1])].Id(),)
                tuple_of_maps=tuple_of_maps+(A.mu(*word[n:n+m+1]),)
                for i in range(1,n+1):
                    # The signs in the A_\infty equation are handled by using
                    # the function sigma here:
                    tuple_of_maps=tuple_of_maps+(A[(word[n-i],word[n-i+1])].sigma(),)
                F=LinearMap.tensor(*tuple_of_maps)
                cut_word=word[0:n+1]+word[n+m:d+1]
                if m>1 and m!=d:
                    return (A.mu(*cut_word).circ(F)).flatten(d-n-m)
                else:
                    return A.mu(*cut_word).circ(F)

            associativity_test=LinearMap(A.hom(*word),A[(word[0],word[d])],3-d)
            for N in range(0,d+1):
                for M in range(1,d-N+1):
                    associativity_test+=term(N,M)
            return (associativity_test==LinearMap(A.hom(*word),A[(word[0],word[d])],3-d))

        truth_dictionary={}
        for word in super_words:
            truth_dictionary[word]=check_word(*word)
        if not all(truth_dictionary):
            raise ValueError('A_\infty equations not satisfied')
        else:
            return True

    def yoneda(self,Q):
        '''Returns the Yoneda A_\infty-module over A associated to the object X.'''
        A=self
        if Q in A.objects:
            modules={X: A[(X,Q)] for X in A.objects if (X,Q) in A.morphisms}
            operations={word[:-1]: A.mu(*word)
                        for word in A.operations
                        if word[-1]==Q}
            return A8Module(A,modules,operations)
        else:
            raise ValueError('Cannot form Yoneda module: ',
                             Q,' is not an object of the category ',A)

    def total_yoneda(self):
        total=A8Module(self,{},{})
        for X in self.objects:
            total=(total.oplus(self.yoneda(X))).simplify()
        return total
        
class DynkinGraph():
    '''The class of Dynkin graphs.

    Input a directed graph in the form:

        V - a set
        A - a dictionary of dictionaries
            {source1: {target1:grading of arrow1,...,targetn:grading of arrow n},
             source2: {target1:grading of arrow1,...,targetn:grading of arrow n},
             ...}

    Example:

        V={'A','B','C'}
        A={'A':{'B':1,'C':1},'B':{'C':2}}
        G=DynkinGraph(V,A)

    encodes the graph with three vertices A,B,C and arrows from A to B
    and C (both degree 1) and one from B to C (degree 2).

    METHODS:

        G.categorify(K,N)

            Produces a Calabi-Yau-N A_\infty-category over the field K
            whose objects are in bijection with G.vertices, whose
            morphisms are determined by G.arrows in the following way:

                hom(X,X) = H^*(S^N;K)
                hom(X,Y) = K in degree G.arrows[X][Y]
                              if Y is in G.arrows[X]
                hom(Y,X) = K in degree N-G.arrows[X][Y]
                              if Y is in G.arrows[X]

            and whose operations are uniquely determined by the
            conditions that the category is CY-N and that a triangle
            X-->Y-->Z, X-->Z in the Dynkin graph gives a nontrivial
            composition

                hom(Y,Z) (x) hom(X,Y) --> hom(X,Z).
    '''

    def __init__(self,V,A):
        self.vertices,self.arrows=V,A
        
    def categorify(self,K,N):
        'Produce the associated A_\infty-category of dimension N'
        
        def sph(K,N):
            '''Returns the graded vector space H^*(S^N;K).'''
            V=VectorSpace(K)
            V.basis.update({0:0,1:N})
            F=LinearMap(V.otimes(V),V,0)
            F.maps.update({(0,0):V[0],
                           (1,0):-(V.sigma())(V[1]),
                           (0,1):V[1]})
            return V,F
        
        def pt(K,D):
            '''Returns the field K in degree D'''
            V=VectorSpace(K)
            V.basis.update({0:D})
            return V
        
        objects=self.vertices
        arrow=self.arrows
        morphisms={}
        operations={}

        A=A8Category(K,objects,morphisms,operations)

        # Initialise endomorphism spaces to be H^*(S^N;K)
        for X in objects:
            A.morphisms[(X,X)],A.operations[(X,X,X)]=sph(K,N)
        # Initialise morphism spaces between objects
        # connected by arrows
        A.morphisms.update({(X,Y): pt(K,arrow[X][Y])
                          for X in objects
                          for Y in arrow[X]})
        A.morphisms.update({(Y,X): pt(K,N-arrow[X][Y])
                          for X in objects
                          for Y in arrow[X]})
        # Add in operations to ensure the CY condition holds
        for X in objects:
            for Y in arrow[X]:
                F={}
                for triple in [(X,X,Y),(X,Y,X),(Y,X,X),(X,Y,Y),(Y,X,Y),(Y,Y,X)]:
                    F[triple]=LinearMap(A.hom(*triple),A[(triple[0],triple[2])],0)
                    V=F[triple].target
                    if triple[0]==triple[2]:
                        F[triple].maps.update({(0,0):V[1]})
                    else:
                        F[triple].maps.update({(0,0):V[0]})
                A.operations.update(F)
        # Add in operations corresponding to triangles in the graph
        for X in objects:
            for Y in arrow[X]:
                for Z in arrow[Y]:
                    if Z in arrow[X]:
                        F=LinearMap(A.hom(X,Y,Z),A[(X,Z)],0)
                        F.maps.update({(0,0):A[(X,Z)][0]})
                        A.operations.update({(X,Y,Z): F})
        
        return A
                        
def BP(p,q,K,N,d):
    '''Generates the Fukaya category of a Brieskorn-Pham Milnor fibre of
    the form x^p+y^q+z_1^2+...+z_{N-1}^2=1. The morphisms between
    spheres are in degree d or 2d.

    '''
    M=(p-1)*(q-1)
    vertices={m for m in range(1,M+1)}
    def nbhd(m):
        # Change your grading conventions here:
        #nbh={m+1:0,m+p:0,m+p+1:0}
        nbh={m+1:d,m+p:d,m+p+1:2*d}
        return {x: nbh[x] for x in nbh if x in range(1,M+1)}
        
        
    arrows={m : nbhd(m) for m in range(1,M+1)}
    G=DynkinGraph(vertices,arrows)
    return G.categorify(K,N)

    
class CochainComplex(AlgebraicStructure):
    '''The class of cochain complexes.

    USAGE:

        Z=CochainComplex(cochains,differential)
    
    ATTRIBUTES:

        cochains [VectorSpace]

        differential [LinearMap]

    METHODS:

        Z.cohomology()

            Returns the cohomology of Z as a dictionary of the form:

                {i: rank of H^i(Z)}

        Z.display()

            Displays the space of cochains and the
            differential. Mostly a debugging tool.

        Z.verify()

            Returns true if d^2=0, otherwise raises exception.

        Z.otimes(M)

            Returns the tensor product of a cochain complex with an
            A_\infty-module M.

    '''

    _required_fields=['cochains','differential']

    def cohomology(self):
        '''Returns the cohomology of the cochain complex as a dictionary of
        the form {i: rank of H^i(Z)}
        '''
        G=self.cochains.graded_pieces # Dictionary of Z_n, graded pieces of Z
        graded_maps={} # Dictionary to store restriction d_n of d to Z_n
        kernels={}     # Dictionary to store kernels of d_n
        images={}      # Dictionary to store images of d_n
        cohom={}       # Dictionary to store cohomology groups
        graded_maps={n: self.differential.restrict(G[n]) for n in G}
        for n in G:
            kernels[n],images[n]=graded_maps[n].ker_im()
        for n in G:
            if n-1 in G:
                cohom[n]=kernels[n]-images[n-1]
            else:
                cohom[n]=kernels[n]
        return cohom

    def display(self):
        '''Prints the space of cochains and the differential.'''
        print(self.cochains)
        self.differential.display()
    
    def verify(self):
        '''Returns true if d^2=0 and raises an exception otherwise.'''
        d=self.differential
        Z=self.cochains
        d.verify()
        zero_map=LinearMap(Z,Z,2)
        if d.circ(d) != zero_map:
            raise TypeError('Not a cochain complex')
        else:
            return True

    def otimes(self,M):
        '''Returns the tensor product of a chain complex with
        with an A_\infty-module M.
        '''
        Z,diff=self.cochains,self.differential
        A=M.cat
        modules={X: Z.otimes(M[X]) for X in M.modules}
        operations={(X,): Z.Id().otimes(M.mu(X))+diff.otimes(M[X].sigma())
                    for X in M.modules}
        for word in M.operations:
            if len(word)>1:
                new_op=Z.Id().otimes(M.mu(*word))
                operations.update({word: new_op.flatten(1).unflatten(0,2)})               
        return A8Module(A,modules,operations)
                
class A8Module(AlgebraicStructure):
    '''The class of A_\infty-modules.

    USAGE:

        M=A8Module(A,modules,operations)

    creates an A_\infty-module over the A_\infty-category A with
    specified module-spaces and operations.

    ATTRIBUTES:

        M.cat [A8Category]

            The category over which this module is defined.

        M.field [fi.Field]

            The field over which the module is defined (inherited from
            underlying category).

        M.modules [dict] {X: VectorSpace}

            A dictionary indexed by A.objects. The entry M.modules[X]
            is the module M at the object X.  This can also be
            accessed using:

                M[X]

            which will return the zero vector space
            if M.modules[(X,Y)] has not been set.

        M.operations [dict] {(X_0,...,X_{d-1}): LinearMap}

            A dictionary indexed by d-tuples from the set
            A.objects. The entry M.operations[(X_0,...,X_{d-1})]
            encodes the A_\infty module-operation

        \mu_M^d: M[X_{d-1}] (x) A[(X_{d-2},X_{d-1})] (x) ... (x) A[(X_0,X_1)]
                                    ----> M[X_0]

            This can also be accessed using

                M.mu(X_0,...,X_{d-1})

            which will return 0 if the operation has not been set.

    METHODS:

        M[X]

            Returns M.modules[X] if defined or else zero.

        M.mu(X_0,...,X_d)

            Returns M.operations[(X_0,...,X_d)] if defined or else
            zero.

        M.display()

        M.verify()

            Returns true if the A_\infty-module equations are
            satisfied by M. Otherwise raises an exception.

        M.simplify()

            Returns the same A_\infty-module but re-indexes the bases
            for its vector spaces to make it more compact.

        M.cpx(X)

            Returns the cochain complex M[X] with differential
            M.mu(X).

        M.total()

            Returns the direct sum of cohomology groups:

                (+)_{X in A.objects} M.cpx(X).cohomology()

        M.width()

            Returns the difference between the maximal and minimal
            grading present in M.total().

        M.shift(m=1)

            Returns M shifted down in degree by m.

        M.twist(X)

            Returns the twist of M around the object X.
    '''

    _required_fields=['cat','modules','operations']

    @lazyproperty
    def field(self):
        '''Returns the field of definition of the module.'''
        return self.cat.field
    
    def __getitem__(self,X):
        '''M[X] returns M.modules[X] if defined and zero otherwise.'''
        M=self
        K=M.field
        if X in M.modules:
            return M.modules[X]
        else:
            return VectorSpace(K)

    def mu(self,*word):
        '''M.mu(X_0,...,X_{d-1}) returns M.operations[(X_0,...,X_{d-1})] if
        defined and zero otherwise.'''
        M=self
        if word in M.operations:
            return M.operations[word]
        else:
            A=M.cat
            d=len(word)
            X,Y=word[d-1],word[0]
            if d>2:
                new_domain=(M[X].otimes(A.hom(*word))).flatten(1)
            elif d==2:
                new_domain=M[X].otimes(A.hom(*word))
            else:
                new_domain=M[X]
            return LinearMap(new_domain,M[Y],2-d)

    def display(self):
        '''Displays the modules and operations of M.'''
        for X in self.modules:
            print('M(',X,') = ',self[X])
        for word in self.operations:
            if len(word)==1:
                print('M(',word[-1],') = ')
            else:
                print('M(',word[-1],') * A.hom(',word,') = ')
            self.operations[word].display()

    def verify(self):
        '''Returns true if the A_\infty-module equations hold for M and raises
        an exception otherwise.'''
        M,A=self,self.cat
        super_words=set()
        for word in M.operations:
            print('Verifying module operations are well-defined: ',word)
            M.mu(*word).verify()
            if M.mu(*word).deg!=2-len(word):
                raise ValueError('Not an A_\infty-module: '+
                                 'degrees of operations are wrong')

        for inner_op in A.operations:
            for outer_op in M.operations:
                pair=(inner_op[0],inner_op[-1])
                for k, val in enumerate(outer_op):
                    if ((k<len(outer_op)-1) and ((val,outer_op[k+1])==pair)):
                        if k+2<len(outer_op):
                            super_words.update({outer_op[:k]+inner_op+outer_op[k+2:]})
                        else:
                            super_words.update({outer_op[:k]+inner_op})

        super_words=set()
        for inner_op in M.operations:
            for outer_op in M.operations:
                if outer_op[-1:]==inner_op[0]:
                    super_words.update({outer_op[:-1]+inner_op})
        
        def check_word(*word):
            d=len(word) # For consistency with Seidel's notation

            def term_1(n,m):
                tuple_of_maps=(M[word[d-1]].Id(),)
                for i in range(2,d-n-m+1):
                    tuple_of_maps=tuple_of_maps+(A[(word[d-i],word[d-i+1])].Id(),)
                tuple_of_maps=tuple_of_maps+(A.mu(*word[n:n+m+1]),)
                for i in range(1,n+1):
                    # The signs in the A_\infty equation are handled by using
                    # the function sigma here:
                    tuple_of_maps=tuple_of_maps+(A[(word[n-i],word[n-i+1])].sigma(),)
                F=LinearMap.tensor(*tuple_of_maps)
                cut_word=word[0:n+1]+word[n+m:]
                if m>1 and m!=d-1:
                    return (M.mu(*cut_word).circ(F)).flatten(d-n-m)
                else:
                    return M.mu(*cut_word).circ(F)
                
            def term_2(n):
                tuple_of_maps=(M.mu(*word[n:]),)
                for i in range(1,n+1):
                    tuple_of_maps=tuple_of_maps+(A[(word[d-i],word[d-i+1])].sigma(),)
                F=LinearMap.tensor(*tuple_of_maps)
                cut_word=word[:n+1]+(word[-1],)
                if n!=0 and n!=d-1:
                    return (M.mu(*cut_word).circ(F)).flatten(0)
                else:
                    return M.mu(*cut_word).circ(F)

            module_source,module_target=M.mu(*word).source,M.mu(*word).target
            associativity_test=LinearMap(module_source,module_target,3-d)
            for N in range(0,d):
                for M in range(1,d-N):
                    associativity_test+=term_1(N,M)
                associativity_test+=term_2(N)
            return (associativity_test==LinearMap(module_source,module_target),3-d)
        
        truth_dictionary={}
        for word in super_words:
            truth_dictionary[word]=check_word(*word)
        if not all(truth_dictionary):
            raise ValueError('A_\infty equations not satisfied')
        else:
            return True

    def simplify(self):
        '''After repeatedly twisting modules, the bases of the resulting
        modules are indexed by (often messy) nested tuples:

            M.simplify()

        will return the A_\infty-module M but with all its vector spaces
        indexed simply by integers.
        '''
        M=self
        N=A8Module(self.cat,{},{})
        translator={}
        for X in self.modules:
            N.modules[X]=VectorSpace(N.field)
            translator[X]=({i: j for j,i in
                            enumerate(list(M[X].basis.keys()))})
            N.modules[X].basis.update({translator[X][i]: M[X].basis[i]
                                       for i in translator[X]})
       
        def tr_bas(Y,*i):
            return (translator[Y][i[0]],)+i[1:]

        def tr_fun(F,*word):
            X=word[-1]
            d=len(word)
            Z=word[0]
            new_maps={}
            for i in F.maps:
                new_sp=N.mu(*word).target
                if d==1:
                    new_idx=translator[X][i]
                else:
                    new_idx=tr_bas(X,*i)
                new_cpts={translator[Z][j]: F.maps[i].components[j]
                          for j in F.maps[i].components}
                new_maps[new_idx]=Vector(new_sp,new_cpts)
            return new_maps
        
        for word in M.operations:
            d=len(word)
            N.operations[word]=LinearMap(N.mu(*word).source,N.mu(*word).target,2-d)
            N.operations[word].maps.update(tr_fun(M.operations[word],*word))

        return N

    def cpx(self,X):
        '''Returns the cochain complex M(X), \mu^1.'''
        cochains=self[X]
        differential=self.mu(X)
        return CochainComplex(cochains,differential)

    def total(self):
        '''Given an A_\infty-module M over a category A,
        M.total() returns the direct sum of Ext-groups
        H^*(M[X],mu^1) over all objects X in A.
        '''
        A,M=self.cat,self
        coh_gps={X: M.cpx(X).cohomology() for X in A.objects}
        total_coh={i: sum(coh_gps[X].get(i,0)
                          for X in coh_gps)
                   for i in ChainMap(*(coh_gps[X] for X in coh_gps))}
        answer={i: total_coh[i] for i in total_coh if total_coh[i]!=0}
        return answer

    def width(self):
        '''Returns the maximal and minimal degree of an element in M.total().'''
        keys=self.total().keys()
        return min(keys),max(keys)

    def shift(self,m=1):
        '''Returns the A_\infty module shifted in degree by m.'''
        M,A=self,self.cat
        modules={X: M[X].shift(m) for X in M.modules}
        operations={(X,): M.mu(X).rejig_2()
                    for X in A.objects if (X,) in M.operations}
        for word in M.operations:
            F=M.operations[word]
            if len(word)!=1:
                if len(word)==2:
                    new_source=M[word[-1]].shift(m).otimes(A.hom(*word))
                else:
                    new_source=M[word[-1]].shift(m).otimes(A.hom(*word)).flatten(1)
                new_target=M[word[0]].shift(m)
                new_op=LinearMap(new_source,new_target,F.deg)
                new_op.maps.update({i: F[i].shift(m) for i in F.maps})
                operations.update({word: new_op})
        return A8Module(A,modules,operations)

    def oplus(self,other):
        '''Returns the module which is the direct sum of self and other.'''
        (M,N)=(self,other)
        A=M.cat
        new_modules={X: M[X].oplus(N[X])
                     for X in ChainMap(M.modules,N.modules)}
        new_operations={}
        all_keys=ChainMap(M.operations,N.operations)
        zeroNM=A8ModuleMap(N,M,0,{})
        zeroMN=A8ModuleMap(M,N,0,{})
        new_operations.update({word:
                               LinearMap.block(
                                   M.mu(*word),zeroNM.cpt(*word),
                                   zeroMN.cpt(*word),N.mu(*word))
                               for word in all_keys
                               if len(word)==1})
        new_operations.update({word:
                               LinearMap.block(
                                   M.mu(*word),zeroNM.cpt(*word),
                                   zeroMN.cpt(*word),N.mu(*word))
                               .flatten(1).unflatten(0,2)
                               for word in all_keys
                               if len(word)>1})
        new_a8mod=A8Module(A,new_modules,new_operations)
        return new_a8mod
        
    
    def twist(self,X):
        '''Returns the twist of the module M around the object X.

        This is obtained in several steps:
        1. We form
            (a) the Yoneda module  Y = yoneda(X)
            (b) the chain complex  Z = (M(X), \mu^1)
            (c) the tensor product T = Z (x) Y
        2. We form the canonical evaluation morphism, ev:

            ev^d(c(x)b,a_{d-1},...,a_1) = \mu^{d+1}_M(c,b,a_{d-1},...,a_1)
        
        3. We return the cone on ev.
        '''
        Y=self.cat.yoneda(X)
        Z=self.cpx(X)
        T=Z.otimes(Y)
        new_cpts={}
        new_cpts.update({word[:-1]: self.mu(*word)
                         for word in self.operations
                         if word[-1]==X and len(word)==2})
        new_cpts.update({word[:-1]: self.mu(*word).unflatten(0,2)
                         for word in self.operations
                         if word[-1]==X and len(word)>2})
        ev=A8ModuleMap(T,self,0,new_cpts)
        return ev.cone().simplify()

class A8ModuleMap(AlgebraicStructure):
    '''The class of A_\infty pre-module homomorphisms.

    USAGE:

        e = A8ModuleMap(source,target,deg,components)

    creates an A_\infty pre-module map of degree deg from the module
    source to the module target, with the given components.

    ATTRIBUTES:

        e.source, e.target [A8Modules]

            The source and target of e.

        e.deg [int]

            The degree of e.

        e.components [dict] {(X_0,...,X_{d-1}): LinearMap}

            A dictionary indexed by d-tuples from the set
            A.objects. The entry e.components[(X_0,...,X_{d-1})]
            encodes the component

     e^d: source(X_{d-1}) (x) A[(X_{d-1},X_d)] (x) ... (x) A[(X_0,X_1)]
                                    ----> target[(X_0,X_d)]

            of the pre-module map.

            This can also be accessed using

                e.cpt(X_0,...,X_{d-1})

            which will return 0 if the component has not been set.
    

    METHODS:

        e.cpt(X_0,...,X_{d-1})

            Returns e.components[(X_0,...,X_{d-1})] if defined and
            zero otherwise.

        e.display()

            Displays the components of e (mostly for debugging
            purposes).

        e.cone()

            Returns the cone on the pre-module homomorphism e; if e is
            a module map (i.e. closed) then this cone is an A_\infty
            module.

    Currently there is no method for verifying if e is closed.
    '''
    _required_fields=['source','target','deg','components']

    def cpt(self,*word):
        '''Returns e.components[word] if defined and zero otherwise.'''        
        if word in self.components:
            return self.components[word]
        else:
            cpt_source=self.source.mu(*word).source
            cpt_target=self.target[word[0]]
            return LinearMap(cpt_source,cpt_target,1+self.deg-len(word))

    def display(self):
        '''Displays the components of e.'''
        for word in self.components:
            print('F^',len(word))
            if len(word)==1:
                print('M(',word[-1],') = ')
            else:
                print('M(',word[-1],') * A.hom(',word,') = ')
            self.components[word].display()
        
    def cone(self):
        '''Returns the cone on an A_\infty pre-module morphism.

        The modules for the cone are
        C(X) = M(X)[1] (+) N(X)
        and the operations are the block matrices

            (\mu_M  0    )
            (F      \mu_N)

        but F and \mu_M need to be suitably shifted
        (using rejig_*) to make sense and some flattening/
        unflattening is needed to get the bases to match up.
        '''
        (M,N)=(self.source,self.target)
        A=M.cat
        new_modules={X: M[X].shift().oplus(N[X])
                     for X in ChainMap(M.modules,N.modules)}
        new_operations={}
        all_keys=ChainMap(self.components,M.operations,N.operations)
        zero=A8ModuleMap(N,M.shift(),1,{})
        new_operations.update({word:
                               LinearMap.block(
                                   M.mu(*word).rejig_2(),zero.cpt(*word),
                                   self.cpt(*word).rejig_3(),N.mu(*word))
                               for word in all_keys
                               if len(word)==1})
        new_operations.update({word:
                               LinearMap.block(
                                   M.mu(*word).rejig_2(),zero.cpt(*word),
                                   self.cpt(*word).rejig_3(),N.mu(*word))
                               .flatten(1).unflatten(0,2)
                               for word in all_keys
                               if len(word)>1})
        new_a8mod=A8Module(A,new_modules,new_operations)
        return new_a8mod
