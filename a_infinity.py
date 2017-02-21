#!/usr/bin/python

import fields as fi
from linear_algebra import *
from itertools import permutations as perm
from itertools import combinations_with_replacement as comb

class A8Category(AlgebraicStructure):
    _required_fields=['field','objects','morphisms','operations']

    def __getitem__(self,i):
        if i in self.morphisms:
            return self.morphisms[i]
        else:
            return VectorSpace(self.field)
        
    def hom(self,*word):
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
        
        '''# We now have a better implementation
        max_eq=1+max([len(word) for word in A.operations])
        super_words=[]
        for i in range(2,max_eq+1):
            for choice in comb(A.objects,i):
                for word in set(perm(choice)):
                    super_words.append(word)'''

        truth_dictionary={}
        for word in super_words:
            truth_dictionary[word]=check_word(*word)
        if not all(truth_dictionary):
            raise ValueError('A_\infty equations not satisfied')
        else:
            print('A_\infty equations hold')

    def yoneda(self,Q):
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
            return None


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
        
        def pt(K,N):
            V=VectorSpace(K)
            V.basis.update({0:N})
            return V
        
        objects=self.vertices
        arrow=self.arrows
        morphisms={}
        operations={}

        A=A8Category(K,objects,morphisms,operations)
        
        for X in objects:
            A.morphisms[(X,X)],A.operations[(X,X,X)]=sph(K,N)
            
        A.morphisms.update({(X,Y): pt(K,0)
                          for X in objects
                          for Y in arrow[X]})
        A.morphisms.update({(Y,X): pt(K,N)
                          for X in objects
                          for Y in arrow[X]})
        
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
                
        for X in objects:
            for Y in arrow[X]:
                for Z in arrow[Y]:
                    if Z in arrow[X]:
                        F=LinearMap(A.hom(X,Y,Z),A[(X,Z)],0)
                        F.maps.update({(0,0):A[(X,Z)][0]})
                        A.operations.update({(X,Y,Z): F})
        
        return A
                        
def BP(p,q,K,N):
    '''Generates the Fukaya category of a Brieskorn-Pham Milnor fibre
        of the form x^p+y^q+z_1^2+...+z_{N-1}^2=1.'''
    M=(p-1)*(q-1)
    vertices={m for m in range(1,M+1)}
    def nbhd(m):
        return {m+1,m+p,m+p+1} & {x for x in range(1,M+1)}
    arrows={m : nbhd(m) for m in range(1,M+1)}
    G=DynkinGraph(vertices,arrows)
    return G.categorify(K,N)

    
class CochainComplex(AlgebraicStructure):
    _required_fields=['cochains','differential']

    def cohomology(self):
        '''Returns the cohomology of the cochain complex (Z,d).'''
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
        print(self.cochains)
        self.differential.display()
    
    def verify(self):
        d=self.differential
        Z=self.cochains
        d.verify()
        zero_map=LinearMap(Z,Z,2)
        if d.circ(d) != zero_map:
            raise TypeError('Not a cochain complex')
        else:
            print('This is a cochain complex')

    def otimes(self,M):
        '''Returns the tensor product of a chain complex with
        with an A_\infty-module
        '''
        Z,diff=self.cochains,self.differential
        A=M.cat
        modules={X: Z.otimes(M[X]) for X in M.modules}
        operations={(X,): Z.Id().otimes(M.mu(X))+diff.otimes(M[X].sigma())
                    for X in A.objects if (X,) in M.operations}
        for word in M.operations:
            if len(word)>1:
                new_op=Z.Id().otimes(M.mu(*word))
                operations.update({word: new_op.flatten(1).unflatten(0,2)})               
        return A8Module(A,modules,operations)
                
class A8Module(AlgebraicStructure):
    _required_fields=['cat','modules','operations']

    @lazyproperty
    def field(self):
        return self.cat.field
    
    def __getitem__(self,X):
        M=self
        K=M.field
        if X in M.modules:
            return M.modules[X]
        else:
            return VectorSpace(K)

    def mod(self,Q):
        if Q in self.modules:
            return self.modules[Q]
        else:
            # If we haven't bothered to define module[X]
            # then it's the zero vector space, which can be created
            # like this:
            return VectorSpace(self.field)

    def mu(self,*word):
        M=self
        if word in M.operations:
            return M.operations[word]
        else:
            A=M.cat
            d=len(word)
            X,Y=word[d-1],word[0]
            if d>1:
                new_domain=(M[X].otimes(A.hom(*word))).flatten(1)
            else:
                new_domain=M[X]
            return LinearMap(new_domain,M[Y],2-d)

    def shift(self,m=1):
        '''Returns the A_\infty module shifted in degree by m.

        This is achieved by tensoring with K[m].
        '''
        M,A=self,self.cat
        modules={X: M[X].shift(m) for X in M.modules}
        operations={(X,): M.mu(X).rejig_2()
                    for X in A.objects if (X,) in M.operations}
        for word in M.operations:
            if len(word)!=1:
                F=M.operations(word)
                new_source=M.mod(word[-1]).shift(m).otimes(A.hom(*word)).flatten(1)
                new_target=M.mod(word[0]).shift(m)
                new_op=LinearMap(new_source,new_target,F.deg)
                new_op.maps.update({i: F[i].shift(m)})
                operations.update({word: new_op})
        return A8Module(A,modules,operations)
        
        
        '''K=self.field
        cochains=VectorSpace(K)
        cochains.basis.update({0:-m})
        differential=LinearMap(cochains,cochains,1)
        Z=CochainComplex(cochains,differential)
        return Z.otimes(self)'''
        
    def verify(self):
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
            print('A_\infty-module equations hold')

    def cpx(self,Q):
        '''Returns the cochain complex M(Q), \mu^1.'''
        cochains=self.mod(Q)
        differential=self.mu(Q)
        return CochainComplex(cochains,differential)
    
    def display(self):
        for X in self.modules:
            print('M(',X,') = ',self[X])
        for word in self.operations:
            print('M(',word[-1],') * A.hom(',word,') = ')
            self.operations[word].display()

    def twist(self,Q):
        '''Returns the twist of the module M around the object Q.

        This is obtained in several steps:
        1. We form
            (a) the Yoneda module  Y = yoneda(Q)
            (b) the chain complex  Z = (M(Q), \mu^1)
            (c) the tensor product T = Z (x) Y
        2. We form the canonical evaluation morphism, ev:

            ev^d(c(x)b,a_{d-1},...,a_1) = \mu^{d+1}_M(c,b,a_{d-1},...,a_1)
        
        3. We return the cone on ev.
        '''
        Y=self.cat.yoneda(Q)
        Z=self.cpx(Q)
        T=Z.otimes(Y)
        new_cpts={}
        new_cpts.update({word[:-1]: self.mu(*word)
                         for word in self.operations
                         if word[-1]==Q and len(word)==2})
        new_cpts.update({word[:-1]: self.mu(*word).unflatten(0,2)
                         for word in self.operations
                         if word[-1]==Q and len(word)!=2})
        print('checking ev')
        for x in new_cpts:
            new_cpts[x].verify()
        print('done')
        ev=A8ModuleMap(T,self,0,new_cpts)
        return ev.cone()

    '''def simplify(self):
        M=A8Module(self.cat,{},{})
        translator={}
        for X in self.modules:
            M.modules[X]=VectorSpace(M.field)
            V=self.mod(X)
            translator[X]=({x: y for x,y in
                            enumerate(list(V.basis.keys()))})
            M.modules[X].basis.update({x: V.basis[translator[x]]
                                       for x in translator[X]})
        for word in self.operations:
            d=len(word)
            M.operations[word]=LinearMap(M.mu(*word).source,M.mu(*word).target,2-d)
            if d==1:
                M.operations[word].maps.update({x: self.operations[word](})
            else:'''
                
                

class A8ModuleMap(AlgebraicStructure):
    _required_fields=['source','target','deg','components']

    def cpt(self,*word):
        if word in self.components:
            return self.components[word]
        else:
            cpt_source=self.source.mu(*word).source
            cpt_target=self.target.mod(word[0])
            return LinearMap(cpt_source,cpt_target,1+self.deg-len(word))

    def display(self):
        for word in self.components:
            if len(word)==1:
                print('d_M(',word[-1],') = ')
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
        new_modules={X: M.mod(X).shift().oplus(N.mod(X))
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
        print('Verifying the cone')
        new_a8mod.verify()
        print('Done')
        return new_a8mod
        '''
        for word in all_keys:
            (new_source,new_target)=(the_cone.mu(*word).source,
                                     the_cone.mu(*word).target)
            d=len(word)
            the_cone.operations[word]=LinearMap(new_source,new_target,2-d)
            new_keys=ChainMap(M.operations[word].maps,self.components[word].maps)
            new_ops={}
            new_ops.update({(('a',i[0]),)+i[1:]:
                            M.mu(*word)[i].oplus(self.cpt(*word)[i])
                            for i in new_keys if len(i)>1})
            new_ops.update({('a',i[0]):
                            M.mu(*word)[i].oplus(self.cpt(*word)[i])
                            for i new_keys if len(i)==1})
            new_ops.update({(('b',j[0]),)+j[1:]:
                            Vector(M[word[-1],{}).oplus(N.mu(*word)[j])
                            for j in N.operations[word].maps
                            if len(j)>1})
            new_ops.update({('b',j[0]):
                            Vector(M[word[-1]],{}).oplus(N.mu(*word)[j])
                            for j in N.operations[word].maps
                            if len(j)==1})
            the_cone.operations[word].maps.update(new_ops)
        '''


