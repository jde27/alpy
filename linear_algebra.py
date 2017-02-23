#!/usr/bin/python

import fields as fi
from collections import Counter, ChainMap

class AlgebraicStructure:
    _required_fields=[]
    _empty_dictionaries=[]
    def __init__(self,*args):
        L=len(self._required_fields)
        if len(args) != L:
            raise TypeError('Expected {} arguments'.format(L))
        for name, value in zip(self._required_fields,args):
            setattr(self,name,value)
        for name in self._empty_dictionaries:
            setattr(self,name,{})

    @staticmethod
    def compat(objects,*args):
        truth=True
        for test in args:
            check_list=[getattr(X,test) for X in objects]
            if check_list.count(check_list[0])!=len(check_list):
                truth=False
                raise TypeError('Objects incompatible for '+
                                'attribute {}'.format(test))
        return truth

    
class lazyproperty:
    def __init__(self,func):
        self.func=func
        
    def __get__(self,instance,cls):
        if instance is None:
            return self
        else:
            value=self.func(instance)
            setattr(instance,self.func.__name__,value)
            return value
            
class VectorSpace(AlgebraicStructure):
    '''Class of vector spaces

    V=VectorSpace(K) creates a vector space over the field K.

    Attributes:
    
        V.field: An instance of fi.Field
        V.basis: A dictionary of vectors.

    Methods:

        V.gr_dim(n): Returns the dimension of the nth graded piece.

        V.oplus(W): Returns the direct sum of V and W.
    
        V.otimes(W): Returns the tensor product of V and W.

        V.Id(): Returns the identity map V--->V.

        V.build(graded_dim): Builds a graded vector space with
                             graded pieces of the specified dimensions;
                             namely, {d:n} will give an n-dimensional
                             piece in degree d.
    '''
    _required_fields=['field']
    _empty_dictionaries=['basis']

    def __getitem__(self,i):
        return Vector(self,{i:self.field(1)})
    
    def __eq__(self,other):
        if type(other) is VectorSpace:
            return (self.basis==other.basis)
        elif other==0:
            return not self.basis

    def __str__(self):
        return "%s" % str(self.basis)

    def flatten(self,k):
        '''Given a vector space V whose basis elements are
        indexed by tuples (i_0,...,i_n) where i_k is itself
        a tuple (j_0,...,j_m), V.flatten(k) returns
        the same vector space, but with its basis indexed by
        tuples (i_0,...,i_{k-1},j_0,...,j_m,i_{k+1},...,i_n).
        '''
        V=VectorSpace(self.field)
        for x in self.basis:
            if k!=0:
                if k+1!=len(x):
                    V.basis.update({x[:k]+x[k]+x[k+1:]: self.basis[x]})
                else:
                    V.basis.update({x[:k]+x[k]: self.basis[x]})
            else:
                if k+1!=len(x):
                    V.basis.update({x[k]+x[k+1:]: self.basis[x]})
                else:
                    V.basis.update({x[k]: self.basis[x]})
        return V

    def unflatten(self,m,n):
        '''Given a vector space V whose basis elements are
        indexed by tuples (i_0,...,i_k) V.unflatten(m,n) returns
        the same vector space, but with its basis indexed by
        tuples (i_0,...,i_{m-1},(i_m,...,i_{n-1}),i_n,...,i_k).
        '''
        V=VectorSpace(self.field)
        for x in self.basis:
            if m!=0:
                if n!=len(x):
                    V.basis.update({x[:m]+(x[m:n],)+x[n:]: self.basis[x]})
                else:
                    V.basis.update({x[:m]+(x[m:n],): self.basis[x]})
            else:
                if n!=len(x):
                    V.basis.update({(x[m:n],)+x[n:]: self.basis[x]})
                else:
                    V.basis.update({(x,): self.basis[x]})
        return V

        
    
    def build(self,graded_dim):
        '''Builds up the basis as specified by the given dictionary;
        for example, V.build({0:1,2:2}) creates a basis with one vector
        in degree 0 and two in degree 2.
        '''
        ctr=0
        for d in graded_dim:
            self.basis.update({i: d for i in range(ctr,ctr+graded_dim[d])})
            ctr+=graded_dim[d]
        return self

    @lazyproperty
    def gr_dim(self):
        '''Returns a dictionary of gradings of V together with the dimensions
        of the graded pieces of V in that grading.'''
        return Counter(self.basis.values())

    @lazyproperty
    def graded_pieces(self):
        '''Returns a dictionary of vector spaces which are the graded
        pieces of V.
        '''
        graded_pieces={}
        for n in self.gr_dim:
            graded_pieces[n]=VectorSpace(self.field)
            graded_pieces[n].basis.update({i: n for i in self.basis
                                          if self.basis[i]==n})
        return graded_pieces
    
    def shift(self,n=1):
        W=VectorSpace(self.field)
        W.basis.update({i: self.basis[i]-n for i in self.basis})
        return W
    
    def oplus(self,other):
        if AlgebraicStructure.compat((self,other),'field'):
            U=VectorSpace(self.field)
            U.basis.update({('a',i): self.basis[i] for i in self.basis})
            U.basis.update({('b',j): other.basis[j] for j in other.basis})
            return U

    def otimes(self,other):
        return VectorSpace.tensor(self,other)

    @staticmethod
    def tensor(*args):
        '''Given a tuple of vector spaces args=(X_1,...,X_n),
        tensor(*args) returns the tensor product X_1*X_2*...*X_n
        whose basis is indexed by tuples (x_1,...,x_n) with
        x_i a basis element of X_i.
        '''
        if len(args)>1:
            if AlgebraicStructure.compat(args,'field'):
                U=VectorSpace(args[0].field)
                build_basis=[()]
                for X in args:
                    build_basis=[x+(y,) for x in build_basis for y in X.basis]
                    
                def grading(elt):
                    '''Returns the grading |elt| of elt=(x_1,...,x_n)
                    where x_i in X_i and |elt|=\sum_{i=1}^n |x_i|.
                    '''
                    ans=0
                    for i, val in enumerate(elt):
                        ans+=args[i].basis[val]
                    return ans
                
                U.basis.update({x: grading(x) for x in build_basis})
                return U
        else:
            return args[0]
        
    def Id(self):
        '''Returns the identity map V-->V.'''
        V=self
        F=LinearMap(V,V,0)
        F.maps.update({i: V[i] for i in V.basis})
        return F

    def sigma(self):
        '''Returns the linear map b |--> (-1)^{|b|-1} b'''
        V=self
        F=self.Id()
        F.maps.update({i: -self[i]
                       for i in self.basis
                       if not (self.basis[i])%2})
        return F


    
class Vector(AlgebraicStructure):
    '''Class of vectors.'''
    _required_fields=['space','components']

    def __getitem__(self,i):
        if i in self.components:
            return self.components[i]
        else:
            return self.space.field(0)

    def __str__(self):
        s=''
        for i in self.components:
            s=s+' %s e_{%s} ' % (self.components[i],i)
        return '%s' % s
        
    def __add__(self,other):
        if AlgebraicStructure.compat((self,other),'space'):
            sum_cpts=ChainMap(self.components,other.components)
            return Vector(self.space,{i: self[i]+other[i] for i in sum_cpts}).chomp()
            
    def __sub__(self,other):
        if AlgebraicStructure.compat((self,other),'space'):
            sum_cpts=ChainMap(self.components,other.components)
            return Vector(self.space,{i: self[i]-other[i]
                                      for i in sum_cpts}).chomp()

    def __mul__(self,t):
        T=self.space.field(t)
        return Vector(self.space,{i: T*(self[i])
                                  for i in self.components}).chomp()
    
    def __rmul__(self,t):
        return self*t

    def __eq__(self,other):
        if type(other) is Vector:
            return not (self-other).components
        elif other==0:
            return not self.components


    def __neg__(self):
        return Vector(self.space,{i: -self[i] for i in self.components})

    def flatten(self,k):
        '''Given a vector v in V, v.flatten(k) returns the
        corresponding vector in V.flatten(k).
        '''
        v=Vector(self.space.flatten(k),{})
        for x in self.components:
            if k!=0:
                if k+1!=len(x):
                    v.components.update({x[:k]+x[k]+x[k+1:]: self[x]})
                else:
                    v.components.update({x[:k]+x[k]: self[x]})
            else:
                if k+1!=len(x):
                    v.components.update({x[k]+x[k+1:]: self[x]})
                else:
                    v.components.update({x[k]: self[x]})
        return v
    
    def unflatten(self,m,n):
        '''Given a vector v in V, v.unflatten(m,n) returns the
        corresponding vector in V.unflatten(m,n).
        '''
        v=Vector(self.space.unflatten(m,n),{})
        for x in self.components:
            if m!=0:
                if n!=len(x):
                    v.components.update({x[:m]+(x[m:n],)+x[n:]: self[x]})
                else:
                    v.components.update({x[:m]+(x[m:n],): self[x]})
            else:
                if n!=len(x):
                    v.components.update({(x[m:n],)+x[n:]: self[x]})
                else:
                    v.components.update({(x,): self[x]})
        return v

    @lazyproperty
    def field(self):
        return self.space.field
    
    @lazyproperty
    def deg(self):
        gradings={self.space.basis[i] for i in self.components}
        if len(gradings) >= 1:
            raise TypeError('Vector of mixed grading')
            return None
        elif len(gradings)==0:
            return 0
        else:
            return list(gradings)[0]            

    def shift(self,m=1):
        '''Returns the same vector in a shifted vector space.'''
        return Vector(self.space.shift(m),self.components)
        
    def otimes(self,other):
        return Vector.tensor(self,other)

    @staticmethod
    def tensor(*args):
        '''Given a tuple of vectors args=(x_1,...,x_n),
        tensor(*args) returns the tensor x_1*...*x_n.'''
        K=args[0].space.field
        if len(args)>1:
            new_cpts=[()]
            for v in args:
                new_cpts=[x+(y,) for x in new_cpts for y in v.components]
                
            def cpt(x):
                ans=K(1)
                for i, val in enumerate(x):
                    ans*=args[i].components[val]
                return ans
            
            tensor_space=VectorSpace.tensor(*tuple(X.space for X in args))
            return Vector(tensor_space,{x: cpt(x) for x in new_cpts}).chomp()
        else:
            return args[0]
    
    def oplus(self,other):
        x=Vector(self.space.oplus(other.space),
                 {('a',i): self[i] for i in self.components})
        x.components.update({('b',j): other[j] for j in other.components})
        return x
    
    def chomp(self):
        null_cpts={x:y for x,y in self.components.items() if y==0}
        for i in null_cpts:
            del self.components[i]
        return self

class LinearMap(AlgebraicStructure):
    '''Class of linear maps

    F=LinearMap(V,W,d) creates a linear map of degree d between V and W

    Attributes:

        F.source,F.target
        F.maps

    Methods:

        F[i]

        
    '''
    _required_fields=['source','target','deg']
    _empty_dictionaries=['maps']

    def __getitem__(self,i):
        '''Given a basis element e_i in the source of F,
        F[i] returns the vector F(e_i).'''
        if i in self.maps:
            return self.maps[i]
        else:
            return Vector(self.target,{})
    
    def __call__(self,other):
        '''Evaluates a function on a vector.'''
        if other.space==self.source:
            w=Vector(self.target,{})
            for i in other.components:
                w+=other[i]*self[i]
            return w.chomp()
        else:
            raise TypeError('Cannot apply this map to this vector')    


    def verify(self):
        for i in self.maps:
            if i not in self.source.basis:
                print("Map failed: source bases don't match")
            if self[i].space!=self.target:
                print("Map failed: target doesn't match value:")
                print(self[i].space)
                print(self.target)
                exit()

            
    @lazyproperty
    def field(self):
        '''Returns the field over which the linear map is linear.'''
        return self.source.field
    
    def flatten(self,k):
        '''Given a linear map V-->W, F.flatten(k) returns
        the corresponding map F:V.flatten(k)-->W.
        '''
        F=LinearMap(self.source.flatten(k),self.target,self.deg)
        for x in self.maps:
            if k!=0:
                if k+1!=len(x):
                    F.maps.update({x[:k]+x[k]+x[k+1:]: self[x]})
                else:
                    F.maps.update({x[:k]+x[k]: self[x]})
            else:
                if k+1!=len(x):
                    F.maps.update({x[k]+x[k+1:]: self[x]})
                else:
                    F.maps.update({x[k]: self[x]})
        return F
        
    def unflatten(self,m,n):
        '''Given a linear map F:V-->W, F.unflatten(m,n) returns
        the corresponding map F:V.unflatten(m,n)-->W.
        '''
        F=LinearMap(self.source.unflatten(m,n),self.target,self.deg)
        for x in self.maps:
            if m!=0:
                if n!=len(x):
                    F.maps.update({x[:m]+(x[m:n],)+x[n:]: self[x]})
                else:
                    F.maps.update({x[:m]+(x[m:n],): self[x]})
            else:
                if n!=len(x):
                    F.maps.update({(x[m:n],)+x[n:]: self[x]})
                else:
                    F.maps.update({(x,): self[x]})
        return F
    
    def restrict(self,V):
        '''Returns the restriction of a linear map to a subspace of its domain.'''
        subbasis_set={i for i in V.basis}
        basis_set={i for i in self.source.basis}
        if subbasis_set<basis_set:
            new_map=LinearMap(V,self.target,self.deg)
            new_map.maps.update({i: self[i] for i in subbasis_set})
            return new_map
        else:
            raise TypeError('Can only restrict linear map to a subspace.')
        

    def circ(self,other):
        '''Returns the composition of F with G.'''
        H=LinearMap(other.source,self.target,self.deg+other.deg)
        if other.target==self.source:
            H.maps.update({i: self(other[i]) for i in other.maps})
            return H.chomp()
        else:
            raise TypeError('Cannot compose these maps')

    def rejig_1(self,m=1):
        '''This implements the canonical isomorphism
        hom(M,N)[m] --> hom(M,N[m])
        '''
        (M,N,d)=(self.source,self.target,self.deg)
        new_map=LinearMap(M,N.shift(m),d-m)
        for i in self.maps:
            new_map.maps[i]=self.maps[i].shift(m)
        return new_map
        
    def rejig_2(self,m=1):
        '''This implements the canonical isomorphism
        hom(M,N)[d]-->hom(M[m],N[m])[d]
        f--------->((-1)**(m(d-1))) f
        '''
        K=self.field
        if not m%2:
            (M,N,d)=(self.source,self.target,self.deg)
            new_map=LinearMap(M.shift(m),N.shift(m),d)
            new_map.maps.update({i: (self.maps[i].shift(m))*K(-1)
                                 for i in self.maps})
            return new_map
        else:
            (M,N,d)=(self.source,self.target,self.deg)
            new_map=LinearMap(M.shift(m),N.shift(m),d)
            new_map.maps.update({i: self.maps[i].shift(m)
                                 for i in self.maps})
            return new_map

    def rejig_3(self,m=1):
        '''This implements the canonical isomorphism
        hom(M,N)-->hom(M[m],N)[m]
        Remember there are two different ways to implement this
        which differ by a sign - the other is rejig_2(rejig_1(self,-n),n)
        '''
        return (self.rejig_2(m)).rejig_1(-m)
        
    def __eq__(self,other):
        '''Tests if the difference of two linear maps is zero.'''
        if type(other) is LinearMap:
            return not (self-other).maps
        elif other==0:
            return not self.maps

    def __add__(self,other):
        '''Adds two linear maps with the same source, target and degree.'''
        if AlgebraicStructure.compat((self,other),'source','target','deg'):
            sum_maps=ChainMap(self.maps,other.maps)
            H=LinearMap(self.source,self.target,self.deg)
            H.maps.update({i: self[i]+other[i] for i in sum_maps})
            return H.chomp()

    def __sub__(self,other):
        '''Subtracts two linear maps with the same source, target and degree.'''
        if AlgebraicStructure.compat((self,other),'source','target','deg'):
            sum_maps=ChainMap(self.maps,other.maps)
            H=LinearMap(self.source,self.target,self.deg+other.deg)
            H.maps.update({i: self[i]-other[i] for i in sum_maps})
            return H.chomp()

    def __mul__(self,t):
        '''Rescales a linear map by a scalar.'''
        if AlgebraicStructure.compat((self.space,t),'field'):
            T=self.source.field(t)
            H=LinearMap(self.source,self.target,self.deg)
            H.maps.update({i: T*(F[i]) for i in F.maps})
            return H.chomp()
        
    def __rmul__(self,t):
        '''Rescales a linear map by a scalar.'''
        return self*t

    def chomp(self):
        '''Linear maps are stored as dictionaries listing their nonzero values
        on basis elements. Sometimes the value ends up being zero and we no
        longer want to store it; chomp() removes all entries with zero value.'''
        null_maps={x:y for x,y in self.maps.items() if y==0}
        for i in null_maps:
            del self.maps[i]
        return self

    def otimes(self,other):
        '''Returns the tensor product of two linear maps,
        so, given:

        F:V-->W, F':V'-->W'

        it returns (F (x) F')(a (x) b)=F(a) (x) F(b).
        '''
        return LinearMap.tensor(self,other)

    @staticmethod
    def tensor(*args):
        '''Given a tuple of linear maps args=(f_1,...,f_n),
        tensor(*args) returns the map f_1*f_2*...*f_n.
        '''
        if len(args)>1:
            new_source=VectorSpace.tensor(*tuple(f.source for f in args))
            new_target=VectorSpace.tensor(*tuple(f.target for f in args))
            new_deg=sum(f.deg for f in args)
            F=LinearMap(new_source,new_target,new_deg)
            new_maps=[()]
            for f in args:
                new_maps=[x+(y,) for x in new_maps for y in f.maps]
                
            def result(x):
                vec_tup=tuple(args[i].maps[val] for i, val in enumerate(x))
                return Vector.tensor(*vec_tup)
            
            F.maps.update({x: result(x) for x in new_maps})
            return F.chomp()
        else:
            return args[0]
    
    @staticmethod
    def block(A,B,C,D):
        '''Given linear maps:
        A:V-->W, D:V'-->W',
        B:V-->W', C:V'-->W
        LinearMap.block(A,B,C,D) returns a linear map
        (A B): V(+)V' --> W(+)W'
        (C D)
        '''
        V1,W1=A.source,A.target
        V2,W2=D.source,D.target
        X2,Y1=B.source,B.target
        X1,Y2=C.source,C.target
        # Raise error if X1\neq V1 etc
        E=LinearMap(V1.oplus(V2),W1.oplus(W2),A.deg)
        first_keys=ChainMap(A.maps,C.maps)
        second_keys=ChainMap(B.maps,D.maps)
        E.maps.update({('a',i):(A[i]).oplus(C[i]) for i in first_keys})
        E.maps.update({('b',j):(B[j]).oplus(D[j]) for j in second_keys})
        return E

    def display(self):
        '''Displays linear map by printing its nonzero values on a basis.'''
        for i in self.maps:
            print("F(",i,") = ",self[i])

    def ker_im(self):
        '''Implements Gaussian elimination to find the nullity and rank
        of a linear map.'''
        F,V=self,self.source
        F.chomp()
        candidates=[V[i] for i in V.basis]
        kernel=[]
        image=[]
        
        def ker_pop():
            casualties=[]
            survivors=[]
            for v in candidates:
                if F(v)==0:
                    casualties.append(v)
                else:
                    survivors.append(v)
            kernel.extend(casualties)
            return survivors
        
        candidates=ker_pop()
        while candidates:
            x=candidates[0]
            image.append(x)
            if len(candidates)!=1:
                m=list(F(x).components.keys())[0]
                new_candidates=candidates[1:]
                candidates=[(y-(F(y)[m]/F(x)[m])*x).chomp()
                            for y in new_candidates]
                candidates=ker_pop()
            else:
                candidates=[]
        return len(kernel),len(image)
