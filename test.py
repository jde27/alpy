#!/usr/bin/python
D={(1,):6,(2,):7,(1,2):10}
Q=2
print({w: D[w+(Q,)]
       for w in D
       if w+(Q,) in D})
