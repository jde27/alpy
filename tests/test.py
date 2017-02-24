#!/usr/bin/python

from collections import ChainMap

dict={}
dict[1]={0:1, 2:2, 3:4}
dict[2]={0:2, 4:1, 2:1}
dict[3]={0:1, 3:1, 5:1}
#keys=set().union(*(set(dict[X]) for X in dict))
total={i: sum(dict[X].get(i,0) for X in dict)
       for i in ChainMap(*(dict[X] for X in dict))}
print(total)
