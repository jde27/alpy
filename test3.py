#!/usr/bin/python

class fancyclass:
    def __init__(self):
        pass
    
    def getargs(self,*args):
        arr={(1,1):
             'String',(1,2):'Other'}
        print(args)
        if args in arr:
            print(arr[args])
P=fancyclass()
P.getargs(1,2)
