#! python
from sys import argv
import numpy as np
import numpy
import scipy
from bitstring import BitArray
filename='inputBinary'
numpy.set_printoptions(threshold='nan')
f = scipy.fromfile(open(filename), dtype=scipy.uint8)
#print(f)
print(f.__len__())
x=np.average(f)
print(x*10**-6)
#a= BitArray(float=0.34, length=32)
#a.bin=f
#print(a.float)
st=""
for item in range(0,f.__len__()):
    if(item!=(f.__len__()-1)):
      st=st+str(f[item])+","
    else:
        st=st+str(f[item])
    
    print f[item],
    if(item!=(f.__len__() - 1)):
     print ","



print(f.__len__())
f=open("inputBinaryIN",'w+')
f.write(st)
f.close()
 
