#! python
from sys import argv
import numpy
import re
import scipy
from bitstring import BitArray
numpy.set_printoptions(threshold='nan')
#print "Hello, Python below is your file"
#filename= raw_input('ENTER YOUR  FIRST Filename: ')
script, filename1, filename2, correlation = argv
print(filename2)
A1 = scipy.fromfile(open(filename1), dtype=scipy.uint8)
#filename2= raw_input('ENTER YOUR  SECOND Filename: ')
A2 = scipy.fromfile(open(filename2), dtype=scipy.uint8)
#print("A1 in file STARTS HERE")
#print(A1)
#print("A2 out file STARTS HERE")
#print(A2)
if A2.__len__()!=A1.__len__():
	print("NOT EQUAL LENGTH")
print( "length of in =" +str(A1.__len__()))
print( "length of out =" +str(A2.__len__()))
test=0
test=int(correlation)
print(" correlation =" +str(test))
count=0;
if(A2.__len__()>A1.__len__()):
    end=A1.__len__()
else: end=A2.__len__()
index=0;
delay=int(correlation)
print("this is the end= "+str(end))
if(delay>=0):
#if correlation is +ve
  for i in range(delay,end-delay):
      if A1[i-delay]!=A2[i]:
            count+=1
            index=i
#if correlation is -ve

else:
    for i in range(0, end+delay):
        if A1[i-delay]!=A2[i]:
            count+=1
            index=i






print("count ="+str(count)+ " index="+ str(index))
print(str(float(count)/float(end-float(delay))) + "    this is end=" +str(end))
po2=re.findall("[/].+[_].+",filename2)
result2=str(po2[0])
SNRin=re.findall("[_].+",result2)
presult=str(SNRin[0])
po= re.findall("\w+[/]",filename2)
result=str(po[0])+ str(po[1])+"BER"
f=open(result,'a')
print("THIS is SNR from BER "+ presult[1:])
data=str(float(count)/float((float(end)-float(delay))))+ ","+ str(presult[1:])
f.write(data+ "\n")
f.close()
#a= BitArray(float=0.34, length=32)
#a.bin=f
#print(a.float)
