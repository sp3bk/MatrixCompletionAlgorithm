import numpy as np
import pandas as pd
import csv
from sklearn import Imputer

def main():
   imp = Imputer(strategy="mean")
   #a=csv.reader(open("x.csv","rb"),delimiter=',')
   a = pd.read_csv('x.csv')
   #rint("done")
   #print (a)
   print("done")
   a = imp.fit_transform(a)
   #print (a)
   #____________________
   c = csv.writer(open("x4.csv", "wb"))
   listA=[]
   for x in range(0,122):
      listA.append(x)  
   c.writerow(sorted(listA))
   for x in range(0,121):
      c.writerow(a[x])
      
main()