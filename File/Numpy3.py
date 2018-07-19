import numpy as np
import pandas as pd
import csv
from sklearn.preprocessing import Imputer

def main():
   imp = Imputer(strategy="most_frequent")
   a = pd.read_csv('x.csv')
   a = imp.fit_transform(a)
   c = csv.writer(open("x6.csv", "wb"))
   listA=[]
   for x in range(0,122):
      listA.append(x)  
   c.writerow(sorted(listA))
   for x in range(0,121):
      c.writerow(a[x])
      
   
main()