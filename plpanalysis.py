import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("zap", header=None, sep=" ",usecols=[1,2,3,4,5,6,7,8,9,10],
                 names=["b4","b5","H","plq","P4r","P4i","P4","P5r","P5i","P5"])

df.b = np.sqrt(df.b4 * df.b5)

hvals = np.unique(df.H)
bvals = np.unique(df.b)

for b in bvals:
#    for h in hvals:
    for h in [0]:
        dfbh = df[(df.H == h) & (df.b == b)]
        plt.plot(dfbh.P5r, dfbh.P5i, "r.", markersize=5, label="b:"+str(b)+" h"+str(h), alpha=0.5 )
        plt.plot(dfbh.P4r, dfbh.P4i, "b.", markersize=5, label="b:"+str(b)+" h"+str(h), alpha=0.5 )


        plt.plot([-1/3, np.cos(np.pi/3)/3, np.cos(-np.pi/3)/3], 
         [0.0, np.sin(np.pi/3)/3, np.sin(-np.pi/3)/3], "ko", alpha=0.3)
        plt.plot([1, np.cos(2*np.pi/3), np.cos(-2*np.pi/3)],
         [0.0, np.sin(2*np.pi/3), np.sin(-2*np.pi/3)], "k^", alpha=0.3)
        plt.plot(0,0, "k*", markersize=12, alpha=0.2)

        plt.axis('equal')
        plt.xlim(-2,2)
        plt.ylim(-2,2)
        plt.legend(numpoints=1,bbox_to_anchor=(0.8, 1.1), loc=2, borderpad=1,markerscale=10)

        plt.show()
