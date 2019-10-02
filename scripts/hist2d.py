import pandas as pd
import matplotlib.pyplot as plt

metrics = pd.read_csv('/Users/massoudmaher/data/sc1935_to_1937_metrics.csv')

fig, ax = plt.subplots()
ax.set_aspect("equal")

hist, xbins, ybins, im = plt.hist2d(metrics["total_mapped_reads_hmmcopy"], 
                                    metrics["quality"])                    
                          
for i in range(len(ybins)-1):
    for j in range(len(xbins)-1):
        ax.text(xbins[j]+0.5,ybins[i]+0.5, hist[i,j], 
                color="w", ha="center", va="center", fontweight="bold")

plt.savefig("/Users/massoudmaher/Desktop/hist2d_annot.png")
