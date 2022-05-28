import matplotlib.pyplot as plt
import numpy as np

#time = [1.16, 1.83, 4.42, 3.43]
louvain_time = [1757, 699, 2372, 163, 104]
#stan_time = [66951, 59103, 37093, 4298, 11959, 3901, 3197, 3070, 3311, 4259]
stan_time = [66951, 59103, 37093, 4298, 3901, 3197, 3070, 3311, 4259, 11959] 
x = np.arange(len(stan_time))
#labels = ["Opt-High", "Opt-Normal", "Opt-Low", "No-Opt"]
louvain_labels = ["s-normal", "s-opt-high", "s-opt-low", "m-normal", "m-opt-high"]
stan_labels = ["s-normal", "s-opt", "s-cut-3", "m-normal","m-cut-1", "m-cut-2", "m-cut-3", "m-cut-4", "m-cut-5", "m-opt"]
width = 0.4

plt.bar(x, stan_time, width=width, align='center', tick_label=stan_labels)
plt.yscale("log")
plt.show()
