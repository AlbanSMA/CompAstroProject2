import pandas as pd
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1,1)
path = "C:\\Users\\User\\Documents\\Documents\\University\\Master_SU\\Year_1\\IIB\\CompAstro\\Projects\\Project2\\Project\\Plots\\"

data = pd.read_csv("data.txt", sep="\\s+", header=None)
data = data.values

get_label = ["time", "p", "n", "D", r"$\gamma$", "T", r"$^3He$", "$^4He$"]
for i in range(1,8):
    ax.plot(data[::10,0], data[::10,1]*data[::10,i], "-", label=get_label[i])
ax.grid(alpha=0.3)
ax.set_xlabel("time (s)")
ax.set_ylabel(r"Density as $n/n_p$")
ax.set_title("Big Bang Nucleosynthesis")
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend()
fig.savefig(path+"0.1_10_5e5.pdf", dpi=300)