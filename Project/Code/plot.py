import pandas as pd
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1,1)

data = pd.read_csv("data.txt", sep="\s+", header=None)
print(data[7])

get_label = ["time", "p", "n", "D", r"$\gamma$", "T", r"$^3He$", "$^4He$"]
for i in (1, 2, 3, 4, 5, 6, 7):
    ax.plot(data[0], data[i], label=get_label[i])
ax.grid(alpha=0.3)
ax.set_xlabel("time (s)")
ax.set_ylabel(r"number density ($g/cm^-3$)")
ax.set_title("Big Bang Nucleosynthesis")
ax.set_yscale("log")
ax.set_xscale("log")
ax.legend()
plt.show()