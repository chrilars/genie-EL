import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False)
x1 = [1, 2, 3, 4, 5, 6]
x2 = [1, 2, 3, 4, 5, 6, 7, 8, 9]
yr = [3, 9, 31, 129, 651, 3913]
yo = [3, 5, 7, 9, 11, 13, 15, 17, 19]

ax1.plot(x1, yr, 'r--', label="Rabin (n pairs)")
ax2.plot(x2, yo, 'g--', label="Parity (2n colors)")
ax1.set_xlabel("n")
ax1.set_ylabel("Number of nodes")
ax2.set_xlabel("n")
ax2.set_ylabel("Number of nodes")
ax1.legend()
ax2.legend()
fig.set_figwidth(10)
plt.savefig("graph_nodes.png", bbox_inches="tight")
plt.show()
