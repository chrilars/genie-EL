import matplotlib.pyplot as plt

fig = plt.figure()
x1 = [1, 2, 3, 4, 5, 6]
x2 = [1, 2, 3, 4, 5, 6, 7, 8, 9]
yr = [3, 9, 31, 129, 651, 3913]
yo = [3, 5, 7, 9, 11, 13, 15, 17, 19]

plt.plot(x1, yr, 'r--', label="Rabin (n pairs)")
plt.plot(x2, yo, 'g--', label="Parity (2n colors)")
plt.xlabel("n")
plt.ylabel("Number of nodes")
plt.legend()
plt.savefig("graph_nodes.png")
plt.show()
