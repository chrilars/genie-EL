# Symbolic solutions of emerson-lei games for reactive synthesis notes

## Sets
- $T_\circ$: Winning Zielonka nodes for circle player
- $T_\square$: Winning Zielonka nodes for square player
- $R(s)$: Set of nodes $\{t\in\ T \mid (s,t)\in\ R\}$, direct successors
- $\gamma(v)$: Set of colors at node $v\in\V$
- $V$: Nodes in game


## Questions
- How does Zielonka trees work when there are a mix of finite and infinite conditions in the winning conditions. Example 2, page 5.

- Explain $anc^s_t$ and $\gamma^{-1}$

- Explain the arguments of Solve

- Are there Zielonka trees built into Maskot or do we need to make our own? Cannot find in docs.


## Algorithm
Solve(s: Vertex of the Zielonka tree, ls: List of subsets from V)

Initialize $X_s$ appropriately according to if the node $s$ is winning for circle or square player in the Zielonka tree.

Make a set $W$ with nodes in $V$ but not in $X_s$.

While $X_s\neq\ W$, do

$W \leftarrow\ X_s$
