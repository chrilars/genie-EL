# Notes meeting

* Negated normal form maybe, just con- and disjunction and negation.
* Positive boolean formulae
* standard options for parity and buchi and so on
    * For such, generate tree from template

* ls : W means append W
* Crucial to do the recursive uions and interesections the way that mascot does, bc otherwise inefficient 
* Basically switch everything to BDDs
* Maybe fairsyn is the part which uses the formulae


## Questions for meeting
* read about BDDs
* Graph nodes, arena gamma and so on [BDDs]
* data structures for the zielonka tree [use own, but for arena and graph use the symbolicset or bdds]
* how to rearrange the fixpoints to make them appropriate for our cause
* live and nonlive vertices [should not matter]
* cube of variables? [ask kaushik, genie include ubdd baseubdd]
* Grid?


## Answers from Daniel
* Dont convert everything to BDDs, use own structures for the ZTree, but arena and game is in BDD form still because they use bdds for nodes and graph
* The options to do buchi or rabin or something would be from our generalized fixpoint function, not just plain mascot
* Maybe use binary strings for the colors, also check how they treat the rabin pairs
* mu är least
* stort sett samma CPre
* generalisera deras Rabinfunktion till en EM-funktion
* Start work on CUDD
* Ändra i FairSyn basefixpoint
* Ändra Rabinfunktionen i basefixpoint till att ta in ett zielonkatree och modifiera för mer generellt beteende
* Ändra sequentialrabin och rabinconstructor
* Gör en conditionparser (Inf c | Fin c | phi&psi | phi|psi | (psi)) [done]


## Nästa möte
* Mascot and ZTree can be totally independent, we dont even need the graph from mascot
* Keep old unecessary things and just expand - dont remove all old functionality
* Arbeta med genie istället för mascot
* Genie är iallafall delen som löser spelen
* Change in genie/include/lib/basefixpoint
* Do not need "live_"
* Can use fairsyn to run the experiments




## Cpre implementation
```cpp
UBDD cpre(const UBDD &Zi) {
    UBDD Z = Zi;
    /* project onto state alphabet */
    Z = Z.existAbstract(CubeNotState());
    /* swap variables */
    Z = Z.permute(preVars_, postVars_);
    /* the controllable system edges */
    UBDD W0 = tr_ & Z & sys_nodes_;
    /* the controllable environment edges */
    UBDD nZ = !Z & nodes_;
    /* the environment nodes having an outgoing edge outside Zi */
    UBDD F = tr_.andAbstract(nZ, CubeNotState()) & env_nodes_;
    /* the other environment nodes are controllable */
    UBDD nF = !F;
    UBDD W1 = tr_ & nF & env_nodes_;
    /* return all the controllable edges */
    return (W0 | W1);
}
```





## Notes to self
BDD subset: S subset of T -> if T & S == S
BDD intersection: S & T
BDD union: S | T

set of BDDs intersection:
for s in S:
    if s in T:
        W.push(s)

set of BDDs union:
for s in S:
    if s not in T:
        (copy of T or something) T.push(s)

set of BDDs diff:
(S-T)
for s in S:
    if s not in T:
        W.push(s)

## Report
* Future work
    * Fairness of EM games



## Improvement
* yac or lex

## Theory
* Zielonkatree
* emerson-lei game
* Arena
* All concepts we are using
* Main results of the paper - just refer to daniels paper about the proofs and so on
* "this is how it works in theory, now to how we implement it"


Given set of colors and convert to BDD that encodes that includes at least one of those colors


## Meeting
* BDD verkar finnas constructor för från set
* Think about how testing should be done
* In the algorithm we pass a whole list of sets, but we want to just pass along the term with the combined BDDs
* So dont pass the W, pass along the term1 which we have discussed, on the leaf - just return the argument


## TODO
* Send mail to yahia
* Maybe prove that the term1, term2, term3 chain is correct instead of strategy ext
