# MEM-Rearrange
<!-- (-   [Overview](#overview) -->
-   [Overview](#overview)
-   [Input Format](#if)
<!--
<a name='overview'>Overview</a>
--------

-->
<a name='overview'>Overview</a>
--------
MEM-Rearrange is an implementation of the algorithms described in the paper **New Algorithms for Structure Informed Genome
Rearrangement**.
The data used in the expirement of the paper, is placed in the folder "input_families".
In the paper we define two problems, *Constrained TreeToString Divergence (CTTSD)* and *TreeToString Divergence (TTSD)*.

The input to *CTTSD* consists of two signed permutations of length $n$, $S_1$ and $S_2$ ; a PQ-tree $T$ ordered as $S_1$ ; and two numbers $\delta^Q_{\mathsf{ord}}$ and $\delta^Q_{\mathsf{flip}}$ indicating the penalty of the events of changing order and flipping, respectively, a Q-node.. We aim to perform actions on $T$ to reorder it as $S_2$. That is, we reorder $T$ as $T'$ so that $F(T')=S_2$ and we return the divergence from $T$ to $S_2$.

Generalizing *CTTSD*, in *TTSD* we do not assume that the input strings are permutations, and we allow deletions. The input to *TTSD* consists of two signed strings, $S_1$ and $S_2$; a PQ-tree $T$ ordered as $S_1$; $d_T\in\mathbb{N} \cup \{0\}$, which specifies the number of allowed deletions from $T$; $d_S\in\mathbb{N}\cup \{0\}$, which specifies the number of allowed deletions from $S_2$; and two numbers $\qpenaltyorder$, $\qpenaltyflip$ indicating the penalty of the events of changing order and flipping, respectively, a Q-node.
In *TTSD* we perform actions on $T$ to reorder it as a subsequence $S_2'$ of $S_2$, allowing $d_T$ deletions from $T$, and so that $S_2'$ is obtained from $S_2$ by using up to $d_S$ deletions from $S_2$. That is, after reordering $T$ as $T'$ (and performing up to $d_T$ deletions), $F(T')=S_2'$ and we return the divergence from $T$ to $S_2'$ corresponding to $d_T$ and $d_S$.


<a name='if'>Input Format</a>
--------
The input files format is as follows.
For each CSB in the cluster, there is a line in the next format:

<CSB_ID><\t><Length><\t><Score><\t>Instance_Count><\t><CSB><\t><Main_Category><\t><Family_ID>
Please see an example in the folder "input_families".


