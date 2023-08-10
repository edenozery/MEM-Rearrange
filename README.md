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
The input to CTTSD consists of two signed permutations of length $n$, $S_1=\sigma_1 \dotso \sigma_n \in \Sigma^n$, $|\Sigma|=n$, such that $\sigma_i \neq \sigma_j$ for all $1 \leq i<j \leq n$, and $S_2=\lambda_1 \dotso \lambda_n \in \Sigma^n$ such that $\lambda_i \neq \lambda_j$ for all $1 \leq i<j \leq n$; a \pqt{} $T$ ordered as $S_1$ with $m_p$ P-nodes and $m_q$ Q-nodes; and two numbers $\qpenaltyorder$ and $\qpenaltyflip$. We aim to perform actions on $T$ to reorder it as $S_2$. That is, we reorder $T$ as $T'$ so that $F(T')=S_2$. We return the divergence from $T$ to $S_2$.


<a name='if'>Input Format</a>
--------
The input files format is as follows.
For each CSB in the cluster, there is a line in the next format:

<CSB_ID><\t><Length><\t><Score><\t>Instance_Count><\t><CSB><\t><Main_Category><\t><Family_ID>
Please see an example in the folder "input_families".


