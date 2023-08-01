# Master-Thesis

Learning Graphical Lyapunov models with best-subset selection methods

Author: Rahul Radhakrishnan, Technical University of Munich, Germany 

Supervisor: Dr. Mathias Drton, Department of Mathematical statistics, Technical University of Munich, Germany

Advisor: Msc. Philipp Dettling, Department of Mathematical statistics, Technical University of Munich, Germany

## Abstract
Cause and effect relationships between variables involved in a dynamical process are of great interest. Existing graphical model approaches allow limited exploration of relationships. In particular, cyclic interactions between different components that represent feedback loops are generally harder to capture and also to interpret using such methods. In this light, graphical Lyapunov models offer a new perspective on modeling causally interpretable dependence structure in multivariate data by treating each independent observation as a one-time cross-sectional snapshot of a multivariate Ornstein-Uhlenbeck processes in equilibrium. The Gaussian equilibrium exists under a stability assumption on the drift matrix, and the equilibrium covariance matrix is determined by the continuous Lyapunov equation. The models assume sparsity on the drift matrix and its support is determined by a directed graph. This directed graph, in particular, can include directed cycles and is one of the main advantage of this approach. Past works have shown one approach to this support recovery by using an $l_1-$ regularized approach (lasso) that seeks to find sparse approximate solution to the continuous Lyapunov equation, given a sample covariance matrix, called the \textit{Direct Lyapunov Lasso}. 

This work focuses on recovering the support of the drift matrix from the Lyapunov equation using best-subset selection (BS). The classical methods of solving best-subset selection, a NP-hard problem, are not scalable to large problem sizes or cannot guarantee global optimality. One of the recent works focuses on converting best-subset selection to a mixed integer optimization (MIO) problem. With advancement in optimization techniques and parallel computing, mixed integer programming has become very powerful. Its practical applicability and feasibility relies on the strength of the formulation, warm starts and their quality. Out of the two formulations presented for best-subset selection, the stronger and computationally more efficient formulation was chosen in this study. Furthermore, the formulations were adapted to recover stable drift matrices. The MIO formulation was solved using Gurobi on its R interface. 

For the warm starts, four different initializations were explored: three variants of projected gradient descent solutions and direct lasso solutions ($\mathrm{BS_{lasso}}$). The projected gradient solutions achieve a $k-$sparse solution using marginal regression coefficients as initializations. Its variants differ in terms of the strategy used to supply these marginal regression coefficients to the algorithm as a starting point. In particular, one can directly choose the $k$ edges with largest coefficients ($\mathrm{BS_{reg}}$), choose only 1 edge ($\mathrm{BS_{1edge}}$), and choose edges based on sparsity pattern of the estimate of inverse covariance matrix ($\mathrm{BS_{glasso}}$). 

All the initialization strategies are compared against each other and to the Direct Lyapunov Lasso on synthetic and real dataset.  Information criterion are used to select drift matrices that are best fitting to the observed data, while balancing the sparsity. The results show that lasso outperforms BS, provided there are enough samples. BS has the potential to perform better than lasso for low sample size, especially for large graph structures. One also finds in this study that lasso solutions seem to provide better bounds to the MIO formulation compared to projected gradient descent solutions. So $\mathrm{BS_{lasso}}$ performs as good as or better than the other BS methods. All the projected gradient descent based methods have similar performance, with $\mathrm{BS_{1edge}}$ and $\mathrm{BS_{glasso}}$ being marginally better. These comparisons are relative. It is shown that all the methods perform very poor beyond certain problem size and below a certain value of sample size.

## Code files
Simulated data: 

Run recoverLyapunovGraphs.R for generating signals and recover them. 


Isoprenoid data: 

Run recoverPathway.R for loading the isoprenoid dataset and recover the support of the drift matrix (on the assumption that the volatility matrix is identity - Future works need to focus on joint estimation of M and C). 
Data was obtained from: Wille, Anja, et al. "Sparse graphical Gaussian modeling of the isoprenoid gene network in Arabidopsis thaliana." Genome biology 5.11 (2004): 1-13.
