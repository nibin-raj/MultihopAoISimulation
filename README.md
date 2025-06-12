
# Sampling and scheduling policies for minimizing age-of-information in multi-hop wireless networks
Submitted for review

---

### Multihop AoI Simulator

This repository contains the source code for the Multihop Age-of-Information (AoI) Index Policy Simulator, designed to validate and reproduce key results from Tripathi et al., "Information Freshness in Multihop Wireless Networks" [[1]](#1).

---

## Age-Difference Policy and Simulator Validation Results

This section presents simulation results obtained using our simulator to reproduce the Age-Difference policy results from [1] and to compare them against our proposed Index policy, denoted as $\pi_{\text{Index}}$.

We consider a network where each link has a transmission success probability $p_t$, and failure probability $1 - p_t$. For the results presented here, we assume reliable links, i.e., $p_t = 1$.

---

### Validation Results

In [[1]](#1), two interference models were analyzed:

1. $K$-hop Interference Model (with $K = 1$): Also known as the 1-hop interference model, where two links can be active simultaneously only if they are at least $K$ hops apart.

2. Complete Interference Model: All nodes interfere with each other, allowing only one transmission in the network per time slot.

Figures 8 and 9 in [[1]](#1) present the Age-of-Information performance under these models using the Age-Difference policy. We have reproduced these results using our simulator and have also compared them with our proposed Index policy $\pi_{\text{Index}}$. 

---

### Simulation Results

#### 1-Hop Interference Model

![Figure A](https://github.com/nibin-raj/Multihop-AoI-IndexPolicy-Simulation/tree/main/figures/AoI_Line_Khop1_p1.png)

*Reproduced Age-Difference results for a single-source scenario under the 1-hop interference constraint ($`K = 1`$), replicating the results presented in Figure 8 of [[1]](#1). The performance of the proposed Index policy $`\pi_{\text{Index}}`$ is also shown for comparison.*


#### Complete Interference Model

![Figure B](https://github.com/nibin-raj/Multihop-AoI-IndexPolicy-Simulation/tree/main/figures/AoI_Line_ComINT_p1.png)

*Reproduced Age-Difference results for a single-source scenario under the Complete Interference model, replicating the results presented in Figure 9 of [[1]](#1). The results also include the proposed Index policy $`\pi_{\text{Index}}`$ for comparison.*

---

## References
<a id="1">[1]</a> 
V. Tripathi et.al., "Information Freshness in Multihop Wireless Networks,", in IEEE/ACM Transactions on Networking, vol. 31, no. 2, pp. 784-799, April 2023, 
 doi: 10.1109/TNET.2022.3201751.

## Acknowledgment

This work is supported by the IIT-Palakkad IPTIF grant IPTIF/HRD/DF/022.

