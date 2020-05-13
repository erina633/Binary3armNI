# Binary3armNI: R codes for "A More Powerful Test for Three-Arm Non-Inferiority via Risk Difference: Frequentist and Bayesian Approaches"

This research is supported by PCORI contract number ME-1409-21410, PI: Samiran Ghosh

The codes are provided to generate sample size and type I error for both Frequentist and Bayesian approaches. The following notations are used in this package:

• parE, parR, and parP: Binary proportion parameters for the arm E, R, and P, respectively

• nP: The sample size in the placebo arm (P)

• nR: The sample size in the reference arm (R)

• nE: The sample size in the experimental arm (E)

• theta: effect retention parameter

• alloc: Allocation vector which can be (1:1:1), (2:2:1) or (3:2:1) for nE:nR:nP

• n: Total sample size

• n_total: Maximum number of n to get the sample size

• theta: effect retention parameter

We give brief description of the R files below:

1. freq_samplesizecode.R

This function calculates the Frequentist sample size for a given value of theta, allocation, alpha, parP, parR, and parE.

Output: Minimum sample size of the arm P satisfying power>=1-beta

2. approxbayes_samplesizecode.R

This function calculates the sample size corresponding to approximate Bayesian method for a given value of theta, allocation, alpha, parP, parR, and parE.

Output: Minimum sample size of the arm P satisfying power>=1-beta

3. exactbayes_samplesizecode.R

This function calculates the sample size corresponding to exact Bayesian method for a given value of theta, allocation, alpha, parP, parR, and parE.

Output: Minimum sample size of the arm P satisfying power>=1-beta

4. type1error_exactbayes.R

This function calculates the estimated type I error under exact Bayesian approach for different values of theta, calculated sample size, allocation, and parE.

Output: Estimated type I error for exact Bayesian approach
