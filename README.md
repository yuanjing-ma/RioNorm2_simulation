# RioNorm2_simulation
simulation/evaluation codes for RioNorm2 approach (both Multinomial and Dirichlet-Multinomial-based model)

Large scale of simulation studies are conducted to compare the performance of RioNorm2 to that of DESeq, DESeq2, metagenomeSeq, RAIDA, Omnibus and ANCOM. Specifically, we focus on evaluating the impact of sample size, library size, and effect size on the performance of different methods. 

We explore two simulation settings which are based on different distributions:

-The first one is adopted from the simulation setting B in McMurdie and Holmes (2014), which is based on Multinomial distribution. The simulation and evaluation codes for RioNorm2 and all other compared methods are saved under Multinomial folder.

  - Main files:

    - Simulation.R: generate multinomial distribution-based data
    - simulation_performance_comparison.R: contains (1) codes for all compared methods, e.g. RioNorm2, RAIDA, DESeq, DESeq2, metagenomeSeq, Omnibus (2) organize results and generate figures

  - Other files:
    - ANCOM_RioNorm2/ANCOM_RioNorm2_subset_multi.R: compare RioNorm2 and ANCOM on subset of simulated data
    - Other_DA_methods/...: combine RioNorm2 common divisors with other down-stream testing methods, such as zero-inflation possion, zero-inflated negative binomial, t-test etc. 
    - RAIDA_RioNorm2_common_divisor/RAIDA_RioNorm2_common_divisor.R: compare common divisors of RAIDA and RioNorm2
    - diff_percentage_TPs/...: generate simulation data with different percentage of DA-OTUs and compare different approaches
    - h_robustness_analysis/Simulation_h_robustness_RioNorm2.R: robustness analysis of h value in RioNorm2 algorithm

-The second simulation setting is based on the Dirichlet- Multinomial distribution. The simulation and evaluation codes for RioNorm2 and all other compared methods are saved under Dirichlet-Multinomial folder.

References:
McMurdie, P.J., and Holmes, S., (2014). Waste not, want not: why rarefying microbiome data is inadmissible. PLoS computational biology, 10(4), p.e1003531.
