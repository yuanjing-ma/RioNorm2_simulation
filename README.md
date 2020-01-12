# RioNorm2_simulation
simulation/evaluation codes for RioNorm2 approach (both Multinomial and Dirichlet-Multinomial-based model)

Large scale of simulation studies are conducted to compare the performance of RioNorm2 to that of DESeq, DESeq2, metagenomeSeq, RAIDA, Omnibus and ANCOM. Specifically, we focus on evaluating the impact of sample size, library size, and effect size on the performance of different methods. 

We explore two simulation settings which are based on different distributions:

-The first one is adopted from the simulation setting B in McMurdie and Holmes (2014), which is based on Multinomial distribution. The simulation and evaluation codes for RioNorm2 and all other compared methods are saved under Multinomial folder.

-The second simulation setting is based on the Dirichlet- Multinomial distribution. The simulation and evaluation codes for RioNorm2 and all other compared methods are saved under Dirichlet-Multinomial folder.

References:
McMurdie, P.J., and Holmes, S., (2014). Waste not, want not: why rarefying microbiome data is inadmissible. PLoS computational biology, 10(4), p.e1003531.
