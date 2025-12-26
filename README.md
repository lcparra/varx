# VARX Granger Analysis

VARX model estimation and statistical analysis based on Granger formalism. 

The methods are described in: 

Parra, L. C., Silvan, A., Nentwich, M., Madsen, J., Parra, V. E., & Babadi, B. (2025). VARX Granger analysis: Models for neuroscience, physiology, sociology and econometrics. PLOS ONE, 20(1), e0313875. https://doi.org/10.1371/journal.pone.0313875

The code is available for matlab, python and R.

This repository contains example code from multiple projects. To generate the figures of the paper follow these instructions:

- Figure 2: run code cell in varx_demo.m
- Figure 3: run code cell in varx_demo.m
- Figure 4: run examples/simulate_gain_adaptation.m and follow instructions on where to get the data. 
- Figure 5: data and code not yet available.
- Figure 6: run examples/varx_union_history.m
- Figure 7: run examples/varx_economy_example.m
- Figure 8: run code cell in varx_demo.m
- Figure 9: run code cell in varx_demo.m
- Figure 10: run examples/granger_simulation.m

Under examples/ are several applications of the methods to neurals signals published here:  

Maximilian Nentwich, Marcin Leszczynski, Charles E Schroeder, Stephan Bickel, Lucas C. Parra, "Intrinsic dynamic shapes responses to external stimulation in the human brain", eLife, 104996, February 27, 2025, https://elifesciences.org/articles/104996

- Figure 2 part of this figure, run examples/innovation_rest_vs_movie_hfa.m	
- Figure 7 run phython and malap bode under  examples/neurolib_brainmodel/
- Appendix 1 run examples/simulate_gain_adaptation.m
- Figure 3b to compute H use varx_trf.m 

Addions of and explained in

Jens Madsen, Aimar Silvan, Behtash Babadi, Lucas C Parra "Brain-body dynamics is asymmetric and stable across cognitive states", bioRxiv, December 2025 , https://www.biorxiv.org/content/10.64898/2025.12.22.696055v1

- Figure 2c-d to compute inv(1-A) use varx_trf.m 
- Figure S6:  Akike Information Criterion to find optimal model parameters na, nb, part of varx.m
- Fiugre S8 a-c:  similarity_shuffle_test.m 
