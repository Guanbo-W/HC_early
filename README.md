# HybridControl_early
A massive simulation for assessing hybrid control methods (existing and their variation) in early phase trials

This repository contains the simulation R and jags code for evaluating hybrid control methods in early phase trials. Survival outcomes are being investigated. Four Frequnetist and eight Bayesian methods were being considered. In the simulation, we assumed two historical trilas are available.

Frequntist methods: 
1. Separated IPTW (AFT), 2. Separated full matching (AFT), 3. Joint IPTW (AFT), 4. Joint full matching (AFT)
5. Separated IPTW (Cox), 6. Separated full matching (Cox), 7. Joint IPTW (Cox), 8. Joint full matching (Cox)

Bayesian Dynamic borrowing using:
1. Informative prior (Same), 2. Weakly-informative prior (Same), 3. Non-informative prior (Same)
4. Informative prior (Distinguish), 5. Weakly-informative prior (Distinguish), 6. Non-informative prior (Distinguish)
7. No borrpwong, 8. No borrpwong

Sample size: 
Conconrrent trial: treated: 16/40, control: 
