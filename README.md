#### HybridControl_early
A massive simulation for assessing hybrid control methods (existing and their variation) in early phase trials

This repository contains the simulation R and jags code for evaluating hybrid control methods in early phase trials. Survival outcomes are being investigated. Four Frequnetist and eight Bayesian methods were being considered. In the simulation, we assumed two historical trilas are available.

# Frequntist methods: 
1. Separated IPTW (AFT), 2. Separated full matching (AFT), 3. Joint IPTW (AFT), 4. Joint full matching (AFT)
5. Separated IPTW (Cox), 6. Separated full matching (Cox), 7. Joint IPTW (Cox), 8. Joint full matching (Cox)

# Bayesian Dynamic borrowing using:
1. Informative prior (Same), 2. Weakly-informative prior (Same), 3. Non-informative prior (Same)
4. Informative prior (Distinguish), 5. Weakly-informative prior (Distinguish), 6. Non-informative prior (Distinguish)
7. No borrpwong, 8. No borrpwong

# Sample size: 
Conconrrent trial: treated (T): 16/40, control (C): 22
Historical control 1 (HC0): 145
Historical control 2 (HC1): 467

# Treatment effects (log harvard ratio scale):
Scenario 1 (All equivalent)
beta1=c(log(1), log(1), log(1))

Scenario 2 (HC0=HC1=C<T)
beta2=c(log(1), log(1), log(2))

Scenario 3 (HC0=HC1<C<T)
beta3=c(log(0.7), log(0.7), log(2))

Scenario 4  (HC0<HC1<C<T)
beta4=c(log(0.5), log(0.7), log(2))
