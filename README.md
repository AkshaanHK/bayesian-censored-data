# Applied Bayesian Statistics — Left-Censored Environmental Data

## Overview

This project applies Bayesian inference to estimate key summary statistics
from **left-censored concentration data** — a common challenge in environmental
and occupational health studies where measurements below a detection limit are
not fully observed.

Two model settings are studied:

- **Single-element model** — one concentration variable with multiple reporting limits
- **Two-element model** — two correlated concentrations, jointly censored

Both Bayesian and frequentist methods are implemented and compared across
a simulation study (Objective A) and a real industrial dataset (Objective B).

> Note: the dataset and industry partner are confidential and not included in this repository.

---

## Research Question

How do Bayesian methods (Gibbs sampler, Metropolis-Hastings) perform relative
to frequentist approaches (substitution, ROS) for estimating summary statistics
from heavily left-censored data — and does this hold on real-world data?

---

## Methods

**Bayesian**
- Metropolis-Hastings sampler with log-normal likelihood and half-Cauchy prior on σ
- Gibbs sampler with data augmentation for censored observations
- Convergence diagnostics: Gelman-Rubin R̂, Effective Sample Size (ESS), ACF

**Frequentist**
- Substitution methods (½ DL, 1/√2 DL, beta-substitution)
- Regression on Order Statistics (ROS)

**Estimands** — for each model: geometric mean (GM), geometric standard deviation (GSD),
arithmetic mean (AM), 95th quantile (Q95)

---

## Objective A — Simulation Study

**Design:** 40 data-generating mechanisms per model (single and two-element),
varying sample size, censoring proportion, and reporting limit structure.
100 Monte Carlo replications per scenario.

Parameters:
- Sample sizes: n ∈ {25, 75, 125, 175, 225}
- Censoring proportions: p ∈ {0.2, 0.4, 0.6, 0.8}
- Single vs. multiple reporting limits

**Key findings:**

- Bayesian methods provide consistent estimates of AM, GM, and GSD across all scenarios,
  with values close to the true ones even under high censoring.
- Bayesian methods tend to underestimate the 95th quantile across all scenarios.
- Frequentist methods are less consistent, especially for AM and GM at high censoring
  proportions — they tend to underestimate AM and overestimate GM as censoring increases.
- Gibbs and Metropolis produce nearly identical results in all configurations.
- Frequentist quantile estimates are typically close to the true values but show
  higher variability.
  
📊 [View full simulation report (HTML)](https://akshaanhk.github.io/bayesian-censored-data/1_Slytherin_report_objective_A.html)

---

## Objective B — Real Dataset

**Data:** Two confidential industrial datasets provided by an industry partner.

- Single-element dataset: 11 534 observations, 77% censored, 61 unique reporting limits
  ranging from 3.0 × 10⁻¹⁰ to 8.3 × 10⁻⁴
- Two-element dataset: 4 828 observations, 96.5% censored for element 1
  and 93.5% for element 2, with 29 reporting limits per element

**Convergence diagnostics** were assessed via Gelman-Rubin R̂, ESS, and ACF plots.
**Model checking** was performed visually by overlapping original and posterior
predictive distributions.

📊 [View full simulation report (HTML)](https://akshaanhk.github.io/bayesian-censored-data/Slytherin_report_objective_B.html)

---

## Repository Structure

```
.
├── functions.R                              Core implementation
│   ├── 0. Helper functions                  Convergence diagnostics, reproducible streams
│   ├── 1. Data-generating mechanisms        Simulate single and two-element censored data
│   ├── 2. Single-element model
│   │   ├── 2.1 Metropolis-Hastings
│   │   ├── 2.2 Gibbs sampler
│   │   └── 2.3 Frequentist (substitution, ROS)
│   └── 3. Two-element model
│       ├── 3.1 Metropolis-Hastings
│       ├── 3.2 Gibbs sampler
│       └── 3.3 Frequentist
│
├── 1_Slytherin_report_objective_A.html      Simulation study report
└── Slytherin_report_objective_B.html        Real data analysis report
```

---

## Reproducibility

The simulation study (Objective A) and real-data analysis (Objective B) were
implemented in R. The core functions are fully self-contained in `functions.R`.

**Dependencies**

```r
install.packages(c("coda", "rstream", "mvtnorm", "ggplot2", "dplyr", "tidyr",
                   "scales", "gridExtra", "parallel"))
```

**Run a single Metropolis chain (example)**

```r
source("functions.R")

set.seed(1)
data <- generate_single_element(n_obs = 100, p_censoring = 0.4)

chain <- metropolis_single(
  data      = data,
  n_iter    = 10000,
  init      = c(mu = 0.5, log_sd = 0),
  mu.prior  = c(0, 10),
  sd.prior.scale = 5
)
```

> The dataset used in Objective B is confidential and not included.

---

## Key Takeaways

- Bayesian methods with proper priors handle extreme censoring (up to 96%) reliably
- Gibbs and Metropolis converge to equivalent posteriors — choice is implementation-driven
- Frequentist substitution methods degrade meaningfully above 60% censoring
- The 95th quantile is the hardest estimand to recover under heavy censoring,
  regardless of method

---

## Authors

Akshaan Murugesu
Huimin Hu
Lorena Coppola
Paul Edward Schumacher
Pooya Sabbagh Savoojbulagh
  
MSc Statistics — University of Geneva
