# Audit Sampling & Evaluation Toolkit

A Streamlit-based calculator for statistical audit sampling. Provides five calculators covering sample size determination, stratified allocation, error rate estimation, and binomial tolerance testing.

## Running the App

```bash
pip install -r requirements.txt
streamlit run audit_sampling_gui.py
```

---

## Calculators & Formulas

### 1. Discovery Sampling

**Purpose:** Determine the minimum sample size needed to detect at least one error with a given confidence, assuming no errors will be found (i.e., a "clean" population test).

**Formula:**

$$n = \left\lceil \frac{\ln(\alpha)}{\ln(1 - p_0)} \right\rceil$$

Where:
- `α = 1 − confidence` — the acceptable probability of missing an error
- `p₀` — the maximum tolerable error rate you want to be able to detect
- `ln` — natural logarithm

**Derivation:** The probability of drawing zero errors in `n` samples from a population with error rate `p₀` is `(1 − p₀)ⁿ`. Setting this equal to `α` (the acceptable miss probability) and solving for `n` gives the formula above.

**Inputs:** Expected max error rate (`p₀`), confidence level  
**Output:** Minimum sample size where observing 0 errors is sufficient to conclude the population error rate is ≤ `p₀`

---

### 2. Sample Size — Cochran's Formula

**Purpose:** Calculate the sample size needed to estimate a proportion (error rate) within a specified margin of error.

**Formula (infinite population):**

$$n_0 = \frac{z^2 \cdot p(1-p)}{e^2}$$

**Finite Population Correction (FPC):**

$$n = \left\lceil \frac{n_0}{1 + \dfrac{n_0 - 1}{N}} \right\rceil$$

Where:
- `z` — z-score for the chosen confidence level (e.g., 1.96 for 95%)
- `p` — estimated error rate (use 0.5 if unknown, maximises sample size)
- `e` — desired margin of error (half-width of the confidence interval)
- `N` — population size (omit or set to 0 to skip FPC)

**Inputs:** Estimated error rate (`p`), margin of error (`e`), confidence level, optional population size (`N`)  
**Output:** Required sample size

---

### 3. Stratified Sampling

**Purpose:** Allocate a total sample across strata (sub-populations) to achieve the overall desired precision.

**Step 1 — Total sample size** uses Cochran's formula with FPC applied to the full population `N_total`:

$$n_{total} = \left\lceil \frac{n_0}{1 + \dfrac{n_0 - 1}{N_{total}}} \right\rceil$$

**Step 2 — Allocation per stratum:**

#### Proportional Allocation
Each stratum receives a share proportional to its size:

$$n_h = \frac{N_h}{N_{total}} \cdot n_{total}$$

Where `N_h` is the size of stratum `h`.

#### Neyman (Optimal) Allocation
Each stratum is weighted by its size *and* variability — strata with higher spread receive more samples:

$$n_h = \frac{N_h \cdot \sigma_h}{\sum_{h} N_h \cdot \sigma_h} \cdot n_{total}$$

Where `σ_h` is the standard deviation of stratum `h` (requires a `std_dev` column in the uploaded CSV).

**CSV format:** `stratum, N [, std_dev]`  
**Inputs:** CSV file, estimated proportion (`p`), margin of error (`e`), confidence level, allocation method  
**Output:** Per-stratum sample sizes (`n_h`) and total sample size

---

### 4. Evaluate Error Rate (Confidence Interval)

**Purpose:** Given errors found in a sample, estimate the true population error rate with a confidence interval.

**Observed error rate:**

$$\hat{p} = \frac{x}{n}$$

**Standard error:**

$$SE = \sqrt{\frac{\hat{p}(1-\hat{p})}{n}}$$

**With Finite Population Correction:**

$$SE_{FPC} = SE \cdot \sqrt{\frac{N - n}{N - 1}}$$

(Applied only when population size `N` is provided and `N > n`)

**Confidence interval:**

$$CI = \hat{p} \pm z \cdot SE$$

Where:
- `x` — number of errors found
- `n` — sample size
- `N` — population size (optional)
- `z` — z-score for the chosen confidence level

**Inputs:** Sample size (`n`), errors found (`x`), optional population size (`N`), confidence level  
**Output:** Point estimate `p̂`, confidence interval lower/upper bounds, and margin of error

---

### 5. Max Allowed Errors (Binomial)

**Purpose:** Given a fixed sample size, determine how many errors can be observed while still concluding (with the desired confidence) that the true error rate does not exceed `p₀`.

**Method:** Finds the smallest integer `x` satisfying:

$$P(X \leq x \mid n, p_0) \geq 1 - \alpha$$

where `X ~ Binomial(n, p₀)` and `α = 1 − confidence`.

This is equivalent to finding the `(1 − α)` quantile of the binomial distribution. If the number of observed errors is ≤ `x`, the sample is consistent with a population error rate ≤ `p₀` at the specified confidence level.

**Inputs:** Sample size (`n`), maximum tolerable error rate (`p₀`), confidence level  
**Output:** Maximum number of errors that can be observed while still passing the test

---

## Method Selection Guide

| Situation | Recommended Calculator |
|---|---|
| Testing for any error in a clean population | Discovery Sampling |
| Estimating how many errors exist | Sample Size (Cochran) |
| Population divided into segments of different risk | Stratified Sampling |
| You have results and want to bound the error rate | Evaluate Error Rate |
| You have a fixed sample and need a pass/fail threshold | Max Allowed Errors (Binomial) |

---

## Dependencies

- `streamlit`
- `pandas`
- `scipy`
