import streamlit as st
import pandas as pd
from math import ceil, log
from scipy.stats import norm, binom, beta

st.set_page_config(page_title="Audit Sampling Calculator", layout="centered")
st.title("📊 Audit Sampling & Evaluation Toolkit")

with st.sidebar:
    st.header("Select Calculator")
    calc_type = st.radio("Choose Function:", [
        "Discovery Sampling",
        "Sample Size (Cochran)",
        "Stratified Sampling",
        "Evaluate Error Rate",
        "Max Allowed Errors (Binomial)",
        "Methodology Statement"
    ])

if calc_type == "Discovery Sampling":
    st.subheader("Discovery Sampling")
    p0 = st.number_input("Expected Max Error Rate (e.g., 0.1 = 10%)", min_value=0.001, max_value=0.99, value=0.1)
    confidence = st.slider("Confidence Level", 0.8, 0.999, 0.95)

    if st.button("Calculate Sample Size"):
        alpha = 1 - confidence
        n = ceil(log(alpha) / log(1 - p0))
        st.success(f"Required sample size (0 errors allowed): **{n}**")

elif calc_type == "Sample Size (Cochran)":
    st.subheader("Sample Size via Cochran's Formula")
    p = st.number_input("Estimated Error Rate (p)", 0.01, 0.99, 0.5)
    e = st.number_input("Margin of Error (e)", 0.005, 0.2, 0.05)
    confidence = st.slider("Confidence Level", 0.8, 0.999, 0.95)
    N = st.number_input("Population Size (optional)", min_value=0, value=0)

    if st.button("Compute Cochran Sample Size"):
        z = norm.ppf(1 - (1 - confidence) / 2)
        n0 = (z**2 * p * (1 - p)) / e**2
        if N > 0:
            n = ceil(n0 / (1 + (n0 - 1) / N))
        else:
            n = ceil(n0)
        st.success(f"Required sample size: **{n}**")

elif calc_type == "Stratified Sampling":
    st.subheader("Stratified Sampling Calculator")
    uploaded = st.file_uploader("Upload CSV with columns: stratum, N [, std_dev]", type="csv")
    p = st.number_input("Estimated Proportion (p)", 0.01, 0.99, 0.5)
    e = st.number_input("Margin of Error (e)", 0.005, 0.2, 0.05)
    confidence = st.slider("Confidence Level", 0.8, 0.999, 0.95)
    allocation = st.selectbox("Allocation Method", ["proportional", "neyman"])

    if uploaded:
        strata_df = pd.read_csv(uploaded)
        z = norm.ppf(1 - (1 - confidence) / 2)
        n0 = (z**2 * p * (1 - p)) / e**2
        N_total = strata_df['N'].sum()
        n_adj = n0 / (1 + (n0 - 1) / N_total)
        n_total = ceil(n_adj)

        if allocation == 'proportional':
            strata_df['n_h'] = strata_df['N'] / N_total * n_total
        else:
            if 'std_dev' not in strata_df.columns:
                st.error("std_dev column required for Neyman allocation")
            else:
                strata_df['weight'] = strata_df['N'] * strata_df['std_dev']
                total_weight = strata_df['weight'].sum()
                strata_df['n_h'] = strata_df['weight'] / total_weight * n_total

        strata_df['n_h'] = strata_df['n_h'].apply(ceil)
        st.dataframe(strata_df[['stratum', 'N', 'n_h']])
        st.success(f"Total sample size: {n_total}")

elif calc_type == "Evaluate Error Rate":
    st.subheader("Error Rate Confidence Interval")
    n = st.number_input("Sample Size (n)", min_value=1)
    x = st.number_input("Number of Errors Found (x)", min_value=0)
    confidence = st.slider("Confidence Level", 0.8, 0.999, 0.95)
    limit = st.radio("Limit", ["Two-sided (lower and upper)", "Upper limit only"])

    if st.button("Estimate Error Rate CI"):
        if x > n:
            st.error("Errors found cannot exceed the sample size.")
        else:
            # Exact (Clopper-Pearson) binomial limits: valid with zero or few errors,
            # where the normal approximation collapses or understates the upper limit.
            p_hat = x / n
            if limit.startswith("Upper"):
                ci_lower = 0.0
                ci_upper = 1.0 if x == n else beta.ppf(confidence, x + 1, n - x)
                label = f"{int(confidence*100)}% upper limit"
            else:
                a = 1 - confidence
                ci_lower = 0.0 if x == 0 else beta.ppf(a / 2, x, n - x + 1)
                ci_upper = 1.0 if x == n else beta.ppf(1 - a / 2, x + 1, n - x)
                label = f"{int(confidence*100)}% CI"
            st.success(f"Observed Error Rate: {p_hat:.2%}")
            if limit.startswith("Upper"):
                st.info(f"{label}: error rate is at most {ci_upper:.2%}")
            else:
                st.info(f"{label}: {ci_lower:.2%} to {ci_upper:.2%}")

elif calc_type == "Max Allowed Errors (Binomial)":
    st.subheader("Max Errors Allowed for Given Sample Size")
    n = st.number_input("Sample Size (n)", min_value=1)
    p0 = st.number_input("Expected Max Error Rate (p0)", 0.001, 0.99, 0.1)
    confidence = st.slider("Confidence Level", 0.8, 0.999, 0.95)
    alpha = 1 - confidence

    if st.button("Compute Max Allowed Errors"):
        # Acceptance number: the largest error count x where, if the true error rate
        # were p0, seeing x or fewer errors would happen no more than alpha of the time.
        # Finding x or fewer errors then supports "error rate <= p0" at this confidence.
        max_errors = None
        for x in range(n + 1):
            if binom.cdf(x, n, p0) <= alpha:
                max_errors = x
            else:
                break
        if max_errors is not None:
            st.success(f"You may observe up to **{max_errors}** errors and still conclude error rate ≤ {p0:.2%} with {int(confidence*100)}% confidence.")
        else:
            min_n = ceil(log(alpha) / log(1 - p0))
            st.warning(f"This sample is too small: even 0 errors would not support that conclusion. Use at least {min_n} items (the discovery sample size).")

elif calc_type == "Methodology Statement":
    st.subheader("Methodology Statement")
    st.caption("Statistical methods used by each calculator in this tool, for citation in audit workpapers.")

    methodology_text = r"""# Audit Sampling Toolkit — Methodology Statement

## Discovery Sampling

**Plain language:** Finds the smallest sample size such that, if a deviation exists in the
population at or above a stated rate, the sample has the chosen probability (confidence level)
of containing at least one instance of it. The test assumes the sample must come back with
**zero** errors in order to conclude the true rate is below the target.

**Equation:**

$$ n = \left\lceil \frac{\ln(\alpha)}{\ln(1 - p_0)} \right\rceil, \qquad \alpha = 1 - \text{confidence} $$

## Sample Size (Cochran)

**Plain language:** Estimates the sample size needed to estimate a proportion within a target
margin of error, using the normal approximation to the binomial distribution. Applies a
finite-population correction when a population size is supplied.

**Equations:**

$$ n_0 = \frac{z^2\, p(1-p)}{e^2} \qquad n = \frac{n_0}{1 + \dfrac{n_0 - 1}{N}}\ \text{(when } N \text{ is given)} $$

where $z$ is the standard normal quantile for the chosen confidence level, $p$ is the estimated
proportion, and $e$ is the target margin of error.

## Stratified Sampling

**Plain language:** Splits a total sample size across subgroups ("strata") of the population,
either in proportion to each stratum's size, or (Neyman allocation) weighted toward strata with
more internal variability, to get the most precision for a given total sample size.

**Equations:**

$$ \text{Proportional: } n_h = \frac{N_h}{N}\, n \qquad \text{Neyman: } n_h = \frac{N_h \sigma_h}{\sum_h N_h \sigma_h}\, n $$

## Evaluate Error Rate

**Plain language:** Given $x$ errors found in a sample of $n$ items, computes the range of
population error rates consistent with that result, using the **exact (Clopper-Pearson)**
binomial confidence interval rather than a normal approximation. The exact method stays valid
with zero or very few errors, where the normal approximation collapses or understates the true
upper limit. The "upper limit only" option reports just the upper bound, for sign-off decisions
that only need to establish a ceiling on the error rate.

**Equations:**

Two-sided interval (with $\alpha = 1 - \text{confidence}$):

$$ \text{Lower} = \begin{cases} 0 & x = 0 \\ B^{-1}\!\left(\tfrac{\alpha}{2};\, x,\, n-x+1\right) & x > 0 \end{cases} \qquad \text{Upper} = \begin{cases} 1 & x = n \\ B^{-1}\!\left(1-\tfrac{\alpha}{2};\, x+1,\, n-x\right) & x < n \end{cases} $$

One-sided upper limit:

$$ \text{Upper} = \begin{cases} 1 & x = n \\ B^{-1}\!\left(\text{confidence};\, x+1,\, n-x\right) & x < n \end{cases} $$

where $B^{-1}(q;\, a,\, b)$ is the inverse CDF (quantile function) of the Beta distribution.

## Max Allowed Errors (Binomial)

**Plain language:** For a given sample size and a tolerable error rate $p_0$, finds the
**acceptance number** — the largest number of errors that could be found in the sample while
still supporting the conclusion, at the chosen confidence level, that the true error rate is at
or below $p_0$. If even zero errors would not support that conclusion, the sample size is too
small for the target rate and confidence level.

**Equation:**

$$ x^* = \max\{\, x : P(X \le x \mid n, p_0) \le \alpha \,\}, \qquad X \sim \text{Binomial}(n, p_0), \quad \alpha = 1 - \text{confidence} $$
"""

    st.markdown(methodology_text)
    st.download_button(
        "Download methodology statement (Markdown)",
        data=methodology_text,
        file_name="audit_sampling_methodology.md",
        mime="text/markdown"
    )
