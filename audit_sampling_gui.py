import streamlit as st
import pandas as pd
from math import ceil, log, sqrt
from scipy.stats import norm, binom

st.set_page_config(page_title="Audit Sampling Calculator", layout="centered")
st.title("📊 Audit Sampling & Evaluation Toolkit")

with st.sidebar:
    st.header("Select Calculator")
    calc_type = st.radio("Choose Function:", [
        "Discovery Sampling",
        "Sample Size (Cochran)",
        "Stratified Sampling",
        "Evaluate Error Rate",
        "Max Allowed Errors (Binomial)"
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
    N = st.number_input("Population Size (optional)", min_value=0, value=0)
    confidence = st.slider("Confidence Level", 0.8, 0.999, 0.95)

    if st.button("Estimate Error Rate CI"):
        p_hat = x / n
        z = norm.ppf(1 - (1 - confidence) / 2)
        se = sqrt(p_hat * (1 - p_hat) / n)
        if N > n:
            se *= sqrt((N - n) / (N - 1))
        margin = z * se
        ci_lower = max(0, round(p_hat - margin, 4))
        ci_upper = min(1, round(p_hat + margin, 4))
        st.success(f"Observed Error Rate: {round(p_hat, 4)}")
        st.info(f"{int(confidence*100)}% CI: ({ci_lower}, {ci_upper}) ± {round(margin, 4)}")

elif calc_type == "Max Allowed Errors (Binomial)":
    st.subheader("Max Errors Allowed for Given Sample Size")
    n = st.number_input("Sample Size (n)", min_value=1)
    p0 = st.number_input("Expected Max Error Rate (p0)", 0.001, 0.99, 0.1)
    confidence = st.slider("Confidence Level", 0.8, 0.999, 0.95)
    alpha = 1 - confidence

    if st.button("Compute Max Allowed Errors"):
        max_errors = None
        for x in range(n + 1):
            if binom.cdf(x, n, p0) >= 1 - alpha:
                max_errors = x
                break
        if max_errors is not None:
            st.success(f"You may observe up to **{max_errors}** errors and still conclude error rate ≤ {p0:.2f} with {int(confidence*100)}% confidence.")
        else:
            st.warning("No valid error threshold found. Try increasing sample size or adjusting error rate.")
