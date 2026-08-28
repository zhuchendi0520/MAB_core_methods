#!/usr/bin/env python3
"""
Estimate pre/post pN/pS switch probabilities from synonym/nonsynonym counts.

Input is the wide pNS table where columns D-I are pre-expansion and
columns J-O are post-expansion:
  D-I: pNS, total_events, observed_syn, observed_nonsyn, expected_syn, expected_nonsyn
  J-O: pNS.1, total_events.1, observed_syn.1, observed_nonsyn.1, expected_syn.1, expected_nonsyn.1

For each gene, the observed nonsynonymous fraction is modeled as:
  p ~ Beta(observed_nonsyn + alpha, observed_syn + beta)

Each sampled p is converted to omega / pN-pS using the expected nonsyn:syn
opportunity ratio:
  omega = (p / (1 - p)) / (expected_nonsyn / expected_syn)

The switch probability is the posterior probability that omega crosses 1
between pre and post.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


PRE_COLS = {
    "pNS_pre": "pNS",
    "total_events_pre": "total_events",
    "observed_syn_pre": "observed_syn",
    "observed_nonsyn_pre": "observed_nonsyn",
    "expected_syn_pre": "expected_syn",
    "expected_nonsyn_pre": "expected_nonsyn",
}

POST_COLS = {
    "pNS_post": "pNS.1",
    "total_events_post": "total_events.1",
    "observed_syn_post": "observed_syn.1",
    "observed_nonsyn_post": "observed_nonsyn.1",
    "expected_syn_post": "expected_syn.1",
    "expected_nonsyn_post": "expected_nonsyn.1",
}


def numeric_series(df: pd.DataFrame, col: str) -> pd.Series:
    return pd.to_numeric(df[col].astype(object).replace("-", np.nan), errors="coerce")


def omega_from_beta_draws(
    rng: np.random.Generator,
    observed_syn: float,
    observed_nonsyn: float,
    expected_syn: float,
    expected_nonsyn: float,
    draws: int,
    prior_alpha: float,
    prior_beta: float,
) -> np.ndarray | None:
    if not np.isfinite([observed_syn, observed_nonsyn, expected_syn, expected_nonsyn]).all():
        return None
    if expected_syn <= 0 or expected_nonsyn <= 0:
        return None

    p = rng.beta(observed_nonsyn + prior_alpha, observed_syn + prior_beta, size=draws)
    p = np.clip(p, 1e-12, 1 - 1e-12)
    opportunity_odds = expected_nonsyn / expected_syn
    return (p / (1 - p)) / opportunity_odds


def point_switch_type(pre: float, post: float) -> str:
    if not np.isfinite(pre) or not np.isfinite(post):
        return "insufficient"
    if pre < 1 and post > 1:
        return "purifying_to_positive"
    if pre > 1 and post < 1:
        return "positive_to_purifying"
    if pre < 1 and post < 1:
        return "both_below_1"
    if pre > 1 and post > 1:
        return "both_above_1"
    return "touches_1"


def confidence_label(prob: float, direction: str) -> str:
    if not np.isfinite(prob):
        return "insufficient"
    if prob >= 0.95:
        return f"high_confidence_{direction}"
    if prob >= 0.90:
        return f"strong_candidate_{direction}"
    if prob >= 0.80:
        return f"moderate_candidate_{direction}"
    return "no_high_confidence_switch"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate Bayesian posterior probabilities for pre/post pN/pS switches."
    )
    parser.add_argument("input_xlsx", type=Path)
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--draws", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--prior-alpha", type=float, default=1.0)
    parser.add_argument("--prior-beta", type=float, default=1.0)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    df = pd.read_excel(args.input_xlsx)

    out = df[["DCC", "locus", "gene", "product"]].copy()
    for new_col, old_col in PRE_COLS.items():
        out[new_col] = numeric_series(df, old_col)
    for new_col, old_col in POST_COLS.items():
        out[new_col] = numeric_series(df, old_col)

    records: list[dict[str, object]] = []
    for row in out.itertuples(index=False):
        pre = omega_from_beta_draws(
            rng,
            row.observed_syn_pre,
            row.observed_nonsyn_pre,
            row.expected_syn_pre,
            row.expected_nonsyn_pre,
            args.draws,
            args.prior_alpha,
            args.prior_beta,
        )
        post = omega_from_beta_draws(
            rng,
            row.observed_syn_post,
            row.observed_nonsyn_post,
            row.expected_syn_post,
            row.expected_nonsyn_post,
            args.draws,
            args.prior_alpha,
            args.prior_beta,
        )

        base = {
            "omega_pre_median": np.nan,
            "omega_pre_ci2.5": np.nan,
            "omega_pre_ci97.5": np.nan,
            "omega_post_median": np.nan,
            "omega_post_ci2.5": np.nan,
            "omega_post_ci97.5": np.nan,
            "prob_pre_gt_1": np.nan,
            "prob_post_gt_1": np.nan,
            "prob_post_gt_pre": np.nan,
            "prob_pre_gt_post": np.nan,
            "prob_purifying_to_positive": np.nan,
            "prob_positive_to_purifying": np.nan,
            "posterior_switch_probability": np.nan,
            "posterior_switch_direction": "insufficient",
            "posterior_switch_call": "insufficient",
        }

        if pre is not None and post is not None:
            prob_p2p = float(np.mean((pre < 1) & (post > 1)))
            prob_pos2pur = float(np.mean((pre > 1) & (post < 1)))
            if prob_p2p >= prob_pos2pur:
                switch_prob = prob_p2p
                switch_direction = "purifying_to_positive"
            else:
                switch_prob = prob_pos2pur
                switch_direction = "positive_to_purifying"

            base.update(
                {
                    "omega_pre_median": float(np.median(pre)),
                    "omega_pre_ci2.5": float(np.quantile(pre, 0.025)),
                    "omega_pre_ci97.5": float(np.quantile(pre, 0.975)),
                    "omega_post_median": float(np.median(post)),
                    "omega_post_ci2.5": float(np.quantile(post, 0.025)),
                    "omega_post_ci97.5": float(np.quantile(post, 0.975)),
                    "prob_pre_gt_1": float(np.mean(pre > 1)),
                    "prob_post_gt_1": float(np.mean(post > 1)),
                    "prob_post_gt_pre": float(np.mean(post > pre)),
                    "prob_pre_gt_post": float(np.mean(pre > post)),
                    "prob_purifying_to_positive": prob_p2p,
                    "prob_positive_to_purifying": prob_pos2pur,
                    "posterior_switch_probability": switch_prob,
                    "posterior_switch_direction": switch_direction,
                    "posterior_switch_call": confidence_label(switch_prob, switch_direction),
                }
            )

        records.append(base)

    stats = pd.DataFrame.from_records(records)
    result = pd.concat([out, stats], axis=1)

    result["point_switch_type"] = [
        point_switch_type(pre, post)
        for pre, post in zip(result["pNS_pre"], result["pNS_post"])
    ]
    result["delta_log2_pNS_post_vs_pre"] = np.log2(
        (result["pNS_post"].clip(lower=1e-4)) / (result["pNS_pre"].clip(lower=1e-4))
    )
    result["min_total_events_pre_post"] = result[
        ["total_events_pre", "total_events_post"]
    ].min(axis=1)
    result["total_events_pre_post"] = result[
        ["total_events_pre", "total_events_post"]
    ].sum(axis=1)

    sort_cols = ["posterior_switch_probability", "total_events_pre_post"]
    result_sorted = result.sort_values(sort_cols, ascending=[False, False])

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_prefix.parent / f"{args.output_prefix.name}.csv"
    xlsx_path = args.output_prefix.parent / f"{args.output_prefix.name}.xlsx"
    summary_path = args.output_prefix.parent / f"{args.output_prefix.name}_summary.csv"
    candidates_path = (
        args.output_prefix.parent / f"{args.output_prefix.name}_posterior_switch_ge0.80.csv"
    )

    result_sorted.to_csv(csv_path, index=False)
    result_sorted.to_excel(xlsx_path, index=False)

    summary = (
        result_sorted.groupby(["posterior_switch_call", "posterior_switch_direction"], dropna=False)
        .size()
        .reset_index(name="gene_count")
        .sort_values("gene_count", ascending=False)
    )
    summary.to_csv(summary_path, index=False)

    candidates = result_sorted[result_sorted["posterior_switch_probability"] >= 0.80]
    candidates.to_csv(candidates_path, index=False)

    print(f"Input rows: {len(df)}")
    print(f"Output: {csv_path}")
    print(f"Output: {xlsx_path}")
    print(f"Summary: {summary_path}")
    print(f"Posterior switch >= 0.80: {len(candidates)} genes")
    print(candidates["posterior_switch_call"].value_counts(dropna=False).to_string())


if __name__ == "__main__":
    main()
