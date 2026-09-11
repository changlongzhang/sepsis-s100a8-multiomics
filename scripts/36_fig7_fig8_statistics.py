"""Reproduce final Fig 7/8 statistics with biological-repeat blocking."""
from pathlib import Path

import pandas as pd
import statsmodels.api as sm
from patsy import build_design_matrices
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "in_vitro_source_data" / "Fig7_Fig8_quantitative_source_data.csv"
OUT = ROOT / "figure_source_data"


def model_contrast(fit, row_a, row_b):
    design = fit.model.data.design_info
    matrix = build_design_matrices([design], pd.DataFrame([row_a, row_b]))[0]
    test = fit.t_test(matrix[0] - matrix[1])
    return float(test.effect.ravel()[0]), float(test.pvalue)


data = pd.read_csv(SOURCE)
anova_rows, contrast_rows = [], []

for panel in list("BCDEGHI"):
    z = data.loc[data.panel.eq(panel)].dropna(subset=["value", "group"]).copy()
    fit = ols("value ~ C(group)", data=z).fit()
    table = sm.stats.anova_lm(fit, typ=2)
    term = table.loc["C(group)"]
    figure_panel = "Fig 7A" if panel in list("BCDE") else "Fig 7C"
    anova_rows.append({"source_panel": panel, "figure_panel": figure_panel,
        "endpoint": z.endpoint.iloc[0], "model": "One-way ANOVA",
        "term": "Treatment group", "df_numerator": int(term["df"]),
        "df_denominator": int(fit.df_resid), "F_statistic": float(term["F"]),
        "p_value": float(term["PR(>F)"])})
    tukey = pairwise_tukeyhsd(z.value, z.group)
    for row in tukey.summary().data[1:]:
        p_adjusted = float(row[3])
        contrast_rows.append({"source_panel": panel, "figure_panel": figure_panel,
            "endpoint": z.endpoint.iloc[0], "contrast": f"{row[0]} vs {row[1]}",
            "estimate": float(row[2]), "p_value_raw": "", "p_adjusted": p_adjusted,
            "method": "Tukey HSD", "significant": p_adjusted < 0.05})

for panel, figure_panel in [("N", "Fig 8B"), ("O", "Fig 8C")]:
    z = data.loc[data.panel.eq(panel)].dropna(subset=["value", "replicate", "siRNA", "LPS", "TPL"]).copy()
    z["si"] = (z.siRNA == "siS100A8").astype(int)
    z["lps"] = z.LPS.astype(int)
    z["tpl"] = z.TPL.astype(int)
    fit = ols("value ~ C(replicate) + C(si)*C(lps)*C(tpl)", data=z).fit()
    table = sm.stats.anova_lm(fit, typ=2)
    labels = {"C(replicate)": "Biological repeat", "C(si)": "siRNA", "C(lps)": "LPS",
        "C(tpl)": "TPL", "C(si):C(lps)": "siRNA × LPS", "C(si):C(tpl)": "siRNA × TPL",
        "C(lps):C(tpl)": "LPS × TPL", "C(si):C(lps):C(tpl)": "siRNA × LPS × TPL"}
    for term_name, term_label in labels.items():
        term = table.loc[term_name]
        anova_rows.append({"source_panel": panel, "figure_panel": figure_panel,
            "endpoint": z.endpoint.iloc[0], "model": "Repeat-blocked 2x2x2 ANOVA",
            "term": term_label, "df_numerator": int(term["df"]),
            "df_denominator": int(fit.df_resid), "F_statistic": float(term["F"]),
            "p_value": float(term["PR(>F)"])})
    specs = [("siNC: LPS vs Control", (0,1,0), (0,0,0)),
        ("siNC: LPS+TPL vs LPS", (0,1,1), (0,1,0)),
        ("siS100A8: LPS vs Control", (1,1,0), (1,0,0)),
        ("siS100A8: LPS+TPL vs LPS", (1,1,1), (1,1,0)),
        ("LPS: siS100A8 vs siNC", (1,1,0), (0,1,0))]
    raw = []
    for name, a, b in specs:
        common = {"replicate": z.replicate.iloc[0]}
        raw.append((name, *model_contrast(fit,
            {**common, "si": a[0], "lps": a[1], "tpl": a[2]},
            {**common, "si": b[0], "lps": b[1], "tpl": b[2]})))
    adjusted = multipletests([item[2] for item in raw], method="holm")[1]
    for (name, estimate, p_value), p_adjusted in zip(raw, adjusted):
        contrast_rows.append({"source_panel": panel, "figure_panel": figure_panel,
            "endpoint": z.endpoint.iloc[0], "contrast": name, "estimate": estimate,
            "p_value_raw": p_value, "p_adjusted": p_adjusted,
            "method": "Repeat-blocked model contrast; Holm",
            "significant": bool(p_adjusted < 0.05)})

for panel, figure_panel in [("Q", "Fig 8E"), ("R", "Fig 8F"), ("S", "Fig 8G")]:
    z = data.loc[data.panel.eq(panel)].dropna(subset=["value", "replicate", "siRNA", "LPS"]).copy()
    z["si"] = (z.siRNA == "siS100A8").astype(int)
    z["lps"] = z.LPS.astype(int)
    fit = ols("value ~ C(replicate) + C(si)*C(lps)", data=z).fit()
    table = sm.stats.anova_lm(fit, typ=2)
    labels = {"C(replicate)": "Biological repeat", "C(si)": "siRNA", "C(lps)": "LPS",
        "C(si):C(lps)": "siRNA × LPS"}
    for term_name, term_label in labels.items():
        term = table.loc[term_name]
        anova_rows.append({"source_panel": panel, "figure_panel": figure_panel,
            "endpoint": z.endpoint.iloc[0], "model": "Repeat-blocked 2x2 ANOVA",
            "term": term_label, "df_numerator": int(term["df"]),
            "df_denominator": int(fit.df_resid), "F_statistic": float(term["F"]),
            "p_value": float(term["PR(>F)"])})
    specs = [("siNC: LPS vs Control", (0,1), (0,0)),
        ("LPS: siS100A8 vs siNC", (1,1), (0,1))]
    raw = []
    for name, a, b in specs:
        common = {"replicate": z.replicate.iloc[0]}
        raw.append((name, *model_contrast(fit,
            {**common, "si": a[0], "lps": a[1]},
            {**common, "si": b[0], "lps": b[1]})))
    adjusted = multipletests([item[2] for item in raw], method="holm")[1]
    for (name, estimate, p_value), p_adjusted in zip(raw, adjusted):
        contrast_rows.append({"source_panel": panel, "figure_panel": figure_panel,
            "endpoint": z.endpoint.iloc[0], "contrast": name, "estimate": estimate,
            "p_value_raw": p_value, "p_adjusted": p_adjusted,
            "method": "Repeat-blocked model contrast; Holm",
            "significant": bool(p_adjusted < 0.05)})

OUT.mkdir(parents=True, exist_ok=True)
pd.DataFrame(anova_rows).to_csv(OUT / "Fig7_Fig8_factorial_ANOVA.csv", index=False)
pd.DataFrame(contrast_rows).to_csv(OUT / "Fig7_Fig8_planned_contrasts.csv", index=False)
print(f"Wrote {len(anova_rows)} ANOVA rows and {len(contrast_rows)} contrast rows.")
