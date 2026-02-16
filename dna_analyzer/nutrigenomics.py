"""Nutrigenomics analysis — nutrient metabolism, dietary sensitivities."""

from .snp_database import NUTRIGENOMICS_SNPS


CATEGORY_LABELS = {
    "caffeine_metabolism": "Caffeine Metabolism",
    "lactose": "Lactose Tolerance",
    "alcohol": "Alcohol Metabolism",
    "vitamin_d": "Vitamin D",
    "folate": "Folate / MTHFR",
    "vitamin_b12": "Vitamin B12",
    "vitamin_a": "Vitamin A / Beta-carotene",
    "fat_metabolism": "Omega-3 / Fat Metabolism",
    "salt_sensitivity": "Salt Sensitivity",
    "choline": "Choline",
}


def _assess_impact(genotype, entry):
    """Assess the nutritional impact of a genotype."""
    effect = entry.get("effect", {})
    effect_text = effect.get(genotype, "")
    recommendation = entry.get("recommendation", {}).get(genotype, "")

    text_lower = effect_text.lower()

    if any(w in text_lower for w in ["reduced", "lower", "impaired", "poor",
                                      "intolerant", "slow metabolizer", "non-secretor",
                                      "non-functional", "salt-sensitive"]):
        if "slightly" in text_lower or "mildly" in text_lower:
            impact = "moderate"
        else:
            impact = "significant"
    elif any(w in text_lower for w in ["intermediate", "moderate", "slightly"]):
        impact = "moderate"
    elif any(w in text_lower for w in ["normal", "tolerant", "fast", "good",
                                        "efficient", "standard"]):
        impact = "minimal"
    else:
        impact = "moderate"

    return impact, effect_text, recommendation


def analyze_nutrigenomics(user_snps):
    """Analyze nutrigenomics SNPs.

    Args:
        user_snps: dict of rsid -> {'genotype': str, ...}

    Returns:
        dict with categorized nutrigenomic findings
    """
    results = {cat: [] for cat in CATEGORY_LABELS}
    dietary_recommendations = []
    supplement_recommendations = []

    for rsid, entry in NUTRIGENOMICS_SNPS.items():
        if rsid in user_snps:
            genotype = user_snps[rsid]["genotype"]
            impact, effect_text, recommendation = _assess_impact(genotype, entry)

            finding = {
                "rsid": rsid,
                "gene": entry["gene"],
                "genotype": genotype,
                "condition": entry["condition"],
                "impact": impact,
                "effect": effect_text,
                "recommendation": recommendation,
                "source": entry.get("source", ""),
            }

            category = entry["category"]
            if category in results:
                results[category].append(finding)

            # Collect actionable recommendations
            if recommendation and impact in ("significant", "moderate"):
                if "supplement" in recommendation.lower() or "mcg" in recommendation.lower() \
                   or "iu" in recommendation.lower() or "mg" in recommendation.lower():
                    supplement_recommendations.append({
                        "gene": entry["gene"],
                        "recommendation": recommendation,
                        "priority": "high" if impact == "significant" else "moderate",
                    })
                else:
                    dietary_recommendations.append({
                        "gene": entry["gene"],
                        "recommendation": recommendation,
                        "priority": "high" if impact == "significant" else "moderate",
                    })

    # Sort by impact
    impact_order = {"significant": 0, "moderate": 1, "minimal": 2}
    for cat in results:
        results[cat].sort(key=lambda x: impact_order.get(x["impact"], 2))

    return {
        "findings": results,
        "dietary_recommendations": dietary_recommendations,
        "supplement_recommendations": supplement_recommendations,
        "category_labels": CATEGORY_LABELS,
    }


def get_mthfr_status(user_snps):
    """Specifically assess MTHFR compound heterozygosity."""
    c677t = user_snps.get("rs1801133", {}).get("genotype", "")
    a1298c = user_snps.get("rs1801131", {}).get("genotype", "")

    if not c677t and not a1298c:
        return None

    status = {
        "c677t": {"genotype": c677t, "rsid": "rs1801133"},
        "a1298c": {"genotype": a1298c, "rsid": "rs1801131"},
    }

    # Assess compound heterozygosity
    c677t_het = c677t in ("AG", "GA", "CT", "TC")
    a1298c_het = a1298c in ("GT", "TG", "AC", "CA")
    c677t_hom = c677t in ("AA", "TT")
    a1298c_hom = a1298c in ("GG", "CC")

    if c677t_hom:
        status["severity"] = "significant"
        status["interpretation"] = (
            "Homozygous MTHFR C677T — ~70% reduced enzyme activity. "
            "Methylfolate supplementation strongly recommended."
        )
    elif c677t_het and a1298c_het:
        status["severity"] = "moderate"
        status["interpretation"] = (
            "Compound heterozygous MTHFR (C677T + A1298C) — "
            "~50% reduced activity. Methylfolate recommended."
        )
    elif c677t_het:
        status["severity"] = "mild"
        status["interpretation"] = (
            "Heterozygous MTHFR C677T — ~35% reduced activity. "
            "Consider methylfolate."
        )
    elif a1298c_hom:
        status["severity"] = "mild"
        status["interpretation"] = (
            "Homozygous MTHFR A1298C — mild reduction. "
            "May affect BH4 (neurotransmitter cofactor)."
        )
    elif a1298c_het:
        status["severity"] = "minimal"
        status["interpretation"] = (
            "Heterozygous MTHFR A1298C — minimal clinical impact alone."
        )
    else:
        status["severity"] = "none"
        status["interpretation"] = "Normal MTHFR function."

    return status
