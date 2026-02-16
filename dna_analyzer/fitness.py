"""Fitness & recovery analysis — muscle type, injury risk, endurance, recovery."""

from .snp_database import FITNESS_SNPS


CATEGORY_LABELS = {
    "muscle_fiber": "Muscle Fiber Type",
    "endurance": "Endurance Capacity",
    "injury_risk": "Injury Susceptibility",
    "recovery": "Recovery & Oxidative Stress",
}


def _assess_fitness_trait(genotype, entry):
    """Assess the fitness implication of a genotype."""
    effect = entry.get("effect", {})
    effect_text = effect.get(genotype, "")
    recommendation = entry.get("recommendation", {}).get(genotype, "")

    text_lower = effect_text.lower()

    if any(w in text_lower for w in ["increased risk", "reduced", "poorer",
                                      "slower", "less efficient"]):
        trait_impact = "attention_needed"
    elif any(w in text_lower for w in ["power", "sprint", "endurance",
                                        "advantage", "good", "efficient"]):
        trait_impact = "favorable"
    elif "mixed" in text_lower or "intermediate" in text_lower or "versatile" in text_lower:
        trait_impact = "neutral"
    else:
        trait_impact = "neutral"

    return trait_impact, effect_text, recommendation


def analyze_fitness(user_snps):
    """Analyze fitness and recovery SNPs.

    Args:
        user_snps: dict of rsid -> {'genotype': str, ...}

    Returns:
        dict with fitness analysis results
    """
    results = {cat: [] for cat in CATEGORY_LABELS}
    training_recommendations = []
    injury_alerts = []

    for rsid, entry in FITNESS_SNPS.items():
        # Handle cross-referenced SNPs with suffixes
        lookup_rsid = rsid.split("_")[0] if "_" in rsid else rsid

        if lookup_rsid in user_snps:
            genotype = user_snps[lookup_rsid]["genotype"]
            trait_impact, effect_text, recommendation = _assess_fitness_trait(genotype, entry)

            finding = {
                "rsid": lookup_rsid,
                "gene": entry["gene"],
                "genotype": genotype,
                "condition": entry["condition"],
                "trait_impact": trait_impact,
                "effect": effect_text,
                "recommendation": recommendation,
                "source": entry.get("source", ""),
            }

            category = entry["category"]
            if category in results:
                results[category].append(finding)

            if recommendation:
                training_recommendations.append({
                    "gene": entry["gene"],
                    "recommendation": recommendation,
                    "category": category,
                })

            if category == "injury_risk" and trait_impact == "attention_needed":
                injury_alerts.append(finding)

    # Determine training profile
    profile = _determine_training_profile(user_snps)

    return {
        "findings": results,
        "training_profile": profile,
        "training_recommendations": training_recommendations,
        "injury_alerts": injury_alerts,
        "category_labels": CATEGORY_LABELS,
    }


def _determine_training_profile(user_snps):
    """Determine overall genetic training profile (power vs endurance)."""
    power_score = 0
    endurance_score = 0

    # ACTN3
    actn3 = user_snps.get("rs1815739", {}).get("genotype", "")
    if actn3 == "CC":
        power_score += 2
    elif actn3 in ("CT", "TC"):
        power_score += 1
        endurance_score += 1
    elif actn3 == "TT":
        endurance_score += 2

    # ACE
    ace = user_snps.get("rs4341", {}).get("genotype", "")
    if ace == "GG":
        power_score += 1
    elif ace == "CC":
        endurance_score += 1

    # PPARGC1A
    pparg = user_snps.get("rs8192678", {}).get("genotype", "")
    if pparg == "GG":
        endurance_score += 1
    elif pparg == "AA":
        pass  # No bonus

    total = power_score + endurance_score
    if total == 0:
        return {
            "type": "unknown",
            "description": "Insufficient data to determine training profile",
            "power_score": 0,
            "endurance_score": 0,
        }

    if power_score > endurance_score + 1:
        profile_type = "power"
        description = (
            "Genetically favored for power/sprint activities. "
            "May excel in explosive movements, sprinting, weightlifting."
        )
    elif endurance_score > power_score + 1:
        profile_type = "endurance"
        description = (
            "Genetically favored for endurance activities. "
            "May excel in long-distance running, cycling, swimming."
        )
    else:
        profile_type = "mixed"
        description = (
            "Balanced genetic profile for both power and endurance. "
            "Versatile training response expected."
        )

    return {
        "type": profile_type,
        "description": description,
        "power_score": power_score,
        "endurance_score": endurance_score,
    }
