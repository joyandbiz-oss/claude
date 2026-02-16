"""Pharmacogenomics analysis — drug metabolism and medication sensitivities."""

from .snp_database import PHARMACOGENOMICS_SNPS


CATEGORY_LABELS = {
    "cyp2d6": "CYP2D6 (Codeine, Tamoxifen, Antidepressants, Beta-blockers)",
    "cyp2c19": "CYP2C19 (Clopidogrel, PPIs, SSRIs)",
    "cyp2c9": "CYP2C9 (Warfarin, NSAIDs, Phenytoin)",
    "cyp3a4": "CYP3A4 (Statins, Immunosuppressants)",
    "cyp3a5": "CYP3A5 (Tacrolimus)",
    "warfarin": "Warfarin Sensitivity (VKORC1)",
    "statins": "Statin Myopathy Risk (SLCO1B1)",
    "chemotherapy": "Chemotherapy Toxicity (DPYD)",
    "thiopurines": "Thiopurine Metabolism (TPMT)",
}


METABOLIZER_STATUS = {
    "poor": "Poor Metabolizer",
    "intermediate": "Intermediate Metabolizer",
    "normal": "Normal (Extensive) Metabolizer",
    "rapid": "Rapid Metabolizer",
    "ultra_rapid": "Ultra-Rapid Metabolizer",
}


def _determine_metabolizer_status(genotype, entry):
    """Determine metabolizer status from genotype."""
    effect = entry.get("effect", {})
    effect_text = effect.get(genotype, "")
    text_lower = effect_text.lower()

    if "poor" in text_lower or "non-functional" in text_lower or "deficient" in text_lower:
        return "poor", effect_text
    if "ultra-rapid" in text_lower or "ultra rapid" in text_lower:
        return "ultra_rapid", effect_text
    if "rapid" in text_lower:
        return "rapid", effect_text
    if "intermediate" in text_lower or "reduced" in text_lower or "one functional" in text_lower:
        return "intermediate", effect_text
    if "normal" in text_lower or "extensive" in text_lower or "reference" in text_lower:
        return "normal", effect_text
    if "highly sensitive" in text_lower:
        return "poor", effect_text
    if "moderately sensitive" in text_lower:
        return "intermediate", effect_text

    return "normal", effect_text


def _assess_clinical_urgency(metabolizer_status, entry):
    """Determine clinical urgency of a pharmacogenomic finding."""
    drugs = entry.get("drugs_affected", [])
    has_critical = any("CRITICAL" in d for d in drugs)

    if has_critical and metabolizer_status in ("poor", "ultra_rapid"):
        return "critical"
    if metabolizer_status in ("poor", "ultra_rapid"):
        return "high"
    if metabolizer_status == "intermediate":
        return "moderate"
    return "informational"


def analyze_pharmacogenomics(user_snps):
    """Analyze pharmacogenomics SNPs.

    Args:
        user_snps: dict of rsid -> {'genotype': str, ...}

    Returns:
        dict with categorized pharmacogenomic findings
    """
    results = {cat: [] for cat in CATEGORY_LABELS}
    critical_alerts = []
    drug_interactions = []

    for rsid, entry in PHARMACOGENOMICS_SNPS.items():
        if rsid in user_snps:
            genotype = user_snps[rsid]["genotype"]
            met_status, effect_text = _determine_metabolizer_status(genotype, entry)
            urgency = _assess_clinical_urgency(met_status, entry)

            finding = {
                "rsid": rsid,
                "gene": entry["gene"],
                "genotype": genotype,
                "condition": entry["condition"],
                "metabolizer_status": met_status,
                "metabolizer_label": METABOLIZER_STATUS.get(met_status, met_status),
                "effect": effect_text,
                "drugs_affected": entry.get("drugs_affected", []),
                "urgency": urgency,
                "source": entry.get("source", ""),
            }

            category = entry["category"]
            if category in results:
                results[category].append(finding)

            if urgency == "critical":
                critical_alerts.append(finding)

            # Collect drug interaction summaries
            if met_status != "normal":
                for drug in entry.get("drugs_affected", []):
                    drug_interactions.append({
                        "drug": drug,
                        "gene": entry["gene"],
                        "status": METABOLIZER_STATUS.get(met_status, met_status),
                        "urgency": urgency,
                    })

    # Sort by urgency
    urgency_order = {"critical": 0, "high": 1, "moderate": 2, "informational": 3}
    for cat in results:
        results[cat].sort(key=lambda x: urgency_order.get(x["urgency"], 3))

    drug_interactions.sort(key=lambda x: urgency_order.get(x["urgency"], 3))

    # Determine CYP2D6 composite status
    cyp2d6_status = _composite_cyp_status(user_snps, "cyp2d6",
                                           ["rs3892097", "rs16947", "rs1065852"])
    cyp2c19_status = _composite_cyp_status(user_snps, "cyp2c19",
                                            ["rs4244285", "rs4986893", "rs12248560"])
    cyp2c9_status = _composite_cyp_status(user_snps, "cyp2c9",
                                           ["rs1799853", "rs1057910"])

    return {
        "findings": results,
        "critical_alerts": critical_alerts,
        "drug_interactions": drug_interactions,
        "cyp_summary": {
            "CYP2D6": cyp2d6_status,
            "CYP2C19": cyp2c19_status,
            "CYP2C9": cyp2c9_status,
        },
        "category_labels": CATEGORY_LABELS,
    }


def _composite_cyp_status(user_snps, category, rsids):
    """Determine composite metabolizer status for a CYP enzyme from multiple SNPs."""
    statuses = []
    for rsid in rsids:
        if rsid in user_snps:
            entry = PHARMACOGENOMICS_SNPS.get(rsid)
            if entry:
                genotype = user_snps[rsid]["genotype"]
                status, _ = _determine_metabolizer_status(genotype, entry)
                statuses.append(status)

    if not statuses:
        return "Not tested"

    # Worst-case status
    status_priority = {"poor": 0, "ultra_rapid": 1, "intermediate": 2, "rapid": 3, "normal": 4}
    worst = min(statuses, key=lambda s: status_priority.get(s, 4))
    return METABOLIZER_STATUS.get(worst, worst)
