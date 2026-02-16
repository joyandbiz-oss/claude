"""Health risk analysis — cardiovascular, diabetes, cancer, autoimmune, neurological."""

from .snp_database import HEALTH_RISK_SNPS


CATEGORY_LABELS = {
    "cardiovascular": "Cardiovascular Health",
    "diabetes": "Type 2 Diabetes",
    "cancer": "Cancer Predisposition",
    "autoimmune": "Autoimmune Conditions",
    "neurological": "Neurological / Mental Health",
}

RISK_LEVELS = {
    "high": "HIGH RISK",
    "moderate": "MODERATE RISK",
    "low": "LOW RISK",
    "average": "AVERAGE RISK",
    "protective": "PROTECTIVE",
    "uncertain": "UNCERTAIN",
}


def _classify_risk(genotype, entry):
    """Classify risk level based on genotype and SNP entry."""
    effect = entry.get("effect", {})
    effect_text = effect.get(genotype, "")
    text_lower = effect_text.lower()

    if not effect_text:
        return "uncertain", effect_text

    if any(w in text_lower for w in ["strongly increased", "high risk", "~12x",
                                      "critical", "non-functional", "no "]):
        return "high", effect_text
    if any(w in text_lower for w in ["increased risk", "~1.6x", "~1.7x",
                                      "reduced", "higher", "elevated",
                                      "homozygous", "poor"]):
        if any(w in text_lower for w in ["slightly", "mildly", "moderate"]):
            return "moderate", effect_text
        return "high", effect_text
    if any(w in text_lower for w in ["slightly increased", "moderately",
                                      "intermediate", "carrier"]):
        return "moderate", effect_text
    if any(w in text_lower for w in ["protective", "reduced risk"]):
        return "protective", effect_text
    if any(w in text_lower for w in ["average", "normal", "common"]):
        return "low", effect_text

    return "moderate", effect_text


def analyze_health_risks(user_snps):
    """Analyze health risk SNPs from user's genotype data.

    Args:
        user_snps: dict of rsid -> {'genotype': str, ...}

    Returns:
        dict with categorized health risk findings
    """
    results = {cat: [] for cat in CATEGORY_LABELS}
    not_found = []

    for rsid, entry in HEALTH_RISK_SNPS.items():
        if rsid in user_snps:
            genotype = user_snps[rsid]["genotype"]
            risk_level, effect_text = _classify_risk(genotype, entry)

            finding = {
                "rsid": rsid,
                "gene": entry["gene"],
                "genotype": genotype,
                "condition": entry["condition"],
                "risk_level": risk_level,
                "effect": effect_text,
                "clinical_significance": entry.get("clinical_significance", ""),
                "source": entry.get("source", ""),
            }
            category = entry["category"]
            if category in results:
                results[category].append(finding)
        else:
            not_found.append(rsid)

    # Sort each category by risk level (high first)
    risk_order = {"high": 0, "moderate": 1, "low": 2, "average": 3, "protective": 4, "uncertain": 5}
    for cat in results:
        results[cat].sort(key=lambda x: risk_order.get(x["risk_level"], 5))

    # Compute APOE haplotype if both SNPs present
    apoe = _determine_apoe(user_snps)

    return {
        "findings": results,
        "apoe_status": apoe,
        "snps_tested": len(HEALTH_RISK_SNPS) - len(not_found),
        "snps_missing": len(not_found),
        "category_labels": CATEGORY_LABELS,
    }


def _determine_apoe(user_snps):
    """Determine APOE haplotype from rs429358 and rs7412."""
    if "rs429358" not in user_snps or "rs7412" not in user_snps:
        return None

    e4_snp = user_snps["rs429358"]["genotype"]  # C = E4
    e2_snp = user_snps["rs7412"]["genotype"]      # T = E2

    # Determine alleles
    # rs429358: T=E3/E2, C=E4
    # rs7412: C=E3/E4, T=E2
    haplotypes = []

    if e4_snp == "TT" and e2_snp == "CC":
        haplotypes = ["E3", "E3"]
    elif e4_snp == "TT" and e2_snp == "CT":
        haplotypes = ["E3", "E2"]
    elif e4_snp == "TT" and e2_snp == "TT":
        haplotypes = ["E2", "E2"]
    elif e4_snp == "CT" and e2_snp == "CC":
        haplotypes = ["E3", "E4"]
    elif e4_snp == "CC" and e2_snp == "CC":
        haplotypes = ["E4", "E4"]
    elif e4_snp == "CT" and e2_snp == "CT":
        haplotypes = ["E2", "E4"]
    else:
        return {"genotype": f"rs429358={e4_snp}, rs7412={e2_snp}", "interpretation": "Unusual combination"}

    haplotype_str = "/".join(sorted(haplotypes))

    risk_map = {
        "E2/E2": "Lowest Alzheimer's risk. Possible hyperlipoproteinemia type III risk.",
        "E2/E3": "Below average Alzheimer's risk. Mildly protective.",
        "E3/E3": "Most common genotype. Average Alzheimer's risk.",
        "E3/E4": "~3x increased Alzheimer's risk. Earlier onset possible.",
        "E4/E4": "~12x increased Alzheimer's risk. Discuss with doctor. Consider lifestyle interventions.",
        "E2/E4": "Risk partially offset. Complex interaction.",
    }

    return {
        "haplotype": haplotype_str,
        "interpretation": risk_map.get(haplotype_str, "Unknown combination"),
        "rs429358": e4_snp,
        "rs7412": e2_snp,
    }


def get_high_risk_findings(analysis_results):
    """Extract only high-risk findings from analysis results."""
    high_risk = []
    for category, findings in analysis_results["findings"].items():
        for f in findings:
            if f["risk_level"] == "high":
                high_risk.append(f)
    return high_risk
