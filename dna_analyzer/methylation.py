"""Methylation & detox pathways analysis — MTHFR, COMT, CBS, SOD2, GST, etc."""

from .snp_database import METHYLATION_DETOX_SNPS


CATEGORY_LABELS = {
    "methylation": "Methylation Cycle",
    "transsulfuration": "Transsulfuration Pathway (CBS)",
    "antioxidant_defense": "Antioxidant Defense (SOD2)",
    "detoxification": "Phase II Detoxification (GST, NQO1, NAT2)",
}


def _assess_pathway_status(genotype, entry):
    """Assess pathway function based on genotype."""
    effect = entry.get("effect", {})
    effect_text = effect.get(genotype, "")
    supplement_rec = entry.get("supplement_recommendation", {}).get(genotype, "")

    text_lower = effect_text.lower()

    if any(w in text_lower for w in ["~70% reduced", "no ", "non-functional",
                                      "significantly reduced", "impaired",
                                      "slow comt", "upregulated"]):
        status = "impaired"
    elif any(w in text_lower for w in ["~35% reduced", "mildly", "moderately",
                                        "intermediate", "reduced"]):
        status = "suboptimal"
    elif any(w in text_lower for w in ["normal", "standard"]):
        status = "normal"
    else:
        status = "suboptimal"

    return status, effect_text, supplement_rec


def analyze_methylation_detox(user_snps):
    """Analyze methylation and detox pathway SNPs.

    Args:
        user_snps: dict of rsid -> {'genotype': str, ...}

    Returns:
        dict with methylation/detox analysis results
    """
    results = {cat: [] for cat in CATEGORY_LABELS}
    supplement_protocol = []
    pathway_warnings = []

    for rsid, entry in METHYLATION_DETOX_SNPS.items():
        lookup_rsid = rsid.split("_")[0] if "_" in rsid else rsid

        if lookup_rsid in user_snps:
            genotype = user_snps[lookup_rsid]["genotype"]
            status, effect_text, supplement_rec = _assess_pathway_status(genotype, entry)

            finding = {
                "rsid": lookup_rsid,
                "db_key": rsid,
                "gene": entry["gene"],
                "genotype": genotype,
                "condition": entry["condition"],
                "pathway_status": status,
                "effect": effect_text,
                "supplement_recommendation": supplement_rec,
                "source": entry.get("source", ""),
            }

            category = entry["category"]
            if category in results:
                results[category].append(finding)

            if supplement_rec and status in ("impaired", "suboptimal"):
                supplement_protocol.append({
                    "gene": entry["gene"],
                    "status": status,
                    "recommendation": supplement_rec,
                    "priority": "high" if status == "impaired" else "moderate",
                })

            if status == "impaired":
                pathway_warnings.append(finding)

    # Check for COMT + MTHFR interaction
    interactions = _check_pathway_interactions(user_snps)

    # Assess overall methylation capacity
    methylation_score = _score_methylation(user_snps)

    return {
        "findings": results,
        "supplement_protocol": supplement_protocol,
        "pathway_warnings": pathway_warnings,
        "interactions": interactions,
        "methylation_score": methylation_score,
        "category_labels": CATEGORY_LABELS,
    }


def _check_pathway_interactions(user_snps):
    """Check for important gene-gene interactions in methylation/detox."""
    interactions = []

    # COMT slow + MTHFR impaired = be careful with methyl donors
    comt = user_snps.get("rs4680", {}).get("genotype", "")
    mthfr_677 = user_snps.get("rs1801133", {}).get("genotype", "")
    mthfr_1298 = user_snps.get("rs1801131", {}).get("genotype", "")

    comt_slow = comt in ("AA",)
    mthfr_impaired = mthfr_677 in ("AA",) or (
        mthfr_677 in ("AG", "GA") and mthfr_1298 in ("GT", "TG")
    )

    if comt_slow and mthfr_impaired:
        interactions.append({
            "genes": "COMT (slow) + MTHFR (impaired)",
            "warning": (
                "Slow COMT with impaired MTHFR creates a methylation paradox: "
                "you need methyl support (MTHFR), but excess methyl donors can "
                "overwhelm slow COMT causing anxiety/overstimulation. "
                "Start methylfolate at LOW doses (100-200mcg) and titrate slowly. "
                "Avoid high-dose SAMe. Support with magnesium glycinate and B6 (P5P)."
            ),
            "severity": "important",
        })

    # CBS upregulated + MTHFR impaired = sulfur sensitivity
    cbs = user_snps.get("rs234706", {}).get("genotype", "")
    if cbs in ("AA",) and mthfr_impaired:
        interactions.append({
            "genes": "CBS (upregulated) + MTHFR (impaired)",
            "warning": (
                "Upregulated CBS pulls homocysteine toward transsulfuration, "
                "while impaired MTHFR can't efficiently recycle it. "
                "This can deplete methionine/SAMe. Limit sulfur supplements "
                "(NAC, MSM) and support methylation with low-dose methylfolate "
                "and methylcobalamin. Add molybdenum for sulfite processing."
            ),
            "severity": "important",
        })

    # SOD2 Val/Val + impaired GST = high oxidative stress
    sod2 = user_snps.get("rs4880", {}).get("genotype", "")
    gstp1 = user_snps.get("rs1695", {}).get("genotype", "")
    if sod2 in ("CC",) and gstp1 in ("GG", "AG"):
        interactions.append({
            "genes": "SOD2 (Val/Val) + GSTP1 (reduced)",
            "warning": (
                "Impaired mitochondrial antioxidant defense (SOD2) combined with "
                "reduced glutathione conjugation (GSTP1) = high oxidative stress. "
                "Priority: CoQ10 (ubiquinol) 200-400mg, NAC 600-1200mg, "
                "selenium 200mcg, vitamin C 1-2g, glutathione support."
            ),
            "severity": "important",
        })

    return interactions


def _score_methylation(user_snps):
    """Score overall methylation capacity (0-10, 10 = best)."""
    score = 10
    details = []

    # MTHFR C677T
    g = user_snps.get("rs1801133", {}).get("genotype", "")
    if g == "AA":
        score -= 3
        details.append("MTHFR C677T homozygous (-3)")
    elif g in ("AG", "GA"):
        score -= 1.5
        details.append("MTHFR C677T heterozygous (-1.5)")

    # MTHFR A1298C
    g = user_snps.get("rs1801131", {}).get("genotype", "")
    if g in ("GG", "CC"):
        score -= 1.5
        details.append("MTHFR A1298C homozygous (-1.5)")
    elif g in ("GT", "TG", "AC", "CA"):
        score -= 0.5
        details.append("MTHFR A1298C heterozygous (-0.5)")

    # MTR
    g = user_snps.get("rs1805087", {}).get("genotype", "")
    if g == "GG":
        score -= 1
        details.append("MTR A2756G homozygous (-1)")
    elif g in ("AG", "GA"):
        score -= 0.5
        details.append("MTR A2756G heterozygous (-0.5)")

    # MTRR
    g = user_snps.get("rs2066470", {}).get("genotype", "")
    if g == "AA":
        score -= 1
        details.append("MTRR homozygous (-1)")
    elif g in ("AG", "GA"):
        score -= 0.5
        details.append("MTRR heterozygous (-0.5)")

    # COMT slow (uses up methyl groups faster)
    g = user_snps.get("rs4680", {}).get("genotype", "")
    if g == "AA":
        score -= 0.5
        details.append("COMT slow — higher SAMe demand (-0.5)")

    score = max(0, round(score, 1))

    if score >= 8:
        interpretation = "Good methylation capacity"
    elif score >= 6:
        interpretation = "Mildly impaired methylation — targeted support recommended"
    elif score >= 4:
        interpretation = "Moderately impaired methylation — supplement protocol recommended"
    else:
        interpretation = "Significantly impaired methylation — comprehensive protocol needed"

    return {
        "score": score,
        "max_score": 10,
        "interpretation": interpretation,
        "details": details,
    }
