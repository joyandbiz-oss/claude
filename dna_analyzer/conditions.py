"""Cross-reference genetic findings with existing health conditions."""


# Map of conditions to relevant genetic categories and SNPs
CONDITION_SNP_RELEVANCE = {
    "cervical_disc_herniation": {
        "label": "Cervical Disc Herniations (C3-C7)",
        "relevant_genes": {
            "COL5A1": "Collagen variant may contribute to disc degeneration susceptibility",
            "COL1A1": "Collagen type I variants affect connective tissue integrity",
            "IL6": "Pro-inflammatory IL-6 contributes to disc inflammation and degeneration",
            "IL10": "Reduced anti-inflammatory IL-10 may worsen disc inflammation",
            "VDR": "Vitamin D receptor variants affect disc health and collagen metabolism",
            "SOD2": "Oxidative stress accelerates disc degeneration",
            "COMT": "Pain sensitivity varies with COMT status — affects pain management approach",
            "OPRM1": "Opioid receptor variant affects pain medication response",
        },
        "recommendations": [
            "Optimize vitamin D levels (especially with VDR variants)",
            "Anti-inflammatory support: omega-3, curcumin, boswellia",
            "Collagen support: vitamin C, proline, glycine, collagen peptides",
            "If COMT slow (Met/Met): may be more pain-sensitive, non-opioid approaches preferred",
            "Address SOD2 status — oxidative stress worsens disc degeneration",
        ],
    },
    "lumbar_disc_herniation": {
        "label": "Lumbar Disc Herniation L5-S1, Protrusion L4-L5",
        "relevant_genes": {
            "COL5A1": "Connective tissue integrity — affects disc and ligament health",
            "COL1A1": "Bone and connective tissue composition",
            "IL6": "Inflammation mediator in disc degeneration",
            "IL10": "Anti-inflammatory capacity",
            "VDR": "Vitamin D metabolism affects disc hydration and repair",
            "SOD2": "Mitochondrial oxidative stress in disc cells",
            "COMT": "Pain perception modulation",
            "MTHFR": "Homocysteine elevation may impair collagen cross-linking",
        },
        "recommendations": [
            "Same collagen/anti-inflammatory support as cervical",
            "MTHFR variants with elevated homocysteine can impair collagen quality",
            "Prioritize methylfolate if MTHFR impaired (supports collagen synthesis)",
            "Weight management — IL6/FTO variants may affect weight control difficulty",
            "Gentle core strengthening adapted to injury risk profile",
        ],
    },
    "seborrheic_dermatitis": {
        "label": "Seborrheic Dermatitis",
        "relevant_genes": {
            "VDR": "Vitamin D is immunomodulatory — deficiency worsens skin conditions",
            "IL6": "Inflammatory mediator — drives skin inflammation",
            "IL10": "Anti-inflammatory — low IL10 may contribute to chronic skin inflammation",
            "MTHFR": "Impaired methylation can affect skin barrier and immune regulation",
            "SOD2": "Oxidative stress contributes to skin inflammation",
            "GSTP1": "Detoxification capacity affects skin health",
            "HLA": "HLA variants may predispose to autoimmune skin manifestations",
            "FUT2": "Non-secretor status affects gut microbiome → skin axis",
        },
        "recommendations": [
            "Optimize vitamin D (target 50-70 ng/mL with VDR variants)",
            "Support gut-skin axis: probiotics (L. rhamnosus, B. longum)",
            "Zinc 15-30mg (antifungal + anti-inflammatory)",
            "Biotin and B-vitamins (especially if MTHFR impaired)",
            "Omega-3 fatty acids for inflammation control",
            "Consider food sensitivity testing (especially if HLA-DQ2/DQ8 positive)",
            "Selenium 200mcg (supports GPx for skin antioxidant defense)",
        ],
    },
    "elevated_monocytes": {
        "label": "Elevated Monocytes (12.4%, ref 2-10%)",
        "relevant_genes": {
            "IL6": "IL-6 stimulates monocyte production — chronic inflammation driver",
            "IL10": "Low IL-10 reduces monocyte regulation",
            "MTHFR": "Elevated homocysteine activates monocytes",
            "SOD2": "Oxidative stress activates innate immune cells",
            "COMT": "Slow COMT → higher catecholamines → immune activation",
            "PTPN22": "Autoimmune predisposition variant",
            "HFE": "Iron overload activates monocytes/macrophages",
        },
        "recommendations": [
            "Check homocysteine levels (MTHFR variants may contribute)",
            "Comprehensive inflammation panel: hs-CRP, IL-6, ferritin, ESR",
            "Check iron studies (ferritin, transferrin sat) if HFE variants present",
            "Anti-inflammatory protocol: omega-3, curcumin, vitamin D",
            "Rule out chronic infection or subclinical autoimmune process",
            "If PTPN22 variant: consider ANA, anti-thyroid antibodies screening",
        ],
    },
    "elevated_plcr": {
        "label": "Elevated P-LCR (48.5%, ref 13-43%)",
        "relevant_genes": {
            "IL6": "IL-6 affects thrombopoiesis and platelet size",
            "IL10": "Anti-inflammatory regulation of platelet production",
            "MTHFR": "Elevated homocysteine affects platelet function",
            "PON1": "Paraoxonase affects platelet oxidative stress",
            "9p21": "Cardiovascular risk locus may affect platelet biology",
        },
        "recommendations": [
            "P-LCR elevation may indicate platelet activation or chronic inflammation",
            "Check MPV (mean platelet volume) and platelet count together",
            "If homocysteine elevated (MTHFR): methylfolate + B12 to normalize",
            "Omega-3 fatty acids have mild antiplatelet effects",
            "Discuss with hematologist if persistently elevated",
            "Monitor cardiovascular risk markers given 9p21 genotype",
        ],
    },
}


def cross_reference_conditions(conditions, all_analysis_results, user_snps):
    """Cross-reference user's existing conditions with genetic findings.

    Args:
        conditions: list of condition keys from CONDITION_SNP_RELEVANCE
        all_analysis_results: dict of all analysis results
        user_snps: dict of rsid -> snp info

    Returns:
        list of condition cross-reference reports
    """
    reports = []

    for condition_key in conditions:
        if condition_key not in CONDITION_SNP_RELEVANCE:
            continue

        condition_info = CONDITION_SNP_RELEVANCE[condition_key]
        relevant_findings = []

        # Find genetic findings relevant to this condition
        for gene_name, relevance_note in condition_info["relevant_genes"].items():
            gene_findings = _find_findings_for_gene(gene_name, all_analysis_results)
            for finding in gene_findings:
                relevant_findings.append({
                    **finding,
                    "relevance_to_condition": relevance_note,
                })

        reports.append({
            "condition_key": condition_key,
            "condition_label": condition_info["label"],
            "relevant_findings": relevant_findings,
            "genetic_recommendations": condition_info["recommendations"],
            "finding_count": len(relevant_findings),
        })

    return reports


def _find_findings_for_gene(gene_name, all_analysis_results):
    """Find all genetic findings related to a gene name across all analyses."""
    findings = []
    gene_lower = gene_name.lower()

    for analysis_name, analysis in all_analysis_results.items():
        if "findings" not in analysis:
            continue
        for category, category_findings in analysis["findings"].items():
            for finding in category_findings:
                finding_gene = finding.get("gene", "").lower()
                if gene_lower in finding_gene or finding_gene.startswith(gene_lower):
                    findings.append({
                        "analysis_type": analysis_name,
                        "gene": finding.get("gene", ""),
                        "rsid": finding.get("rsid", ""),
                        "genotype": finding.get("genotype", ""),
                        "effect": finding.get("effect", ""),
                    })

    return findings
