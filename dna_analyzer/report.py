"""Generate comprehensive DNA analysis report."""

import json
from datetime import datetime


def generate_report(parsed_data, health_results, nutri_results,
                    pharma_results, fitness_results, methyl_results,
                    condition_reports, output_format="text"):
    """Generate a comprehensive analysis report.

    Args:
        parsed_data: Raw parsing results
        health_results: Health risk analysis
        nutri_results: Nutrigenomics analysis
        pharma_results: Pharmacogenomics analysis
        fitness_results: Fitness analysis
        methyl_results: Methylation/detox analysis
        condition_reports: Condition cross-reference reports
        output_format: 'text' or 'json'

    Returns:
        Formatted report string or dict
    """
    if output_format == "json":
        return _generate_json_report(
            parsed_data, health_results, nutri_results,
            pharma_results, fitness_results, methyl_results,
            condition_reports
        )

    return _generate_text_report(
        parsed_data, health_results, nutri_results,
        pharma_results, fitness_results, methyl_results,
        condition_reports
    )


def _generate_text_report(parsed_data, health_results, nutri_results,
                           pharma_results, fitness_results, methyl_results,
                           condition_reports):
    """Generate human-readable text report."""
    lines = []
    sep = "=" * 80
    subsep = "-" * 60

    lines.append(sep)
    lines.append("       DNA HEALTH RISK ANALYSIS REPORT")
    lines.append(f"       Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append(sep)
    lines.append("")
    lines.append("DISCLAIMER: This report is for EDUCATIONAL and RESEARCH purposes only.")
    lines.append("It is NOT a substitute for clinical genetic testing, genetic counseling,")
    lines.append("or medical advice. Always consult qualified healthcare professionals")
    lines.append("before making medical decisions based on genetic information.")
    lines.append("")

    # --- OVERVIEW ---
    lines.append(sep)
    lines.append("  1. DATA OVERVIEW")
    lines.append(sep)
    lines.append(f"  Provider detected: {parsed_data['provider']}")
    lines.append(f"  Total SNPs in file: {parsed_data['total_snps']:,}")
    lines.append(f"  Clinically relevant SNPs analyzed: {health_results['snps_tested']}"
                 f" (health) + nutrigenomics + pharmacogenomics + fitness + methylation")
    lines.append("")

    # --- CRITICAL ALERTS ---
    critical = pharma_results.get("critical_alerts", [])
    high_risk = []
    for cat, findings in health_results["findings"].items():
        for f in findings:
            if f["risk_level"] == "high":
                high_risk.append(f)

    if critical or high_risk:
        lines.append(sep)
        lines.append("  *** CRITICAL FINDINGS — DISCUSS WITH YOUR DOCTOR ***")
        lines.append(sep)
        lines.append("")

        if critical:
            lines.append("  PHARMACOGENOMIC ALERTS:")
            for alert in critical:
                lines.append(f"    ! {alert['gene']} ({alert['rsid']}): {alert['genotype']}")
                lines.append(f"      {alert['effect']}")
                for drug in alert.get("drugs_affected", []):
                    lines.append(f"      - {drug}")
                lines.append("")

        if high_risk:
            lines.append("  HIGH-RISK HEALTH FINDINGS:")
            for finding in high_risk:
                lines.append(f"    ! {finding['gene']} ({finding['rsid']}): {finding['genotype']}")
                lines.append(f"      {finding['condition']}: {finding['effect']}")
            lines.append("")

    # --- HEALTH RISKS ---
    lines.append(sep)
    lines.append("  2. HEALTH RISK ANALYSIS")
    lines.append(sep)

    if health_results.get("apoe_status"):
        apoe = health_results["apoe_status"]
        lines.append("")
        lines.append(f"  APOE Status: {apoe.get('haplotype', 'Unknown')}")
        lines.append(f"  {apoe.get('interpretation', '')}")
        lines.append("")

    for category, label in health_results["category_labels"].items():
        findings = health_results["findings"].get(category, [])
        if not findings:
            continue
        lines.append(f"  {label}:")
        lines.append(subsep)
        for f in findings:
            risk_marker = _risk_marker(f["risk_level"])
            lines.append(f"  {risk_marker} {f['gene']} ({f['rsid']}): {f['genotype']}")
            lines.append(f"     {f['condition']}")
            lines.append(f"     {f['effect']}")
            lines.append("")

    # --- NUTRIGENOMICS ---
    lines.append(sep)
    lines.append("  3. NUTRIGENOMICS")
    lines.append(sep)

    for category, label in nutri_results["category_labels"].items():
        findings = nutri_results["findings"].get(category, [])
        if not findings:
            continue
        lines.append(f"\n  {label}:")
        lines.append(subsep)
        for f in findings:
            impact_marker = _impact_marker(f["impact"])
            lines.append(f"  {impact_marker} {f['gene']} ({f['rsid']}): {f['genotype']}")
            lines.append(f"     {f['effect']}")
            if f.get("recommendation"):
                lines.append(f"     >> {f['recommendation']}")
            lines.append("")

    # --- PHARMACOGENOMICS ---
    lines.append(sep)
    lines.append("  4. PHARMACOGENOMICS")
    lines.append(sep)

    # CYP Summary
    cyp_summary = pharma_results.get("cyp_summary", {})
    if cyp_summary:
        lines.append("\n  CYP Enzyme Summary:")
        lines.append(subsep)
        for enzyme, status in cyp_summary.items():
            lines.append(f"    {enzyme}: {status}")
        lines.append("")

    for category, label in pharma_results["category_labels"].items():
        findings = pharma_results["findings"].get(category, [])
        if not findings:
            continue
        lines.append(f"\n  {label}:")
        lines.append(subsep)
        for f in findings:
            urgency_marker = _urgency_marker(f["urgency"])
            lines.append(f"  {urgency_marker} {f['gene']} ({f['rsid']}): {f['genotype']}")
            lines.append(f"     Status: {f['metabolizer_label']}")
            lines.append(f"     {f['effect']}")
            if f.get("drugs_affected"):
                lines.append("     Drugs affected:")
                for drug in f["drugs_affected"]:
                    lines.append(f"       - {drug}")
            lines.append("")

    # --- FITNESS ---
    lines.append(sep)
    lines.append("  5. FITNESS & RECOVERY")
    lines.append(sep)

    profile = fitness_results.get("training_profile", {})
    if profile.get("type") != "unknown":
        lines.append(f"\n  Training Profile: {profile.get('type', '').upper()}")
        lines.append(f"  {profile.get('description', '')}")
        lines.append(f"  Power score: {profile.get('power_score', 0)} | "
                     f"Endurance score: {profile.get('endurance_score', 0)}")
        lines.append("")

    for category, label in fitness_results["category_labels"].items():
        findings = fitness_results["findings"].get(category, [])
        if not findings:
            continue
        lines.append(f"\n  {label}:")
        lines.append(subsep)
        for f in findings:
            lines.append(f"  {f['gene']} ({f['rsid']}): {f['genotype']}")
            lines.append(f"     {f['effect']}")
            if f.get("recommendation"):
                lines.append(f"     >> {f['recommendation']}")
            lines.append("")

    if fitness_results.get("injury_alerts"):
        lines.append("  INJURY RISK ALERTS:")
        for alert in fitness_results["injury_alerts"]:
            lines.append(f"    ! {alert['gene']}: {alert['effect']}")
            if alert.get("recommendation"):
                lines.append(f"      >> {alert['recommendation']}")
        lines.append("")

    # --- METHYLATION & DETOX ---
    lines.append(sep)
    lines.append("  6. METHYLATION & DETOX PATHWAYS")
    lines.append(sep)

    methyl_score = methyl_results.get("methylation_score", {})
    if methyl_score:
        lines.append(f"\n  Methylation Capacity Score: "
                     f"{methyl_score.get('score', '?')}/{methyl_score.get('max_score', 10)}")
        lines.append(f"  {methyl_score.get('interpretation', '')}")
        if methyl_score.get("details"):
            for detail in methyl_score["details"]:
                lines.append(f"    - {detail}")
        lines.append("")

    for category, label in methyl_results["category_labels"].items():
        findings = methyl_results["findings"].get(category, [])
        if not findings:
            continue
        lines.append(f"\n  {label}:")
        lines.append(subsep)
        for f in findings:
            status_marker = _pathway_marker(f["pathway_status"])
            lines.append(f"  {status_marker} {f['gene']} ({f['rsid']}): {f['genotype']}")
            lines.append(f"     {f['effect']}")
            if f.get("supplement_recommendation"):
                lines.append(f"     >> {f['supplement_recommendation']}")
            lines.append("")

    # Pathway interactions
    interactions = methyl_results.get("interactions", [])
    if interactions:
        lines.append("  IMPORTANT GENE-GENE INTERACTIONS:")
        lines.append(subsep)
        for interaction in interactions:
            lines.append(f"  ! {interaction['genes']}")
            lines.append(f"    {interaction['warning']}")
            lines.append("")

    # --- CONDITION CROSS-REFERENCE ---
    if condition_reports:
        lines.append(sep)
        lines.append("  7. CROSS-REFERENCE WITH EXISTING CONDITIONS")
        lines.append(sep)

        for report in condition_reports:
            lines.append(f"\n  {report['condition_label']}:")
            lines.append(subsep)

            if report["relevant_findings"]:
                lines.append("  Relevant genetic findings:")
                for rf in report["relevant_findings"]:
                    lines.append(f"    - {rf['gene']} ({rf['rsid']}): {rf['genotype']}")
                    lines.append(f"      Effect: {rf['effect']}")
                    lines.append(f"      Relevance: {rf['relevance_to_condition']}")
                lines.append("")

            lines.append("  Genetic-informed recommendations:")
            for rec in report["genetic_recommendations"]:
                lines.append(f"    * {rec}")
            lines.append("")

    # --- SUPPLEMENT PROTOCOL ---
    lines.append(sep)
    lines.append("  8. CONSOLIDATED SUPPLEMENT RECOMMENDATIONS")
    lines.append(sep)
    lines.append("")
    lines.append(_build_supplement_section(nutri_results, methyl_results, fitness_results))

    # --- WHAT TO DISCUSS WITH DOCTOR ---
    lines.append(sep)
    lines.append("  9. WHAT TO DISCUSS WITH YOUR DOCTOR")
    lines.append(sep)
    lines.append("")
    lines.append(_build_doctor_discussion(health_results, pharma_results,
                                           methyl_results, condition_reports))

    # --- FOOTER ---
    lines.append("")
    lines.append(sep)
    lines.append("  END OF REPORT")
    lines.append(sep)
    lines.append("")
    lines.append("  Sources: ClinVar, SNPedia, PharmGKB, CPIC Guidelines,")
    lines.append("  published GWAS studies. See individual entries for specific citations.")
    lines.append("")
    lines.append("  This analysis is NOT a diagnostic tool. Genetic risk factors are")
    lines.append("  probabilistic, not deterministic. Environmental factors, lifestyle,")
    lines.append("  and gene-gene interactions all modify actual risk.")
    lines.append("")

    return "\n".join(lines)


def _build_supplement_section(nutri_results, methyl_results, fitness_results):
    """Build consolidated supplement recommendations organized by priority."""
    lines = []
    high_priority = []
    moderate_priority = []

    # Collect from all sources
    for src in [nutri_results.get("supplement_recommendations", []),
                methyl_results.get("supplement_protocol", [])]:
        for rec in src:
            if rec.get("priority") == "high":
                high_priority.append(rec)
            else:
                moderate_priority.append(rec)

    if high_priority:
        lines.append("  HIGH PRIORITY:")
        for rec in high_priority:
            lines.append(f"    [{rec['gene']}] {rec['recommendation']}")
        lines.append("")

    if moderate_priority:
        lines.append("  MODERATE PRIORITY:")
        for rec in moderate_priority:
            lines.append(f"    [{rec['gene']}] {rec['recommendation']}")
        lines.append("")

    # Training-specific
    training_recs = fitness_results.get("training_recommendations", [])
    if training_recs:
        lines.append("  FITNESS-SPECIFIC:")
        for rec in training_recs:
            lines.append(f"    [{rec['gene']}] {rec['recommendation']}")
        lines.append("")

    if not lines:
        lines.append("  No specific supplement recommendations based on available data.")
        lines.append("")

    return "\n".join(lines)


def _build_doctor_discussion(health_results, pharma_results, methyl_results,
                              condition_reports):
    """Build list of items to discuss with doctor."""
    lines = []
    items = []

    # Critical pharma alerts
    for alert in pharma_results.get("critical_alerts", []):
        items.append(f"PHARMACOGENOMICS: {alert['gene']} — {alert['metabolizer_label']}. "
                     f"Affects: {', '.join(alert.get('drugs_affected', [])[:3])}")

    # High-risk health findings
    for cat, findings in health_results["findings"].items():
        for f in findings:
            if f["risk_level"] == "high":
                items.append(f"HEALTH RISK: {f['gene']} — {f['condition']}: {f['effect']}")

    # APOE
    apoe = health_results.get("apoe_status")
    if apoe and "E4" in apoe.get("haplotype", ""):
        items.append(f"APOE STATUS: {apoe['haplotype']} — {apoe['interpretation']}")

    # Methylation warnings
    methyl_score = methyl_results.get("methylation_score", {})
    if methyl_score.get("score", 10) < 6:
        items.append(f"METHYLATION: Score {methyl_score['score']}/10 — "
                     f"{methyl_score['interpretation']}. Request homocysteine test.")

    # Pathway interactions
    for interaction in methyl_results.get("interactions", []):
        items.append(f"GENE INTERACTION: {interaction['genes']} — {interaction['warning'][:100]}...")

    # Condition-specific
    if condition_reports:
        items.append("EXISTING CONDITIONS: Review genetic cross-reference findings for "
                     "personalized treatment optimization")

    if items:
        for i, item in enumerate(items, 1):
            lines.append(f"  {i}. {item}")
    else:
        lines.append("  No critical items identified. Standard preventive care recommended.")

    lines.append("")
    lines.append("  RECOMMENDED LAB TESTS:")
    lines.append("    - Homocysteine (fasting)")
    lines.append("    - 25(OH) Vitamin D")
    lines.append("    - hs-CRP (high-sensitivity C-reactive protein)")
    lines.append("    - Complete iron panel (ferritin, iron, TIBC, transferrin sat)")
    lines.append("    - Vitamin B12 and folate levels")
    lines.append("    - Comprehensive metabolic panel")
    lines.append("")

    return "\n".join(lines)


def _risk_marker(level):
    """Return a text marker for risk level."""
    markers = {
        "high": "[!!!]",
        "moderate": "[!!]",
        "low": "[OK]",
        "average": "[OK]",
        "protective": "[+]",
        "uncertain": "[?]",
    }
    return markers.get(level, "[?]")


def _impact_marker(impact):
    """Return a text marker for nutritional impact."""
    markers = {
        "significant": "[!!!]",
        "moderate": "[!!]",
        "minimal": "[OK]",
    }
    return markers.get(impact, "[?]")


def _urgency_marker(urgency):
    """Return a text marker for clinical urgency."""
    markers = {
        "critical": "[CRITICAL]",
        "high": "[HIGH]",
        "moderate": "[MOD]",
        "informational": "[INFO]",
    }
    return markers.get(urgency, "[?]")


def _pathway_marker(status):
    """Return a text marker for pathway status."""
    markers = {
        "impaired": "[!!!]",
        "suboptimal": "[!!]",
        "normal": "[OK]",
    }
    return markers.get(status, "[?]")


def _generate_json_report(parsed_data, health_results, nutri_results,
                           pharma_results, fitness_results, methyl_results,
                           condition_reports):
    """Generate JSON report."""
    report = {
        "generated_at": datetime.now().isoformat(),
        "disclaimer": (
            "This report is for EDUCATIONAL and RESEARCH purposes only. "
            "Not a substitute for clinical genetic testing or medical advice."
        ),
        "data_overview": {
            "provider": parsed_data["provider"],
            "total_snps": parsed_data["total_snps"],
        },
        "health_risks": health_results,
        "nutrigenomics": nutri_results,
        "pharmacogenomics": pharma_results,
        "fitness": fitness_results,
        "methylation_detox": methyl_results,
        "condition_cross_reference": condition_reports,
    }
    return json.dumps(report, indent=2, default=str)
