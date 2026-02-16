#!/usr/bin/env python3
"""DNA Health Risk Analysis Tool — Main CLI entry point.

Parses raw DNA genotyping data and generates a comprehensive health report
covering health risks, nutrigenomics, pharmacogenomics, fitness/recovery,
and methylation/detox pathways.

Usage:
    python analyze_dna.py <path_to_dna_file> [options]

Examples:
    python analyze_dna.py ~/Downloads/MyHeritage_raw_dna_data.zip
    python analyze_dna.py ~/Downloads/23andme_data.txt --format json --output report.json
    python analyze_dna.py data.csv --conditions cervical_disc,elevated_monocytes
"""

import argparse
import os
import sys

from dna_analyzer.parser import parse_raw_dna_file, filter_snps
from dna_analyzer.snp_database import get_all_tracked_rsids
from dna_analyzer.health_risks import analyze_health_risks
from dna_analyzer.nutrigenomics import analyze_nutrigenomics, get_mthfr_status
from dna_analyzer.pharmacogenomics import analyze_pharmacogenomics
from dna_analyzer.fitness import analyze_fitness
from dna_analyzer.methylation import analyze_methylation_detox
from dna_analyzer.conditions import cross_reference_conditions, CONDITION_SNP_RELEVANCE
from dna_analyzer.report import generate_report


DEFAULT_CONDITIONS = [
    "cervical_disc_herniation",
    "lumbar_disc_herniation",
    "seborrheic_dermatitis",
    "elevated_monocytes",
    "elevated_plcr",
]


def main():
    parser = argparse.ArgumentParser(
        description="DNA Health Risk Analysis Tool — Analyze raw DNA genotyping data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Supported providers: MyHeritage, 23andMe, AncestryDNA, FTDNA
Supported formats: .txt, .csv, .zip (containing .txt or .csv)

Available conditions for cross-referencing:
  cervical_disc_herniation  — Cervical disc herniations
  lumbar_disc_herniation    — Lumbar disc herniation
  seborrheic_dermatitis     — Seborrheic dermatitis
  elevated_monocytes        — Elevated monocytes
  elevated_plcr             — Elevated P-LCR

DISCLAIMER: This tool is for educational/research purposes only.
Not a substitute for clinical genetic testing or medical advice.
        """,
    )

    parser.add_argument(
        "dna_file",
        help="Path to raw DNA data file (.txt, .csv, or .zip)",
    )
    parser.add_argument(
        "--provider",
        choices=["auto", "myheritage", "23andme", "ancestrydna", "ftdna"],
        default="auto",
        help="DNA data provider (default: auto-detect)",
    )
    parser.add_argument(
        "--format",
        choices=["text", "json"],
        default="text",
        help="Output format (default: text)",
    )
    parser.add_argument(
        "--output", "-o",
        help="Output file path (default: stdout)",
    )
    parser.add_argument(
        "--conditions",
        help="Comma-separated list of conditions to cross-reference "
             "(default: all predefined conditions). Use 'none' to skip.",
    )
    parser.add_argument(
        "--filter-only",
        action="store_true",
        help="Only extract clinically relevant SNPs and output as CSV "
             "(skip full analysis)",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed progress",
    )

    args = parser.parse_args()

    # Validate input file
    if not os.path.exists(args.dna_file):
        print(f"Error: File not found: {args.dna_file}", file=sys.stderr)
        sys.exit(1)

    # Parse DNA file
    _log(args.verbose, f"Parsing DNA file: {args.dna_file}")
    parsed_data = parse_raw_dna_file(args.dna_file, provider=args.provider)
    _log(args.verbose, f"Provider: {parsed_data['provider']}")
    _log(args.verbose, f"Total SNPs loaded: {parsed_data['total_snps']:,}")

    # Filter to clinically relevant SNPs
    tracked_rsids = get_all_tracked_rsids()
    relevant_snps = filter_snps(parsed_data, tracked_rsids)
    _log(args.verbose, f"Clinically relevant SNPs found: {len(relevant_snps)}/{len(tracked_rsids)}")

    # Filter-only mode: just output the relevant SNPs as CSV
    if args.filter_only:
        _output_filtered_csv(relevant_snps, args.output)
        return

    # Run all analyses
    _log(args.verbose, "Running health risk analysis...")
    health_results = analyze_health_risks(relevant_snps)

    _log(args.verbose, "Running nutrigenomics analysis...")
    nutri_results = analyze_nutrigenomics(relevant_snps)
    mthfr_status = get_mthfr_status(relevant_snps)

    _log(args.verbose, "Running pharmacogenomics analysis...")
    pharma_results = analyze_pharmacogenomics(relevant_snps)

    _log(args.verbose, "Running fitness & recovery analysis...")
    fitness_results = analyze_fitness(relevant_snps)

    _log(args.verbose, "Running methylation & detox analysis...")
    methyl_results = analyze_methylation_detox(relevant_snps)

    # Determine conditions to cross-reference
    if args.conditions == "none":
        conditions = []
    elif args.conditions:
        conditions = [c.strip() for c in args.conditions.split(",")]
    else:
        conditions = DEFAULT_CONDITIONS

    _log(args.verbose, f"Cross-referencing with {len(conditions)} conditions...")
    all_analysis = {
        "health_risks": health_results,
        "nutrigenomics": nutri_results,
        "pharmacogenomics": pharma_results,
        "fitness": fitness_results,
        "methylation_detox": methyl_results,
    }
    condition_reports = cross_reference_conditions(conditions, all_analysis, relevant_snps)

    # Generate report
    _log(args.verbose, "Generating report...")
    report = generate_report(
        parsed_data, health_results, nutri_results,
        pharma_results, fitness_results, methyl_results,
        condition_reports, output_format=args.format,
    )

    # Output
    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write(report)
        print(f"Report written to: {args.output}")
    else:
        print(report)


def _output_filtered_csv(relevant_snps, output_path):
    """Output filtered SNPs as CSV for further analysis."""
    lines = ["rsid,chromosome,position,genotype"]
    for rsid, info in sorted(relevant_snps.items()):
        lines.append(f"{rsid},{info['chromosome']},{info['position']},{info['genotype']}")

    output = "\n".join(lines) + "\n"
    if output_path:
        with open(output_path, "w") as f:
            f.write(output)
        print(f"Filtered SNPs written to: {output_path}")
        print(f"Total clinically relevant SNPs: {len(relevant_snps)}")
    else:
        print(output)


def _log(verbose, message):
    """Print progress message if verbose mode enabled."""
    if verbose:
        print(f"[*] {message}", file=sys.stderr)


if __name__ == "__main__":
    main()
