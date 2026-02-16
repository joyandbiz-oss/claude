"""Parse raw DNA genotyping data from various providers."""

import csv
import io
import os
import zipfile


SUPPORTED_PROVIDERS = ["myheritage", "23andme", "ancestrydna", "ftdna", "auto"]


def detect_provider(header_lines):
    """Detect the DNA data provider from file header lines."""
    joined = "\n".join(header_lines[:20]).lower()
    if "myheritage" in joined:
        return "myheritage"
    if "23andme" in joined:
        return "23andme"
    if "ancestrydna" in joined or "ancestry" in joined:
        return "ancestrydna"
    if "ftdna" in joined or "family tree dna" in joined:
        return "ftdna"
    return "unknown"


def _parse_line_fields(line, provider):
    """Parse a single data line into (rsid, chromosome, position, genotype)."""
    line = line.strip()
    if not line or line.startswith("#") or line.startswith("RSID") or line.startswith("rsid"):
        return None

    # Try tab-separated first, then comma
    if "\t" in line:
        parts = line.split("\t")
    elif "," in line:
        parts = line.split(",")
    else:
        parts = line.split()

    if len(parts) < 4:
        return None

    rsid = parts[0].strip().strip('"')
    chrom = parts[1].strip().strip('"')
    pos = parts[2].strip().strip('"')
    genotype = parts[3].strip().strip('"')

    # Some formats have allele1 and allele2 in separate columns
    if len(parts) >= 5 and len(parts[3].strip().strip('"')) == 1:
        allele1 = parts[3].strip().strip('"')
        allele2 = parts[4].strip().strip('"')
        genotype = allele1 + allele2

    if not rsid.startswith("rs") and not rsid.startswith("i"):
        return None

    return (rsid, chrom, pos, genotype)


def parse_raw_dna_file(filepath, provider="auto"):
    """Parse a raw DNA data file and return a dict of rsid -> genotype info.

    Args:
        filepath: Path to raw DNA file (can be .txt, .csv, or .zip)
        provider: One of 'myheritage', '23andme', 'ancestrydna', 'ftdna', 'auto'

    Returns:
        dict with keys:
            'provider': detected provider name
            'snps': dict of rsid -> {'chromosome': str, 'position': str, 'genotype': str}
            'total_snps': int
    """
    lines = _read_file_lines(filepath)

    if provider == "auto":
        provider = detect_provider(lines[:20])

    snps = {}
    skipped = 0

    for line in lines:
        result = _parse_line_fields(line, provider)
        if result is None:
            continue
        rsid, chrom, pos, genotype = result
        # Normalize genotype
        genotype = genotype.upper().replace("-", "").replace("0", "")
        if genotype and rsid:
            snps[rsid] = {
                "chromosome": chrom,
                "position": pos,
                "genotype": genotype,
            }
        else:
            skipped += 1

    return {
        "provider": provider,
        "snps": snps,
        "total_snps": len(snps),
        "skipped": skipped,
    }


def _read_file_lines(filepath):
    """Read lines from a file, handling zip archives."""
    if filepath.endswith(".zip"):
        return _read_zip_lines(filepath)
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        return f.readlines()


def _read_zip_lines(filepath):
    """Extract and read lines from a zip archive containing DNA data."""
    with zipfile.ZipFile(filepath, "r") as zf:
        # Find the data file inside the zip
        data_file = None
        for name in zf.namelist():
            lower = name.lower()
            if lower.endswith(".csv") or lower.endswith(".txt"):
                data_file = name
                break
        if data_file is None:
            # Just use the first file
            data_file = zf.namelist()[0]

        with zf.open(data_file) as f:
            content = f.read().decode("utf-8", errors="replace")
            return content.splitlines()


def filter_snps(parsed_data, rsid_list):
    """Filter parsed SNP data to only include specified rsIDs.

    Args:
        parsed_data: Output from parse_raw_dna_file()
        rsid_list: List of rsIDs to keep

    Returns:
        dict of rsid -> snp info for matching SNPs
    """
    rsid_set = set(rsid_list)
    return {
        rsid: info
        for rsid, info in parsed_data["snps"].items()
        if rsid in rsid_set
    }
