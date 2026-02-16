"""Clinically significant SNP reference database.

Sources: ClinVar, SNPedia, PharmGKB, published GWAS studies.
Each entry includes rsID, gene, risk allele(s), category, clinical significance,
and relevant references.

DISCLAIMER: This database is for educational/research purposes only.
Not a substitute for clinical genetic testing or medical advice.
"""


# --- HEALTH RISK SNPs ---
HEALTH_RISK_SNPS = {
    # === CARDIOVASCULAR ===
    "rs1333049": {
        "gene": "9p21.3 (CDKN2A/B)",
        "risk_allele": "C",
        "category": "cardiovascular",
        "condition": "Coronary artery disease",
        "effect": {
            "CC": "~1.6x increased risk of CAD",
            "CG": "~1.3x increased risk of CAD",
            "GG": "Average risk",
        },
        "clinical_significance": "pathogenic_risk",
        "source": "ClinVar, Wellcome Trust GWAS",
    },
    "rs10757278": {
        "gene": "9p21.3 (CDKN2A/B)",
        "risk_allele": "G",
        "category": "cardiovascular",
        "condition": "Myocardial infarction / CAD",
        "effect": {
            "GG": "~1.6x increased risk",
            "AG": "~1.3x increased risk",
            "AA": "Average risk",
        },
        "clinical_significance": "pathogenic_risk",
        "source": "McPherson et al. Science 2007",
    },
    "rs6922269": {
        "gene": "MTHFD1L",
        "risk_allele": "A",
        "category": "cardiovascular",
        "condition": "Coronary heart disease",
        "effect": {
            "AA": "Increased risk",
            "AG": "Slightly increased risk",
            "GG": "Average risk",
        },
        "clinical_significance": "risk_factor",
        "source": "GWAS meta-analysis",
    },
    "rs1746048": {
        "gene": "CXCL12",
        "risk_allele": "C",
        "category": "cardiovascular",
        "condition": "Coronary artery disease",
        "effect": {
            "CC": "Increased risk of CAD",
            "CT": "Slightly increased risk",
            "TT": "Average risk",
        },
        "clinical_significance": "risk_factor",
        "source": "Samani et al. NEJM 2007",
    },
    "rs9818870": {
        "gene": "MRAS",
        "risk_allele": "T",
        "category": "cardiovascular",
        "condition": "Coronary artery disease",
        "effect": {
            "TT": "Increased risk",
            "CT": "Slightly increased risk",
            "CC": "Average risk",
        },
        "clinical_significance": "risk_factor",
        "source": "Erdmann et al. 2009",
    },
    "rs1799945": {
        "gene": "HFE (H63D)",
        "risk_allele": "G",
        "category": "cardiovascular",
        "condition": "Hemochromatosis / iron overload",
        "effect": {
            "GG": "Homozygous H63D — iron overload risk",
            "CG": "Carrier — mild increased iron absorption",
            "CC": "Normal",
        },
        "clinical_significance": "pathogenic_risk",
        "source": "ClinVar",
    },
    "rs1800562": {
        "gene": "HFE (C282Y)",
        "risk_allele": "A",
        "category": "cardiovascular",
        "condition": "Hereditary hemochromatosis",
        "effect": {
            "AA": "Homozygous — high risk of hemochromatosis",
            "GA": "Carrier — compound heterozygote risk with H63D",
            "GG": "Normal",
        },
        "clinical_significance": "pathogenic",
        "source": "ClinVar",
    },
    "rs5082": {
        "gene": "APOA2",
        "risk_allele": "G",
        "category": "cardiovascular",
        "condition": "Saturated fat sensitivity / obesity",
        "effect": {
            "GG": "Increased BMI with high saturated fat intake",
            "AG": "Moderate sensitivity",
            "AA": "Normal response to saturated fat",
        },
        "clinical_significance": "risk_factor",
        "source": "Corella et al. 2009",
    },
    "rs662": {
        "gene": "PON1",
        "risk_allele": "G",
        "category": "cardiovascular",
        "condition": "Reduced antioxidant protection (paraoxonase)",
        "effect": {
            "GG": "Lower PON1 activity — reduced HDL protection",
            "AG": "Intermediate PON1 activity",
            "AA": "Normal PON1 activity — better HDL antioxidant function",
        },
        "clinical_significance": "risk_factor",
        "source": "SNPedia, Mackness et al.",
    },
    "rs1800795": {
        "gene": "IL6",
        "risk_allele": "C",
        "category": "cardiovascular",
        "condition": "Chronic inflammation / cardiovascular risk",
        "effect": {
            "CC": "Higher IL-6 production — increased inflammation",
            "CG": "Intermediate IL-6 levels",
            "GG": "Lower IL-6 production",
        },
        "clinical_significance": "risk_factor",
        "source": "Fishman et al. 1998",
    },
    "rs1800896": {
        "gene": "IL10",
        "risk_allele": "A",
        "category": "cardiovascular",
        "condition": "Reduced anti-inflammatory response",
        "effect": {
            "AA": "Lower IL-10 — reduced anti-inflammatory capacity",
            "AG": "Intermediate IL-10 levels",
            "GG": "Higher IL-10 — better anti-inflammatory response",
        },
        "clinical_significance": "risk_factor",
        "source": "Turner et al. 1997",
    },

    # === TYPE 2 DIABETES ===
    "rs7903146": {
        "gene": "TCF7L2",
        "risk_allele": "T",
        "category": "diabetes",
        "condition": "Type 2 diabetes",
        "effect": {
            "TT": "~1.7x increased risk of T2D (strongest known genetic risk)",
            "CT": "~1.4x increased risk",
            "CC": "Average risk",
        },
        "clinical_significance": "pathogenic_risk",
        "source": "Grant et al. Nature Genetics 2006",
    },
    "rs1801282": {
        "gene": "PPARG (Pro12Ala)",
        "risk_allele": "C",
        "category": "diabetes",
        "condition": "Type 2 diabetes",
        "effect": {
            "CC": "~1.25x increased risk of T2D",
            "CG": "Average risk",
            "GG": "Slightly reduced risk (protective Ala allele)",
        },
        "clinical_significance": "risk_factor",
        "source": "Altshuler et al. 2000",
    },
    "rs5219": {
        "gene": "KCNJ11 (E23K)",
        "risk_allele": "T",
        "category": "diabetes",
        "condition": "Type 2 diabetes",
        "effect": {
            "TT": "Increased risk — affects potassium channel in beta cells",
            "CT": "Slightly increased risk",
            "CC": "Average risk",
        },
        "clinical_significance": "risk_factor",
        "source": "Gloyn et al. 2003",
    },
    "rs13266634": {
        "gene": "SLC30A8",
        "risk_allele": "C",
        "category": "diabetes",
        "condition": "Type 2 diabetes",
        "effect": {
            "CC": "Increased risk — zinc transporter variant",
            "CT": "Slightly increased risk",
            "TT": "Reduced risk (protective)",
        },
        "clinical_significance": "risk_factor",
        "source": "Sladek et al. 2007",
    },
    "rs4402960": {
        "gene": "IGF2BP2",
        "risk_allele": "T",
        "category": "diabetes",
        "condition": "Type 2 diabetes",
        "effect": {
            "TT": "Increased risk of T2D",
            "GT": "Slightly increased risk",
            "GG": "Average risk",
        },
        "clinical_significance": "risk_factor",
        "source": "GWAS diabetes consortium",
    },
    "rs10811661": {
        "gene": "CDKN2A/B",
        "risk_allele": "T",
        "category": "diabetes",
        "condition": "Type 2 diabetes",
        "effect": {
            "TT": "Increased risk of T2D",
            "CT": "Slightly increased risk",
            "CC": "Average risk",
        },
        "clinical_significance": "risk_factor",
        "source": "DIAGRAM consortium",
    },

    # === CANCER PREDISPOSITION ===
    "rs1042522": {
        "gene": "TP53 (Arg72Pro)",
        "risk_allele": "C",
        "category": "cancer",
        "condition": "Various cancers (modified apoptosis)",
        "effect": {
            "CC": "Pro/Pro — altered apoptosis efficiency",
            "CG": "Arg/Pro — heterozygous",
            "GG": "Arg/Arg — more efficient apoptosis",
        },
        "clinical_significance": "risk_factor",
        "source": "Dumont et al. 2003, SNPedia",
    },
    "rs1799966": {
        "gene": "BRCA1",
        "risk_allele": "T",
        "category": "cancer",
        "condition": "Breast/ovarian cancer susceptibility",
        "effect": {
            "TT": "Potentially increased risk — discuss with genetic counselor",
            "CT": "Carrier — discuss with genetic counselor",
            "CC": "Common variant",
        },
        "clinical_significance": "risk_factor",
        "source": "ClinVar",
    },
    "rs766173": {
        "gene": "BRCA2",
        "risk_allele": "T",
        "category": "cancer",
        "condition": "Breast cancer susceptibility",
        "effect": {
            "TT": "Potentially increased risk",
            "CT": "Carrier variant",
            "CC": "Common variant",
        },
        "clinical_significance": "risk_factor",
        "source": "ClinVar",
    },
    "rs1800566": {
        "gene": "NQO1",
        "risk_allele": "T",
        "category": "cancer",
        "condition": "Benzene toxicity / cancer susceptibility",
        "effect": {
            "TT": "No NQO1 enzyme activity — increased susceptibility to carcinogens",
            "CT": "Reduced NQO1 activity",
            "CC": "Normal NQO1 activity",
        },
        "clinical_significance": "risk_factor",
        "source": "Ross et al. 2000",
    },
    "rs4986850": {
        "gene": "BRCA1 (missense)",
        "risk_allele": "A",
        "category": "cancer",
        "condition": "Breast cancer",
        "effect": {
            "AA": "Variant — consult genetic counselor",
            "AG": "Heterozygous",
            "GG": "Common genotype",
        },
        "clinical_significance": "uncertain_significance",
        "source": "ClinVar",
    },
    "rs1801516": {
        "gene": "ATM (D1853N)",
        "risk_allele": "A",
        "category": "cancer",
        "condition": "Breast cancer / radiation sensitivity",
        "effect": {
            "AA": "Increased radiation sensitivity — discuss screening",
            "AG": "Carrier — some increased sensitivity",
            "GG": "Normal ATM function",
        },
        "clinical_significance": "risk_factor",
        "source": "Renwick et al. 2006",
    },
    "rs2981582": {
        "gene": "FGFR2",
        "risk_allele": "T",
        "category": "cancer",
        "condition": "Breast cancer",
        "effect": {
            "TT": "~1.6x increased breast cancer risk",
            "CT": "~1.3x increased risk",
            "CC": "Average risk",
        },
        "clinical_significance": "risk_factor",
        "source": "Easton et al. Nature 2007",
    },
    "rs6983267": {
        "gene": "8q24 (MYC region)",
        "risk_allele": "G",
        "category": "cancer",
        "condition": "Colorectal cancer / prostate cancer",
        "effect": {
            "GG": "Increased risk of colorectal and prostate cancer",
            "GT": "Moderately increased risk",
            "TT": "Average risk",
        },
        "clinical_significance": "risk_factor",
        "source": "Tomlinson et al. 2007",
    },

    # === AUTOIMMUNE ===
    "rs2476601": {
        "gene": "PTPN22 (R620W)",
        "risk_allele": "A",
        "category": "autoimmune",
        "condition": "Multiple autoimmune diseases (RA, T1D, lupus, thyroid)",
        "effect": {
            "AA": "Strongly increased autoimmune risk",
            "AG": "~1.5-2x increased autoimmune risk",
            "GG": "Average risk",
        },
        "clinical_significance": "pathogenic_risk",
        "source": "Begovich et al. 2004",
    },
    "rs3087243": {
        "gene": "CTLA4",
        "risk_allele": "G",
        "category": "autoimmune",
        "condition": "Autoimmune thyroid disease, T1D",
        "effect": {
            "GG": "Increased autoimmune risk",
            "AG": "Moderately increased risk",
            "AA": "Average risk",
        },
        "clinical_significance": "risk_factor",
        "source": "Ueda et al. Nature 2003",
    },
    "rs2104286": {
        "gene": "IL2RA",
        "risk_allele": "A",
        "category": "autoimmune",
        "condition": "Multiple sclerosis, T1D",
        "effect": {
            "AA": "Increased autoimmune risk",
            "AG": "Slightly increased risk",
            "GG": "Average risk",
        },
        "clinical_significance": "risk_factor",
        "source": "GWAS studies",
    },
    "rs3184504": {
        "gene": "SH2B3 (LNK)",
        "risk_allele": "T",
        "category": "autoimmune",
        "condition": "Celiac disease, T1D, cardiovascular",
        "effect": {
            "TT": "Increased risk for celiac, T1D, and hypertension",
            "CT": "Moderately increased risk",
            "CC": "Average risk",
        },
        "clinical_significance": "risk_factor",
        "source": "Hunt et al. 2008",
    },
    "rs2187668": {
        "gene": "HLA-DQ2.5",
        "risk_allele": "T",
        "category": "autoimmune",
        "condition": "Celiac disease",
        "effect": {
            "TT": "Homozygous HLA-DQ2.5 — high celiac risk (~10x)",
            "CT": "Heterozygous HLA-DQ2.5 — moderate celiac risk",
            "CC": "Low celiac risk (from DQ2)",
        },
        "clinical_significance": "pathogenic_risk",
        "source": "van Heel et al. 2007",
    },
    "rs7454108": {
        "gene": "HLA-DQ8",
        "risk_allele": "C",
        "category": "autoimmune",
        "condition": "Celiac disease",
        "effect": {
            "CC": "HLA-DQ8 positive — celiac risk",
            "CT": "HLA-DQ8 carrier",
            "TT": "HLA-DQ8 negative",
        },
        "clinical_significance": "risk_factor",
        "source": "HLA typing studies",
    },
    "rs10484554": {
        "gene": "HLA-C",
        "risk_allele": "T",
        "category": "autoimmune",
        "condition": "Psoriasis",
        "effect": {
            "TT": "Increased psoriasis risk (HLA-Cw6)",
            "CT": "Moderately increased risk",
            "CC": "Average risk",
        },
        "clinical_significance": "risk_factor",
        "source": "de Cid et al. 2009",
    },

    # === NEUROLOGICAL / MENTAL HEALTH ===
    "rs429358": {
        "gene": "APOE (determines E4 vs E3)",
        "risk_allele": "C",
        "category": "neurological",
        "condition": "Alzheimer's disease",
        "effect": {
            "CC": "APOE4/E4 (with rs7412 TT) — ~12x Alzheimer's risk",
            "CT": "One APOE4 allele — ~3x risk (check rs7412)",
            "TT": "No APOE4 allele (E3 or E2 depending on rs7412)",
        },
        "clinical_significance": "pathogenic_risk",
        "source": "Corder et al. Science 1993",
    },
    "rs7412": {
        "gene": "APOE (determines E2 vs E3)",
        "risk_allele": "T",
        "category": "neurological",
        "condition": "Alzheimer's disease (protective if APOE2)",
        "effect": {
            "TT": "APOE2/E2 — protective against Alzheimer's (if rs429358=TT)",
            "CT": "One APOE2 allele — somewhat protective",
            "CC": "No APOE2 (E3 or E4 depending on rs429358)",
        },
        "clinical_significance": "protective/risk",
        "source": "SNPedia",
    },
    "rs6265": {
        "gene": "BDNF (Val66Met)",
        "risk_allele": "T",
        "category": "neurological",
        "condition": "Depression, anxiety, reduced neuroplasticity",
        "effect": {
            "TT": "Met/Met — reduced BDNF secretion, higher depression/anxiety risk",
            "CT": "Val/Met — slightly reduced BDNF",
            "CC": "Val/Val — normal BDNF secretion",
        },
        "clinical_significance": "risk_factor",
        "source": "Egan et al. Cell 2003",
    },
    "rs4680": {
        "gene": "COMT (Val158Met)",
        "risk_allele": "A",
        "category": "neurological",
        "condition": "Stress response, pain sensitivity, dopamine metabolism",
        "effect": {
            "AA": "Met/Met — slow COMT, higher dopamine, more stress-sensitive ('worrier')",
            "AG": "Val/Met — intermediate COMT activity",
            "GG": "Val/Val — fast COMT, lower dopamine, stress-resilient ('warrior')",
        },
        "clinical_significance": "risk_factor",
        "source": "Zubieta et al. 2003",
    },
    "rs1800497": {
        "gene": "ANKK1/DRD2 (Taq1A)",
        "risk_allele": "T",
        "category": "neurological",
        "condition": "Dopamine receptor density, addiction risk",
        "effect": {
            "TT": "A1/A1 — fewer D2 receptors, higher addiction risk",
            "CT": "A1/A2 — moderately reduced D2 density",
            "CC": "A2/A2 — normal D2 receptor density",
        },
        "clinical_significance": "risk_factor",
        "source": "Noble et al. 1991",
    },
    "rs25531": {
        "gene": "SLC6A4 (5-HTTLPR)",
        "risk_allele": "A",
        "category": "neurological",
        "condition": "Serotonin transporter expression, SSRI response",
        "effect": {
            "AA": "Lower serotonin reuptake — may affect SSRI response",
            "AG": "Intermediate",
            "GG": "Higher serotonin reuptake",
        },
        "clinical_significance": "risk_factor",
        "source": "Hu et al. 2006",
    },
    "rs4570625": {
        "gene": "TPH2",
        "risk_allele": "G",
        "category": "neurological",
        "condition": "Serotonin synthesis, depression risk",
        "effect": {
            "GG": "Altered serotonin synthesis — depression risk",
            "GT": "Intermediate",
            "TT": "Normal serotonin synthesis",
        },
        "clinical_significance": "risk_factor",
        "source": "Zhang et al. 2005",
    },
    "rs1799971": {
        "gene": "OPRM1 (A118G)",
        "risk_allele": "G",
        "category": "neurological",
        "condition": "Pain sensitivity, opioid response, addiction risk",
        "effect": {
            "GG": "Reduced opioid receptor binding — higher pain sensitivity, altered addiction risk",
            "AG": "Intermediate opioid sensitivity",
            "AA": "Normal opioid receptor function",
        },
        "clinical_significance": "risk_factor",
        "source": "Bond et al. PNAS 1998",
    },
}


# --- NUTRIGENOMICS SNPs ---
NUTRIGENOMICS_SNPS = {
    # === CAFFEINE ===
    "rs762551": {
        "gene": "CYP1A2",
        "risk_allele": "C",
        "category": "caffeine_metabolism",
        "condition": "Caffeine metabolism speed",
        "effect": {
            "AA": "Fast metabolizer — caffeine cleared quickly, lower CV risk from coffee",
            "AC": "Slow metabolizer — caffeine stays longer, may increase CV risk with >2 cups/day",
            "CC": "Slow metabolizer — significantly slower clearance, limit caffeine",
        },
        "recommendation": {
            "AA": "Moderate coffee intake (3-4 cups/day) may be protective",
            "AC": "Limit to 1-2 cups/day, avoid afternoon caffeine",
            "CC": "Limit to 1 cup/day or less, consider decaf",
        },
        "source": "Cornelis et al. JAMA 2006",
    },

    # === LACTOSE ===
    "rs4988235": {
        "gene": "MCM6/LCT",
        "risk_allele": "T",
        "category": "lactose",
        "condition": "Lactose tolerance/intolerance",
        "effect": {
            "TT": "Lactose tolerant — produces lactase into adulthood",
            "CT": "Likely lactose tolerant",
            "CC": "Lactose intolerant — lactase production decreases after childhood",
        },
        "recommendation": {
            "TT": "Dairy is well tolerated",
            "CT": "Dairy is likely tolerated, monitor symptoms",
            "CC": "Consider lactose-free dairy, fermented dairy (yogurt, kefir), or calcium supplements",
        },
        "source": "Enattah et al. 2002",
    },

    # === ALCOHOL ===
    "rs671": {
        "gene": "ALDH2",
        "risk_allele": "A",
        "category": "alcohol",
        "condition": "Alcohol metabolism (acetaldehyde detox)",
        "effect": {
            "AA": "Non-functional ALDH2 — severe alcohol intolerance, high esophageal cancer risk",
            "AG": "Reduced ALDH2 — alcohol flush, increased cancer risk from drinking",
            "GG": "Normal ALDH2 — standard alcohol metabolism",
        },
        "recommendation": {
            "AA": "Avoid alcohol completely — high acetaldehyde accumulation",
            "AG": "Minimize alcohol — increased esophageal cancer risk",
            "GG": "Standard alcohol guidelines apply",
        },
        "source": "Brooks et al. 2009",
    },
    "rs1229984": {
        "gene": "ADH1B",
        "risk_allele": "T",
        "category": "alcohol",
        "condition": "Alcohol metabolism speed",
        "effect": {
            "TT": "Ultra-rapid alcohol → acetaldehyde conversion (protective against alcoholism)",
            "CT": "Fast conversion",
            "CC": "Normal conversion speed",
        },
        "source": "Li et al. 2011",
    },

    # === VITAMIN D ===
    "rs2228570": {
        "gene": "VDR (FokI)",
        "risk_allele": "T",
        "category": "vitamin_d",
        "condition": "Vitamin D receptor function",
        "effect": {
            "TT": "Less active VDR — may need more vitamin D",
            "CT": "Intermediate VDR activity",
            "CC": "More active VDR — better vitamin D utilization",
        },
        "recommendation": {
            "TT": "Higher vitamin D supplementation may be needed (4000-5000 IU/day)",
            "CT": "Moderate supplementation (2000-4000 IU/day)",
            "CC": "Standard supplementation (1000-2000 IU/day)",
        },
        "source": "Uitterlinden et al. 2004",
    },
    "rs10741657": {
        "gene": "CYP2R1",
        "risk_allele": "A",
        "category": "vitamin_d",
        "condition": "Vitamin D 25-hydroxylation",
        "effect": {
            "AA": "Reduced conversion of vitamin D to 25(OH)D",
            "AG": "Slightly reduced conversion",
            "GG": "Normal vitamin D activation",
        },
        "source": "Wang et al. Lancet 2010",
    },
    "rs12794714": {
        "gene": "CYP2R1",
        "risk_allele": "A",
        "category": "vitamin_d",
        "condition": "Vitamin D levels",
        "effect": {
            "AA": "Lower circulating 25(OH)D levels",
            "AG": "Slightly lower levels",
            "GG": "Normal levels",
        },
        "source": "Ahn et al. 2010",
    },
    "rs2282679": {
        "gene": "GC (vitamin D binding protein)",
        "risk_allele": "C",
        "category": "vitamin_d",
        "condition": "Vitamin D transport",
        "effect": {
            "CC": "Lower vitamin D binding protein — lower total 25(OH)D",
            "AC": "Intermediate levels",
            "AA": "Normal vitamin D transport",
        },
        "source": "Wang et al. Lancet 2010",
    },

    # === FOLATE / MTHFR ===
    "rs1801133": {
        "gene": "MTHFR (C677T)",
        "risk_allele": "A",
        "category": "folate",
        "condition": "Folate metabolism, homocysteine levels",
        "effect": {
            "AA": "~70% reduced MTHFR activity — elevated homocysteine, impaired folate metabolism",
            "AG": "~35% reduced MTHFR activity",
            "GG": "Normal MTHFR activity",
        },
        "recommendation": {
            "AA": "Take methylfolate (L-5-MTHF) instead of folic acid, 800-1000 mcg/day. Monitor homocysteine.",
            "AG": "Consider methylfolate 400-800 mcg/day",
            "GG": "Standard folate intake sufficient",
        },
        "source": "Frosst et al. 1995, ClinVar",
    },
    "rs1801131": {
        "gene": "MTHFR (A1298C)",
        "risk_allele": "G",
        "category": "folate",
        "condition": "Folate metabolism (complementary to C677T)",
        "effect": {
            "GG": "Reduced MTHFR activity (compounds with C677T)",
            "GT": "Slightly reduced MTHFR activity",
            "TT": "Normal for this variant",
        },
        "recommendation": {
            "GG": "Use methylfolate; especially important if also C677T heterozygous",
            "GT": "Consider methylfolate if combined with C677T variant",
            "TT": "No action needed for this variant alone",
        },
        "source": "van der Put et al. 1998",
    },

    # === VITAMIN B12 ===
    "rs602662": {
        "gene": "FUT2",
        "risk_allele": "A",
        "category": "vitamin_b12",
        "condition": "Vitamin B12 absorption",
        "effect": {
            "AA": "Non-secretor — lower B12 absorption from food",
            "AG": "Intermediate B12 absorption",
            "GG": "Secretor — normal B12 absorption",
        },
        "recommendation": {
            "AA": "Monitor B12 levels, consider methylcobalamin supplementation",
            "AG": "Monitor B12 levels periodically",
            "GG": "Standard B12 intake usually sufficient",
        },
        "source": "Hazra et al. 2008",
    },
    "rs1801198": {
        "gene": "TCN2 (transcobalamin)",
        "risk_allele": "G",
        "category": "vitamin_b12",
        "condition": "Vitamin B12 cellular delivery",
        "effect": {
            "GG": "Reduced B12 delivery to cells",
            "CG": "Intermediate B12 transport",
            "CC": "Normal B12 cellular delivery",
        },
        "source": "Hazra et al. 2009",
    },

    # === VITAMIN A ===
    "rs7501331": {
        "gene": "BCMO1",
        "risk_allele": "T",
        "category": "vitamin_a",
        "condition": "Beta-carotene to vitamin A conversion",
        "effect": {
            "TT": "~50% reduced conversion — poor converter of beta-carotene",
            "CT": "~30% reduced conversion",
            "CC": "Normal conversion",
        },
        "recommendation": {
            "TT": "Get preformed vitamin A (retinol) from animal sources or supplements",
            "CT": "Include some preformed vitamin A sources",
            "CC": "Beta-carotene from plants is efficiently converted",
        },
        "source": "Leung et al. 2009",
    },
    "rs12934922": {
        "gene": "BCMO1",
        "risk_allele": "T",
        "category": "vitamin_a",
        "condition": "Beta-carotene conversion (secondary variant)",
        "effect": {
            "TT": "Further reduced beta-carotene conversion",
            "AT": "Slightly reduced conversion",
            "AA": "Normal conversion",
        },
        "source": "Leung et al. 2009",
    },

    # === OMEGA-3 / FAT METABOLISM ===
    "rs174547": {
        "gene": "FADS1",
        "risk_allele": "C",
        "category": "fat_metabolism",
        "condition": "Omega-3/6 fatty acid conversion",
        "effect": {
            "CC": "Reduced conversion of ALA to EPA/DHA — need direct omega-3 sources",
            "CT": "Intermediate conversion",
            "TT": "Good conversion of plant omega-3 to EPA/DHA",
        },
        "recommendation": {
            "CC": "Supplement with fish oil or algal DHA/EPA directly",
            "CT": "Include fish 2-3x/week or consider fish oil",
            "TT": "Plant-based omega-3 (flax, chia, walnuts) may be sufficient",
        },
        "source": "Tanaka et al. 2009",
    },

    # === GLUTEN SENSITIVITY ===
    # rs2187668 and rs7454108 already in autoimmune section — cross-referenced

    # === SALT SENSITIVITY ===
    "rs4961": {
        "gene": "ADD1 (Gly460Trp)",
        "risk_allele": "T",
        "category": "salt_sensitivity",
        "condition": "Salt-sensitive hypertension",
        "effect": {
            "TT": "Salt-sensitive — higher BP response to sodium",
            "GT": "Moderately salt-sensitive",
            "GG": "Normal salt sensitivity",
        },
        "recommendation": {
            "TT": "Limit sodium to <1500mg/day, increase potassium",
            "GT": "Moderate sodium intake <2000mg/day",
            "GG": "Standard sodium guidelines apply (<2300mg/day)",
        },
        "source": "Cusi et al. 1997",
    },
    "rs1799998": {
        "gene": "CYP11B2 (aldosterone synthase)",
        "risk_allele": "T",
        "category": "salt_sensitivity",
        "condition": "Aldosterone production / blood pressure",
        "effect": {
            "TT": "Higher aldosterone — may increase sodium retention and BP",
            "CT": "Intermediate",
            "CC": "Normal aldosterone regulation",
        },
        "source": "Davies et al. 1999",
    },

    # === CHOLINE ===
    "rs12325817": {
        "gene": "PEMT",
        "risk_allele": "C",
        "category": "choline",
        "condition": "Choline synthesis",
        "effect": {
            "CC": "Reduced endogenous choline production — dietary choline essential",
            "CG": "Slightly reduced choline synthesis",
            "GG": "Normal choline synthesis",
        },
        "recommendation": {
            "CC": "Increase dietary choline (eggs, liver) or supplement",
            "CG": "Ensure adequate dietary choline",
            "GG": "Standard diet usually provides sufficient choline",
        },
        "source": "da Costa et al. 2006",
    },
}


# --- PHARMACOGENOMICS SNPs ---
PHARMACOGENOMICS_SNPS = {
    # === CYP2D6 ===
    "rs3892097": {
        "gene": "CYP2D6 (*4)",
        "risk_allele": "A",
        "category": "cyp2d6",
        "condition": "CYP2D6 poor metabolizer",
        "effect": {
            "AA": "Poor metabolizer — non-functional CYP2D6",
            "GA": "Intermediate metabolizer — one functional copy",
            "GG": "Extensive (normal) metabolizer",
        },
        "drugs_affected": [
            "codeine (no conversion to morphine — ineffective)",
            "tramadol (reduced efficacy)",
            "tamoxifen (reduced activation to endoxifen)",
            "atomoxetine (increased side effects)",
            "many antidepressants (SSRIs, TCAs)",
            "metoprolol (increased levels)",
            "ondansetron (reduced efficacy)",
        ],
        "source": "PharmGKB, CPIC guidelines",
    },
    "rs16947": {
        "gene": "CYP2D6 (*2)",
        "risk_allele": "A",
        "category": "cyp2d6",
        "condition": "CYP2D6 activity variant",
        "effect": {
            "AA": "CYP2D6*2 homozygous — normal to increased activity",
            "GA": "Heterozygous",
            "GG": "Reference allele",
        },
        "source": "PharmGKB",
    },
    "rs1065852": {
        "gene": "CYP2D6 (*10)",
        "risk_allele": "A",
        "category": "cyp2d6",
        "condition": "CYP2D6 reduced function",
        "effect": {
            "AA": "Reduced CYP2D6 activity — *10 homozygous",
            "GA": "Intermediate activity",
            "GG": "Normal activity",
        },
        "source": "PharmGKB",
    },

    # === CYP2C19 ===
    "rs4244285": {
        "gene": "CYP2C19 (*2)",
        "risk_allele": "A",
        "category": "cyp2c19",
        "condition": "CYP2C19 poor metabolizer",
        "effect": {
            "AA": "Poor metabolizer — no CYP2C19 activity from this allele",
            "GA": "Intermediate metabolizer",
            "GG": "Extensive (normal) metabolizer",
        },
        "drugs_affected": [
            "clopidogrel (Plavix) — CRITICAL: reduced activation, increased CV events",
            "omeprazole/esomeprazole (increased efficacy)",
            "citalopram/escitalopram (increased levels)",
            "voriconazole (increased levels)",
            "diazepam (slower clearance)",
        ],
        "source": "PharmGKB, CPIC, FDA label",
    },
    "rs4986893": {
        "gene": "CYP2C19 (*3)",
        "risk_allele": "A",
        "category": "cyp2c19",
        "condition": "CYP2C19 non-functional",
        "effect": {
            "AA": "Non-functional CYP2C19",
            "GA": "One non-functional allele",
            "GG": "Normal for this variant",
        },
        "source": "PharmGKB",
    },
    "rs12248560": {
        "gene": "CYP2C19 (*17)",
        "risk_allele": "T",
        "category": "cyp2c19",
        "condition": "CYP2C19 ultra-rapid metabolizer",
        "effect": {
            "TT": "Ultra-rapid metabolizer — faster drug clearance",
            "CT": "Rapid metabolizer",
            "CC": "Normal metabolizer for this variant",
        },
        "drugs_affected": [
            "clopidogrel (increased activation — good)",
            "PPIs (reduced efficacy — faster clearance)",
            "escitalopram (lower levels — may need higher dose)",
            "voriconazole (reduced efficacy)",
        ],
        "source": "PharmGKB, CPIC",
    },

    # === CYP2C9 ===
    "rs1799853": {
        "gene": "CYP2C9 (*2)",
        "risk_allele": "T",
        "category": "cyp2c9",
        "condition": "CYP2C9 reduced function",
        "effect": {
            "TT": "Significantly reduced CYP2C9 — *2 homozygous",
            "CT": "Intermediate metabolizer",
            "CC": "Normal metabolizer",
        },
        "drugs_affected": [
            "warfarin (CRITICAL: lower dose needed, bleeding risk)",
            "NSAIDs (ibuprofen, celecoxib — slower clearance)",
            "phenytoin (increased levels)",
            "losartan (reduced activation)",
        ],
        "source": "PharmGKB, CPIC, FDA label",
    },
    "rs1057910": {
        "gene": "CYP2C9 (*3)",
        "risk_allele": "C",
        "category": "cyp2c9",
        "condition": "CYP2C9 poor function",
        "effect": {
            "CC": "Poor CYP2C9 function — *3 homozygous",
            "AC": "Intermediate metabolizer",
            "AA": "Normal metabolizer",
        },
        "drugs_affected": [
            "warfarin (CRITICAL: significantly lower dose needed)",
            "NSAIDs (accumulation risk)",
            "sulfonylureas (hypoglycemia risk)",
        ],
        "source": "PharmGKB, CPIC",
    },

    # === CYP3A4/5 ===
    "rs35599367": {
        "gene": "CYP3A4 (*22)",
        "risk_allele": "T",
        "category": "cyp3a4",
        "condition": "CYP3A4 reduced expression",
        "effect": {
            "TT": "Significantly reduced CYP3A4 activity",
            "CT": "Intermediate CYP3A4 activity",
            "CC": "Normal CYP3A4 activity",
        },
        "drugs_affected": [
            "statins (atorvastatin, simvastatin — increased levels)",
            "tacrolimus (lower dose needed)",
            "cyclosporine (increased levels)",
            "many chemotherapy drugs",
        ],
        "source": "PharmGKB",
    },
    "rs776746": {
        "gene": "CYP3A5 (*3)",
        "risk_allele": "C",
        "category": "cyp3a5",
        "condition": "CYP3A5 non-expressor",
        "effect": {
            "CC": "CYP3A5 non-expressor (*3/*3) — common in Europeans",
            "CT": "Intermediate expressor",
            "TT": "CYP3A5 expressor (*1/*1)",
        },
        "drugs_affected": [
            "tacrolimus (non-expressors need lower dose)",
        ],
        "source": "PharmGKB, CPIC",
    },

    # === WARFARIN SENSITIVITY ===
    "rs9923231": {
        "gene": "VKORC1",
        "risk_allele": "T",
        "category": "warfarin",
        "condition": "Warfarin sensitivity",
        "effect": {
            "TT": "Highly sensitive — need ~50% lower warfarin dose",
            "CT": "Moderately sensitive — need ~25% lower dose",
            "CC": "Normal sensitivity — standard dosing",
        },
        "drugs_affected": ["warfarin (CRITICAL: dose adjustment required)"],
        "source": "CPIC, FDA label",
    },

    # === STATIN MYOPATHY ===
    "rs4149056": {
        "gene": "SLCO1B1 (*5)",
        "risk_allele": "C",
        "category": "statins",
        "condition": "Statin-induced myopathy risk",
        "effect": {
            "CC": "High risk of simvastatin myopathy (~17x risk)",
            "CT": "Moderate risk (~4x risk)",
            "TT": "Normal risk",
        },
        "drugs_affected": [
            "simvastatin (CRITICAL: consider alternative statin)",
            "atorvastatin (moderate risk)",
            "rosuvastatin (lower risk but monitor)",
        ],
        "source": "CPIC, SEARCH Collaborative 2008",
    },

    # === 5-FU TOXICITY ===
    "rs3918290": {
        "gene": "DPYD (*2A)",
        "risk_allele": "A",
        "category": "chemotherapy",
        "condition": "5-fluorouracil/capecitabine toxicity",
        "effect": {
            "AA": "CRITICAL: DPD deficient — life-threatening 5-FU toxicity",
            "GA": "Intermediate DPD — reduce 5-FU dose by 50%",
            "GG": "Normal DPD activity",
        },
        "drugs_affected": [
            "5-fluorouracil (CRITICAL: life-threatening toxicity if DPD deficient)",
            "capecitabine (same risk)",
        ],
        "source": "CPIC, FDA label",
    },

    # === METHOTREXATE ===
    # rs1801133 (MTHFR) also affects methotrexate — cross-referenced from folate

    # === OPIOID RESPONSE ===
    # rs1799971 (OPRM1) in neurological section — also affects opioid dosing

    # === THIOPURINES ===
    "rs1800462": {
        "gene": "TPMT (*2)",
        "risk_allele": "G",
        "category": "thiopurines",
        "condition": "Thiopurine toxicity",
        "effect": {
            "GG": "TPMT deficient — high risk of myelosuppression",
            "AG": "Intermediate TPMT",
            "AA": "Normal TPMT",
        },
        "drugs_affected": [
            "azathioprine (dose reduction needed)",
            "6-mercaptopurine (dose reduction needed)",
            "thioguanine",
        ],
        "source": "CPIC",
    },
    "rs1142345": {
        "gene": "TPMT (*3C)",
        "risk_allele": "C",
        "category": "thiopurines",
        "condition": "Thiopurine metabolism",
        "effect": {
            "CC": "TPMT deficient",
            "AC": "Intermediate TPMT",
            "AA": "Normal TPMT",
        },
        "source": "CPIC",
    },
}


# --- FITNESS & RECOVERY SNPs ---
FITNESS_SNPS = {
    "rs1815739": {
        "gene": "ACTN3 (R577X)",
        "risk_allele": "T",
        "category": "muscle_fiber",
        "condition": "Muscle fiber type composition",
        "effect": {
            "TT": "XX — no alpha-actinin-3 in fast-twitch fibers, endurance-oriented",
            "CT": "RX — mixed fiber type, good all-around",
            "CC": "RR — full alpha-actinin-3, power/sprint-oriented",
        },
        "recommendation": {
            "TT": "Focus on endurance training (running, cycling, swimming). May excel at marathon distance.",
            "CT": "Versatile — can develop both power and endurance effectively",
            "CC": "May have advantage in power/sprint events. Include explosive training.",
        },
        "source": "Yang et al. AJHG 2003",
    },
    "rs4341": {
        "gene": "ACE (I/D proxy)",
        "risk_allele": "G",
        "category": "endurance",
        "condition": "ACE activity / endurance capacity",
        "effect": {
            "GG": "DD-like — higher ACE, power/strength oriented",
            "CG": "ID-like — mixed endurance/power",
            "CC": "II-like — lower ACE, endurance oriented",
        },
        "recommendation": {
            "GG": "Strength/power training may yield better results",
            "CG": "Balanced training approach recommended",
            "CC": "Endurance training may yield better results, good VO2max response",
        },
        "source": "Rigat et al. 1990, Montgomery et al. 1998",
    },
    "rs12722": {
        "gene": "COL5A1",
        "risk_allele": "T",
        "category": "injury_risk",
        "condition": "Tendon/ligament injury susceptibility",
        "effect": {
            "TT": "Increased risk of tendon/ligament injuries (ACL, Achilles)",
            "CT": "Moderate injury risk",
            "CC": "Lower injury susceptibility — more flexible connective tissue",
        },
        "recommendation": {
            "TT": "Prioritize warm-up, flexibility work, and gradual load increases. Consider collagen supplementation.",
            "CT": "Standard injury prevention measures",
            "CC": "Normal injury risk — maintain standard prevention",
        },
        "source": "Posthumus et al. 2009",
    },
    "rs8192678": {
        "gene": "PPARGC1A (Gly482Ser)",
        "risk_allele": "A",
        "category": "endurance",
        "condition": "Mitochondrial biogenesis / VO2max",
        "effect": {
            "AA": "Ser/Ser — reduced PGC-1α activity, lower aerobic capacity baseline",
            "GA": "Intermediate",
            "GG": "Gly/Gly — better mitochondrial biogenesis, higher VO2max response",
        },
        "recommendation": {
            "AA": "May need more training volume for aerobic gains. Focus on Zone 2 training.",
            "GA": "Standard endurance training approach",
            "GG": "Good responder to endurance training",
        },
        "source": "Ling et al. 2004",
    },
    "rs1049434": {
        "gene": "MCT1 (SLC16A1)",
        "risk_allele": "T",
        "category": "recovery",
        "condition": "Lactate clearance",
        "effect": {
            "TT": "Reduced lactate transport — slower clearance during exercise",
            "AT": "Intermediate lactate handling",
            "AA": "Normal lactate transport — efficient clearance",
        },
        "recommendation": {
            "TT": "Allow longer recovery between high-intensity intervals. Lactate threshold training important.",
            "AT": "Standard interval recovery",
            "AA": "Good lactate clearance — can handle shorter rest intervals",
        },
        "source": "Merezhinskaya et al. 2000",
    },
    "rs7460": {
        "gene": "COL1A1 (Sp1 binding site)",
        "risk_allele": "T",
        "category": "injury_risk",
        "condition": "Bone density / stress fracture risk",
        "effect": {
            "TT": "Reduced bone mineral density — higher stress fracture risk",
            "GT": "Intermediate bone density",
            "GG": "Normal bone mineral density",
        },
        "recommendation": {
            "TT": "Ensure adequate calcium (1200mg) and vitamin D. Gradual load progression.",
            "GT": "Standard bone health measures",
            "GG": "Normal fracture risk",
        },
        "source": "Mann et al. 2001",
    },
    "rs2228570_fitness": {
        "gene": "VDR (FokI) — fitness context",
        "risk_allele": "T",
        "category": "recovery",
        "condition": "Bone health and muscle recovery",
        "effect": {
            "TT": "Poorer vitamin D utilization affecting bone/muscle recovery",
            "CT": "Intermediate",
            "CC": "Good vitamin D utilization for recovery",
        },
        "source": "Uitterlinden et al. 2004",
    },
    "rs4880_fitness": {
        "gene": "SOD2 (fitness context)",
        "risk_allele": "T",
        "category": "recovery",
        "condition": "Exercise-induced oxidative stress",
        "effect": {
            "TT": "Ala/Ala — more SOD2 in mitochondria, but may overproduce H2O2 with exercise",
            "CT": "Val/Ala — intermediate",
            "CC": "Val/Val — less efficient mitochondrial SOD2, more oxidative stress",
        },
        "recommendation": {
            "TT": "May benefit from moderate antioxidant intake (not excessive)",
            "CT": "Standard antioxidant intake",
            "CC": "Higher antioxidant needs — vitamin C, E, selenium, CoQ10 post-exercise",
        },
        "source": "Shimoda-Matsubayashi et al. 1996",
    },
}


# --- METHYLATION & DETOX SNPs ---
METHYLATION_DETOX_SNPS = {
    # MTHFR already in NUTRIGENOMICS — cross-referenced
    "rs1801133_methyl": {
        "gene": "MTHFR (C677T)",
        "risk_allele": "A",
        "category": "methylation",
        "condition": "Methylation cycle — rate-limiting step",
        "effect": {
            "AA": "~70% reduced — impaired methylation, elevated homocysteine",
            "AG": "~35% reduced — mildly impaired methylation",
            "GG": "Normal methylation capacity",
        },
        "supplement_recommendation": {
            "AA": "Methylfolate (L-5-MTHF) 800-1000mcg, methylcobalamin (B12) 1000mcg, "
                  "riboflavin (B2) 25-50mg, TMG (betaine) 500-1000mg",
            "AG": "Methylfolate 400-800mcg, methylcobalamin 500-1000mcg, B2 25mg",
            "GG": "Standard B-complex sufficient",
        },
        "source": "Frosst et al. 1995",
    },
    "rs1801131_methyl": {
        "gene": "MTHFR (A1298C)",
        "risk_allele": "G",
        "category": "methylation",
        "condition": "BH4 production (neurotransmitter cofactor)",
        "effect": {
            "GG": "Reduced BH4 — may affect serotonin/dopamine synthesis",
            "GT": "Mildly reduced BH4",
            "TT": "Normal BH4 production",
        },
        "supplement_recommendation": {
            "GG": "Methylfolate, SAMe (if not bipolar), consider BH4 cofactors",
            "GT": "Methylfolate may help",
            "TT": "No specific action needed",
        },
        "source": "van der Put et al. 1998",
    },
    "rs4680_methyl": {
        "gene": "COMT (Val158Met)",
        "risk_allele": "A",
        "category": "methylation",
        "condition": "Catechol-O-methyltransferase — dopamine/estrogen clearance",
        "effect": {
            "AA": "Met/Met — slow COMT: higher dopamine, estrogen, catecholamines. "
                  "Better focus but more stress-sensitive. SAMe may overstimulate.",
            "AG": "Val/Met — intermediate COMT",
            "GG": "Val/Val — fast COMT: rapid dopamine clearance. More stress-resilient "
                  "but may need dopamine support. SAMe well-tolerated.",
        },
        "supplement_recommendation": {
            "AA": "AVOID high-dose SAMe, methyl donors may worsen anxiety. "
                  "Support with magnesium glycinate, NAC. Avoid excess catechol-containing foods.",
            "AG": "Moderate methyl donor support. Magnesium helpful.",
            "GG": "Can tolerate SAMe, methyl donors. May benefit from dopamine precursors "
                  "(tyrosine, phenylalanine) if needed.",
        },
        "source": "Lotta et al. 1995",
    },
    "rs234706": {
        "gene": "CBS (C699T)",
        "risk_allele": "A",
        "category": "transsulfuration",
        "condition": "Cystathionine beta-synthase — homocysteine to cysteine",
        "effect": {
            "AA": "Upregulated CBS — faster homocysteine → taurine/sulfate. "
                  "May deplete homocysteine too fast, affecting methylation.",
            "AG": "Mildly upregulated CBS",
            "GG": "Normal CBS activity",
        },
        "supplement_recommendation": {
            "AA": "Limit sulfur-containing supplements (NAC, MSM, alpha-lipoic acid). "
                  "Support methylation with methylfolate and B12. Molybdenum 150-300mcg "
                  "to help process excess sulfites.",
            "AG": "Monitor sulfur tolerance",
            "GG": "Standard supplement protocol",
        },
        "source": "Kraus et al. 1999",
    },
    "rs4880": {
        "gene": "SOD2 (Ala16Val)",
        "risk_allele": "T",
        "category": "antioxidant_defense",
        "condition": "Mitochondrial superoxide dismutase",
        "effect": {
            "TT": "Ala/Ala — SOD2 efficiently imported to mitochondria. "
                  "More H2O2 produced → need adequate catalase/GPx.",
            "CT": "Val/Ala — intermediate",
            "CC": "Val/Val — less efficient SOD2 import. More superoxide damage. "
                  "Higher oxidative stress.",
        },
        "supplement_recommendation": {
            "TT": "Manganese (cofactor), avoid excessive iron/copper (Fenton reaction). "
                  "Glutathione support (NAC, selenium for GPx).",
            "CT": "Balanced antioxidant approach",
            "CC": "CoQ10 200-400mg, MitoQ, NAC 600-1200mg, selenium 200mcg, "
                  "vitamin C, vitamin E (mixed tocopherols). High antioxidant diet.",
        },
        "source": "Shimoda-Matsubayashi et al. 1996",
    },
    "rs1695": {
        "gene": "GSTP1 (Ile105Val)",
        "risk_allele": "G",
        "category": "detoxification",
        "condition": "Glutathione S-transferase Pi — Phase II detox",
        "effect": {
            "GG": "Val/Val — reduced GST-Pi activity, impaired detoxification",
            "AG": "Ile/Val — moderately reduced detox capacity",
            "AA": "Ile/Ile — normal GST-Pi detox function",
        },
        "supplement_recommendation": {
            "GG": "NAC 600-1200mg, liposomal glutathione, cruciferous vegetables "
                  "(broccoli sprouts for sulforaphane), milk thistle, selenium",
            "AG": "Support glutathione: NAC, cruciferous vegetables, selenium",
            "AA": "Standard detox support sufficient",
        },
        "source": "Watson et al. 1998",
    },
    "rs1800566_methyl": {
        "gene": "NQO1 (Pro187Ser)",
        "risk_allele": "T",
        "category": "detoxification",
        "condition": "NAD(P)H quinone dehydrogenase — benzene/quinone detox",
        "effect": {
            "TT": "No NQO1 activity — cannot detoxify quinones, increased oxidative damage",
            "CT": "Reduced NQO1 activity",
            "CC": "Normal NQO1 activity",
        },
        "supplement_recommendation": {
            "TT": "CoQ10 (ubiquinol form) 200-400mg essential. Avoid benzene exposure. "
                  "Riboflavin (B2) 25-50mg as NQO1 cofactor.",
            "CT": "CoQ10 100-200mg, riboflavin 25mg",
            "CC": "Standard CoQ10 supplementation optional",
        },
        "source": "Ross et al. 2000",
    },
    "rs1801280": {
        "gene": "NAT2 (slow acetylator)",
        "risk_allele": "A",
        "category": "detoxification",
        "condition": "N-acetyltransferase 2 — Phase II acetylation",
        "effect": {
            "AA": "Slow acetylator — slower processing of aromatic amines, some drugs",
            "AG": "Intermediate acetylator",
            "GG": "Rapid acetylator",
        },
        "supplement_recommendation": {
            "AA": "Reduce exposure to heterocyclic amines (charred meat). "
                  "B-vitamins support acetylation. Cruciferous vegetables helpful.",
            "AG": "Moderate charred food intake",
            "GG": "Normal detox capacity for amines",
        },
        "source": "Hein et al. 2000",
    },
    "rs1056806": {
        "gene": "CBS (N212N)",
        "risk_allele": "T",
        "category": "transsulfuration",
        "condition": "CBS expression regulation",
        "effect": {
            "TT": "Possible CBS upregulation",
            "CT": "Intermediate",
            "CC": "Normal expression",
        },
        "source": "Kraus et al. 1999",
    },
    "rs2066470": {
        "gene": "MTRR (Ile22Met)",
        "risk_allele": "A",
        "category": "methylation",
        "condition": "Methionine synthase reductase — B12 recycling",
        "effect": {
            "AA": "Reduced MTRR — impaired B12 recycling for methylation",
            "AG": "Mildly reduced",
            "GG": "Normal MTRR function",
        },
        "supplement_recommendation": {
            "AA": "Methylcobalamin or hydroxocobalamin (not cyanocobalamin). "
                  "Riboflavin as MTRR cofactor.",
            "AG": "Consider active B12 forms",
            "GG": "Standard B12 supplementation",
        },
        "source": "Wilson et al. 1999",
    },
    "rs1805087": {
        "gene": "MTR (A2756G)",
        "risk_allele": "G",
        "category": "methylation",
        "condition": "Methionine synthase — homocysteine remethylation",
        "effect": {
            "GG": "Increased MTR activity — may deplete methyl-B12 faster",
            "AG": "Mildly increased activity",
            "AA": "Normal MTR activity",
        },
        "supplement_recommendation": {
            "GG": "Higher methylcobalamin needs. Methylfolate support.",
            "AG": "Standard methylcobalamin",
            "AA": "Standard B12",
        },
        "source": "Leclerc et al. 1998",
    },
    "rs7946": {
        "gene": "PEMT (G5465A)",
        "risk_allele": "A",
        "category": "methylation",
        "condition": "Phosphatidylcholine synthesis via methylation",
        "effect": {
            "AA": "Reduced PEMT — lower endogenous choline/phosphatidylcholine synthesis. "
                  "Uses more SAMe for choline → affects methylation.",
            "GA": "Mildly reduced",
            "GG": "Normal PEMT activity",
        },
        "supplement_recommendation": {
            "AA": "Dietary choline essential (eggs, liver). Phosphatidylcholine supplement. "
                  "Betaine (TMG) to spare methyl groups.",
            "GA": "Ensure adequate dietary choline",
            "GG": "Standard diet usually sufficient",
        },
        "source": "da Costa et al. 2006",
    },
}


def get_all_tracked_rsids():
    """Return a set of all rsIDs tracked across all databases."""
    all_rsids = set()
    for db in [HEALTH_RISK_SNPS, NUTRIGENOMICS_SNPS, PHARMACOGENOMICS_SNPS,
               FITNESS_SNPS, METHYLATION_DETOX_SNPS]:
        for rsid in db:
            # Strip suffixes like _methyl, _fitness for actual rsID lookup
            clean_rsid = rsid.split("_")[0] if "_" in rsid and rsid.split("_")[0].startswith("rs") else rsid
            all_rsids.add(clean_rsid)
    return all_rsids


def get_all_databases():
    """Return all SNP databases as a dict."""
    return {
        "health_risks": HEALTH_RISK_SNPS,
        "nutrigenomics": NUTRIGENOMICS_SNPS,
        "pharmacogenomics": PHARMACOGENOMICS_SNPS,
        "fitness": FITNESS_SNPS,
        "methylation_detox": METHYLATION_DETOX_SNPS,
    }


def lookup_snp(rsid):
    """Look up an rsID across all databases. Returns list of (database_name, entry) tuples."""
    results = []
    databases = get_all_databases()
    for db_name, db in databases.items():
        # Check exact match
        if rsid in db:
            results.append((db_name, db[rsid]))
        # Check suffixed versions (e.g., rs4680_methyl)
        for key, entry in db.items():
            if key.startswith(rsid + "_"):
                results.append((db_name, entry))
    return results
