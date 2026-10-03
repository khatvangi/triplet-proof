#!/usr/bin/env python
"""
assemble POSITION_RATES.json from the per-source outputs:
  katoh2023_position_summary.json   (scripts/katoh_position_rates.py)
  landerer2024_position_summary.json (scripts/landerer_position_rates.py)
  landerer2024_bootstrap.json       (scripts/landerer_dataset_bootstrap.py)
plus Kramer & Farabaugh 2007 Table 1, WT column, transcribed verbatim from
raw/kramer2007/87tbl1.jpg (values x 1e-4, mean +/- error as printed).
ratios are normalised to eps2 = 1.
"""
import json

kat = json.load(open("katoh2023_position_summary.json"))
lan = json.load(open("landerer2024_position_summary.json"))
boot = json.load(open("landerer2024_bootstrap.json"))


def katohEntry(key):
    r = kat[key]
    p = r["positions"]
    return dict(
        definition=r["definition"],
        per_position={k: {x: p[k][x] for x in ("n_tested", "n_misdecoded", "freq", "freq_wilson95",
                                              "mean_intensity", "share_of_events")} for k in p},
        eps1_eps2=r["ratio_freq_eps1_eps2"], eps1_eps2_katz95=r["ratio_freq_eps1_eps2_katz95"],
        eps3_eps2=r["ratio_freq_eps3_eps2"], eps3_eps2_katz95=r["ratio_freq_eps3_eps2_katz95"],
        eps1_eps2_by_mean_intensity=r["ratio_meanint_eps1_eps2"],
        eps3_eps2_by_mean_intensity=r["ratio_meanint_eps3_eps2"],
    )


def landererEntry(org):
    r = lan[org]
    q = r["per_opportunity_all_pairs"]
    qo = r["per_opportunity_observed_pairs"]
    b = boot[org]
    return dict(
        n_psm=r["n_psm"], n_unique_sites=r["n_unique_sites"], n_datasets=r["n_datasets"],
        n_psm_dropped_codon_origin_mismatch=r["n_psm_dropped_codon_origin_mismatch"],
        class_counts_psm=r["class_counts_psm"],
        share_single_mismatch_psm=r["share_single_mismatch_psm"],
        share_single_mismatch_sites=r["share_single_mismatch_sites"],
        per_position={k: q[k] for k in "123"},
        eps1_eps2=q["eps1_over_eps2"], eps3_eps2=q["eps3_over_eps2"],
        dataset_bootstrap_2p5_50_97p5=b["boot_2p5_50_97p5"],
        leave_largest_dataset_out=dict(dataset=b["largest_dataset"],
                                       its_share_of_single_mismatch_psm=b["largest_dataset_share_of_single_mismatch_psm"],
                                       **b["leave_largest_out"]),
        sensitivity_observed_pairs=dict(eps1_eps2=qo["eps1_over_eps2"], eps3_eps2=qo["eps3_over_eps2"]),
        sensitivity_no_cys=dict(eps1_eps2=lan[org + "_no_cys"]["per_opportunity_all_pairs"]["eps1_over_eps2"],
                                eps3_eps2=lan[org + "_no_cys"]["per_opportunity_all_pairs"]["eps3_over_eps2"]),
        top_pairs_by_position=r["top_pairs_by_position"],
    )


# kramer & farabaugh 2007 table 1, WT column (x1e-4), transcribed from raw/kramer2007/87tbl1.jpg
kramer_wt = [
    dict(codon="UAA", aa="Ter", position_misread=1, freq_e4=4.1, err_e4=0.51),
    dict(codon="UAG", aa="Ter", position_misread=1, freq_e4=14, err_e4=0.66),
    dict(codon="AUA", aa="Ile", position_misread=2, freq_e4=3.5, err_e4=0.44),
    dict(codon="AGA", aa="Arg", position_misread=2, freq_e4=36, err_e4=2.7),
    dict(codon="AGG", aa="Arg", position_misread=2, freq_e4=31, err_e4=1.5),
    dict(codon="AAU", aa="Asn", position_misread=3, freq_e4=16, err_e4=0.78),
    dict(codon="AAC", aa="Asn", position_misread=3, freq_e4=3.8, err_e4=0.52),
    dict(codon="Others", aa="-", position_misread=None, freq_e4=3.1, err_e4=0.74,
         note="table footnote b: 'Average of all remaining tested codons'"),
]

out = dict(
    note=("eps_k = misreading attributable to a single codon-anticodon mismatch at mRNA codon position k "
          "(k=3 is the wobble position), normalised to eps2 = 1. only amino-acid-changing misreading is "
          "observable in any of these assays; synonymous position-3 misreading is invisible."),
    sources=[
        dict(id="katoh2023", citation="Katoh T, Suga H. 2023. Nucleic Acids Res 51(19):10642-10652",
             doi="10.1093/nar/gkad707", pmcid="PMC10602915", organism="E. coli translation machinery, in vitro",
             method=("reconstituted E. coli cell-free translation (FIT system) of reporter mRNAs with the "
                     "cognate amino acid of the test codon omitted; misincorporation of each of 19 amino "
                     "acids scored by LC/MS as relative peak intensity 0-5 (61 codons x 19 aa = 1159)"),
             used="Supplementary Figure 8 (pdf pp. 12-15; per-position match/mismatch and intensity per "
                  "combination), parsed from vector text; validated against main-text Figure 4",
             status="computed",
             strict=katohEntry("strict"), authors_class=katohEntry("authors"),
             source_inconsistencies=kat["source_inconsistencies"]),
        dict(id="landerer2024_ecoli", citation="Landerer C, Poehls J, Toth-Petroczy A. 2024. Mol Biol Evol 41(3):msae048",
             doi="10.1093/molbev/msae048", pmcid="PMC10939442", organism="E. coli",
             method=("re-analysis of public PRIDE shotgun proteomes by open search (MSFragger); mass shifts "
                     "equal to an amino-acid difference are called substitutions; I/L merged; known-PTM "
                     "matches removed"),
             used="Data S4 (ecoli_filtered/*_substitution_errors.csv, *_codon_counts.csv); position "
                  "assignment is ours (hamming-1 near-cognate rule), not the authors'",
             status="computed", **landererEntry("ecoli")),
        dict(id="landerer2024_yeast", citation="Landerer C, Poehls J, Toth-Petroczy A. 2024. Mol Biol Evol 41(3):msae048",
             doi="10.1093/molbev/msae048", pmcid="PMC10939442", organism="S. cerevisiae",
             method="as landerer2024_ecoli",
             used="Data S4 (yeast_filtered/*); position assignment is ours",
             status="computed", **landererEntry("yeast")),
        dict(id="landerer2024_figS5", citation="Landerer et al. 2024, Supplementary Figure S5c,d",
             doi="10.1093/molbev/msae048", organism="E. coli, S. cerevisiae",
             used="Fig. S5c,d: 'Percentage of observed amino acid misincorporation caused by a given codon position'",
             status="figure only (not digitised); authors' text: 'Misincorporations caused by a single "
                    "mismatch are most frequently observed at the second codon position'"),
        dict(id="kramer2007", citation="Kramer EB, Farabaugh PJ. 2007. RNA 13(1):87-96",
             doi="10.1261/rna.294907", pmcid="PMC1705757", organism="E. coli, in vivo",
             method=("firefly luciferase K529 reporter; each near-cognate codon of AAA/AAG substituted at "
                     "K529; residual activity / wild-type activity = frequency of misreading by tRNA-Lys(UUU)"),
             used="Table 1, WT column (transcribed from the typeset table image, raw/kramer2007/87tbl1.jpg)",
             status=("partial: Table 1 lists only the 7 near-cognates elevated above background plus an "
                     "'Others' average; per-codon values for the remaining near-cognates (incl. all four "
                     "sense position-1 codons CAA, CAG, GAA, GAG) are in Figure 2 (bar chart) only. "
                     "eps ratios not computed."),
             table1_wt=kramer_wt),
        dict(id="mordret2019", citation="Mordret E et al. 2019. Mol Cell 75(3):427-441.e5",
             doi="10.1016/j.molcel.2019.06.041", pmid="31353208", organism="E. coli, S. cerevisiae",
             status=("not accessed: cell.com and sciencedirect return HTTP 403 (Cloudflare) to curl and to a "
                     "headless browser; bioRxiv 10.1101/255943 returns 429; no PMC copy; github "
                     "ernestmordret/substitutions has pipeline code only. Per-position mismatch data are in "
                     "Figure 3A (4x4 confusion matrices shaded by log count) = figure only."),
             ),
        dict(id="manickam2016", citation="Manickam N, Joshi K, Bhatt MJ, Farabaugh PJ. 2016. Nucleic Acids Res 44(4):1871 (end page not verified)",
             doi="10.1093/nar/gkv1506", pmcid="PMC4770228", organism="E. coli, in vivo",
             status="figure only: misreading frequencies are shown in bar-chart figures; supplement holds tRNA modification/abundance methods only"),
    ],
)
json.dump(out, open("POSITION_RATES.json", "w"), indent=2)
print("wrote POSITION_RATES.json")
