gold_standard_tests = [
    {
        "index": 0,
        "pmid": "PMC5065698",
        "title": "Evidence that the 5p12 Variant rs10941679 Confers Susceptibility to Estrogen-Receptor-Positive Breast Cancer through FGF10 and MRPS30 Regulation",
        "chrom": "5",
        "chromBand": "5p12",
        "disease": "Estrogen-Receptor-Positive Breast Cancer",
        "efoIds": ["EFO_0000305"],
        "mshIds": ["D001943"],
        "indexSnps": ["rs7703618", "rs981782", "rs4866929"],
        "ldSnps": ["Refer to table 1: many"],
        "causativeSnps": ["rs10941679"],
        "putativeTargetGenes": ["RP11-473L15.2", "RP11-473L15.3", "RN7SL383P", "MRPS30", "BRCAT54", "RP11-357F12.1",
                                "HCN1"],
        "realTargetGenes": ["FGF10", "MRPS30"],
        "targetIsNearest": "Multiple",
        "pmidAssignmentConfidence": "HIGH"
    },
    {
        "index": 1,
        "pmid": "PMC4321900",
        "title": "Evidence that breast cancer risk at the 2q35 locus is mediated through IGFBP5 regulation",
        "chrom": "2",
        "chromBand": "2q35",
        "disease": "Breast Cancer",
        "efoIds": ["EFO_0000305"],
        "mshIds": ["D001943"],
        "indexSnps": ["rs13387042"],
        "ldSnps": ["rs4442975", "rs6721996", "rs13387042", "rs13412666", "rs13426489", "rs12465506", "rs12613955",
                   "rs12620095"],
        "causativeSnps": ["rs4442975"],
        "putativeTargetGenes": ["TNP1 (181Kb)", "DIRC3 (243Kb)", "IGFBP5 (345 kb proximal)",
                                "IGFBP2 (376 kb proximal)"],
        "realTargetGenes": ["IGFBP5"],
        "targetIsNearest": "No",
        "pmidAssignmentConfidence": "HIGH"
    },
    {
        "index": 2,
        "pmid": "PMC3140991",
        "title": "A Functional Variant at a Prostate Cancer Predisposition Locus at 8q24 Is Associated with PVT1 Expression",
        "chrom": "8",
        "chromBand": "8q24",
        "disease": "Prostate Cancer",
        "efoIds": ["EFO_0001663"],
        "mshIds": ["D011471"],
        "indexSnps": ["rs620861"],
        "ldSnps": ["Not presented"],
        "causativeSnps": ["rs378854"],
        "putativeTargetGenes": ["FAM84B", "POU5F1P1", "MYC", "PVT1"],
        "realTargetGenes": ["MYC", "PVT1"],
        "targetIsNearest": "Multiple",
        "pmidAssignmentConfidence": "HIGH"
    },
    {
        "index": 3,
        "pmid": "PMC3617380",
        "title": "Functional Variants at the 11q13 Risk Locus for Breast Cancer Regulate Cyclin D1 Expression through Long-Range Enhancers",
        "chrom": "11",
        "chromBand": "11q13",
        "disease": "Breast Cancer",
        "efoIds": ["EFO_0000305"],
        "mshIds": ["D001943"],
        "indexSnps": ["rs614367"],
        "dSnps": ["Not presented"],
        "causativeSnps": ["rs554219", "rs78540526", "rs75915166"],
        "putativeTargetGenes": ["MYEOV1", "CCND1", "ORAOV1", "FGF3", "FGF4", "FGF19"],
        "realTargetGenes": ["CCND1"],
        "targetIsNearest": "Yes",
        "pmidAssignmentConfidence": "HIGH"
    }
]


def main():
    print("GoldStandardTests.main()")
    print_examples()


def print_examples():
    print("print_examples()")
    for i in range(0, len(gold_standard_tests)):
        gold_standard_test = gold_standard_tests[i]
        print("{}. test: {}".format(i, gold_standard_test))
        chrom = gold_standard_test['chrom']
        print("chrom: '{}'".format(chrom))


if __name__ == '__main__':
    main()
