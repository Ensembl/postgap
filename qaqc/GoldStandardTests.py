array_of_dict = [
    {
        "a": [1, 2, 3],
        "b": [2, 3, 4],
        "c": [3, 4, 5]
    },
    {
        "x": 10,
        "y": 11,
        "z": 12
    }
]

gold_standard_tests = [
    {
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
    }
]

def main():
    print("GoldStandardTests.main()")

if __name__ == '__main__':
    main()
