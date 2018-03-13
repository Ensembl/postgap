__author__ = 'wnewell'

gold_standard_examples = [
    dict(
        index=0,
        pmid="PMC5065698",
        title="Evidence that the 5p12 Variant rs10941679 Confers Susceptibility to Estrogen-Receptor-Positive Breast Cancer through FGF10 and MRPS30 Regulation",
        chrom="5",
        chromBand="5p12",
        disease="Estrogen-Receptor-Positive Breast Cancer",
        efoIds=["EFO_0000305"],
        mshIds=["D001943"],
        indexSnps=["rs7703618", "rs981782", "rs4866929"],
        ldSnps=["Refer to table 1: many"],
        causativeSnps=["rs10941679"],
        putativeTargetGenes=["RP11-473L15.2", "RP11-473L15.3", "RN7SL383P", "MRPS30", "BRCAT54", "RP11-357F12.1",
                             "HCN1"],
        realTargetGenes=["FGF10", "MRPS30"],
        targetIsNearest="Multiple",
        pmidAssignmentConfidence="HIGH"),
    dict(
        index=1,
        pmid="PMC4321900",
        title="Evidence that breast cancer risk at the 2q35 locus is mediated through IGFBP5 regulation",
        chrom="2",
        chromBand="2q35",
        disease="Breast Cancer",
        efoIds=["EFO_0000305"],
        mshIds=["D001943"],
        indexSnps=["rs13387042"],
        ldSnps=["rs4442975", "rs6721996", "rs13387042", "rs13412666", "rs13426489", "rs12465506", "rs12613955",
                "rs12620095"],
        causativeSnps=["rs4442975"],
        putativeTargetGenes=["TNP1", "DIRC3", "IGFBP5", "IGFBP2"],
        putativeTargetGenesLong=["TNP1 (181Kb)", "DIRC3 (243Kb)", "IGFBP5 (345 kb proximal)",
                                 "IGFBP2 (376 kb proximal)"],
        realTargetGenes=["IGFBP5"],
        targetIsNearest="No",
        pmidAssignmentConfidence="HIGH"),
    dict(
        index=2,
        pmid="PMC3140991",
        title="A Functional Variant at a Prostate Cancer Predisposition Locus at 8q24 Is Associated with PVT1 Expression",
        chrom="8",
        chromBand="8q24",
        disease="Prostate Cancer",
        efoIds=["EFO_0001663"],
        mshIds=["D011471"],
        indexSnps=["rs620861"],
        ldSnps=["Not presented"],
        causativeSnps=["rs378854"],
        putativeTargetGenes=["FAM84B", "POU5F1P1", "MYC", "PVT1"],
        realTargetGenes=["MYC", "PVT1"],
        targetIsNearest="Multiple",
        pmidAssignmentConfidence="HIGH"),
    dict(
        index=3,
        pmid="PMC3617380",
        title="Functional Variants at the 11q13 Risk Locus for Breast Cancer Regulate Cyclin D1 Expression through Long-Range Enhancers",
        chrom="11",
        chromBand="11q13",
        disease="Breast Cancer",
        efoIds=["EFO_0000305"],
        mshIds=["D001943"],
        indexSnps=["rs614367"],
        ldSnps=["Not presented"],
        causativeSnps=["rs554219", "rs78540526", "rs75915166"],
        putativeTargetGenes=["MYEOV1", "CCND1", "ORAOV1", "FGF3", "FGF4", "FGF19"],
        realTargetGenes=["CCND1"],
        targetIsNearest="Yes",
        pmidAssignmentConfidence="HIGH"),
    dict(
        index=4,
        pmid="PMC4113484",
        title="Obesity-associated variants within FTO form long-range functional connections with IRX3.",
        chrom="16",
        chromBand="16q12",
        disease="Obesity/T2D",
        efoIds=["EFO_0001073"],
        mshIds=["D009765", "D063766"],
        indexSnps=["rs9939609", "rs9930506", "rs17817449", "rs3751812", "rs1421085"],
        ldSnps=["Not presented"],
        causativeSnps=["rs1421085"],
        putativeTargetGenes=["MM9", "FTO", "RPGRIP11", "IRX3"],
        realTargetGenes=["IRX3"],
        targetIsNearest="Single",
        pmidAssignmentConfidence="HIGH"),
    dict(
        index=5,
        pmid="PMC4959911",
        title="(IRX3/5, ARID5B (TF)). FTO Obesity Variant Circuitry and Adipocyte Browning in Humans",
        chrom="16",
        chromBand="16q12",
        disease="Obesity",
        efoIds=["EFO_0001073"],
        mshIds=["D009765", "D063766"],
        indexSnps=["rs9930506"],
        ldSnps=["Not presented"],
        causativeSnps=["rs1421085"],
        putativeTargetGenesLong=["FTO", "CHD9", "RBL2", "IRX3", "RP6RIP1L", "AKTIP", "LOC643802", "IRX3", "CRNDE",
                                 "IRX5", "IRX6 (2.5 Mb window)"],
        putativeTargetGenes=["FTO", "CHD9", "RBL2", "IRX3", "RP6RIP1L", "AKTIP", "LOC643802", "IRX3", "CRNDE",
                             "IRX5", "IRX6"],
        realTargetGenes=["IRX3", "IRX5", "ARID5B"],
        targetIsNearest="Multiple",
        pmidAssignmentConfidence="HIGH"),
    dict(
        index=6,
        pmid="PMC3062476",
        title="From noncoding variant to phenotype via SORT1 at the 1p13 cholesterol locus", chrom="1",
        chromBand="1p13",
        disease="Serum low-density lipoprotein cholesterol (LDL-C) and myocardial infarction",
        efoIds=["EFO_0000612"],
        mshIds=["D006937", "D009203"],
        indexSnps=["rs646776", "rs599839", "rs12740374", "rs629301"],
        ldSnps=["Not presented"],
        causativeSnps=["rs12740374"],
        putativeTargetGenes=["SARS", "CELSR2", "PSRC1", "MYBPHL", "SORT1", "PSMA5", "SYPL2"],
        realTargetGenes=["SORT1"],
        targetIsNearest="Single",
        pmidAssignmentConfidence="HIGH"),
    dict(
        index=7,
        pmid="PMC5023955",
        title="(ABHD8 (++), ANKLE1 (potential)). Functional mechanisms underlying pleiotropic risk alleles at the 19p13.1 breast–ovarian cancer susceptibility locus",
        chrom="19",
        chromBand="19p13.1",
        disease="Breast Cancer",
        efoIds=["EFO_0000305"],
        mshIds=["D001943"],
        indexSnps=["rs56069439", "rs111961716", "rs4808075", "rs10419397", "rs113299211", "rs4808616", "rs13343778"],
        ldSnps=["Not presented"],
        causativeSnps=["rs56069439", "rs111961716"],
        putativeTargetGenes=["ANKLE1", "BABAM1", "ABHD8", "USHBP1"],
        realTargetGenes=["ABHD8", "ANKLE1"],
        targetIsNearest="Multiple",
        pmidAssignmentConfidence="HIGH"),
    dict(
        index=8,
        pmid="PMC_unknown_8q24_PCAT1",
        title="Not available (Long noncoding RNA PCAT1)",
        chrom="8",
        chromBand="8q24",
        disease="Prostate Cancer",
        efoIds=["EFO_0001663"],
        mshIds=["D011471"],
        indexSnps=["rs13254738", "rs6983561"],
        ldSnps=["Not presented"],
        causativeSnps=["rs7463708"],
        putativeTargetGenes=["Not presented"],
        realTargetGenes=["PCAT1"],
        targetIsNearest="Single",
        pmidAssignmentConfidence="HIGH"),
    dict(
        index=9,
        pmid="PMC_unknown_TBR1",
        title="Not available",
        chrom="?",
        chromBand="?",
        disease="intellectual disability / autism",
        efoIds=["EFO_0003758"],
        mshIds=[],
        indexSnps=[],
        ldSnps=[],
        causativeSnps=["rs4500960"],
        putativeTargetGenes=[],
        realTargetGenes=["TBR1"],
        targetIsNearest="No",
        pmidAssignmentConfidence="?"
    ),
    dict(
        index=10,
        pmid="PMC_unknown_SNCA",
        title="Not available",
        chrom="?",
        chromBand="?",
        disease="Parkinson's disease",
        efoIds=["EFO_0002508", "EFO_0004720"],
        mshIds=["D010300"],
        indexSnps=[],
        ldSnps=[],
        causativeSnps=["rs356168"],
        putativeTargetGenes=[],
        realTargetGenes=["SNCA"],
        targetIsNearest="No",
        pmidAssignmentConfidence="?"
    ),
    dict(
        index=11,
        pmid="PMC_unknown_rs9926296",
        title="rs9926296",
        chrom="?",
        chromBand="?",
        disease="?",
        efoIds=[],
        mshIds=[],
        indexSnps=["rs9926296"],
        ldSnps=[],
        causativeSnps=[],
        putativeTargetGenes=[],
        realTargetGenes=[],
        targetIsNearest="",
        pmidAssignmentConfidence=""
    ),
    dict(
        index=12,
        pmid="PMC5476422",
        title="Large, Diverse Population Cohorts of hiPSCs and Derived Hepatocyte-like Cells Reveal Functional Genetic Variation at Blood Lipid-Associated Loci",
        chrom="20",
        chromBand="20q11",
        disease="Blood lipid (HDL, LDL, TG, PP, …)",
        efoIds=[],
        mshIds=[],
        indexSnps=["rs9926296"],
        ldSnps=[],
        causativeSnps=["rs2277862"],
        putativeTargetGenes=[],
        realTargetGenes=["CPNE1"],
        targetIsNearest="",
        pmidAssignmentConfidence=""
    ),
    dict(
        index=13,
        pmid="PMC5476422",
        title="Large, Diverse Population Cohorts of hiPSCs and Derived Hepatocyte-like Cells Reveal Functional Genetic Variation at Blood Lipid-Associated Loci",
        chrom="1",
        chromBand="1p31",
        disease="Blood lipid (HDL, LDL, TG, PP, …)",
        efoIds=[],
        mshIds=[],
        indexSnps=[],
        ldSnps=[],
        causativeSnps=["rs10889356"],
        putativeTargetGenes=[],
        realTargetGenes=["DOCK7, ANGPTL3"],
        targetIsNearest="",
        pmidAssignmentConfidence=""
    ),
    dict(
        index=14,
        pmid="PMC5476422",
        title="Large, Diverse Population Cohorts of hiPSCs and Derived Hepatocyte-like Cells Reveal Functional Genetic Variation at Blood Lipid-Associated Loci",
        chrom="6",
        chromBand="6q22",
        disease="Blood lipid (HDL, LDL, TG, PP, …)",
        efoIds=[],
        mshIds=[],
        indexSnps=[],
        ldSnps=[],
        causativeSnps=["rs10872142"],
        putativeTargetGenes=[],
        realTargetGenes=["FRK"],
        targetIsNearest="",
        pmidAssignmentConfidence=""
    ),
    dict(
        index=15,
        pmid="PMID28753427",
        title="A Genetic Variant Associated with Five Vascular Diseases Is a Distal Regulator of Endothelin-1 Gene Expression",
        chrom="6",
        chromBand="6p24",
        disease="Coronary artery disease,  coronary artery calcification , migraine headache , cervical artery dissection  fi- bromuscular dysplasia  and hypertension ",
        efoIds=[],
        mshIds=[],
        indexSnps=[],
        ldSnps=[],
        causativeSnps=["rs9926296"],
        putativeTargetGenes=[],
        realTargetGenes=["EDN1"],
        targetIsNearest="",
        pmidAssignmentConfidence=""
    ),
    dict(
        index=16,
        pmid="PMID24115442",
        title="An erythroid enhancer of BCL11A subject to genetic variation determines fetal hemoglobin level",
        chrom="2",
        chromBand="2p16",
        disease="Fetal hemoglobin levels",
        efoIds=[],
        mshIds=[],
        indexSnps=[],
        ldSnps=[],
        causativeSnps=["rs1427407", "rs7606173"],
        putativeTargetGenes=[],
        realTargetGenes=["BCL11A"],
        targetIsNearest="",
        pmidAssignmentConfidence=""
    )

]

import pprint
import os


def main():
    print("GoldStandardTests.main()")
    # print_examples()
    postgap_gold_standard_dataset = PostgapGoldStandardDataset()


def print_examples():
    print("print_examples()")
    for i in range(0, len(gold_standard_examples)):
        gold_standard_test = gold_standard_examples[i]
        print("{}. test: {}".format(i, gold_standard_test))
        chrom = gold_standard_test['chrom']
        print("chrom: '{}'".format(chrom))
    # filter tests with no "realTargetGene"
    non_empty_real_target_examples = []
    for i, example in enumerate(gold_standard_examples):
        print("{}, {}".format(i, example))


class PostgapGoldStandardDataset:
    '''
    Extract data from full POSTGAP data
    containing all SNPs and all genes from the gold standard examples.
    Write to file Postagap.20180209.goldstandard.txt
    for testing assignments made by POSTGAP
    '''

    def __init__(self):
        print("\n\nnew PostgapGoldStandardDataset")
        self.non_empty_real_target_gene_examples = self.fetch_non_empty_real_target_gene_examples()
        self.snp_to_gene_dict = self.extract_gold_standard_associations(self.non_empty_real_target_gene_examples)
        self.print_associations(self.snp_to_gene_dict)
        self.extract_postgap_data(self.non_empty_real_target_gene_examples)

    def fetch_non_empty_real_target_gene_examples(self):
        print("fetch_non_empty_real_target_gene_examples(), orig: {}".format(len(gold_standard_examples)))
        gold_standard_examples_filtered = []
        for i in range(0, len(gold_standard_examples)):
            gold_standard_example = gold_standard_examples[i]
            if (gold_standard_example['realTargetGenes'] == []):
                print("no real target genes")
                continue
            else:
                gold_standard_examples_filtered.append(gold_standard_example)
        print("fetch_non_empty_real_target_gene_examples(), filtered: {}".format(len(gold_standard_examples_filtered)))
        return gold_standard_examples_filtered

    def extract_gold_standard_associations(self, examples):
        ''' make snp([genes]) dict'''
        print("extract_gold_standard_associations")
        snp_to_gene_dict = dict()
        for i in range(0, len(examples)):
            example = examples[i]
            causative_snps = example['causativeSnps']
            real_target_genes = example['realTargetGenes']
            print("{}. causative SNPs: {}. real target genes: {}".format(i, causative_snps, real_target_genes))
            for j in range(0, len(causative_snps)):
                causative_snp = causative_snps[j]
                print("\t add causative SNP: {}".format(causative_snp))
                genes = snp_to_gene_dict.get(causative_snp, set())
                genes.update(real_target_genes)
                snp_to_gene_dict[causative_snp] = genes
            print("snp_to_gene_dict: {}".format(snp_to_gene_dict))
        return snp_to_gene_dict

    def print_associations(self, snp_to_gene_dict):
        print("print_associations()")
        pp = pprint.PrettyPrinter(indent=4)
        pp.pprint(snp_to_gene_dict)

    def extract_postgap_data(self, examples):
        '''
        Extract all SNPs and genes from Gold Standard examples.
        Extract all POSTGAP lines containing GS SNPs and genes
        :return:
        '''
        print("\nextract_postgap_data()")
        # extract all SNPs
        (all_snps, all_genes) = self.extract_gold_standard_snps_and_genes(examples)
        self.extract_postgap_lines(all_snps, all_genes)

    def extract_gold_standard_snps_and_genes(self, examples):
        print("\nextract_gold_standard_snps, examples: {}".format(len(examples)))
        snp_set = set()
        gene_set = set()
        for i, example in enumerate(examples):
            snp_set.update(example['indexSnps'])
            snp_set.update(example['ldSnps'])
            snp_set.update(example['causativeSnps'])
            gene_set.update(example['putativeTargetGenes'])
            gene_set.update(example['realTargetGenes'])
        print("{}. snps: {}, {}".format(i, len(snp_set), snp_set))
        print("{}. genes: {}, {}".format(i, len(gene_set), gene_set))
        snps = sorted(snp_set)
        genes = sorted(gene_set)
        print("{}. snps: {}, {}".format(i, len(snps), snps))
        print("{}. genes: {}, {}".format(i, len(genes), genes))
        return snps, genes

    def extract_postgap_lines(self, snps, genes):
        print("\nextract_postgap_lines, snps: {}, genes: {}".format(len(snps), len(genes)))
        cwd = os.getcwd()
        print("cwd: {}".format(cwd))
        out_dir = cwd + "/../tests/sample_data/qaqc_data"
        print("out_dir: {}".format(out_dir))
        out_filename = "postgap.20180209.gold_standard_snps.txt"
        out_filepath = out_dir + "/" + out_filename
        out_fp = open(out_filepath, "w")
        data_dir = "/Users/otvisitor/Documents/biodata/s2g_comparison_data/postgap_20180209"
        filename = "postgap.20180209.txt"
        filepath = data_dir + "/" + filename
        n_lines_read = 0
        n_snp_lines_found = 0
        n_gene_lines_found = 0
        with open(filepath, "r") as fp:
            for line in fp:
                n_lines_read += 1
                if (n_lines_read == 1):
                    out_fp.write(line)
                if (n_lines_read % 100000 == 0):
                    print("nLines_read: {}, n_snp_lines: {}, n_gene_lines: {},".format(n_lines_read, n_snp_lines_found,
                                                                                       n_gene_lines_found))
                tokens = line.split("\t")
                if (len(tokens) != 44):
                    continue
                rsid = tokens[0]
                if rsid in snps:
                    n_snp_lines_found += 1
                    out_fp.write(line)
                gene_symbol = tokens[10]
                if gene_symbol in genes:
                    n_gene_lines_found += 1
            print("nLines_read: {}, n_snp_lines: {}, n_gene_lines: {},".format(n_lines_read, n_snp_lines_found, n_gene_lines_found))


if __name__ == '__main__':
    main()
