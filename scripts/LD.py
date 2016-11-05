#! /usr/bin/env python

"""

Copyright [1999-2016] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License")
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

		 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

"""

"""

	Please email comments or questions to the public Ensembl
	developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

	Questions may also be sent to the Ensembl help desk at
	<http://www.ensembl.org/Help/Contact>.

"""
import os
import os.path
import sys
import re
import subprocess
import tempfile

def calculate_window(snp, window_len=500000,population='CEPH',cutoff=0.5,db=0):

    """
    Given a SNP id, calculate the pairwise LD between all SNPs within window_size base pairs.
    """

    SNP_id = snp.rsID
    ### Get the SNP location from ENSEMBL
    position = int(snp.pos)

    ### Define the necessary region.
    from_pos = position - (window_len / 2)
    to_pos = position + (window_len / 2)
    chromosome = snp.chrom
    region = '{}:{}-{}'.format(chromosome,from_pos,to_pos)

    chrom_file = "ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf" % (chromosome)
    if not os.path.isfile(chrom_file):
	return dict([(snp, 1)])

    ### Extract this region out from the 1000 genomes BCF

    region_file, region_file_name = tempfile.mkstemp()
    extract_region_comm = "bcftools view -r %s %s -O v -o %s" % (region, os.path.join(DATABASES_DIR, '1000Genomes', population, chrom_file), region_file_name)

    subprocess.call(extract_region_comm.split(" "))
    region_file = open(region_file_name,'r')
    region_vcf = region_file.read()


    ### Find the order of SNPs in the VCF
    SNPs_order = re.findall('rs[0-9]+', region_vcf)


    ### Calculate the pairwise LD using plink2
    ld_file, ld_file_name = tempfile.mkstemp()
    plinkcomm = "plink --vcf %s --r2 --ld-snp %s --inter-chr --out %s" % (region_file_name, SNP_id, ld_file_name)
    print plinkcomm
    plinkcomm_list = plinkcomm.split(" ")

    try:
        subprocess.call(plinkcomm_list)
    except:
        raise Exception("Plink2 needs to be installed")
        sys.exit()

    ### Remove intermediate LD file
    f = open(ld_file_name + '.ld','r')
    g = f.read().splitlines()

    r2_dict = {}
    for s in [l.split() for l in g][1:]:
        ld = s[-1]
        the_snp = SNP(s[-2].split(';')[0], s[-4], int(s[-3]))
        r2_dict[the_snp] = float(ld)


    pruned_ld_snps = dict([x for x in r2_dict.items() if x[1] > cutoff])

    if db != 1:
        subprocess.call(['rm', region_file_name])
        subprocess.call(['rm', ld_file_name + '.ld', ld_file_name + '.log', ld_file_name + '.nosex'])

    return pruned_ld_snps


def get_lds_from_top_gwas(gwas_snp, ld_snps, population='CEPH', region=None,db=0, cutoff=0.5):
    """
    For large numbers of SNPs, best to specify SNP region with chrom:to-from, e.g. 1:7654947-8155562
    For small numbers (<10), regions are extracted from ENSEMBL REST API.
    SNPs can be inputted in a list or from a file with one SNP id per line.
    """


    if gwas_snp.rsID not in [x.rsID for x in ld_snps]:
        ld_snps.append(gwas_snp)


    snp_position_map = {}
    for s in ld_snps:
        snp_position_map[s.rsID] = s.pos



    assert [x.chrom for x in ld_snps] == ([ld_snps[0].chrom] * len(ld_snps))
    positions = [x.pos for x in ld_snps]
    region = '{}:{}-{}'.format(ld_snps[0].chrom, min(positions), max(positions))

    chromosome = ld_snps[0].chrom


    rsID_file, rsID_file_name = tempfile.mkstemp()
    h = open(rsID_file_name, 'w')
    for x in ld_snps:
        h.write('{}\n'.format(str(x.rsID)))
    h.close()


    ### Extract the required region from the VCF
    chrom_file = 'ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf' % (chromosome)
    if not os.path.isfile(chrom_file):
	return dict((snp, 1) for snp in ld_snps)

    region_file, region_file_name = tempfile.mkstemp()
    extract_region_comm = "bcftools view -r %s %s -O z -o %s" % (region, os.path.join(DATABASES_DIR, '1000Genomes', population, chrom_file), region_file_name)
    subprocess.call(extract_region_comm.split(" "))


    ### Extract the list of SNP ids from this region
    vcfcomm = "vcftools --gzvcf %s --snps %s --recode --stdout" % (region_file_name, rsID_file_name)
    print vcfcomm
    vcf = subprocess.check_output(vcfcomm.split(" "))


    ### Remove intermediate region VCF file


    snp_file, snp_file_name = tempfile.mkstemp()
    f = open(snp_file_name, 'w')
    f.write(vcf)
    f.close()

    ### Extract out the order of SNPs
    SNPs_order = re.findall('rs[0-9]+', vcf)


    ### Use plink2 to calculate pairwise LD between these SNPs.
    ld_file, ld_file_name = tempfile.mkstemp()
    plinkcomm = "plink --vcf %s --r2 square --out %s" % (snp_file_name, ld_file_name)
    print plinkcomm
    plinkcomm_list = plinkcomm.split(" ")
    subprocess.call(plinkcomm_list)



    ### Read from the generated results file and output an array.
    LD_file = open(ld_file_name + '.ld','r')
    g = LD_file.read()
    LD_array = [x.split('\t') for x in g.splitlines()]
    LD_file.close


    snp_index = SNPs_order.index(gwas_snp.rsID)
    ld_vector = LD_array[snp_index]


    SNPs = [SNP(snp_id, chromosome, snp_position_map[snp_id]) for snp_id in SNPs_order]
    ld_scores = [(SNPs[i],float(ld_vector[i])) for i in range(len(ld_vector))]
    ld_scores_cutoff = dict([x for x in ld_scores if x[1] > cutoff])

    ### Remove intermediate files
    if db != 1:
        subprocess.call(['rm', ld_file_name + '.ld', ld_file_name + '.log', ld_file_name + '.nosex'])
        subprocess.call(['rm', snp_file_name])
        subprocess.call(['rm', region_file_name])
        subprocess.call(['rm', rsID_file_name])

    return ld_scores_cutoff

