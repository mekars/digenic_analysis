import os
import subprocess
from itertools import combinations
import pandas as pd
from rpy2.robjects import r, pandas2ri

pandas2ri.activate()

# Define input files
variantIDs_genes_file = "variantIDs_genes.txt"
plink_input_file = "input"
fam_file = "input.fam"
phenotypes_with_pcs_file = "phenotypes_with_pcs.txt"

# Step 2: Convert to raw file format
subprocess.run(["plink", "--bfile", "input", "--recode", "A", "--out", "input_raw"])

# Step 3: Read the raw data and create a new file with the required format
input_raw = "input_raw.raw"
output_file = "input.txt"

with open(input_raw, "r") as raw_file:
    lines = raw_file.readlines()

# Extract sample ID and genotype data, remove the header line, and add a new header line containing sample ID, variant IDs, and gene names
header = lines[0].strip().split()
sample_data = [line.strip().split() for line in lines[1:]]

variant_names = [line.split()[0] + "_" + line.split()[1] for line in open(variantIDs_genes_file)]
header_line = ["ID"] + variant_names

# Write the new header and data to the output file
with open(output_file, "w") as out_file:
    out_file.write(" ".join(header_line) + "\n")
    for data in sample_data:
        out_file.write(" ".join([data[1]] + data[6:]) + "\n")

# Step 4: Extract variant pairs
def binary_combinations(input_array):
    binary_combos = []
    for combo in combinations(input_array, 2):
        binary_combos.append(list(combo))
    return binary_combos

input_file_path = "input.txt"
output_file_path = "variant_pairs.txt"
with open(input_file_path, "r") as file:
    first_line = file.readline().strip()
    columns = first_line.split()
    variant_names = columns[1:]

    with open(output_file_path, "w") as ofile:
        for line in file:
            variant_indexes_and_values = []
            values = line.strip().split()
            sample_id = values[0]
            for idx, value in enumerate(values[1:]):
                if value in ["1", "2"]:
                    variant_indexes_and_values.append((idx, value))
            combs = binary_combinations(variant_indexes_and_values)

            for binary_combination in combs:
                ofile.write(
                    f"{sample_id}\t{'HET' if binary_combination[0][1]=='1' else 'HOM'}_{variant_names[binary_combination[0][0]]}@"
                    f"{'HET' if binary_combination[1][1]=='1' else 'HOM'}_{variant_names[binary_combination[1][0]]}\n"
                )

# Step 5: Prepare the fam file and sort variant_pairs.txt in the order of "proband", "father" and "mother"
fam_df = pd.read_csv(fam_file, delim_whitespace=True, header=None)
fam_df = fam_df.iloc[:, [1, 2, 3]]
fam_df = fam_df[(fam_df[1] != "0") & (fam_df[2] != "0")]
fam_df.to_csv("fam.txt", sep="\t", index=False, header=False)

command1 = "awk '{for(i=1;i<=NF;i++){print $i \"\\t\" ++count}}' fam.txt > tmp1.txt"
command2 = "awk 'FNR==NR{a[$1]=$0;next}{if($1 in a){print $0 \"\\t\" a[$1]}}' tmp1.txt variant_pairs.txt > tmp2.txt"
command3 = "sort -g -k 4 tmp2.txt | cut -f -2 > variant_pairs.txt"
subprocess.run(command1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
subprocess.run(command2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
subprocess.run(command3, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
os.remove("tmp1.txt")
os.remove("tmp2.txt")

# Step 6: Filter variant pairs based on family relationships
does_print = False
family_input_file_name = "fam.txt"
variant_pairs_file_name = "variant_pairs.txt"
output_file_name = "filtered_variant_pairs.txt"

try:
    os.remove(output_file_name)
except Exception as e:
    pass

variants_of_proband = set()
previous_proband = ""
variants_of_previous_proband = set()

with open(variant_pairs_file_name, "r") as variant_pairs_file:
    with open(family_input_file_name, "r") as family_file:
        for line in family_file:
            if line.isspace():
                continue
            family_values = line.strip().split("\t")
            proband_id = family_values[0]
            parent1_id = family_values[1]
            parent2_id = family_values[2]
            if does_print:
                print(f"{proband_id} {parent1_id} {parent2_id}")
                print("----------------------------------")
            if len(variants_of_previous_proband) != 0:
                if previous_proband != proband_id:
                    variants_of_previous_proband.clear()
                else:
                    variants_of_proband.update(variants_of_previous_proband)
                    variants_of_previous_proband.clear()
            else:
                variant_pairs_line = variant_pairs_file.readline()
            while variant_pairs_line:
                variant_pairs_line = variant_pairs_line.strip()
                variant_pairs_line_values = variant_pairs_line.split('\t')
                sample_id = variant_pairs_line_values[0]
                variants = variant_pairs_line_values[1]
                variant_pairs_line = variant_pairs_file.readline()
                if does_print:
                    print(f"{sample_id} {variants}")
                if sample_id not in [proband_id, parent1_id, parent2_id]:
                    previous_proband = sample_id
                    variants_of_previous_proband.clear()
                    variants_of_previous_proband.add(variants)
                    break
                if sample_id == proband_id:
                    variants_of_proband.add(variants)
                else:
                    if variants in variants_of_proband:
                        variants_of_proband.remove(variants)
            with open(output_file_name, "a") as ofile:
                for variants in variants_of_proband:
                    ofile.write(f"{proband_id} {variants}\n")
                variants_of_proband.clear()


# Step 7: Select gene pairs to test and perform burden testing
header = ["ID", "PC1", "PC2", "PC3", "phenotype", "carrier"]
phenotypes_df = pd.read_csv(phenotypes_with_pcs_file, sep="\t", header=None)
phenotypes_df.columns = ["ID", "PC1", "PC2", "PC3", "phenotype"]

filtered_variants_df = pd.read_csv("filtered_variant_pairs.txt", delim_whitespace=True, header=None)
filtered_variants_df.columns = ["SampleID", "VariantPair"]

# Extract gene pairs
filtered_variants_df['GenePair'] = filtered_variants_df['VariantPair'].apply(
    lambda vp: '@'.join([v.split('_')[-1] for v in vp.split('@')])
)

# Drop duplicate gene pairs within the same sample
unique_gene_pairs_df = filtered_variants_df.drop_duplicates(subset=["SampleID", "GenePair"])

# Filter gene pairs with more than 5 observations
gene_pair_counts = unique_gene_pairs_df['GenePair'].value_counts()
gene_pairs_to_test = gene_pair_counts[gene_pair_counts > 5].index

# Iterate over filtered gene pairs and perform Firth's logistic regression
burden_output_file = "final_burden_testing_result.txt"
with open(burden_output_file, "w") as burden_file:
    burden_file.write("GenePair\tP\tOR\tL95\tU95\tCaseCount\tControlCount\n")

for gene_pair in gene_pairs_to_test:
    current_pair = unique_gene_pairs_df[unique_gene_pairs_df['GenePair'] == gene_pair]
    sample_ids = current_pair['SampleID'].tolist()
    case_count = phenotypes_df[phenotypes_df['ID'].isin(sample_ids) & (phenotypes_df['phenotype'] == 1)].shape[0]
    control_count = phenotypes_df[phenotypes_df['ID'].isin(sample_ids) & (phenotypes_df['phenotype'] == 0)].shape[0]
    current_phenotypes_df = phenotypes_df.copy()
    current_phenotypes_df['carrier'] = current_phenotypes_df['ID'].apply(lambda x: 1 if x in sample_ids else 0)
    current_phenotypes_df.to_csv("currentinput.txt", sep="\t", index=False)
    r_script = """
    library(logistf)
    data <- read.csv("currentinput.txt", header = TRUE, sep = "\t")
    model <- logistf(phenotype ~ carrier + PC1 + PC2 + PC3, data = data)
    p_value <- model$prob[2]
    odds_ratio <- exp(coef(model)["carrier"])
    ci <- exp(confint(model)["carrier", ])
    result <- data.frame(GenePair="%s", P=p_value, OR=odds_ratio, L95=ci[1], U95=ci[2], CaseCount=%d, ControlCount=%d)
    write.table(result, "final_burden_testing_result.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append=TRUE)
    """ % (gene_pair, case_count, control_count)
    r(r_script)

    os.remove("currentinput.txt")
