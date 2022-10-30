from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv


# takes a str input of your protein sequence and returns the protein properties including AA #, MW, e, and # of Cys residues
def protein_features(sequence):
    protein_analysis = ProteinAnalysis(sequence)
    amino_acid_number = len(sequence)
    molecular_weight = protein_analysis.molecular_weight()/1000

    # [reduced, oxidized]
    extinction_coef = protein_analysis.molar_extinction_coefficient()

    print("Number of amino acids: " + str(amino_acid_number))
    print("Molecular weight: " + str(molecular_weight) + " kDa")
    print("Number of Cysteines: " +
          str(protein_analysis.amino_acid_number['C']))
    print("Extinction coefficient with disulfid bridges: " +
          str(extinction_coef[1]))
    print("Extinction coefficient with reduced cysteines: " +
          str(extinction_coef[0]))
    print("To calculate molar concentration of protein stock --> Optical density at 280/extinction coefficient = M")

    return [amino_acid_number, molecular_weight, extinction_coef[0], extinction_coef[1]]


# input value of optical density at 280 nm, and "y" or "n" for reducing agent, and str of protein sequence to get concentration of sample
def calculate_protein_conc(OD280, reducing_agent, sequence):
    if reducing_agent == "y":
        extinction_coef = protein_features(sequence)[2]
    else:
        extinction_coef = protein_features(sequence)[3]

    molecular_weight = protein_features(sequence)[1]
    molar_concentration = OD280/extinction_coef
    total_mg = molar_concentration * (molecular_weight*1000)
    print("[mol/L, mg/ml]")
    return [molar_concentration, total_mg]  # [mol/L, mg/mL]


# this function takes your protein sequence as input and proposes a purifcation method
def protein_prep_protocol(sequence, protein_name):
    protein = ProteinAnalysis(sequence)
    instability_index = protein.instability_index()
    isoelectric_point = protein.isoelectric_point()
    cys_count = protein.count_amino_acids()['C']

    if instability_index < 40:
        glycerol = ""
        salt = "200 mM"
        solubility_tag = "solubility tag is likely not needed"
    else:
        glycerol = "5-10% glycerol"
        salt = "300 - 500 mM"
        solubility_tag = "solubility tag recommended"

    buffer_pH = isoelectric_point + 1

    if buffer_pH <= 6.5:
        buffer = "25 mM MES"
    elif buffer_pH <= 8:
        buffer = "25 mM HEPES"
    elif buffer_pH <= 9:
        buffer = "25 mM Tris"
    else:
        buffer = "25 mM CAPSO"

    if cys_count > 1:
        reducing_agent = "2 mM BME"
    else:
        reducing_agent = ""

    # Write a file
    title = "Recommended Protein Purification for " + protein_name
    header_1 = "Buffer: " + buffer + " " + salt + \
        " " + reducing_agent + " " + glycerol
    header_2 = "Protein tags: " + solubility_tag

    return [title, header_1, header_2]


materials = [
    "acrylamide gel for SDS-PAGE",
    "10X Running buffer",
    "Coomassie Blue stain",
    "LB broth",
    "500 mL plastic centrifuge bottles",
    "1000X antibiotic stock",
    "1M IPTG (-20C)",
    "1M stock of buffer at appropriate pH",
    "5M NaCl stock",
    "1M DTT stock(-20C)",
    "BME",
    "EDTA-free Protease Inhibitor",
    "DNase",
    "Lysozyme",
    "2L plastic buckets",
    "1 cm plastic tubing",
    "10 mL syringe",
    "Purification column",
    "Nickel or cobalt resin",
    "FPLC size exclusion column",
    "Dialysis tubing"
]

protocol_1 = [
    "Transform competent E coli with plasmid",
    "1. Thaw competent cells over ice ~10 min (50 ul per transformation).",
    "2. Add 30-100 ng of DNA before electroporation.",
    "3. Rescue shocked cells with 300 ul SOB/SOC media.",
    "4. Incubate cells at 37C, 900 rpm for 1 hour.",
    "5. Plate entire transformation on an LB agar plate with selective antibiotic."
    "6. Incubate plate at 37C overnight."
]
