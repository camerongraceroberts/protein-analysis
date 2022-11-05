from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
import pandas as pd
import requests
import io
import urllib.parse


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
          str(protein_analysis.count_amino_acids()['C']))
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

# this function takes a protein accession number as an input (int) and returns a pandas dataframe with information from Uniprot


def uniprot_search(accession):
    request_params = {
        "query": accession,
        "format": "tsv",
        "fields": "id,protein_name,gene_primary,organism_name,cc_catalytic_activity,protein_families,ft_motif,ft_domain,xref_alphafolddb"
        # "columns": "id,entry name,protein names,gene names,organism"
    }

    response = requests.get(
        f"https://rest.uniprot.org/uniprotkb/search?{urllib.parse.urlencode(request_params)}"
    )

    return pd.read_csv(io.StringIO(response.text), sep="\t")


def download_uniprot_search(accession):
    protein_df = uniprot_search(accession)

    features = []
    for (column_name, column_data) in protein_df.items():
        features.append(column_name + ": " + str(column_data[0]))

    p = "\n"
    output_file = p.join(features)

    # to save the file as a .txt
    with open("protein_features.txt", "w") as output:
        output.write(str(output_file))

    print("To view AlphaFold Model visit: " +
          "https://alphafold.ebi.ac.uk/entry/" + accession)

    return


# this function takes your protein sequence as input and proposes a purifcation method
# designed for soluble proteins
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
    header_1 = "Recommended Buffer: " + buffer + ", " + \
        salt + " NaCl, " + reducing_agent + ", " + glycerol
    header_2 = "Protein tags?: " + solubility_tag

    recommendations = [title, header_1, header_2]

    # if you want to generate a list that is separated by returns
    # p = "\n"
    # file1 = p.join(recommendations)

    # # to save the file as a .txt
    # with open("protein_purification_recommendations.txt", "w") as output:
    #     output.write(str(file1))

    # file2 = p.join(materials)
    # with open("protein_purification_material.txt", "w") as output:
    #     output.write(str(file2))

    # file3 = p.join(general_protocol)
    # with open("protein_purification_protocol.txt", "w") as output:
    #     output.write(str(file3))

    print(title)
    print(header_1)
    print(header_2)

    return


materials = [
    "Acrylamide gel for SDS-PAGE",
    "10X Running buffer",
    "Coomassie Blue stain",
    "LB broth or 2X YT",
    "500 mL plastic centrifuge bottles",
    "1000X antibiotic stock",
    "1M IPTG (-20C)",
    "500 mL of 1M stock of buffer at appropriate pH",
    "500 mL of 5M NaCl",
    "1M DTT stock(-20C)",
    "BME",
    "EDTA-free Protease Inhibitor",
    "DNase",
    "Lysozyme",
    "2L plastic buckets",
    "1 cm long plastic tubing",
    "10 mL syringe",
    "Purification column",
    "Affinity resin",
    "Dialysis tubing",
]

general_protocol = [
    "Transform competent E. coli with plasmid",
    "1. Thaw competent cells over ice ~10 min (50 ul per transformation).",
    "2. Add 30-100 ng of DNA before electroporation.",
    "3. Rescue shocked cells with 300 ul SOB/SOC media.",
    "4. Incubate cells at 37C, 900 rpm for 1 hour.",
    "5. Plate entire transformation on an LB agar plate with selective antibiotic."
    "6. Incubate plate at 37C overnight.",
    "Growth and protein expression",
    "1. Scape plate (add 2 mL LB per plate, scrape with bent sterile p200 tip, add 1 mL bacteria into 1L culture)",
    "2. Add 1 mL of 1000X antibiotic stock to 1L LB.",
    "3. Incubate at 37C shaking at 220 rpm for ~4 hours or until the optical density at 600 nm reaches 0.5-0.8.",
    "4. Decrease temperature of shaker to 18C and cool cells on ice for ~10 min.",
    "5. Add 1 mL 1000X stock (0.2M) IPTG to each L of cells.",
    "6. Induce at 18C O/N at 180 rpm.",
    "Harvest cells and enrich protein",
    "1. Transfer cultures to 500 mL or 1L centrifuge bottles (can be non-sterile).",
    "2. Pellet cells at 3500 rpm for 30 min at 4C.",
    "3. Prepare buffer.",
    "4. Resuspend each liter in 10 mL of buffer.",
    "5. Transfer to 50 mL conical tube (freeze at -80C O/N or longer). If needed day of you can flash freeze in liquid nitrogen.",
    "6. Thaw culture at room temperature.",
    "7. Add DNAse, Lysozyme, protease inhibitor and incubate rocking on ice for 30-60 min.",
    "8. Transfer to glass beaker on ice and lyse by sonication or fresh press.",
    "9. Conditions for sonication: 70% amplitude, pulse on 10 sec, off 10 sec for up to 3 min.",
    "10. Centrifuge lysate at 35,000xg for 30 min at 4C.",
    "11. Prepare 1L of buffer.",
    "12. Add resin to glass column (200 - 1000 ul) and wash with molecular grade water to remove storage buffer and then with 3 passes of buffer to equilibrate.",
    "13. Transfer clarified lysate to column containing pre-equilibrated resin",
    "14. Allow protein to bind to resin for appropriate amount of time (refer to product manual)",
    "15. Wash resin bound with protein with 1L of buffer.",
    "16. Raise 1L bucket containing buffer above the column, clamp thin (1 cm) plastic tubing to bucket and use a syringe to pull buffer through, allow buffer to drip into column at steady state with 2 mL head volume.",
    "17. After washing, elute protein with recommended concentration of a competitive agent (5 passes).",
    "18. Run SDS-PAGE of elution as well as samples of other stages of the prep.",
    "19. Dialyze protein overnight at 4C to remove elution reagent and buffer exchange to desired conditions for further use. (10 mL sample to 2L dialysis buffer).",
    "20. Depending on purity of elution, one can proceed with additional purification methods such as size exclusion chromatography.",
    "21. After purification is complete, concentrate protein and store in -80C."
]
