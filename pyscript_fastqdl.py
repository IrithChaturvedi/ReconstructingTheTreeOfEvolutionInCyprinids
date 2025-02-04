import os

# Define the list of species with their corresponding SRR numbers
species_srr_list = [
    ("Poropuntiinae", "Discherodontus_ashmeadi", "SRR4556475"),
    ("Poropuntiinae", "Mystacoleucus_obtusirostris", "SRR4556519"),
    ("Poropuntiinae", "Poropuntius_normani", "SRR4578777"),
    ("Poropuntiinae", "Sawbwa_resplendens", "SRR4556474"),
    ("Poropuntiinae", "Barbonymus_schwanenfeldii", "SRR4556460"),
    ("Poropuntiinae", "Albulichthys_albuloides", "SRR4556487"),
    ("Poropuntiinae", "Cosmochilus_harmandi", "SRR4556481"),
    ("Barbinae", "Cyprinion_semiplotum", "SRR4556501"),
    ("Barbinae", "Capoeta_aculeata", "SRR4556449"),
    ("Schizothoracinae", "Oreinus_dolungensis", "SRR4578800"),
    ("Schizothoracinae", "Percocypris_tchangi", "SRR4556516"),
    ("Smiliogastrinae", "Chagunius_chagunio", "SRR4578766"),
    ("Spinibarbinae", "Spinibarbus_caldwelli", "SRR4578785"),
    ("Acrossocheilinae", "Acrossocheilus_monticola", "SRR4556507"),
    ("Schizopygopsinae", "Gymnodiptychus_integrigymnatus", "SRR4556505")
]

# Iterate through the list of species and run fastq-dl command for each
for genus, species, srr_number in species_srr_list:
    # Create folder for genus_species if it doesn't exist
    folder_name = f'{genus}_{species}'
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    # Run fastq-dl command
    cmd = f'fastq-dl --accession {srr_number} --provider sra -o {folder_name}'
    os.system(cmd)
    