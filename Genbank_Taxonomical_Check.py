#!/usr/bin/env python3
"""
Author: Arshia Farajollahi
Purpose: Groups GenBank records by taxonomy wordlist match.
Usage: python3 Genbank_Taxonomical_Check.py <input.gb> | -a | --all


The wordlist may be edited for customized matching.
"""

import glob
import os
import sys

from Bio import SeqIO

if len(sys.argv) != 2 or sys.argv[1] in ("-h", "--help"):
    print("Usage: python3 Genbank_Taxonomical_Check.py <input.gb> | -a | --all")
    sys.exit(0 if len(sys.argv) == 2 else 1)

# Check for the -a or --all flag
if sys.argv[1] in ("-a", "--all"):
    target_files = sorted(glob.glob("*.gb"))
    if not target_files:
        print("Error: No .gb files found in the current directory.")
        sys.exit(1)
else:
    if not os.path.isfile(sys.argv[1]):
        print(f"Error: file not found: {sys.argv[1]}")
        sys.exit(1)
    target_files = [sys.argv[1]]

# ─────────────────────────────────────────────────────────────────────────────
# ONLY EDIT taxonomy_wordlist = [] below
# Plain string   → "Bacteria"
#   Searches for "Bacteria" in the lineage.
#   Output file: inputname_Bacteria.gb
#
# Prefixed string → "Infraorder;taxonname"
#   Searches for the part AFTER the semicolon ("taxonname") in the lineage.
#   Output file: inputname_Infraorder_taxonname.gb  (semicolon becomes _ in filename)
#
# Example:
#   taxonomy_wordlist = [
#       "Bacteria",
#       "Infraorder1;taxonname1",
#       "TagID2;taxonname2",
#   ]
# ─────────────────────────────────────────────────────────────────────────────
taxonomy_wordlist = [
    "Holocephali",
    "Selachii",
    "Batoidea",
    "Cladistia",
    "Chondrostei",
    "Holostei",
    "Teleostei",
    "Coelacanthiformes",
    "Dipnomorpha",
    "Anura",
    "Caudata",
    "Gymnophiona",
    "Testudines",
    "Crocodylia",
    "Aves",
    "Sphenodontia",
    "Gekkota",
    "Anguimorpha;Anguidae",  # Anguimorpha
    "Anguimorpha;Anniellidae",  # Anguimorpha
    "Anguimorpha;Diploglossidae",  # Anguimorpha
    "Anguimorpha;Helodermatidae",  # Anguimorpha
    "Anguimorpha;Lanthanotidae",  # Anguimorpha
    "Anguimorpha;Shinisauridae",  # Anguimorpha
    "Anguimorpha;Varanidae",  # Anguimorpha
    "Anguimorpha;Xenosauridae",  # Anguimorpha
    "Iguania;Agamidae",  # Iguania
    "Iguania;Chamaeleonidae",  # Iguania
    "Iguania;Dactyloidae",  # Iguania
    "Iguania;Iguanidae",  # Iguania
    "Iguania;Leiosauridae",  # Iguania
    "Iguania;Phrynosomatidae",  # Iguania
    "Iguania;Polychrotidae",  # Iguania
    "Serpentes;Acrochordidae",  # Serpentes
    "Serpentes;Aniliidae",  # Serpentes
    "Serpentes;Anomalepididae",  # Serpentes
    "Serpentes;Anomochilidae",  # Serpentes
    "Serpentes;Boidae",  # Serpentes
    "Serpentes;Bolyeridae",  # Serpentes
    "Serpentes;Colubridae",  # Serpentes
    "Serpentes;Cyclocoridae",  # Serpentes
    "Serpentes;Cylindrophiidae",  # Serpentes
    "Serpentes;Dipsadidae",  # Serpentes
    "Serpentes;Elapidae",  # Serpentes
    "Serpentes;Gerrhopilidae",  # Serpentes
    "Serpentes;Homalopsidae",  # Serpentes
    "Serpentes;Lamprophiidae",  # Serpentes
    "Serpentes;Leptotyphlopidae",  # Serpentes
    "Serpentes;Loxocemidae",  # Serpentes
    "Serpentes;Pareatidae",  # Serpentes
    "Serpentes;Pythonidae",  # Serpentes
    "Serpentes;Tropidophiidae",  # Serpentes
    "Serpentes;Typhlopidae",  # Serpentes
    "Serpentes;Uropeltidae",  # Serpentes
    "Serpentes;Viperidae",  # Serpentes
    "Serpentes;Xenodermatidae",  # Serpentes
    "Serpentes;Xenopeltidae",  # Serpentes
    "Monotremata",
    "Metatheria",
    "Eutheria",
]

if not taxonomy_wordlist:
    print("Error: taxonomy_wordlist is empty. Please add terms before running.")
    sys.exit(1)


def clarity_Y(entry_output):
    if "Y :" in entry_output:
        return f"\033[1m{entry_output}\033[0m"
    return entry_output


# Loop through all targeted files
for file_path in target_files:
    print(f"\n\033[1m{file_path}\033[0m")

    # Reset groups and no_match for each file
    groups = {word: [] for word in taxonomy_wordlist}
    no_match = []

    base_name = os.path.splitext(os.path.basename(file_path))[0]
    output_folder = f"{base_name}_taxonsort"
    os.makedirs(output_folder, exist_ok=True)

    for record in SeqIO.parse(file_path, "genbank"):
        lineage = record.annotations.get("taxonomy", []) + [
            record.annotations.get("organism", "")
        ]
        lineage_lower = {t.strip().lower() for t in lineage}
        raw = record.format("genbank")

        matched = False
        for word in taxonomy_wordlist:
            search_term = word.split(";", 1)[-1].strip().lower()
            if search_term in lineage_lower:
                groups[word].append(raw)
                matched = True
                break

        if not matched:
            no_match.append(raw)

    report_lines = []

    for word in taxonomy_wordlist:
        count = len(groups[word])
        out_path = os.path.join(
            output_folder, f"{base_name}_{word.replace(';', '_')}.gb"
        )
        with open(out_path, "w") as f:
            if count:
                f.writelines(groups[word])
        per_output = f"{'Y' if count else '-'} : {word:35} {count if count > 0 else ''}"
        print(clarity_Y(per_output))
        report_lines.append(per_output)

    no_match_count = len(no_match)
    if no_match:
        with open(os.path.join(output_folder, f"{base_name}_NO-MATCH.gb"), "w") as f:
            f.writelines(no_match)

    per_output = f"{'Y' if no_match_count else '-'} : {'NO-MATCH':35} {no_match_count if no_match_count > 0 else ''}"
    print(clarity_Y(per_output))
    report_lines.append(per_output)

    with open(os.path.join(output_folder, f"result_report_{base_name}.txt"), "w") as f:
        f.write("\n".join(report_lines) + "\n")
