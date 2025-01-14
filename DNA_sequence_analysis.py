# Required packages
from Bio import SeqIO
from itertools import product

# Dictionary mapping ambiguous codes to possible bases
ambiguous_codes = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "K": "GT",
    "M": "AC",
    "N": "ACGT",
    "R": "AG",
    "S": "CG",
    "V": "ACG",
    "W": "AT",
    "Y": "CT"
}

# Function to expand ambiguous DNA sequences into all possible DNA sequences they may represent. This is necessary to
# handle restriction enzyme recognition sites with ambiguous codes.
def expand_recognition_sequence(sequence, code_map):
    bases = [code_map[char] for char in sequence]
    possible_sequences = ["".join(comb) for comb in product(*bases)]
    return possible_sequences


# Function to find ORFs
def find_orfs(sequence, start_codon="ATG", stop_codons=("TAA", "TAG", "TGA")):
    orfs = []
    for frame in range(3):
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i + 3]
            if codon == start_codon:
                for j in range(i, len(sequence) - 2, 3):
                    stop_codon = sequence[j:j + 3]
                    if stop_codon in stop_codons:
                        orfs.append(sequence[i:j + 3])
                        break
    return orfs


# Function to find restriction sites without overlapping matches
def find_restriction_sites(sequence, enzyme_dict):
    restriction_sites = {}
    for enzyme, sites in enzyme_dict.items():
        positions = []
        for site in sites:
            start = 0
            while True:
                start = sequence.find(site, start)
                if start == -1:
                    break
                positions.append(start + 1)
                start += 1

        if positions:
            restriction_sites[enzyme] = (len(positions), positions)
        else:
            restriction_sites[enzyme] = (0, [])
    return restriction_sites


# Function to calculate GC content
def calculate_gc_content(sequence):
    g = sequence.count("G")
    c = sequence.count("C")
    return round((g + c) / len(sequence) * 100, 2)


# Function to parse FASTA files
def parse_fasta(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences


# Restriction enzyme dictionary with recognition sequences including ambiguous codes, translated to exact codes using
# expand_recognition_sequence function
restriction_enzymes = {
    "AatII": expand_recognition_sequence("GACGTC", ambiguous_codes),
    "Acc65I": expand_recognition_sequence("GGTACC", ambiguous_codes),
    "AccI": expand_recognition_sequence("GTMKAC", ambiguous_codes),
    "AclI": expand_recognition_sequence("AACGTT", ambiguous_codes),
    "AfeI": expand_recognition_sequence("AGCGCT", ambiguous_codes),
    "AflII": expand_recognition_sequence("CTTAAG", ambiguous_codes),
    "AflIII": expand_recognition_sequence("ACRYGT", ambiguous_codes),
    "AgeI-HF®": expand_recognition_sequence("ACCGGT", ambiguous_codes),
    "AhdI": expand_recognition_sequence("GACNNNNNGTC", ambiguous_codes),
    "AleI-v2": expand_recognition_sequence("CACNNNNGTG", ambiguous_codes),
    "AluI": expand_recognition_sequence("AGCT", ambiguous_codes),
    "AlwNI": expand_recognition_sequence("CAGNNNCTG", ambiguous_codes),
    "ApaI": expand_recognition_sequence("GGGCCC", ambiguous_codes),
    "ApaLI": expand_recognition_sequence("GTGCAC", ambiguous_codes),
    "ApoI-HF": expand_recognition_sequence("RAATTY", ambiguous_codes),
    "AscI": expand_recognition_sequence("GGCGCGCC", ambiguous_codes),
    "AseI": expand_recognition_sequence("ATTAAT", ambiguous_codes),
    "AsiSI": expand_recognition_sequence("GCGATCGC", ambiguous_codes),
    "AvaI BsoBI": expand_recognition_sequence("CYCGRG", ambiguous_codes),
    "AvaII": expand_recognition_sequence("GGWCC", ambiguous_codes),
    "AvrII": expand_recognition_sequence("CCTAGG", ambiguous_codes),
    "BaeGI": expand_recognition_sequence("GKGCMC", ambiguous_codes),
    "BamHI, BamHI-HF®": expand_recognition_sequence("GGATCC", ambiguous_codes),
    "BanI": expand_recognition_sequence("GGYRCC", ambiguous_codes),
    "BanII": expand_recognition_sequence("GRGCYC", ambiguous_codes),
    "BclI, BclI-HF": expand_recognition_sequence("TGATCA", ambiguous_codes),
    "BfaI": expand_recognition_sequence("CTAG", ambiguous_codes),
    "BglI": expand_recognition_sequence("GCCNNNNNGGC", ambiguous_codes),
    "BglII": expand_recognition_sequence("AGATCT", ambiguous_codes),
    "BlpI": expand_recognition_sequence("GCTNAGC", ambiguous_codes),
    "BmtI-HF®": expand_recognition_sequence("GCTAGC", ambiguous_codes),
    "BsaAI": expand_recognition_sequence("YACGTR", ambiguous_codes),
    "BsaBI": expand_recognition_sequence("GATNNNNATC", ambiguous_codes),
    "BsaHI": expand_recognition_sequence("GRCGYC", ambiguous_codes),
    "BsaJI": expand_recognition_sequence("CCNNGG", ambiguous_codes),
    "BsaWI": expand_recognition_sequence("WCCGGW", ambiguous_codes),
    "BsiEI": expand_recognition_sequence("CGRYCG", ambiguous_codes),
    "BsiHKAI": expand_recognition_sequence("GWGCWC", ambiguous_codes),
    "BsiWI, BsiWI-HF®": expand_recognition_sequence("CGTACG", ambiguous_codes),
    "BslI": expand_recognition_sequence("CCNNNNNNNGG", ambiguous_codes),
    "BsmBI-v2": expand_recognition_sequence("CGTCTC", ambiguous_codes),
    "Bsp1286I": expand_recognition_sequence("GDGCHC", ambiguous_codes),
    "BspDI, ClaI": expand_recognition_sequence("ATCGAT", ambiguous_codes),
    "BspEI": expand_recognition_sequence("TCCGGA", ambiguous_codes),
    "BspHI": expand_recognition_sequence("TCATGA", ambiguous_codes),
    "BsrFI-v2": expand_recognition_sequence("RCCGGY", ambiguous_codes),
    "BsrGI-HF®": expand_recognition_sequence("TGTACA", ambiguous_codes),
    "BssHII": expand_recognition_sequence("GCGCGC", ambiguous_codes),
    "BstAPI": expand_recognition_sequence("GCANNNNNTGC", ambiguous_codes),
    "BstBI": expand_recognition_sequence("TTCGAA", ambiguous_codes),
    "BstEII-HF®": expand_recognition_sequence("GGTNACC", ambiguous_codes),
    "BstNI": expand_recognition_sequence("CCWGG", ambiguous_codes),
    "BstUI": expand_recognition_sequence("CGCG", ambiguous_codes),
    "BstXI": expand_recognition_sequence("CCANNNNNNTGG", ambiguous_codes),
    "BstYI": expand_recognition_sequence("RGATCY", ambiguous_codes),
    "Bsu36I": expand_recognition_sequence("CCTNAGG", ambiguous_codes),
    "BtgI": expand_recognition_sequence("CCRYGG", ambiguous_codes),
    "Cac8I": expand_recognition_sequence("GCNNGC", ambiguous_codes),
    "CviKI-1": expand_recognition_sequence("RGCY", ambiguous_codes),
    "CviQI": expand_recognition_sequence("GTAC", ambiguous_codes),
    "DdeI": expand_recognition_sequence("CTNAG", ambiguous_codes),
    "DpnI": expand_recognition_sequence("GATC", ambiguous_codes),
    "DraI": expand_recognition_sequence("TTTAAA", ambiguous_codes),
    "DraIII-HF®": expand_recognition_sequence("CACNNNGTG", ambiguous_codes),
    "DrdI": expand_recognition_sequence("GACNNNNNNGTC", ambiguous_codes),
    "EaeI": expand_recognition_sequence("YGGCCR", ambiguous_codes),
    "EagI-HF®": expand_recognition_sequence("CGGCCG", ambiguous_codes),
    "Eco53kI": expand_recognition_sequence("GAGCTC", ambiguous_codes),
    "EcoNI": expand_recognition_sequence("CCTNNNNNAGG", ambiguous_codes),
    "EcoO109I": expand_recognition_sequence("RGGNCCY", ambiguous_codes),
    "EcoRI, EcoRI-HF®": expand_recognition_sequence("GAATTC", ambiguous_codes),
    "EcoRV-HF®, EcoRV": expand_recognition_sequence("ACRYGT", ambiguous_codes),
    "FatI": expand_recognition_sequence("CATG", ambiguous_codes),
    "Fnu4HI": expand_recognition_sequence("GCNGC", ambiguous_codes),
    "FseI": expand_recognition_sequence("GGCCGGCC", ambiguous_codes),
    "FspI": expand_recognition_sequence("TGCGCA", ambiguous_codes),
    "HaeII": expand_recognition_sequence("RGCGCY", ambiguous_codes),
    "HaeIII": expand_recognition_sequence("GGCC", ambiguous_codes),
    "HhaI": expand_recognition_sequence("GCGC", ambiguous_codes),
    "HincII": expand_recognition_sequence("GTYRAC", ambiguous_codes),
    "HindIII, HindIII-HF®": expand_recognition_sequence("AAGCTT", ambiguous_codes),
    "HinfI": expand_recognition_sequence("GANTC", ambiguous_codes),
    "HinP1I": expand_recognition_sequence("GCGC", ambiguous_codes),
    "HpaI": expand_recognition_sequence("GTTAAC", ambiguous_codes),
    "Hpy166II": expand_recognition_sequence("GTNNAC", ambiguous_codes),
    "Hpy188I": expand_recognition_sequence("TCNGA", ambiguous_codes),
    "Hpy188III": expand_recognition_sequence("TCNNGA", ambiguous_codes),
    "Hpy99I": expand_recognition_sequence("CGWCG", ambiguous_codes),
    "HpyCH4III": expand_recognition_sequence("ACNGT", ambiguous_codes),
    "HpyCH4IV": expand_recognition_sequence("ACGT", ambiguous_codes),
    "HpyCH4V": expand_recognition_sequence("TGCA", ambiguous_codes),
    "KasI": expand_recognition_sequence("GGCGCC", ambiguous_codes),
    "KpnI-HF®": expand_recognition_sequence("GGTACC", ambiguous_codes),
    "MfeI-HF®": expand_recognition_sequence("CAATTG", ambiguous_codes),
    "MluCI": expand_recognition_sequence("AATT", ambiguous_codes),
    "MluI-HF®": expand_recognition_sequence("ACGCGT", ambiguous_codes),
    "MscI": expand_recognition_sequence("TGGCCA", ambiguous_codes),
    "MseI": expand_recognition_sequence("TTAA", ambiguous_codes),
    "MslI": expand_recognition_sequence("CAYNNNNRTG", ambiguous_codes),
    "MspA1I": expand_recognition_sequence("CMGCKG", ambiguous_codes),
    "MspI, HpaII": expand_recognition_sequence("CCGG", ambiguous_codes),
    "MwoI": expand_recognition_sequence("GCNNNNNNNGC", ambiguous_codes),
    "NaeI": expand_recognition_sequence("GCCGGC", ambiguous_codes),
    "NarI": expand_recognition_sequence("GGCGCC", ambiguous_codes),
    "NciI": expand_recognition_sequence("CCSGG", ambiguous_codes),
    "NcoI, NcoI-HF®": expand_recognition_sequence("CCATGG", ambiguous_codes),
    "NdeI": expand_recognition_sequence("CATATG", ambiguous_codes),
    "NheI-HF®": expand_recognition_sequence("GCTAGC", ambiguous_codes),
    "NlaIII": expand_recognition_sequence("CATG", ambiguous_codes),
    "NlaIV": expand_recognition_sequence("GGNNCC", ambiguous_codes),
    "NotI-HF®, NotI": expand_recognition_sequence("GCGGCCGC", ambiguous_codes),
    "NruI-HF®": expand_recognition_sequence("TCGCGA", ambiguous_codes),
    "NsiI, NsiI-HF®": expand_recognition_sequence("ATGCAT", ambiguous_codes),
    "NspI": expand_recognition_sequence("RCATGY", ambiguous_codes),
    "PacI": expand_recognition_sequence("TTAATTAA", ambiguous_codes),
    "PciI": expand_recognition_sequence("ACATGT", ambiguous_codes),
    "PflMI": expand_recognition_sequence("CCANNNNNTGG", ambiguous_codes),
    "PluTI": expand_recognition_sequence("GGCGCC", ambiguous_codes),
    "PmeI": expand_recognition_sequence("GTTTAAAC", ambiguous_codes),
    "PmlI": expand_recognition_sequence("CACGTG", ambiguous_codes),
    "PpuMI": expand_recognition_sequence("RGGWCCY", ambiguous_codes),
    "PshAI": expand_recognition_sequence("GACNNNNGTC", ambiguous_codes),
    "PsiI-v2": expand_recognition_sequence("TTATAA", ambiguous_codes),
    "PspGI": expand_recognition_sequence("CCWGG", ambiguous_codes),
    "PspOMI": expand_recognition_sequence("GGGCCC", ambiguous_codes),
    "PspXI": expand_recognition_sequence("VCTCGAGB", ambiguous_codes),
    "PstI, PstI-HF®": expand_recognition_sequence("CTGCAG", ambiguous_codes),
    "PvuI-HF®": expand_recognition_sequence("CGATCG", ambiguous_codes),
    "PvuII-HF®, PvuII": expand_recognition_sequence("CAGCTG", ambiguous_codes),
    "RsaI": expand_recognition_sequence("GTAC", ambiguous_codes),
    "RsrII": expand_recognition_sequence("CGGWCCG", ambiguous_codes),
    "SacI-HF®": expand_recognition_sequence("GAGCTC", ambiguous_codes),
    "SacII": expand_recognition_sequence("CCGCGG", ambiguous_codes),
    "SalI-HF®, SalI": expand_recognition_sequence("GTCGAC", ambiguous_codes),
    "Sau3AI, DpnII, MboI": expand_recognition_sequence("GATC", ambiguous_codes),
    "Sau96I": expand_recognition_sequence("GGNCC", ambiguous_codes),
    "SbfI-HF®": expand_recognition_sequence("CCTGCAGG", ambiguous_codes),
    "ScaI-HF®": expand_recognition_sequence("AGTACT", ambiguous_codes),
    "ScrFI": expand_recognition_sequence("CCNGG", ambiguous_codes),
    "SexAI": expand_recognition_sequence("ACCWGGT", ambiguous_codes),
    "SfcI": expand_recognition_sequence("CTRYAG", ambiguous_codes),
    "SfiI": expand_recognition_sequence("GGCCNNNNNGGCC", ambiguous_codes),
    "SfoI": expand_recognition_sequence("GGCGCC", ambiguous_codes),
    "SgrAI": expand_recognition_sequence("CRCCGGYG", ambiguous_codes),
    "SmaI": expand_recognition_sequence("CCCGGG", ambiguous_codes),
    "SmlI": expand_recognition_sequence("CTYRAG", ambiguous_codes),
    "SnaBI": expand_recognition_sequence("TACGTA", ambiguous_codes),
    "SpeI-HF®": expand_recognition_sequence("ACTAGT", ambiguous_codes),
    "SphI-HF®, SphI": expand_recognition_sequence("GCATGC", ambiguous_codes),
    "SrfI": expand_recognition_sequence("GCCCGGGC", ambiguous_codes),
    "SspI-HF®": expand_recognition_sequence("AATATT", ambiguous_codes),
    "StuI": expand_recognition_sequence("AGGCCT", ambiguous_codes),
    "StyD4I": expand_recognition_sequence("CCNGG", ambiguous_codes),
    "StyI-HF®": expand_recognition_sequence("CCWWGG", ambiguous_codes),
    "SwaI": expand_recognition_sequence("ATTTAAAT", ambiguous_codes),
    "TaqI-v2": expand_recognition_sequence("TCGA", ambiguous_codes),
    "TfiI": expand_recognition_sequence("GAWTC", ambiguous_codes),
    "TseI, ApeKI": expand_recognition_sequence("GCWGC", ambiguous_codes),
    "Tsp45I": expand_recognition_sequence("GTSAC", ambiguous_codes),
    "TspMI XmaI": expand_recognition_sequence("CCCGGG", ambiguous_codes),
    "TspRI": expand_recognition_sequence("NNCASTGNN", ambiguous_codes),
    "Tth111I, PflFI": expand_recognition_sequence("GACNNNGTC", ambiguous_codes),
    "XbaI": expand_recognition_sequence("TCTAGA", ambiguous_codes),
    "XcmI": expand_recognition_sequence("CCANNNNNNNNNTGG", ambiguous_codes),
    "XhoI, PaeR7I": expand_recognition_sequence("CTCGAG", ambiguous_codes),
    "XmnI": expand_recognition_sequence("GAANNNNTTC", ambiguous_codes),
    "ZraI": expand_recognition_sequence("GACGTC", ambiguous_codes),

}


# Function to write results to a .txt file
def write_results_to_file(output_file, seq_id, sequence, orfs, gc_content, sorted_sites):
    with open(output_file, "a") as f:
        f.write(f"> {seq_id}\n")
        f.write(f"Sequence: {sequence[:100]}...\n")  # Print the first 100 bases for inspection
        f.write(f"  GC Content: {gc_content}%\n")

        f.write(f"  ORFs (total: {len(orfs)}):\n")
        for orf in orfs:
            f.write(f"    - {orf}\n")

        f.write("  Restriction Sites:\n")
        for enzyme, count, positions in sorted_sites:
            f.write(f"    - {enzyme}: {count} site(s) at positions {positions}\n")

        f.write("\n")  # Add a newline between results


# Function to save the sequence as a .fasta file
def save_sequence_to_fasta(sequence, sequence_name):
    with open(f"{sequence_name}.fasta", "w") as f:
        f.write(f">{sequence_name}\n")
        f.write(f"{sequence}\n")


# Main function to process the pipeline
def main():
    # Prompt user for input
    input_choice = input(
        "Would you like to provide a FASTA file or a DNA sequence string? (Enter 'file' or 'string'): ").strip().lower()

    output_file = "output_results.txt"  # Output file for results
    with open(output_file, "w") as f:  # Clear the previous content of the output file
        f.write("Sequence Analysis Results\n")

    if input_choice == 'file':
        # Get the file path from user
        fasta_file = input("Enter the path to the FASTA file: ").strip()

        try:
            # Parse the FASTA file
            sequences = parse_fasta(fasta_file)
        except Exception as e:
            print(f"Error reading FASTA file: {e}")
            return
    elif input_choice == 'string':
        # Get the DNA sequence from user
        sequence = input("Enter the DNA sequence: ").strip().upper()
        sequence_name = input("Enter a name for your sequence: ").strip()
        sequences = {sequence_name: sequence}  # Store the string as a dictionary
        save_sequence_to_fasta(sequence, sequence_name)  # Save the sequence as a FASTA file
    else:
        print("Invalid input. Please enter either 'file' or 'string'.")
        return

    # Process each sequence
    for seq_id, sequence in sequences.items():
        print(f"> {seq_id}")
        print(f"Sequence: {sequence[:100]}...")  # Print the first 100 bases for inspection

        # 1. Find ORFs
        orfs = find_orfs(sequence)
        print(f"  ORFs (total: {len(orfs)}):")
        for orf in orfs:
            print(f"    - {orf}")

        # 2. Calculate GC content
        gc_content = calculate_gc_content(sequence)
        print(f"  GC Content: {gc_content}%")

        # 3. Find and sort restriction sites
        print("  Restriction Sites:")
        restriction_sites = find_restriction_sites(sequence, restriction_enzymes)

        # Filter and sort restriction sites by count
        sorted_sites = sorted(
            [(enzyme, count, positions) for enzyme, (count, positions) in restriction_sites.items() if count > 0],
            key=lambda x: x[1]
        )

        # Print sorted restriction sites
        for enzyme, count, positions in sorted_sites:
            print(f"    - {enzyme}: {count} site(s) at positions {positions}")

        # Write results to the output file
        write_results_to_file(output_file, seq_id, sequence, orfs, gc_content, sorted_sites)

        print("\n")

    print(f"Results have been written to {output_file}")

main()