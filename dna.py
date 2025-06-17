#DNA ANALYZER


def count_nucleotides(sequence):
    sequence = sequence.upper()
    out_dict = {}
    for s in sequence:
        if s in out_dict:
            out_dict[s] += 1
        else:
            out_dict[s] = 1

    return out_dict

def gc_content(sequence):

    sequence = sequence.upper()
    count=0
    g_count=0
    c_count=0
    for s in sequence:
        count +=1
        if s == 'G':
            g_count +=1
        if s == 'C':
            c_count +=1

    if sequence == '':
        count = 1
    return ((g_count + c_count) / count) * 100

def reverse_complement(sequence):
    sequence = sequence.upper()
    complement = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    new_sequence = ""
    for s in sequence:
        new_sequence += complement[s]

    return new_sequence[::-1]


def transcribe_to_rna(sequence):
    sequence = sequence.upper()
    return sequence.upper().replace('T','U')

def translate_rna_to_protein(rna_sequence):
    rna_sequence = rna_sequence.upper()
    codon_table = {
        'AUG': 'M',
        'UUU': 'F', 'UUC': 'F',
        'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y',
        'CAU': 'H', 'CAC': 'H',
        'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N',
        'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D',
        'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C',
        'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S',
        'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
    }

    protein = ""
    start = False

    for i in range(0, len(rna_sequence)-2,3):
        try:
            codon = rna_sequence[i:i+3]
        except:
            break
        if start:
            if codon_table[codon]=='Stop':
                break
            elif codon in codon_table:
                protein += codon_table[codon]
            else:
                break
        else:
            if codon == "AUG":
                start = True
                protein += codon_table[codon]


    return protein



if __name__ == "__main__":

    test_dna_1 = "ATGCAGGTCCGAT"
    print(f"\nTesting DNA Sequence 1: {test_dna_1}")
    print(f"  Nucleotide Counts: {count_nucleotides(test_dna_1)}")
    print(f"  GC Content: {gc_content(test_dna_1):.2f}%")
    print(f"  Reverse Complement: {reverse_complement(test_dna_1)}")
    rna_seq_1 = transcribe_to_rna(test_dna_1)
    print(f"  Transcribed RNA: {rna_seq_1}")
    protein_1 = translate_rna_to_protein(rna_seq_1)
    print(f"  Translated Protein: {protein_1}")

    test_dna_2 = "GATTACA"
    print(f"\nTesting DNA Sequence 2: {test_dna_2}")
    print(f"  Nucleotide Counts: {count_nucleotides(test_dna_2)}")
    print(f"  GC Content: {gc_content(test_dna_2):.2f}%")
    print(f"  Reverse Complement: {reverse_complement(test_dna_2)}")
    rna_seq_2 = transcribe_to_rna(test_dna_2)
    print(f"  Transcribed RNA: {rna_seq_2}")
    protein_2 = translate_rna_to_protein(rna_seq_2)
    print(f"  Translated Protein: {protein_2}") # Should be empty as no AUG

    test_dna_3 = ""
    print(f"\nTesting Empty Sequence: '{test_dna_3}'")
    print(f"  GC Content: {gc_content(test_dna_3):.2f}%")

    test_dna_4 = "agcttgacgt" # Test mixed case
    print(f"\nTesting Mixed Case Sequence: {test_dna_4}")
    print(f"  Nucleotide Counts: {count_nucleotides(test_dna_4)}")
    print(f"  Reverse Complement: {reverse_complement(test_dna_4)}")

    test_rna_for_translation = "AUGGCCAUGUAAAGCCUGA" # Should translate M A M *
    print(f"\nTesting RNA for Translation: {test_rna_for_translation}")
    protein_translation_test = translate_rna_to_protein(test_rna_for_translation)
    print(f"  Translated Protein: {protein_translation_test}")

    test_rna_no_start = "UUGGCAUGCAGUAA" # No AUG at the beginning
    print(f"\nTesting RNA with no leading AUG: {test_rna_no_start}")
    protein_no_start = translate_rna_to_protein(test_rna_no_start)
    print(f"  Translated Protein: {protein_no_start}") # Should be empty

    test_rna_incomplete_codon = "AUGGC" # Incomplete final codon
    print(f"\nTesting RNA with incomplete final codon: {test_rna_incomplete_codon}")
    protein_incomplete = translate_rna_to_protein(test_rna_incomplete_codon)
    print(f"  Translated Protein: {protein_incomplete}") # Should be 'M'




