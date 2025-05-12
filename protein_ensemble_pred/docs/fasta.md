## Supported Input Formats

This application supports the following input formats:

1.  **FASTA files (`.fasta`, `.fa`, etc.):** Ideal for simple inputs, especially when you want the application to handle Multiple Sequence Alignment (MSA) generation automatically. See detailed FASTA formatting guidelines below.
2.  **AlphaFold3 JSON (`.json`):** You can provide a pre-existing AlphaFold3 JSON file. This is useful if you have already run part of the AlphaFold3 pipeline (e.g., obtained an MSA-enriched JSON) or have a complex setup defined in AF3's native format. The system will use the sequences and any available MSA information from this file and will also convert it to a suitable format for Boltz-1 predictions.
3.  **Boltz-1 YAML (`.yaml`, `.yml`):** You can provide a Boltz-1 YAML input file. This is useful if you have specific constraints, modified residues, or other advanced features defined in Boltz-1's native format. The system will use the sequences and any available MSA/constraint information from this file and will also convert it to a suitable format for AlphaFold3 predictions.

The system is designed to internally manage and convert these different input types as needed to run predictions with both AlphaFold3 and Boltz-1. 

## FASTA Input Formatting Guidelines

When providing input as a FASTA file, the system will parse each entry and attempt to determine the molecule type based on the sequence content. Here are some important guidelines to ensure your FASTA file is processed correctly:


### Explicitly Specifying Molecule Types and MSAs in FASTA Headers

In addition to type guessing based on sequence content, you can provide more explicit information directly in the FASTA header using a pipe-separated (`|`) format. This is particularly useful for ensuring correct molecule typing, assigning specific chain IDs, and providing pre-computed Multiple Sequence Alignments (MSAs) for protein chains.

The recognized format is:

`>CHAIN_ID|ENTITY_TYPE|MSA_PATH`

Where:

*   **`CHAIN_ID`**: Your desired unique identifier for this chain (e.g., `PROTA`, `LIG1`, `RNA_R`). This ID will be used in the generated configuration files. If left empty but pipes are present (e.g. `>|protein|...`), a chain ID will be automatically generated.
*   **`ENTITY_TYPE`**: Specifies the type of the molecule. Supported values are (case-insensitive):
    *   `protein`: For a protein sequence.
    *   `rna`: For an RNA sequence.
    *   `dna`: For a DNA sequence.
    *   `smiles`: For a ligand defined by a SMILES string. The sequence line should contain the SMILES string.
    *   `ccd`: For a ligand defined by a PDB Chemical Component Dictionary (CCD) code. The sequence line should contain the CCD code (e.g., `ATP`).
*   **`MSA_PATH`**: (Optional, primarily for `protein` entity type)
    *   The full or relative path to a pre-computed Multiple Sequence Alignment file for this protein (e.g., in `.a3m` format).
    *   If you want to run this protein explicitly without an MSA (single sequence mode), use the special keyword `empty`.
    *   If this field is omitted or left blank for a protein, the system will attempt to generate an MSA automatically.
    *   This field is ignored if the determined `ENTITY_TYPE` is not `protein`.

**How it works with this format:**

*   If a FASTA header matches this pipe-separated format (i.e., contains `|` characters):
    *   The system will attempt to parse `CHAIN_ID`, `ENTITY_TYPE`, and `MSA_PATH`.
    *   The `CHAIN_ID` from the header will be used if provided and non-empty; otherwise, a new chain ID will be generated.
    *   If the parsed `ENTITY_TYPE` is one of the recognized values listed above (case-insensitive), it will be used directly. Sequence-based type guessing will be skipped for this entry.
    *   If the parsed `ENTITY_TYPE` is *not* one of the recognized values (or is missing), the system will fall back to sequence-based type guessing (as described in "Sequence Interpretation and Type Guessing") to determine the molecule type. However, the `CHAIN_ID` (either parsed from the header or generated if the header's `CHAIN_ID` part was empty) will still be used.
    *   If an `MSA_PATH` is provided and the *explicitly specified* `ENTITY_TYPE` in the header was `protein`, that MSA will be used. If type guessing was invoked for an entry (even if it guessed 'protein'), an `MSA_PATH` from the header for that entry will not be used automatically; MSA generation would proceed based on the guessed type.

**Examples using the explicit header format:**

1.  **Protein with a pre-computed MSA:**
    ```fasta
    >ChainA|protein|/path/to/my_msas/prota.a3m
    MPEPTIDESEQUENCE...
    ```

2.  **Protein to be run in single-sequence mode:**
    ```fasta
    >ChainB|protein|
    ANOTHERPEPTIDE...
    ```

3.  **RNA molecule:**
    ```fasta
    >RNASegment1|rna|
    AGCUAGCUAGCU...
    ```

4.  **Ligand specified by a SMILES string:**
    ```fasta
    >MyDrug|smiles|
    CC(=O)Oc1ccccc1C(=O)OH
    ```

5.  **Ligand specified by a CCD code:**
    ```fasta
    >CofactorX|ccd|
    NAD
    ```

By using this explicit header format, you gain precise control over how each sequence in your FASTA file is interpreted and prepared for the prediction models. If the header does not follow this pipe-separated format (i.e. no `|` characters), the system will fall back to the general type guessing and automatic chain ID generation described in "Sequence Interpretation and Type Guessing".

### Sequence Interpretation and Type Guessing

Each sequence in the FASTA file will be assigned a molecule type. The type guessing prioritizes as follows:

1.  **RNA:** Sequences composed strictly of `A`, `C`, `G`, `U` .
2.  **DNA:** Sequences composed strictly of `A`, `C`, `G`, `T` .
3.  **Protein:** Sequences composed of the standard 20 amino acid one-letter codes.
4.  **Ligand (CCD Code):** If a sequence is 1-3 characters long and consists of alphanumeric characters, it will be interpreted as a potential PDB Chemical Component Dictionary (CCD) code for a ligand (e.g., `ATP`, `HEM`).
5.  **Ligand (SMILES String):** If a sequence does not match any of the above and contains characters commonly found in SMILES strings (e.g., `(`, `)`, `=`, `#`, `[`, `]`, numbers, `+`, `-`, `@`), it will be heuristically interpreted as a SMILES string for a ligand.
6.  **Unknown:** If a sequence does not fit any of the above categories, it will be marked as `unknown`. These sequences wil be ignored and cause warnings during processing.

**Note:** The type guessing for ligands is heuristic. For precise ligand definition, especially for novel ligands not in standard CCDs, providing an AlphaFold3 JSON input with explicit SMILES or custom CCD definitions is recommended.

### Header Lines (`>`)

*   The entire line following the `>` symbol is captured as the original name or description of the sequence. This name is currently used for informational purposes and does **not** influence molecule type determination or stoichiometry.
*   Special header notations like `#protein`, `#dna`, `#ligand`, or `#3` (for counts) found in some other FASTA utilities are **not** currently used by this parser for type or stoichiometry. 

### Complexes

*   If your FASTA file contains multiple sequences, they will be treated as individual components of a single complex. Each sequence will be assigned a unique chain ID (e.g., A, B, C, ...).

### Homomeric Chains (Multiple Copies of the Same Sequence)

*   The current FASTA parser assigns a unique chain ID to each entry in the FASTA file. To model a homomer (e.g., a homotrimer where three copies of the same protein chain are present):
    *   **Option 1 (Recommended for FASTA input):** Include the sequence multiple times in your FASTA file. Each will be treated as a separate chain that happens to have the same sequence, and they will be assigned different chain IDs (e.g., A, B, C for three identical sequences).
    *   **Option 2 (Manual JSON edit):** After the AlphaFold3 JSON is generated from your FASTA, you would need to manually edit the JSON to represent homomeric relationships if you listed the sequence only once. For example, for a protein, you would change `"id": "A"` to `"id": ["A", "B", "C"]` for a trimer, ensuring these IDs are unique across the whole system. 