# ProbeMaker

A Python tool to generate complementary sequences to mRNA for given genes. This program can take gene sequences (DNA or RNA) as input and output complementary sequences that could be used as probes, primers, or for other molecular biology applications.

## Features

- **DNA to RNA Transcription**: Automatically transcribes DNA sequences to mRNA
- **Gene Name Input**: Accepts gene names and automatically fetches sequences from NCBI
- **Reverse Complement Generation**: Generates reverse complementary sequences to mRNA
- **Multiple Probe Design**: Generates 3 different 50-base single-stranded DNA probes per gene
- **LHS/RHS Probe Splitting**: Automatically splits each 50-base probe into two 25-base probes
- **Probe Design Constraints**:
  - Position 25 must be T (complementing to A in input RNA)
  - GC content must be between 44% and 72%
  - No homopolymer repeats >4 of the same nucleotide in a row
- **Multiple Input Methods**: 
  - Single sequence input
  - File-based batch processing
  - Gene name input (NCBI integration)
  - Interactive mode
- **Flexible Output**: Console output or save to file
- **Input Validation**: Ensures sequences contain valid nucleotides and are long enough
- **Support for Both DNA and RNA**: Can handle both input types

## Installation

1. Clone or download this repository
2. Ensure you have Python 3.6+ installed
3. Install required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Command Line Interface

#### Single Sequence Mode
```bash
# Generate complement for a DNA sequence
python probe_maker.py -s "ATGCGATCGATCG"

# Generate complement for an RNA sequence
python probe_maker.py -s "AUGCGAUCGAUCG" --rna

# Save output to file
python probe_maker.py -s "ATGCGATCGATCG" -o results.txt
```

#### File Mode
```bash
# Process multiple sequences from a file
python probe_maker.py -f genes.txt

# Process RNA sequences from file
python probe_maker.py -f genes.txt --rna
```

#### Gene Name Mode (NCBI Integration)
```bash
# Process gene names and fetch sequences from NCBI
python probe_maker.py -g gene_names.txt

# Save gene name results to file
python probe_maker.py -g gene_names.txt -o results.txt
```

#### Interactive Mode
```bash
# Run in interactive mode
python probe_maker.py --interactive
```

### Interactive Mode Commands

When running in interactive mode, you can use these commands:
- `help` - Show available commands
- `process` - Process all entered sequences
- `clear` - Clear all entered sequences
- `list` - Show all entered sequences
- `quit` - Exit the program

## Input Format

### DNA Sequences
- Use standard DNA nucleotides: A, T, G, C
- Case insensitive (a, t, g, c also accepted)
- Whitespace is automatically removed
- Example: `ATGCGATCGATCG`

### RNA Sequences
- Use standard RNA nucleotides: A, U, G, C
- Case insensitive (a, u, g, c also accepted)
- Use `--rna` flag to indicate RNA input
- Example: `AUGCGAUCGAUCG`

### File Input
- One sequence per line
- Empty lines are ignored
- Example file content:
```
ATGCGATCGATCG
GCTAGCTAGCTAG
TATATATATATA
```

### Gene Name Input
- One gene name per line
- Gene names are automatically searched in NCBI
- mRNA sequences are fetched and processed
- Example file content:
```
BRCA1
TP53
EGFR
KRAS
BRAF
```

## Output Format

The program outputs:
1. **Input Sequence**: The original sequence provided
2. **mRNA Sequence**: The transcribed mRNA (if DNA input) or cleaned RNA
3. **Reverse Complementary Sequence**: The reverse complementary sequence to the mRNA
4. **50-Base Probe**: Single-stranded DNA probe (exactly 50 bases, U converted to T)
5. **LHS Probe (25 bases)**: Left-hand side probe (first 25 bases of the 50-base probe)
6. **RHS Probe (25 bases)**: Right-hand side probe (last 25 bases of the 50-base probe)
7. **Probe Start Position**: Starting position in the reverse complement sequence
8. **GC Content**: Percentage of G and C nucleotides

**For Gene Name Input**: The program generates 3 different probe pairs per gene, each with different starting positions to provide variety in probe design.

Example output:
```
Gene 1:
  Input (DNA): ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
  mRNA: AUGCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCG
  Reverse Complementary: CGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGCAU
  50-Base Probe: TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
  LHS Probe (25 bases): TCGATCGATCGATCGATCGATCGAT
  RHS Probe (25 bases): CGATCGATCGATCGATCGATCGATC
  Probe Start Position: 3
  GC Content: 50.0%
```

## How It Works

1. **DNA Input**: 
   - Transcribes DNA to RNA (T → U)
   - Generates reverse complementary sequence to mRNA
   - Creates 50-base single-stranded DNA probe (U → T)

2. **RNA Input**:
   - Cleans and validates RNA sequence
   - Generates reverse complementary sequence directly
   - Creates 50-base single-stranded DNA probe (U → T)

3. **Reverse Complement Generation**:
   - A ↔ U (Adenine ↔ Uracil)
   - G ↔ C (Guanine ↔ Cytosine)
   - Sequence is then reversed

4. **50-Base Probe Creation**:
   - Takes first 50 bases of reverse complement
   - Converts U to T for single-stranded DNA
   - Ensures probe is exactly 50 bases long
   - **Constraint Validation**:
     - Position 25 must be T (complementing to A in input RNA)
     - GC content must be between 44% and 72%
     - No homopolymer repeats >4 of the same nucleotide in a row
     - If constraints cannot be met, tries different starting positions
     - Fails gracefully if no valid probe can be designed

5. **LHS/RHS Probe Splitting**:
   - Automatically splits the 50-base probe into two 25-base probes
   - **LHS Probe**: First 25 bases (positions 1-25)
   - **RHS Probe**: Last 25 bases (positions 26-50)
   - Both probes maintain all constraint validations
   - Useful for experimental design requiring shorter probe segments

## Examples

### Example 1: Single DNA Sequence
```bash
python probe_maker.py -s "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
```

Output:
```
Gene N/A:
  Input (DNA): ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
  mRNA: AUGCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCG
  Reverse Complementary: CGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGCAU
  50-Base Probe: CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
```

### Example 2: Batch Processing
Create a file `genes.txt`:
```
ATGCGATCGATCG
GCTAGCTAGCTAG
```

Run:
```bash
python probe_maker.py -f genes.txt
```

### Example 3: Gene Name Input (NCBI Integration)
Create a file `gene_names.txt`:
```
BRCA1
TP53
EGFR
```

Run:
```bash
python probe_maker.py -g gene_names.txt
```

This will automatically fetch sequences from NCBI and generate 3 probe pairs for each gene.

### Example 4: Interactive Mode
```bash
python probe_maker.py --interactive
```

Then enter sequences one by one and use commands to process them.

## Error Handling

The program includes comprehensive error handling for:
- Invalid nucleotide characters
- File not found errors
- Empty input files
- Invalid sequences

## Use Cases

- **Probe Design**: Generate 50-base single-stranded DNA probes for hybridization
- **LHS/RHS Probe Design**: Create two 25-base probes for flexible experimental design
- **Primer Design**: Create reverse primers for PCR
- **Research**: Analyze gene sequences and their reverse complements
- **Education**: Learn about DNA/RNA transcription, complementarity, and probe design
- **Molecular Biology**: Design oligonucleotide probes for various applications
- **Experimental Flexibility**: Use full 50-base probe or individual 25-base segments as needed

## NCBI Integration

When using gene name input (`-g` flag), the program automatically:
- Searches NCBI databases for gene sequences
- Fetches mRNA sequences for the specified genes
- Processes sequences to generate multiple probe pairs
- Implements rate limiting to respect NCBI server policies

**Note**: NCBI servers have rate limits. The program automatically adds delays between requests to avoid overwhelming the servers.

## Probe Design Constraints

The program enforces three critical constraints for optimal probe design:

1. **Position 25 Constraint**: The 25th nucleotide of the probe must be T
   - This ensures the probe has an A to complement to in the input RNA
   - Critical for proper hybridization and binding specificity

2. **GC Content Constraint**: GC content must be between 44% and 72%
   - Prevents probes that are too AT-rich (poor binding) or too GC-rich (high melting temperature)
   - Ensures optimal hybridization conditions

3. **Homopolymer Repeat Constraint**: No more than 4 of the same nucleotide in a row
   - Prevents probes with long homopolymer stretches (e.g., AAAAA, TTTTT, GGGGG, CCCCC)
   - Avoids poor hybridization and sequencing issues
   - Ensures probe quality and reliability

If a sequence cannot meet these constraints, the program:
- Tries different starting positions in the reverse complement
- Provides clear error messages for sequences that fail
- Continues processing other valid sequences

## Technical Details

- **Language**: Python 3.6+
- **Dependencies**: requests library for NCBI integration
- **Input Validation**: Regex-based nucleotide validation
- **Memory Efficient**: Processes sequences one at a time
- **Error Recovery**: Continues processing even if individual sequences fail
- **NCBI Integration**: Uses NCBI E-utilities for sequence fetching

## Contributing

Feel free to submit issues, feature requests, or pull requests to improve this tool.

## License

This project is open source and available under the MIT License.
