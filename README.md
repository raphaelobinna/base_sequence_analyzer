# FASTA Base Composition Analyzer

A Python tool for detecting sequences with unusual base compositions in FASTA files. This analyzer identifies sequences that deviate significantly from the dataset mean, making it useful for quality control, contamination detection, and identifying sequences with atypical nucleotide content.

## Features

- **Adaptive thresholds**: Uses percentage-based deviations from dataset mean rather than fixed cutoffs
- **Multiple detection criteria**: Identifies unusual GC%, AT%, high N content, and base composition skew
- **Genus-agnostic**: Automatically adapts to the base composition characteristics of your dataset
- **Comprehensive output**: Provides both summary statistics and detailed per-sequence analysis
- **Export capability**: Save results to TSV format for further analysis

## Installation

### Requirements

- Python 3.6 or higher
- No external dependencies (uses only Python standard library)

### Setup

```bash
# Clone or download the script
wget https://raw.githubusercontent.com/your-repo/fasta_analyzer.py
# or
git clone https://github.com/your-repo/fasta-analyzer.git

# Make executable (optional)
chmod +x fasta_analyzer.py
```

## Usage

### Basic Usage

```bash
python fasta_analyzer.py sequences.fasta
```

### Command Line Options

```bash
python fasta_analyzer.py [FASTA_FILE] [OPTIONS]

Positional Arguments:
  FASTA_FILE              Input FASTA file to analyze

Optional Arguments:
  -h, --help              Show help message and exit
  --gc-threshold FLOAT    Fractional deviation from mean GC% to flag as unusual
                          (default: 0.2 = ±20%)
  --at-threshold FLOAT    Fractional deviation from mean AT% to flag as unusual
                          (default: 0.2 = ±20%)
  --n-threshold FLOAT     Percentage of N bases to flag as unusual (default: 5.0)
  -o, --output FILE       Output file for detailed results in TSV format
```

### Examples

#### Standard analysis with default thresholds

```bash
python fasta_analyzer.py my_sequences.fasta
```

#### More stringent detection (10% deviation threshold)

```bash
python fasta_analyzer.py my_sequences.fasta --gc-threshold 0.1 --at-threshold 0.1
```

#### More lenient detection with high N-content tolerance

```bash
python fasta_analyzer.py my_sequences.fasta --gc-threshold 0.3 --n-threshold 15
```

#### Save detailed results to file

```bash
python fasta_analyzer.py my_sequences.fasta --output analysis_results.tsv
```

#### Complete analysis with custom parameters

```bash
python fasta_analyzer.py bacterial_genomes.fasta \
  --gc-threshold 0.15 \
  --at-threshold 0.15 \
  --n-threshold 2.0 \
  --output bacterial_analysis.tsv
```

## Detection Criteria

The analyzer flags sequences as unusual based on the following criteria:

### 1. GC Content Deviation

- **Method**: Percentage deviation from dataset mean
- **Default threshold**: ±20% from mean
- **Example**: If mean GC% = 45%, flags sequences <36% or >54%

### 2. AT Content Deviation

- **Method**: Percentage deviation from dataset mean
- **Default threshold**: ±20% from mean
- **Example**: If mean AT% = 55%, flags sequences <44% or >66%

### 3. High N Content

- **Method**: Direct percentage threshold
- **Default threshold**: >5% N bases
- **Purpose**: Identifies low-quality or ambiguous sequences

### 4. Base Composition Skew

- **Method**: Single base dominance
- **Threshold**: Any single base (A, T, G, or C) >60%
- **Purpose**: Detects highly biased sequences (e.g., poly-A tails, repeats)

## Output Format

### Console Output

```
FASTA BASE COMPOSITION ANALYSIS SUMMARY
========================================

Dataset Statistics:
  Total sequences analyzed: 1,247
  Mean GC content: 44.32% (±3.15)
  Mean AT content: 55.68% (±3.15)

Unusual Sequences Found: 23

Detailed Results:
------------------------------

1. sequence_001_suspicious
   Length: 2,847 bp
   Composition: A=15.2% T=18.9% G=31.1% C=34.8%
   GC content: 65.9%
   Flags: GC%: 65.9% (+48.7% from mean 44.3%)
   Preview: ATGGCGCGCGCGATCGATCGCGCGCGATATCGATCGCGCGCGC...
```

### TSV Output (--output option)

```
Header	Length	A%	T%	G%	C%	GC%	AT%	N%	Flags
seq1	1250	25.12	24.88	24.96	25.04	50.00	50.00	0.00	Normal
seq2	987	15.20	18.94	31.12	34.74	65.86	34.14	0.00	GC%: 65.9% (+48.7% from mean 44.3%)
```

## Understanding Thresholds

### GC/AT Thresholds

The `--gc-threshold` and `--at-threshold` parameters accept fractional values:

- `0.1` = ±10% deviation from mean
- `0.2` = ±20% deviation from mean (default)
- `0.3` = ±30% deviation from mean

**Example calculation:**

- Dataset mean GC% = 40%
- Threshold = 0.2 (20%)
- Lower bound = 40% × (1 - 0.2) = 32%
- Upper bound = 40% × (1 + 0.2) = 48%
- Sequences with GC% <32% or >48% are flagged

### Choosing Appropriate Thresholds

| Organism Type     | Suggested GC Threshold | Reasoning                                   |
| ----------------- | ---------------------- | ------------------------------------------- |
| Bacterial genomes | 0.15-0.20              | Relatively stable GC content within species |
| Viral sequences   | 0.25-0.30              | More variable, especially RNA viruses       |
| Mixed datasets    | 0.30-0.40              | Account for inter-species variation         |
| Quality control   | 0.10-0.15              | Strict detection of anomalies               |

## Use Cases

### 1. Quality Control

Identify potentially problematic sequences in genomic datasets:

```bash
python fasta_analyzer.py raw_sequences.fasta --gc-threshold 0.1 --n-threshold 2
```

### 2. Contamination Detection

Find sequences that don't match expected organism characteristics:

```bash
python fasta_analyzer.py bacterial_assembly.fasta --gc-threshold 0.2
```

### 3. Sequence Classification

Pre-screen sequences before phylogenetic analysis:

```bash
python fasta_analyzer.py 16S_sequences.fasta --output classification_prep.tsv
```

### 4. Dataset Characterization

Understand the composition diversity in your dataset:

```bash
python fasta_analyzer.py metagenome_contigs.fasta --gc-threshold 0.5
```

## Troubleshooting

### Common Issues

#### Empty or No Results

```
No sequences found in the FASTA file.
```

**Solution**: Check file format, ensure sequences are properly formatted with `>` headers

#### File Not Found

```
Error: Input file 'sequences.fasta' does not exist.
```

**Solution**: Verify file path and name, ensure file exists in specified location

#### All Sequences Flagged

If most sequences are flagged as unusual, consider:

- Increasing thresholds (e.g., `--gc-threshold 0.3`)
- Checking if dataset contains mixed organism types
- Verifying sequence quality

#### No Sequences Flagged

If no sequences are flagged but you expect some:

- Decrease thresholds (e.g., `--gc-threshold 0.1`)
- Check N-content threshold
- Verify input sequences are of expected quality

### Performance Notes

- Memory usage: ~1KB per sequence (scales linearly)
- Processing speed: ~10,000-50,000 sequences/second (depends on sequence length)
- File size limit: No hard limit, but very large files (>1M sequences) may take several minutes

## Output Interpretation

### Normal Sequences

Sequences within expected parameters show:

- `Flags: Normal` in TSV output
- Not listed in console "Unusual Sequences" section

### Flagged Sequences

Each flag indicates a specific issue:

- `GC%: X% (+Y% from mean Z%)`: GC content deviation
- `AT%: X% (+Y% from mean Z%)`: AT content deviation
- `High N content: X%`: Too many ambiguous bases
- `Base skew detected (max: X%)`: Single base dominance

### Statistical Context

The dataset statistics help interpret individual sequence flags:

- **Mean values**: Expected composition for your dataset
- **Standard deviations**: Natural variation in your dataset
- **Total sequences**: Sample size for statistical confidence

## Contributing

Suggestions and improvements welcome! Common enhancement requests:

- Additional base composition metrics
- Support for amino acid sequences
- Integration with sequence databases
- Batch processing capabilities

## License

This tool is provided as-is for research and educational purposes. Feel free to modify and distribute according to your needs.

## Citation

If you use this tool in your research, please cite:

```
FASTA Base Composition Analyzer [Computer software].
Retrieved from https://github.com/raphaelobinna/fasta-analyzer
```
