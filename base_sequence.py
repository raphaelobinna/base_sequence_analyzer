#!/usr/bin/env python3
"""
FASTA Base Composition Analyzer
Detects sequences with unusual base compositions in FASTA files.
"""

import argparse
import sys
from collections import Counter
from pathlib import Path
import statistics

class BaseCompositionAnalyzer:
    def __init__(self, gc_threshold=0.2, at_threshold=0.2, n_threshold=5.0):
        """
        Initialize the analyzer with thresholds for unusual composition detection.
        
        Args:
            gc_threshold: Fractional deviation from mean GC% to flag as unusual (e.g., 0.2 = ±20%)
            at_threshold: Fractional deviation from mean AT% to flag as unusual (e.g., 0.2 = ±20%)
            n_threshold: Percentage of N's to flag as unusual
        """
        self.gc_threshold = gc_threshold
        self.at_threshold = at_threshold
        self.n_threshold = n_threshold
        self.sequences = []
        
    def parse_fasta(self, filepath):
        """Parse FASTA file and store sequences."""
        current_header = None
        current_seq = []
        
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # Save previous sequence if exists
                        if current_header:
                            self.sequences.append({
                                'header': current_header,
                                'sequence': ''.join(current_seq).upper()
                            })
                        # Start new sequence
                        current_header = line[1:]  # Remove '>'
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                # Don't forget the last sequence
                if current_header:
                    self.sequences.append({
                        'header': current_header,
                        'sequence': ''.join(current_seq).upper()
                    })
                    
        except FileNotFoundError:
            print(f"Error: File '{filepath}' not found.")
            sys.exit(1)
        except Exception as e:
            print(f"Error reading file: {e}")
            sys.exit(1)
            
        print(f"Parsed {len(self.sequences)} sequences from {filepath}")
        
    def calculate_base_composition(self, sequence):
        """Calculate base composition statistics for a sequence."""
        # Count bases
        base_counts = Counter(sequence)
        total_length = len(sequence)
        
        if total_length == 0:
            return None
            
        # Calculate percentages
        composition = {
            'A': (base_counts.get('A', 0) / total_length) * 100,
            'T': (base_counts.get('T', 0) / total_length) * 100,
            'G': (base_counts.get('G', 0) / total_length) * 100,
            'C': (base_counts.get('C', 0) / total_length) * 100,
            'N': (base_counts.get('N', 0) / total_length) * 100,
            'other': ((total_length - sum(base_counts[base] for base in 'ATGCN' if base in base_counts)) / total_length) * 100
        }
        
        # Calculate derived metrics
        composition['GC'] = composition['G'] + composition['C']
        composition['AT'] = composition['A'] + composition['T']
        composition['length'] = total_length
        composition['counts'] = dict(base_counts)
        
        return composition
        
    def analyze_all_sequences(self):
        """Analyze all sequences and identify unusual compositions."""
        results = []
        gc_percentages = []
        at_percentages = []
        
        # First pass: calculate compositions
        for seq_data in self.sequences:
            comp = self.calculate_base_composition(seq_data['sequence'])
            if comp:
                results.append({
                    'header': seq_data['header'],
                    'composition': comp,
                    'sequence_preview': seq_data['sequence'][:50] + ('...' if len(seq_data['sequence']) > 50 else '')
                })
                gc_percentages.append(comp['GC'])
                at_percentages.append(comp['AT'])
        
        # Calculate statistics for the dataset
        if gc_percentages:
            gc_mean = statistics.mean(gc_percentages)
            gc_stdev = statistics.stdev(gc_percentages) if len(gc_percentages) > 1 else 0
            at_mean = statistics.mean(at_percentages)
            at_stdev = statistics.stdev(at_percentages) if len(at_percentages) > 1 else 0
        else:
            gc_mean = gc_stdev = at_mean = at_stdev = 0
            
        # Second pass: identify unusual sequences
        unusual_sequences = []
        for result in results:
            comp = result['composition']
            flags = []
            
            # Check GC content using percentage deviation from mean
            if gc_mean > 0:  # Avoid division by zero
                gc_lower_bound = gc_mean * (1 - self.gc_threshold)
                gc_upper_bound = gc_mean * (1 + self.gc_threshold)
                if comp['GC'] < gc_lower_bound or comp['GC'] > gc_upper_bound:
                    deviation_pct = ((comp['GC'] - gc_mean) / gc_mean) * 100
                    flags.append(f"GC%: {comp['GC']:.1f}% ({deviation_pct:+.1f}% from mean {gc_mean:.1f}%)")
            
            # Check AT content using percentage deviation from mean
            if at_mean > 0:  # Avoid division by zero
                at_lower_bound = at_mean * (1 - self.at_threshold)
                at_upper_bound = at_mean * (1 + self.at_threshold)
                if comp['AT'] < at_lower_bound or comp['AT'] > at_upper_bound:
                    deviation_pct = ((comp['AT'] - at_mean) / at_mean) * 100
                    flags.append(f"AT%: {comp['AT']:.1f}% ({deviation_pct:+.1f}% from mean {at_mean:.1f}%)")
            
            # Check N content
            if comp['N'] > self.n_threshold:
                flags.append(f"High N content: {comp['N']:.1f}%")
                
            # Check for very skewed individual bases
            max_single_base = max(comp['A'], comp['T'], comp['G'], comp['C'])
            if max_single_base > 60:  # If any single base is >60%
                flags.append(f"Base skew detected (max: {max_single_base:.1f}%)")
            
            if flags:
                result['flags'] = flags
                unusual_sequences.append(result)
                
        return {
            'all_results': results,
            'unusual_sequences': unusual_sequences,
            'dataset_stats': {
                'total_sequences': len(results),
                'gc_mean': gc_mean,
                'gc_stdev': gc_stdev,
                'at_mean': at_mean,
                'at_stdev': at_stdev
            }
        }
    
    def print_summary(self, analysis_results):
        """Print a summary of the analysis."""
        stats = analysis_results['dataset_stats']
        unusual = analysis_results['unusual_sequences']
        
        print("\n" + "="*60)
        print("FASTA BASE COMPOSITION ANALYSIS SUMMARY")
        print("="*60)
        
        print(f"\nDataset Statistics:")
        print(f"  Total sequences analyzed: {stats['total_sequences']}")
        print(f"  Mean GC content: {stats['gc_mean']:.2f}% (±{stats['gc_stdev']:.2f})")
        print(f"  Mean AT content: {stats['at_mean']:.2f}% (±{stats['at_stdev']:.2f})")
        
        print(f"\nUnusual Sequences Found: {len(unusual)}")
        
        if unusual:
            print(f"\nDetailed Results:")
            print("-" * 60)
            
            for i, seq in enumerate(unusual, 1):
                comp = seq['composition']
                print(f"\n{i}. {seq['header'][:80]}{'...' if len(seq['header']) > 80 else ''}")
                print(f"   Length: {comp['length']} bp")
                print(f"   Composition: A={comp['A']:.1f}% T={comp['T']:.1f}% G={comp['G']:.1f}% C={comp['C']:.1f}%")
                print(f"   GC content: {comp['GC']:.1f}%")
                if comp['N'] > 0:
                    print(f"   N content: {comp['N']:.1f}%")
                print(f"   Flags: {'; '.join(seq['flags'])}")
                print(f"   Preview: {seq['sequence_preview']}")
        else:
            print("No sequences with unusual base composition detected.")
            
    def save_results(self, analysis_results, output_file):
        """Save detailed results to a file."""
        with open(output_file, 'w') as f:
            f.write("Header\tLength\tA%\tT%\tG%\tC%\tGC%\tAT%\tN%\tFlags\n")
            
            for seq in analysis_results['all_results']:
                comp = seq['composition']
                flags = '; '.join(seq.get('flags', ['Normal']))
                
                f.write(f"{seq['header']}\t"
                       f"{comp['length']}\t"
                       f"{comp['A']:.2f}\t"
                       f"{comp['T']:.2f}\t"
                       f"{comp['G']:.2f}\t"
                       f"{comp['C']:.2f}\t"
                       f"{comp['GC']:.2f}\t"
                       f"{comp['AT']:.2f}\t"
                       f"{comp['N']:.2f}\t"
                       f"{flags}\n")
        
        print(f"\nDetailed results saved to: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Analyze FASTA files for unusual base compositions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python fasta_analyzer.py sequences.fasta
  python fasta_analyzer.py sequences.fasta --gc-threshold 1.5 --output results.tsv
  python fasta_analyzer.py sequences.fasta --n-threshold 10 --at-threshold 3.0
        """)
    
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('--gc-threshold', type=float, default=0.2,
                       help='Fractional deviation from mean GC%% to flag as unusual (default: 0.2 = ±20%%)')
    parser.add_argument('--at-threshold', type=float, default=0.2,
                       help='Fractional deviation from mean AT%% to flag as unusual (default: 0.2 = ±20%%)')
    parser.add_argument('--n-threshold', type=float, default=5.0,
                       help='Percentage of N bases to flag as unusual (default: 5.0)')
    parser.add_argument('--output', '-o', 
                       help='Output file for detailed results (TSV format)')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not Path(args.fasta_file).exists():
        print(f"Error: Input file '{args.fasta_file}' does not exist.")
        sys.exit(1)
    
    # Initialize analyzer
    analyzer = BaseCompositionAnalyzer(
        gc_threshold=args.gc_threshold,
        at_threshold=args.at_threshold,
        n_threshold=args.n_threshold
    )
    
    # Parse FASTA and analyze
    analyzer.parse_fasta(args.fasta_file)
    
    if not analyzer.sequences:
        print("No sequences found in the FASTA file.")
        sys.exit(1)
        
    results = analyzer.analyze_all_sequences()
    
    # Print summary
    analyzer.print_summary(results)
    
    # Save detailed results if requested
    if args.output:
        analyzer.save_results(results, args.output)

if __name__ == "__main__":
    main()