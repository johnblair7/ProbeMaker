#!/usr/bin/env python3
"""
ProbeMaker - A tool to generate complementary sequences to mRNA for given genes.
This program can take gene sequences or gene names and output complementary sequences
that could be used as probes or primers.
"""

import argparse
import sys
import requests
import time
from typing import List, Dict, Optional
import re


class GeneSequenceFetcher:
    """Class for fetching gene sequences from online databases."""
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'ProbeMaker/1.0 (Research Tool)'
        })
    
    def fetch_gene_sequence(self, gene_name: str) -> Optional[str]:
        """Fetch mRNA sequence for a given gene name from NCBI."""
        try:
            # Search for the gene in NCBI
            search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            search_params = {
                'db': 'nucleotide',
                'term': f"{gene_name}[Gene Name] AND (transcript[Title] OR mRNA[Title] OR RNA[Title])",
                'retmode': 'json',
                'retmax': 10
            }
            
            response = self.session.get(search_url, params=search_params)
            response.raise_for_status()
            
            # Parse the search results
            search_data = response.json()
            if 'esearchresult' not in search_data or 'idlist' not in search_data['esearchresult']:
                return None
            
            gene_ids = search_data['esearchresult']['idlist']
            if not gene_ids:
                return None
            
            # Try multiple results to find a good sequence
            for gene_id in gene_ids[:5]:  # Try up to 5 results
                try:
                    fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                    fetch_params = {
                        'db': 'nucleotide',
                        'id': gene_id,
                        'rettype': 'fasta',
                        'retmode': 'text'
                    }
                    
                    response = self.session.get(fetch_url, params=fetch_params)
                    response.raise_for_status()
                    
                    # Parse FASTA format and extract sequence
                    lines = response.text.strip().split('\n')
                    sequence = ''.join(lines[1:])  # Skip header line
                    
                    # Clean the sequence (remove any non-nucleotide characters)
                    sequence = re.sub(r'[^ATGCUatgcu]', '', sequence)
                    
                    # Check if this sequence is long enough and looks like a good candidate
                    if len(sequence) >= 50 and self._is_good_sequence_candidate(sequence):
                        return sequence
                        
                except Exception as e:
                    # If this result fails, try the next one
                    continue
            
            # If no good sequence found, try a broader search without transcript filter
            print(f"  Trying broader search for {gene_name}...")
            fallback_search_params = {
                'db': 'nucleotide',
                'term': f"{gene_name}[Gene Name]",
                'retmode': 'json',
                'retmax': 5
            }
            
            response = self.session.get(search_url, params=fallback_search_params)
            response.raise_for_status()
            
            fallback_data = response.json()
            if 'esearchresult' in fallback_data and 'idlist' in fallback_data['esearchresult']:
                fallback_ids = fallback_data['esearchresult']['idlist']
                
                for gene_id in fallback_ids[:3]:
                    try:
                        fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                        fetch_params = {
                            'db': 'nucleotide',
                            'id': gene_id,
                            'rettype': 'fasta',
                            'retmode': 'text'
                        }
                        
                        response = self.session.get(fetch_url, params=fetch_params)
                        response.raise_for_status()
                        
                        lines = response.text.strip().split('\n')
                        sequence = ''.join(lines[1:])
                        sequence = re.sub(r'[^ATGCUatgcu]', '', sequence)
                        
                        if len(sequence) >= 50 and self._is_good_sequence_candidate(sequence):
                            return sequence
                            
                    except Exception as e:
                        continue
            
            return None
            
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 429:
                print(f"Rate limited for {gene_name}. Waiting longer...")
                time.sleep(2.0)  # Wait longer for rate limited requests
                return None
            else:
                print(f"HTTP error fetching sequence for {gene_name}: {e}")
                return None
        except Exception as e:
            print(f"Error fetching sequence for {gene_name}: {e}")
            return None
    
    def _is_good_sequence_candidate(self, sequence: str) -> bool:
        """Check if a sequence is a good candidate for probe design."""
        # Must be at least 50 bases long
        if len(sequence) < 50:
            return False
        
        # Reject extremely long sequences (likely genomic DNA)
        if len(sequence) > 100000:  # 100kb max
            return False
        
        # Should have a reasonable GC content (not too extreme)
        gc_count = sequence.upper().count('G') + sequence.upper().count('C')
        gc_percentage = (gc_count / len(sequence)) * 100
        
        # Accept sequences with GC content between 20% and 80%
        if gc_percentage < 20 or gc_percentage > 80:
            return False
        
        # Should not be mostly N's or other ambiguous characters
        ambiguous_count = sequence.upper().count('N') + sequence.upper().count('X')
        if ambiguous_count > len(sequence) * 0.1:  # More than 10% ambiguous
            return False
        
        return True
    
    def close(self):
        """Close the session."""
        self.session.close()


class ProbeMaker:
    """Main class for generating complementary sequences to mRNA."""
    
    def __init__(self):
        # DNA to RNA transcription rules
        self.dna_to_rna = {
            'A': 'U',  # Adenine -> Uracil
            'T': 'A',  # Thymine -> Adenine
            'G': 'C',  # Guanine -> Cytosine
            'C': 'G',  # Cytosine -> Guanine
            'a': 'u',  # Lowercase versions
            't': 'a',
            'g': 'c',
            'c': 'g'
        }
        
        # RNA to complementary RNA rules
        self.rna_complement = {
            'A': 'U',  # Adenine -> Uracil
            'U': 'A',  # Uracil -> Adenine
            'G': 'C',  # Guanine -> Cytosine
            'C': 'G',  # Cytosine -> Guanine
            'a': 'u',  # Lowercase versions
            'u': 'a',
            'g': 'c',
            'c': 'g'
        }
    
    def clean_sequence(self, sequence: str) -> str:
        """Clean and validate a DNA/RNA sequence."""
        # Remove whitespace and convert to uppercase
        cleaned = re.sub(r'\s+', '', sequence.upper())
        
        # Check for empty sequence
        if not cleaned:
            raise ValueError("Empty sequence provided")
        
        # Validate that sequence only contains valid nucleotides
        valid_nucleotides = set('ATGCU')
        if not all(nuc in valid_nucleotides for nuc in cleaned):
            raise ValueError(f"Invalid nucleotide found in sequence: {sequence}")
        
        return cleaned
    
    def dna_to_rna_transcribe(self, dna_sequence: str) -> str:
        """Transcribe DNA sequence to RNA sequence."""
        cleaned_dna = self.clean_sequence(dna_sequence)
        rna_sequence = ''
        
        for nucleotide in cleaned_dna:
            if nucleotide == 'T':
                rna_sequence += 'U'
            else:
                rna_sequence += nucleotide
        
        return rna_sequence
    
    def generate_reverse_complementary_sequence(self, sequence: str) -> str:
        """Generate reverse complementary sequence to the given RNA sequence."""
        cleaned_seq = self.clean_sequence(sequence)
        complementary = ''
        
        # Generate complementary sequence
        for nucleotide in cleaned_seq:
            if nucleotide in self.rna_complement:
                complementary += self.rna_complement[nucleotide]
            else:
                # Handle DNA nucleotides if present
                if nucleotide == 'T':
                    complementary += 'A'
                else:
                    complementary += nucleotide
        
        # Reverse the complementary sequence
        reverse_complement = complementary[::-1]
        
        return reverse_complement
    
    def create_50_base_probe(self, sequence: str) -> str:
        """Create a 50-base probe from the given sequence."""
        if len(sequence) < 50:
            raise ValueError(f"Sequence too short ({len(sequence)} bases) to create 50-base probe")
        
        # Take the first 50 bases of the sequence
        probe_50 = sequence[:50]
        
        # Convert U to T to make it single-stranded DNA
        probe_50_dna = probe_50.replace('U', 'T').replace('u', 't')
        
        return probe_50_dna
    
    def generate_multiple_probe_pairs(self, input_sequence: str, is_dna: bool = True, num_pairs: int = 3) -> List[Dict[str, str]]:
        """Generate multiple probe pairs for a single gene sequence."""
        try:
            if is_dna:
                # Transcribe DNA to RNA first
                mrna = self.dna_to_rna_transcribe(input_sequence)
            else:
                # Assume input is already RNA
                mrna = self.clean_sequence(input_sequence)
            
            # Generate reverse complementary sequence
            reverse_complement = self.generate_reverse_complementary_sequence(mrna)
            
            probe_pairs = []
            
            # Try to find non-overlapping probes
            used_positions = set()  # Track used starting positions to avoid overlap
            min_spacing = 50  # Minimum spacing between probe starts to ensure no overlap
            
            for pair_num in range(num_pairs):
                try:
                    # Find valid 50-base probe that doesn't overlap with previous ones
                    probe_50, probe_start_pos = self.find_non_overlapping_probe_start(
                        reverse_complement, mrna, used_positions, min_spacing
                    )
                    
                    # Mark this position as used
                    used_positions.add(probe_start_pos)
                    
                    # Split the 50-base probe into two 25-base probes
                    lhs_probe = probe_50[:25]
                    rhs_probe = probe_50[25:]
                    
                    # Calculate GC content for the probe
                    gc_count = probe_50.count('G') + probe_50.count('C')
                    gc_percentage = (gc_count / len(probe_50)) * 100
                    
                    probe_pair = {
                        'pair_number': pair_num + 1,
                        'input_sequence': input_sequence,
                        'mrna_sequence': mrna,
                        'reverse_complementary_sequence': reverse_complement,
                        'probe_50_base': probe_50,
                        'lhs_probe': lhs_probe,
                        'rhs_probe': rhs_probe,
                        'probe_start_position': probe_start_pos,
                        'gc_content_percentage': gc_percentage,
                        'input_type': 'DNA' if is_dna else 'RNA'
                    }
                    
                    probe_pairs.append(probe_pair)
                    
                except ValueError as e:
                    # If we can't find a non-overlapping probe, try with reduced spacing
                    min_spacing = max(25, min_spacing - 10)  # Reduce spacing requirement
                    try:
                        probe_50, probe_start_pos = self.find_non_overlapping_probe_start(
                            reverse_complement, mrna, used_positions, min_spacing
                        )
                        used_positions.add(probe_start_pos)
                        
                        lhs_probe = probe_50[:25]
                        rhs_probe = probe_50[25:]
                        gc_count = probe_50.count('G') + probe_50.count('C')
                        gc_percentage = (gc_count / len(probe_50)) * 100
                        
                        probe_pair = {
                            'pair_number': pair_num + 1,
                            'input_sequence': input_sequence,
                            'mrna_sequence': mrna,
                            'reverse_complementary_sequence': reverse_complement,
                            'probe_50_base': probe_50,
                            'lhs_probe': lhs_probe,
                            'rhs_probe': rhs_probe,
                            'probe_start_position': probe_start_pos,
                            'gc_content_percentage': gc_percentage,
                            'input_type': 'DNA' if is_dna else 'RNA'
                        }
                        
                        probe_pairs.append(probe_pair)
                        
                    except ValueError:
                        # If still can't find one, skip this pair
                        continue
            
            if not probe_pairs:
                raise ValueError("Could not generate any valid probe pairs")
            
            return probe_pairs
            
        except Exception as e:
            raise ValueError(f"Error generating probe pairs: {e}")
    
    def validate_probe_constraints(self, probe: str, original_mrna: str) -> tuple[bool, str]:
        """Validate that the probe meets all constraints."""
        # Check position 25 constraint (0-indexed, so position 24)
        if len(probe) < 25:
            return False, "Probe too short to check position 25 constraint"
        
        if probe[24] != 'T':
            return False, f"Position 25 must be T, but got {probe[24]}"
        
        # Check GC content constraint
        gc_count = probe.count('G') + probe.count('C')
        gc_percentage = (gc_count / len(probe)) * 100
        
        if gc_percentage < 44 or gc_percentage > 72:
            return False, f"GC content {gc_percentage:.1f}% is outside allowed range (44%-72%)"
        
        # Check homopolymer repeat constraint (no more than 4 of the same nucleotide in a row)
        for i in range(len(probe) - 3):  # Check 4-nucleotide windows
            window = probe[i:i+4]
            if len(set(window)) == 1:  # All nucleotides in window are the same
                return False, f"Homopolymer repeat found: {window} at position {i+1}-{i+4}"
        
        return True, f"Valid probe: GC content {gc_percentage:.1f}%, position 25 is T, no homopolymer repeats"
    
    def find_valid_probe_start(self, reverse_complement: str, original_mrna: str, start_offset: int = 0) -> tuple[str, int]:
        """Find a valid starting position for a 50-base probe that meets all constraints."""
        if len(reverse_complement) < 50:
            raise ValueError(f"Reverse complement too short ({len(reverse_complement)} bases) to create 50-base probe")
        
        # Start from the offset position to find different probes
        start_range = range(start_offset, len(reverse_complement) - 49)
        
        # Try different starting positions to find a valid probe
        for start_pos in start_range:
            if start_pos < 0 or start_pos + 50 > len(reverse_complement):
                continue
                
            probe_candidate = reverse_complement[start_pos:start_pos + 50]
            probe_dna = probe_candidate.replace('U', 'T').replace('u', 't')
            
            is_valid, message = self.validate_probe_constraints(probe_dna, original_mrna)
            if is_valid:
                return probe_dna, start_pos
        
        # If no valid probe found, raise an error
        raise ValueError("No valid 50-base probe found that meets all constraints")
    
    def find_non_overlapping_probe_start(self, reverse_complement: str, original_mrna: str, used_positions: set, min_spacing: int) -> tuple[str, int]:
        """Find a valid starting position for a 50-base probe that doesn't overlap with previous probes."""
        if len(reverse_complement) < 50:
            raise ValueError(f"Reverse complement too short ({len(reverse_complement)} bases) to create 50-base probe")
        
        # Try different starting positions to find a valid, non-overlapping probe
        for start_pos in range(0, len(reverse_complement) - 49):
            # Check if this position is too close to any previously used position
            too_close = False
            for used_pos in used_positions:
                if abs(start_pos - used_pos) < min_spacing:
                    too_close = True
                    break
            
            if too_close:
                continue
                
            probe_candidate = reverse_complement[start_pos:start_pos + 50]
            probe_dna = probe_candidate.replace('U', 'T').replace('u', 't')
            
            is_valid, message = self.validate_probe_constraints(probe_dna, original_mrna)
            if is_valid:
                return probe_dna, start_pos
        
        # If no non-overlapping probe found, raise an error
        raise ValueError(f"No valid non-overlapping 50-base probe found with minimum spacing {min_spacing}")
    
    def generate_mrna_complement(self, input_sequence: str, is_dna: bool = True) -> Dict[str, str]:
        """Generate complementary sequence to mRNA from DNA or RNA input."""
        try:
            if is_dna:
                # Transcribe DNA to RNA first
                mrna = self.dna_to_rna_transcribe(input_sequence)
            else:
                # Assume input is already RNA
                mrna = self.clean_sequence(input_sequence)
            
            # Generate reverse complementary sequence
            reverse_complement = self.generate_reverse_complementary_sequence(mrna)
            
            # Find valid 50-base probe that meets all constraints
            probe_50, probe_start_pos = self.find_valid_probe_start(reverse_complement, mrna)
            
            # Split the 50-base probe into two 25-base probes
            lhs_probe = probe_50[:25]
            rhs_probe = probe_50[25:]
            
            # Calculate GC content for the probe
            gc_count = probe_50.count('G') + probe_50.count('C')
            gc_percentage = (gc_count / len(probe_50)) * 100
            
            return {
                'input_sequence': input_sequence,
                'mrna_sequence': mrna,
                'reverse_complementary_sequence': reverse_complement,
                'probe_50_base': probe_50,
                'lhs_probe': lhs_probe,
                'rhs_probe': rhs_probe,
                'probe_start_position': probe_start_pos,
                'gc_content_percentage': gc_percentage,
                'input_type': 'DNA' if is_dna else 'RNA'
            }
        
        except ValueError as e:
            raise ValueError(f"Error processing sequence: {e}")
    
    def process_gene_list(self, gene_list: List[str], is_dna: bool = True) -> List[Dict[str, str]]:
        """Process a list of genes and generate complementary sequences."""
        results = []
        
        for i, gene in enumerate(gene_list, 1):
            try:
                result = self.generate_mrna_complement(gene, is_dna)
                result['gene_number'] = i
                results.append(result)
            except ValueError as e:
                print(f"Warning: Could not process gene {i} ('{gene}'): {e}")
                continue
        
        return results
    
    def process_gene_names(self, gene_names: List[str], sequence_fetcher: GeneSequenceFetcher) -> List[Dict[str, str]]:
        """Process a list of gene names, fetch sequences, and generate multiple probe pairs."""
        all_results = []
        
        for i, gene_name in enumerate(gene_names, 1):
            print(f"Processing gene {i}: {gene_name}")
            
            try:
                # Fetch the sequence for this gene
                sequence = sequence_fetcher.fetch_gene_sequence(gene_name)
                if not sequence:
                    print(f"Warning: Could not fetch sequence for gene '{gene_name}'")
                    continue
                
                print(f"  Fetched sequence: {len(sequence)} bases")
                
                # Generate multiple probe pairs for this gene
                probe_pairs = self.generate_multiple_probe_pairs(sequence, is_dna=True, num_pairs=3)
                
                # Add gene information to each probe pair
                for probe_pair in probe_pairs:
                    probe_pair['gene_name'] = gene_name
                    probe_pair['gene_number'] = i
                
                all_results.extend(probe_pairs)
                
                # Add a longer delay to be respectful to the NCBI servers and avoid rate limiting
                time.sleep(1.0)
                
            except Exception as e:
                print(f"Error processing gene '{gene_name}': {e}")
                continue
        
        return all_results
    
    def print_results(self, results: List[Dict[str, str]], output_file: Optional[str] = None):
        """Print results in a formatted way."""
        output_lines = []
        
        # Add header for the three-column format
        output_lines.append("LHS Probe (25 bases)\tRHS Probe (25 bases)\tGene Name")
        output_lines.append("-" * 50 + "\t" + "-" * 50 + "\t" + "-" * 20)
        
        for result in results:
            # Show the probe sequences and gene name in three columns
            lhs_probe = result['lhs_probe']
            rhs_probe = result['rhs_probe']
            gene_name = result.get('gene_name', 'Unknown')
            
            # Format as tab-separated columns
            output_lines.append(f"{lhs_probe}\t{rhs_probe}\t{gene_name}")
        
        # Print to console
        for line in output_lines:
            print(line.rstrip())
        
        # Write to file if specified
        if output_file:
            try:
                with open(output_file, 'w') as f:
                    for line in output_lines:
                        f.write(line + '\n')
                print(f"Results saved to {output_file}")
            except IOError as e:
                print(f"Error writing to file {output_file}: {e}")


def main():
    """Main function to handle command line interface."""
    parser = argparse.ArgumentParser(
        description="Generate complementary sequences to mRNA for given genes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate complement for a single DNA sequence
  python probe_maker.py -s "ATGCGATCGATCG"
  
  # Generate complement for a single RNA sequence
  python probe_maker.py -s "AUGCGAUCGAUCG" --rna
  
  # Process multiple genes from a file
  python probe_maker.py -f genes.txt
  
  # Process gene names and fetch sequences from NCBI
  python probe_maker.py -g top20_genes.txt
  
  # Save output to file
  python probe_maker.py -s "ATGCGATCGATCG" -o results.txt
        """
    )
    
    parser.add_argument(
        '-s', '--sequence',
        help='Single gene sequence (DNA or RNA)'
    )
    
    parser.add_argument(
        '-f', '--file',
        help='File containing gene sequences (one per line)'
    )
    
    parser.add_argument(
        '-g', '--genes',
        help='File containing gene names (one per line) to fetch from NCBI'
    )
    
    parser.add_argument(
        '--rna',
        action='store_true',
        help='Treat input sequence as RNA instead of DNA'
    )
    
    parser.add_argument(
        '-o', '--output',
        help='Output file to save results'
    )
    
    parser.add_argument(
        '--interactive',
        action='store_true',
        help='Run in interactive mode'
    )
    
    args = parser.parse_args()
    
    # Initialize ProbeMaker
    probe_maker = ProbeMaker()
    
    if args.interactive:
        run_interactive_mode(probe_maker)
    elif args.sequence:
        # Single sequence mode
        try:
            result = probe_maker.generate_mrna_complement(args.sequence, not args.rna)
            probe_maker.print_results([result], args.output)
        except ValueError as e:
            print(f"Error: {e}")
            sys.exit(1)
    
    elif args.file:
        # File mode
        try:
            with open(args.file, 'r') as f:
                gene_list = [line.strip() for line in f if line.strip()]
            
            if not gene_list:
                print("Error: File is empty or contains no valid sequences")
                sys.exit(1)
            
            results = probe_maker.process_gene_list(gene_list, not args.rna)
            if results:
                probe_maker.print_results(results, args.output)
            else:
                print("Error: No valid sequences could be processed")
                sys.exit(1)
                
        except FileNotFoundError:
            print(f"Error: File '{args.file}' not found")
            sys.exit(1)
        except IOError as e:
            print(f"Error reading file: {e}")
            sys.exit(1)
    
    elif args.genes:
        # Gene names mode - fetch sequences from NCBI
        try:
            with open(args.genes, 'r') as f:
                gene_names = [line.strip() for line in f if line.strip()]
            
            if not gene_names:
                print("Error: File is empty or contains no valid gene names")
                sys.exit(1)
            
            print(f"Fetching sequences for {len(gene_names)} genes from NCBI...")
            print("This may take a few moments...")
            
            # Initialize sequence fetcher
            sequence_fetcher = GeneSequenceFetcher()
            
            try:
                # Process gene names and generate multiple probe pairs
                results = probe_maker.process_gene_names(gene_names, sequence_fetcher)
                
                if results:
                    print(f"\nSuccessfully generated {len(results)} probe pairs for {len(gene_names)} genes")
                    probe_maker.print_results(results, args.output)
                else:
                    print("Error: No valid probe pairs could be generated")
                    sys.exit(1)
                    
            finally:
                sequence_fetcher.close()
                
        except FileNotFoundError:
            print(f"Error: File '{args.genes}' not found")
            sys.exit(1)
        except IOError as e:
            print(f"Error reading file: {e}")
            sys.exit(1)
    
    else:
        # No arguments provided, show help
        parser.print_help()
        print("\nNo input provided. Use -s for single sequence, -f for file, or --interactive for interactive mode.")


def run_interactive_mode(probe_maker: ProbeMaker):
    """Run the program in interactive mode."""
    print("=== ProbeMaker Interactive Mode ===")
    print("Enter gene sequences (one per line). Type 'quit' to exit.")
    print("Type 'help' for usage information.\n")
    
    gene_list = []
    
    while True:
        try:
            user_input = input("Enter gene sequence (or command): ").strip()
            
            if user_input.lower() == 'quit':
                break
            elif user_input.lower() == 'help':
                print("\nCommands:")
                print("  quit - Exit the program")
                print("  help - Show this help")
                print("  process - Process all entered sequences")
                print("  clear - Clear all entered sequences")
                print("  list - Show all entered sequences")
                print("\nEnter DNA or RNA sequences directly.\n")
                continue
            elif user_input.lower() == 'process':
                if gene_list:
                    print(f"\nProcessing {len(gene_list)} sequences...\n")
                    results = probe_maker.process_gene_list(gene_list, True)
                    probe_maker.print_results(results)
                    print("\n" + "="*50 + "\n")
                else:
                    print("No sequences to process. Enter some sequences first.\n")
                continue
            elif user_input.lower() == 'clear':
                gene_list.clear()
                print("All sequences cleared.\n")
                continue
            elif user_input.lower() == 'list':
                if gene_list:
                    print("\nEntered sequences:")
                    for i, seq in enumerate(gene_list, 1):
                        print(f"  {i}. {seq}")
                    print()
                else:
                    print("No sequences entered yet.\n")
                continue
            elif user_input:
                # Validate sequence
                try:
                    probe_maker.clean_sequence(user_input)
                    gene_list.append(user_input)
                    print(f"Sequence added: {user_input}")
                except ValueError as e:
                    print(f"Invalid sequence: {e}")
                print()
        
        except KeyboardInterrupt:
            print("\n\nExiting...")
            break
        except EOFError:
            print("\n\nExiting...")
            break


if __name__ == "__main__":
    main()
