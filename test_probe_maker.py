#!/usr/bin/env python3
"""
Test script for ProbeMaker functionality.
Run this to verify that the program works correctly.
"""

import sys
import os

# Add current directory to path so we can import probe_maker
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from probe_maker import ProbeMaker

def test_basic_functionality():
    """Test basic DNA to RNA transcription and complement generation."""
    print("Testing basic functionality...")
    
    probe_maker = ProbeMaker()
    
    # Test DNA sequence (longer to test 50-base probe)
    test_dna = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    result = probe_maker.generate_mrna_complement(test_dna, is_dna=True)
    
    print(f"Input DNA: {result['input_sequence']}")
    print(f"mRNA: {result['mrna_sequence']}")
    print(f"Reverse Complementary: {result['reverse_complementary_sequence']}")
    print(f"50-Base Probe: {result['probe_50_base']}")
    print(f"LHS Probe (25 bases): {result['lhs_probe']}")
    print(f"RHS Probe (25 bases): {result['rhs_probe']}")
    print(f"Probe Start Position: {result['probe_start_position']}")
    print(f"GC Content: {result['gc_content_percentage']:.1f}%")
    
    # Verify transcription (T should become U)
    expected_mrna = "AUGCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCG"
    if result['mrna_sequence'] == expected_mrna:
        print("✓ DNA to RNA transcription: PASSED")
    else:
        print(f"✗ DNA to RNA transcription: FAILED. Expected {expected_mrna[:20]}..., got {result['mrna_sequence'][:20]}...")
    
    # Verify reverse complement generation (first 20 bases for readability)
    expected_reverse_complement = "CGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGCAU"
    if result['reverse_complementary_sequence'] == expected_reverse_complement:
        print("✓ Reverse complement generation: PASSED")
    else:
        print(f"✗ Reverse complement generation: FAILED. Expected {expected_reverse_complement[:20]}..., got {result['reverse_complementary_sequence'][:20]}...")
    
    # Verify 50-base probe (should be exactly 50 bases, converted to DNA)
    if len(result['probe_50_base']) == 50:
        print("✓ 50-base probe generation: PASSED (exactly 50 bases)")
        # Check that U's are converted to T's for DNA probe
        if 'U' not in result['probe_50_base'] and 'u' not in result['probe_50_base']:
            print("✓ DNA probe conversion: PASSED (U's converted to T's)")
        else:
            print(f"✗ DNA probe conversion: FAILED. Probe still contains U: {result['probe_50_base']}")
    else:
        print(f"✗ 50-base probe generation: FAILED. Expected length 50, got {len(result['probe_50_base'])}")
    
    # Verify position 25 constraint
    if result['probe_50_base'][24] == 'T':
        print("✓ Position 25 constraint: PASSED (T at position 25)")
    else:
        print(f"✗ Position 25 constraint: FAILED. Expected T at position 25, got {result['probe_50_base'][24]}")
    
    # Verify GC content constraint
    gc_percentage = result['gc_content_percentage']
    if 44 <= gc_percentage <= 72:
        print(f"✓ GC content constraint: PASSED ({gc_percentage:.1f}% within 44%-72% range)")
    else:
        print(f"✗ GC content constraint: FAILED. GC content {gc_percentage:.1f}% outside 44%-72% range")
    
    # Verify homopolymer repeat constraint
    probe_sequence = result['probe_50_base']
    has_homopolymer = False
    for i in range(len(probe_sequence) - 3):
        window = probe_sequence[i:i+4]
        if len(set(window)) == 1:  # All nucleotides in window are the same
            has_homopolymer = True
            break
    
    if not has_homopolymer:
        print("✓ Homopolymer repeat constraint: PASSED (no repeats >4 in a row)")
    else:
        print(f"✗ Homopolymer repeat constraint: FAILED. Found homopolymer repeat: {window}")
    
    # Verify LHS and RHS probe splitting
    lhs_probe = result['lhs_probe']
    rhs_probe = result['rhs_probe']
    
    if len(lhs_probe) == 25 and len(rhs_probe) == 25:
        print("✓ LHS/RHS probe splitting: PASSED (both exactly 25 bases)")
        # Verify that LHS + RHS = full probe
        if lhs_probe + rhs_probe == result['probe_50_base']:
            print("✓ LHS/RHS probe concatenation: PASSED (LHS + RHS = full probe)")
        else:
            print("✗ LHS/RHS probe concatenation: FAILED (LHS + RHS ≠ full probe)")
    else:
        print(f"✗ LHS/RHS probe splitting: FAILED. LHS: {len(lhs_probe)} bases, RHS: {len(rhs_probe)} bases")
    
    print()

def test_rna_input():
    """Test RNA input handling."""
    print("Testing RNA input...")
    
    probe_maker = ProbeMaker()
    
    # Test RNA sequence (longer to test 50-base probe)
    test_rna = "AUGCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCG"
    result = probe_maker.generate_mrna_complement(test_rna, is_dna=False)
    
    print(f"Input RNA: {result['input_sequence']}")
    print(f"mRNA: {result['mrna_sequence']}")
    print(f"Reverse Complementary: {result['reverse_complementary_sequence']}")
    print(f"50-Base Probe: {result['probe_50_base']}")
    print(f"LHS Probe (25 bases): {result['lhs_probe']}")
    print(f"RHS Probe (25 bases): {result['rhs_probe']}")
    print(f"Probe Start Position: {result['probe_start_position']}")
    print(f"GC Content: {result['gc_content_percentage']:.1f}%")
    
    # Verify RNA is handled correctly
    if result['mrna_sequence'] == test_rna:
        print("✓ RNA input handling: PASSED")
    else:
        print(f"✗ RNA input handling: FAILED. Expected {test_rna}, got {result['mrna_sequence']}")
    
    # Verify reverse complement generation for RNA
    expected_reverse_complement = "CGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGAUCGCAU"
    if result['reverse_complementary_sequence'] == expected_reverse_complement:
        print("✓ RNA reverse complement generation: PASSED")
    else:
        print(f"✗ RNA reverse complement generation: FAILED. Expected {expected_reverse_complement[:20]}..., got {result['reverse_complementary_sequence'][:20]}...")
    
    # Verify 50-base probe for RNA input
    if len(result['probe_50_base']) == 50:
        print("✓ RNA 50-base probe generation: PASSED (exactly 50 bases)")
        # Check that U's are converted to T's for DNA probe
        if 'U' not in result['probe_50_base'] and 'u' not in result['probe_50_base']:
            print("✓ RNA DNA probe conversion: PASSED (U's converted to T's)")
        else:
            print(f"✗ RNA DNA probe conversion: FAILED. Probe still contains U: {result['probe_50_base']}")
    else:
        print(f"✗ RNA 50-base probe generation: FAILED. Expected length 50, got {len(result['probe_50_base'])}")
    
    # Verify position 25 constraint for RNA
    if result['probe_50_base'][24] == 'T':
        print("✓ RNA position 25 constraint: PASSED (T at position 25)")
    else:
        print(f"✗ RNA position 25 constraint: FAILED. Expected T at position 25, got {result['probe_50_base'][24]}")
    
    # Verify GC content constraint for RNA
    gc_percentage = result['gc_content_percentage']
    if 44 <= gc_percentage <= 72:
        print(f"✓ RNA GC content constraint: PASSED ({gc_percentage:.1f}% within 44%-72% range)")
    else:
        print(f"✗ RNA GC content constraint: FAILED. GC content {gc_percentage:.1f}% outside 44%-72% range")
    
    # Verify homopolymer repeat constraint for RNA
    probe_sequence = result['probe_50_base']
    has_homopolymer = False
    for i in range(len(probe_sequence) - 3):
        window = probe_sequence[i:i+4]
        if len(set(window)) == 1:  # All nucleotides in window are the same
            has_homopolymer = True
            break
    
    if not has_homopolymer:
        print("✓ RNA homopolymer repeat constraint: PASSED (no repeats >4 in a row)")
    else:
        print(f"✗ RNA homopolymer repeat constraint: FAILED. Found homopolymer repeat: {window}")
    
    # Verify LHS and RHS probe splitting for RNA
    lhs_probe = result['lhs_probe']
    rhs_probe = result['rhs_probe']
    
    if len(lhs_probe) == 25 and len(rhs_probe) == 25:
        print("✓ RNA LHS/RHS probe splitting: PASSED (both exactly 25 bases)")
        # Verify that LHS + RHS = full probe
        if lhs_probe + rhs_probe == result['probe_50_base']:
            print("✓ RNA LHS/RHS probe concatenation: PASSED (LHS + RHS = full probe)")
        else:
            print("✗ RNA LHS/RHS probe concatenation: FAILED (LHS + RHS ≠ full probe)")
    else:
        print(f"✗ RNA LHS/RHS probe splitting: FAILED. LHS: {len(lhs_probe)} bases, RHS: {len(rhs_probe)} bases")
    
    print()

def test_sequence_cleaning():
    """Test sequence cleaning and validation."""
    print("Testing sequence cleaning...")
    
    probe_maker = ProbeMaker()
    
    # Test sequences with whitespace and mixed case
    test_cases = [
        ("  atgcgatcgatcg  ", "ATGCGATCGATCG"),
        ("AtGcGaTcGaTcG", "ATGCGATCGATCG"),
        ("ATG\nCGA\tTCG", "ATGCGATCG")
    ]
    
    for input_seq, expected in test_cases:
        try:
            cleaned = probe_maker.clean_sequence(input_seq)
            if cleaned == expected:
                print(f"✓ Sequence cleaning: PASSED for '{input_seq}'")
            else:
                print(f"✗ Sequence cleaning: FAILED for '{input_seq}'. Expected {expected}, got {cleaned}")
        except ValueError as e:
            print(f"✗ Sequence cleaning: FAILED for '{input_seq}' with error: {e}")
    
    print()

def test_batch_processing():
    """Test batch processing functionality."""
    print("Testing batch processing...")
    
    probe_maker = ProbeMaker()
    
    # Test multiple sequences (longer to test 50-base probes)
    test_sequences = [
        "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG",
        "TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA"
    ]
    
    results = probe_maker.process_gene_list(test_sequences, is_dna=True)
    
    if len(results) == len(test_sequences):
        print(f"✓ Batch processing: PASSED. Processed {len(results)} sequences.")
        for i, result in enumerate(results):
            print(f"  {i+1}. Input: {len(result['input_sequence'])} bases -> Probe: {len(result['probe_50_base'])} bases")
            # Verify each probe is exactly 50 bases and contains no U's
            if len(result['probe_50_base']) == 50 and 'U' not in result['probe_50_base'] and 'u' not in result['probe_50_base']:
                print(f"     ✓ Probe {i+1}: Valid 50-base DNA probe")
            else:
                print(f"     ✗ Probe {i+1}: Invalid probe")
    else:
        print(f"✗ Batch processing: FAILED. Expected {len(test_sequences)} results, got {len(results)}")
    
    print()

def test_probe_constraints():
    """Test probe constraint validation."""
    print("Testing probe constraints...")
    
    probe_maker = ProbeMaker()
    
    # Test valid probe (with T at position 25)
    valid_probe = "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    # Ensure position 25 (index 24) is T
    valid_probe = valid_probe[:24] + "T" + valid_probe[25:]
    is_valid, message = probe_maker.validate_probe_constraints(valid_probe, "test_mrna")
    if is_valid:
        print("✓ Valid probe validation: PASSED")
    else:
        print(f"✗ Valid probe validation: FAILED: {message}")
    
    # Test position 25 constraint (should be T)
    invalid_pos25_probe = "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    # Change position 25 (index 24) to something other than T
    invalid_pos25_probe = invalid_pos25_probe[:24] + "A" + invalid_pos25_probe[25:]
    is_valid, message = probe_maker.validate_probe_constraints(invalid_pos25_probe, "test_mrna")
    if not is_valid and "Position 25 must be T" in message:
        print("✓ Position 25 constraint validation: PASSED")
    else:
        print(f"✗ Position 25 constraint validation: FAILED: {message}")
    
    # Test GC content constraint (too low)
    low_gc_probe = "TATATATATATATATATATATATATATATATATATATATATATATATATA"
    is_valid, message = probe_maker.validate_probe_constraints(low_gc_probe, "test_mrna")
    if not is_valid and "GC content" in message and "outside allowed range" in message:
        print("✓ Low GC content validation: PASSED")
    else:
        print(f"✗ Low GC content validation: FAILED: {message}")
    
    # Test GC content constraint (too high)
    high_gc_probe = "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG"
    # Ensure position 25 (index 24) is T
    high_gc_probe = high_gc_probe[:24] + "T" + high_gc_probe[25:]
    is_valid, message = probe_maker.validate_probe_constraints(high_gc_probe, "test_mrna")
    if not is_valid and "GC content" in message and "outside allowed range" in message:
        print("✓ High GC content validation: PASSED")
    else:
        print(f"✗ High GC content validation: FAILED: {message}")
    
    # Test homopolymer repeat constraint
    homopolymer_probe = "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    # Ensure position 25 (index 24) is T
    homopolymer_probe = homopolymer_probe[:24] + "T" + homopolymer_probe[25:]
    # Add a homopolymer repeat (5 A's in a row)
    homopolymer_probe = homopolymer_probe[:10] + "AAAAA" + homopolymer_probe[15:]
    is_valid, message = probe_maker.validate_probe_constraints(homopolymer_probe, "test_mrna")
    if not is_valid and "Homopolymer repeat found" in message:
        print("✓ Homopolymer repeat validation: PASSED")
    else:
        print(f"✗ Homopolymer repeat validation: FAILED: {message}")
    
    print()

def test_error_handling():
    """Test error handling for invalid sequences."""
    print("Testing error handling...")
    
    probe_maker = ProbeMaker()
    
    # Test invalid sequences
    invalid_sequences = [
        "ATGCGATCGATCGX",  # Invalid nucleotide X
        "ATGCGATCGATCG123",  # Numbers
        "ATGCGATCGATCG!",  # Special characters
        "",  # Empty string
    ]
    
    for invalid_seq in invalid_sequences:
        try:
            probe_maker.clean_sequence(invalid_seq)
            print(f"✗ Error handling: FAILED. Should have rejected '{invalid_seq}'")
        except ValueError:
            print(f"✓ Error handling: PASSED. Correctly rejected '{invalid_seq}'")
    
    print()

def main():
    """Run all tests."""
    print("=" * 50)
    print("ProbeMaker Test Suite")
    print("=" * 50)
    print()
    
    try:
        test_basic_functionality()
        test_rna_input()
        test_sequence_cleaning()
        test_batch_processing()
        test_probe_constraints()
        test_error_handling()
        
        print("=" * 50)
        print("All tests completed!")
        print("=" * 50)
        
    except Exception as e:
        print(f"Test suite failed with error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
