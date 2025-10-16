#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  create_functional_regions.py — create functional region coordinates     ║
# ║                                                                          ║
# ║  Based on the old R script strategy but with configurable parameters    ║
# ║  from the user's pipeline configuration.                                ║
# ║                                                                          ║
# ║  Inputs  : genes.tsv file (from gtf_to_catalog.py)                      ║
# ║  Outputs : functional_regions.tsv with all region coordinates           ║
# ║                                                                          ║
# ║  Strategy (matching old R script):                                      ║
# ║   • Calculate functional regions for each gene based on TSS/TES         ║
# ║   • Use configurable parameters instead of hardcoded values             ║
# ║   • Output comprehensive coordinate file for all functional regions     ║
# ╚══════════════════════════════════════════════════════════════════════════╝

from __future__ import annotations
import argparse
import sys
import os
from typing import Dict, List, Tuple

def _clamp_pair(a: int, b: int) -> Tuple[int, int]:
    """Ensure start/end are ordered and non-negative (BED-style)."""
    if a > b:
        a, b = b, a
    if a < 0:
        b = b - a  # preserve length when shifting to 0
        a = 0
    return a, b
def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Create functional region coordinates from gene list",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("--genes", required=True, help="Input genes.tsv file")
    parser.add_argument("--output", required=True, help="Output functional_regions.tsv file")
    
    # Functional region parameters (with defaults matching old R script)
    parser.add_argument("--prom-up", type=int, default=250, 
                       help="Promoter upstream flank (bp)")
    parser.add_argument("--prom-down", type=int, default=250, 
                       help="Promoter downstream flank (bp)")
    parser.add_argument("--div-inner", type=int, default=251, 
                       help="Divergent transcription inner boundary (bp)")
    parser.add_argument("--div-outer", type=int, default=750, 
                       help="Divergent transcription outer boundary (bp)")
    parser.add_argument("--cps-offset", type=int, default=500, 
                       help="CPS (cleavage/polyA site) offset (bp)")
    parser.add_argument("--tw-length", type=int, default=10000, 
                       help="Termination window length (bp)")
    
    return parser.parse_args()

def read_genes(genes_file: str) -> List[Dict[str, any]]:
    """Read genes.tsv file and return list of gene dictionaries."""
    genes = []
    
    with open(genes_file, 'r') as f:
        header = f.readline().strip().split('\t')
        
        # Expected columns: gene_id, gene_name, chr, strand, start, end, tss, tes, biotype
        expected_cols = ['gene_id', 'gene_name', 'chr', 'strand', 'start', 'end', 'tss', 'tes', 'biotype']
        
        for i, col in enumerate(header):
            if col.lower() in expected_cols:
                expected_cols[expected_cols.index(col.lower())] = col
        
        for line in f:
            if not line.strip():
                continue
                
            fields = line.strip().split('\t')
            if len(fields) != len(header):
                continue
                
            gene = {}
            for i, field in enumerate(fields):
                if i < len(header):
                    gene[header[i]] = field
            
            # Convert numeric fields
            try:
                gene['start'] = int(gene['start'])
                gene['end'] = int(gene['end'])
                gene['tss'] = int(gene['tss'])
                gene['tes'] = int(gene['tes'])
                
                # Include all genes (user can control filtering via parameters)
                genes.append(gene)
                    
            except (ValueError, KeyError):
                continue
    
    return genes

def calculate_regions(gene: Dict[str, any], args: argparse.Namespace) -> Dict[str, any]:
    """Calculate functional regions for a single gene based on strand."""
    
    strand = gene['strand']
    tss = gene['tss']
    tes = gene['tes']
    
    # Initialize result with gene info
    result = gene.copy()
    
    if strand == "+":
        # Plus strand calculations
        result['TSS'] = tss
        result['CPS'] = tes
        
        # Divergent transcription (opposite strand)
        divs, dive = tss - args.div_outer, tss - args.div_inner
        result['DIVs'], result['DIVe'] = _clamp_pair(divs, dive)
        
        # Promoter (same strand)
        pps, ppe = tss - args.prom_up, tss + args.prom_down
        result['PPs'], result['PPe'] = _clamp_pair(pps, ppe)
        
        # Gene body
        gbs, gbe = tss + args.prom_down, tes - args.cps_offset - 1
        result['GBs'], result['GBe'] = _clamp_pair(gbs, gbe)
        
        # CPS (cleavage/polyA site)
        cpss, cpse = tes - args.cps_offset, tes + args.cps_offset - 1
        result['CPSs'], result['CPSe'] = _clamp_pair(cpss, cpse)
        
        # Termination window
        tws, twe = tes + args.cps_offset, tes + args.tw_length
        result['TWs'], result['TWe'] = _clamp_pair(tws, twe)
        
    else:  # strand == "-"
        # Minus strand calculations: use TSS at the gene end (higher coordinate)
        # Inputs from genes.tsv already encode: for minus, tss = gene end, tes = gene start
        result['TSS'] = max(0, tss)
        result['CPS'] = max(0, tes)

        # Divergent transcription (opposite strand) upstream of TSS in transcriptional sense
        # For minus, upstream lies toward increasing genomic coordinates (+)
        divs, dive = tss + args.div_inner, tss + args.div_outer
        result['DIVs'], result['DIVe'] = _clamp_pair(divs, dive)

        # Promoter (same strand) centered on TSS
        pps, ppe = tss - args.prom_up, tss + args.prom_down
        result['PPs'], result['PPe'] = _clamp_pair(pps, ppe)

        # Gene body: from CPS end toward TSS
        gbs, gbe = tes + args.cps_offset + 1, tss - args.prom_down
        result['GBs'], result['GBe'] = _clamp_pair(gbs, gbe)

        # CPS (cleavage/polyA site) around CPS coordinate
        cpss, cpse = tes - args.cps_offset + 1, tes + args.cps_offset
        result['CPSs'], result['CPSe'] = _clamp_pair(cpss, cpse)

        # Termination window downstream of CPS (in transcriptional sense)
        # For minus, this is toward decreasing genomic coordinates (-)
        tws, twe = tss - args.tw_length, tss - args.cps_offset
        result['TWs'], result['TWe'] = _clamp_pair(tws, twe)
    
    # Add promoter coordinates (±500 bp around TSS for compatibility)
    pc1, pc2 = result['TSS'] - 500, result['TSS'] + 500
    result['promC1'], result['promC2'] = _clamp_pair(pc1, pc2)
    
    return result

def write_functional_regions(genes_with_regions: List[Dict[str, any]], output_file: str) -> None:
    """Write functional regions to output file."""
    
    # Define column order (matching old R script output)
    columns = [
        'chr', 'start', 'end', 'strand', 'gene_name',  # Original gene info
        'TSS', 'CPS',  # Key coordinates
        'DIVs', 'DIVe',  # Divergent transcription
        'PPs', 'PPe',    # Promoter
        'GBs', 'GBe',    # Gene body
        'CPSs', 'CPSe',  # CPS region
        'TWs', 'TWe',    # Termination window
        'promC1', 'promC2'  # Promoter coordinates
    ]
    
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(columns) + '\n')
        
        # Write data
        for gene in genes_with_regions:
            row = []
            for col in columns:
                if col in gene:
                    row.append(str(gene[col]))
                else:
                    row.append('')
            f.write('\t'.join(row) + '\n')

def main() -> int:
    """Main function."""
    args = parse_args()
    
    if not os.path.exists(args.genes):
        print(f"ERROR: Input file {args.genes} does not exist", file=sys.stderr)
        return 1
    
    print(f"INFO: Reading genes from {args.genes}", file=sys.stderr)
    genes = read_genes(args.genes)
    
    if not genes:
        print("ERROR: No valid genes found in input file", file=sys.stderr)
        return 1
    
    print(f"INFO: Processing {len(genes)} genes", file=sys.stderr)
    
    # Calculate functional regions for each gene
    genes_with_regions = []
    for gene in genes:
        try:
            regions = calculate_regions(gene, args)
            genes_with_regions.append(regions)
        except Exception as e:
            print(f"WARNING: Failed to process gene {gene.get('gene_name', 'unknown')}: {e}", file=sys.stderr)
            continue
    
    print(f"INFO: Writing functional regions to {args.output}", file=sys.stderr)
    write_functional_regions(genes_with_regions, args.output)
    
    print(f"INFO: Successfully created functional regions for {len(genes_with_regions)} genes", file=sys.stderr)
    return 0

if __name__ == "__main__":
    sys.exit(main())
