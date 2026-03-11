#!/usr/bin/env python3

import os
import sys
import argparse
import time
from io import StringIO
import pandas as pd
from tqdm import tqdm

API_KEY = 'AIzaSyBlunrEyxdUwcWHjoMA7dlMOE_FSee_MXg'

def parse_args():
    parser = argparse.ArgumentParser(description='AlphaGenome variant score')
    parser.add_argument('-i', '--input', type=str, default=None, help='inpute VCF')
    parser.add_argument('-o', '--output', type=str, default='variant_scores.csv', help='output')
    return parser.parse_args()

try:
    from alphagenome import colab_utils
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers
except ImportError as e:
    print(f"error AlphaGenome model: {e}")
    print("conda 'alphagenome'")
    print("conda activate alphagenome")
    sys.exit(1)

def create_dna_model(api_key):
    try:
        model = dna_client.create(api_key)
        print("✓ AlphaGenome API")
        return model
    except Exception as e:
        print(f"✗AlphaGenome API faile: {e}")
        print("please check API")
        return None

def load_vcf_data(input_file=None):
        vcf = pd.read_csv(StringIO(vcf_file), sep='\t')
        print(f"✓  {len(vcf)} variants test")
        return vcf
    else:
        try:
            vcf = pd.read_csv(input_file, sep='\t')
            print(f"✓ {len(vcf)} variants: {input_file}")
            return vcf
        except Exception as e:
            print(f"✗ faile VCF: {e}")
            return None

def validate_vcf(vcf):
    required_columns = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']
    missing_columns = [col for col in required_columns if col not in vcf.columns]
    
    if missing_columns:
        print(f"✗ VCF error: {missing_columns}")
        return False
    
    print("✓ VCF OK")
    return True

def get_scorers(organism='human'):
    scorer_selections = {
        'rna_seq': True,
        'cage': True,
        'procap': True,
        'atac': True,
        'dnase': True,
        'chip_histone': True,
        'chip_tf': True,
        'polyadenylation': True,
        'splice_sites': True,
        'splice_site_usage': True,
        'splice_junctions': True,
    }

    organism_map = {
        'human': dna_client.Organism.HOMO_SAPIENS,
        'mouse': dna_client.Organism.MUS_MUSCULUS,
    }
    organism_enum = organism_map.get(organism.lower(), dna_client.Organism.HOMO_SAPIENS)

    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    selected_scorers = [
        all_scorers[key]
        for key in all_scorers
        if scorer_selections.get(key.lower(), False)
    ]
    
    unsupported_scorers = [
        scorer
        for scorer in selected_scorers
        if (
            organism_enum.value
            not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
        )
        | (
            (scorer.requested_output == dna_client.OutputType.PROCAP)
            & (organism_enum == dna_client.Organism.MUS_MUSCULUS)
        )
    ]
    
    if unsupported_scorers:
        print(f"⚠ : {len(unsupported_scorers)}")
        for unsupported_scorer in unsupported_scorers:
            selected_scorers.remove(unsupported_scorer)
    
    return selected_scorers, organism_enum

def score_variants(dna_model, vcf, selected_scorers, organism, sequence_length='1MB'):

    sequence_length_enum = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f'SEQUENCE_LENGTH_{sequence_length}'
    ]
    
    results = []
    total_variants = len(vcf)
    successful_variants = 0
    failed_variants = 0
    
    print(f"start {total_variants} variants..")
    print(f"length: {sequence_length}")
    print(f"score: {len(selected_scorers)}")
    print("-" * 50)

    start_time = time.time()
    
    with tqdm(total=total_variants, desc="", unit="") as pbar:
        for i, vcf_row in vcf.iterrows():
            try:
                pbar.set_description(f": {vcf_row.variant_id}")
                
                variant = genome.Variant(
                    chromosome=str(vcf_row.CHROM),
                    position=int(vcf_row.POS),
                    reference_bases=vcf_row.REF,
                    alternate_bases=vcf_row.ALT,
                    name=vcf_row.variant_id,
                )
                interval = variant.reference_interval.resize(sequence_length_enum)

                variant_scores = dna_model.score_variant(
                    interval=interval,
                    variant=variant,
                    variant_scorers=selected_scorers,
                    organism=organism,
                )
                results.append(variant_scores)
                successful_variants += 1
                
                if successful_variants % 10 == 0:
                    progress = (i + 1) / total_variants * 100
                    elapsed_time = time.time() - start_time
                    avg_time_per_variant = elapsed_time / (i + 1)
                    remaining_variants = total_variants - (i + 1)
                    estimated_remaining_time = remaining_variants * avg_time_per_variant
                    
                    def format_time(seconds):
                        if seconds < 60:
                            return f"{seconds:.0f}s"
                        elif seconds < 3600:
                            minutes = seconds / 60
                            return f"{minutes:.1f}min"
                        else:
                            hours = seconds / 3600
                            return f"{hours:.1f}h"
                    
                    pbar.set_postfix({
                        'success': successful_variants,
                        'faile': failed_variants,
                        'process': f'{progress:.1f}%',
                        'remains': format_time(estimated_remaining_time)
                    })
                
            except Exception as e:
                failed_variants += 1
                print(f"\n⚠ {vcf_row.variant_id}: {e}")
                pbar.set_postfix({
                    'success': successful_variants,
                    'faile': failed_variants,
                    'process': vcf_row.variant_id
                })
                continue
            
            pbar.update(1)
    
    total_time = time.time() - start_time

    return results

def main():
    args = parse_args()
    print("=" * 50)
    print("AlphaGenome")
    print("=" * 50)
 
    dna_model = create_dna_model(API_KEY)
    if not dna_model:
        sys.exit(1)

    vcf = load_vcf_data(args.input)
    if vcf is None:
        sys.exit(1)

    if not validate_vcf(vcf):
        sys.exit(1)

    selected_scorers, organism = get_scorers('human')
    print(f"✓ {len(selected_scorers)} used")
    
    results = score_variants(dna_model, vcf, selected_scorers, organism)
    
    if not results:
        sys.exit(1)
    

    try:
        df_scores = variant_scorers.tidy_scores(results)
        print("✓ complete!")
        print("\nresults:")
        print(df_scores)
        

        df_scores.to_csv(args.output, index=False)
        print(f"\nresults sorted at: {args.output}")
        
    except Exception as e:
        print(f"✗ results failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 
