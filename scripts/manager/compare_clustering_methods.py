#!/usr/bin/env python3
"""
Comparison and Evaluation Module for GO Term Clustering Methods

This script compares different clustering approaches by evaluating:
- Runtime and efficiency metrics
- Cluster statistics (size, density, p-value distribution)
- Overlap between methods
- Label quality and biological interpretability
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple
import logging

import pandas as pd
import numpy as np

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ClusteringComparator:
    """Compare and evaluate clustering methods."""
    
    def __init__(self):
        self.methods = {}
        self.metrics = {}
    
    def load_cluster_file(self, filepath: str, method_name: str) -> pd.DataFrame:
        """Load clustering results from TSV file."""
        try:
            df = pd.read_csv(filepath, sep="\t")
            self.methods[method_name] = df
            logger.info(f"Loaded {method_name}: {len(df)} clusters")
            return df
        except Exception as e:
            logger.error(f"Failed to load {filepath}: {e}")
            return None
    
    def compute_cluster_statistics(self, df: pd.DataFrame, method_name: str) -> Dict:
        """Compute statistics for a clustering method."""
        stats = {
            'method': method_name,
            'num_clusters': len(df),
            'min_cluster_size': df.get('cluster_size', [0]).min() if 'cluster_size' in df else np.nan,
            'max_cluster_size': df.get('cluster_size', [0]).max() if 'cluster_size' in df else np.nan,
            'mean_cluster_size': df.get('cluster_size', [0]).mean() if 'cluster_size' in df else np.nan,
            'median_cluster_size': df.get('cluster_size', [0]).median() if 'cluster_size' in df else np.nan,
            'std_cluster_size': df.get('cluster_size', [0]).std() if 'cluster_size' in df else np.nan,
        }
        
        if 'p.adjust' in df.columns:
            stats['mean_padjust'] = df['p.adjust'].mean()
            stats['median_padjust'] = df['p.adjust'].median()
            stats['min_padjust'] = df['p.adjust'].min()
        elif 'p_adjust' in df.columns:
            stats['mean_padjust'] = df['p_adjust'].mean()
            stats['median_padjust'] = df['p_adjust'].median()
            stats['min_padjust'] = df['p_adjust'].min()
        
        # Cluster density score (lower is more cohesive)
        if 'cluster_size' in df.columns:
            stats['mean_size_normalized'] = df['cluster_size'].mean() / df['cluster_size'].max()
        
        return stats
    
    def compute_method_overlap(self, method1: str, method2: str) -> Dict:
        """Compute overlap metrics between two methods."""
        df1 = self.methods.get(method1)
        df2 = self.methods.get(method2)
        
        if df1 is None or df2 is None:
            return {}
        
        # Extract GO IDs from representative_id or rep_go_id columns
        id_col1 = 'representative_id' if 'representative_id' in df1.columns else 'rep_go_id'
        id_col2 = 'representative_id' if 'representative_id' in df2.columns else 'rep_go_id'
        
        ids1 = set(df1[id_col1].unique()) if id_col1 in df1.columns else set()
        ids2 = set(df2[id_col2].unique()) if id_col2 in df2.columns else set()
        
        if not ids1 or not ids2:
            return {'jaccard': np.nan, 'overlap_count': 0}
        
        overlap = ids1.intersection(ids2)
        union = ids1.union(ids2)
        
        return {
            'jaccard_similarity': len(overlap) / len(union) if union else 0,
            'overlap_count': len(overlap),
            'total_unique_ids': len(union)
        }
    
    def generate_efficiency_report(self, output_file: str) -> pd.DataFrame:
        """Generate efficiency comparison report."""
        report_data = []
        
        for method_name, df in self.methods.items():
            stats = self.compute_cluster_statistics(df, method_name)
            report_data.append(stats)
        
        report_df = pd.DataFrame(report_data)
        report_df.to_csv(output_file, sep="\t", index=False)
        logger.info(f"Efficiency report saved to {output_file}")
        
        return report_df
    
    def generate_pairwise_comparison(self, output_file: str) -> pd.DataFrame:
        """Generate pairwise overlap comparison."""
        method_names = list(self.methods.keys())
        comparison_data = []
        
        for i, m1 in enumerate(method_names):
            for m2 in method_names[i+1:]:
                overlap = self.compute_method_overlap(m1, m2)
                overlap['method_1'] = m1
                overlap['method_2'] = m2
                comparison_data.append(overlap)
        
        if comparison_data:
            comparison_df = pd.DataFrame(comparison_data)
            comparison_df.to_csv(output_file, sep="\t", index=False)
            logger.info(f"Pairwise comparison saved to {output_file}")
            return comparison_df
        
        return pd.DataFrame()
    
    def generate_json_report(self, output_file: str) -> Dict:
        """Generate comprehensive JSON report."""
        report = {
            'methods': {},
            'pairwise_comparisons': {}
        }
        
        # Method-specific statistics
        for method_name in self.methods.keys():
            stats = self.compute_cluster_statistics(self.methods[method_name], method_name)
            report['methods'][method_name] = stats
        
        # Pairwise comparisons
        method_names = list(self.methods.keys())
        for i, m1 in enumerate(method_names):
            for m2 in method_names[i+1:]:
                key = f"{m1}_vs_{m2}"
                report['pairwise_comparisons'][key] = self.compute_method_overlap(m1, m2)
        
        with open(output_file, 'w') as f:
            json.dump(report, f, indent=2, default=str)
        
        logger.info(f"JSON report saved to {output_file}")
        return report


def main():
    parser = argparse.ArgumentParser(
        description="Compare GO term clustering methods"
    )
    parser.add_argument("--results", nargs="+", required=True,
                        help="Cluster result TSV files")
    parser.add_argument("--method-names", nargs="+",
                        help="Method names (default: auto-detect from filenames)")
    parser.add_argument("--output", default="comparison_metrics.json",
                        help="Output JSON metrics file")
    parser.add_argument("--report", default="efficiency_report.tsv",
                        help="Output efficiency report file")
    parser.add_argument("--pairwise", default="pairwise_comparison.tsv",
                        help="Output pairwise comparison file")
    
    args = parser.parse_args()
    
    comparator = ClusteringComparator()
    
    # Load cluster files
    if args.method_names:
        if len(args.method_names) != len(args.results):
            logger.error("Number of method names must match number of result files")
            sys.exit(1)
        method_name_pairs = zip(args.results, args.method_names)
    else:
        # Auto-detect from filenames
        method_name_pairs = [
            (f, Path(f).stem.rsplit('_', 1)[0]) for f in args.results
        ]
    
    for filepath, method_name in method_name_pairs:
        comparator.load_cluster_file(filepath, method_name)
    
    # Generate reports
    comparator.generate_efficiency_report(args.report)
    comparator.generate_pairwise_comparison(args.pairwise)
    comparator.generate_json_report(args.output)
    
    logger.info("Comparison complete!")


if __name__ == "__main__":
    main()
