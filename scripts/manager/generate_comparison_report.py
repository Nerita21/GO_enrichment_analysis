#!/usr/bin/env python3
"""
Generate comprehensive HTML/PDF comparison report for GO clustering methods.
"""

import argparse
import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, Any

import pandas as pd

try:
    from jinja2 import Template
    HAS_JINJA2 = True
except ImportError:
    HAS_JINJA2 = False


HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GO Enrichment Analysis: Clustering Methods Comparison</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 20px;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.1);
            overflow: hidden;
        }
        
        header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px 20px;
            text-align: center;
        }
        
        header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
        }
        
        header p {
            font-size: 1.1em;
            opacity: 0.9;
        }
        
        main {
            padding: 40px;
        }
        
        section {
            margin-bottom: 40px;
        }
        
        h2 {
            color: #667eea;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }
        
        h3 {
            color: #764ba2;
            margin-top: 20px;
            margin-bottom: 10px;
        }
        
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }
        
        th {
            background: #667eea;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }
        
        td {
            padding: 12px;
            border-bottom: 1px solid #ddd;
        }
        
        tr:hover {
            background-color: #f5f7fa;
        }
        
        .metric-card {
            background: #f8f9fa;
            border-left: 4px solid #667eea;
            padding: 20px;
            margin: 15px 0;
            border-radius: 5px;
        }
        
        .metric-card strong {
            color: #667eea;
        }
        
        .summary {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }
        
        .summary-item {
            background: #f0f4ff;
            border: 2px solid #667eea;
            border-radius: 8px;
            padding: 20px;
            text-align: center;
        }
        
        .summary-item h4 {
            color: #764ba2;
            margin-bottom: 10px;
        }
        
        .summary-item .value {
            font-size: 1.8em;
            font-weight: bold;
            color: #667eea;
        }
        
        footer {
            background: #f5f7fa;
            padding: 20px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
            border-top: 1px solid #ddd;
        }
        
        .warning {
            background: #fff3cd;
            border: 1px solid #ffc107;
            color: #856404;
            padding: 15px;
            border-radius: 5px;
            margin: 15px 0;
        }
        
        .success {
            background: #d4edda;
            border: 1px solid #28a745;
            color: #155724;
            padding: 15px;
            border-radius: 5px;
            margin: 15px 0;
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>GO Enrichment Analysis</h1>
            <p>Clustering Methods Comparison Report</p>
            <p>Generated: {{ timestamp }}</p>
        </header>
        
        <main>
            <!-- Executive Summary -->
            <section>
                <h2>Executive Summary</h2>
                <p>This report compares multiple approaches for clustering Gene Ontology (GO) enrichment results, 
                   including R-based (Wang, Lin) and Python-based (Lin) methods, along with custom clustering variants.</p>
                
                <div class="summary">
                    <div class="summary-item">
                        <h4>Total Methods Compared</h4>
                        <div class="value">{{ num_methods }}</div>
                    </div>
                    <div class="summary-item">
                        <h4>Total Clusters</h4>
                        <div class="value">{{ total_clusters }}</div>
                    </div>
                    <div class="summary-item">
                        <h4>Average Cluster Size</h4>
                        <div class="value">{{ avg_cluster_size }}</div>
                    </div>
                </div>
            </section>
            
            <!-- Method Statistics -->
            <section>
                <h2>Method Statistics</h2>
                <p>Summary statistics for each clustering method:</p>
                
                <table>
                    <thead>
                        <tr>
                            <th>Method</th>
                            <th>Num Clusters</th>
                            <th>Min Size</th>
                            <th>Max Size</th>
                            <th>Mean Size</th>
                            <th>Median Size</th>
                            <th>Std Dev</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for method, stats in methods.items() %}
                        <tr>
                            <td><strong>{{ method }}</strong></td>
                            <td>{{ stats.num_clusters }}</td>
                            <td>{{ "%.2f"|format(stats.min_cluster_size) if stats.min_cluster_size == stats.min_cluster_size else "N/A" }}</td>
                            <td>{{ "%.2f"|format(stats.max_cluster_size) if stats.max_cluster_size == stats.max_cluster_size else "N/A" }}</td>
                            <td>{{ "%.2f"|format(stats.mean_cluster_size) if stats.mean_cluster_size == stats.mean_cluster_size else "N/A" }}</td>
                            <td>{{ "%.2f"|format(stats.median_cluster_size) if stats.median_cluster_size == stats.median_cluster_size else "N/A" }}</td>
                            <td>{{ "%.2f"|format(stats.std_cluster_size) if stats.std_cluster_size == stats.std_cluster_size else "N/A" }}</td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
                
                <div class="metric-card">
                    <strong>Interpretation:</strong> 
                    <ul style="margin-left: 20px; margin-top: 10px;">
                        <li>Methods with <strong>fewer clusters</strong> tend toward more aggressive merging (higher similarity thresholds)</li>
                        <li>Methods with <strong>smaller average size</strong> may be more sensitive to subtle term differences</li>
                        <li><strong>High standard deviation</strong> indicates variable cluster cohesion across the term space</li>
                    </ul>
                </div>
            </section>
            
            <!-- P-value Distribution -->
            <section>
                <h2>Significance Metrics</h2>
                <p>Distribution of adjusted p-values across methods:</p>
                
                <table>
                    <thead>
                        <tr>
                            <th>Method</th>
                            <th>Mean p.adjust</th>
                            <th>Median p.adjust</th>
                            <th>Min p.adjust</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for method, stats in methods.items() %}
                        {% if stats.mean_padjust == stats.mean_padjust %}
                        <tr>
                            <td><strong>{{ method }}</strong></td>
                            <td>{{ "%.2e"|format(stats.mean_padjust) }}</td>
                            <td>{{ "%.2e"|format(stats.median_padjust) }}</td>
                            <td>{{ "%.2e"|format(stats.min_padjust) }}</td>
                        </tr>
                        {% endif %}
                        {% endfor %}
                    </tbody>
                </table>
            </section>
            
            <!-- Method Overlap -->
            <section>
                <h2>Method Overlap Analysis</h2>
                <p>Similarity between different methods based on representative GO terms:</p>
                
                <table>
                    <thead>
                        <tr>
                            <th>Method 1</th>
                            <th>Method 2</th>
                            <th>Jaccard Similarity</th>
                            <th>Overlap Count</th>
                            <th>Total Unique</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for pair, metrics in pairwise_comparisons.items() %}
                        {% if metrics.jaccard_similarity == metrics.jaccard_similarity %}
                        <tr>
                            <td><strong>{{ metrics.method_1 }}</strong></td>
                            <td><strong>{{ metrics.method_2 }}</strong></td>
                            <td>{{ "%.3f"|format(metrics.jaccard_similarity) }}</td>
                            <td>{{ metrics.overlap_count }}</td>
                            <td>{{ metrics.total_unique_ids }}</td>
                        </tr>
                        {% endif %}
                        {% endfor %}
                    </tbody>
                </table>
                
                <div class="metric-card">
                    <strong>Interpretation:</strong> 
                    <ul style="margin-left: 20px; margin-top: 10px;">
                        <li><strong>High Jaccard similarity (>0.7)</strong> indicates methods produce similar clusterings</li>
                        <li><strong>Low similarity (<0.3)</strong> suggests fundamental differences in clustering philosophy</li>
                        <li>High overlap with <strong>fewer total unique terms</strong> indicates convergence on key GO terms</li>
                    </ul>
                </div>
            </section>
            
            <!-- Recommendations -->
            <section>
                <h2>Recommendations</h2>
                
                <div class="success">
                    <h3>Best Practices for Method Selection:</h3>
                    <ul style="margin-left: 20px; margin-top: 10px;">
                        <li><strong>For fine-grained analysis:</strong> Use methods with higher cluster counts</li>
                        <li><strong>For summarization:</strong> Use methods with more aggressive merging</li>
                        <li><strong>For robustness:</strong> Combine results from multiple methods</li>
                        <li><strong>For interpretation:</strong> Validate top clusters with literature review</li>
                    </ul>
                </div>
                
                <div class="warning">
                    <h3>Considerations:</h3>
                    <ul style="margin-left: 20px; margin-top: 10px;">
                        <li>Results depend heavily on similarity threshold and clustering parameters</li>
                        <li>Background gene annotations affect IC computation (where applicable)</li>
                        <li>Built-in methods may be more computationally efficient for large datasets</li>
                    </ul>
                </div>
            </section>
        </main>
        
        <footer>
            <p>GO Enrichment Analysis - Clustering Comparison Report</p>
            <p>For questions or improvements, please contact the analysis team.</p>
        </footer>
    </div>
</body>
</html>
"""


def load_metrics(metrics_file: str) -> Dict[str, Any]:
    """Load comparison metrics from JSON file."""
    with open(metrics_file, 'r') as f:
        return json.load(f)


def generate_html_report(metrics: Dict, output_file: str):
    """Generate HTML report from metrics."""
    if not HAS_JINJA2:
        print("Warning: jinja2 not installed. Generating simple HTML report.")
        # Fallback to simple HTML generation
        with open(output_file, 'w') as f:
            f.write(f"""
            <html>
            <head><title>GO Enrichment Report</title></head>
            <body>
            <h1>GO Enrichment Analysis Report</h1>
            <p>Generated: {datetime.now().isoformat()}</p>
            <pre>{json.dumps(metrics, indent=2)}</pre>
            </body>
            </html>
            """)
        return
    
    # Calculate summary statistics
    methods = metrics.get('methods', {})
    pairwise = metrics.get('pairwise_comparisons', {})
    
    total_clusters = sum(m.get('num_clusters', 0) for m in methods.values())
    num_methods = len(methods)
    avg_cluster_size = np.mean([
        m.get('mean_cluster_size', 0) for m in methods.values()
    ]) if methods else 0
    
    template_vars = {
        'timestamp': datetime.now().isoformat(),
        'num_methods': num_methods,
        'total_clusters': total_clusters,
        'avg_cluster_size': f"{avg_cluster_size:.2f}",
        'methods': methods,
        'pairwise_comparisons': pairwise
    }
    
    template = Template(HTML_TEMPLATE)
    html_content = template.render(**template_vars)
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"HTML report saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate comprehensive comparison report"
    )
    parser.add_argument("--metrics", required=True,
                        help="JSON metrics file from compare_clustering_methods.py")
    parser.add_argument("--output", default="final_report",
                        help="Output filename (without extension)")
    
    args = parser.parse_args()
    
    # Load metrics
    metrics = load_metrics(args.metrics)
    
    # Generate HTML
    html_file = f"{args.output}.html"
    generate_html_report(metrics, html_file)
    
    print(f"Report generation complete!")


# Fallback for numpy if not available
import numpy as np


if __name__ == "__main__":
    main()
