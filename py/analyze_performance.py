#!/usr/bin/env python3
"""
Performance analysis script for mechanosynthesis calculations.

Provides utilities to analyze and visualize calculation performance metrics
from the CSV tracking database.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import sys
from datetime import datetime, timedelta
import argparse

sys.path.append(str(Path(__file__).parent))
import settings
from kernel.utils.performance_tracker import PerformanceTracker


def load_performance_data(csv_path=None):
    """Load performance data from CSV into pandas DataFrame."""
    if csv_path is None:
        csv_path = settings.PROJECT_ROOT / "performance_logs.csv"
    
    if not csv_path.exists():
        print(f"No performance data found at {csv_path}")
        return pd.DataFrame()
    
    try:
        df = pd.read_csv(csv_path)
        
        # Convert timestamp to datetime
        df['timestamp'] = pd.to_datetime(df['timestamp'])
        
        # Convert numeric columns
        numeric_cols = ['n_atoms', 'n_electrons', 'n_cores', 'wall_time_sec', 
                       'memory_peak_mb', 'convergence_steps', 'final_energy', 'quality_rating']
        
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
                
        # Convert boolean column
        if 'success' in df.columns:
            df['success'] = df['success'].astype(bool)
            
        return df
    except Exception as e:
        print(f"Error loading performance data: {e}")
        return pd.DataFrame()


def print_summary_stats(df):
    """Print summary statistics."""
    if df.empty:
        print("No data to analyze")
        return
        
    print("=== PERFORMANCE SUMMARY ===\n")
    
    # Basic stats
    print(f"Total calculations: {len(df)}")
    print(f"Date range: {df['timestamp'].min().strftime('%Y-%m-%d')} to {df['timestamp'].max().strftime('%Y-%m-%d')}")
    print(f"Success rate: {df['success'].mean():.1%}")
    print()
    
    # Method breakdown
    if 'method' in df.columns:
        print("Method usage:")
        method_counts = df['method'].value_counts()
        for method, count in method_counts.items():
            success_rate = df[df['method'] == method]['success'].mean()
            print(f"  {method}: {count} calculations ({success_rate:.1%} success)")
        print()
    
    # Performance stats
    if 'wall_time_sec' in df.columns and not df['wall_time_sec'].isna().all():
        successful = df[df['success'] == True]
        if not successful.empty:
            print("Timing statistics (successful calculations only):")
            print(f"  Average time: {successful['wall_time_sec'].mean():.1f}s")
            print(f"  Median time: {successful['wall_time_sec'].median():.1f}s")
            print(f"  Fastest: {successful['wall_time_sec'].min():.1f}s")
            print(f"  Slowest: {successful['wall_time_sec'].max():.1f}s")
            print()
    
    # Molecule size stats
    if 'n_atoms' in df.columns and not df['n_atoms'].isna().all():
        print("Molecule size statistics:")
        print(f"  Average atoms: {df['n_atoms'].mean():.1f}")
        print(f"  Size range: {df['n_atoms'].min():.0f} - {df['n_atoms'].max():.0f} atoms")
        print()
        
    # Quality ratings
    if 'quality_rating' in df.columns and not df['quality_rating'].isna().all():
        print("Quality ratings:")
        quality_counts = df['quality_rating'].value_counts().sort_index()
        for rating, count in quality_counts.items():
            print(f"  Rating {rating}: {count} calculations")
        avg_quality = df['quality_rating'].mean()
        print(f"  Average quality: {avg_quality:.1f}/5")
        print()


def plot_method_comparison(df, save_path=None):
    """Plot method comparison charts."""
    if df.empty or 'method' not in df.columns:
        return
        
    successful = df[df['success'] == True].copy()
    if successful.empty:
        print("No successful calculations to plot")
        return
        
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Method Performance Comparison', fontsize=16, fontweight='bold')
    
    # 1. Success rate by method
    if len(df) > 0:
        success_by_method = df.groupby('method')['success'].agg(['count', 'sum']).reset_index()
        success_by_method['success_rate'] = success_by_method['sum'] / success_by_method['count']
        
        axes[0,0].bar(success_by_method['method'], success_by_method['success_rate'])
        axes[0,0].set_title('Success Rate by Method')
        axes[0,0].set_ylabel('Success Rate')
        axes[0,0].set_ylim(0, 1)
        for i, v in enumerate(success_by_method['success_rate']):
            axes[0,0].text(i, v + 0.02, f'{v:.1%}', ha='center')
    
    # 2. Time vs atoms scatter
    if 'wall_time_sec' in successful.columns and 'n_atoms' in successful.columns:
        methods = successful['method'].unique()
        colors = plt.cm.Set1(np.linspace(0, 1, len(methods)))
        
        for method, color in zip(methods, colors):
            method_data = successful[successful['method'] == method]
            if not method_data.empty:
                axes[0,1].scatter(method_data['n_atoms'], method_data['wall_time_sec'], 
                                label=method, alpha=0.7, color=color)
        
        axes[0,1].set_xlabel('Number of Atoms')
        axes[0,1].set_ylabel('Wall Time (seconds)')
        axes[0,1].set_title('Performance Scaling')
        axes[0,1].legend()
        axes[0,1].set_yscale('log')
    
    # 3. Average time by method
    if 'wall_time_sec' in successful.columns:
        time_by_method = successful.groupby('method')['wall_time_sec'].agg(['mean', 'std']).reset_index()
        
        axes[1,0].bar(time_by_method['method'], time_by_method['mean'], 
                     yerr=time_by_method['std'], capsize=5)
        axes[1,0].set_title('Average Calculation Time')
        axes[1,0].set_ylabel('Time (seconds)')
        axes[1,0].set_yscale('log')
    
    # 4. Quality distribution
    if 'quality_rating' in successful.columns:
        quality_data = []
        methods = []
        for method in successful['method'].unique():
            method_data = successful[successful['method'] == method]['quality_rating'].dropna()
            if len(method_data) > 0:
                quality_data.append(method_data.tolist())
                methods.append(method)
        
        if quality_data:
            axes[1,1].boxplot(quality_data, labels=methods)
            axes[1,1].set_title('Quality Rating Distribution')
            axes[1,1].set_ylabel('Quality Rating (1-5)')
            axes[1,1].set_ylim(0.5, 5.5)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Method comparison plot saved to {save_path}")
    else:
        plt.show()


def plot_timeline(df, save_path=None):
    """Plot calculation timeline."""
    if df.empty:
        return
        
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    fig.suptitle('Calculation Timeline', fontsize=16, fontweight='bold')
    
    # 1. Calculations per day
    df_daily = df.set_index('timestamp').resample('D').size()
    axes[0].plot(df_daily.index, df_daily.values, marker='o')
    axes[0].set_title('Calculations per Day')
    axes[0].set_ylabel('Number of Calculations')
    axes[0].grid(True, alpha=0.3)
    
    # 2. Cumulative time spent
    successful = df[df['success'] == True].copy()
    if 'wall_time_sec' in successful.columns and not successful['wall_time_sec'].isna().all():
        successful = successful.sort_values('timestamp')
        successful['cumulative_time'] = successful['wall_time_sec'].cumsum() / 3600  # Convert to hours
        
        axes[1].plot(successful['timestamp'], successful['cumulative_time'], marker='.')
        axes[1].set_title('Cumulative Computation Time')
        axes[1].set_ylabel('Total Hours')
        axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Timeline plot saved to {save_path}")
    else:
        plt.show()


def generate_report(output_dir=None):
    """Generate comprehensive performance report."""
    if output_dir is None:
        output_dir = settings.PROJECT_ROOT / "performance_reports"
    else:
        output_dir = Path(output_dir)
        
    output_dir.mkdir(exist_ok=True)
    
    df = load_performance_data()
    
    if df.empty:
        print("No performance data available for report generation")
        return
    
    # Generate timestamp for report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    
    # Summary statistics
    print_summary_stats(df)
    
    # Save plots
    method_plot_path = output_dir / f"method_comparison_{timestamp}.png"
    timeline_plot_path = output_dir / f"timeline_{timestamp}.png"
    
    plot_method_comparison(df, method_plot_path)
    plot_timeline(df, timeline_plot_path)
    
    # Save detailed CSV
    detailed_csv = output_dir / f"performance_data_{timestamp}.csv"
    df.to_csv(detailed_csv, index=False)
    print(f"Detailed data saved to {detailed_csv}")
    
    print(f"\n📊 Performance report generated in {output_dir}")


def filter_data(df, **filters):
    """Filter performance data based on criteria."""
    filtered = df.copy()
    
    for key, value in filters.items():
        if key in filtered.columns:
            if isinstance(value, list):
                filtered = filtered[filtered[key].isin(value)]
            else:
                filtered = filtered[filtered[key] == value]
                
    return filtered


def main():
    """Main CLI interface."""
    parser = argparse.ArgumentParser(description='Analyze mechanosynthesis calculation performance')
    parser.add_argument('--summary', action='store_true', help='Print summary statistics')
    parser.add_argument('--plot-methods', action='store_true', help='Plot method comparison')
    parser.add_argument('--plot-timeline', action='store_true', help='Plot timeline')
    parser.add_argument('--report', action='store_true', help='Generate full report')
    parser.add_argument('--method', help='Filter by calculation method')
    parser.add_argument('--molecule', help='Filter by molecule name')
    parser.add_argument('--successful-only', action='store_true', help='Show only successful calculations')
    parser.add_argument('--output-dir', help='Output directory for reports and plots')
    
    args = parser.parse_args()
    
    # Load data
    df = load_performance_data()
    
    if df.empty:
        print("No performance data found. Run some calculations with performance tracking first.")
        return
    
    # Apply filters
    filters = {}
    if args.method:
        filters['method'] = args.method
    if args.molecule:
        filters['molecule_name'] = args.molecule
    if args.successful_only:
        filters['success'] = True
        
    if filters:
        df = filter_data(df, **filters)
        print(f"Applied filters: {filters}")
        print(f"Filtered data: {len(df)} calculations\n")
    
    # Execute requested actions
    if args.summary or not any([args.plot_methods, args.plot_timeline, args.report]):
        print_summary_stats(df)
    
    if args.plot_methods:
        plot_method_comparison(df)
    
    if args.plot_timeline:
        plot_timeline(df)
        
    if args.report:
        generate_report(args.output_dir)


if __name__ == "__main__":
    main()