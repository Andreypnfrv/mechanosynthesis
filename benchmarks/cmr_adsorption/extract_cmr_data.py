#!/usr/bin/env python3
"""
Script to extract CMR adsorption data from the web interface.
Since direct database download is not available, we'll scrape the web interface.
"""

import requests
import pandas as pd
from bs4 import BeautifulSoup
import yaml
from pathlib import Path
import time

# CMR Adsorption database URL
BASE_URL = "https://cmrdb.fysik.dtu.dk/adsorption/"

# Metals of interest for mechanosynthesis
PRIORITY_METALS = {
    'Au': 'HIGH',  # Gold - electrodes, platforms
    'W': 'HIGH',   # Tungsten - STM tips
    'Cu': 'MEDIUM', # Copper - electrodes  
    'Pt': 'MEDIUM', # Platinum - alternative tips
    'Ag': 'MEDIUM', # Silver
    'Fe': 'MEDIUM', # Iron - some mechanosynthesis applications
    'Ni': 'MEDIUM', # Nickel
}

# DFT functionals to extract
DFT_FUNCTIONALS = ['LDA', 'PBE', 'RPBE', 'BEEF-vdW']

def extract_adsorption_data():
    """Extract adsorption data from CMR web interface."""
    
    print("Extracting CMR adsorption data...")
    
    try:
        response = requests.get(BASE_URL)
        response.raise_for_status()
        
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Find the data table
        table = soup.find('table')
        if not table:
            print("No table found on the page")
            return None
            
        # Extract table headers and data
        headers = []
        header_row = table.find('tr')
        if header_row:
            for th in header_row.find_all(['th', 'td']):
                headers.append(th.get_text(strip=True))
        
        # Extract data rows
        data_rows = []
        for row in table.find_all('tr')[1:]:  # Skip header row
            row_data = []
            for cell in row.find_all(['td', 'th']):
                row_data.append(cell.get_text(strip=True))
            if row_data:
                data_rows.append(row_data)
        
        # Convert to DataFrame
        df = pd.DataFrame(data_rows, columns=headers)
        
        print(f"Extracted {len(df)} data points")
        print(f"Columns: {list(df.columns)}")
        print(f"Available metals: {df.iloc[:, 0].unique() if len(df) > 0 else 'None'}")
        
        return df
        
    except Exception as e:
        print(f"Error extracting data: {e}")
        return None

def filter_mechanosynthesis_data(df):
    """Filter data for mechanosynthesis-relevant metals."""
    
    if df is None or len(df) == 0:
        return None
        
    # Assume first column contains metal information
    metal_column = df.columns[0]
    
    # Filter for priority metals
    priority_data = df[df[metal_column].isin(PRIORITY_METALS.keys())].copy()
    
    # Add relevance column
    priority_data['mechanosynthesis_relevance'] = priority_data[metal_column].map(PRIORITY_METALS)
    
    print(f"Filtered to {len(priority_data)} mechanosynthesis-relevant data points")
    
    return priority_data

def save_extracted_data(df, output_dir):
    """Save extracted data in multiple formats."""
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    if df is None or len(df) == 0:
        print("No data to save")
        return
    
    # Save as CSV
    csv_path = output_dir / "cmr_adsorption_data.csv"
    df.to_csv(csv_path, index=False)
    print(f"Saved CSV: {csv_path}")
    
    # Save as YAML for integration with benchmark system
    yaml_path = output_dir / "cmr_adsorption_data.yaml"
    
    # Convert DataFrame to list of dictionaries
    data_list = df.to_dict('records')
    
    # Create benchmark-compatible format
    benchmark_data = {
        'dataset': 'cmr_adsorption',
        'source': 'DTU Computational Materials Repository',
        'url': BASE_URL,
        'description': 'Adsorption energies on transition metal surfaces',
        'data_points': len(df),
        'metals': list(df.iloc[:, 0].unique()) if len(df) > 0 else [],
        'mechanosynthesis_metals': list(PRIORITY_METALS.keys()),
        'data': data_list
    }
    
    with open(yaml_path, 'w') as f:
        yaml.dump(benchmark_data, f, default_flow_style=False)
    print(f"Saved YAML: {yaml_path}")
    
    # Save summary statistics
    summary_path = output_dir / "summary.txt"
    with open(summary_path, 'w') as f:
        f.write("CMR Adsorption Database Summary\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Total data points: {len(df)}\n")
        f.write(f"Columns: {list(df.columns)}\n\n")
        
        if len(df) > 0:
            f.write("Available metals:\n")
            for metal in sorted(df.iloc[:, 0].unique()):
                relevance = PRIORITY_METALS.get(metal, 'LOW')
                f.write(f"  {metal}: {relevance} relevance\n")
            
            f.write(f"\nMechanosynthesis priority metals found:\n")
            for metal, relevance in PRIORITY_METALS.items():
                if metal in df.iloc[:, 0].values:
                    count = len(df[df.iloc[:, 0] == metal])
                    f.write(f"  {metal}: {count} data points ({relevance} priority)\n")
    
    print(f"Saved summary: {summary_path}")

def main():
    """Main extraction workflow."""
    
    print("Starting CMR adsorption data extraction...")
    
    # Extract data from web interface
    raw_data = extract_adsorption_data()
    
    if raw_data is None:
        print("Failed to extract data")
        return
    
    # Filter for mechanosynthesis relevance
    mech_data = filter_mechanosynthesis_data(raw_data)
    
    # Save both raw and filtered data
    output_dir = Path(__file__).parent
    
    print("\nSaving raw data...")
    save_extracted_data(raw_data, output_dir / "raw")
    
    print("\nSaving mechanosynthesis-filtered data...")
    save_extracted_data(mech_data, output_dir / "mechanosynthesis")
    
    print("\nExtraction complete!")
    print(f"Data saved in: {output_dir}")

if __name__ == "__main__":
    main()