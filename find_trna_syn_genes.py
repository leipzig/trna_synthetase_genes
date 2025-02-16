import requests
import pandas as pd
import time
from typing import List, Dict, Optional

class EnsemblAPI:
    def __init__(self, species="homo_sapiens"):
        self.base_url = "https://rest.ensembl.org"
        self.species = species
        self.headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }

    def search_genes(self, query: str) -> Optional[List[Dict]]:
        """
        Search for genes using Ensembl REST API
        
        Parameters:
        query (str): Search term
        
        Returns:
        List[Dict]: List of matching genes
        """
        endpoint = f"/xrefs/symbol/{self.species}/{query}?"
        
        try:
            response = requests.get(
                f"{self.base_url}{endpoint}",
                headers=self.headers
            )
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error querying Ensembl: {e}")
            return None

    def get_gene_info(self, gene_id: str) -> Optional[Dict]:
        """
        Get detailed information for a specific gene
        
        Parameters:
        gene_id (str): Ensembl gene ID
        
        Returns:
        Dict: Gene information
        """
        endpoint = f"/lookup/id/{gene_id}?expand=1"
        
        try:
            response = requests.get(
                f"{self.base_url}{endpoint}",
                headers=self.headers
            )
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error fetching gene info: {e}")
            return None

def search_trna_synthetases() -> Optional[pd.DataFrame]:
    """
    Search for tRNA synthetase genes using Ensembl API
    
    Returns:
    pandas.DataFrame: Table of tRNA synthetase genes
    """
    ensembl = EnsemblAPI()
    
    # Search patterns for tRNA synthetases
    patterns = [
        "AARS1", "AARS2", "CARS1", "CARS2", "DARS1", "DARS2", "EARS2",
        "EPRS1", "FARSA", "FARSB", "FARS2", "GARS1", "HARS1", "HARS2",
        "IARS1", "IARS2", "KARS1", "LARS1", "LARS2", "MARS1", "MARS2",
        "NARS1", "NARS2", "PARS2", "QARS1", "RARS1", "RARS2", "SARS1",
        "SARS2", "TARS1", "TARS2", "VARS1", "VARS2", "WARS1", "WARS2",
        "YARS1", "YARS2"
    ]
    
    all_genes = []
    
    for pattern in patterns:
        print(f"Searching for {pattern}...")
        genes = ensembl.search_genes(pattern)
        
        if genes:
            for gene in genes:
                # Get detailed information for each gene
                gene_info = ensembl.get_gene_info(gene.get('id'))
                if gene_info:
                    # Only include genes on standard chromosomes
                    chromosome = gene_info.get('seq_region_name')
                    # only include if the symbol appears in the genes list
                    if pattern in gene_info.get('display_name'):
                        if chromosome in [str(i) for i in range(1,23)] + ['X', 'Y', 'MT']:
                            all_genes.append({
                                'gene_id': gene_info.get('id'),
                                'symbol': gene_info.get('display_name'),
                                'description': gene_info.get('description'),
                                'chromosome': chromosome,
                                'start': gene_info.get('start'),
                                'end': gene_info.get('end'),
                                'strand': gene_info.get('strand'),
                                'biotype': gene_info.get('biotype')
                            })
                
                # Respect API rate limits
                time.sleep(0.1)
    
    if all_genes:
        # remove duplicates
        df = pd.DataFrame(all_genes)
        df = df.drop_duplicates(subset=['gene_id'])
        return df
    
    return None

def main():
    print("Searching for tRNA synthetase genes in Ensembl...")
    results = search_trna_synthetases()
    
    if results is not None:
        print(f"\nFound {len(results)} tRNA synthetase genes")
        print("\nFirst few results:")
        print(results[['symbol', 'chromosome', 'start', 'end', 'description']].head())
        
        # Save to file
        output_file = "trna_synthetase_genes_ensembl.csv"
        results.to_csv(output_file, index=False)
        print(f"\nResults saved to {output_file}")
    else:
        print("No results found or error occurred")

if __name__ == "__main__":
    main()