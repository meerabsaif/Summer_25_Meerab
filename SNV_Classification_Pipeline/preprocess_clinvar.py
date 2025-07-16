
import pandas as pd
import gzip
import os
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_info_field(info_string):
    info_dict = {}
    if info_string:
        pairs = info_string.split(';')
        for pair in pairs:
            if '=' in pair:
                key, value = pair.split('=', 1)
                info_dict[key] = value
    return info_dict

def parse_vcf_to_csv(vcf_path, output_csv_path):
    logging.info(f"Starting VCF parsing for: {vcf_path}")
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    processed_variants = []

    try:
        with open_func(vcf_path, 'rt') as f:
            header_cols = []
            for line in f:
                line = line.strip()
                if line.startswith('##'):
                    continue
                elif line.startswith('#CHROM'):
                    header_cols = line[1:].split('\t')
                    continue
                
                if not header_cols:
                    continue
                
                parts = line.split('\t')
                if len(parts) != len(header_cols):
                    continue

                variant_data = dict(zip(header_cols, parts))
                info_field = variant_data.get('INFO', '')
                info_dict = parse_info_field(info_field)
                clnsig = info_dict.get('CLNSIG', '')

                if 'Pathogenic' in clnsig or 'Likely_pathogenic' in clnsig:
                    label = 1
                elif 'Benign' in clnsig or 'Likely_benign' in clnsig:
                    label = 0
                else:
                    continue

                processed_variants.append({
                    'CHROM': variant_data.get('CHROM'),
                    'POS': variant_data.get('POS'),
                    'REF': variant_data.get('REF'),
                    'ALT': variant_data.get('ALT'),
                    'INFO': info_field,
                    'label': label
                })

        if processed_variants:
            df = pd.DataFrame(processed_variants)
            logging.info(f"Created DataFrame with {len(df)} variants and columns: {df.columns.tolist()}")
            df.to_csv(output_csv_path, index=False)
            logging.info(f"Successfully processed {len(processed_variants)} variants and saved to {output_csv_path}")
        else:
            logging.warning("No suitable variants found or processed from the VCF file.")

    except FileNotFoundError:
        logging.error(f"VCF file not found: {vcf_path}")
        raise
    except Exception as e:
        logging.error(f"An error occurred during VCF parsing: {e}")
        raise

if __name__ == "__main__":
    input_vcf_file = 'clinvar.vcf.gz'
    output_csv_file = 'data/clinvar_processed.csv'
    
    os.makedirs('data', exist_ok=True)
    parse_vcf_to_csv(input_vcf_file, output_csv_file)
    logging.info("VCF parsing script finished.")
