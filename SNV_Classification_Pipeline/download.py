import os
import requests
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def download_clinvar(url="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz", output_file="clinvar.vcf.gz"):
    if os.path.exists(output_file):
        logging.info(f"'{output_file}' already exists. Skipping download.")
        return
    try:
        logging.info(f"Downloading {url} ...")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(output_file, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        logging.info(f"Downloaded and saved as '{output_file}'")
    except Exception as e:
        logging.error(f"Error downloading file: {e}")
        if os.path.exists(output_file):
            os.remove(output_file)
        raise

if __name__ == "__main__":
    download_clinvar()
