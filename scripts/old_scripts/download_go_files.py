import sys, os
import urllib

main_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
toolbox_dir = os.path.join(main_dir, 'diana/toolbox')

# Download files
print("  DIANA INFO:\t.... Downloading go-basic.obo ....\n")
urllib.urlretrieve('http://geneontology.org/ontology/go-basic.obo', os.path.join(toolbox_dir, 'go-basic.obo'))

print("  DIANA INFO:\t.... Downloading gene2go.gz ....\n")
compressed_file = os.path.join(toolbox_dir, 'gene2go.gz')
urllib.urlretrieve('ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz', compressed_file)

# Decompress gene2go.gz
command = 'gzip -d {}'.format(compressed_file)
os.system(command)

print("  DIANA INFO:\t.... Download finished. Thank you for your patience!\n")