# Start by moving into your paired folder
# cd ~/paired
# Run the prepDE.py script from the command line
# The -l flag allows us to specify that the average read length in our alignment was 127 bp (average is 126.5 but rounded up)
# The -i flag tells the script where to look for the input files. Our stringtie outputs for each sample were written to a sample-specific folder in ~/stringtie_results. 
# The -g and -t flags allow us to specify where we want the gene and transcript result files to be output to, respectively. 
python prepDE.py -l 127 -i ~/stringtie_results/ -g ~/gene_count_matrix.csv -t ~/transcript_count_matrix.csv
