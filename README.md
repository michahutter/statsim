# computing the statistical similarity 
Chemical Similarity derived from Statistical Evaluation of Pharmaceutical Drugs
(M.C. Hutter, to be published)

statsim.pl
Compute the similarities between a query molecule and a list of other molecules (SMILES format)

usage: statsim.pl [-options] query_smiles.smi list_of_smiles.smi [alternative_freq_mtx]

options: -v  verbose output 
         -d  debugging output. Also enables verbose output
         -m  print atom to atom assignment matrix between query compound
            and the respective other compound 
If alternative_freq_mtx is provided, those will be used instead of atpairfrq{1-4}.mtx
Output is written to standard out in space separated format:
SMILES compound-name harmonic-mean geometric-mean arithmetic-mean divided-by-smaller-molecule

limitations:
SMILES containing lower case letters (e.g. c1ccccc1) are not allowed
isotopes (e.g. [13C]), except [2H] should be converted to their normal element
inorganic elements and counter ions (e.g. Fe, Mg, Na, K) should be removed


recsim.pl 

Compute the inter-similarities between a molecules given as a list of SMILES.
To determine suitable query molecules it performs also a 10 nearest neighbor 
approach as sum of ranks/scores, as well as DBSCAN cluster algorithm.
Count of top ranks for each molecule is given.
Outlier compounds are marked as "N"

usage: recsim.pl [-options] list_of_smiles.smi [alternative_freq_mtx]

options: -v  verbose output 
         -d  debugging output. Also enables verbose output
if alternative_freq_mtx is provided, those will be used instead of atpairfrq{1-4}.mtx
Output is written to standard out in space separated format:
(symmetic) matrix of simililarities between the compounds
10 NN compounds
DBSCAN results

Generating alternative frequency matrices for your set of compounds:
Convert SMILES to individual .hin files that contain the required atom-types

smi2hin.pl

convert SMILES from a .smi file to individual .hin files

usage: smi2hin.pl -mm infile.smi

limitations: aromaticity defined by lower case letters does not work!

distrib.pl

Compute the distributions/frequencies of atom pair interactions

usage: distrib.pl [-options] list_of_hin_files.txt

options: -<matrix_name>  override "atpairfrq" as matrix name 
          <matrix_name>.mtx  matrix files will be generated
