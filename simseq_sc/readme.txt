# install cmake in advance
conda install cmake
conda install -c conda-forge cxx-compiler

# install fast simulator
bash install.sh

# run examples
## fasta references are *.fa files which are located in ref/ folder

### 1) an example to generate sequencing with uniform distribution
./build/simseq_sc ref/ -o simUni -b 0 
./build/simseq_sc ref/ -o simUni -b 0 -rc readcount.txt

### 2) an example to generate sequencing with 3 prime bias
./build/simseq_sc ref/ -o sim3 -b 3
./build/simseq_sc ref/ -o sim3 -b 0 -rc readcount.txt


### 3) an example to generate sequencing with 5 prime bias
./build/simseq_sc ref/ -o sim5 -b 5
./build/simseq_sc ref/ -o sim5 -b 0 -rc readcount.txt
