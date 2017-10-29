[![Documentation Status](https://readthedocs.org/projects/epigbs/badge/?version=latest)](http://epigbs.readthedocs.io/?badge=latest)

# epiGBS

Code for working with epiGBS data.

# Install requirements

1.Tools and dependencies

```SH
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install build-essential git zlib1g zlib1g-dev lbzip2 \
bzip2 dh-autoreconf python-pip python-dev pigz ncurses-dev \
libpng-dev libfreetype6-dev pkg-config gfortran \
libopenblas-dev liblapack-dev
```

2.[PEAR: a fast and accurate Illumina Paired-End reAd mergeR](https://dx.doi.org/10.1093/bioinformatics/btt593)

```SH
wget http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.10-bin-64.tar.gz
tar -xf pear-0.9.10-bin-64.tar.gz
sudo cp pear-0.9.10-bin-64/pear-0.9.10-bin-64 /usr/local/bin
sudo ln -s /usr/local/bin/pear-0.9.10-bin-64 /usr/local/bin/pear
```

3.[vsearch](https://github.com/torognes/vsearch)

```SH
git clone https://github.com/torognes/vsearch.git
cd vsearch
./autogen.sh
./configure
make
sudo make install
cd
```

4.[Seqtk](https://github.com/lh3/seqtk.git)

```SH
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
sudo cp seqtk /usr/local/bin/
cd
```

5.[samtools](http://github.com/samtools/)

```SH
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar -xf samtools-1.3.1.tar.bz2
cd samtools-1.3.1
./configure
make
sudo make install
cd
```

6.[bcftools](http://samtools.github.io/bcftools/) 

```SH
wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2
tar -xf bcftools-1.3.1.tar.bz2
cd bcftools-1.3.1
make
sudo make install
cd
```

7.[bwa-mem](https://github.com/lh3/bwa)

```SH
git clone https://github.com/lh3/bwa.git
cd bwa
make
find ./ -maxdepth 1 -type f -perm /a+x -exec sudo cp {} /usr/local/bin \;
cd
```

8.[seqtk](https://github.com/lh3/seqtk)

```SH
git clone https://github.com/lh3/seqtk)
cd seqtk
make
sudo cp seqtk /usr/local/bin/
```
9.[sambamba](https://github.com/lomereiter/sambamba)
 follow instructions on [website](https://github.com/lomereiter/sambamba).
 
```SH
git clone https://github.com/lomereiter/sambamba
cd sambamba
make

10. Install remaining requirements with pip

```SH
git clone https://github.com/thomasvangurp/epiGBS
sudo pip install -r epiGBS/requirements.txt
```









10. [usearch (C) Copyright 2013 Robert C. Edgar, all rights reserved.](http://drive5.com/usearch)

Click [here](http://drive5.com/usearch) and follow the instructions to get the 32bit binary.

```SH
sudo cp /path/to/usearch9.2.64_i86linux32 /usr/local/bin/
sudo cp usearch9.2.64_i86linux32 /usr/local/bin/
sudo chmod +x /usr/local/bin/usearch9.2.64_i86linux32
sudo ln -s /usr/local/bin/usearch9.2.64_i86linux32 /usr/local/bin/usearch
```




