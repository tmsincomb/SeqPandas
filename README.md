# BioPandas
Import genomic data to get a custom Pandas &amp; Biopython hybrid class object with fancy shortcuts to make Machine Learning preprocessing easy!

# Install
### Python3.7+
```python3
> git clone git@github.com:tmsincomb/BioPandas.git
> pip3 install pandas
> pip3 install biopython
```

# Usage
```python3
from BioPandas import pandas as pd

# Direct File Path
df = pd.read_seq('file.fasta', format='fasta')

# Already Opened BioPython Handle
from Bio import SeqIO
seqrecords = SeqIO.parse('file.fasta', format='fasta')
df = pd.DataFrame.from_seqrecords(seqrecords)
```

# Tutorial
For a complete walkthrough and to use it for a machine learning pipeline please follow the [tutorial notebook](https://github.com/tmsincomb/BioPandas/blob/master/tutorial.ipynb)
