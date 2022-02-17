=========
SeqPandas 
=========
Import genomic data to get a custom Pandas & Biopython hybrid class object with fancy shortcuts to make Machine Learning preprocessing easy!

* Free software: MIT license
* Documentation: https://seqpandas.readthedocs.io.


Installation
------------

.. code:: bash
    
    pip install seqpandas


Usage
-----

.. code:: python

    import seqpandas as spd

    # Direct File Path
    df = spd.read_seq('file.fasta', format='fasta')
    df = spd.read_seq('file.sam', format='sam')
    df = spd.read_vcf('file.vcf', format='vcf')
    df = spd.read_bed('file.bed', format='bed')

    # Just need BioPython Seqs? No problem!
    seqrecords = spd.read('file.fasta', format='fasta')

    # Already Opened BioPython Handle
    from Bio import SeqIO
    seqrecords = SeqIO.parse('file.fasta', format='fasta')
    df = spd.BioDataFrame.from_seqrecords(seqrecords)


Tutorial
--------
For a complete walkthrough and to use it for a machine learning pipeline please follow the `tutorial notebook <https://github.com/tmsincomb/SeqPandas/blob/master/tutorial.ipynb>`_.


Credits
-------

This package was created with Cookiecutter and the audreyr/cookiecutter-pypackage project template.
