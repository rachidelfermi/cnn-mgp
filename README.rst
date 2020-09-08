
    
Convolutional Neural Networks for Metagenomics Gene Prediction(CNN-MGP)
--------------------------------------------

CNN-MGP is a metagenomic ORF finding tool for the prediction of protein coding genes in short, environmental DNA sequences with unknown phylogenetic origin. CNN-MGP is based on 10  models trained on 10 pre-defined GC content ranges. The scripts extract all ORFS, one hot ecoded them afterward feeds the ORFS to the right model. 
CNN-MGP analyses can be performed via the CNN-MGP website `Coming soon <https://cnnmgp.herokuapp.com/>`_, or alternatively you can run the script from the commande line. The instructions below discuss use of CNN-MGP at the command line, following a general overview of how CNN-MGP works.


Predicting genes
-----------------------------------------------------------------------

If metagenomic fragements are submitted, CNN-MGP first extracts all complete and incomplete open reading frames (ORFs) using our integrated orf finder (casting ORFS less than 80 bp), The user customize the desirable parametre for the script and the output is eather the predicted nucliotide CDS, protein sequences or both.

Table of Contents
-------------------------------------

- `License`_
- `Citation`_
- `Support & Bug Reports`_
- `Requirements`_
- `Install Dependencies`_
- `Help Menu and Usage`_

License
--------

Use or reproduction of these materials, in whole or in part, by any commercial organization whether or not for non-commercial (including research) or commercial purposes is allowed.


Citation
--------

Al-Ajlan, A., El Allali, A. CNN-MGP: Convolutional Neural Networks for Metagenomics Gene Prediction. [`CNN-MGP <https://doi.org/10.1007/s12539-018-0313-4>`_]

Support & Bug Reports
----------------------

Please log an issue on `github issue <https://github.com/rachidelfermi/cnn-mgp/issues>`_.

You can email the CARD curators or developers directly at `rachid.elfermi@gmail.com <rachid.elfermi@gmail.com>`_.

Python version
--------------------

-Install python 3.7  (64bit) or higher from the officiel website `Python 3.7 <https://www.python.org/downloads/release/python-377/>`_.

Requirements
--------------------

- tqdm
- tensorflow==2.3.0


Install Dependencies
--------------------

- pip3 install -r requirements.txt
- or
- pip3 install tqdm tensorflow==2.3.0


Help Menu and Usage
----------------------

The following command will bring up CNN-MGP's main help menu:

.. code-block:: sh

   CNN-MGP --help

.. code-block:: sh

      usage: CNN-MGP <command> [<args>]
            commands are:
               ---------------------------------------------------------------------------------------
               -i
               ---------------------------------------------------------------------------------------

               load the input file(fasta)

               ---------------------------------------------------------------------------------------
               -o
               ---------------------------------------------------------------------------------------

              Specify the output file name 

               ---------------------------------------------------------------------------------------
               -min 
               ---------------------------------------------------------------------------------------
               
               The minimun orf lenght, Default 80
               
               ---------------------------------------------------------------------------------------
               -u
               ---------------------------------------------------------------------------------------
               
               Type 1 for unresolved start codons(ie, ATG, CTG, GTG, TTG) recommanded 
               Type 0 for start codon(ATG)

               ---------------------------------------------------------------------------------------
               -st
               ---------------------------------------------------------------------------------------
               
               Type nucl for the output file to be nucleotide CDS
               Type Prot for the output file to be protein translated gene
               Type Both for two output files 
               
               ---------------------------------------------------------------------------------------           
