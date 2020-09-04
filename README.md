

Convolutional Neural Networks for Metagenomics Gene Prediction(CNN-MGP)
--------------------------------------------

CNN-MGP is a metagenomic ORF finding tool for the prediction of protein coding genes in short, environmental DNA sequences with unknown phylogenetic origin. CNN-MGP is based on 10  models trained on 10 pre-defined GC content ranges. The scripts extract all ORFS, one hot ecoded them afterward feeds the ORFS to the right model. 
CNN-MGP analyses can be performed via the CNN-MGP website [coming soon](), or alternatively you can run the script from the commande line. The instructions below discuss use of CNN-MGP at the command line, following a general overview of how CNN-MGP works.


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

Al-Ajlan, A., El Allali, A. CNN-MGP: Convolutional Neural Networks for Metagenomics Gene Prediction. [Epub ahead of print](https://doi.org/10.1007/s12539-018-0313-4)

Support & Bug Reports
----------------------

Please log an issue on [github issue]().

You can email the CARD curators or developers directly at [rachid.elfermi@gail.com](rachid.elfermi@gail.com)

Python version
--------------------

-Install python 3.7.7 (64bit) from the officiel website [Python 3.7.7](https://www.python.org/downloads/release/python-377/) 

Requirements
--------------------

- absl-py==0.10.0
- astunparse==1.6.3
- cachetools==4.1.1
- certifi==2020.6.20
- chardet==3.0.4
- gast==0.3.3
- google-auth==1.21.0
- google-auth-oauthlib==0.4.1
- google-pasta==0.2.0
- grpcio==1.31.0
- h5py==2.10.0
- idna==2.10
- importlib-metadata==1.7.0
- Keras-Preprocessing==1.1.2
- Markdown==3.2.2
- numpy==1.18.5
- oauthlib==3.1.0
- opt-einsum==3.3.0
- protobuf==3.13.0
- pyasn1==0.4.8
- pyasn1-modules==0.2.8
- requests==2.24.0
- requests-oauthlib==1.3.0
- rsa==4.6
- scipy==1.4.1
- six==1.15.0
- tensorboard==2.3.0
- tensorboard-plugin-wit==1.7.0
- tensorflow==2.3.0
- tensorflow-estimator==2.3.0
- termcolor==1.1.0
- tqdm==4.48.2
- urllib3==1.25.10
- Werkzeug==1.0.1
- wrapt==1.12.1
- zipp==3.1.0

Install Dependencies
--------------------

- pip3 install -r requirements.txt
- or
- pip3 install -r requirements.txt --user

Help Menu and Usage
----------------------

The following command will bring up CNN-MGP's main help menu:


   `CNN-MGP --help`


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

The scripts  support interactive mode as follow 
* First excute the script 
* submit each message prompt in the terminal

             
             
             
            
