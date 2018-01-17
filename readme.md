# MHC - Peptide binding predictor using RNN
## Setup the system
All code is written in Python3. To install all libraries using pip:

~~~sh 
pip install -r requirements.txt
~~~

## Collect the data
*This step is optional*  
Download the binding data file mhc\_ligand\_full.zip from
<http://www.iedb.org/database\_export\_v3.php>. Due to its large size, this
file is not included in the project by default.

Download the full MHC data (A\_prot.fasta, B\_prot.fasta and C\_prot.fasta) from <ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta>

Store these files in 'data/original/' directory. If it does not exist, create
the directory.

## Preprocess the data
*This step is optional*  
(If this step is executed, the previous step is required as well)  

The experimental data containing the peptide sequences and the MHC details have
to be parsed, cleaned and combined. This can be done by running
src/Preprocessor.py. It is currently configured to be executed as module. This
can be done by loading or installing the necessary libraries (found in
'requirements.txt') and executing the following command from the current
directory:

~~~sh
python3 src/Preprocessor.py
~~~

This process may take a while and results in the creation of the file TrainingDataAll.csv.

## Patient data
*This step is optional*
Original (anonymous) patient data can be found in 'data/prior\_df.txt'. To
parse the data run the script in src/patients/ExtractData.py. It will write the
data to a file called 'data.txt'. The models use this file to predict the
binding preferences for the patients data.

~~~sh 
python3 src/patients/ExtractData.py
~~~

## Train the model
The Wide and Shallow model (as described in the thesis for which this code was
written) is already set up in the file src/Predictor.py. It is set to train and
test only on peptides of length 9. This was added to the current script as an
illustration of how conditions w.r.t. the peptides used can be enforced. To train it and view
its performance, execute the following command from the current directory:

~~~sh 
python3 src/Predictor.py
~~~

Of course other models can be added by calling them in a similar way as was
done using the WS model.

The prediction on patient data will automatically be executed after the training is
complete.
