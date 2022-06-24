# Deep Evolutionary Forecasting

This README documents the creation of potential variants, using files in the xtal_pompd_inputs folder. Other README.md files are available in other folders. 

This git repository includes submodules. Please clone recursively:

> git clone --recurse-submodules git@github.com:deep-evo4cast/deep-evo4cast.git

## Running pompd for sequence enumeration (xtal_pompd_inputs folder)

We assume you have installed PyRosetta (release r425).

Parameters for pompd are included in the xtal_pompd_inputs/pars directory. They match the pompd defaults except for the following ones:

* scorefunction : beta_genpot
* nbenum : ALL
* threshold : 8
* rosetta_flags : design_flags

Enter the xtal_pompd_inputs folder and run:

> make complex.enum

This command will generate several files, including 'complex.enum' which contains the result of the sequence enumeration.

For more information about pompd and how to use it, please refer to the documentation in pompd submodule folder.


## Running ToulBar2 in parallel mode to compute energies on monomers

Clone toulbar2 "cpd" branch from its repository (https://github.com/toulbar2/toulbar2) and compile it with MPI option turned on, then go to the xtal_pompd_inputs directory and run:

> make monomer.cfn.gz

This command will generate a file 'monomer.cfn.gz' representing the monomer energy matrix. 

Then, all energies can be computed on the sequences enumerated on the complex using the command:

> path_to_toulbar2/build/bin/Linux/toulbar2 monomer.cfn.gz  -dee: -O=-3 -s --cpd --jobs=monomer.bow --negative-sequences=complex.enum --diffneg

Where 'monomer.bow' is a file containing the string "monomer.cfn.gz" and 'complex.enum' is the output file from the enumeration on the complex.

This will generate a report file ('complex_negative.txt') containing all sequences, their energy on the complex form, pn the monomeric form, and the difference between the two (dG)
