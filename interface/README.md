Assuming PyRosetta is installed on your system, the interface between RBD and ACE2 can be computed with the following command:

		python3 ./interface.py -i ../xtal_pompd_inputs/complex.pdb -p get_ace2rbd_interface.xml -o interface.pdb
		
This command will generate a new PDB file with interface residue numbers at the end
