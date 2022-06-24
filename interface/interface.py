from pyrosetta import *
from pyrosetta.rosetta.protocols.rosetta_scripts import *

import argparse


parser = argparse.ArgumentParser(
    description='Run rosetta protocol')
parser.add_argument('-p',  '--protocol',
                    help='Protocol in XML format')
parser.add_argument('-i',  '--input',
                    help='Input pdb')
parser.add_argument('-f',  '--flags',
                    help='flags file')

parser.add_argument('-o', '--output', help='Output pdb file')

args = parser.parse_args()
inpdb = args.input
flags = args.flags
protocol_file = args.protocol
outfile = args.output

if flags:
    init('-beta @'+flags)
else:
    init('-beta')
    
pose = pose_from_pdb(inpdb)
parser = RosettaScriptsParser()
protocol = parser.generate_mover(protocol_file)
protocol.apply(pose)
pose.dump_pdb(outfile)
