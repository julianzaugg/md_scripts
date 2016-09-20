#!/bin/bash

"""
Read in a directory of pdbs and minimise them using a provided em.mdp file.
Assume GROMACS double precision installed to path (gmx_d).
em.mdp file should be in the output directory.
Also assumes that the pdbs have been pre-processed, i.e., waters and ligands removed.

Edit this script to be specific to your system as required, i.e., grep # # #.
"""



"""
For pdb in input directory
	Create folder for system/em
	Generate topology
	Editconf
	Solvate
	Grompp using Use base em.mdp
	MDrun
	Generate pdb
"""

import sys, os
import argparse
import subprocess
from collections import defaultdict

def _create_pdb2gmx_params(structure_filename):
	with open(structure_filename, 'r') as fh:
		data = [line.strip().split() for line in fh.readlines()]
		counts = defaultdict(int)
		residues_seen = set()
		for line in data:
			# print line
			if not line or line[0] != "ATOM" or line[4] in residues_seen: continue
			counts[line[3]] += 1
			residues_seen.add(line[4])
		# print structure_filename.split("/")[-1], counts
		# sys.exit()
	param_list = [] # 14 is the forcefield number
	for aa in ["LYS", "ARG", "GLN", "ASP", "GLU", "HIS"]:
		if aa in ["LYS", "ARG"]:
			param_list += ["1" for i in xrange(counts[aa])]
		elif aa in ["GLN", "ASP", "GLU"]:
			param_list += ["0" for i in xrange(counts[aa])]
		elif aa in ["HIS"]:
			param_list += ["2" for i in xrange(counts[aa])]
	param_list += ["0", "0"]
	return " ".join(param_list)

def _process_arguments(myparser, myargs):
	base_out_dirname =  os.path.abspath(myargs.output)
	for pdb_filename in [filename for filename in os.listdir(myargs.input) if filename.endswith(".pdb")]:
		structure_name = pdb_filename.split(".")[-2]
		if structure_name == "B_1": continue
		pdb_out_dirname =  os.path.abspath(base_out_dirname + "/" + structure_name)
		input_pdb_full_filename = os.path.abspath("/".join([myargs.input, pdb_filename]))
		if not os.path.exists(pdb_out_dirname): os.mkdir(pdb_out_dirname)
		os.chdir(pdb_out_dirname)
		pdb2gmx_params = _create_pdb2gmx_params(input_pdb_full_filename)
		# print "echo %s | gmx_d pdb2gmx -f %s -o %s.pdb -p %s.top -i %s.itp -water spc -inter -ignh" % (pdb2gmx_params, input_pdb_full_filename,
		# 																				structure_name,
		# 																				structure_name,
		# 																				structure_name)
		subprocess.call("echo %s | gmx_d pdb2gmx -f %s -o %s.pdb -p %s.top -i %s.itp -water spc -inter -ignh -ff gromos54a7" % (pdb2gmx_params, input_pdb_full_filename,
																						structure_name,
																						structure_name,
																						structure_name), shell = True)
		# subprocess.call("gmx_d pdb2gmx -v -f %s -o %s.pdb -p %s.top -i %s.itp -water spc -inter -ignh" % (input_pdb_full_filename,
		# 																				structure_name,
		# 																				structure_name,
		# 																				structure_name), shell = True)
		subprocess.call("gmx_d editconf -f %s.pdb -o %s_box.gro -bt octahedron -d 1.4" % (structure_name, structure_name), shell = True)
		subprocess.call("gmx_d solvate -cp %s_box.gro -cs spc216.gro -o %s_wat.gro -p %s.top" % (structure_name, structure_name, structure_name), shell = True)
		# Minimisation
		subprocess.call("gmx_d grompp -f %s.mdp -c %s_wat.gro -p %s.top -o em.tpr" % (base_out_dirname + "/" + "em", structure_name, structure_name), shell = True)
		subprocess.call("gmx_d mdrun -s em.tpr -o em.trr -c em.gro -e em.edr -g em.log -nt 6", shell = True)
		# subprocess.call("gmx_d grompp -f %s_cg.mdp -c em.tpr -p %s.top -o em_cg.tpr" % (base_out_dirname + "/" + "em", structure_name), shell = True)
		# subprocess.call("gmx_d mdrun -s em_cg.tpr -o em_cg.trr -c em_cg.gro -e em_cg.edr -g em_cg.log", shell = True)
		subprocess.call("echo 14 | gmx_d trjconv -f em.trr -s em.tpr -o %s_minimised.pdb -pbc mol" % (structure_name), shell = True)
		subprocess.call("echo 7 | gmx_d energy -f em.edr -s em.tpr -o em.xvg", shell = True)
		sys.exit()
		# subprocess.call("echo gmx_d pdb2gmx -f %s -o %s -p %s.top -i %s.itp -water spc -inter -ignh" % (), shell = True)

		# pdb2gmx -f 3G0I_245A_variant_chain_a_repaired_res-renum.pdb -o 3G0I_protein.pdb -p 3G0I_protein.top -i 3G0I_protein_posre.itp -water spc -inter -ignh
		os.chdir("../")
		# sys.exit()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Energy minimise multiple pdbs')
	parser.add_argument('-i', '--input', help='Input directory', required=True)
	parser.add_argument('-o', '--output', help='Output directory', required=True)
	args = parser.parse_args(["-i", "../renamed_pdbs", "-o", "../test"])
	# args = parser.parse_args()
	_process_arguments(parser, args)

