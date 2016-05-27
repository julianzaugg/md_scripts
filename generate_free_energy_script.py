

import sys, os

import argparse



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-s', '--start', help='Start lambda', required=False, default=0.0)
	parser.add_argument('-e', '--end', help='End lambda', required=False, default=1.0)
	parser.add_argument('-st', '--stepsize', help='Step size', required=False, default=0.1)
	parser.add_argument('-o', '--output', help='Output filename', required=False, default="free_energy_script.sh")
	parser.add_argument('-top', '--topology_name', help='Topology filename', required=True)
	args = parser.parse_args()
	
	with open(args.output, 'w') as fh:
		current_step = args.start
		while current_step < args.end:
			if current_step == 0.0:
				print >> fh, "cd l_%0.2f" % (current_step)
				print >> fh, "gmx_d grompp -f MD.mdp -c ../em.gro -p ../%s.top -o ./md1_1.tpr -maxwarn 2" % args.topology_name
				print >> fh, "gmx_d mdrun -s ./md1_1.tpr -o ./md1_1.trr -c ./md1_1.gro -e md1_1.edr -g ./md1_1.log"
				current_step += args.stepsize
				continue

			print >> fh, "cd ../l_%0.2f" % (current_step)
			print >> fh, "gmx_d grompp -f MD.mdp -c ../l_%0.2f/md1_1.gro -p ../%s.top -o ./md1_1.tpr -maxwarn 2" % (current_step - args.stepsize,args.topology_name)
			print >> fh, "gmx_d mdrun -s ./md1_1.tpr -o ./md1_1.trr -c ./md1_1.gro -e md1_1.edr -g ./md1_1.log"
			current_step += args.stepsize
