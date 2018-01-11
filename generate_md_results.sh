# Currently assumes single trajectory. Specifically for simulations where we have added ions.


################################
# Run from MD directory
STRUCTURE_NAME="3G02"
LIGAND_NAME="S_GPE"

# Generate index file for catalytic atoms. 13 should be the ligand, 22-28 the groups made.
echo -e "r 192 & a OD1
r 192 & a OD2
r 251 & a OH
r 314 & a OH
13 & a C8 
13 & a C9
13 & a O2
22 | 26
22 | 27
23 | 26
23 | 27
24 | 28
25 | 28
del 22
del 22
del 22
del 22
del 22
del 22
del 22
q" | gmx_d make_ndx -f PR/pr.gro -o ${STRUCTURE_NAME}_ligand_protein_catalytic_atoms_ions.ndx


# Generate index file for secondary structure residues. Combine with index file in EM_ions to retain defined groups
echo -e "r 24 | r 25 | r 26 | r 27 | r 28 | r 29 | r 30 | r 31 | r 32 | r 33 | r 34 | r 35 | r 36 | r 37 | r 38 | r 45 | r 46 | r 47 | r 48 | r 49 | r 57 | r 58 | r 59 | r 60 | r 61 | r 62 | r 63 | r 64 | r 65 | r 66 | r 67 | r 68 | r 69 | r 70 | r 71 | r 73 | r 74 | r 75 | r 76 | r 77 | r 78 | r 79 | r 80 | r 81 | r 82 | r 86 | r 87 | r 88 | r 89 | r 90 | r 91 | r 94 | r 95 | r 96 | r 97 | r 98 | r 99 | r 100 | r 101 | r 109 | r 110 | r 111 | r 112 | r 113 | r 114 | r 120 | r 121 | r 122 | r 123 | r 124 | r 125 | r 126 | r 127 | r 128 | r 129 | r 130 | r 131 | r 132 | r 133 | r 134 | r 135 | r 142 | r 143 | r 144 | r 145 | r 146 | r 147 | r 165 | r 166 | r 167 | r 168 | r 169 | r 170 | r 171 | r 172 | r 173 | r 174 | r 175 | r 176 | r 177 | r 178 | r 179 | r 180 | r 186 | r 187 | r 188 | r 189 | r 190 | r 192 | r 193 | r 194 | r 195 | r 196 | r 197 | r 198 | r 199 | r 200 | r 201 | r 202 | r 203 | r 204 | r 205 | r 208 | r 209 | r 210 | r 211 | r 212 | r 213 | r 226 | r 227 | r 228 | r 229 | r 230 | r 231 | r 232 | r 233 | r 234 | r 235 | r 236 | r 237 | r 238 | r 239 | r 240 | r 241 | r 242 | r 243 | r 244 | r 245 | r 246 | r 247 | r 248 | r 249 | r 250 | r 251 | r 252 | r 253 | r 254 | r 255 | r 256 | r 257 | r 258 | r 259 | r 260 | r 261 | r 262 | r 263 | r 264 | r 265 | r 266 | r 267 | r 268 | r 269 | r 270 | r 271 | r 272 | r 273 | r 274 | r 275 | r 276 | r 277 | r 278 | r 279 | r 280 | r 281 | r 282 | r 283 | r 284 | r 290 | r 291 | r 292 | r 293 | r 294 | r 295 | r 296 | r 297 | r 298 | r 299 | r 300 | r 301 | r 302 | r 303 | r 304 | r 305 | r 306 | r 307 | r 308 | r 309 | r 310 | r 311 | r 313 | r 314 | r 315 | r 316 | r 317 | r 318 | r 328 | r 329 | r 330 | r 331 | r 332 | r 336 | r 337 | r 338 | r 339 | r 340 | r 341 | r 342 | r 343 | r 344 | r 353 | r 354 | r 355 | r 356 | r 357 | r 358 | r 361 | r 362 | r 363 | r 364 | r 365 | r 366 | r 367 | r 368 | r 375 | r 376 | r 377 | r 378 | r 379 | r 380 | r 381 | r 382 | r 383 | r 384 | r 385 | r 386 | r 387 | r 388 | r 389 | r 390 | r 391 | r 392 | r 393 | r 394 | r 395 | r 396 & 4
q"| gmx_d make_ndx -f PR/pr.gro -o ${STRUCTURE_NAME}_SS_residues.ndx -n EM_ions/index.ndx

# Generate index file for residues not fluctuating >2.5 angstroms in original simulations
echo -e "4 & ! r 3-57 & ! r 66-73 & ! r 156-162 & ! r 220-236 & ! r 242-258 & ! r 318-327
q"| gmx_d make_ndx -f PR/pr.gro -o ${STRUCTURE_NAME}_not_fluctating_residues_removed_noncontinguous_ions.ndx

SKIP=false
if [ "$SKIP" = false ]; then
	for dir in ./MD[0-9]; do
		echo $dir
		cd $dir
		CURRENT_REPLICATE=$(echo $dir | grep -Eo '[0-9]+$')
		#echo $CURRENT_REPLICATE
		# NOT NEEDED IF ONLY SINGLE TRAJECTORY
		gmx_d trjcat -f md${CURRENT_REPLICATE}_[0-9].xtc -tu ns -xvg xmgrace -o md${CURRENT_REPLICATE}_1_out.xtc
		#cp md${CURRENT_REPLICATE}_1.xtc md${CURRENT_REPLICATE}_1_out.xtc
		# Fitting and centering system
		echo 22 | gmx_d trjconv -f md${CURRENT_REPLICATE}_1_out.xtc -s md${CURRENT_REPLICATE}_1.tpr -o md${CURRENT_REPLICATE}_1_out_noPBC_whole.xtc -pbc whole -ur compact -skip 1 -n ../EM_ions/index.ndx

		echo 22 | gmx_d trjconv -f md${CURRENT_REPLICATE}_1_out_noPBC_whole.xtc -s md${CURRENT_REPLICATE}_1.tpr -o md${CURRENT_REPLICATE}_1_out_noPBC_nojump.xtc -pbc nojump -ur compact -skip 1 -n ../EM_ions/index.ndx

		echo 1 22 | gmx_d trjconv -f md${CURRENT_REPLICATE}_1_out_noPBC_nojump.xtc -s md${CURRENT_REPLICATE}_1.tpr -o md${CURRENT_REPLICATE}_1_out_noPBC_center.xtc -center -ur compact -skip 1 -n ../EM_ions/index.ndx

		# Fit to backbone 
		#echo 4 22 | gmx_d trjconv -f md${CURRENT_REPLICATE}_1_out_noPBC_center.xtc -s md${CURRENT_REPLICATE}_1.tpr -o md${CURRENT_REPLICATE}_1_out_noPBC_fitted.xtc -ur compact -fit progressive -n ../EM_ions/index.ndx
		
		# Fit to backbone of secondary structure residues 
		echo 23 22 | gmx_d trjconv -f md${CURRENT_REPLICATE}_1_out_noPBC_center.xtc -s md${CURRENT_REPLICATE}_1.tpr -o md${CURRENT_REPLICATE}_1_out_noPBC_fitted.xtc -ur compact -fit progressive -n ../${STRUCTURE_NAME}_SS_residues.ndx
		

		# Simulation specific cluster files
		mkdir cluster_files
		echo 13 22 | gmx_d cluster -f md${CURRENT_REPLICATE}_1_out_noPBC_fitted.xtc -s md${CURRENT_REPLICATE}_1.tpr -cl cluster_files/md${CURRENT_REPLICATE}_1_out_clusters.pdb -g cluster_files/cluster.log -o cluster_files/md${CURRENT_REPLICATE}_1_out_clusters.xpm -dist cluster_files/${STRUCTURE_NAME}_${LIGAND_NAME}_md${CURRENT_REPLICATE}_1-1_out_cluster_dists.xvg -ntr cluster_files/${STRUCTURE_NAME}_${LIGAND_NAME}_md${CURRENT_REPLICATE}_1-1_out_cluster_trans.xvg -cutoff 0.3 -fit no -n ../EM_ions/index.ndx

		echo 22 | gmx_d trjconv -f md${CURRENT_REPLICATE}_1_out_noPBC_fitted.xtc -s md${CURRENT_REPLICATE}_1.tpr -o ${STRUCTURE_NAME}_${LIGAND_NAME}_md${CURRENT_REPLICATE}_1_out.pdb -n ../EM_ions/index.ndx
		# Extract last frame from PDB
		if [ "$(uname -s)" == "Darwin" ] # If the system is Mac OS X
		then
				mysed="gsed"
				mygrep="grep -FE"
		else
				mysed="sed"
				mygrep="grep -P"
		fi
		$mysed -ne '/^TITLE/,/^ENDMDL/{/^TITLE/{x;d};H}' -e '${x;p}' ${STRUCTURE_NAME}_${LIGAND_NAME}_md${CURRENT_REPLICATE}_1_out.pdb > ${STRUCTURE_NAME}_${LIGAND_NAME}_md${CURRENT_REPLICATE}_1_last_frame.pdb

		# Calculate the RMSD of backbone in the whole structure or defined secondary structure
		# 22 = backbone for residues with low fluctuation
		echo 22 22 | gmx_d rms -s ../PR/pr.gro -f md${CURRENT_REPLICATE}_1_out_noPBC_fitted.xtc -o ${STRUCTURE_NAME}_${LIGAND_NAME}_md${CURRENT_REPLICATE}_1-1_out_rmsd.xvg -n ../${STRUCTURE_NAME}_not_fluctating_residues_removed_noncontinguous_ions.ndx -tu ns
		#echo 4 4 | gmx_d rms -s ../PR/pr.gro -f md${CURRENT_REPLICATE}_1_out_noPBC_fitted.xtc -o ${STRUCTURE_NAME}_${LIGAND_NAME}_md${CURRENT_REPLICATE}_1-1_out_rmsd.xvg -tu ns


		# Calculate the RMSD of ligand based on whole structure backbone fit or defined secondary structure
		echo 4 13 | gmx_d rms -s ../PR/pr.gro -f md${CURRENT_REPLICATE}_1_out_noPBC_fitted.xtc -o ${STRUCTURE_NAME}_${LIGAND_NAME}_md${CURRENT_REPLICATE}_1-1_out_ligand_rmsd.xvg -tu ns

		# Generates the average distances for the ligand-protein reaction critical atoms
		gmx_d distance -f md${CURRENT_REPLICATE}_1_out_noPBC_fitted.xtc -s md${CURRENT_REPLICATE}_1.tpr -oav ${STRUCTURE_NAME}_${LIGAND_NAME}_md${CURRENT_REPLICATE}_1-1_out_ligand_protein_distances.xvg -tu ns -n ../${STRUCTURE_NAME}_ligand_protein_catalytic_atoms_ions.ndx  -select 22 23 24 25 26 27
		cd ../
	done
fi
echo $PWD
################################
# Generate single cluster for all simulation trajectories. Run from base structure directory, i.e. S_GPE
mkdir ${STRUCTURE_NAME}_cluster_files
XTC_FILES=$(find . -name 'md*_1_out_noPBC_fitted.xtc' | tr '\n' ' ')
echo $XTC_FILES
gmx_d trjcat -f $XTC_FILES -o ${STRUCTURE_NAME}_cluster_files/${STRUCTURE_NAME}_noPBC_fitted_for_clustering.xtc -cat
echo 13 22 | gmx_d cluster -f ${STRUCTURE_NAME}_cluster_files/${STRUCTURE_NAME}_noPBC_fitted_for_clustering.xtc -s PR/pr.gro -cl ${STRUCTURE_NAME}_cluster_files/${STRUCTURE_NAME}_clusters.pdb -g ${STRUCTURE_NAME}_cluster_files/${STRUCTURE_NAME}_clusters.log -o ${STRUCTURE_NAME}_cluster_files/${STRUCTURE_NAME}_clusters.xpm -dist ${STRUCTURE_NAME}_cluster_files/${STRUCTURE_NAME}_clusters_dists.xvg -ntr ${STRUCTURE_NAME}_cluster_files/${STRUCTURE_NAME}_clusters_trans.xvg -cutoff 0.3 -method gromos -fit no -n EM_ions/index.ndx

# RMSD of backbone of residues not fluctuating >2.5 angstroms in original simulations
echo 22 22 | gmx_d rms -s PR/pr.gro -f ${STRUCTURE_NAME}_cluster_files/${STRUCTURE_NAME}_noPBC_fitted_for_clustering.xtc -o ${STRUCTURE_NAME}_rmsd.xvg -tu ns -n ${STRUCTURE_NAME}_not_fluctating_residues_removed_noncontinguous_ions.ndx
