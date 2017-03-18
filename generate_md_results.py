"""
Script for my own work. Goes through molecular dynamics simulation results and performs further processing and 
generation of plots using xmgrace. Assumes specific naming scheme:
e.g., md1_1.trr, md1_2.trr....etc - md1 = md for replicate 1, _2 = time step.

This was originally written for liganding binding MD analysis (hence the `S/R GPE').
Assumes a number of index files were available.
Folder heirarchy was - `script_base_directory -> Structure_directory -> Ligand_directory -> MD_replicate_directory'

I would just use this script for `inspiration' in writing your own.
"""

import os, sys

import subprocess
import math


# OVERWRITE_IMAGES = False
STARTING_FRAME_RMSF = 30000
SKIP_AMOUNT = 1

OVERWRITE_RMSF = False
OVERWRITE_RMSD = False
OVERWRITE_LIGAND_RMSD = False
OVERWRITE_BINDING_POCKET_RMSD = False
OVERWRITE_LIGAND_PROTEIN_DISTANCES = False
OVERWRITE_PDB = False
OVERWRITE_LASTFRAME = False
OVERWRITE_NOPBC = False
OVERWRITE_CLUSTERING = False
SS_BACKBONE = True # use secondary structure backbone residues only. Requires defined index file

# OVERWRITE_RESULTS_FOR_STRUCTURES = ["poly3_FDVFFFVV", "poly3_FDVFVGDV", "poly3_FVVFFCLV"]
OVERWRITE_RESULTS_FOR_STRUCTURES = [] #["3G0I", "3G02"]

# OVERWRITE_LIG_DISTANCES_FOR_STRUCTURES = ["poly3_FDVFFFVV", "poly3_FDVFVGDV", "poly3_FVVFFCLV"]
OVERWRITE_LIG_DISTANCES_FOR_STRUCTURES = []
OVERWRITE_CLUSTERING_FOR_STRUCTURES = ["3G02"]


TAG = "SS_fitted"

for filename in os.listdir("."):
    if os.path.isdir(filename) and "free" not in filename:
        # TODO : collect individual rmsd filenames for binding pocket and full protein. Generate combined plots for these and write to file
        binding_pocket_rmsd_filenames = []
        full_structure_rmsd_filenames = []
        full_structure_rmsf_filenames = []
        base_dir = filename + "/"
        structure_name = filename.split("/")[-1]
        if structure_name != "3G0I" and structure_name != "3G02": continue
        #if structure_name != "3G02": continue
        #for ligand_name in ["S_GPE", "R_GPE"]:
        for ligand_name in ["S_GPE_crystal_water_gromacs_333", "R_GPE_crystal_water_gromacs_333"]:
            if structure_name != "3G02" and ligand_name != "R_GPE_crystal_water_gromacs_333": continue
            cur_dir = base_dir + ligand_name
            xtc_list_for_clustering = []
            for MDdir in [dir for dir in os.listdir(cur_dir) if dir.startswith("MD")]:
                md_dir = cur_dir + "/" + MDdir + "/"
                print cur_dir + "/" + MDdir #, [(rep,step) for trr in os.listdir(md_dir) for rep,step in trr.split("_") if trr.endswith("trr")]
                # Get trajectory filenames
                trr_filenames = [mdfile for mdfile in os.listdir(md_dir) if mdfile.endswith(".trr")]
                full_trr_filesnames = [md_dir + trrfn for trrfn in trr_filenames]
                max_step = 1
                if not os.path.exists("%scluster_files/" % md_dir):
                    subprocess.call("mkdir %scluster_files/" % md_dir, shell = True)
                for trr_file in trr_filenames:
                    split_trr_filename = trr_file.split("_")
                    if "backup" in trr_file: continue
                    replicate, step = split_trr_filename[0][-1], split_trr_filename[1].split(".")[0]
                    if int(step) > max_step: max_step = int(step)
                time_threshold = max_step # Currently limiting max x-axis value. 
                if max_step > 50:
                    time_threshold = 50
                
                # Join the separate trajectory files together
                # We remove the trr -> xtc file for space reasons. If the noPBC xtc does not exist, we will re-concatenate
                if not os.path.exists("%smd%s_1-%s.xtc" % (md_dir, replicate, max_step)) and \
                not os.path.exists("%smd%s_1-%s_noPBC.xtc" % (md_dir, replicate, max_step)) \
                    or structure_name in OVERWRITE_RESULTS_FOR_STRUCTURES:
                    subprocess.call("gmx_d trjcat -f %s -tu ns -xvg xmgrace -o %smd%s_1-%s.xtc" % (" ".join(full_trr_filesnames), md_dir, replicate, max_step), shell=True)

                if not os.path.exists("%smd%s_1-%s_noPBC_fitted.xtc" % (md_dir, replicate, max_step)) or OVERWRITE_NOPBC or structure_name in OVERWRITE_RESULTS_FOR_STRUCTURES:
                    subprocess.call("echo 16 | gmx_d trjconv -f %smd%s_1-%s.xtc -s %smd%s_1.tpr -o %smd%s_1-%s_noPBC_whole.xtc -pbc whole -ur compact -skip %i" % (md_dir,replicate,max_step, md_dir, replicate, md_dir,replicate,max_step, SKIP_AMOUNT), shell = True)
                    subprocess.call("echo 16 | gmx_d trjconv -f %smd%s_1-%s_noPBC_whole.xtc -s %smd%s_1.tpr -o %smd%s_1-%s_noPBC_nojump.xtc -ur compact -pbc nojump" % (md_dir,replicate,max_step, md_dir, replicate, md_dir,replicate,max_step), shell = True)
                    # 1 = center to protein, 16 - output non-water
                    subprocess.call("echo 1 16 | gmx_d trjconv -f %smd%s_1-%s_noPBC_nojump.xtc -s %smd%s_1.tpr -o %smd%s_1-%s_noPBC_center.xtc  -ur compact -center" % (md_dir,replicate,max_step, md_dir, replicate, md_dir,replicate,max_step), shell = True)
                    # 4/17 = fit to backbone/secondary structure backbone, 16 - output non-water
                    if SS_BACKBONE:
                        subprocess.call("echo 17 16 | gmx_d trjconv -f %smd%s_1-%s_noPBC_center.xtc -s %smd%s_1.tpr -o %smd%s_1-%s_noPBC_fitted.xtc -ur compact -fit progressive -n %s_SS_residues.ndx" % (md_dir,replicate,max_step, md_dir, replicate, md_dir,replicate,max_step, structure_name), shell = True)
                    else:
                        subprocess.call("echo 4 16 | gmx_d trjconv -f %smd%s_1-%s_noPBC_center.xtc -s %smd%s_1.tpr -o %smd%s_1-%s_noPBC_fitted.xtc -ur compact -fit progressive" % (md_dir,replicate,max_step, md_dir, replicate, md_dir,replicate,max_step), shell = True)
                    # subprocess.call("rm %smd%s_1-%s.xtc" % (md_dir, replicate, max_step), shell = True)
                if structure_name == "3G02" and replicate == "1" and "R_GPE" in ligand_name:
                    subprocess.call("echo 17 16 | gmx_d trjconv -f %smd%s_1-%s_noPBC_center.xtc -s %smd%s_1.tpr -o %smd%s_1-%s_noPBC_fitted.xtc -ur compact -fit progressive -n %s_SS_residues.ndx -e 35000" % (md_dir,replicate,max_step, md_dir, replicate, md_dir,replicate,max_step, structure_name), shell = True)
                xtc_list_for_clustering.append(("%smd%s_1-%s_noPBC_fitted.xtc" % (md_dir, replicate, max_step), max_step))

                if not os.path.exists("%scluster_files/md%s_1-%s_clusters.pdb" % (md_dir,replicate,max_step)) or OVERWRITE_CLUSTERING or structure_name in OVERWRITE_CLUSTERING_FOR_STRUCTURES:
                    subprocess.call("echo 13 16 | gmx_d cluster -f %smd%s_1-%s_noPBC_fitted.xtc -s %smd%s_1.tpr -cl %scluster_files/md%s_1-%s_clusters.pdb -g %scluster_files/cluster.log -o %scluster_files/md%s_1-%s_clusters.xpm -dist %scluster_files/md%s_1-%s_cluster_dists.xvg -ntr %scluster_files/md%s_1-%s_cluster_trans.xvg -cutoff 0.75 -method gromos -fit no"  % (md_dir, replicate, max_step, md_dir, replicate, md_dir, replicate, max_step, md_dir, md_dir, replicate, max_step, md_dir, replicate, max_step, md_dir, replicate, max_step), shell = True)

                # Generate pdb from trajectory
                if not os.path.exists("%smd%s_1-%s.pdb" % (md_dir, replicate, max_step)) or OVERWRITE_PDB or structure_name in OVERWRITE_RESULTS_FOR_STRUCTURES:
                    subprocess.call("echo 16 | gmx_d trjconv -f %smd%s_1-%s_noPBC_fitted.xtc -s %smd%s_1.tpr -o %smd%s_1-%s.pdb" % (md_dir,replicate,max_step, md_dir, replicate, md_dir,replicate, max_step), shell = True)
                # Extracts the last frame
                try:
                    last_frame_filename = "%s%s_%s_md%s_1-%s_last_frame.pdb" % (md_dir, structure_name, ligand_name, replicate, max_step)
                    if not os.path.exists(last_frame_filename) or OVERWRITE_LASTFRAME or structure_name in OVERWRITE_RESULTS_FOR_STRUCTURES:
                        print "Extracting last frame from PDB"
                        subprocess.call("sed -ne '/^TITLE/,/^ENDMDL/{/^TITLE/{x;d};H}' -e '${x;p}' %smd%s_1-%s.pdb > %s" % (md_dir, replicate, max_step, last_frame_filename), shell = True)
                except:
                    print "Failed to extract last frame"
                
                try:
                    # Calculate the RMSD of backbone in the whole structure or defined secondary structure
                    structure_rmsd_filename = "%s%s_%s_md%s_1-%s_rmsd.xvg" % (md_dir, structure_name, ligand_name, replicate, max_step)
                    full_structure_rmsd_filenames.append(structure_rmsd_filename)
                    if not os.path.exists(structure_rmsd_filename) or OVERWRITE_RMSD or structure_name in OVERWRITE_RESULTS_FOR_STRUCTURES:
                        if SS_BACKBONE:
                            subprocess.call("echo 17 17 | gmx_d rms -s %s../PR/pr_corrected.gro -f %smd%s_1-%s_noPBC_fitted.xtc -o %s -tu ns -e %i -n %s_SS_residues.ndx" % (md_dir,  md_dir, replicate, max_step, structure_rmsd_filename, time_threshold * 1000, structure_name), shell = True)
                        else:
                            subprocess.call("echo 4 4 | gmx_d rms -s %s../PR/pr_corrected.gro -f %smd%s_1-%s_noPBC_fitted.xtc -o %s -tu ns -e %i" % (md_dir,  md_dir, replicate, max_step, structure_rmsd_filename, time_threshold * 1000), shell = True)
                except Exception, e:
                    print str(e)
                    print "Failed to generate structure rmsd"

                try:
                    # Calculate the RMSD of ligand based on whole structure backbone fit or defined secondary structure
                    structure_ligand_rmsd_filename = "%s%s_%s_md%s_1-%s_ligand_rmsd.xvg" % (md_dir, structure_name, ligand_name, replicate, max_step)
                    full_structure_rmsd_filenames.append(structure_rmsd_filename)
                    if not os.path.exists(structure_ligand_rmsd_filename) or OVERWRITE_LIGAND_RMSD or structure_name in OVERWRITE_RESULTS_FOR_STRUCTURES:
                        if SS_BACKBONE:
                            subprocess.call("echo 17 13 | gmx_d rms -s %s../PR/pr_corrected.gro -f %smd%s_1-%s_noPBC_fitted.xtc -o %s -tu ns -e %i -n %s_SS_residues.ndx" % (md_dir,  md_dir, replicate, max_step, structure_ligand_rmsd_filename, time_threshold * 1000, structure_name), shell = True)
                        else:
                            subprocess.call("echo 4 13 | gmx_d rms -s %s../PR/pr_corrected.gro -f %smd%s_1-%s_noPBC_fitted.xtc -o %s -tu ns -e %i" % (md_dir,  md_dir, replicate, max_step, structure_ligand_rmsd_filename, time_threshold * 1000), shell = True)
                except Exception, e:
                    print str(e)
                    print "Failed to generate ligand rmsd"

                # try:
                #     # Calculate residue fluctuations (RMSF)
                #     structure_rmsf_filename = "%s%s_%s_md%s_1-%s_rmsf.xvg" % (md_dir, structure_name, ligand_name, replicate, max_step)
                #     full_structure_rmsf_filenames.append(structure_rmsf_filename)
                #     if not os.path.exists(structure_rmsf_filename) or OVERWRITE_RMSF or structure_name in OVERWRITE_RESULTS_FOR_STRUCTURES:
                #         # Carbon alpha rmsf makes more sense for rmsf
                #         subprocess.call("echo 3 | gmx_d rmsf -s %s../PR/pr_corrected.gro -f %smd%s_1-%s_noPBC_fitted.xtc -o %s -e %i -res yes -b %i" % (md_dir,  md_dir, replicate, max_step, structure_rmsf_filename, time_threshold * 1000, STARTING_FRAME_RMSF), shell = True)
                # except Exception, e:
                #     print str(e)
                #     print "Failed to generate structure rmsf"

                # try:
                #     # Calculate the RMSD of residues around the binding pocket. Assumes the appropriate index file exists in base directory.
                #     # binding_pocket_rmsd_filename = "%smd%s_1-%s_binding_pocket_rmsd.xvg" % (md_dir, replicate, max_step)
                #     binding_pocket_rmsd_filename = "%s%s_%s_md%s_1-%s_binding_pocket_rmsd.xvg" % (md_dir, structure_name, ligand_name, replicate, max_step)
                #     binding_pocket_rmsd_filenames.append(binding_pocket_rmsd_filename)
                #     if not os.path.exists(binding_pocket_rmsd_filename) or OVERWRITE_BINDING_POCKET_RMSD or structure_name in OVERWRITE_RESULTS_FOR_STRUCTURES:
                #         subprocess.call("echo 17 17 | gmx_d rms -s %s../PR/pr_corrected.gro -f %smd%s_1-%s_noPBC_fitted.xtc -o %s -tu ns -e %s -n %s_binding_pocket_ndx.ndx" % (md_dir,  md_dir, replicate, max_step, binding_pocket_rmsd_filename, time_threshold * 1000, structure_name), shell = True)
                # except Exception, e:
                #     print str(e)

                try:
                    # Generates the average distances for the ligand-protein reaction critical atoms
                    distances_xvg_filename = "%s%s_%s_md%s_1-%s_ligand_protein_distances.xvg" % (md_dir, structure_name, ligand_name, replicate, max_step)
                    if not os.path.exists(distances_xvg_filename) or structure_name in OVERWRITE_RESULTS_FOR_STRUCTURES or structure_name in OVERWRITE_LIG_DISTANCES_FOR_STRUCTURES:
                        subprocess.call("gmx_d distance -f %smd%s_1-%s_noPBC_fitted.xtc -s %smd%s_%s.tpr -n %s_ligand_protein_catalytic_atoms.ndx -oav %s -tu ns -select 17 18 19 20 21 22 -e %s" % \
                                                                (md_dir, replicate, max_step, md_dir, replicate, max_step, structure_name, distances_xvg_filename, time_threshold * 1000), shell = True)
                except Exception, e:
                    print str(e)
                    print "Failed to generate ligand-protein distances / distance plot"
            # Generate single cluster for all trajectories
            if not os.path.exists("%s/%s_cluster_files/" % (cur_dir, structure_name)):
                subprocess.call("mkdir %s/%s_cluster_files/" % (cur_dir, structure_name), shell = True)
            str_xtc_list = " ".join([x[0] for x in xtc_list_for_clustering])
            str_xtc_frame_steps = " ".join([str(x[1] * 1000) for x in xtc_list_for_clustering])
            print "clustering these files - ", str_xtc_list
            if not os.path.exists("%s/%s_cluster_files/%s_noPBC_fitted_for_clustering.xtc" % (cur_dir, structure_name, structure_name)) or OVERWRITE_CLUSTERING or structure_name in OVERWRITE_CLUSTERING_FOR_STRUCTURES:
                frame_steps = []
                subprocess.call("gmx_d trjcat -f %s -o %s/%s_cluster_files/%s_noPBC_fitted_for_clustering.xtc -cat" % (str_xtc_list, cur_dir, structure_name,structure_name), shell = True)

            if not os.path.exists("%s/%s_cluster_files/%s_clusters.pdb" % (cur_dir,structure_name,structure_name)) or OVERWRITE_CLUSTERING:
                subprocess.call("echo 13 16 | gmx_d cluster -f %s/%s_cluster_files/%s_noPBC_fitted_for_clustering.xtc  -s %s/PR/pr_corrected.gro -cl %s/%s_cluster_files/%s_clusters.pdb -g %s/%s_cluster_files/%s_cluster.log -o %s/%s_cluster_files/%s_clusters.xpm -dist %s/%s_cluster_files/%s_cluster_dists.xvg -ntr %s/%s_cluster_files/%s_cluster_trans.xvg -cutoff 0.75 -method gromos -fit no" %  (cur_dir,structure_name,structure_name, cur_dir, cur_dir, structure_name, structure_name,cur_dir, structure_name, structure_name,cur_dir, structure_name, structure_name, cur_dir, structure_name, structure_name, cur_dir, structure_name, structure_name), shell = True)




