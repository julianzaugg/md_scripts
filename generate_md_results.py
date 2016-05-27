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


def _generate_grace_param_file(structure_name, ligand_name, replicate, xthreshold, ythreshold):
    with open("grace_params.txt", 'w') as fh:
        print >> fh, "WORLD XMAX %s" % xthreshold
        print >> fh, "WORLD YMAX %s" % ythreshold
        # print >> fh, "WORLD YMAX 5.0"
        print >> fh, "WORLD YMIN 0.0"
        print >> fh, "VIEW XMAX 1.2"
        print >> fh, "PAGE RESIZE 3000, 2000"
        print >> fh, 'TITLE "%s, %s, MD replicate - %s"' % (structure_name, ligand_name, replicate)
        # print >> fh, 'SUBTITLE "Ligand-Protein atom distances"'
        print >> fh, "LEGEND LOCTYPE WORLD"
        print >> fh, "LEGEND ON"
        print >> fh, "LEGEND 0.0, %s" % ythreshold
        print >> fh, "LEGEND CHAR SIZE .7"
        print >> fh, 'DEVICE "JPEG" DPI 200'
        print >> fh, 'HARDCOPY DEVICE "JPEG"'

def _get_max_value(filename):
    """
    Read an xvg file and find the maximum y value. Use this to automatically set a max y-axes value.
    """
    max_value = 0.0
    with open(filename, 'r') as fh:
        data = [line.strip().split()[1:] for line in fh.readlines()][50:] # assume the file is more than 50 lines
        for line in data:
            if max_value < max(map(float, line)):
                max_value = max(map(float, line))
    return max_value



OVERWRITE_IMAGES = True

for filename in os.listdir("."):
    if os.path.isdir(filename) and "free" not in filename:
        # TODO : collect individual rmsd filenames for binding pocket and full protein. Generate combined plots for these and write to file
        binding_pocket_rmsd_filenames = []
        full_structure_rmsd_filenames = []
        full_structure_rmsf_filenames = []
        base_dir = filename + "/"
        structure_name = filename.split("/")[-1]
        for ligand_name in ["S_GPE", "R_GPE"]:
            cur_dir = base_dir + ligand_name
            for MDdir in [dir for dir in os.listdir(cur_dir) if dir.startswith("MD")]:
                md_dir = cur_dir + "/" + MDdir + "/"
                print cur_dir + "/" + MDdir #, [(rep,step) for trr in os.listdir(md_dir) for rep,step in trr.split("_") if trr.endswith("trr")]
                trr_filenames = [mdfile for mdfile in os.listdir(md_dir) if mdfile.endswith(".trr")]
                full_trr_filesnames = [md_dir + trrfn for trrfn in trr_filenames]
                max_step = 1
                for trr_file in trr_filenames:
                    split_trr_filename = trr_file.split("_")
                    replicate, step = split_trr_filename[0][-1], split_trr_filename[1].split(".")[0]
                    if int(step) > max_step: max_step = int(step)
                time_threshold = max_step # Currently limiting max x-axis value. 
                if max_step > 40:
                    time_threshold = 40
                if not os.path.exists("%smd%s_1-%s.xtc" % (md_dir, replicate, max_step)):
                    subprocess.call("gmx_d trjcat -f %s -tu ns -xvg xmgrace -o %smd%s_1-%s.xtc" % (" ".join(full_trr_filesnames), md_dir, replicate, max_step), shell=True)
                
                if not os.path.exists("%smd%s_1-%s_noPBC.xtc" % (md_dir, replicate, max_step)):
                    subprocess.call("echo 16 | gmx_d trjconv -f %smd%s_1-%s.xtc -s %smd%s_1.tpr -o %smd%s_1-%s_noPBC.xtc -pbc mol -ur compact" % (md_dir,replicate,max_step, md_dir, replicate, md_dir,replicate,max_step), shell = True)
                
                if not os.path.exists("%smd%s_1-%s.pdb" % (md_dir, replicate, max_step)):
                    subprocess.call("echo 16 | gmx_d trjconv -f %smd%s_1-%s.xtc -s %smd%s_1.tpr -o %smd%s_1-%s.pdb -pbc mol -ur compact" % (md_dir,replicate,max_step, md_dir, replicate, md_dir,replicate,max_step), shell = True)
                
                # Extracts the last frame
                try:
                    last_frame_filename = "%s%s_%s_md%s_1-%s_last_frame.pdb" % (md_dir, structure_name, ligand_name, replicate, max_step)
                    if not os.path.exists(last_frame_filename):
                        print "Extracting last frame from PDB"
                        subprocess.call("sed -ne '/^TITLE/,/^ENDMDL/{/^TITLE/{x;d};H}' -e '${x;p}' %smd%s_1-%s.pdb >> %s" % (md_dir, replicate, max_step, last_frame_filename), shell = True)
                except:
                    print "Failed to extract last frame"
                
                try:
                    # Calculate the RMSD of carbon alpha in the whole structure
                    structure_rmsd_filename = "%s%s_%s_md%s_1-%s_rmsd.xvg" % (md_dir, structure_name, ligand_name, replicate, max_step)
                    full_structure_rmsd_filenames.append(structure_rmsd_filename)
                    if not os.path.exists(structure_rmsd_filename):
                        subprocess.call("echo 3 3 | gmx_d rms -s %smd%s_1.tpr -f %smd%s_1-%s_noPBC.xtc -o %s -tu ns -e %s" % (md_dir, replicate,  md_dir, replicate, max_step, structure_rmsd_filename, time_threshold), shell = True)
                except Exception, e:
                    print str(e)
                    print "Failed to generate structure rmsd"

                try:
                    # Calculate residue fluctuations (RMSF)
                    # structure_rmsf_filename = "%smd%s_1-%s_rmsf.xvg" % (md_dir, replicate, max_step)
                    structure_rmsf_filename = "%s%s_%s_md%s_1-%s_rmsf.xvg" % (md_dir, structure_name, ligand_name, replicate, max_step)
                    full_structure_rmsf_filenames.append(structure_rmsf_filename)
                    if not os.path.exists(structure_rmsf_filename):
                        # Carbon alpha rmsf
                        subprocess.call("echo 3 | gmx_d rmsf -s %smd%s_1.tpr -f %smd%s_1-%s_noPBC.xtc -o %s -e %s -res yes" % (md_dir, replicate,  md_dir, replicate, max_step, structure_rmsf_filename, time_threshold), shell = True)
                    # Generate RMSF plot
                    structure_rmsf_jpg_filename = "%s%s_%s_md%s_1-%s_rmsf.jpg" % (md_dir, structure_name, ligand_name, replicate, max_step)

                    # if not os.path.exists(structure_rmsf_jpg_filename) or OVERWRITE_IMAGES:
                    #     # Generate grace parameter file
                    #     y_threshold = str(_get_max_value(structure_rmsf_filename) * 1.1)
                    #     _generate_grace_param_file(structure_name, ligand_name, replicate, time_threshold, y_threshold)
                    #     subprocess.call("grace -nxy %s -autoscale xy -printfile %s -hardcopy -param grace_params.txt" % (structure_rmsf_filename, structure_rmsf_jpg_filename), shell = True)
                except Exception, e:
                    print str(e)
                    print "Failed to generate structure rmsf"

                try:
                    # Calculate the RMSD of residues around the binding pocket. Assumes the appropriate index file exists in base directory.
                    # binding_pocket_rmsd_filename = "%smd%s_1-%s_binding_pocket_rmsd.xvg" % (md_dir, replicate, max_step)
                    binding_pocket_rmsd_filename = "%s%s_%s_md%s_1-%s_binding_pocket_rmsd.xvg" % (md_dir, structure_name, ligand_name, replicate, max_step)
                    binding_pocket_rmsd_filenames.append(binding_pocket_rmsd_filename)
                    if not os.path.exists(binding_pocket_rmsd_filename):
                        # subprocess.call("echo 17 17 | gmx_d rms -s %smd%s_1.tpr -f %smd%s_1-%s_noPBC.xtc -o %smd%s_1-%s_binding_pocket_rmsd.xvg -tu ns -e %s -n binding_pocket_ndx.ndx" % (md_dir, replicate,  md_dir, replicate, max_step, md_dir, replicate,max_step, time_threshold), shell = True)
                        subprocess.call("echo 17 17 | gmx_d rms -s %smd%s_1.tpr -f %smd%s_1-%s_noPBC.xtc -o %s -tu ns -e %s -n binding_pocket_ndx.ndx" % (md_dir, replicate,  md_dir, replicate, max_step, binding_pocket_rmsd_filename, time_threshold), shell = True)

                    # Set maximum y axis value
                    y_threshold = str(_get_max_value(binding_pocket_rmsd_filename) * 1.1)
                    binding_pocket_jpg_filename = "%s%s_%s_md%s_1-%s_binding_pocket_rmsd.jpg" % (md_dir, structure_name, ligand_name, replicate, max_step)

                    # Grace must be installed to path for this to work
                    # if not os.path.exists(binding_pocket_jpg_filename) or OVERWRITE_IMAGES:
                        #Generate grace parameter file
                        # _generate_grace_param_file(structure_name, ligand_name, replicate, time_threshold, y_threshold)
                        # subprocess.call("grace -nxy %s -autoscale xy -printfile %s -hardcopy -param grace_params.txt" % (binding_pocket_rmsd_filename, binding_pocket_jpg_filename), shell = True)
                except Exception, e:
                    print str(e)

                try:
                    # Generates the average distances for the ligand-protein reaction critical atoms
                    # View with `xmgrace -nxy outputfilename.xvg
                    # distances_xvg_filename = "%smd%s_1-%s_ligand_protein_distances.xvg" % (md_dir, replicate, max_step)
                    distances_xvg_filename = "%s%s_%s_md%s_1-%s_ligand_protein_distances.xvg" % (md_dir, structure_name, ligand_name, replicate, max_step)
                    if not os.path.exists(distances_xvg_filename):
                        subprocess.call("gmx_d distance -f %smd%s_1-%s_noPBC.xtc -s %smd%s_%s.tpr -n ligand_protein_catalytic_atoms.ndx -oav %s -tu ns -select 17 18 19 20 -e %s" % \
                                                                (md_dir, replicate, max_step, md_dir, replicate, max_step, distances_xvg_filename, time_threshold), shell = True)
                    # Set maximum y axis value
                    y_threshold = str(math.ceil(_get_max_value(distances_xvg_filename)) * 1.1)

                    # Generate grace parameter file
                    _generate_grace_param_file(structure_name, ligand_name, replicate, 40, y_threshold)

                    # Generates jpgs of distance rmsd plots
                    distances_jpg_filename = "%s%s_%s_md%s_1-%s_distances.jpg" % (md_dir, structure_name, ligand_name, replicate, max_step)

                    # Grace must be installed to path for this to work
                    if not os.path.exists(distances_jpg_filename) or OVERWRITE_IMAGES:
                        subprocess.call("grace -nxy %s -autoscale xy -printfile %s -hardcopy -param grace_params.txt" % (distances_xvg_filename, distances_jpg_filename), shell = True)
                except Exception, e:
                    print str(e)
                    print "Failed to generate ligand-protein distances / distance plot"
        try:
            # Generate joint rmsd plots for binding pocket and full structure.
            # Generate joint rsmf plots for full structure
            # One plot for each pair of ligands and their duplicates
            binding_pocket_joint_jpg_filename = "%s%s_binding_pocket_rmsd.jpg" % (base_dir, structure_name)
            structure_joint_jpg_filename = "%s%s_rmsd.jpg" % (base_dir, structure_name)
            structure_joint_rmsf_jpg_filename = "%s%s_rmsf.jpg" % (base_dir, structure_name)
            
            if not os.path.exists(binding_pocket_joint_jpg_filename) or not os.path.exists(structure_joint_jpg_filename) or OVERWRITE_IMAGES:
                # Find the max value over all rmsd files
                binding_pocket_max_y = 0.0
                structure_max_y = 0.0
                structure_rmsf_max_y = 0.0
                for bp_filename in binding_pocket_rmsd_filenames:
                    value = _get_max_value(bp_filename)
                    if value > binding_pocket_max_y:
                        binding_pocket_max_y = value

                for struc_filename in full_structure_rmsd_filenames:
                    value = _get_max_value(struc_filename)
                    if value > structure_max_y:
                        structure_max_y = value

                # for struc_filename in full_structure_rmsf_filenames:
                #     value = _get_max_value(struc_filename)
                #     if value > structure_rmsf_max_y:
                #         structure_rmsf_max_y = value
                # binding_pocket_max_y = str((math.ceil(binding_pocket_max_y)) * 1.1)
                # structure_max_y = str((math.ceil(structure_max_y)) * 1.1)

                # Generate grace parameter file for binding pocket results
                _generate_grace_param_file(structure_name, "S/R-GPE", "1,2", 40, binding_pocket_max_y)

                subprocess.call("grace -nxy %s -autoscale xy -printfile %s -hardcopy -param grace_params.txt" % (" ".join(binding_pocket_rmsd_filenames), binding_pocket_joint_jpg_filename), shell = True)

                # Generate grace parameter file for binding pocket results
                _generate_grace_param_file(structure_name, "S/R-GPE", "1,2", 40, structure_max_y)
                subprocess.call("grace -nxy %s -autoscale xy -printfile %s -hardcopy -param grace_params.txt" % (" ".join(full_structure_rmsd_filenames), structure_joint_jpg_filename), shell = True)

                # TODO: joint rmsf plots
                # _generate_grace_param_file(structure_name, "S/R-GPE", "1,2", 40, structure_rmsf_max_y)
                subprocess.call("grace %s -printfile %s -hardcopy -hdevice JPEG" % (" ".join(full_structure_rmsf_filenames), structure_joint_rmsf_jpg_filename), shell = True)                


        except Exception, e:
            print str(e)
            print "Failed to generate joint rmsd plots"





