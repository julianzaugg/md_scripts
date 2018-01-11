
mkdir -p ../free_energy_results
rm -r ../free_energy_results/*
mkdir ../free_energy_results/dhdl_files
mkdir ../free_energy_results/last_frames

find . -wholename './3G0*_free_energy/[SR]_[SR]*/last_frames' -exec cp -r --parents \{\} ../free_energy_results/last_frames \;
find . -wholename './3G0*_free_energy/[SR]_[SR]*/dhdl_files' -exec cp -r --parents \{\} ../free_energy_results/dhdl_files \;
