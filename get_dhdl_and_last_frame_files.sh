# Generate last frame pdbs and collect dhdl files from free energy calculations

mkdir last_frames
mkdir dhdl_files
rm -r dhdl_files/*
rm -r last_frames/*

for d in l_*[0-9]/ ; do
    cd $d
    mkdir ../dhdl_files/${d}
    find . -name "dhdl*.xvg" -exec cp \{\} ../dhdl_files/${d} \;
    LAMBDA_NAME=$(pwd | rev | cut -d / -f1 | rev | cut -d _ -f2 | tr -d .)
    MAXSTEP=$(ls md*.xtc | cut -d_ -f2 | cut -d. -f1)
    REPLICATE=$(ls md*.xtc | cut -d_ -f1 | cut -dd -f2)
    gmxcheck_d -f md${REPLICATE}_${MAXSTEP}.xtc 2> temp.txt
    LAST_FRAME=$(grep -Po "Last frame.+" temp.txt | grep -Po "(?<=time\s)[0-9].+" | cut -d . -f1)
    echo 13 17 |trjconv_d -f md${REPLICATE}_${MAXSTEP}.xtc -s md${REPLICATE}_${MAXSTEP}.tpr -o ../last_frames/md${REPLICATE}_${MAXSTEP}_l_${LAMBDA_NAME}_${LAST_FRAME}.pdb -dump ${LAST_FRAME} -center -ur compact -pbc mol
    cd ../
done