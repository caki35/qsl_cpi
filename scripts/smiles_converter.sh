#!bash
mkdir smiles_output
cwd=$(pwd)
cd $1
counter=1
for i in $(ls |grep mol)
do
    echo $('C:\Program Files\ChemAxon\JChemSuite\bin/molconvert.exe' -Y smiles:a ${i}) >> ${cwd}/smiles_output/${i%.*}.smi
    echo -en "\r$counter. compound is converted"
    (( counter += 1 ))
done
