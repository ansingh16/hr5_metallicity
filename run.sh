# Sort the files
source ~/miniconda3/bin/activate ICL  

trap 'break' ERR

for i in $(ls *.dat | tac); do
    [ -f "$i" ] || break
    echo "Sorting $i"
    sort -n "$i" -o "$i"
    snap=$(echo "${i%.*}")
    if [ "$snap" -le 300 ]; then
        #python ICL_work.py "$snap"
        python Background.py "$snap"
    else
        echo "Skipping $snap"
    fi 
done


