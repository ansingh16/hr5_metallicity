# Sort the files
source ~/miniconda3/bin/activate ICL  

trap 'break' ERR

cd ./snapfiles/

file_list=$(ls *.dat | tac)
sorted_list=$(echo "$file_list" | sort -t'.' -k1 -n -r | tr '\n' ' ')

for i in $sorted_list; do
    [ -f "$i" ] || break
    echo "Sorting $i"
    sort -n "$i" -o "$i"
    snap=$(echo "${i%.*}")
    if [ "$snap" -le 300 ]; then
        python ../Code/ICL_work.py "$snap"
        #python ../Code/Background.py "$snap"
    else
        echo "Skipping $snap"
    fi 
done

# source ~/miniconda3/bin/activate ICL  

# cd ./Code/
# trap 'break' ERR

# find ../snapfiles/ -name "*.dat" -print0 | sort -zrV | xargs -0 -n 1 -P 4 sh -c '
#     # Run the command for each file
#     stem=$(basename "$1" .dat)
#     echo Processing $stem
#     python ICL_work.py "$stem" >> "$stem.log"
# ' sh