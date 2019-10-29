count=0
mypath=$1
for files in $mypath/*; do
    ./runJetFinding -hard $files -output output/$2-$count.root -nev $3 -hardtype $4 -tag $5 >&log/$2-$count.log &
    let count+=1
    [[ $((count%10)) -eq 0 ]] && wait
done
