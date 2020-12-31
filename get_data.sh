mkdir data_files
mkdir data_files/tsv
mkdir data_files/rds

while IFS= read -r link; do
    
    if [[ $link == *".RDS"* ]]
    then
        wget -P data_files/rds $link
    else
        wget -P data_files/tsv $link
    fi

done < links.txt
