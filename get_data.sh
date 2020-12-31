mkdir data_files
mkdir data_files/tsv
mkdir data_files/rds

while IFS= read -r links; do
    
    if [[ $link == *".RDS"* ]]
    then
        wget -P data_files/tsv $link
    else
        wget -P data_file/rds $link
    fi

done < links.txt
