while IFS= read -r link; do
  wget -P data_files/ $link
done
