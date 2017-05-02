set -u
set -e
dstd=tmp_input_data
mkdir -p "$dstd"
while read url ; do
    filename=$( echo "$url" | perl -MURI -le 'chomp($url = <>); print URI->new($url)->path')
    filename=$(basename "$filename")
    storage_path="$dstd/$filename"
    if [ -f $storage_path ] ; then
        echo "File exists. Skipping $storage_path"
    else
        echo "Downloading $storage_path"
        curl -o "$storage_path" "$url"
    fi
done < external_input_data.txt

