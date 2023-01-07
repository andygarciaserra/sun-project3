start=$(date +%s)
python3 plots.py
end=$(date +%s)
echo "Elapsed Time: $(($end-$start)) seconds"
