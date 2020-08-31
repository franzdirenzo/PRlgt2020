touch MYtest
for X in $(seq 48)
do
    echo "seed $RANDOM$RANDOM" >> MYtest
done
