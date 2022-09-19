PI=$(echo "scale=5; 4*a(1)" | bc -l)
NUM_STEPS=250
STEP_SIZE=$(echo "scale=5; 4*a(1)/$NUM_STEPS" | bc -l)

echo "Taking $NUM_STEPS steps of size $STEP_SIZE from 0 to $PI"
for angle in $(seq 0.00 $STEP_SIZE $PI); do
    echo $angle
    python3 "./generate_angle.py" "0.0" "$angle"
    ./Constrained-MD > ./output.txt
    python3 "./plot_excitement.py"
    precision_angle=$(printf "%.3f" $angle)
    mv "./output.xyz" "./output_$precision_angle.xyz"
    mv "./fft.pdf" "./fft_$precision_angle.pdf"
    mv "./frequencies.txt" "./frequencies_$precision_angle.txt"
done
