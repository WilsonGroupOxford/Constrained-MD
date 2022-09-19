COUNTER=0
for file in fft_*.pdf; do
    PADDEDCOUNTER=$(printf "%02d\n" $COUNTER)
    cp $file ./output/fft_$PADDEDCOUNTER.pdf
    let COUNTER++
    /usr/bin/convert -density 300 ./output/fft_$PADDEDCOUNTER.pdf -resize 1920x1080 ./output/fft_$PADDEDCOUNTER.png
    rm ./output/fft_$PADDEDCOUNTER.pdf
done
rm ./output/test.mp4
ffmpeg -r 2 -i ./output/fft_%02d.png ./output/test.mp4
