if [ -z ${1+x} ]; then OUT="../pic/waelderDesOstalpenraumes.pdf";
else OUT=$1;
fi

if [ ! -f waelderDesOstalpenraumesORIG.jpg ]; then
    curl -o waelderDesOstalpenraumesORIG.jpg "https://anno.onb.ac.at/cgi-content/annoshow-plus?call=cbf|1977|0004|00000152|00000001|jpg||O|"
fi
    
tmp=$(mktemp -d)

convert waelderDesOstalpenraumesORIG.jpg -crop 5000x3500+170+240 $tmp/crop.jpg
gdal_translate -of GTiff -gcp 599.982 1332.5 254558 5.28499e+06 -gcp 2245.71 1147.59 464694 5.31029e+06 -gcp 3179.02 445.799 583416 5.39663e+06 -gcp 4803.47 1132.69 787250 5.31032e+06 -gcp 2765.35 3067.82 533530 5.0687e+06 -a_nodata 255 $tmp/crop.jpg $tmp/cropProj.jpg
gdalwarp -r bilinear -tps -co COMPRESS=DEFLATE  -t_srs EPSG:31286 $tmp/cropProj.jpg $tmp/cropProj.tif
jbig2 -s -p -T 150 -v $tmp/cropProj.tif && python2 /usr/local/bin/pdf.py output >$tmp/waelderDesOstalpenraumesBW.pdf
rm output.0000 output.sym

if [ ! -f gadm40_AUT_1.kmz ]; then
    wget https://geodata.ucdavis.edu/gadm/gadm4.0/kmz/gadm40_AUT_1.kmz
fi
ogr2ogr -t_srs EPSG:31286 -f csv -lco GEOMETRY=AS_WKT /vsistdout/ gadm40_AUT_1.kmz | sed -r '1d; s/,/\n/g; s/[a-zA-Z"(]//g; s/\).*/\nNA NA/g' >$tmp/border.txt

Rscript --vanilla plotBordersAt.r $tmp

pdfcrop $tmp/border.pdf
pdftk $tmp/border-crop.pdf background $tmp/waelderDesOstalpenraumesBW.pdf output $OUT

rm -r $tmp

