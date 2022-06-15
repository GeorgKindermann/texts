if [ -z ${1+x} ]; then OUT="../pic/fichtenHerkuenfte.jpg";
else OUT=$1;
fi

if [ ! -f fichtenHerkuenfteORIG.jpg ]; then
    curl -o fichtenHerkuenfteORIG.jpg "https://bibdigital.rjb.csic.es/viewer/16441/download?file=HEG_Ill_Fl_Mitt_Eur_1_000276&crop=180.77735124760076,604.6951219512194,1242.4856046065258,765.1829268292682&CNT=0&type=img"
fi

convert -sampling-factor 4:2:0 -strip -quality 20 -interlace JPEG -colorspace Gray fichtenHerkuenfteORIG.jpg $OUT

