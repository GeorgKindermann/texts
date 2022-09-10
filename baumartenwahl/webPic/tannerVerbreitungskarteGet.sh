if [ -z ${1+x} ]; then OUT="../pic/tannerVerbreitungskarte.jpg";
else OUT=$1;
fi

if [ ! -f tannerVerbreitungskarteORIG.png ]; then
    curl -o tannerVerbreitungskarteORIG.png "https://upload.wikimedia.org/wikipedia/commons/thumb/9/9c/Abies_alba_distribution_map.svg/800px-Abies_alba_distribution_map.svg.png"
fi

convert -sampling-factor 4:2:0 -strip -quality 30 -interlace JPEG -colorspace sRGB tannerVerbreitungskarteORIG.png $OUT
