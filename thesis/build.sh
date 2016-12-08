#!/bin/bash
set -e

texname="main"

if [ ! -e $texname".tex" ]; then
    echo -e "\e[31m$texname.tex doesn't exist. Exiting...\e[0m"
    exit 1
fi

echo -n "Cleaning last build... "
rm -rf .build
echo -e "\e[32mOK\e[0m"

set +e

echo -n "Recreating build directory and copying files... "
mkdir -p .build
for file in *; do
    [ $file == "build.sh" ] || cp -r $file .build
	[ $? -ne 0 ] && echo -e "\e[31mFAILED\e[0m" && exit 1
done
echo -e "\e[32mOK\e[0m"

cd .build

# generate all plots
if ls gen/* &>/dev/null; then
    for f in *.plt; do
        echo -n "Generating $f... "
        chmod u+x $f
        ./$f 1>/dev/null 2>/dev/null && echo -e "\e[32mOK\e[0m" || echo -e "\e[31mFAILED\e[0m" && exit 1
    done
fi


echo -n "Generating $texname... "
pdflatex -interaction=nonstopmode $texname.tex 1>/dev/null
if [ $? -ne 0 ]; then
    echo -e "\e[31mFAILED\e[0m"
    rm -f $texname.aux
    pdflatex $texname.tex
    exit 1
else
    echo -e "\e[32mOK\e[0m"
fi

set -e

echo -n "Second pass... "
pdflatex -interaction=nonstopmode $texname.tex 1>/dev/null
echo -e "\e[32mOK\e[0m"

cd ..
cp .build/$texname.pdf ./$texname.pdf
echo -e "Output copied to \e[34m"`pwd`"/$texname.pdf\e[0m"
