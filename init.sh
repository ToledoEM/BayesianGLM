#!bin/sh

dir=./OUTPUT; [ ! -e $dir ] && mkdir -p $dir
dir=./DATA; [ ! -e $dir ] && mkdir -p $dir

cd DATA

wget https://github.com/linnarsson-lab/ipynb-lamanno2016/archive/master.zip
unzip master.zip
mv ipynb-lamanno2016-master/data/*.cef .
mv ipynb-lamanno2016-master/data/*.tsv .
rm master.zip
rm -r ipynb-lamanno2016-master

