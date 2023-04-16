#!/usr/bin/env bash

export TEMP=/tmp
for i in data
do
    if [ ! -d ${i} ]; then
        mkdir -p ${i}
    fi
done

if [ ! -d venv ]; then
    echo Set up python3 virtual environment
    python3 -m venv venv
    source venv/bin/activate
    pip3 install --upgrade pip
    pip3 install -r requirements.txt
else
    source venv/bin/activate
fi

if [ ! -h momepy2 ]; then
    git clone https://github.com/anisotropi4/momepy.git
    ln -s momepy/momepy momepy2
fi

if [ ! -d jay ]; then
    git clone https://github.com/anisotropi4/jay.git
fi

FILESTUB=network-model
if [ ! -s data/${FILESTUB}-simple.gpkg ]; then
    URI=https://github.com/openraildata/network-rail-gis/releases/download/20230317-01
    if [ ! -s data/${FILESTUB}.gpkg ]; then
        curl -Lo data/${FILESTUB}.gpkg ${URI}/${FILESTUB}.gpkg
    fi
    (cd jay; ./simplify.sh ../data/${FILESTUB}.gpkg)
    ln jay/output/${FILESTUB}-simple.gpkg data/${FILESTUB}-simple.gpkg
fi

FILESTUB=great-britain-rail
if [ ! -s data/${FILESTUB}-simple.gpkg ]; then
    URI=https://github.com/anisotropi4/magpie/blob/master/great-britain-rail.gpkg?raw=true
    if [ ! -s data/${FILESTUB}.gpkg ]; then
        curl -Lo data/${FILESTUB}.gpkg ${URI}/${FILESTUB}.gpkg
    fi
    (cd jay; ./simplify.sh ../data/${FILESTUB}.gpkg)
    ln jay/output/${FILESTUB}-simple.gpkg data/${FILESTUB}-simple.gpkg
fi

if [ ! -s linetrack.gpkg ]; then
    ./trackcheck.py
fi

if [ ! -s linetrack.gpkg ]; then
    ./trackcheck.py
fi

if [ ! -s outputx.gpkg ]; then
    ./gettrack.py
fi

