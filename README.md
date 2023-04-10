# Woodpecker
This is a way to overlay and combine two linear geographic network models.

## Overview
The centre-line track-model is a dataset published via an Environmental Information Regulations (EIR) [request](https://www.whatdotheyknow.com/request/geospatial_data) by [Network Rail](https://www.networkrail.co.uk/who-we-are/transparency-and-ethics/freedom-of-information-foi) under the Open Government Licence:

Open Street Map contains current and historic mainland Britain railway data [here](https://www.openstreetmap.org/) published under the [Open Database License](https://opendatacommons.org/licenses/odbl/1-0/).

Run the following script:

    $ ./run.sh

This will downloads dependencies, centre-line and OSM data, creates the network overlay and combines this in a `outputx.gpkg` repository.

## Overlay and combine approach
The approach used to overlay and combine data is:

1. create a cropped OSM and centre-line track-model LineString segments within scaled hexagons
2. create line-extended centre-line track-model rectangular polygons within a given hexagon
3. test if the OSM LineStrings are within these rectangular polygons
4. any OSM LineString segments that fall within a centre-line track-model rectangular polygon are discarded.
5. hexagons are classified as to whether they contain no rail, only OSM, only centre-line or a mixture of elements.
6. OSM data in hexagons only containing OSM elements are retained.
7. hexagons with mixed OSM and centre-line data are reassessed but at the next granular hexagon level.
8. the combined OSM only data and centre-line model then seems to give the required result.
