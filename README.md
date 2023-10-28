# Woodpecker
This is a way to overlay and combine two linear geographic network models using the rail centre-line track model and OpenStreetMap rail data as examples. The resulting network model is segmented using a rail waymark point layer and assigns approximate operation reporting points (STANOX) and planning timing point (TIPLOC) locations snapped to the network model segments

## Overview
The centre-line track-model is a dataset published via an Environmental Information Regulations (EIR) [request](https://www.whatdotheyknow.com/request/geospatial_data) by [Network Rail](https://www.networkrail.co.uk/who-we-are/transparency-and-ethics/freedom-of-information-foi) under the Open Government Licence:

Open Street Map contains current and historic mainland Britain railway data [here](https://www.openstreetmap.org/) published under the [Open Database License](https://opendatacommons.org/licenses/odbl/1-0/).

Run the following script:

    $ ./run.sh

This will downloads dependencies, centre-line and OSM data, creates the network overlay and combines this in a [outputx.gpkg](https://github.com/anisotropi4/woodpecker/blob/main/outputx.gpkg) GeoPKG repository. This combined network is then segmented using the waymark point layer and switch, end-point and locations in [segmentx.gpkg](https://github.com/anisotropi4/woodpecker/blob/main/segmentx.gpkg) and assign approximate operation reporting points (STANOX) and planning timing point (TIPLOC) locations snapped to OpenStreetMap and network model segments in [tiploc-location.gpkg](https://github.com/anisotropi4/woodpecker/blob/main/tiploc-location.gpkg)

### Overlay and combine approach
The approach used to overlay and combine data is:

1. create a cropped OSM and centre-line track-model LineString segments within scaled hexagons
2. create line-extended centre-line track-model rectangular polygons within a given hexagon
3. test if the OSM LineStrings are within these rectangular polygons
4. any OSM LineString segments that fall within a centre-line track-model rectangular polygon are discarded.
5. hexagons are classified as to whether they contain no rail, only OSM, only centre-line or a mixture of elements.
6. OSM data in hexagons only containing OSM elements are retained.
7. hexagons with mixed OSM and centre-line data are reassessed but at the next granular hexagon level.
8. the combined OSM only data and centre-line model then seems to give the required result.

### Segmentation approach
The approach used to segment the combined network is:

1. identify waymarks within a rectangle 12 track centre-to-centre separation of the network
2. split linear network features at the near nearest point the the waymark

The centre-to-centre track separation on the rail network is 3.26m

### Rail operation and planning location approach
The approach used to identify rail operation STANOX and planning TIPLOC locations is:

1. snap approximate STANOX locations to nearest OSM linear feature (`ox` layer in [tiploc-location.gpkg](https://github.com/anisotropi4/woodpecker/blob/main/tiploc-location.gpkg) 
2. snap approximate STANOX to combined network-model segement (`sx` layer in [tiploc-location.gpkg](https://github.com/anisotropi4/woodpecker/blob/main/tiploc-location.gpkg)
3. use the Network Rail [CORPUS](https://wiki.openraildata.com/index.php/Reference_Data#Downloading_CORPUS_Data) open rail reference data to link TIPLOC to STANOX locatation 

## Notes
