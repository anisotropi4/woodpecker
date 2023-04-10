#!/usr/bin/env python

import os
from functools import partial
import datetime as dt
import pandas as pd

os.environ["USE_PYGEOS"] = "0"
import geopandas as gp
from pyogrio import read_dataframe, write_dataframe, list_layers

START = dt.datetime.now()
pd.set_option("display.max_columns", None)
CRS = "EPSG:32630"
WGS84 = "EPSG:4326"
CENTRE2CENTRE = 3.26

INFILE = "linetrack.gpkg"
OUTFILE = "outputx.gpkg"

NETWORK = read_dataframe("data/network-model-simple.gpkg", layer="NetworkLinks")
NETWORK = NETWORK.dropna(how="all", axis=1)
OSMNX = read_dataframe("data/great-britain-rail-simple.gpkg", layer="lines")
OSMNX = OSMNX[OSMNX["location"] == "GB"].reset_index(drop=True)
write_dataframe(OSMNX, OUTFILE, layer="osmnx")

HXOSM = read_dataframe(INFILE, layer="osmnx")
HEXAGON = read_dataframe("hexagon4.gpkg", layer="hexagon4-00")
write_dataframe(HEXAGON, OUTFILE, layer='hex')


def main():
    ix = OSMNX["railway"].isin(["rail", "light_rail", "tram"])
    osmnx = OSMNX[ix].reset_index(drop=True)
    data = []
    layers = {k for k, _ in list_layers(INFILE)}
    ix = pd.Index([])
    for n in range(16):
        layer = f"OSM{str(n).zfill(2)}"
        if layer in layers:
            print(f"layer: {n}", dt.datetime.now() - START)
            osmx = read_dataframe(INFILE, layer=layer)
            idx = osmx[osmx["m"] == 0].index
            data.append(osmx.loc[idx])
    fields = [
        "osmid",
        "maxspeed",
        "name",
        "ref",
        "electrified",
        "railway",
        "tunnel",
        "bridge",
        "oneway",
        "frequency",
        "voltage",
        "ref:tiploc",
        "landuse",
        "type",
        "location",
        "geometry",
        "hex_id",
    ]
    basenx = pd.concat(data).reset_index(drop=True)
    basenx["length"] = basenx.length
    basenx = basenx.sjoin_nearest(OSMNX, max_distance=1.0, distance_col="d")
    keys = ["osmx2", "hex_id", "OSM"]
    basenx = basenx.sort_values("d").drop_duplicates(subset=keys)
    write_dataframe(basenx[fields], OUTFILE, layer="basex")
    fullnx = basenx[["osmid", "geometry"]].copy()
    fullnx["length"] = fullnx.length
    fullnx = fullnx.sort_values("length", ascending=False)
    fullnx = fullnx.drop_duplicates(subset="osmid")
    ix = pd.Index(fullnx["osmid"])
    fullnx = OSMNX.set_index("osmid").loc[ix].reset_index()
    write_dataframe(fullnx, OUTFILE, layer="fullnx")
    overlap = NETWORK.sjoin_nearest(
        osmnx, max_distance=2 * CENTRE2CENTRE, distance_col="d"
    )
    overlap = overlap.sort_values("d").drop_duplicates(subset="ASSETID")
    overlap = overlap.reset_index(drop=True)
    idx = pd.Index(NETWORK["ASSETID"]).difference(pd.Index(overlap["ASSETID"]))
    clipped = NETWORK.set_index("ASSETID").loc[idx].reset_index()
    clipped[["railway", "location"]] = ["rail", "GB"]
    network = pd.concat([overlap, clipped]).fillna("-").reset_index(drop=True)
    network = network.sort_values("ASSETID")
    fields = [
        "OBJECTID",
        "ASSETID",
        "L_LINK_ID",
        "L_SYSTEM",
        "L_VAL",
        "L_QUALITY",
        "ELR",
        "TRID",
        "TRCODE",
        "L_M_FROM",
        "L_M_TO",
        "TRACK_STAT",
        "SHAPE_Leng",
        "geometry",
        "osmid",
        "maxspeed",
        "name",
        "ref",
        "electrified",
        "railway",
        "tunnel",
        "bridge",
        "oneway",
        "frequency",
        "voltage",
        "ref:tiploc",
        "landuse",
        "type",
        "location",
    ]
    write_dataframe(network[fields], OUTFILE, layer="network")
    # IDX = NETWORK[['geometry']].sjoin(GDATA[['OSM', 'geometry']], how='inner')
    # IDX = pd.Index(IDX['OSM']).drop_duplicates()

    # GDATA = HXOSM.set_index('OSM').loc[OKEY].reset_index(names='OSM')
    # write_dataframe(GDATA, 'xy.gpkg', layer='hxosm')


if __name__ == "__main__":
    main()
