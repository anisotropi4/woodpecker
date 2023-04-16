#!/usr/bin/env python

import os
import datetime as dt
import pandas as pd

os.environ["USE_PYGEOS"] = "0"
import geopandas as gp
from pyogrio import read_dataframe, write_dataframe, list_layers
from shapely import line_locate_point, Point, LineString
from shapely.ops import nearest_points, split, snap

import momepy2 as mp2

START = dt.datetime.now()
pd.set_option("display.max_columns", None)
CRS = "EPSG:32630"
WGS84 = "EPSG:4326"
CENTRE2CENTRE = 3.26

INFILE = "linetrack.gpkg"

NETWORK = read_dataframe("data/network-model-simple.gpkg", layer="NetworkLinks")
NETWORK = NETWORK.dropna(how="all", axis=1)
OSMNX = read_dataframe("data/great-britain-rail-simple.gpkg", layer="lines")
OSMNX = OSMNX[OSMNX["location"] == "GB"].reset_index(drop=True)
OSMNX = OSMNX.drop_duplicates()

HEXAGON = read_dataframe("hexagon4.gpkg", layer="hexagon4-00")

WAYMARKS = read_dataframe("data/network-model-simple.gpkg", layer="WaymarksShape")

def get_buffer(gs, width=4, length=0):
    # style = {"cap_style": "square", "join_style": "mitre", "mitre_limit": length}
    style = {"cap_style": "flat"}
    gs = gs["geometry"].copy()
    return gs.buffer(width, **style)


def get_ds(this_array, name=None, axis=None):
    r = pd.Series(this_array.T[1], this_array.T[0])
    return r.rename_axis(axis).rename(name)


def get_overlap(gf1, gf2):
    data = []
    gf1 = gf1.copy()
    gf1["n"] = gf1.reset_index().index
    for k, v in gf1["geometry"].items():
        n = gf1.loc[k, "n"]
        if (n % 1024) == 0:
            print(f"{k}\t{n}\t{round(n / gf1.shape[0], 3):.3f}")
        ix = gf2["geometry"].within(v)
        s = pd.Series({i: k for i in ix[ix].index})
        data.append(pd.Series(s.index, index=s))
    return pd.concat(data).rename_axis(gf1.index.name).rename(gf2.index.name)


def get_npoint(gs1, gs2):
    r = [nearest_points(*i)[1] for i in zip(gs1, gs2)]
    return gp.GeoSeries(r, crs=CRS).rename("geometry")


def get_nearest(gs1, gs2):
    r = [LineString(nearest_points(*i)) for i in zip(gs1, gs2)]
    return gp.GeoSeries(r, crs=CRS).rename("geometry")


def get_split(gs1, gs2):
    r = [split(*i) for i in zip(gs1, gs2)]
    return gp.GeoSeries(r, crs=CRS).rename("geometry")


def get_extend(gs, length=1000.0):
    r = mp2.extend_lines(gs.to_frame(), 10.0, extension=length)
    return r.set_crs(CRS)


def get_offset(line, point):
    return line_locate_point(line, point)


def get_splits(v, separation=1.0e-4):
    line, point = v["line"], v["geometry"]
    return list(split(snap(line, point, separation), point).geoms)


def get_end_point(v):
    return Point(v.reverse().coords[0])


def get_end_segment(this_gf):
    r = this_gf.set_index("ASSETID").drop_duplicates(subset="line")
    r["offset"] = r["line"].length
    r["M_POST_ID"] = 0
    r["geometry"] = r["line"].apply(get_end_point)
    return r.reset_index()


def get_gs(this_array, this_index):
    return gp.GeoSeries(this_array, index=this_index).dropna()


def get_segments(this_network):
    data = []
    gf = this_network.drop_duplicates("ASSETID")
    ix = gf.index
    end_ix = this_network[this_network["M_POST_ID"] == 0].index
    s, rest = gf.apply(get_splits, axis=1).apply(pd.Series).to_numpy().T
    data.append(gp.GeoSeries(s, index=ix, name="segment"))

    while not (gs := get_gs(rest, (ix + 1))).empty:
        print(f'{gs.shape[0]}\t{gf.shape[0]}')
        ix = gs.index.intersection(end_ix)
        if not ix.empty:
            data.append(gs.loc[ix])
        ix = gs.index.difference(end_ix)
        if ix.empty:
            break
        gf = this_network.loc[ix]
        gf["line"] = gs
        s, rest = gf.apply(get_splits, axis=1).apply(pd.Series).to_numpy().T
        data.append(gp.GeoSeries(s, index=ix, name="segment"))

    return gp.GeoSeries(pd.concat(data)).sort_index()


def get_basenx():
    data = []
    layers = {k for k, _ in list_layers(INFILE)}
    for n in range(16):
        layer = f"OSM{str(n).zfill(2)}"
        if layer in layers:
            print(f"layer: {n}", dt.datetime.now() - START)
            osmx = read_dataframe(INFILE, layer=layer)
            idx = osmx[osmx["m"] == 0].index
            data.append(osmx.loc[idx])
    basenx = pd.concat(data).reset_index(drop=True)
    basenx["length"] = basenx.length
    basenx = basenx.sjoin_nearest(OSMNX, max_distance=1.0, distance_col="d")
    keys = ["osmx2", "hex_id", "OSM"]
    basenx = basenx.sort_values("d").drop_duplicates(subset=keys)
    return basenx


def write_basenx(basenx, outfile, layer):
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
    write_dataframe(basenx[fields], outfile, layer=layer)


def set_fullnx(basenx, outfile, layer):
    fullnx = basenx[["osmid", "geometry"]].copy()
    fullnx["length"] = fullnx.length
    fullnx = fullnx.sort_values("length", ascending=False)
    fullnx = fullnx.drop_duplicates(subset="osmid")
    ix = pd.Index(fullnx["osmid"])
    fullnx = OSMNX.set_index("osmid").loc[ix].reset_index()
    write_dataframe(fullnx, outfile, layer=layer)


def get_network(osmnx):
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
    network["ELD"] = network["ELR"].str[:3]
    return network


def write_network(network, outfile, layer):
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
        "ELD",
    ]
    write_dataframe(network[fields], outfile, layer=layer)


def set_elrnx(osmnx):
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
        "ELD",
    ]

    elrnx = osmnx.sjoin_nearest(
        NETWORK, max_distance=2 * CENTRE2CENTRE, distance_col="d"
    )
    elrnx = (
        elrnx.sort_values("d").drop_duplicates(subset="osmid").reset_index(drop=True)
    )
    elrnx["elr_match"] = "differ"
    elrnx.loc[elrnx["ref"].isna(), "elr_match"] = "no osm"
    elrnx.loc[elrnx["ref"] == elrnx["ELR"], "elr_match"] = "match"
    elrnx["ELD"] = elrnx["ELR"].str[:3]
    if "elr_match" not in fields:
        fields += ["elr_match"]
    write_dataframe(elrnx[fields], OUTFILE, layer="ELR")


def combine_network(waymarks, network, width=CENTRE2CENTRE):
    gf1 = waymarks[["M_POST_ID", "ELR", "geometry"]]
    gf2 = network[["ASSETID", "ELR", "geometry"]]

    gf = get_buffer(gf2.set_index("ASSETID"), width=width)
    gf = gf.to_frame("geometry")
    ix = get_overlap(gf, gf1.set_index("M_POST_ID"))
    ix = ix.reset_index().drop_duplicates()
    ix = ix.sort_values("ASSETID").reset_index(drop=True)

    gs1 = gf1.set_index("M_POST_ID").loc[ix["M_POST_ID"]]
    gs2 = gf2.set_index("ASSETID").loc[ix["ASSETID"]]
    r = get_npoint(gs1["geometry"], gs2["geometry"]).to_frame()
    r["M_POST_ID"] = gs1.index
    r.index = gs2.index
    r["line"] = gf2.set_index("ASSETID").loc[r.index, "geometry"]
    r["offset"] = get_offset(r["line"], r["geometry"])
    return r.sort_values(["ASSETID", "offset"]).reset_index()


def get_segmented_nx(network):
    gf1 = combine_network(WAYMARKS, network, 10 * CENTRE2CENTRE)
    gf2 = get_end_segment(gf1)
    r = pd.concat([gf1, gf2]).drop_duplicates(subset=["ASSETID", "geometry"])
    r["segment"] = LineString([])
    r = r.sort_values(["ASSETID", "offset"]).reset_index(drop=True)

    fields = ["ASSETID", "M_POST_ID", "offset", "geometry"]
    post = r[fields]

    fields = ["ASSETID", "M_POST_ID", "geometry", "line"]
    gs = get_segments(r[fields])
    r.loc[gs.index, "segment"] = gs

    fields = ["ASSETID", "M_POST_ID", "offset", "segment"]
    segment = gp.GeoDataFrame(r[fields].rename(columns={"segment": "geometry"}))
    segment["length"] = segment.length
    segment = segment.set_crs(CRS)

    nx = network.set_index("ASSETID")
    fields = [
        "ELR",
        "TRID",
        "TRCODE",
        "L_M_FROM",
        "L_M_TO",
        "L_SYSTEM",
        "L_LINK_ID",
        "SHAPE_LEN",
    ]
    ix = segment.set_index("ASSETID").index
    segment[fields] = nx[fields].loc[ix].values

    ix = nx.index.difference(segment['ASSETID'])
    remainder = nx.loc[ix, fields + ["geometry"]].reset_index()
    remainder[["M_POST_ID", "offset"]] = [0, 0.0]
    remainder["length"] = remainder.length

    segment = pd.concat([segment, remainder])
    segment = segment.drop_duplicates(["ASSETID", "geometry"])
    segment = segment.sort_values(["ASSETID", "offset"]).reset_index(drop=True)

    return (post, segment)


def main():
    outfile = "outputx.gpkg"
    write_dataframe(HEXAGON, outfile, layer="hex")

    basenx = get_basenx()
    write_basenx(basenx, outfile, "basenx")
    set_fullnx(basenx, outfile, "fullnx")

    ix = OSMNX["railway"].isin(["rail", "light_rail", "tram"])
    osmnx = OSMNX[ix].reset_index(drop=True)
    network = get_network(osmnx)
    write_network(network, outfile, layer="network")

    outfile = "segmentx.gpkg"
    post, segment = get_segmented_nx(network)
    write_dataframe(post, outfile, layer="post")
    write_dataframe(segment, outfile, layer="segment")


if __name__ == "__main__":
    main()
