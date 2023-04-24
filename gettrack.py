#!/usr/bin/env python

import os
import warnings
import datetime as dt
import pandas as pd

os.environ["USE_PYGEOS"] = "0"
import geopandas as gp
from pyogrio import read_dataframe, write_dataframe, list_layers
from shapely import line_locate_point, line_interpolate_point, Point, LineString
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
NODE = read_dataframe("data/network-model-simple.gpkg", layer="VectorNodes")
NODE = NODE.dropna(how="all", axis=1)

OSMNX = read_dataframe("data/great-britain-rail-simple.gpkg", layer="lines")
OSMNX = OSMNX[OSMNX["location"] == "GB"].reset_index(drop=True)
OSMNX = OSMNX.drop_duplicates()

HEXAGON = read_dataframe("hexagon4.gpkg", layer="hexagon4-00")

WAYMARK = read_dataframe("data/network-model-simple.gpkg", layer="Waymarks")


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
        if (n % 4096) == 0:
            print(f"{k}\t{n}\t{round(n / gf1.shape[0], 3):.3f}")
        ix = gf2["geometry"].within(v)
        s = pd.Series({i: k for i in ix[ix].index})
        data.append(pd.Series(s.index, index=s))
    return pd.concat(data).rename_axis(gf1.index.name).rename(gf2.index.name)


def get_npoint(gs1, gs2):
    r = [nearest_points(*i)[1] for i in zip(gs1, gs2)]
    return gp.GeoSeries(r, crs=CRS).rename("geometry")


def get_offset(line, point):
    return line_locate_point(line, point)


def get_nearest(gs1, gs2):
    r = [LineString(nearest_points(*i)) for i in zip(gs1, gs2)]
    return gp.GeoSeries(r, crs=CRS).rename("geometry")


def get_split(gs1, gs2):
    r = [split(*i) for i in zip(gs1, gs2)]
    return gp.GeoSeries(r, crs=CRS).rename("geometry")


def get_extend(gs, length=1000.0):
    r = mp2.extend_lines(gs.to_frame(), 10.0, extension=length)
    return r.set_crs(CRS)


def get_splits(v, separation=1.0e-4):
    line, point = v["line"], v["geometry"]
    return list(split(snap(line, point, separation), point).geoms)


def get_extrema(v):
    return (Point(v.coords[0]), Point(v.coords[-1]))


def get_start_point(v):
    return Point(v.coords[0])


def get_end_point(v):
    return Point(v.coords[-1])


def get_end_segment(this_gf):
    r = this_gf.set_index("ASSETID").drop_duplicates(subset="line")
    r["offset"] = r["line"].length
    r["M_POST_ID"] = 0
    r["geometry"] = r["line"].apply(get_end_point)
    return r.reset_index()


def get_gs(this_array, this_index):
    return gp.GeoSeries(this_array, index=this_index).dropna()


def get_sections(this_gf, separation=1.0e-4):
    line, offset = this_gf["line"], this_gf["offset"]
    point = line_interpolate_point(line, offset).set_crs(CRS).rename("geometry")
    r = pd.concat([line, point], axis=1).apply(get_splits, axis=1)
    r = r.apply(pd.Series).rename(columns={0: "segment"})
    return r["segment"]


def get_segments(this_network):
    warnings.warn("ignore", RuntimeWarning)
    data = []
    gf = this_network.drop_duplicates("ASSETID")
    ix = gf.index
    end_ix = this_network[this_network["M_POST_ID"] == 0].index
    s, rest = gf.apply(get_splits, axis=1).apply(pd.Series).to_numpy().T
    data.append(gp.GeoSeries(s, index=ix, name="segment"))
    n = gf.shape[0]
    while not (gs := get_gs(rest, (ix + 1))).empty:
        if gs.shape[0] > 128:
            percent = round(n / this_network.shape[0], 3)
            print(f"{n}\t{this_network.shape[0]}\t{percent}")
        n += gs.shape[0]
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
    percent = round(n / this_network.shape[0], 3)
    print(f"{n}\t{this_network.shape[0]}\t{percent}")
    warnings.warn("default", RuntimeWarning)
    return gp.GeoSeries(pd.concat(data)).sort_index()


def get_basenx(osmnx):
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
    basenx = basenx.sjoin_nearest(osmnx, max_distance=1.0, distance_col="d")
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


def set_fullnx(basenx, osmnx, outfile, layer):
    fullnx = basenx.copy()
    fullnx["length"] = fullnx.length
    fullnx = basenx[["osmid", "length"]].sort_values("length", ascending=False)
    ix = pd.Index(fullnx["osmid"].drop_duplicates())
    fullnx = osmnx.set_index("osmid").loc[ix].reset_index()
    write_dataframe(fullnx, outfile, layer=layer)


def get_network(network, osmnx, distance=CENTRE2CENTRE):
    overlap = network.sjoin_nearest(osmnx, max_distance=distance, distance_col="d")
    overlap = overlap.sort_values("d").drop_duplicates(subset="ASSETID")
    overlap = overlap.reset_index(drop=True)
    idx = pd.Index(network["ASSETID"]).difference(pd.Index(overlap["ASSETID"]))
    clipped = network.set_index("ASSETID").loc[idx].reset_index()
    clipped[["railway", "location"]] = ["rail", "GB"]
    network = pd.concat([overlap, clipped]).fillna("-").reset_index(drop=True)
    network = network.sort_values(["ASSETID", "L_M_FROM"])
    network["ELD"] = network["ELR"].str[:3]
    network = network.rename(columns={"SHAPE_Leng": "SHAPE_LEN"})
    return network.reset_index(drop=True)


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


def set_elrnx(osmnx, network, outfile, distance=2 * CENTRE2CENTRE):
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
    elrnx = osmnx.sjoin_nearest(network, max_distance=distance, distance_col="d")
    elrnx = (
        elrnx.sort_values("d").drop_duplicates(subset="osmid").reset_index(drop=True)
    )
    elrnx["elr_match"] = "differ"
    elrnx.loc[elrnx["ref"].isna(), "elr_match"] = "no osm"
    elrnx.loc[elrnx["ref"] == elrnx["ELR"], "elr_match"] = "match"
    elrnx["ELD"] = elrnx["ELR"].str[:3]
    if "elr_match" not in fields:
        fields += ["elr_match"]
    write_dataframe(elrnx[fields], outfile, layer="ELR")


def combine_network(waymarks, network, width=CENTRE2CENTRE):
    gf1 = waymarks[["M_POST_ID", "geometry"]]
    gf2 = network[["ASSETID", "geometry"]]

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
    gf2 = gf2.set_index("ASSETID")
    r["line"] = gf2.loc[r.index, "geometry"]
    warnings.warn("ignore", RuntimeWarning)
    r["offset"] = get_offset(r["line"], r["geometry"])
    warnings.warn("default", RuntimeWarning)
    return r.sort_values(["ASSETID", "offset"]).reset_index()


def get_post(this_gf, node):
    fields = ["ASSETID", "M_POST_ID", "ELR", "offset", "geometry"]
    post = this_gf[fields]
    fields = ["M_POST_ID", "M_SYSTEM", "WAYMARK_VA"]
    post = post.join(waymark[fields].set_index("M_POST_ID"), on="M_POST_ID")
    post["M_SYSTEM"] = post["M_SYSTEM"].fillna("-")
    post = post.fillna(0.0)
    post["id"] = post["M_POST_ID"]
    post["waymark"] = post["id"] > 0
    return post


def get_waynode(waymark, node):
    fields = ["ASSETID", "geometry"]
    r = pd.concat([waymark[fields], node[fields]]).set_index("ASSETID")
    r["waymark"] = False
    ix = r.index.difference(node["ASSETID"])
    r.loc[ix, "waymark"] = True
    return r.sort_index()


def get_point_key(this_gf, key):
    this_gf = this_gf.set_index(key)
    r = pd.concat([this_gf["source"], this_gf["target"]]).rename("point_id")
    r = r.reset_index().drop_duplicates()
    r = r.drop_duplicates(subset="point_id", keep=False)
    return r.set_index("point_id").sort_index()


def get_point_count(this_gf):
    r = pd.concat([this_gf["source"], this_gf["target"]]).to_frame("point_id")
    r["#"] = 1
    return r.groupby("point_id").count()


def get_segmented_nx(this_nx, this_waymark, this_node):
    gf1 = combine_network(this_waymark, this_nx, 12 * CENTRE2CENTRE)
    gf2 = get_end_segment(gf1)
    r = pd.concat([gf1, gf2]).drop_duplicates(subset=["ASSETID", "geometry"])
    r["segment"] = LineString([])
    r = r.sort_values(["ASSETID", "offset"]).reset_index(drop=True)

    fields = ["ASSETID", "M_POST_ID", "geometry", "line"]
    segments = get_segments(r[fields])
    r.loc[segments.index, "segment"] = segments
    ix = r.index.difference(segments.index)
    sections = get_sections(r.loc[ix, ["ASSETID", "line", "offset"]])
    r.loc[sections.index, "segment"] = sections

    nx = this_nx.set_index("ASSETID")
    ix = nx.index.difference(r["ASSETID"])
    remainder = nx.loc[ix, "geometry"].reset_index()
    remainder[["M_POST_ID", "offset"]] = [0, 0.0]
    for c in ["segment", "line"]:
        remainder[c] = remainder["geometry"]
    remainder["geometry"] = remainder["segment"].apply(get_end_point)
    r = pd.concat([r, remainder]).reset_index(drop=True)
    r = r.sort_values(["ASSETID", "offset"]).reset_index(drop=True)
    r["segment_id"] = "S" + (1 + r.index).astype(str).str.zfill(7)

    post = r[["M_POST_ID", "geometry"]].drop_duplicates().reset_index(drop=True)
    wx = this_waymark.set_index("M_POST_ID")
    post = post.join(wx["ASSETID"], on="M_POST_ID").dropna()
    nx = this_node[["ASSETID", "geometry"]].copy()
    nx["M_POST_ID"] = 0
    post = pd.concat([post, nx]).sort_values("ASSETID").reset_index(drop=True)
    post = post.reset_index(drop=True)
    post["point_id"] = "P" + (1 + post.index).astype(str).str.zfill(7)

    rx = r.set_index("segment_id")
    px = post[["geometry", "point_id"]]
    gs = rx["geometry"]
    target = gs.to_frame("geometry").sjoin_nearest(px, distance_col="d")
    rx["target"] = target["point_id"]
    gs = r[["segment_id", "segment"]].set_index("segment_id")
    gs = gs["segment"].apply(get_start_point).rename("geometry")
    gs = gp.GeoSeries(gs, crs=CRS)
    source = gs.to_frame("geometry").sjoin_nearest(px, distance_col="d")
    rx["source"] = source["point_id"]
    fields = [
        "source",
        "target",
        "ASSETID",
        "M_POST_ID",
        "offset",
        "segment",
    ]
    r = rx[fields].reset_index().rename(columns={"segment": "geometry"})

    segment = gp.GeoDataFrame(r, crs=CRS)
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
    ix = r["ASSETID"]
    segment[fields] = this_nx.set_index("ASSETID").loc[ix, fields].values
    segment = segment.drop_duplicates(["ASSETID", "geometry"])
    segment = segment.sort_values(["ASSETID", "offset"]).reset_index(drop=True)

    nx = this_node.set_index("ASSETID")
    post["CONNECTED_"] = 2
    ix = post.loc[post["M_POST_ID"] == 0, "ASSETID"]
    post.loc[ix.index, "CONNECTED_"] = nx.loc[ix, "CONNECTED_"].values
    fields = [
        "M_SYSTEM",
        "ELR",
        "WAYMARK_VA",
    ]
    post[fields] = ["U", "-", 0.0]

    wx = this_waymark.set_index("ASSETID")
    px = post.set_index("point_id")
    ix = pd.Index(post["ASSETID"]).difference(ix)
    ix = post.set_index("ASSETID").loc[ix, "point_id"]
    px.loc[ix, fields] = wx.loc[ix.index, fields].values
    gs = get_point_key(segment, "ELR")
    px.loc[gs.index, "ELR"] = gs
    gs = get_point_key(segment, "L_SYSTEM")
    px.loc[gs.index, "M_SYSTEM"] = gs.values

    ix = px[px["M_POST_ID"] > 0].index
    gs = get_point_count(segment)
    px.loc[ix, "CONNECTED_"] = gs.loc[ix].values
    post = px.reset_index()

    return (post, segment)


def main():
    outfile = "outputx.gpkg"
    write_dataframe(HEXAGON, outfile, layer="hex")

    basenx = get_basenx(OSMNX)
    write_basenx(basenx, outfile, "basenx")
    set_fullnx(basenx, OSMNX, outfile, "fullnx")

    ix = OSMNX["railway"].isin(["rail", "light_rail", "tram"])
    osmnx = OSMNX[ix].reset_index(drop=True)
    write_dataframe(osmnx, outfile, layer="osmnx")

    network = get_network(NETWORK, osmnx, 2 * CENTRE2CENTRE)
    write_network(network, outfile, layer="network")

    outfile = "segmentx.gpkg"
    post, segment = get_segmented_nx(network, WAYMARK, NODE)
    write_dataframe(post, outfile, layer="post")
    write_dataframe(segment, outfile, layer="segment")
    print(f"segemented", dt.datetime.now() - START)


if __name__ == "__main__":
    main()
