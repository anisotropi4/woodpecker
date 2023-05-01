#!/usr/bin/env python

import os
import datetime as dt
import pandas as pd
from functools import partial

os.environ["USE_PYGEOS"] = "0"
import geopandas as gp
from pyogrio import read_dataframe, write_dataframe, list_layers
from shapely import line_locate_point, line_interpolate_point, Point, LineString
from shapely.ops import nearest_points, split, snap

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


def get_nearest_point(gs1, gs2):
    r = [nearest_points(*i)[0] for i in zip(gs1, gs2)]
    return gp.GeoSeries(r, crs=CRS).rename("geometry")


def get_offset(line, point):
    return line_locate_point(line, point)


def get_split(v, separation=1.0e-4):
    line, point = v["line"], v["geometry"]
    return list(split(snap(line, point, separation), point).geoms)


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


def get_section(this_gf, ix, separation=1.0e-4):
    get_splits = partial(get_split, separation=separation)
    this_gf = this_gf.loc[this_gf.index.difference(ix)]
    line, offset = this_gf["line"].rename("line"), this_gf["offset"].values
    point = line_interpolate_point(line, offset).rename("geometry")
    section = pd.concat([line, point], axis=1).apply(get_splits, axis=1)
    section = section.apply(pd.Series).rename(columns={0: "segment"})
    return section["segment"]


def get_segment(this_network, separation=1.0e-4):
    get_splits = partial(get_split, separation=separation)
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


def overlay_nx_waymark(
    network, waymarks, width=CENTRE2CENTRE, nxkey="ASSETID", wxkey="M_POST_ID"
):
    """overlay_nx_waymark: create a rectangular polygon buffer of width on network line elements, identify waymarks within the buffer. Return a GeoDataFrame of all line waymark"""
    gf1 = network[[nxkey, "geometry"]]
    gf2 = waymarks[[wxkey, "geometry"]]

    this_buffer = get_buffer(gf1.set_index(nxkey), width=width)
    this_buffer = this_buffer.to_frame("geometry")
    ix = get_overlap(this_buffer, gf2.set_index(wxkey))
    ix = ix.reset_index().drop_duplicates()
    ix = ix.sort_values(nxkey).reset_index(drop=True)

    gs1 = gf1.set_index(nxkey).loc[ix[nxkey]]
    gs2 = gf2.set_index("M_POST_ID").loc[ix[wxkey]]
    r = get_nearest_point(gs1["geometry"], gs2["geometry"]).to_frame()
    r.index = gs1.index
    r[wxkey] = gs2.index
    gf1 = gf1.set_index(nxkey)
    r["line"] = gf1.loc[r.index, "geometry"]
    r["offset"] = get_offset(r["line"], r["geometry"])
    return r.sort_values([nxkey, "offset"]).reset_index()


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


def get_point_m_value(this_segment, ix):
    source = this_segment.set_index("source")["L_M_FROM"]
    source = source.loc[ix.intersection(source.index)]
    target = this_segment.set_index("target")["L_M_TO"]
    target = target.loc[ix.intersection(target.index)]
    r = pd.concat([source, target]).rename("WAYMARK_VA").sort_index()
    r = r.reset_index().drop_duplicates().rename(columns={"index": "point_id"})
    r = r.drop_duplicates(subset="point_id", keep=False)
    return r.set_index("point_id")


def get_base_sx(this_nx, this_waymark):
    overlay_nx = overlay_nx_waymark(this_nx, this_waymark, 12 * CENTRE2CENTRE)
    end_nx = get_end_segment(overlay_nx)
    r = pd.concat([overlay_nx, end_nx])
    r = r.drop_duplicates(subset=["ASSETID", "geometry"])
    r["segment"] = LineString([])
    return r.sort_values(["ASSETID", "offset"]).reset_index(drop=True)


def get_remainder_sx(this_nx, this_sx):
    nx = this_nx.set_index("ASSETID")
    ix = nx.index.difference(this_sx["ASSETID"])
    remainder = nx.loc[ix, "geometry"].reset_index()
    remainder[["M_POST_ID", "offset"]] = [0, 0.0]
    for i in ["segment", "line"]:
        remainder[i] = remainder["geometry"]
    remainder["geometry"] = remainder["segment"].apply(get_end_point)
    return remainder


def get_full_sx(this_nx, this_waymark):
    sx = get_base_sx(this_nx, this_waymark)
    fields = ["ASSETID", "M_POST_ID", "geometry", "line"]
    gs = get_segment(sx[fields])
    sx.loc[gs.index, "segment"] = gs
    section = get_section(sx[["ASSETID", "line", "offset"]], gs.index)
    sx.loc[section.index, "segment"] = section

    remainder = get_remainder_sx(this_nx, sx)
    sx = pd.concat([sx, remainder]).reset_index(drop=True)
    sx = sx.sort_values(["ASSETID", "offset"]).reset_index(drop=True)

    fields = [
        "L_LINK_ID",
        "L_SYSTEM",
        "ELR",
        "TRID",
        "TRCODE",
        "L_M_FROM",
        "L_M_TO",
        "SHAPE_LEN",
    ]
    nx = this_nx.set_index("ASSETID")
    sx[fields] = nx.loc[sx["ASSETID"], fields].values
    sx["segment_id"] = "S" + (1 + sx.index).astype(str).str.zfill(7)
    return sx


def get_full_px(this_sx, this_waymark, this_node):
    px = this_sx[["M_POST_ID", "geometry"]].copy().drop_duplicates()
    px = px.set_index("M_POST_ID")
    ix = px.index[px.index > 0]
    wx = this_waymark.sort_values("WAYMARK_VA").set_index("M_POST_ID")
    fields = ["ASSETID", "M_SYSTEM", "ELR", "WAYMARK_VA"]
    px = px.loc[ix].join(wx[fields]).reset_index()
    nx = this_node[["ASSETID", "geometry"]].copy()
    fields = [
        "M_POST_ID",
        "M_SYSTEM",
        "ELR",
        "WAYMARK_VA",
    ]
    nx[fields] = [0, "U", "-", 0.0]
    px = pd.concat([px, nx]).sort_values("ASSETID").reset_index(drop=True)
    px["CONNECTED"] = 2
    px = px.sort_values("WAYMARK_VA").drop_duplicates(subset="geometry")
    px = px.sort_values(["M_POST_ID", "ASSETID"]).reset_index(drop=True)
    ix = px.loc[px["M_POST_ID"] == 0, "ASSETID"]
    cx = this_node[["ASSETID", "CONNECTED_"]].set_index("ASSETID")
    px.loc[ix.index, "CONNECTED"] = cx.loc[ix, "CONNECTED_"].values
    px["point_id"] = "P" + (1 + px.index).astype(str).str.zfill(7)
    return px


def match_segment_point(this_segment, this_post):
    sx = this_segment.set_index("segment_id")
    px = this_post[["geometry", "point_id", "ASSETID"]]
    gs = sx["segment"].apply(get_start_point).rename("geometry")
    gs = gp.GeoSeries(gs, crs=CRS)
    source = gs.to_frame("geometry").sjoin_nearest(px, distance_col="d")
    sx[["source", "source_asset"]] = source[["point_id", "ASSETID"]]
    gs = sx["geometry"]
    target = gs.to_frame("geometry").sjoin_nearest(px, distance_col="d")
    target = target.sort_index()
    sx[["target", "target_asset"]] = target[["point_id", "ASSETID"]]
    r = sx.reset_index().drop(columns=["geometry", "line"])
    r = r.rename(columns={"segment": "geometry"})
    return gp.GeoDataFrame(r, crs=CRS)


def match_point_segment(this_post, this_segment):
    px = this_post.set_index("point_id")
    gs = get_point_key(this_segment, "ELR")
    px["ELR2"] = "-"
    px.loc[gs.index, "ELR2"] = gs.values

    gs = get_point_key(this_segment, "L_SYSTEM")
    px["L_SYSTEM"] = "-"
    px.loc[gs.index, "L_SYSTEM"] = gs.values

    gs = get_point_m_value(this_segment, px[px["M_POST_ID"] == 0].index)
    px.loc[gs.index, "WAYMARK_VA"] = gs
    px.loc[gs.index, "M_SYSTEM"] = px.loc[gs.index, "L_SYSTEM"]
    return px.reset_index()


def get_segmented_nx(this_nx, this_waymark, this_node):
    segment = get_full_sx(this_nx, this_waymark)
    post = get_full_px(segment, this_waymark, this_node)
    segment = match_segment_point(segment, post)
    post = match_point_segment(post, segment)
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
    # this_nx, this_waymark, this_node = network, WAYMARK, NODE
    post, segment = get_segmented_nx(network, WAYMARK, NODE)
    write_dataframe(post, outfile, layer="post")
    write_dataframe(segment, outfile, layer="segment")
    print(f"segemented {dt.datetime.now() - START}")


if __name__ == "__main__":
    main()
