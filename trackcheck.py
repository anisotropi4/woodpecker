#!/usr/bin/env python3

import os
from functools import partial
import datetime as dt
import pandas as pd

os.environ["USE_PYGEOS"] = "0"
import geopandas as gp
import numpy as np


from shapely.geometry import MultiLineString, Polygon
from shapely import hausdorff_distance
from shapely.ops import linemerge
from shapely import oriented_envelope
import h3

from pyogrio import read_dataframe, write_dataframe
from tobler.util import h3fy

import networkx as nx
import momepy2 as mp2

# from mapclassify import greedy

pd.set_option("display.max_columns", None)
CRS = "EPSG:32630"
WGS84 = "EPSG:4326"
CENTRE2CENTRE = 3.26

START = dt.datetime.now()


def merge_network(this_gs):
    r = gp.GeoSeries(linemerge(MultiLineString(this_gs.values)))
    r = gp.GeoDataFrame(geometry=r.explode(ignore_index=True), crs=CRS)
    return r.drop_duplicates().reset_index(drop=True)


def get_hexagon(hex_id):
    return Polygon(h3.h3_to_geo_boundary(hex_id, geo_json=True))


def get_hexagons(hex_df):
    r = gp.GeoDataFrame(geometry=hex_df.apply(get_hexagon), crs=WGS84)
    r.name = "geometry"
    r.index = hex_df.values
    return r.reset_index(names="hex_id").to_crs(CRS)


def get_subhexagon(gf, resolution):
    gethex = partial(h3fy, resolution=resolution, return_geoms=False)
    df = gf[["geometry"]].apply(gethex)
    r = get_hexagons(df["geometry"])
    k = f"idx{str(resolution - 1).zfill(2)}"
    r[k] = r.index // 7
    return r


OUTFILE = "linetrack.gpkg"
NETWORK = read_dataframe("data/network-model-simple.gpkg", layer="NetworkLinks")
OSMNX = read_dataframe("data/great-britain-rail-simple.gpkg", layer="lines")
OSMNX = OSMNX[OSMNX["location"] == "GB"].reset_index(drop=True)
OSMNX = OSMNX.drop_duplicates()

print(1, dt.datetime.now() - START)

if 'HXNX' not in globals():
    HXNX = mp2.remove_false_nodes(NETWORK["geometry"]).set_crs(CRS)
    HXNX = HXNX.to_frame(name="geometry")
    HXNX = merge_network(NETWORK["geometry"])
    HXNX["NX"] = HXNX.index
    # HXNX['colour'] = greedy(HXNX)
    HXNX["class"] = -1
    NXNX = mp2.gdf_to_nx(HXNX)
    for n, G in enumerate(nx.connected_components(NXNX)):
        _, edges = mp2.nx_to_gdf(NXNX.subgraph(G).copy())
        HXNX.loc[edges.set_index("NX").index, "class"] = n

write_dataframe(HXNX, OUTFILE, layer="linenx")

print(2, dt.datetime.now() - START)
if 'HXOSM' not in globals():
    IDX = OSMNX[OSMNX["railway"] != "rail"].index
    HXNOSM = OSMNX.loc[IDX, "geometry"].drop_duplicates()
    NXNOSM = mp2.remove_false_nodes(HXNOSM).set_crs(CRS)
    HXNOSM = HXNOSM.drop_duplicates().reset_index(drop=True)
    HXNOSM = HXNOSM.to_frame(name="geometry")
    # HXOSM = merge_network(HXOSM['geometry'])
    HXNOSM["OSN"] = HXNOSM.index
    write_dataframe(HXNOSM, OUTFILE, layer="nosmx")

    IDX = OSMNX[OSMNX["railway"] == "rail"].index
    HXOSM = OSMNX.loc[IDX, "geometry"].drop_duplicates()
    HXOSM = mp2.remove_false_nodes(HXOSM).set_crs(CRS)
    HXOSM = HXOSM.drop_duplicates().reset_index(drop=True)
    HXOSM = HXOSM.to_frame(name="geometry")
    # HXOSM = merge_network(HXOSM['geometry'])
    HXOSM["OSM"] = HXOSM.index
    # HXOSM['colour'] = greedy(HXOSM)
    HXOSM["class"] = -1
    NXOSM = mp2.gdf_to_nx(HXOSM)
    for n, G in enumerate(nx.connected_components(NXOSM)):
        _, edges = mp2.nx_to_gdf(NXOSM.subgraph(G).copy())
        HXOSM.loc[edges.set_index("OSM").index, "class"] = n
    write_dataframe(HXOSM, OUTFILE, layer="osmnx")

print(3, dt.datetime.now() - START)
HEXAGON = read_dataframe("hexagon4.gpkg", layer="hexagon4-00")
HX = HEXAGON.dropna().reset_index(drop=True).rename(columns={"index": "hx"})
HX = HX.reset_index().rename(columns={"index": "hx"})
HX = HX[["hx", "geometry"]]
write_dataframe(HX.reset_index(), OUTFILE, layer="HX")

def get_hexlayer(hexagon, resolution):
    s = hexagon.copy()
    edge = hexagon.length.mean() / np.sqrt(7) / 6.0
    height = 2.0 * edge / np.sqrt(3)
    print(f"h0: {round(edge)}m {round(height)}m")
    s["geometry"] = s.buffer(height, cap_style="square", join_style="mitre")
    print("h1", dt.datetime.now() - START)
    r = get_subhexagon(s, resolution)
    r = r.reset_index(names=f"hx{str(resolution).zfill(2)}")
    return r.set_crs(CRS)


def get_overlay(hexagon, net_mx):
    print("h2", dt.datetime.now() - START)
    r = net_mx.overlay(hexagon, how="union", keep_geom_type=True)
    r = r.dropna().explode(ignore_index=True)
    for column in net_mx.columns.difference(["geometry"]):
        try:
            r[column] = r[column].apply(int)
        except ValueError:
            pass
    return r


def get_counts(hexagon, networkx, osmnetx):
    print("h4", dt.datetime.now() - START)
    r = hexagon.copy()
    r = r.set_index("hex_id")
    r[["OSM", "NX", "#"]] = 0, 0, 0
    ds = networkx["hex_id"].value_counts()
    r.loc[ds.index, "NX"] = ds
    ds = osmnetx["hex_id"].value_counts()
    r.loc[ds.index, "OSM"] = ds
    r.loc[r["NX"] > 0, "#"] = 1
    r.loc[r["OSM"] > 0, "#"] += 2
    return r.reset_index().set_crs(CRS)


def get_distance(vf1, vf2):
    r = pd.DataFrame(index=vf2.index)
    r[["distance", "key"]] = 65536.0, -1
    for k, v in vf1.iterrows():
        # h_fn = partial(hausdorff_distance, v['geometry'], densify=0.1)
        h_fn = partial(hausdorff_distance, v["geometry"], densify=1.0)
        d = vf2["geometry"].apply(h_fn)
        ix = d[d < r["distance"]].index
        r.loc[ix, "distance"] = d[idx]
        r.loc[ix, "key"] = k
    return r


def invert_dict(this_dict):
    r = {}
    for k, v in this_dict.items():
        for i in v:
            try:
                r[i].append(k)
            except KeyError:
                r[i] = [k]
    return r


def get_orientedbuffer(gf, width=4):
    r = gf["geometry"].apply(oriented_envelope)
    return r.buffer(width, cap_style="square", join_style="mitre", mitre_limit=0.0)


def get_buffer(gs, width=4, length=1.0e4):
    style = {"cap_style": "square", "join_style": "mitre", "mitre_limit": length}
    gs = gs["geometry"].copy()
    r = mp2.extend_lines(gs.to_frame(), 10.0, extension=length)
    r.index = gs.index
    return r.buffer(width, **style)



def get_overlap(gf1, gf2):
    r, s = [], pd.Index([])
    for k, v in gf1["geometry"].items():
        ix = gf2["geometry"].within(v)
        if ix.any():
            r.append(k)
        s = s.union(ix[ix].index)
    return (pd.Index(r), s)


def get_overlaps(r1, r2, width=4):
    print("h5", dt.datetime.now() - START)
    buffered = get_buffer(r1, width).to_frame("geometry")
    s = get_overlay(buffered, r2)
    buffered["hex_id"] = r1["hex_id"]
    keys = r1["hex_id"].drop_duplicates().values
    r, ix = pd.Index([]), pd.Index([])
    print(f"{len(keys)} {dt.datetime.now() - START}")
    for _, hkey in enumerate(keys):
        # print(f'{i} {dt.datetime.now() - START}')
        ix1 = buffered["hex_id"] == hkey
        ix2 = s["hex_id"] == hkey
        v = get_overlap(buffered[ix1], s[ix2])
        r = r.union(v[0])
        ix = ix.union(v[1])
    s.loc[ix, "m"] = 2
    print("h6", dt.datetime.now() - START)
    return (r, s)


def set_framekv(gf, ix, k, v):
    r = gf.copy()
    r.loc[ix, k] = v
    return r.reset_index()


print(4, dt.datetime.now() - START)
HX0 = HX
netx = HXNX
osmx = HXOSM
for n in range(5, 16):
    print(f"layer: {n}", dt.datetime.now() - START)
    hx = get_hexlayer(HX0, n)
    netx = get_overlay(hx[["hex_id", "geometry"]], netx)
    netx["nx2"] = netx.index
    osmx = get_overlay(hx[["hex_id", "geometry"]], osmx)
    osmx["osmx2"] = osmx.index
    m = f"{str(n).zfill(2)}"
    HX0 = get_counts(hx, netx, osmx)
    netx["m"], osmx["m"] = 3, 3
    idx = HX0.loc[HX0["#"] == 1, ["hex_id"]].set_index("hex_id").index
    netx = set_framekv(netx.set_index("hex_id"), idx, "m", 0)
    idx = HX0.loc[HX0["#"] == 2, ["hex_id"]].set_index("hex_id").index
    osmx = set_framekv(osmx.set_index("hex_id"), idx, "m", 0)
    idx, overlap = get_overlaps(
        netx[netx["m"] == 3], osmx[osmx["m"] == 3], width=2 * CENTRE2CENTRE
    )
    # HX0 key: {0: no rail, 1: centre-line only, 2: OSM only, 3: overlaps,
    #           4: all centre-line overlay}
    HX0["m"] = HX0["#"]
    idx = overlap[["osmx2", "m"]].drop_duplicates()
    idx = idx.drop_duplicates(subset="osmx2", keep=False).set_index("m")
    idx = pd.Index(idx.loc[2, "osmx2"]).drop_duplicates()
    osmx = set_framekv(osmx.set_index("osmx2"), idx, "m", 1)
    idx = overlap[["hex_id", "m"]].drop_duplicates()
    idx = idx.drop_duplicates(subset="hex_id", keep=False).set_index("m")
    idx = pd.Index(idx.loc[2, "hex_id"]).drop_duplicates()
    HX0 = set_framekv(HX0.set_index("hex_id"), idx, "m", 4)
    netx = set_framekv(netx.set_index("hex_id"), idx, "m", 2)
    write_dataframe(HX0, OUTFILE, layer=f"HX{m}")
    write_dataframe(netx, OUTFILE, layer=f"NTX{m}")
    write_dataframe(osmx, OUTFILE, layer=f"OSM{m}")

    HX0 = HX0[HX0["m"] == 3]
    idx = netx[netx["m"] == 3].set_index("NX").index.drop_duplicates()
    netx = HXNX.set_index("NX").loc[idx].reset_index()
    idx = osmx[osmx["m"] == 3].set_index("OSM").index.drop_duplicates()
    osmx = HXOSM.set_index("OSM").loc[idx].reset_index()
    if osmx.empty:
        break

print(5, dt.datetime.now() - START)
