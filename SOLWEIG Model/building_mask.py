import fiona
import rasterio
import rasterio.mask

def adjustWithMask(maskshpfn, rasterfn, newrasterfn):
    with fiona.open(maskshpfn, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]

    geoms2d = [
        {
            "type": g["type"],
            "coordinates": [[xyz[0:2] for xyz in p] for p in g["coordinates"]],
        }
        for g in shapes
    ]

    with rasterio.open(rasterfn) as src:
        out_image, out_transform = rasterio.mask.mask(src, geoms2d, nodata=600, invert=True) #Maybe have to change invert to false
        out_meta = src.meta

    out_meta.update({"driver": "GTiff",
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform})

    with rasterio.open(newrasterfn, "w", **out_meta) as dest:
        dest.write(out_image)


if __name__ == "__main__":
    maskshpfn = 'Maps/Tempe_buildings_NAD.shp'
    rasterfn = 'Maps/2100202328_2100_mrt22.tif'
    newrasterfn = 'Maps/MASKED_NAD_c2100202328_2100_mrt22.tif'

    adjustWithMask(maskshpfn, rasterfn, newrasterfn)