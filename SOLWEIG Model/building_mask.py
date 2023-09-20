import fiona
import rasterio
import rasterio.mask


def adjustWithMask(mask_shape_file_name, rasterfn, newrasterfn):
    with fiona.open(mask_shape_file_name, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]

    geoms2d = [
        {
            "type": g["type"],
            "coordinates": [[xyz[0:2] for xyz in p] for p in g["coordinates"]],
        }
        for g in shapes
    ]

    with rasterio.open(rasterfn, 'r') as src:
        out_image, out_transform = rasterio.mask.mask(src, geoms2d, nodata=10000,
                                                      invert=True)  # Maybe have to change invert to false
        out_meta = src.meta

    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform})

    with rasterio.open(newrasterfn, "w", **out_meta) as dest:
        dest.write(out_image)


if __name__ == "__main__":
    maskshpfn = 'Maps/Tempe_MaskedBuildingsRoads.shp'
    rasterfn = 'output/2023-4-8-2100_mrt.tif'
    newrasterfn = 'output/masked-2023-4-8-2100-mrt.tif'

    adjustWithMask(maskshpfn, rasterfn, newrasterfn)
