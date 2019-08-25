## Install Libraries on Ubuntu
sudo apt install gdal-bin gdal-data libalgorithms1 libbase1 libfileclasses1 libgdal-dev libgdal-doc libgdal-grass libgdal-java libgdal-perl libgdal20 libgdal-perl-doc libimageclasses1 liblasclasses1 libotbgdaladapters-6.6-1 libotbiogdal-6.6-1 node-srs pktools pktools-dev
python3-csvkit python3-rasterio r-cran-sf r-cran-stars rasterio

 - grass
 - gdal
 - node-srs which is spatial reference for node.js.
 - rasterio - command line for geospatial raster layers
 - pktools GDAL add-on for raster processing
 - libgeographic-dev - transforms coordinates

## Coordinate systems
Uganda Lat-long bounding box
* Latitude range: 4.5, -1.6
* Longitude range: 29.4, 35.1

## High-resolution settlement layer
The HRSL is in 30m-ish resolution,
where the -ish is because it's on a regular lat-long grid,
not in projection.
```
$ gdalinfo hrsl_uga_pop.tif
Driver: GTiff/GeoTIFF
Files: hrsl_uga_pop.tif
       hrsl_uga_pop.tif.ovr
       hrsl_uga_pop.tif.aux.xml
Size is 19536, 20540
Coordinate System is:
GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0],
    UNIT["degree",0.0174532925199433],
    AUTHORITY["EPSG","4326"]]
Origin = (29.573476945700001,4.228134261810000)
Pixel Size = (0.000277777780000,-0.000277777780000)
Metadata:
  AREA_OR_POINT=Area
  DataType=Generic
Image Structure Metadata:
  COMPRESSION=LZW
  INTERLEAVE=BAND
Corner Coordinates:
Upper Left  (  29.5734769,   4.2281343) ( 29d34'24.52"E,  4d13'41.28"N)
Lower Left  (  29.5734769,  -1.4774213) ( 29d34'24.52"E,  1d28'38.72"S)
Upper Right (  35.0001437,   4.2281343) ( 35d 0' 0.52"E,  4d13'41.28"N)
Lower Right (  35.0001437,  -1.4774213) ( 35d 0' 0.52"E,  1d28'38.72"S)
Center      (  32.2868103,   1.3753565) ( 32d17'12.52"E,  1d22'31.28"N)
Band 1 Block=128x128 Type=Float64, ColorInterp=Gray
  Min=0.461 Max=265.642 
  Minimum=0.461, Maximum=265.642, Mean=8.357, StdDev=6.256
  NoData Value=-1.79769300000000005e+308
  Overviews: 9768x10270, 4884x5135, 2442x2568, 1221x1284, 611x642, 306x321, 153x161
  Metadata:
    RepresentationType=ATHEMATIC
    STATISTICS_COVARIANCES=39.13179870426735
    STATISTICS_MAXIMUM=265.64248129217
    STATISTICS_MEAN=8.3569307275022
    STATISTICS_MINIMUM=0.46121201315056
    STATISTICS_SKIPFACTORX=1
    STATISTICS_SKIPFACTORY=1
    STATISTICS_STDDEV=6.255541439737
```

### PfPR from Malaria Atlas Project
The PfPR is also on a lat-long grid, not in projection,
this time in about 5km tiles.
```
$ gdalinfo 2019_Global_PfPR_2015.tif
   Driver: GTiff/GeoTIFF
   Files: 2019_Global_PfPR_2015.tif
          2019_Global_PfPR_2015.tif.aux.xml
   Size is 6927, 2607
   Coordinate System is:
   GEOGCS["WGS 84",
       DATUM["WGS_1984",
           SPHEROID["WGS 84",6378137,298.257223563,
               AUTHORITY["EPSG","7030"]],
           AUTHORITY["EPSG","6326"]],
       PRIMEM["Greenwich",0],
       UNIT["degree",0.0174532925199433],
       AUTHORITY["EPSG","4326"]]
   Origin = (-118.375000000000000,53.541623217000001)
   Pixel Size = (0.041666650000000,-0.041666650000000)
   Metadata:
     AREA_OR_POINT=Area
   Image Structure Metadata:
     COMPRESSION=LZW
     INTERLEAVE=BAND
   Corner Coordinates:
   Upper Left  (-118.3750000,  53.5416232) (118d22'30.00"W, 53d32'29.84"N)
   Lower Left  (-118.3750000, -55.0833333) (118d22'30.00"W, 55d 5' 0.00"S)
   Upper Right ( 170.2498845,  53.5416232) (170d14'59.58"E, 53d32'29.84"N)
   Lower Right ( 170.2498845, -55.0833333) (170d14'59.58"E, 55d 5' 0.00"S)
   Center      (  25.9374423,  -0.7708551) ( 25d56'14.79"E,  0d46'15.08"S)
   Band 1 Block=256x256 Type=Float64, ColorInterp=Gray
     Min=0.000 Max=0.875 
     Minimum=0.000, Maximum=0.875, Mean=0.048, StdDev=0.121
     NoData Value=-9999
     Overviews: 3464x1304, 1732x652, 866x326, 433x163
     Metadata:
       STATISTICS_MAXIMUM=0.87466576814651
       STATISTICS_MEAN=0.047807271344601
       STATISTICS_MINIMUM=0
       STATISTICS_STDDEV=0.12072164533226
       STATISTICS_VALID_PERCENT=20.82
```

## Districts from Local Burden of Disease
The LBD group creates admin files for the whole world. Unfortunately, they
stop at Admin 2 units. They are available at /home/j/WORK/11_geospatial/admin_shapefiles/current.
Looks like WGS84 and polygons.

```
$ ogrinfo lbd_standard_admin_2.shp lbd_standard_admin_2 | head -100
INFO: Open of `lbd_standard_admin_2.shp'
      using driver `ESRI Shapefile' successful.

Layer name: lbd_standard_admin_2
Metadata:
  DBF_DATE_LAST_UPDATE=2019-08-23
Geometry: Polygon
Feature Count: 47358
Extent: (-180.000000, -90.000000) - (180.000000, 83.658333)
Layer SRS WKT:
GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,
        AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.0174532925199433,
        AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]]
NAME_0: String (80.0)
NAME_1: String (80.0)
NAME_2: String (80.0)
geo_id: Real (19.11)
ad2_id: Integer64 (10.0)
ad0_parent: Integer64 (10.0)
ad1_parent: Integer64 (10.0)
ADM2_CODE: Integer64 (10.0)
ADM2_NAME: String (100.0)
ADM1_CODE: Integer64 (10.0)
ADM1_NAME: String (100.0)
ADM0_CODE: Integer64 (10.0)
ADM0_NAME: String (100.0)
OGRFeature(lbd_standard_admin_2):0
  NAME_0 (String) = Mexico
  NAME_1 (String) = Aguascalientes
  NAME_2 (String) = Aguascalientes
  geo_id (Real) = 1001143.00000000000
  ad2_id (Integer64) = 1
  ad0_parent (Integer64) = 143
  ad1_parent (Integer64) = 1
  ADM2_CODE (Integer64) = 1001143
  ADM2_NAME (String) = Aguascalientes
  ADM1_CODE (Integer64) = 1143
  ADM1_NAME (String) = Aguascalientes
  ADM0_CODE (Integer64) = 143
  ADM0_NAME (String) = Mexico
  POLYGON(...
```

## Districts from humdata
The features tell you the enclosing admin units.
```
layer uga_admbnda_adm3_UBOS_v5
	feature ADM0_EN String
	feature ADM0_PCODE String
	feature ADM1_EN String
	feature ADM1_PCODE String
	feature ADM2_EN String
	feature ADM2_PCODE String
	feature ADM3_EN String
	feature ADM3_PCODE String
ADM0_EN=Uganda,ADM0_PCODE=UG,ADM1_EN=ABIM,ADM1_PCODE=UG314,ADM2_EN=LABWOR,ADM2_PCODE=UG3141,ADM3_EN=ABIM,ADM3_PCODE=UG314101
ADM0_EN=Uganda,ADM0_PCODE=UG,ADM1_EN=ABIM,ADM1_PCODE=UG314,ADM2_EN=LABWOR,ADM2_PCODE=UG3141,ADM3_EN=ABIM TOWN COUNCIL,ADM3_PCODE=UG314102
ADM0_EN=Uganda,ADM0_PCODE=UG,ADM1_EN=ABIM,ADM1_PCODE=UG314,ADM2_EN=LABWOR,ADM2_PCODE=UG3141,ADM3_EN=ALEREK,ADM3_PCODE=UG314103
```

The Humdata level 3 admins have polygons and multipolygons, meaning the multipolygons have
more than one exterior ring and the polygons can have interior rings.

The proj.4 is
```
PROJ.4 : +proj=longlat +datum=WGS84 +no_defs
```

## Maplibrary Shapefile for Level 1 districts
In lat-long format.
There are four attributes

 * INDEX - unique within the file, not keyed to anything outside.
 * LEFT - Can be -1. Otherwise is an index.
 * RIGHT - Can be the same as this index.
 * LEVEL - 0-3


## The Plan

This should make patches to use for simulation. The target data structure is

 * A vector shapefile
 * with polygons that describe each patch
 * feature for which admin unit contains each patch
 * feature for population count
 * feature for PfPR for the patch

There are a few possible intermediate formats that might be helpful.
Let's list them and choose.

 1. A raster where each point says
    * how much population is there
    * how much PfPR
    * percentage of the point in each containing admin
    
 2. A raster where each point says
    * how much population is there
    * how much PfPR
    * a single containing admin

 3. A vector where each polygon says
    * how much population
    * what PfPR
    * containing admin unit
