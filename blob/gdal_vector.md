# GDAL Vector Geometry
The GDAL documentation is a huge list of classes, so this is my summary of
it so that I can figure out how it works.

The class hierarchy:

 * OGRGeometryFactory
 * OGREnvelope
 * OGRGeometry
   * OGRPoint
   * OGRCurve
     * OGRSimpleCurve
       * OGRLineString - a possibly closed ring
         * OGRLinearRing - a closed ring
       * OGRCircularString - Circular string from several arc circles 
   * OGRSurface
     * OGRCurvePolygon
       * OGRPolygon
         * OGRTriangle
     * OGRPolyhedralSurface
       * OGRTriangulatedSurface
   * OGRGeometryCollection
     * OGRMultiPoint
     * OGRMultiCurve
       * OGRMultiLineString
     * OGRMultiSurface
       * OGRMultiPolygon
     


## OGRGeometryFactory Class

Makes and destroys geometry, but also has some goodies like
an algorithm to make polygons from linestrings and calculations
of geodesic curve approximations. Random.

 * Makes new geometry from various file formats, like WKB and GeoJson.

 * void destroyGeometry(OGRGeometry*) - because free doesn't work
   across library boundaries in all systems.
 
 * OGRGeometry* createGeometry(geometry_type) - To use library's malloc.
 
 * OGRGeometry* forceToPolygon(*geometry) - Changes multipolygons.
   Destroys its input. Similarly, forceToLineString, forceToMultiPolygon,
   forceToMultiPoint, forceTo(*geometry, geometry_type).
   
 * organizePolygons - take many polygons and make a polygon with inner rings
   or a multipolygon.
   
 * transformWithOptions, for instance wrapping the date line.
 
 * approximateArcAngles - Stroke the arc of a circle to a linestring.
 
 * GetCurveParmeters - for the arc of a circle. The 'a' is missing in the header too.
 
 * curveToLineString - arc circle to an approximate line string.
   and curveFromLineString.

## OGREnvelope

Umm, bounding box?

 * MinX, MaxX, MinY, MaxY; All doubles.
 * Merge(another envelope&)
 * Intersect(another envelope&)
 * bool Intersects(other&)
 * bool Contains(other&)

## OGRGeometry

 * getDimension() and getCoordinateDimension(), which is 3 for XYZ and 4 for XYZM.
 
 * IsEmpty() and IsValid(), but even a OGRGeometry* MakeValid, that tries
   to fix things. IsSimple() to look for intersections and self-tangency.
   
 * empty() - Not a check. It clears the geometry. clone() makes a copy.
 
 * setCoordinateDimension(int) - Will squash z if you set to 2. or add it if you set3D().
   setMeasured(bool) - add or remove M coordinate.
  
 * assignSpatialReference(OGRSpatialReference *poSR), and getSpatialReference
   transform(OGRCoordinateTransformation *poCT), transformTo(OGRSpatialReference *poSR)
 
 * getEnvelope(OGREnvelope*) - bounding envelope.
 
 * getGeometryType() - listed in ogr_core.h.
 
 * getGeometryName() - from input. Don't free what's returned.
 
 * void flattenTo2D() - Converts all Z to 0.0.
 
 * Lots of import/export for well-known formats: GML, KML, Json.
 
 * bool hasCurveGeometry() - Tells you if any lines are curved.
 
 * OGRGeometry* getCurveGeometry(options) and getLinearGeometry.
   If you ask for Linear, and it was curved, you get an approximation.
 
 * closeRings() - Force rings to be closed.
 
 * segmentize(double max_length) - Turn long segments into short ones!
 
 * bool Intersects(OGRGeometry*), Equals(), Disjoing(OGRGeometry*),
   Touches(), Crosses(), Within(), Contains(), Overlaps()
   
 * OGRGeometry* boundary() - makes a new geometry object. Wonder what this
   does for a multipolygon.
   
 * double Distance(OGRGeometry* other), Distance3D(other)
 
 * OGRGeometry* ConvexHull()
 
 * OGRGeometry* Buffer(double dist, int quad_segments) - Make an outset version.
 
 * OGRGeometry* Intersection(OGRGeometry*), Union(), UnionCascaded(), Difference(), SymDifference()
 
 * Centroid(*point)
 
 * OGRGeometry* Simplify(double tolerance), SimplifyPreserveTopology()
 
 * OGRGeometry* DelaunayTriangulation(tolerance, only_edges) - triangulate the vertices.
   Uses GEOS library. Returns GeometryCollection of triangular polygons, or, if only edges,
   a multilinestring.
   
 * OGRGeometry* Polygonize() - turns sparse edges into polygons.
 
 * swapXY() - swap coordinates
 
 * toPoint(), toCurve(), toSimpleCurve, toLineString, toLinearRing, ...
 
 * createGEOSContext(), for working with GEOS.
 
## OGRPoint

 * getX(), getY(), getZ(), getM()
 
## OGRCurve

 * child elements are OGRPoint, for begin(), end()
 
 * get_Length(), get_Area()
 
 * StartPoint(*point), EndPoint(*point), getNumPoints(), getPointIterator()
 
 * get_IsClosed(), IsConvex() 
 
 * Value(distance, *point) - get a point at a distance along the curve
 
 * CurveToLine()
 
 * toSimpleCurve() - a downcast, CastToLineString, CastToCompoundCurve, CastToLinearRing
 
## OGRLineString

 * CurveToLine(double max_angle_step_size_degrees, options) - Turns curve geometry into
   a line string. And vice versa: getCurveGeometry()
 
 * get_Area() - It shows up here, not above in OGRGeometry().
 
## OGRLinearRing

 * isClockwise(), reverseWindingOrder()
 
 * isPointInRing(), isPointOnRingBoundary()
 
##  OGRCircularString

One or several arc-circles that make a circular string.

 * get_Length() length of the curve.
 
 * CurveToLine - makes a line string.
 
 * Value(double distance, point*) - fetch point at a distance along the curve
 
 * segmentize(double max_length) - create smaller segments.
 
## OGRSurface

 * get_Area()
 
 * OGRErr PointOnSurface(*point). Returns a point that is internal to the surface?

## OGRCurvePolygon

 * Child type is OGRCurve for begin() and end().
 
 * addRing(OGRCurve*), addRingDirectly(addRing), creates external or internal rings.
   removeRing(int)
 
 * OGRCurve* getExteriorRingCurve(), getNumInteriorRings()

## OGRPolygon

One outer ring, zero or more inner rings.

 * Child type is linear ring for begin() and end().
 
 * CurvePolyToPoly, and curve geometry translation.
 
 * getExteriorRing() and getInteriorRing(int) versus stealExteriorRing(int).
 
 * IsPointOnSurface(*point)
 
 * closeRings() - Force rings to be closed.

## OGRTriangle



## OGRPolyHedralSurface

 * Child type is OGRPolygon
 
 * get_Area()
 
 * addGeometry(), addGeometryDirectly()
 
 * getNumGeometries() getGeometryRef()
 
## OGRTriangulatedSurface

 * child is OGRTriangle
 
 * CastToPolyhedralSurface
 
 ## OGRGeometryCollection
 
  * Children are OGRGeometry for begin(), end().
  
  * get_Area(), get_Length(), getDimension(), getEnvelope()
  
  * getNumGeometries(), getGeometryRef(int), addGeometry(OGRGeometry),
    addGeometryDirectly(OGRGeometry*), removeGeometry()
    
  * CastToGeometryCollection - converts a derived class of geometry collection
    to a plain geometry collection... ah.
 
 ## OGRMultiPoint

 * Child type is Point for begin() and end(). 
 
 ## OGRMultiCurve
 
 * Child type is OGRCurve for begin(), end().
 
 ## OGRMultiLineString
 
 * Child type is OGRSurface
 
 * PointOnSurface() - returns a point on any surface?
 
 * CastToMultiPolygon

## OGRMultiSurface

 * Child type is OGRSurface
 
 * CastToMultiPolygon

## OGRMultiPolygon

 * Child type is OGRPolygon. They contained polygons are non-overlapping.
 
 * CastToMultiSurface
 
 
## Other functions in Geometry

```
const char CPL_DLL * OGRGeometryTypeToName( OGRwkbGeometryType eType );
OGRwkbGeometryType CPL_DLL OGRMergeGeometryTypes( OGRwkbGeometryType eMain,
                                                  OGRwkbGeometryType eExtra );
OGRwkbGeometryType CPL_DLL OGRMergeGeometryTypesEx( OGRwkbGeometryType eMain,
                                                    OGRwkbGeometryType eExtra,
                                                    int bAllowPromotingToCurves );
OGRwkbGeometryType CPL_DLL OGR_GT_Flatten( OGRwkbGeometryType eType );
OGRwkbGeometryType CPL_DLL OGR_GT_SetZ( OGRwkbGeometryType eType );
OGRwkbGeometryType CPL_DLL OGR_GT_SetM( OGRwkbGeometryType eType );
OGRwkbGeometryType CPL_DLL OGR_GT_SetModifier( OGRwkbGeometryType eType, int bSetZ, int bSetM );
int                CPL_DLL OGR_GT_HasZ( OGRwkbGeometryType eType );
int                CPL_DLL OGR_GT_HasM( OGRwkbGeometryType eType );
int                CPL_DLL OGR_GT_IsSubClassOf( OGRwkbGeometryType eType,
                                                OGRwkbGeometryType eSuperType );
int                CPL_DLL OGR_GT_IsCurve( OGRwkbGeometryType );
int                CPL_DLL OGR_GT_IsSurface( OGRwkbGeometryType );
int                CPL_DLL OGR_GT_IsNonLinear( OGRwkbGeometryType );
OGRwkbGeometryType CPL_DLL OGR_GT_GetCollection( OGRwkbGeometryType eType );
OGRwkbGeometryType CPL_DLL OGR_GT_GetCurve( OGRwkbGeometryType eType );
OGRwkbGeometryType CPL_DLL OGR_GT_GetLinear( OGRwkbGeometryType eType );

GDALVersionInfo
GDALCheckVersion
```

## Other functions in GDAL.h

 * GDALComputeMedianCutPCT
 * GDALDitherRGB2PCT
 * GDALComputeProximity
 * GDALFillNodata
 * GDALPolygonize
 * GDALFPolygonize
 * GDALSieveFilter
 * GDALCreateGenImgProjTransformer
 * GDALCreateGenImgProjTransformer2
 * GDALSetGenImgProjTransformerDstGeoTransform
 * GDALSetTransformerDstGeoTransform
 * GDALCreateReprojectionTransformer
 * GDALGCPTransform... lots more projection stuff.
 * Generate contour lines
 * Rasterize vector data
 * Grid vector data
 * Delaunay triangulation interface.  in gdal_alg.h
 