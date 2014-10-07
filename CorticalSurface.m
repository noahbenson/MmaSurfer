(* CorticalSurface.m
 *
 * Utility functions for dealing with (mostly FreeSurfer) cortical surface data in Mathematica.
 *
 * Copyright (C) 2013-2014 by Noah C. Benson.
 * This file is part of the MmaSurfer library, which is provided under the terms of the Eclipse
 * Public License, version 1.0. See the accompanying LICENSE file for more information.
 *)

(**************************************************************************************************)
BeginPackage["CorticalSurface`", {"ComputationalGeometry`"}];
Unprotect["CorticalSurface`*", "CorticalSurface`Private`*"];
ClearAll[ "CorticalSurface`*", "CorticalSurface`Private`*"];

EdgeLengths::usage = "EdgeLengths[s] yields a list of the lengths (Euclidean norm) of each edge in the EdgeList[s] where s may be a surface or a map.";
FaceList::usage = "FaceList[s] yields a list of all faces in the surface or map s, if any, as lists of indices into VertexList[s].";
Faces::usage = "Faces[s] is an alias for FaceList[s].";
FaceCount::usage = "FaceCount[s] yields the count of all faces in the surface or map s.";
FaceAngles::usage = "FaceAngles[s] yields a list of the angles (in radians) of each edge in the FaceList[s] where s may be a surface or a map. FaceAngles[s, X] yields the angles for s if the vertices in s had coordinates equal to those in the list X.";
FacesIndex::usage = "FacesIndex[s] yields a list of the indices into FaceList[s] such that the i'th element in FacesIndex[s] is the list of indices at which the i'th vertex in VertexList[s] appears in FaceList[s].";
Longitude::usage = "Longitude is a keyword that represents the longitude (identical to azimuth angle) in spherical coordinates.";
Latitude::usage = "Latitude is a keyword that represents the latitude in spherical coordinates.";
Radius::usage = "Radius is a keyword that represents the radius of a surface-to-map projection (see SurfaceProjection).";

CartesianToSpherical::usage = "CartesianToSpherical[data] yields an equivalent dataset to data except that the parts of data corresponding to (x,y,z) coordinates are replaced with (\[Theta],\[Phi]) coordinates. The data parameter must be a list of either lists or rules OR a single list or a single rule. In data is a list of lists, then the first three elements of each sublist is taken to be the (x,y,z) part; if data is a list of rules or a rule, then the first element of the rule is taken to be the (x,y,z) part. CartesianToSpherical yields data in the same format as it is given, but with the appropriate replacements."
CartesianToSpherical::badfmt = "Bad data format given to CartesianToSpherical";
CartesianToSphercial::badopt = "Bad option `1` given to CartesianToSpherical: `2`";

SphericalToCartesian::usage = "SphericalToCartesian is identical to CartesianToSpherical, but performs the inverse transform."
SphericalToCartesian::badfmt = "Bad data format given to SphericalToCartesian";
SphercialToCartesian::badopt = "Bad option `1` given to SphericalToCartesian: `2`";

SphericalCoordinateStyle::usage = "SphericalCoordinateStyle is an option for functions such as Surface that accept spherical coordinates. SphericalCoordinateStyle describes coordinates of the form {a, b} such that a SphericalCoordinateStyle value of {x, y} indicates that a represents measurement x and b represents measurement y. The possible measurements are Azimuth/Longitude (which are identical), Latitude, and PolarAngle. At least one of x and y must be Azimuth or Longitude and one of them must be either Latitude or PolarAngle. The ConvertCoordinates function can be used to convert between representations.";
SphericalPolarAngle::usage = "SphericalPolarAngle is a keyword that either represents polar angle data of the visual field or the angle from the pole (z-axis) in spherical coordinate data.";
SphericalAzimuth::usage = "SphericalAzimuth is a keyword that represents the azimuth angle (identical to longitude) in spherical coordinates.";

ConvertCoordinates::usage = "ConvertCoordinates[data, fromStyle -> toStyle] yields a transformation of data such that its coordinates are converted from fromStyle to toStyle; fromStyle and toStyle may be Cartesian or any argument appropriate for SphericalCoordinateStyle. The data parameter may be a list of points or a list of rules whose heads are points.";
ConvertCoordinates::badfmt = "Bad input to ConvertCoordinates: `1`";

Surface::usage = "Surface[data] yields a cortical surface representation of the cortical surface data given in data. The cortical surface data must be a list of rules, which may be in spherical or cartesian coordinates. If the data is in spherical coordinates, then the SphericalCoordinateStyle option specifies how the coordinates should be interpreted; if the coordinates for each point are {a, b}, then the SphericalCoordinateStyle option specifies {a-style, b-style} where each of a-style and b-style must come from these lists and must not come from the same list: {SphericalAzimuth, Longitude} and {Latitude, SphericalPolarAngle}. See ?SphericalCoordinateStyle for more information.
The cortical surface representation yielded by this function is actually a symbol with many tagged values. These each have their own additional documentation and are listed here:
  VertexList[surf],
  FaceList[surf],
  Field[surf],
  SurfaceQ[surf]
  Normal[surf]
  VertexFilter[surf]
  FaceFilter[surf]";
Surface::badarg = "Bad argument `1` to Surface: `2`";

SurfaceFromVTK::usage = "SurfaceFromVTK[filename] yields a cortical surface representation of the cortical surface data given in the VTK file specified; see Surface and ReadVTK for more information.";

VertexFilter::usage = "VertexFilter is an option to Surface and SurfaceProjection; in both cases, it must be a function that takes as arguments a vertex position (in cartesian or map coordinates) and a field value (which may be None). Only those vertices for which the filter function yields True are included in the constructed object. In the case of maps, vertices that are filtered are still included in the domain, but not in the result. VertexFilter[surface] and VertexFilter[map] yield the surface or map's filter.";
FaceFilter::usage = "FaceFilter is an option to Surface and SurfaceProjection; in both cases, it must be a function that takes as arguments a list of vertex indices which together form a face. Only those faces for which the filter function yields True are included in the constructed object. In the case of maps, faces that are filtered are still included in the domain, but not in the result. FaceFilter[surface] and FaceFilter[map] yield the surface or map's filter.";
Polygons::usage = "Polygons[surf] yields a list of all polygons in the surface surf, if any, by cross-referencing VertexList[surf] and Faces[surf]; these polygons are listed in the same order, with the same vertex order, as in Faces[surf].";
Field::usage = "Field[surf] yields a list of all field values associated with cortical surface surf.";
FieldQ::usage = "FieldQ[obj] yields true if and only if obj has a field such that Field[obj] yields a list.";
ProjectedSurface::usage = "ProjectedSurface[map] yields the surface object from which the map was projected.";
SurfaceQ::usage = "SurfaceQ[s] yields true if and only if s is a surface object, as created by Surface[].";
SurfaceName::usage = "SurfaceName[x] yields the symbol to which all the given surface's geometric data is attached. Note that while SurfaceName[s] is always a map with the same vertices, faces, edges, etc. as s, it is not necessarily s, and does not necessarily have the same field as s.";
VertexIndexDispatch::usage = "VertexIndexDispatch[s] yields a dispatch that converts vertex indices from the vertex list given to the Surface function originally for surface s to the index used in the final version of the surface that was produced by the Surface function. The main use case for this dispatch is that the indices in the face list are automatically converted to new indices when a VertexFilter is added to a surface.";

Domain::usage = "Domain[map] yields the Normal of the section of the surface from which map was derived.";
DomainIndices::usage = "DomainIndices[map] yields a list of the indices into Normal[ProjectedSurface[map]] for the points that compose the domain of map.";
ProjectionDispatch::usage = "ProjectionDispatch[map] yields a dispatch list (ie, for use with Replace and /.) that, when given a cartesian point on a surface projects it to the map.";
ProjectionShear::usage = "ProjectionShear[map] yields the shearing matrix associated with map.";
ProjectionRotation::usage = "ProjectionRotation[map] yields the rotation matrix associated with map.";
ProjectionTransform::usage = "ProjectionTransform[map] yields a function that converts cartesian coordinates into map coordinates for the projection that creates map.";
InverseProjectionDispatch::usage = "InverseDispatch[map] yields the dispatch function that, when used to replace something, inverses the transformation from the surface to the map.";
InverseProjectionTransform::usage = "InverseProjectionTransform[map] yields a function that inverts the transform from map's surface to map."
MapQ::usage = "MapQ[x] yields true if and only if x is a map created as a projection of a surface.";
MapName::usage = "MapName[x] yields the symbol to which all the given map's geometric data is attached. Note that while MapName[m] is always a map with the same vertices, faces, edges, etc. as m, it is not necessarily m, and does not necessarily have the same field as m.";
MapConvexHull::usage = "MapConvexHull[m] yields the indices of the convex hull of the vertices in the given map m.";
MapHull::usage = "MapConvexHull[m] yields the indices of the general hull of the vertices in the given map m. The general hull differs from the convex hull in that all vertices with faces ";

MergeSurfaces::usage = "MergeSurfaces[s1, s2, f] yields a surface object such that any point p appearing in either s1 or s2 is represented in the new surface and has field value equal to
 f[f1, f2] where f1 and f2 are the field values at p in surface s1 and s2 respectively (or None if no field is defined at that point).";
MergeSurfaces::dup = "Surface given to MergeSurfaces has duplicate vertices";

SurfaceProjection::usage = "SurfaceProjection[surf] yields a surface projection of the region within Pi/3 radians of the occipital pole, with V1 oriented along the positive x-axis. To alter the location of the map projection constructed, the following options may be used:
  Duplicate may be a surface projection whose projection should be duplicated (other arguments overwrite the duplicated surface's options).
  Center must be a 3D cartesian point representing the center of the map to be constructed.
  Radius must be a number greater than 0 that indicates the distance from the center that the map should encompass.
  OrientPoint must be a point or a point -> an angle (which by default is 0) to use in orienting the map; the point is oriented so that it lies at the given angle in a polar-coordinate version of the map (i.e., x -> Pi/2 indicates that point x should, in the resulting map, lie on the y-axis).
  ProjectionShear must be a shear matrix that is applied to the final map after orientation.
  Scale must be a real number greater than 0 that instructs the projection to scale the resulting map such that the average edge length is equal to the given argument. By default, this is 1; None may be provided to indicate that no scaling should be used (in which case the result is in terms of radians).";
SurfaceProjection::badarg = "Bad argument `1` to SurfaceProjection: `2`";
SurfaceProjection::badfmt = "Transform argument given in unrecognized format";
SurfaceRotation::usage = "SurfaceRotation[map] yields the rotation matrix that produces the map.";

OrientPoint::usage = "OrientPoint is a keyword that represents the ultimate orientation of a surface-to-map projection (see SurfaceProjection). OrientPoint -> (x -> 0) is identical to OrientPoint -> x; otherwise, OrientPoint -> (x -> y) ensures that in the final map projection, the point x will be projected to point {qx,qy} such that the ArcTan[qx, qy] == y.";
OrientMatrix::usage = "OrientMatrix[map] yields the orientation matrix applied to the transformation that produced the surface projection, map.";
Duplicate::usage = "Duplicate is an option to SurfaceProjection which specifies that the projection should, unless otherwise specified in the options list, use the transformation arguments used by the given map.";

MapPlot::usage = "MapPlot[map] yields a 2D map plot of the given map; all options for Graphics[] may be given.";
MapMeshPlot::usage = "MapMeshPlot[map] is identical to MapPlot except that it plots the map as a mesh instead of as a set of polygons.";
MapMeshPlot::badarg = "Bad arguent given to MapMeshPlot: `1`";

WithVertexFilter::usage = "WithVertexFilter[map, filt] and WithVertexFilter[surf, filt] yield a map or surface object that is identical to map or surf except that it has an additional filter applied to it.";
WithFaceFilter::usage = "WithFaceFilter[map, filt] and WithFaceFilter[surf, filt] yield a map or surface object that is identical to map or surf except that it has an additional filter applied to it.";
WithField::usage = "WithField[surf, field] yields a surface object that is identical to surf except that its field has been replaced with the given field. Rule (field -> surf) may be used as a shorthand for WithField[surf, field].";
WithField::incompat = "WithField given surface and field with a different numbers of vertices.";
WithVertexList::usage = "WithVertexList[surf, X] yields a new surface identical to surf, but with the given list of vertices X instead in the place of surf's vertices.";
WithVertexList::incompat = "WithVertexList given a surface and a vertex list with different sizes.";

ReadVTK::usage = "ReadVTK[filename] reads the given VTK file, including fields, and yields a list of rules, {x,y,z} -> field. The filename may be a URL.";
ReadVTK::nofile = "No such file: `1`";
ReadVTK::flerr = "Unspecified error while reading file: `1`";

SurfacePlot::usage = "SurfacePlot[surf] yields a plot of the surface object surf. All options available to Graphics3D may be given,";

SurfaceResample::usage = "SurfaceResample[surf1, surf2] yields a surface equivalent to surf2 but such that the field of the surface has been resampled from surf1. All options except the Field option of Surface are accepted and passed verbatim; a Method option may also specify Nearest (default) for nearest-neighbor interpolation, Interpolation, or List interpolation, for their respective functions. In the latter two cases, A list may be given instead of the argument such that the first argument is Interpolation or LitInterpolation and the remaining elements of the list are options to pass to these functions; e.g. Method -> {Interpolation, InterpolationOrder -> 4}. Note that surf1 -> surf2 is equivalent to SurfaceResample[surf1, surf2] but self-hashes.";

NeighborhoodList::usage = "NeighborhoodList[s] yields a list of length N (where N is the number of vertices in s) of the neighboring vertices of each vertex; each entry k of the list is a list of the integer id's of the neighbors of the kth vertex. The neighbor id's are sorted such that they are listed in a counter-clockwise order around vertex k starting from the x-axis. The argument s may be either a map or a surface.";
NeighborhoodAngles::usage = "NeighborhoodAngles[s] yields a list of the angles of each angle in in the surface or map s such that for the i'th entry in the VertexList[s] of the surface or map s, the i'th entry of the NeighborhoodAngles[s] list is a list of the angles around the i'th vertex such that the first angle is the angle between the virst two entries of the i'th element of the NeighborhoodList[s]. NeighborhoodAngles[s, X] yields the neighborhood angles for s if the vertices of s were replaced with the vertices in the list X.";
NeighborhoodBisectors::usage = "NeighborhoodBisectors[s] yeilds a list of the vectors that bisect each angle in the neighborhood of each vertex in the surface or map s. The result is equivalent to that produced by NeighborhoodAngles[s], but each angle is replaced by the unit vector that bisects that angle. NeighborhoodBisectors[s, X] yields the neighborhood bisecting vectors for the points given by the coordinate matrix X.";
NeighborhoodEdgeLengths::usage = "NeighborhoodEdgeLengths[s] yields a list of the edge lengths for the neighborhood of each vertex in the surface or map s. The results are in the same order as NnighborhoodList[s] such that for the neighborhood list L, and the neighborhood edge length list G, G[[i,j]] is the length of the edge from vertex i to the vertex L[[i,j]].";

EdgesEnergy::usage = "EdgesEnergy[s, X] yeilds the scalar energy for the given surface or map s using the vertex coordinates given in X that is the result of deformation of the edge lengths in s. EdgesEnergy[s] is equivalent to EdgesEnergy[s, VertexList[s]], which is always 0. The energy of edge deformation is the sum of the squares of the deviation between the edge lengths in s and the edge lengths in X divided by the number of edges total.";
EdgesGradient::usage = "EdgesGradient[s, X] yields the gradient matrix for the given surface or map s using the vertex coordinates given in X that is the result of deformation of the edge lengths in s. The result is the gradient of the EdgesEnergy[s,X]. EdgesGradient[s] is equivalent to EdgesGradient[s,VertexList[s]], which is a list of zero-vectors.";
AnglesEnergy::usage = "AnglesEnergy[s, X] yeilds the scalar energy for the given surface or map s using the vertex coordinates given in X that is the result of deformation of the angles in s. AnglesEnergy[s] is equivalent to AnglesEnergy[s, VertexList[s]], which is always 0. The energy of angle deformation is the sum of the squares of the deviation between the angles in s and the angles in X divided by the number of angles total.";
AnglesGradient::usage = "AnglesGradient[s, X] yields the gradient matrix for the given surface or map s using the vertex coordinates given in X that is the result of deformation of the angles in s. The result is the gradient of the AnglesEnergy[s,X]. AnglesGradient[s] is equivalent to AnglesGradient[s,VertexList[s]], which is a list of zero-vectors.";
TrianglesEnergy::usage = "TrianglesEnergy[s, X] is an alternate verision of AnglesEnergy[s, X], which should be better at considering the energy of flipped triangles.";
TrianglesGradient::usage = "TrianglesGradient[s, X] is an alternate verision of AnglesGradient[s, X], which should be better at considering the energy of flipped triangles.";
CorticalPotentialField::usage = "CorticalPotentialField[s, options...] yields a symbol f which, when evaluated as f[X] for a numeric list X with dieensios equial to those of VertexList[s], yields the potential energy of the vertex configuration X according to the options given. Similarly, Gradient[f, X] yields the flattened gradient vector for the conformation X and is appropriate for passing to optimization functions that require a gradient such as FindArgMin.
The following options are accepted:
  EdgesConstant must be a number >= 0 and specifies the relative strength of the edge forces in the potential function; any occurance of Automatic is replaced by 1/n where n is the number of edges in s (default: Automatic).
  AnglesConstant must be a number >= 0 and specifies the relative strength of the angle forces in the potential function; any occurance of Automatic is replaced by 1/n where n is the number of angles in s (default: Automatic).
  TrianglesConstant must be a number >= 0 and specifies the relative strength of the triangles force in the potential function; any occurante of Automatic is replaced by 1/n where n is the number of angles in s (default: 0).
  Other options may be specified if they are defined via the CorticalPotentialTerm interface.";
CorticalPotentialField::badterm = "Bad term given to CorticalPotentialField: `1`";
EdgesConstant::usage = "EdgesConstant is an option to CorticalPotentialField that specifies the raltive strength of the edge forces.";
AnglesConstant::usage = "AnglesConstant is an option to CorticalPotentialField that specifies the raltive strength of the angle forces.";
TrianglesConstant::usage = "TrianglesConstant is an option to CorticalPotentialField that specifies the raltive strength of the triangle forces.";
CorticalPotentialTerm::usage = "CorticalPotentialTerm[s, name -> options] yields the pair {energyFunction, gradientFunction} for the given name with options and surface or map s. This form is partially protected and new values can be defined for it. Any name and option that is defined is a valid argument for the CorticalPotentialField function.";
CorticalPotentialTerm::badarg = "Bad argument given to CorticalPotentialTerm: `1`";

ColorCortex::usage = "ColorCortex[instructions...] yields a color function for a surface or map that follows the instructions given. Each instruction should be of the form <column> -> <color-instruction> where the column is a colum index in the field matrix of the surface or map that is to be colorized. An instruction can be given without the column rule (ie, just <color-instruction>) to indicate that the entire row (ie, when the field is just a vector). Color instructions may be PolarAngle, Eccentricity, Curvature, or a function that takes an argument and yields a color. No field row or cell that is either None or $Failed will ever pass a match, and any function that yields Indeterminate or $Failed will be skipped in the coloring. Instructions are attempted one at a time until there is a match, and if there is no match, then Gray is used.
New cortical colors can be added by interfacing with the CorticalColor form.";
CorticalColor::usage = "CorticalColor[tag] yields the cortical color instruction for the given tag. The CorticalColor form is only partially protected and may be assigned single-argument values such as CorticalColor[\"MyCustomColorScheme\"] = {{0,10}, {Red,Yellow,Green,Cyan}}. In this example, the color instruction \"MyCustomColorScheme\" would then be valid, and would blend the given colors over the range 0 to 10. A Function may also be given as the value, in which case the function is given the field value and expected to yield a color.
See also ColorCortex.";

Curvature::usage = "Curvature is a keyword that can be used to refer to curvature values; it is automatically defined by the CorticalSurface package to include a CorticalColor function as well.";

(**************************************************************************************************)
Begin["`Private`"];

Unprotect[VertexList, EdgeList, EdgeCount, VertexCount];

(* #ReadVTK ***************************************************************************************)
ReadVTK[filename_String] := Catch[
  With[
    {fl = StringToStream[
       (* we do all this instead of a normal stream because only import correctly handles URLs *)
       Apply[
         StringJoin,
         (* The stream for this particular file *)
         Check[
           Import[filename, "Character8"],
           (Message[ReadVTK::nofile, filename];
            Throw[$Failed])]]]},
    Check[
      With[
        {dat = Reap[
           While[ 
             Find[fl, "FIELD"] =!= EndOfFile,
             With[
               {tmp = Read[fl, {Word, Number, Number, Word}]},
               Sow[Partition[ReadList[fl, Number, tmp[[2]]*tmp[[3]]], tmp[[2]]]]]];
           Close[fl];
          ][[2, 1]]},
        (* At this point, we just need to put dat in the form of a list of {x,y,z} -> value rules *)
        MapThread[
          (#1 -> If[Length[#2] == 1, #2[[1]], #2]) &,
          {Import[filename, "VertexData"], Map[Flatten, Transpose[dat]]}]],
      (Message[ReadVTK::flerr, filename]; Close[fl]; $Failed)]]];

(* #CartesianToSpherical **************************************************************************)
Options[CartesianToSpherical] = {SphericalCoordinateStyle -> {Longitude, Latitude}};
CartesianToSpherical[dat_, OptionsPattern[]] := (Message[CartesianToSpherical::badfmt]; $Failed);
CartesianToSpherical[dat : {{_, _, _} ..}, OptionsPattern[]] := With[
  {scFn = Replace[
     OptionValue[SphericalCoordinateStyle],
     {{Latitude, Longitude|SphericalAzimuth} -> Function[
        {ArcSin[#[[3]]/Norm[#]],
         If[#[[1]] == 0 && #[[2]] == 0, 0, ArcTan[#[[1]], #[[2]]]]}],
      {Longitude|SphericalAzimuth, Latitude} -> Function[
        {If[#[[1]] == 0 && #[[2]] == 0, 0, ArcTan[#[[1]], #[[2]]]],
         ArcSin[#[[3]]/Norm[#]]}],
      {SphericalPolarAngle, Longitude|SphericalAzimuth} -> Function[
        {ArcCos[#[[3]]/Norm[#]],
         If[#[[1]] == 0 && #[[2]] == 0, 0, ArcTan[#[[1]], #[[2]]]]}],
      {Longitude|SphericalAzimuth, SphericalPolarAngle} -> Function[
        {If[#[[1]] == 0 && #[[2]] == 0, 0, ArcTan[#[[1]], #[[2]]]],
         ArcCos[#[[3]]/Norm[#]]}],
      _ :> (
        Message[CartesianToSpherical::badopt, SphericalCoordinateStyle, "Style not recognized"];
        $Failed)}]},
  If[scFn === $Failed, $Failed, Map[scFn, dat]]];
CartesianToSpherical[dat : {Rule[{_, _, _}, _] ..}, OptionsPattern[]] := With[
  {scFn = Replace[
     OptionValue[SphericalCoordinateStyle],
     {{Latitude, Longitude|SphericalAzimuth} -> Function[
        Rule[
          {ArcSin[#[[1,3]]/Norm[#[[1]]]],
           If[#[[1,1]] == 0 && #[[1,2]] == 0, 0, ArcTan[#[[1,1]], #[[1,2]]]]},
          #[[2]]]],
      {Longitude|SphericalAzimuth, Latitude} -> Function[
        Rule[
          {If[#[[1,1]] == 0 && #[[1,2]] == 0, 0, ArcTan[#[[1,1]], #[[1,2]]]],
           ArcSin[#[[1,3]]/Norm[#[[1]]]]},
          #[[2]]]],
      {SphericalPolarAngle, Longitude|SphericalAzimuth} -> Function[
        Rule[
          {ArcCos[#[[1,3]]/Norm[#[[1]]]],
           If[#[[1,1]] == 0 && #[[1,2]] == 0, 0, ArcTan[#[[1,1]], #[[1,2]]]]},
          #[[2]]]],
      {Longitude|SphericalAzimuth, SphericalPolarAngle} -> Function[
        Rule[
          {If[#[[1,1]] == 0 && #[[1,2]] == 0, 0, ArcTan[#[[1,1]], #[[1,2]]]],
           ArcCos[#[[1,3]]/Norm[#[[1]]]]},
          #[[2]]]],
      _ :> (
        Message[CartesianToSpherical::badopt, SphericalCoordinateStyle, "Style not recognized"];
        $Failed)}]},
  If[scFn === $Failed, $Failed, Map[scFn, dat]]];
CartesianToSpherical[dat : {_?AtomQ, _?AtomQ, _?AtomQ}, OptionsPattern[]] := First[
  CartesianToSpherical[
    {dat},
    SphericalCoordinateStyle -> OptionValue[SphericalCoordinateStyle]]];
CartesianToSpherical[dat : Rule[{_?AtomQ, _?AtomQ, _?AtomQ}, _], OptionsPattern[]] := First[
  CartesianToSpherical[
    {dat},
    SphericalCoordinateStyle -> OptionValue[SphericalCoordinateStyle]]];

(* #SphericalToCartesian **************************************************************************)
Options[SphericalToCartesian] = {
  SphericalCoordinateStyle -> {Longitude, Latitude}};
SphericalToCartesian[dat_, OptionsPatterm[]] := (Message[SphericalToCartesian::badfmt]; $Failed);
SphericalToCartesian[dat : {{_, _} ..}, OptionsPattern[]] := With[
  {scFn = Replace[
     OptionValue[SphericalCoordinateStyle],
     {{Latitude, Longitude|SphericalAzimuth} -> Function[
        Append[Cos[#[[1]]] * {Cos[#[[2]]], Sin[#[[2]]]}, Sin[#[[1]]]]],
      {Longitude|SphericalAzimuth, Latitude} -> Function[
        Append[Cos[#[[2]]] * {Cos[#[[1]]], Sin[#[[1]]]}, Sin[#[[2]]]]],
      {SphericalPolarAngle, Longitude|SphericalAzimuth} -> Function[
        Append[Sin[#[[1]]] * {Cos[#[[2]]], Sin[#[[2]]]}, Cos[#[[1]]]]],
      {Longitude|SphericalAzimuth, SphericalPolarAngle} -> Function[
        Append[Sin[#[[2]]] * {Cos[#[[1]]], Sin[#[[1]]]}, Cos[#[[2]]]]],
      _ :> (
        Message[SphericalToCartesian::badopt, SphericalCoordinateStyle, "Style not recognized"];
        $Failed)}]},
  If[scFn === $Failed, $Failed, Map[scFn, dat]]];
SphericalToCartesian[dat : {Rule[{_, _}, _] ..}, OptionsPattern[]] := With[
  {scFn = Replace[
     OptionValue[SphericalCoordinateStyle],
     {{Latitude, Longitude|SphericalAzimuth} -> Function[
        Rule[
          Append[Cos[#[[1,1]]] * {Cos[#[[1,2]]], Sin[#[[1,2]]]}, Sin[#[[1,1]]]],
          #[[2]]]],
      {Longitude|SphericalAzimuth, Latitude} -> Function[
        Rule[
          Append[Cos[#[[1,2]]] * {Cos[#[[1,1]]], Sin[#[[1,1]]]}, Sin[#[[1,2]]]],
          #[[2]]]],
      {SphericalPolarAngle, Longitude|SphericalAzimuth} -> Function[
        Rule[
          Append[Sin[#[[1,1]]] * {Cos[#[[1,2]]], Sin[#[[1,2]]]}, Cos[#[[1,1]]]],
          #[[2]]]],
      {Longitude|SphericalAzimuth, SphericalPolarAngle} -> Function[
        Rule[
          Append[Sin[#[[1,2]]] * {Cos[#[[1,1]]], Sin[#[[1,1]]]}, Cos[#[[1,2]]]],
          #[[2]]]],
      _ :> (
        Message[SphericalToCartesian::badopt, SphericalCoordinateStyle, "Style not recognized"];
        $Failed)}]},
  If[scFn === $Failed, $Failed, Map[scFn, dat]]];
SphericalToCartesian[dat : {_?AtomQ, _?AtomQ}, OptionsPattern[]] := SphericalToCartesian[
  {dat}, 
  SphericalCoordinateStyle -> OptionValue[SphericalCoordinateStyle]];
SphericalToCartesian[dat : Rule[{_, _}, _], OptionsPattern[]] := SphericalToCartesian[
  {dat}, 
  SphericalCoordinateStyle -> OptionValue[SphericalCoordinateStyle]];

(* #ConvertCoordinates ****************************************************************************)
ConvertCoordinates[data_List, _] := (
  Message[ConvertCoordinates::badfmt, "Cannot deduce conversion of data"];
  $Failed);
ConvertCoordinates[data:{Rule[{_,_,_},_]..}, Rule[Cartesian, toStyle_]] := If[
  toStyle === Cartesian, data,
  CartesianToSpherical[data, SphericalCoordinateStyle -> toStyle]];
ConvertCoordinates[data:{{_,_,_}..}, Rule[Cartesian, toStyle_]] := If[
  toStyle === Cartesian, data,
  CartesianToSpherical[data, SphericalCoordinateStyle -> toStyle]];
ConvertCoordinates[data:{Rule[{_,_},_]..}, Rule[fromStyle_, Cartesian]] := SphericalToCartesian[
  data,
  SphericalCoordinateStyle -> fromStyle];
ConvertCoordinates[data:{{_,_}..}, Rule[Cartesian, toStyle_]] := SphericalToCartesian[
  data,
  SphericalCoordinateStyle -> fromStyle];
ConvertCoordinates[data:{Rule[{_,_},_]..}, Rule[fromStyle_, toStyle_]] := CartesianToSpherical[
  SphericalToCartesian[data, SphericalCoordinateStyle -> fromStyle],
  SphericalCoordinateStyle -> toStyle];
ConvertCoordinates[data:{{_,_}..}, Rule[Cartesian, toStyle_]] := CartesianToSpherical[
  SphericalToCartesian[data, SphericalCoordinateStyle -> fromStyle],
  SphericalCoordinateStyle -> toStyle];

(* #FieldQ ****************************************************************************************)
FieldQ[obj_] := (Head[Field[obj]] =!= Field);

(* #Surface ***************************************************************************************)
Options[Surface] = {
  SphericalCoordinateStyle -> {Longitude, Latitude}, 
  VertexFilter -> None,
  FaceFilter -> None,
  Field -> None,
  Faces -> None,
  FaceList -> None,
  Polygons -> None};
Surface[data: {{_,_,_}..}, OptionsPattern[]] := Catch[
  With[
    {sym = (SetAttributes[Evaluate[#], Temporary]; #)&@Unique["surf"],
     field = Replace[
       OptionValue[Field],
       {l_List /; Length[l] == Length[data] :> l,
        None -> Table[None, {Length[data]}],
        _ :> (
          Message[Surface::badarg, Field, "must have same length as vertices"];
          Throw[$Failed])}],
     faces = Replace[
       Replace[OptionValue[FaceList], None :> OptionValue[Faces]],
       {(l : {{_Integer, _Integer, _Integer..}..} 
            /; Min[Flatten[l]] >= 1 && Max[Flatten[l]] <= Length[data]) :> l,
        None -> None,
        _ :> (
          Message[Surface::badarg, FaceList, "must consist of index lists"];
          Throw[$Failed])}],
     polygons = Replace[
       OptionValue[Polygons],
       {l_List /; TrueQ[Union[data] == Union[Join[data, Flatten[l, 2]]]] :> l,
        None -> None,
        _ :> (
          Message[Surface::badarg, FaceList, "must consist of index lists"];
          Throw[$Failed])}],
     uFiltFn = Replace[
       OptionValue[VertexFilter], 
       (False|None) :> True],
     fFiltFn = Replace[
       OptionValue[FaceFilter], 
       (False|None) :> True]},
    sym /: SurfaceQ[sym] = True;
    sym /: SurfaceName[sym] = sym;
    sym /: FaceFilter[sym] = OptionValue[FaceFilter];
    sym /: VertexFilter[sym] = OptionValue[VertexFilter];
    With[
      {idcs = If[uFiltFn === True,
         All,
         Select[Range[Length[data]], TrueQ[uFiltFn[data[[#]], field[[#]]]]&]]},
      sym /: VertexIndexDispatch[sym] := TagSet[
        sym,
        VertexIndexDispatch[sym],
        Dispatch[
          Append[
            MapThread[Rule, {idcs, Range[Length@idcs]}],
            _ -> Undefined]]];
      sym /: VertexList[sym] = data[[idcs]];
      sym /: VertexList[sym, toStyle_] := TagSet[
        sym, 
        VertexList[sym, toStyle],
        ConvertCoordinates[data[[idcs]], Cartesian -> toStyle]];
      sym /: Field[sym] := TagSet[
        sym, 
        Field[sym],
        field[[idcs]]];
      sym /: Normal[sym] := TagSet[
        sym,
        Normal[sym],
        MapThread[Rule, {data[[idcs]], field[[idcs]]}]];
      sym /: Normal[sym, toStyle_] := TagSet[
        sym,
        Normal[sym, toStyle],
        MapThread[Rule, {ConvertCoordinates[data[[idcs]], Cartesian -> toStyle], field[[idcs]]}]];
      sym /: FaceList[sym] := TagSet[
        sym, 
        FaceList[sym],
        With[
          {faceList = Which[
             faces =!= None, faces,
             polygons =!= None, Replace[
               polygons,
               Dispatch[MapThread[Rule, {data, Range[Length[data]]}]],
               {2}],
             True, None]},
          With[
            {faceListUncor = Which[
               faceList === None, None,
               idcs === All && fFiltFn === True, faceList,
               idcs === All, Select[faceList, fFiltFn],
               True, With[
                 {disp = Dispatch[Append[Map[(# -> True)&, idcs], _ -> False]],
                  filt = If[TrueQ[fFiltFn], (True&), fFiltFn]},
                 Select[faceList, If[filt[#], And@@Replace[#, disp, {1}], False]&]]]},
            If[idcs === All || Length[idcs] == Length[data],
              faceListUncor,
              Replace[faceListUncor, VertexIndexDispatch[sym], {2}]]]]];
      sym /: Polygons[sym] := TagSet[
        sym,
        Polygons[sym],
        data[[#]]& /@ FaceList[sym]];
      sym /: Polygons[sym, toStyle_] := TagSet[
        sym, 
        Polygons[sym, toStyle],
        Map[ConvertCoordinates[#, Cartesian -> toStyle]&, Polygons[sym]]];
      sym]]];
Surface[data: {Rule[{_,_,_},_]..}, opts:OptionsPattern[]] := Surface[
  data[[All,1]],
  opts,
  Field -> data[[All,2]]];
Surface[data: {{_,_}..}, opts:OptionsPattern[]] := With[
  {sc = OptionValue[SphericalCoordinateStyle]},
  Surface[
    ConvertCoordinates[data, sc -> Cartesian],
    opts]];
Surface[data: {Rule[{_,_},_]..}, opts:OptionsPattern[]] := Surface[
  data[[All, 1]],
  opts,
  Field -> data[[All, 2]]];
Surface[surf_?SurfaceQ, opts:OptionsPattern[]] := Surface[
  VertexList[surf],
  opts,
  Field -> Field[surf],
  FaceList -> FaceList[surf],
  VertexFilter -> VertexFilter[surf],
  FaceFilter -> FaceFilter[surf]];

(* #SurfaceQ **************************************************************************************)
SurfaceQ[x_] = False;

(* #SurfaceFromVTK ********************************************************************************)
Options[SurfaceFromVTK] = Options[Surface];
SurfaceFromVTK[filename_String, opts:OptionsPattern[]] := Surface[
  ReadVTK[filename],
  opts];

(* #MergeSurfaces *********************************************************************************)
Options[MergeSurfaces] = {
  Function -> Function[#2],
  VertexFilter -> None,
  FaceFilter -> None};
MergeSurfaces[surfs:{_?SurfaceQ, _?SurfaceQ..}, OptionsPattern[]] := With[
  {n = Length[surfs],
   ideal = Range[Length[surfs]],
   mergeFn = OptionValue[Function],
   uFilt = Replace[OptionValue[VertexFilter], None|False -> (True&)],
   fFilt = Replace[OptionValue[FaceFilter], None|False -> (True&)]},
  Surface[
    Map[
      Function[
        If[Length[Union[#[[All,3]]]] != Length[#],
          (Message[MergeSurfaces::dup]; #[[1,1]] -> Indeterminate),
          #[[1,1]] -> mergeFn[
            #[[1,1]],
            Replace[Range[n], Append[Map[(#[[3]] -> #[[2]])&, #], (_ -> None)], {1}]]]],
      GatherBy[
        Select[
          Transpose[
            {Flatten[VertexList /@ surfs], 
             Flatten[Field /@ surfs],
             Flatten[Table[k, {k,1,Length[surfs]}, {Length[VertexList[surf[[k]]]]}]]}],
          uFilt[#[[1]]]&],
        First]],
    FaceList -> Select[FaceList[surfs[[1]]], fFilt],
    VertexFilter -> OptionValue[VertexFilter],
    FaceFilter -> OptionValue[FaceFilter]]];

(* #SurfaceProjection *****************************************************************************)
Options[SurfaceProjection] = {
  Center -> Automatic,
  ProjectionShear -> Automatic,
  Radius -> Automatic,
  OrientPoint -> Automatic,
  Duplicate -> None,
  VertexFilter -> None,
  FaceFilter -> None,
  Scale -> 1.0};
SurfaceProjection[surf_, OptionsPattern[]] := Catch[
  With[
    {dflt = Replace[
       OptionValue[Duplicate],
       {(None|False) -> {
          OrientPoint -> ({0,1,0} -> 0),
          Center -> {0,0,1},
          Radius -> (Pi / 3.0),
          ProjectionShear -> {{1,0},{0,1}},
          VertexFilter -> None,
          FaceFilter -> None},
        map_?MapQ :> {
          OrientPoint -> OrientPoint[map],
          Center -> Center[map],
          Radius -> Radius[map],
          ProjectionShear -> ProjectionShear[map],
          VertexFilter -> VertexFilter[map],
          FaceFilter -> FaceFilter[map]},
        _ :> (
          Message[SurfaceProjection::badarg, Duplicate, "must be None or a projection map"];
          Throw[$Failed])}]},
    With[
      {mu = Replace[
       OptionValue[Center],
         {Automatic :> (Center /. dflt),
          Except[{_?NumberQ,_?NumberQ,_?NumberQ}] :> (
            Message[SurfaceProjection::badarg, Center, "must be of the form {x, y, z}"];
            Throw[$Failed])}],
       shear = Replace[
         OptionValue[ProjectionShear],
         {Automatic :> (ProjectionShear /. dflt),
          Except[{{(1|1.), _?AtomQ}, {_?AtomQ, (1|1.)}}] :> (
            Message[
              SurfaceProjection::badarg,
              ProjectionShear,
              "must be of the form {{1, Sx}, {Sy, 1}}"];
            Throw[$Failed])}],
       rad = Replace[
         N[OptionValue[Radius]],
         {Automatic :> (Radius /. dflt),
          x_ /; Or[Not[NumberQ[x]], x <= 0] :> (
            Message[SurfaceProjection::badarg, Radius, "must be a positive number"];
            Throw[$Failed])}],
       orient = Replace[
         OptionValue[OrientPoint],
         {Automatic :> (OrientPoint /. dflt),
          x_Rule :> x,
          x_List :> (x -> 0),
          _ :> (
            Message[SurfaceProjection::badarg, OrientPoint, "must be a rule or a point"];
            Throw[$Failed])}],
       uFilt = Replace[OptionValue[VertexFilter], None|False|Automatic -> (True&)],
       fFilt = Replace[OptionValue[FaceFilter], None|False|Automatic -> (True&)],
       scale = Replace[
         OptionValue[Scale],
         {Automatic -> 1.0,
          None -> None,
          r_?NumericQ /; r > 0 :> r,
          _ :> (
            Message[SurfaceProjection::badarg, Scale, "must be a real number > 0 or None"];
            Throw[$Failed])}],
       sym = (SetAttributes[Evaluate[#], Temporary]; #)&@Unique["map"],
       scFn = Function[CartesianToSpherical[#, SphericalCoordinateStyle->{Longitude, Latitude}]],
       iscFn = Function[SphericalToCartesian[#, SphericalCoordinateStyle->{Longitude, Latitude}]]},
      sym /: SurfaceRotation[sym] = RotationMatrix[{mu, {1,0,0}}];
      With[
        {op = scFn[SurfaceRotation[sym] . orient[[1]]]},
        sym /: OrientMatrix[sym] = If[op == {0,0},
          {{1,0},{0,1}},
          RotationMatrix[{op, {Cos[orient[[2]]],Sin[orient[[2]]]}}]];
        With[
          {opT = Transpose[shear . OrientMatrix[sym]],
           rotT = Transpose[SurfaceRotation[sym]],
           iopT = Transpose[Inverse[shear . OrientMatrix[sym]]],
           irotT = Transpose[Inverse[SurfaceRotation[sym]]]},
          sym /: DomainIndices[sym] = With[
            {cosrad = Cos[rad],
             normmu = Normalize[mu]},
            Flatten[
              Position[
                Normal[surf],
                Rule[xyz_,_] /; Dot[Normalize[xyz], normmu] >= cosrad]]];
          sym /: Domain[sym] = Part[Normal[surf], DomainIndices[sym]];
          With[
            {scaleMult = If[scale === None,
               1.0,
               With[
                 {tmp = Dot[scFn[Dot[#, rotT]], opT]& /@ Domain[sym][[All,1]],
                  idxDisp = Dispatch[
                    MapThread[Rule, {DomainIndices[sym], Range[Length[Domain[sym]]]}]],
                  idxSet = Dispatch[Append[Map[#->True&, DomainIndices[sym]], _->False]]},
                 With[
                   {edges = Transpose[
                      tmp[[#]]& /@ ReplaceAll[
                        Select[
                          EdgeList[surf],
                          And[#[[1]] /. idxSet, #[[2]] /. idxSet]&],
                        idxDisp]]},
                   scale / Mean[Norm /@ (edges[[1]] - edges[[2]])]]]]},
            With[
              {fn = Function[scaleMult*Dot[scFn[Dot[#, rotT]], opT]],
               ifn = Function[Dot[iscFn[Dot[#/scaleMult, iopT]], irotT]]},
              sym /: MapQ[sym] = True;
              sym /: MapName[sym] = sym;
              sym /: VertexFilter[sym] = OptionValue[VertexFilter];
              sym /: FaceFilter[sym] = fFilt;
              sym /: ProjectedSurface[sym] = surf;
              sym /: Center[sym] = mu;
              sym /: ProjectionShear[sym] = shear;
              sym /: Radius[sym] = rad;
              sym /: OrientPoint[sym] = orient;
              sym /: SurfaceRotation[sym] = RotationMatrix[{mu, {1, 0, 0}}];
              sym /: ProjectionDispatch[sym] = Dispatch[
                {pat:{Rule[{_,_,_},_]..} :> MapThread[Rule, {fn[pat[[All,1]]], pat[[All,2]]}],
                 pat:{{_,_,_}..} :> fn[pat],
                 Rule[a:{_,_,_},b_] :> Rule[First[fn[{a}]], b],
                 a:{_?AtomQ,_?AtomQ,_?AtomQ} :> First[fn[{a}]],
                 s_?SurfaceQ :> MapThread[Rule, {fn[VertexList[s]], Field[s]}],
                 _ :> (Message[SurfaceProjection::badfmt]; Indeterminate)}];
              sym /: InverseProjectionDispatch[sym] = Dispatch[
                {pat:{Rule[{_,_},_]..} :> MapThread[Rule, {ifn[pat[[All,1]]], pat[[All,2]]}],
                 pat:{{_,_}..} :> ifn[pat],
                 Rule[a:{_,_},b_] :> Rule[First[ifn[{a}]], b],
                 a:{_?AtomQ,_?AtomQ} :> First[ifn[{a}]],
                 s_?SurfaceQ :> MapThread[Rule, {ifn[VertexList[s]], Field[s]}],
                 _ :> (Message[SurfaceProjection::badfmt]; Indeterminate)}];
              sym /: ProjectionTransform[sym] = Function[Replace[#, ProjectionDispatch[sym]]];
              sym /: InverseProjectionTransform[sym] = Function[
                Replace[#, InverseProjectionDispatch[sym]]];
              sym /: Normal[sym] := TagSet[
                sym,
                Normal[sym],
                Replace[
                  Select[Domain[sym], TrueQ[uFilt[#[[1]], #[[2]]]]&],
                  ProjectionDispatch[sym],
                  {1}]];
              sym /: VertexList[sym] := TagSet[sym, VertexList[sym], Normal[sym][[All,1]]];
              sym /: VertexIndexDispatch[sym] := TagSet[
                sym,
                VertexIndexDispatch[sym],
                Dispatch[
                  Append[
                    MapThread[
                      Rule, 
                      {DomainIndices[sym], Range[Length@DomainIndices[sym]]}],
                    _ -> Undefined]]];
              sym /: Field[sym] := TagSet[sym, Field[sym], Normal[sym][[All,2]]];
              sym /: FaceList[sym] := TagSet[
                sym,
                FaceList[sym],
                With[
                  {dom = Dispatch[
                     Append[
                       Map[(# -> True)&, DomainIndices[sym]],
                       _Integer -> False]],
                   tr = VertexIndexDispatch[sym]},
                  Replace[
                    Select[
                      FaceList[surf],
                      If[fFilt[#], Apply[And, # /. dom], False]&],
                    tr,
                    {2}]]];
              sym /: Polygons[sym] := TagSet[
                sym,
                Polygons[sym],
                Map[VertexList[sym][[#]]&, FaceList[sym]]];
              sym]]]]]]];

(* #StructuredRescale *****************************************************************************)
(* This is used be MapPlot and SurfacePlot when normalizing oddly structured data sets in the color
   function scaling process *)
StructuredRescale[dat_List] := If[ArrayDepth[dat] == 1,
   Rescale[dat],
   Transpose[
     Map[
      StructuredRescale,
      Transpose[dat]]]];

(* #MapPlot ***************************************************************************************)
Options[MapPlot] = Join[
  Options[Graphics],
  {ColorFunction -> Automatic,
   ColorFunctionScaling -> False}];
MapPlot[map_?MapQ, opts:OptionsPattern[]] := Graphics[
  With[
    {V = Normal[VertexList[map]],
     F = Normal[FaceList[map]],
     Z = Normal[Field[map]],
     cfn = Replace[
       OptionValue[ColorFunction],
       Automatic -> Function[Blend[{Blue,Cyan,Gray,Yellow,Red},#]]],
     cfnsc = Replace[
       OptionValue[ColorFunctionScaling],
       {None|False -> Identity,
        Ordering -> Function[
          With[
            {ord = Ordering[#]},
            (Normal[SparseArray[Thread[ord -> Range[Length[ord]]]]] - 1) / (Length[#] - 1)]],
        True|Automatic -> StructuredRescale}]},
    With[
      {ZZ = cfnsc[Z]},
      {EdgeForm[None],
       If[!ListQ[F],
         MapThread[
           Function[{cfn[#1], Point[#2]}],
           {ZZ, V}],
         Map[
           Function[Polygon[V[[#]], VertexColors -> Map[cfn, ZZ[[#]]]]],
           F]]}]],
  Sequence@@FilterRules[
    {opts},
    Except[ColorFunction|ColorFunctionScaling]]];

(* #MapMeshPlot ***********************************************************************************)
Options[MapMeshPlot] = Join[
  Options[Graphics],
  {ColorFunction -> Automatic,
   ColorFunctionScaling -> False}];
MapMeshPlot[map_?MapQ, opts:OptionsPattern[]] := Graphics[
  With[
    {V = Normal[VertexList[map]],
     F = Normal[EdgeList[map]],
     Z = Normal[Field[map]],
     cfn = Replace[
       OptionValue[ColorFunction],
       Automatic -> Function[Blend[{Blue,Cyan,Gray,Yellow,Red},#]]],
     cfnsc = Replace[
       OptionValue[ColorFunctionScaling],
       {None|False -> Identity,
        Ordering -> Function[
          With[
            {ord = Ordering[#]},
            (Normal[SparseArray[Thread[ord -> Range[Length[ord]]]]] - 1) / (Length[#] - 1)]],
        True|Automatic -> StructuredRescale}]},
    With[
      {ZZ = cfnsc[Z]},
      {If[ListQ[F],
         Map[
           Function[Line[V[[#]], VertexColors -> cfn /@ Z[[#]]]],
           F],
         (Message[MapMeshPlot::badarg, "map must have edge data"];
          {})]}]],
  Sequence@@FilterRules[
    {opts},
    Except[ColorFunction|ColorFunctionScaling]]];

(* #MapQ ******************************************************************************************)
MapQ[x_] = False;

(* #WithVertexFilter ******************************************************************************)
WithVertexFilter[surf_?SurfaceQ, filt_] := Surface[surf, VertexFilter -> filt];
WithVertexFilter[map_?MapQ, filt_] := With[
  {newFilt = Replace[
     VertexFilter[map],
     {None -> filt,
      f_ :> (And[f[#1,#2], filt[#1,#2]]&)}]},
  SurfaceProjection[Domain[map], Duplicate -> map, VertexFilter -> newFilt]];

(* #WithFaceFilter ********************************************************************************)
WithFaceFilter[surf_?SurfaceQ, filt_] := Surface[surf, FaceFilter -> filt];
WithFaceFilter[map_?MapQ, filt_] := With[
  {newFilt = Replace[
     FaceFilter[map],
     {None -> filt,
      f_ :> (And[f[#], filt[#]]&)}]},
  SurfaceProjection[Domain[map], Duplicate -> map, FaceFilter -> newFilt]];

(* #MergePolygons *********************************************************************************)
(* Used by SurfacePlot below *)
MergePolygons[polygons_List, X_List, categories_List] := With[
  {nomerge = Select[
     polygons,
     Function[{pg}, Length[Union[categories[[pg]]]] > 1]],
   canmerge = GatherBy[
     Select[
       polygons,
       Function[{pg}, Length[Union[categories[[pg]]]] == 1]],
     categories[[#[[1]]]]&]},
  Join[
    Map[
     Function[{pgs},
       {categories[[pgs[[1, 1]]]],
        Polygon[X[[#]] & /@ pgs]}],
     canmerge],
    Polygon[X[[#]], VertexColors -> categories[[#]]] & /@ nomerge]];

(* #SurfacePlot ***********************************************************************************)
Options[SurfacePlot] = Join[
  FilterRules[
    Options[Graphics3D],
    Except[ColorFunction|ColorFunctionScaling|Lighting|Boxed]],
  {ColorFunction -> Automatic,
   ColorFunctionScaling -> False,
   Lighting -> "Neutral",
   Boxed -> False}];
SurfacePlot[surf_?SurfaceQ, opts:OptionsPattern[]] := Graphics3D[
  With[
    {V = Normal[VertexList[surf]],
     F = Normal[FaceList[surf]],
     Z = Normal[Field[surf]],
     cfn = Replace[
       OptionValue[ColorFunction],
       Automatic -> Function[
         Blend[
           {Blue, Darker[Cyan, 1/6], Darker[Green], Darker[Yellow, 1/6], Red},
           #]]]},
    With[
      {ZZ = Replace[
         OptionValue[ColorFunctionScaling],
         {None|False :> Z,
          Ordering :> With[
            {ord = Ordering[Z]},
            (Normal[SparseArray[Thread[ord -> Range[Length[ord]]]]] - 1) / (Length[Z] - 1)],
          True|Automatic :> StructuredRescale[Z]}]},
      {EdgeForm[None],
       If[!ListQ[F],
         MapThread[
           Function[{cfn[#1], Point[#2]}],
           {ZZ, V}],
         MergePolygons[F, V, cfn /@ ZZ]]}]],
  Sequence@@FilterRules[
    Join[{opts}, Options[SurfacePlot]],
    Except[ColorFunction|ColorFunctionScaling]]];

(* #WithVertexList **********************************************************************************)
Options[WithVertexList] = {SphericalCoordinateStyle -> None};
WithVertexList[surf_?SurfaceQ, X_List] /; Length[VertexList[surf]] != Length[X] := (
  Message[WithVertexList::incompat];
  $Failed);
WithVertexList[surf_?SurfaceQ, X_List /; MatchQ[Dimensions[X], {_,3}]] := Surface[
    X, 
    Field -> Field[surf],
    FaceList -> FaceList[surf],
    FaceFilter -> FaceFilter[surf],
    VertexFilter -> VertexFilter[surf]];
WithVertexList[surf_?MapQ, X_List /; MatchQ[Dimensions[X], {_,2}]] := SurfaceProjection[
    X, 
    Field -> Field[surf],
    FaceList -> FaceList[surf],
    FaceFilter -> FaceFilter[surf],
    VertexFilter -> VertexFilter[surf],
    SphericalCoordinateStyle -> OptionValue[SphericalCoordinateStyle]];

(* #SurfaceResample *******************************************************************************)
Options[SurfaceResample] = Prepend[
  FilterRules[Options[Surface], Except[Field]],
  Method -> Nearest];
SurfaceResample[a_?SurfaceQ, b_?SurfaceQ, opts:OptionsPattern[]] := Catch[
  If[Length[Field[a]] == Length[Field[b]] && VertexList[a] == VertexList[b],
    WithField[b, Field[a]],
    Surface[
      VertexList[b],
      Sequence@@FilterRules[{opts}, Except[Method]],
      FaceList -> FaceList[b],
      Field -> With[
        {method = OptionValue[Method]},
        Switch[
          method,
          Nearest, With[
            {nearest = Nearest[VertexList[a] -> Automatic],
             F = Field[a]},
            Map[Function[F[[First[nearest[#]]]]], VertexList[b]]],
          (Interpolation|{Interpolation, ___}), With[
            {interp = Interpolation[
               MapThread[List, {VertexList[a], Field[a]}],
               Sequence@@If[ListQ[method], Rest[methods], {}]]},
            Map[interp, VertexList[b]]],
          (ListInterpolation|{ListInterpolation, ___}), With[
            {interp = ListInterpolation[
               MapThread[List, {VertexList[a], Field[a]}],
               Sequence@@If[ListQ[method], Rest[methods], {}]]},
            Map[interp, VertexList[b]]],
          _, (
            Message[
              SurfResample::badarg,
              Method,
              "must be Nearest/Interpolation/ListInterpolation"];
            Throw[$Failed])]]]]];

(* #MapConvexHull *********************************************************************************)
MapConvexHull[map_?MapQ] := If[MapName[map] =!= map,
  MapHull[MapName[map]],
  With[
    {res = Check[ConvexHull[VertexList[map]], $Failed]},
    If[res === $Failed, res, (map /: MapConvexHull[map] = res)]]];

(* #MapHull ***************************************************************************************)
MapHull[map_?MapQ] := If[MapName[map] =!= map,
  MapHull[MapName[map]],
  With[
    {res = Check[
       With[
         {hull = Cases[
            Tally[Flatten[Subsets[#,{2}]& /@ FaceList[map], 1]],
            {edge_List, 1} :> edge,
            {1}],
          X = VertexList[map]},
         (* any edge that appears only once is on the hull *)
         SortBy[
           Union[Flatten[hull]],
           Pi + ArcTan[X[[#, 1]], X[[#, 2]]]]],
       $Failed]},
    If[res === $Failed, res, (map /: MapHull[map] = res)]]];

(* #EdgeList **************************************************************************************)
EdgeList[surf_?SurfaceQ] := If[SurfaceName[surf] =!= surf,
  EdgeList[SurfaceName[surf]],
  With[
    {res = Check[
       Union@Map[
         Sort,
         Flatten[Map[Subsets[#,{2}]&, FaceList[surf]], 1]],
       $Failed]},
    If[res === $Failed, res, (surf /: EdgeList[surf] = res)]]];
EdgeList[map_?MapQ] := If[MapName[map] =!= map,
  EdgeList[MapName[map]],
  With[
    {res = Union@Map[
       Sort,
       Flatten[Map[Subsets[#,{2}]&, FaceList[map]], 1]]},
    If[res === $Failed, res, (map /: EdgeList[map] = res)]]];

(* #EdgeLengths ***********************************************************************************)
EdgeLengths[surf_ /; SurfaceQ[surf] || MapQ[surf], X_] := With[
  {EL = Transpose[EdgeList[surf]]},
  With[
    {X1 = Transpose[X[[EL[[1]]]]],
     X2 = Transpose[X[[EL[[2]]]]]},
    Sqrt[
      Plus[
        (X1[[1]] - X2[[1]])^2,
        (X1[[2]] - X2[[2]])^2,
        If[MapQ[surf], 0, (X1[[3]] - X2[[3]])^2]]]]];
EdgeLengths[surf_] := Which[
  SurfaceQ[surf] && surf =!= SurfaceName[surf], EdgeLengths[SurfaceName@surf],
  MapQ[surf] && surf =!= MapName[surf], EdgeLengths[MapName@surf],
  True, With[
    {res = Check[EdgeLengths[surf, VertexList[surf]], $Failed]},
    If[res === $Failed || !ListQ[res], $Failed, (surf /: EdgeLengths[surf] = res)]]];

(* #EdgesEnergy ***********************************************************************************)
EdgesEnergy[surf_ /; SurfaceQ[surf] || MapQ[surf], X_] := Total[
  Log[EdgeLengths[surf] / EdgeLengths[surf, X]]^2];
EdgesEnergy[surf_ /; SurfaceQ[surf] || MapQ[surf], X_, idcs_List] := Total[
  Subtract[
    Flatten[NeighborhoodEdgeLengths[surf][[idcs]]],
    Flatten[NeighborhoodEdgeLengths[surf, X, idcs]]]^2];
EdgesEnergy[surf_ /; SurfaceQ[surf] || MapQ[surf]] := 0;

(*
EdgesEnergy[surf_ /; SurfaceQ[surf] || MapQ[surf], X_] := Total[
  (EdgeLengths[surf] - EdgeLengths[surf, X])^2];
EdgesEnergy[surf_ /; SurfaceQ[surf] || MapQ[surf], X_, idcs_List] := Total[
  Subtract[
    Flatten[NeighborhoodEdgeLengths[surf][[idcs]]],
    Flatten[NeighborhoodEdgeLengths[surf, X, idcs]]]^2];
EdgesEnergy[surf_ /; SurfaceQ[surf] || MapQ[surf]] := 0;
*)

(* #EdgesGradient *********************************************************************************)
EdgesGradientCompiled = Compile[{{x0, _Real, 1}, {x, _Real, 2}, {d0, _Real, 1}},
  With[
    {dx = MapThread[Subtract, {Transpose[x], x0}]},
    With[
      {norms = Sqrt[Total[dx^2]]},
      Dot[dx, -2.0 * Log[norms / d0] / norms^2]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
EdgesGradient[surf_ /; SurfaceQ[surf] || MapQ[surf], X_] := MapThread[
  EdgesGradientCompiled,
  {X, X[[#]]& /@ NeighborhoodList[surf], NeighborhoodEdgeLengths[surf]}];
EdgesGradient[surf_ /; SurfaceQ[surf] || MapQ[surf], X_, idcs_List] := MapThread[
  EdgesGradientCompiled,
  {X[[idcs]],
   X[[#]]& /@ Part[NeighborhoodList[surf], idcs],
   Part[NeighborhoodEdgeLengths[surf], idcs]}];
EdgesGradient[surf_ /; SurfaceQ[surf] || MapQ[surf]] := ConstantArray[
  0,
  Dimensions[VertexList[surf]]];

(*
EdgesGradientCompiled = Compile[{{x0, _Real, 1}, {x, _Real, 2}, {d0, _Real, 1}},
  With[
    {dx = MapThread[Subtract, {Transpose[x], x0}]},
    With[
      {norms = Sqrt[Total[dx^2]]},
      Dot[dx, 2.0 * (d0 - norms) / norms]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
EdgesGradient[surf_ /; SurfaceQ[surf] || MapQ[surf], X_] := MapThread[
  If[Length[#2]<1, Table[0, {Length[#1]}], EdgesGradientCompiled[#1,#2,#3]]&,
  {X, X[[#]]& /@ NeighborhoodList[surf], NeighborhoodEdgeLengths[surf]}];
EdgesGradient[surf_ /; SurfaceQ[surf] || MapQ[surf], X_, idcs_List] := MapThread[
  EdgesGradientCompiled,
  {X[[idcs]],
   X[[#]]& /@ Part[NeighborhoodList[surf], idcs],
   Part[NeighborhoodEdgeLengths[surf], idcs]}];
EdgesGradient[surf_ /; SurfaceQ[surf] || MapQ[surf]] := ConstantArray[
  0,
  Dimensions[VertexList[surf]]];
*)     

(* #FaceAngles ************************************************************************************)
FaceAngles[surf_?SurfaceQ[surf], X_] := With[
  {F = Transpose[FaceList[surf]]},
  With[
    {Xf = X[[#]]& /@ F},
    With[
      {dX = Transpose /@ {
         (Xf[[2]] - Xf[[1]]),
         (Xf[[3]] - Xf[[2]]),
         (Xf[[1]] - Xf[[3]])}},
      With[
        {normed = Map[
           With[{lens = Sqrt[Total[#^2]]}, (#/lens)& /@ #]&,
           dX]},
        Transpose[
          ArcCos[
            {Total[normed[[1]] * -normed[[3]]],
             Total[normed[[2]] * -normed[[1]]],
             Total[normed[[3]] * -normed[[2]]]}]]]]]];
FaceAngles[surf_?MapQ[surf], X_] := With[
  {F = Transpose[FaceList[surf]]},
  With[
    {Xf = X[[#]]& /@ F},
    With[
      {dX = Transpose /@ {
         (Xf[[2]] - Xf[[1]]),
         (Xf[[3]] - Xf[[2]])}},
      With[
        {normed = Map[
           With[{lens = Sqrt[Total[#^2]]}, (#/lens)& /@ #]&,
           dX]},
        Transpose[
          ArcCos[
            {Total[normed[[1]] * -normed[[3]]],
             Total[normed[[2]] * -normed[[1]]]}]]]]]];
FaceAngles[surf_] := Which[
  SurfaceQ[surf] && surf =!= SurfaceName[surf], FaceAngles[SurfaceName@surf],
  MapQ[surf] && surf =!= MapName[surf], FaceAngles[MapName@surf],
  True, With[
    {res = Check[FaceAngles[surf, VertexList[surf]], $Failed]},
    If[res === $Failed || !ListQ[res], $Failed, (surf /: FaceAngles[surf] = res)]]];

(* #NeighborhoodList ******************************************************************************)
NeighborhoodList[surf_?SurfaceQ] := If[SurfaceName[surf] =!= surf,
  NeighborhoodList[SurfaceName[surf]],
  With[
    {res = Check[
       With[
         {V = VertexList[surf]},
         Last[
           Reap[
             Scan[
               Function[{face},
                 Scan[
                   Function[{pair},
                     Sow[pair[[1]], pair[[2]]];
                     Sow[pair[[2]], pair[[1]]]],
                   Subsets[face, {2}]]],
               FaceList[surf]],
             Range[Length[V]],
             Function[{id, neighbors},
               Apply[
                 Sequence,
                 With[
                   {neis = Union[neighbors]},
                   With[
                     {U = Dot[
                       V[[neis]], 
                       Transpose[RotationMatrix[{V[[id]], {0,0,1}}]]
                      ][[All, 1;;2]]},
                     SortBy[Thread[neis -> U], ArcTan[#[[2,1]], #[[2,2]]]&][[All,1]]]]]]]]],
       $Failed]},
    If[res === $Failed, $Failed, (surf /: NeighborhoodList[surf] = res)]]];
NeighborhoodList[map_?MapQ] := If[MapName[map] =!= map,
  NeighborhoodList[MapName[map]],
  With[
    {res = Check[
       With[
         {V = VertexList[map]},
         Last[
           Reap[
             Scan[
               Function[{face},
                 Scan[
                   Function[{pair},
                     Sow[pair[[1]], pair[[2]]];
                     Sow[pair[[2]], pair[[1]]]],
                   Subsets[face, {2}]]],
               FaceList[map]],
             Range[Length[V]],
             Function[{id, neighbors},
               Apply[
                 Sequence,
                 With[
                   {neis = Union[neighbors]},
                   With[
                     {U = V[[neis]]},
                     SortBy[
                       Thread[neis -> U], 
                       ArcTan[#[[2,1]] - V[[id,1]], #[[2,2]] - V[[id,2]]]&
                      ][[All,1]]]]]]]]],
      $Failed]},
    If[res === $Failed, res, (map /: NeighborhoodList[map] = res)]]];

(* #FacesIndex ************************************************************************************)
FacesIndex[surf_ /; SurfaceQ[surf] || MapQ[surf]] := Which[
  MapQ[surf] && MapName[surf] =!= surf, FacesIndex[MapName[surf]],
  SurfaceQ[surf] && SurfaceName[surf] =!= surf, FacesIndex[SurfaceName[surf]],
  True, With[
    {res = Check[
       Last@Reap[
         MapIndexed[Sow[#2,#1]&, FaceList[surf], {2}],
         Range[Length@VertexList[surf]],
         Function[Sequence@@#2]],
       $Failed]},
    If[res === $Failed || !ListQ[res], $Failed, (surf /: FacesIndex[surf] = res)]]];

(* #NeighborhoodAngles ****************************************************************************)
NeighborhoodAngleCompiled = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}},
  With[
    {x = Transpose[xnei]},
    With[
      {dx = Table[x[[i]]-x0[[i]], {i,1,Length[x0]}]},
      With[
        {norms = Sqrt[Total[dx^2]]},
        With[
          {normed = dx/Table[norms, {Length[x0]}]},
          ArcCos[Total[normed*(RotateLeft /@ normed)]]]]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
Protect[NeighborhoodAngleCompiled];
NeighborhoodAngles[surf_ /; SurfaceQ[surf] || MapQ[surf], X_] := MapThread[
  If[Length[#2]<2, {}, NeighborhoodAngleCompiled[#1,#2]]&,
  {X, X[[#]]& /@ NeighborhoodList[surf]}];
NeighborhoodAngles[surf_ /; SurfaceQ[surf] || MapQ[surf], X_, idcs_List] := MapThread[
  If[Length[#2]<2, {}, NeighborhoodAngleCompiled[#1,#2]]&,
  {X[[idcs]], X[[#]]& /@ Part[NeighborhoodList[surf], idcs]}];
NeighborhoodAngles[surf_] := Which[
  MapQ[surf] && MapName[surf] =!= surf, NeighborhoodAngles[MapName[surf]],
  SurfaceQ[surf] && SurfaceName[surf] =!= surf, NeighborhoodAngles[SurfaceName[surf]],
  True, With[
    {res = Check[NeighborhoodAngles[surf, VertexList[surf]], $Failed]},
    If[res === $Failed || !ListQ[res], $Failed, (surf /: NeighborhoodAngles[surf] = res)]]];

(* #NeighborhoodBisectors *************************************************************************)
NeighborhoodBisectorsCompiled = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}},
  With[
    {x = Transpose[xnei]},
    With[
      {dx = Table[x[[i]]-x0[[i]], {i,1,Length[x0]}]},
      With[
        {norms = Sqrt[Total[dx^2]]},
        With[
          {normed = dx/Table[norms, {Length[x0]}]},
          With[
            {means = 0.5*(normed + RotateLeft /@ normed)},
            With[
              {mnorms = Sqrt[Total[means^2]]},
              Transpose[means/Table[mnorms, {Length[x0]}]]]]]]]],
   RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
   Parallelization -> True];
Protect[NeighborhoodBisectorsCompiled];
NeighborhoodBisectors[surf_ /; SurfaceQ[surf] || MapQ[surf], X_] := MapThread[
  NeighborhoodBisectorsCompiled,
  {X, X[[#]] & /@ NeighborhoodList[surf]}];
NeighborhoodBisectors[surf_] := Which[
  MapQ[surf] && MapName[surf] =!= surf, NeighborhoodBisectors[MapName[surf]],
  SurfaceQ[surf] && SurfaceName[surf] =!= surf, NeighborhoodBisectors[SurfaceName[surf]],
  True, With[
    {res = Check[NeighborhoodBisectors[surf, VertexList[surf]], $Failed]},
    If[res === $Failed || !ListQ[res], $Failed, (surf /: NeighborhoodBisectors[surf] = res)]]];

(* #NeighborhoodEdgeLengths ***********************************************************************)
NeighborhoodEdgeLengthsCompiled = Compile[{{x0, _Real, 1}, {x, _Real, 2}},
  Sqrt[Total[MapThread[Subtract, {x0, Transpose[x]}]^2]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
NeighborhoodEdgeLengths[surf_ /; SurfaceQ[surf] || MapQ[surf], X_] := MapThread[
  If[Length[#2]<1, {}, NeighborhoodEdgeLengthsCompiled[#1,#2]]&,
  {X, X[[#]]& /@ NeighborhoodList[surf]}];
NeighborhoodEdgeLengths[surf_ /; SurfaceQ[surf] || MapQ[surf], X_, idcs_List] := MapThread[
  If[Length[#2]<1, {}, NeighborhoodEdgeLengthsCompiled[#1,#2]]&,
  {X[[idcs]], X[[#]]& /@ Part[NeighborhoodList[surf], idcs]}];
NeighborhoodEdgeLengths[surf_ /; SurfaceQ[surf] || MapQ[surf]] := Which[
  MapQ[surf] && MapName[surf] =!= surf, NeighborhoodEdgeLengths[MapName[surf]],
  SurfaceQ[surf] && SurfaceName[surf] =!= surf, NeighborhoodEdgeLengths[SurfaceName[surf]],
  True, With[
    {nel = Check[NeighborhoodEdgeLengths[surf, VertexList[surf]], $Failed]},
    If[nel === $Failed || !ListQ[nel],
      $Failed,
      (surf /: NeighborhoodEdgeLengths[surf] = nel)]]];

(* #AnglesGradient ********************************************************************************)
(* This creates a function that can calculate the gradient of a (2d) angle in terms of x, and y *)
AnglesGradientCompiled2D = Block[{x0, y0, x1, y1, x2, y2},
  Compile @@ List[
    Map[{#, _Real}&, {x0, y0, x1, y1, x2, y2}],
    Map[
      Function[
        Simplify[
          D[Simplify[
              VectorAngle[{x1-x0, y1-y0}, {x2-x0, y2-y0}], 
              Element[{x0, y0, x1, y1, x2, y2}, Reals]],
            #],
          Element[{x0, y0, x1, y1, x2, y2}, Reals]]],
      {x0, y0}],
    RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
    RuntimeAttributes -> {Listable},
    Parallelization -> True]];
AnglesGradientCompiled3D = Block[{x0, y0, z0, x1, y1, z1, x2, y2, z2},
  Compile @@ List[
    Map[{#, _Real}&, {x0, y0, z0, x1, y1, z1, x2, y2, z2}],
    Map[
      Function[
        Simplify[
          D[Simplify[
              VectorAngle[{x1-x0, y1-y0, z1-z0}, {x2-x0, y2-y0, z2-z0}], 
              Element[{x0, y0, z0, x1, y1, z1, x2, y2, z2}, Reals]],
            #],
          Element[{x0, y0, z0, x1, y1, z1, x2, y2, z2}, Reals]]],
      {x0, y0}],
    RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
    RuntimeAttributes -> {Listable},
    Parallelization -> True]];
Protect[AnglesGradientCompiled2D, AnglesGradientCompiled3D];
AnglesGradient[surf_?SurfaceQ, X_] := MapThread[
  Function[
    If[Length[#2] < 2, Table[0,{Length[#1]}],
      With[
        {tr = Transpose[#2]},
        With[
          {rt = RotateLeft /@ tr},
          Total @ AnglesGradienCompiledt3D[
            #1[[1]], #1[[2]], #1[[3]],
            tr[[1]], tr[[2]], tr[[3]],
            rt[[1]], rt[[2]], rt[[3]]]]]]],
  {X, X[[#]] & /@ NeighborhoodList[surf], NeighborhoodAngles[surf]}];
AnglesGradient[surf_?MapQ, X_] := MapThread[
  Function[
    If[Length[#2] < 2, Table[0,{Length[#1]}],
      With[
        {tr = Transpose[#2]},
        With[
          {rt = RotateLeft /@ tr},
          Total @ AnglesGradientCompiled2D[
            #1[[1]], #1[[2]],
            tr[[1]], tr[[2]],
            rt[[1]], rt[[2]]]]]]],
  {X, X[[#]] & /@ NeighborhoodList[surf], NeighborhoodAngles[surf]}];
AnglesGradient[surf_?SurfaceQ, X_, idcs_List] := MapThread[
  If[Length[#2] < 2, Table[0,{Length[#1]}],
    With[
      {tr = Transpose[#2]},
      With[
        {rt = RotateLeft /@ tr},
        Total @ AnglesGradientCompiled3D[
          #1[[1]], #1[[2]], #1[[3]],
          tr[[1]], tr[[2]], tr[[3]],
          rt[[1]], rt[[2]], rt[[3]]]]]],
  {X[[idcs]], X[[#]] & /@ Part[NeighborhoodList[surf], idcs], NeighborhoodAngles[surf][[idcs]]}];
AnglesGradient[surf_?MapQ, X_, idcs_List] := MapThread[
  If[Length[#2] < 2, Table[0,{Length[#1]}],
    With[
      {tr = Transpose[#2]},
      With[
        {rt = RotateLeft /@ tr},
        Total @ AnglesGradientCompiled2D[
          #1[[1]], #1[[2]],
          tr[[1]], tr[[2]],
          rt[[1]], rt[[2]]]]]],
  {X[[idcs]], X[[#]] & /@ Part[NeighborhoodList[surf], idcs], NeighborhoodAngles[surf][[idcs]]}];
AnglesGradient[surf_ /; SurfaceQ[surf] || MapQ[surf]] := ConstantArray[
    0,
    Dimensions[VertexList[surf]]];

(* #TrianglesGradient *****************************************************************************)
TriangleNormalsCompiled = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}},
  With[
    {x = Transpose[xnei]},
    With[
      {u = (RotateLeft /@ x) - x,
       dx = MapThread[Subtract, {x0, x}]},
      With[
        {n = {
          u[[2]]*dx[[3]] - u[[3]]*dx[[2]],
          u[[3]]*dx[[1]] - u[[1]]*dx[[3]], 
          u[[1]]*dx[[2]] - u[[2]]*dx[[1]]}},
        With[
          {d = Sqrt[Total[n^2]]},
          (#/d)& /@ n]]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
Protect[TriangleNormalsCompiled];
TriangleNormals[surf_ /; SurfaceQ[surf], X_] := MapThread[
  TriangleNormalsCompiled,
  {X, X[[#]] & /@ NeighborhoodList[surf]}];
TriangleNormals[surf_ /; SurfaceQ[surf], X_, idcs_List] := MapThread[
  TriangleNormalsCompiled,
  {X[[idcs]], X[[#]] & /@ Part[NeighborhoodList[surf], idcs]}];
TriangleNormals[surf_ /; SurfaceQ[surf]] := Which[
  SurfaceQ[surf] && SurfaceName[surf] =!= surf, TriangleNormals[SurfaceName[surf]],
  True, With[
    {res = Check[TriangleNormals[surf, VertexList[surf]], $Failed]},
    If[!ListQ[res], $Failed, (TriangleNormals[surf] = res)]]];

TriangleAxes2DCompiled = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}},
  With[
    {x = Transpose[xnei]},
    With[
      {u = (RotateLeft /@ x) - x},
      With[
        {n = Sqrt[Total[u^2]]},
        {{u[[1]]/n, u[[2]]/n}, {-u[[2]]/n, u[[1]]/n}}]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
TriangleAxes3DCompiled = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}, {tnorm, _Real, 2}},
  With[
    {x = Transpose[xnei]},
    With[
      {u = (RotateLeft /@ x) - x,
       dx = MapThread[Subtract, {x0, x}]},
      With[
        {n = Sqrt[Total[u^2]]},
        {{u[[1]]/n, u[[2]]/n, u[[3]]/n},
         {tnorm[[2]]*u[[3]]/n - tnorm[[3]]*u[[2]]/n, 
          tnorm[[3]]*u[[1]]/n - tnorm[[1]]*u[[3]]/n,
          tnorm[[1]]*u[[2]]/n - tnorm[[2]]*u[[1]]/n}}]]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
Protect[TriangleAxes2DCompiled, TriangleAxes3DCompiled];
TriangleAxes[surf_ /; SurfaceQ[surf], X_] := With[
  {Xnei = X[[#]] & /@ NeighborhoodList[surf]},
  MapThread[
    TriangleAxes3DCompiled,
    {X, Xnei, MapThread[TriangleNormalsCompiled, {X, Xnei}]}]];
TriangleAxes[surf_ /; SurfaceQ[surf], X_, idcs_List] := With[
  {Xnei = X[[#]] & /@ Part[NeighborhoodList[surf], idcs],
   Xsub = X[[idcs]]},
  MapThread[
    TriangleAxes3DCompiled,
    {Xsub, Xnei, MapThread[TriangleNormalsCompiled, {Xsub, Xnei}]}]];
TriangleAxes[surf_ /; MapQ[surf], X_] := MapThread[
  TriangleAxes2DCompiled,
  {X, X[[#]] & /@ NeighborhoodList[surf]}];
TriangleAxes[surf_ /; MapQ[surf], X_, idcs_List] := MapThread[
  TriangleAxes2DCompiled,
  {X[[idcs]], X[[#]] & /@ Part[NeighborhoodList[surf], idcs]}];
TriangleAxes[surf_ /; SurfaceQ[surf] || MapQ[surf]] := Which[
  MapQ[surf] && MapName[surf] =!= surf, TriangleAxes[MapName[surf]],
  SurfaceQ[surf] && SurfaceName[surf] =!= surf, TriangleAxes[SurfaceName[surf]],
  True, With[
    {res = Check[TriangleAxes[surf, VertexList[surf]], $Failed]},
    If[!ListQ[res], $Failed, (TriangleAxes[surf] = res)]]];

TriangleCoordinatesCompiled = Compile[{{x0, _Real, 1}, {xnei, _Real, 2}, {axes, _Real, 3}},
  With[
    {x = Transpose[xnei]},
    With[
      {dx = MapThread[Subtract, {x0, x}]},
      {Total[axes[[1]] * dx], Total[axes[[2]] * dx]}]],
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
Protect[TriangleCoordinatesCompiled];
TriangleCoordinates[surf_ /; SurfaceQ[surf] || MapQ[surf], X_] := With[
  {Xnei = X[[#]] & /@ NeighborhoodList[surf]},
  MapThread[
    TriangleCoordinatesCompiled,
    {X, Xnei, 
     If[SurfaceQ[surf],
       MapThread[
         TriangleAxes3DCompiled,
         {X, Xnei, MapThread[TriangleNormalsCompiled, {X, Xnei}]}],
       MapThread[TriangleAxes2DCompiled, {X, Xnei}]]}]];
TriangleCoordinates[surf_ /; SurfaceQ[surf] || MapQ[surf], Xall_, idcs_List] := With[
  {Xnei = Xall[[#]] & /@ Part[NeighborhoodList[surf], idcs],
   X = Xall[[idcs]]},
  MapThread[
    TriangleCoordinatesCompiled,
    {X, Xnei, 
     If[SurfaceQ[surf],
       MapThread[
         TriangleAxes3DCompiled,
         {X, Xnei, MapThread[TriangleNormalsCompiled, {X, Xnei}]}],
       MapThread[TriangleAxes2DCompiled, {X, Xnei}]]}]];
TriangleCoordinates[surf_ /; SurfaceQ[surf] || MapQ[surf]] := Which[
  MapQ[surf] && MapName[surf] =!= surf, TriangleCoordinates[MapName[surf]],
  SurfaceQ[surf] && SurfaceName[surf] =!= surf, TriangleCoordinates[SurfaceName[surf]],
  True, With[
    {res = Check[TriangleCoordinates[surf, VertexList[surf]], $Failed]},
    If[!ListQ[res], $Failed, (TriangleCoordinates[surf] = res)]]];

TrianglesGradientCompiled = Compile[
  {{x0, _Real, 1}, {xnei, _Real, 2}, {coords, _Real, 2}, {axes, _Real, 3}, {scale, _Real}},
  With[
    {tmp = -scale * Exp[-scale * coords[[2]]]},
    Total[# * tmp]& /@ axes[[2]]], 
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
Protect[TrianglesGradientCompiled];
TrianglesGradient[surf_ /; SurfaceQ[surf] || MapQ[surf], scale_, X_] := With[
  {Xnei = X[[#]] & /@ NeighborhoodList[surf]},
  With[
    {axes = If[SurfaceQ[surf],
       MapThread[
         If[Length[#2]>0, TriangleAxes3DCompiled[#1,#2,#3], {}],
         {X, Xnei, MapThread[TriangleNormalsCompiled, {X, Xnei}]}],
       MapThread[If[Length[#2]>0, TriangleAxes2DCompiled[#1,#2], {}]&, {X, Xnei}]]},
    With[
      {coords = MapThread[
         If[Length[#2]>0, TriangleCoordinatesCompiled[#1,#2,#3], {}]&,
         {X, Xnei, axes}]},
      MapThread[
        If[Length[#2]>0, TrianglesGradientCompiled[#1,#2,#3,#4,#5], 0*#1]&,
        {X, Xnei, coords, axes, ConstantArray[scale, Length@X]}]]]];
TrianglesGradient[surf_ /; SurfaceQ[surf] || MapQ[surf], scale_, Xall_, idcs_List] := With[
  {Xnei = X[[#]] & /@ Part[NeighborhoodList[surf], idcs],
   X = Xall[[idcs]]},
  With[
    {axes = If[SurfaceQ[surf],
       MapThread[
         If[Length[#2]>0, TriangleAxes3DCompiled[#1,#2,#3], {}],
         {X, Xnei, MapThread[TriangleNormalsCompiled, {X, Xnei}]}],
       MapThread[If[Length[#2]>0, TriangleAxes2DCompiled[#1,#2], {}]&, {X, Xnei}]]},
    With[
      {coords = MapThread[
         If[Length[#2]>0, TriangleCoordinatesCompiled[#1,#2,#3], {}]&,
         {X, Xnei, axes}]},
      MapThread[
        If[Length[#2]>0, TrianglesGradientCompiled[#1,#2,#3,#4,#5], 0*#1]&,
        {X, Xnei, coords, axes, ConstantArray[scale, Length@X]}]]]];
TrianglesGradient[surf_ /; SurfaceQ[surf] || MapQ[surf], scale_] := Which[
  SurfaceQ[surf] && SurfaceName[surf] =!= surf, TrianglesGradient[SurfaceName[surf], scale],
  MapQ[surf] && MapName[surf] =!= surf, TrianglesGradient[MapName[surf], scale],
  True, With[
    {res = Check[TrianglesGradient[surf, scale, VertexList[surf]], $Failed]},
    If[ListQ[res], (TrianglesGradient[surf, scale] = res), $Failed]]];

(* #TrianglesEnergy *******************************************************************************)
TrianglesEnergyCompiled = Compile[
  {{x0, _Real, 1}, {xnei, _Real, 2}, {coords, _Real, 2}, {axes, _Real, 3}, {scale, _Real}},
  Total[Exp[-scale * coords[[2]]]], 
  RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
  Parallelization -> True];
Protect[TrianglesEnergyCompiled];
TrianglesEnergy[surf_ /; SurfaceQ[surf] || MapQ[surf], scale_, X_] := With[
  {Xnei = X[[#]] & /@ NeighborhoodList[surf]},
  With[
    {axes = If[SurfaceQ[surf],
       MapThread[
         If[Length[#2]>0, TriangleAxes3DCompiled[#1,#2,#3], {}],
         {X, Xnei, MapThread[TriangleNormalsCompiled, {X, Xnei}]}],
       MapThread[If[Length[#2]>0, TriangleAxes2DCompiled[#1,#2], {}]&, {X, Xnei}]]},
    With[
      {coords = MapThread[
         If[Length[#2]>0, TriangleCoordinatesCompiled[#1,#2,#3], {}]&,
         {X, Xnei, axes}]},
      Total@MapThread[
        If[Length[#2]>0, TrianglesEnergyCompiled[#1,#2,#3,#4,#5], 0]&,
        {X, Xnei, coords, axes, ConstantArray[scale, Length@X]}]]]];
TrianglesEnergy[surf_ /; SurfaceQ[surf] || MapQ[surf], scale_, Xall_, idcs_List] := With[
  {Xnei = Xall[[#]] & /@ Part[NeighborhoodList[surf], idcs],
   X = Xall[[#]]},
  With[
    {axes = If[SurfaceQ[surf],
       MapThread[
         If[Length[#2]>0, TriangleAxes3DCompiled[#1,#2,#3], {}],
         {X, Xnei, MapThread[TriangleNormalsCompiled, {X, Xnei}]}],
       MapThread[If[Length[#2]>0, TriangleAxes2DCompiled[#1,#2], {}]&, {X, Xnei}]]},
    With[
      {coords = MapThread[
         If[Length[#2]>0, TriangleCoordinatesCompiled[#1,#2,#3], {}]&,
         {X, Xnei, axes}]},
      Total@MapThread[
        If[Length[#2]>0, TrianglesEnergyCompiled[#1,#2,#3,#4,#5], 0]&,
        {X, Xnei, coords, axes, ConstantArray[scale, Length@X]}]]]];
TrianglesEnergy[surf_ /; SurfaceQ[surf] || MapQ[surf], scale_] := Which[
  SurfaceQ[surf] && SurfaceName[surf] =!= surf, TrianglesEnergy[SurfaceName[surf], scale],
  MapQ[surf] && MapName[surf] =!= surf, TrianglesEnergy[MapName[surf], scale],
  True, With[
    {res = Check[TrianglesEnergy[surf, scale, VertexList[surf]], $Failed]},
    If[ListQ[res], (TrianglesEnergy[surf, scale] = res), $Failed]]];

(* #CorticalPotentialField ************************************************************************)
CorticalPotentialField[surf_ /; SurfaceQ[surf] || MapQ[surf], optseq___Rule] := Catch[
  With[
    {opts = With[
       {tmp = {optseq}},
       Join[
         tmp,
         If[FilterRules[tmp, EdgesConstant] == {}, {EdgesConstant -> Automatic}, {}],
         If[FilterRules[tmp, AnglesConstant] == {}, {AnglesConstant -> Automatic}, {}]]]},
    With[
      {parts = Check[
         Select[
           Map[
             Function[
               With[
                 {term = CorticalPotentialTerm[surf, #]},
                 If[term =!= None && !ListQ[term],
                   Message[
                     CorticalPotentialField::badterm, 
                     "term " <> #[[1]] <> " did not yield a 2-element list"];
                   Throw[$Failed]];
                 term]], 
             opts],
           ListQ],
         Throw[$Failed]],
       f = Unique["potential"]},
      SetAttributes[Evaluate[f], Temporary];
      If[Length[parts] == 0,
        Block[{X},
          f[X_ /; ArrayQ[X, 2, NumericQ]] := 0;
          f /: Gradient[f, X_ /; ArrayQ[X, 2, NumericQ]] := ConstantArray[
            0,
            Dimensions@VertexList[s]]],
        Block[{X},
          SetDelayed@@Replace[
            Hold[
              f[X_ /; ArrayQ[X, 2, NumericQ]],
              Evaluate[
                Replace[
                  Hold@@Map[
                    Hold[#[X]]&,
                    parts[[All, 1]]],
                  Hold[f_] :> f,
                  {1}]]],
            Hold[f__] :> Plus[f],
            {1}];
          SetDelayed@@Replace[
            Hold[
              f[X_ /; ArrayQ[X, 2, NumericQ], idcs_List],
              Evaluate[
                Replace[
                  Hold@@Map[
                    Hold[#[X, idcs]]&,
                    parts[[All, 1]]],
                  Hold[f_] :> f,
                  {1}]]],
            Hold[f__] :> Plus[f],
            {1}];
          TagSetDelayed@@Replace[
            Hold[
              f,
              Gradient[f, X_ /; ArrayQ[X, 2, NumericQ]],
              Evaluate[
                Replace[
                  Hold@@Map[
                    Hold[#[X]]&,
                    parts[[All, 2]]],
                  Hold[f_] :> f,
                  {1}]]],
            Hold[f__] :> Flatten[Plus[f]],
            {1}];
          TagSetDelayed@@Replace[
            Hold[
              f,
              Gradient[f, X_ /; ArrayQ[X, 2, NumericQ], idcs_List],
              Evaluate[
                Replace[
                  Hold@@Map[
                    Hold[#[X, idcs]]&,
                    parts[[All, 2]]],
                  Hold[f_] :> f,
                  {1}]]],
            Hold[f__] :> Flatten[Plus[f]],
            {1}]]];
      f]]];

(* #CorticalPotentialTerm *************************************************************************)
CorticalPotentialTerm[s_, opts___] := Message[
  CorticalPotentialTerm::badarg,
  "Unrecognized term: " <> ToString[{opts}]];
CorticalPotentialTerm[s_ /; SurfaceQ[s] || MapQ[s], EdgesConstant -> e_] := With[
  {const = N[e /. Automatic -> (1.0 / Length[EdgeList[s]])]},
  Which[
    !NumericQ[const], (
      Message[CorticalPotentialTerm::badarg, "EdgesConstant must be a number"];
      $Failed),
    const < 0, (
      Message[CorticalPotentialTerm::badarg, "EdgesConstant must be >= 0"];
      $Failed),
    const == 0, None,
    True, {
      (const * If[Length[{##}] == 1, EdgesEnergy[s,#], EdgesEnergy[s,#,{##}[[2]]]])&, 
      (const * If[Length[{##}] == 1, EdgesGradient[s,#], EdgesGradient[s,#,{##}[[2]]]])&}]];
CorticalPotentialTerm[s_ /; SurfaceQ[s] || MapQ[s], AnglesConstant -> e_] := With[
  {const = N[e /. Automatic -> (1.0 / Length[Flatten@NeighborhoodList[s]])]},
  Which[
    !NumericQ[const], (
      Message[CorticalPotentialTerm::badarg, "AnglesConstant must be a number"];
      $Failed),
    const < 0, (
      Message[CorticalPotentialTerm::badarg, "AnglesConstant must be >= 0"];
      $Failed),
    const == 0, None,
    True, {
      (const * If[Length[{##}] == 1, AnglesEnergy[s,#], AnglesEnergy[s,#, {##}[[2]]]])&,
      (const * If[Length[{##}] == 1, AnglesGradient[s,#], AnglesGradient[s,#, {##}[[2]]]])&}]];
CorticalPotentialTerm[s_ /; SurfaceQ[s] || MapQ[s], TrianglesConstant -> e_] := With[
  {const = N[(e /. {c_, sc_} :> c) /. Automatic -> (1.0 / Length[Flatten@NeighborhoodList[s]])],
   scale = N[(e /. {{c_, sc_} :> sc, _ :> 4.0}) /. Automatic -> 4.0]},
  Which[
    !NumericQ[const], (
      Message[CorticalPotentialTerm::badarg, "TrianglesConstant must be a number"];
      $Failed),
    const < 0, (
      Message[CorticalPotentialTerm::badarg, "TrianglesConstant must be >= 0"];
      $Failed),
    const == 0, None,
    True, {
      (const * If[Length[{##}]==1, TrianglesEnergy[s,scale,#], TrianglesEnergy[s,scale,#, {##}[[2]]]])&,
      (const * If[Length[{##}]==1, TrianglesGradient[s,scale,#], TrianglesGradient[s,scale,#, {##}[[2]]]])&}]];
(* Machinery for setting new cortical potential terms *)
SetCorticalPotentialTerm[surf_, rule_, val_] := (
  Unprotect[CorticalPotentialTerm];
  CorticalPotentialTerm /: Set[CorticalPotentialTerm[s_, r_], v_] =.;
  Set[CorticalPotentialTerm[surf, rule], val];
  CorticalPotentialTerm /: Set[CorticalPotentialTerm[s_, r_], v_] := 
    SetCorticalPotentialTerm[s, r, v];
  Protect[CorticalPotentialTerm];
  CorticalPotentialTerm[arg]);
Attributes[SetCorticalPotentialTerm] = Attributes[Set];
CorticalPotentialTerm /: Set[CorticalPotentialTerm[s_, a_],v_] := SetCorticalPotentialTerm[s, a, v];
SetDelayedCorticalPotentialTerm[surf_, rule_, val_] := (
  Unprotect[CorticalPotentialTerm];
  CorticalPotentialTerm /: SetDelayed[CorticalPotentialTerm[s_, r_], v_] =.;
  SetDelayed[CorticalPotentialTerm[surf, rule], val];
  CorticalPotentialTerm /: SetDelayed[CorticalPotentialTerm[s_, r_], v_] := 
    SetDelayedCorticalPotentialTerm[s, r, v];
  Protect[CorticalPotentialTerm];
  Null);
Attributes[SetDelayedCorticalPotentialTerm] = Attributes[SetDelayed];
CorticalPotentialTerm /: SetDelayed[CorticalPotentialTerm[s_, a_],v_] := 
  SetDelayedCorticalPotentialTerm[s, a, v];
UnsetCorticalPotentialTerm[surf_, rule_] := (
  Unprotect[CorticalPotentialTerm];
  CorticalPotentialTerm /: Unset[CorticalPotentialTerm[s_,a_]] =.;
  CorticalPotentialTerm[surf, rule] =.;
  CorticalPotentialTerm /: Unset[CorticalPotentialTerm[s_,a_]] := UnsetCorticalPotentialTerm[s,a];
  Protect[CorticalPotentialTerm];
  Null);
Attributes[UnsetCorticalPotentialTerm] = Attributes[Unset];
CorticalPotentialTerm /: Unset[CorticalPotentialTerm[s_, a_]] := UnsetCorticalPotentialTerm[s,a];


(* #AnglesEnergy **********************************************************************************)
AnglesEnergy[surf_ /; SurfaceQ[surf] || MapQ[surf], X_] := Total[
  Sin[0.5*(Flatten[NeighborhoodAngles[surf]] - Flatten[NeighborhoodAngles[surf, X]])]^2];
(*  (Flatten[NeighborhoodAngles[surf] - NeighborhoodAngles[surf, X]])^2];*)
AnglesEnergy[surf_ /; SurfaceQ[surf] || MapQ[surf], X_, idcs_List] := Total[
  Sin[0.5*(Flatten[NeighborhoodAngles[surf][[idcs]]] - Flatten[NeighborhoodAngles[surf, X, idcs]])^2]];
(*  (Flatten[NeighborhoodAngles[surf][[idcs]] - NeighborhoodAngles[surf, X, idcs]])^2];*)
AnglesEnergy[surf_ /; SurfaceQ[surf] || MapQ[surf]] := 0;

(* #WithField and Machinery for specifying Cortical Surfaces as field -> surface ******************)
WithField[s_?SurfaceQ, f_?FieldQ] := Rule[f, s];
VertexList[Rule[f_?FieldQ, s_?SurfaceQ]] := VertexList[s];
VertexFilter[Rule[f_?FieldQ, s_?SurfaceQ]] := VertexFilter[s];
FaceFilter[Rule[f_?FieldQ, s_?SurfaceQ]] := FaceFilter[s];
Field[Rule[f_?FieldQ, s_?SurfaceQ]] := Field[f];
FaceList[Rule[f_?FieldQ, s_?SurfaceQ]] := FaceList[s];
Polygons[Rule[f_?FieldQ, s_?SurfaceQ]] := Polygons[s];
SurfaceQ[Rule[f_?FieldQ, s_?SurfaceQ]] := True;
SurfaceName[Rule[f_?FieldQ, s_?SurfaceQ]] := SurfaceName[s];

WithField[s_?SurfaceQ, f_List] := Rule[f, s];
VertexList[Rule[f_List, s_?SurfaceQ]] := VertexList[s];
VertexFilter[Rule[f_List, s_?SurfaceQ]] := VertexFilter[s];
FaceFilter[Rule[f_List, s_?SurfaceQ]] := FaceFilter[s];
Field[Rule[f_List, s_?SurfaceQ]] := f;
FaceList[Rule[f_List, s_?SurfaceQ]] := FaceList[s];
Polygons[Rule[f_List, s_?SurfaceQ]] := Polygons[s];
SurfaceQ[Rule[f_List, s_?SurfaceQ]] := True;
SurfaceName[Rule[f_List, s_?SurfaceQ]] := SurfaceName[s];

Radius[Rule[f_?FieldQ, m_?MapQ]] := Radius[m];
ProjectedSurface[Rule[f_?FieldQ, m_?MapQ]] := ProjectedSurface[m];
OrientPoint[Rule[f_?FieldQ, m_?MapQ]] := OrientPoint[m];
InverseProjectionDispatch[Rule[f_?FieldQ, m_?MapQ]] := InverseProjectionDispatch[m];
InverseProjectionTransform[Rule[f_?FieldQ, m_?MapQ]] := InverseProjectionTransform[m];
ProjectionDispatch[Rule[f_?FieldQ, m_?MapQ]] := ProjectionDispatch[m];
ProjectionShear[Rule[f_?FieldQ, m_?MapQ]] := ProjectionShear[m];
SurfaceRotation[Rule[f_?FieldQ, m_?MapQ]] := SurfaceRotation[m];
ProjectionTransform[Rule[f_?FieldQ, m_?MapQ]] := ProjectionTransform[m];
DomainIndices[Rule[f_?FieldQ, m_?MapQ]] := DomainIndices[m];
Domain[Rule[f_?FieldQ, m_?MapQ]] := Domain[m];
Field[Rule[f_?FieldQ, m_?MapQ]] := Part[Field[f], DomainIndices[m]];
MapQ[Rule[f_?FieldQ, m_?MapQ]] := True;
MapName[Rule[f_?FieldQ, m_?MapQ]] := MapName[m];
VertexFilter[Rule[f_?FieldQ, m_?MapQ]] := VertexFilter[m];
FaceFilter[Rule[f_?FieldQ, m_?MapQ]] := FaceFilter[m];
VertexList[Rule[f_?FieldQ, m_?MapQ]] := VertexList[m];
FaceList[Rule[f_?FieldQ, m_?MapQ]] := FaceList[m];
WithField[m_?MapQ, f_?FieldQ] := Rule[f, m];

WithField[m_?MapQ, f_List] := Rule[f, m];
VertexList[Rule[f_List, m_?MapQ]] := VertexList[m];
VertexFilter[Rule[f_List, m_?MapQ]] := VertexFilter[m];
FaceFilter[Rule[f_List, m_?MapQ]] := FaceFilter[m];
Radius[Rule[f_List, m_?MapQ]] := Radius[m];
ProjectedSurface[Rule[f_List, m_?MapQ]] := ProjectedSurface[m];
OrientPoint[Rule[f_List, m_?MapQ]] := OrientPoint[m];
InverseProjectionDispatch[Rule[f_List, m_?MapQ]] := InverseProjectionDispatch[m];
InverseProjectionTransform[Rule[f_List, m_?MapQ]] := InverseProjectionTransform[m];
ProjectionDispatch[Rule[f_List, m_?MapQ]] := ProjectionDispatch[m];
ProjectionShear[Rule[f_List, m_?MapQ]] := ProjectionShear[m];
SurfaceRotation[Rule[f_List, m_?MapQ]] := SurfaceRotation[m];
ProjectionTransform[Rule[f_List, m_?MapQ]] := ProjectionTransform[m];
DomainIndices[Rule[f_List, m_?MapQ]] := DomainIndices[m];
Domain[Rule[f_List, m_?MapQ]] := Domain[m];
Field[Rule[f_List, m_?MapQ]] /; Length[f] == Length[VertexList[m]] := f;
Field[Rule[f_List, m_?MapQ]] /; 
  (Length[f] == Length[VertexList[ProjectedSurface[m]]]) := f[[DomainIndices[m]]];
FaceList[Rule[f_List, m_?MapQ]] := FaceList[m];
MapQ[Rule[f_List, m_?MapQ]] := True;
MapName[Rule[f_List, m_?MapQ]] := MapName[m];

Unprotect[Rule];

Rule /: Center[Rule[f_?FieldQ, m_?MapQ]] := Center[m];
Rule /: Center[Rule[a_List, m_?MapQ]] := Center[m];

Rule /: Normal[Rule[f_?FieldQ, s_?SurfaceQ]] := MapThread[Rule, {VertexList[s], Field[f]}];
Rule /: Normal[Rule[f_?FieldQ, m_?MapQ]] := MapThread[
  Rule,
  {VertexList[m], Part[Field[f], DomainIndices[m]]}];
Rule /: Normal[Rule[a_List, s_?SurfaceQ]] := MapThread[Rule, {VertexList[s], a}];
Rule /: Normal[Rule[a_List, m_?MapQ]] /; Length[a] == Length[DomainIndices[m]] := MapThread[
  Rule,
  {VertexList[m], a}];
Rule /: Normal[Rule[a_List, m_?MapQ]] /; 
  (Length[a] == Length[VertexList[ProjectedSurface[m]]]) := MapThread[
    Rule,
    {VertexList[m], a[[DomainIndices[m]]]}];

Rule[f_?FieldQ, s_?SurfaceQ] /; Length[Field[f]] != Length[Field[s]] := (
  Message[WithField::incompat];
  $Failed);
Rule[f_?FieldQ, m_?MapQ] /; Length[Field[f]] != Length[Field[ProjectedSurface[m]]] := (
  Message[WithField::incompat];
  $Failed);
Rule[a_List, b_?SurfaceQ] /; Length[a] != Length[VertexList[b]] := (
  Message[WithField::incompat];
  $Failed);
Rule[a_List, b_?MapQ] /; (
  Length[a] != Length[VertexList[b]] && Length[a] != Length[VertexList[ProjectedSurface[b]]]
  ) := (
    Message[WithField::incompat];
    $Failed);

Protect[Rule];

ToField[dat_List] := (
  Unprotect[Field];
  Field[dat] = dat;
  Protect[Field];
  dat);

Faces[s_] := FaceList[s];
FaceCount[s_] := Length[FaceList[s]];
VertexCount[s_] := Length[VertexList[s]];
EdgeCount[s_] := Length[EdgeList[s]];

Protect[SphericalAzimuth, Cartesian, CartesianToSpherical,
        ConvertCoordinates, Surface, SurfaceFromVTK, AnglesConstant,
        AnglesEnergy, AnglesGradient, CorticalPotentialField,
        CorticalPotentialTerm, Domain, DomainIndices, Duplicate,
        EdgesConstant, EdgesEnergy, EdgesGradient, EdgeList,
        FaceAngles, FaceList, Faces, FacesIndex, Field, FieldQ,
        FaceFilter, VertexFilter, VertexList,
        InverseProjectionDispatch, InverseProjectionTransform,
        Latitude, Longitude, MapHull, MapName, NeighborhoodAngles,
        NeighborhoodEdgeLengths, NeighborhoodList, MapMeshPlot,
        MapPlot, MapQ, MergeSurfaces, OrientMatrix, OrientPoint,
        SphericalPolarAngle, Polygons, ProjectionDispatch,
        ProjectionRotation, ProjectionShear, ProjectionTransform,
        Radius, ReadVTK, SphericalCoordinateStyle,
        SphericalToCartesian, ProjectedSurface, SurfaeAnglesGradient,
        SurfacePlot, SurfaceProjection, SurfaceQ, SurfaceRotation,
        SurfaceResample, ToField, TrianglesConstant, TrianglesEnergy,
        TrianglesGradient, VertexIndexDispatch, WithField,
        WithFaceFilter, WithVertexFilter, WithVertexList];

ColorCortex[instructions___] := Block[{tmp},
   With[
    {ord = Join @@ Append[
      Map[
        Function[
          Replace[
            #,
            Hold[param_, call_] :> Hold[             
              param =!= None && param =!= $Failed && (tmp = call) =!= None && tmp =!= $Failed,
              tmp]]],
        Replace[
          Hold[instructions],
          {Rule[idx_Integer, instr_] :> With[
             {cc = CorticalColor[instr]},
             If[ListQ[cc],
               With[
                 {min = cc[[1, 1]], max = cc[[1, 2]], clrs = cc[[2]]},
                 Hold[
                   Part[Slot[1], idx],
                   If[ListQ[Slot[1]] && Length[Slot[1]] >= idx && NumericQ[Part[Slot[1], idx]],
                     Blend[clrs, Rescale[Part[Slot[1], idx], {min, max}]],
                    $Failed]]],
               Hold[
                 Part[Slot[1], idx],
                 If[ListQ[Slot[1]] && Length[Slot[1]] >= idx,
                   cc[Part[Slot[1], idx]],
                   $Failed]]]],
           instr_ :> With[
             {cc = CorticalColor[instr]},
             If[ListQ[cc],
               With[
                 {min = cc[[1, 1]], max = cc[[1, 2]], clrs = cc[[2]]},
                 Hold[
                   Slot[1],
                   If[NumericQ[Slot[1]],
                     Blend[clrs, Rescale[Slot[1], {min, max}]],
                     $Failed]]],
               Hold[
                 Slot[1],
                 cc[Slot[1]]]]]},
          {1}]],
      Hold[True, Gray]]},
    Function @@ Replace[
      ord,
      Hold[body__] :> Hold[Module[{tmp}, Which[body]]]]]];

Curvature = Curvature;
CorticalColor[Curvature] = Function[If[# < -0.02, GrayLevel[0.55], GrayLevel[0.2]]];

(* We have a sneaky way of setting cortical colors that allows it to be still protected. *)
SetCorticalColor[arg_, val_] := (
  Unprotect[CorticalColor];
  CorticalColor /: Set[CorticalColor[a_], v_] =.;
  Set[CorticalColor[arg], val];
  CorticalColor /: Set[CorticalColor[a_], v_] := SetCorticalColor[a, v];
  Protect[CorticalColor];
  CorticalColor[arg]);
Attributes[SetCorticalColor] = Attributes[Set];
CorticalColor /: Set[CorticalColor[a_], v_] := SetCorticalColor[a, v];
SetDelayedCorticalColor[arg_, val_] := (
  Unprotect[CorticalColor];
  CorticalColor /: SetDelayed[CorticalColor[a_], v_] =.;
  SetDelayed[CorticalColor[arg], val];
  CorticalColor /: SetDelayed[CorticalColor[a_], v_] := SetDelayedCorticalColor[a, v];
  Protect[CorticalColor];
  Null);
Attributes[SetDelayedCorticalColor] = Attributes[SetDelayed];
CorticalColor /: SetDelayed[CorticalColor[a_], v_] := SetDelayedCorticalColor[a, v];
UnsetCorticalColor[arg_] := (
  Unprotect[CorticalColor];
  CorticalColor /: Unset[CorticalColor[a_]] =.;
  Unset[CorticalColor[arg]];
  CorticalColor /: Unset[CorticalColor[a_]] := UnsetCorticalColor[a];
  Protect[CorticalColor];
  Null);
Attributes[UnsetCorticalColor] = Attributes[Unset];
CorticalColor /: Unset[CorticalColor[a_]] := UnsetCorticalColor[a];

Protect[CorticalColor, ColorCortex, Curvature];

End[];
EndPackage[];
