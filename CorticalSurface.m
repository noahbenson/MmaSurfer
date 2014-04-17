(* CorticalSurface.m
 *
 * Utility functions for dealing with (mostly FreeSurfer) cortical surface data in Mathematica.
 *
 * Copyright (C) 2013-2014 by Noah C. Benson.
 * This file is part of the MmaSurfer library, which is provided under the terms of the Eclipse
 * Public License, version 1.0. See the accompanying LICENSE file for more information.
 *)

(**************************************************************************************************)
BeginPackage["CorticalSurface`"];

Unprotect[Azimuth, Cartesian, CartesianToSpherical,
          ConvertCoordinates, CorticalSurface, CorticalSurfaceFromVTK,
          Domain, DomainIndices, Duplicate, Faces, Field, FieldQ,
          Filter, InverseProjectionDispatch,
          InverseProjectionTransform, Latitude, Longitude, MapPlot,
          MapQ, MergeSurfaces, OrientMatrix, OrientPoint, PolarAngle,
          Polygons, ProjectionDispatch, ProjectionRotation,
          ProjectionShear, ProjectionTransform, Radius, ReadVTK,
          SphericalCoordinateStyle, SphericalToCartesian, Surface,
          SurfaceCases, SurfaceMap, SurfacePlot, SurfaceProjection,
          SurfaceQ, SurfaceReplace, SurfaceResample, SurfaceRotation,
          SurfaceSelect, ToField, Vertices, WithField, WithFilter];

ClearAll[Azimuth, Cartesian, CartesianToSpherical, ConvertCoordinates,
         CorticalSurface, CorticalSurfaceFromVTK, Domain,
         DomainIndices, Duplicate, Faces, Field, FieldQ, Filter,
         InverseProjectionDispatch, InverseProjectionTransform,
         MapPlot, MapQ, MergeSurfaces, OrientMatrix, OrientPoint,
         PolarAngle, Polygons, ProjectionDispatch, ProjectionRotation,
         ProjectionShear, ProjectionTransform, ReadVTK,
         SphericalCoordinateStyle, SphericalToCartesian, Surface,
         SurfaceCases, SurfaceMap, SurfacePlot, SurfaceProjection,
         SurfaceQ, SurfaceReplace, SurfaceResample, SurfaceRotation,
         SurfaceSelect, ToField, WithField, WithFilter];

(*
Faces::usage = "Faces[surf] yields a list of all faces in the surface surf, if any, as lists of
 indices into Vertices[surf].";
Longitude::usage = "Longitude is a keyword that represents the longitude (identical to azimuth
 angle) in spherical coordinates.";
Latitude::usage = "Latitude is a keyword that represents the latitude in spherical coordinates.";
Radius::usage = "Radius is a keyword that represents the radius of a surface-to-map projection (see
 SurfaceProjection).";
Vertices::usage = "Vertices[surf] yields a list of all vertices associated with cortical surface
 surf in cartesian coordinates.  Vertices[surf, rep] yields a list of all vertices associated with
 cortical surface surf in the coordinate representation given by rep; rep can be in any format
 accepted by ConvertCoordinates.";
*)

CartesianToSpherical::usage = "CartesianToSpherical[data] yields an equivalent dataset to data
 except that the parts of data corresponding to (x,y,z) coordinates are replaced with
 (\[Theta],\[Phi]) coordinates. The data parameter must be a list of either lists or rules OR a
 single list or a single rule. In data is a list of lists, then the first three elements of each
 sublist is taken to be the (x,y,z) part; if data is a list of rules or a rule, then the first
 element of the rule is taken to be the (x,y,z) part. CartesianToSpherical yields data in the same
 format as it is given, but with the appropriate replacements."
CartesianToSpherical::badfmt = "Bad data format given to CartesianToSpherical";
CartesianToSphercial::badopt = "Bad option `1` given to CartesianToSpherical: `2`";

SphericalToCartesian::usage = "SphericalToCartesian is identical to CartesianToSpherical, but
 performs the inverse transform."
SphericalToCartesian::badfmt = "Bad data format given to SphericalToCartesian";
SphercialToCartesian::badopt = "Bad option `1` given to SphericalToCartesian: `2`";

SphericalCoordinateStyle::usage = "SphericalCoordinateStyle is an option for functions such as
 CorticalSurface that accept spherical coordinates. SphericalCoordinateStyle describes coordinates
 of the form {a, b} such that a SphericalCoordinateStyle value of {x, y} indicates that a represents
 measurement x and b represents measurement y. The possible measurements are Azimuth/Longitude
 (which are identical), Latitude, and PolarAngle. At least one of x and y must be Azimuth or
 Longitude and one of them must be either Latitude or PolarAngle. The ConvertCoordinates function
 can be used to convert between representations.";
PolarAngle::usage = "PolarAngle is a keyword that either represents polar angle data of the visual
 field or the angle from the pole (z-axis) in spherical coordinate data.";
Azimuth::usage = "Azimuth is a keyword that represents the azimuth angle (identical to longitude)
 in spherical coordinates.";

ConvertCoordinates::usage = "ConvertCoordinates[data, fromStyle -> toStyle] yields a transformation
 of data such that its coordinates are converted from fromStyle to toStyle; fromStyle and toStyle
 may be Cartesian or any argument appropriate for SphericalCoordinateStyle. The data parameter may
 be a list of points or a list of rules whose heads are points.";
ConvertCoordinates::badfmt = "Bad input to ConvertCoordinates: `1`";

CorticalSurface::usage = "CorticalSurface[data] yields a cortical surface representation of the
 cortical surface data given in data. The cortical surface data must be a list of rules, which may
 be in spherical or cartesian coordinates. If the data is in spherical coordinates, then the
 SphericalCoordinateStyle option specifies how the coordinates should be interpreted; if the
 coordinates for each point are {a, b}, then the SphericalCoordinateStyle option specifies {a-style,
 b-style} where each of a-style and b-style must come from these lists and must not come from the
 same list: {Azimuth, Longitude} and {Latitude, PolarAngle}. See ?SphericalCoordinateStyle for more
 information.
The cortical surface representation yielded by this function is actually a symbol with many tagged
 values. These each have their own additional documentation and are listed here:
  Vertices[surf],
  Faces[surf]
  Field[surf],
  SurfaceQ[surf]
  Normal[surf]
  Filter[surf]";
CorticalSurface::badarg = "Bad argument `1` to CorticalSurface: `2`";

CorticalSurfaceFromVTK::usage = "CorticalSurfaceFromVTK[filename] yields a cortical surface
 representation of the cortical surface data given in the VTK file specified; see CorticalSurface
 and ReadVTK for more information.";

Filter::usage = "Filter is an option to CorticalSurface and SurfaceProjection; in both cases, it
 must be a function that takes as arguments a vertex position (in cartesian or map coordinates) and
 a field value (which may be None). Only those vertices for which the filter function yields True
 are included in the constructed object. In the case of maps, vertices that are filtered are still
 included in the domain, but not in the result. Filter[surface] and Filter[map] yield the surface
 or map's filter.";
Polygons::usage = "Polygons[surf] yields a list of all polygons in the surface surf, if any, by
 cross-referencing Vertices[surf] and Faces[surf]; these polygons are listed in the same order, with
 the same vertex order, as in Faces[surf].";
Field::usage = "Field[surf] yields a list of all field values associated with cortical surface
 surf.";
FieldQ::usage = "FieldQ[obj] yields true if and only if obj has a field such that Field[obj] yields
 a list.";
Surface::usage = "Surface[map] yields the surface object from which the map was projected.";
SurfaceQ::usage = "SurfaceQ[s] yields true if and only if s is a surface object, as created by
 CorticalSurface[].";

Domain::usage = "Domain[map] yields the Normal of the section of the surface from which map was
 derived.";
DomainIndices::usage = "DomainIndices[map] yields a list of the indices into Normal[Surface[map]]
 for the points that compose the domain of map.";
ProjectionDispatch::usage = "ProjectionDispatch[map] yields a dispatch list (ie, for use with
 Replace and /.) that, when given a cartesian point on a surface projects it to the map.";
ProjectionShear::usage = "ProjectionShear[map] yields the shearing matrix associated with map.";
ProjectionRotation::usage = "ProjectionRotation[map] yields the rotation matrix associated with
 map.";
ProjectionTransform::usage = "ProjectionTransform[map] yields a function that converts cartesian
 coordinates into map coordinates for the projection that creates map.";
InverseProjectionDispatch::usage = "InverseDispatch[map] yields the dispatch function that, when
 used to replace something, inverses the transformation from the surface to the map.";
InverseProjectionTransform::usage = "InverseProjectionTransform[map] yields a function that inverts
 the transform from map's surface to map."
MapQ::usage = "MapQ[x] yields true if and only if x is a map created as a projection of a surface.";

MergeSurfaces::usage = "MergeSurfaces[s1, s2, f] yields a surface object such that any point p
 appearing in either s1 or s2 is represented in the new surface and has field value equal to
 f[f1, f2] where f1 and f2 are the field values at p in surface s1 and s2 respectively (or None if
 no field is defined at that point).";
MergeSurfaces::dup = "Surface given to MergeSurfaces has duplicate vertices";

SurfaceProjection::usage = "SurfaceProjection[surf] yields a surface projection of the region within
 Pi/3 radians of the occipital pole, with V1 oriented along the positive x-axis. To alter the
 location of the map projection constructed, the following options may be used:
  Duplicate may be a surface projection whose projection should be duplicated (other arguments
    overwrite this argument).
  Center must be a 3D cartesian point representing the center of the map to be constructed.
  Radius must be a number greater than 0 that indicates the distance from the center that the map
    should encompass.
  OrientPoint must be a point or a point -> an angle (which be default is 0) to use in orienting the
    map; the point is oriented so that it lies at the given angle in a polar-coordinate version of
    the map (i.e., x -> Pi/4 indicates that point x should, in the resulting map, lie on the
    y-axis).
  ProjectionShear must be a shear matrix that is applied to the final map after orientation.";
SurfaceProjection::badarg = "Bad argument `1` to SurfaceProjection: `2`";
SurfaceProjection::badfmt = "Transform argument given in unrecognized format";
SurfaceRotation::usage = "SurfaceRotation[map] yields the rotation matrix that produces the map.";

OrientPoint::usage = "OrientPoint is a keyword that represents the ultimate orientation of a
 surface-to-map projection (see SurfaceProjection). OrientPoint -> (x -> 0) is identical to
 OrientPoint -> x; otherwise, OrientPoint -> (x -> y) ensures that in the final map projection, the
 point x will be projected to point {qx,qy} such that the ArcTan[qx, qy] == y.";
OrientMatrix::usage = "OrientMatrix[map] yields the orientation matrix applied to the transformation
 that produced the surface projection, map.";
Duplicate::usage = "Duplicate is an option to SurfaceProjection which specifies that the projection
 should, unless otherwise specified in the options list, use the transformation arguments used by
 the given map.";

MapPlot::usage = "MapPlot[map] yields a 2D map plot of the given map; all options for Graphics[]
 may be given.";

WithFilter::usage = "WithFilter[map, filt] and WithFilter[surf, filt] yield a map or surface object
 that is identical to map or surf except that it has an additional filter applied to it.";
WithField::usage = "WithField[surf, field] yields a surface object that is identical to surf except
 that its field has been replaced with the given field. Rule (field -> surf) may be used as a
 shorthand for WithField[surf, field].";
WithField::incompat = "WithField given surface and field with a different numbers of vertices.";

ToField::usage = "ToField[list] yields a version of list that can be passed as a field.";

ReadVTK::usage = "ReadVTK[filename] reads the given VTK file, including fields, and yields a list
 of rules, {x,y,z} -> field. The filename may be a URL.";
ReadVTK::nofile = "No such file: `1`";
ReadVTK::flerr = "Unspecified error while reading file: `1`";

SurfacePlot::usage = "SurfacePlot[surf] yields a plot of the surface object surf. All options
 available to Graphics3D may be given,";
SurfaceSelect::usage = "SurfaceSelect[surf, fn] yields a surface object which consists of a subset
 of the surface surf containing nodes U for which fn[u -> f[u]] yields true for all u in U and where
 f[u] is the field value of vertex u.";
SurfaceCases::usage =  "SurfaceCases[surf, form] yields a surface object which consists of a subset
 of the surface surf containing nodes U for which u -> f[u] matches form for all u in U and where
 f[u] is the field value of vertex u.";
SurfaceMap::usage = "SurfaceMap[fn, surf] yields a new surface equal to the surface surf but with
 the new surface's vertex/field values determined by fn[u -> f[u]] where u is in Vertices[surf] and
 f[u] is the field of vertex u. SurfaceMap[fn, {surf1, surf2...}] is like SurfaceMap[fn, surfn] but
 calls fn with the vertex to field rules of all the given surfaces.";
SurfaceReplace::usage = "SurfaceReplace[surf, rules] yields a new surface identical to the surface
 surf after having each vertex -> field element of the surface replaced by the list of rules.";

SurfaceResample::usage = "SurfaceResample[surf1, surf2] yields a surface equivalent to surf2 but
 such that the field of the surface has been resampled from surf1. All options except the Field
 option of CorticalSurface are accepted and passed verbatim; a Method option may also specify
 Nearest (default) for nearest-neighbor interpolation, Interpolation, or List interpolation, for
 their respective functions. In the latter two cases, A list may be given instead of the argument
 such that the first argument is Interpolation or LitInterpolation and the remaining elements of
 the list are options to pass to these functions; e.g. Method -> {Interpolation, 
 InterpolationOrder -> 4}. Note that surf1 -> surf2 is equivalent to SurfaceResample[surf1, surf2]
 but self-hashes.";

(**************************************************************************************************)
Begin["`Private`"];

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
        (* At this point, we just need to put dat in the form of a list of {x,y,z} -> value rules. *)
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
     {{Latitude, Longitude|Azimuth} -> Function[
        {ArcSin[#[[3]]/Norm[#]],
         If[#[[1]] == 0 && #[[2]] == 0, 0, ArcTan[#[[1]], #[[2]]]]}],
      {Longitude|Azimuth, Latitude} -> Function[
        {If[#[[1]] == 0 && #[[2]] == 0, 0, ArcTan[#[[1]], #[[2]]]],
         ArcSin[#[[3]]/Norm[#]]}],
      {PolarAngle, Longitude|Azimuth} -> Function[
        {ArcCos[#[[3]]/Norm[#]],
         If[#[[1]] == 0 && #[[2]] == 0, 0, ArcTan[#[[1]], #[[2]]]]}],
      {Longitude|Azimuth, PolarAngle} -> Function[
        {If[#[[1]] == 0 && #[[2]] == 0, 0, ArcTan[#[[1]], #[[2]]]],
         ArcCos[#[[3]]/Norm[#]]}],
      _ :> (
        Message[CartesianToSpherical::badopt, SphericalCoordinateStyle, "Style not recognized"];
        $Failed)}]},
  If[scFn === $Failed, $Failed, Map[scFn, dat]]];
CartesianToSpherical[dat : {Rule[{_, _, _}, _] ..}, OptionsPattern[]] := With[
  {scFn = Replace[
     OptionValue[SphericalCoordinateStyle],
     {{Latitude, Longitude|Azimuth} -> Function[
        Rule[
          {ArcSin[#[[1,3]]/Norm[#[[1]]]],
           If[#[[1,1]] == 0 && #[[1,2]] == 0, 0, ArcTan[#[[1,1]], #[[1,2]]]]},
          #[[2]]]],
      {Longitude|Azimuth, Latitude} -> Function[
        Rule[
          {If[#[[1,1]] == 0 && #[[1,2]] == 0, 0, ArcTan[#[[1,1]], #[[1,2]]]],
           ArcSin[#[[1,3]]/Norm[#[[1]]]]},
          #[[2]]]],
      {PolarAngle, Longitude|Azimuth} -> Function[
        Rule[
          {ArcCos[#[[1,3]]/Norm[#[[1]]]],
           If[#[[1,1]] == 0 && #[[1,2]] == 0, 0, ArcTan[#[[1,1]], #[[1,2]]]]},
          #[[2]]]],
      {Longitude|Azimuth, PolarAngle} -> Function[
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
     {{Latitude, Longitude|Azimuth} -> Function[
        Append[Cos[#[[1]]] * {Cos[#[[2]]], Sin[#[[2]]]}, Sin[#[[1]]]]],
      {Longitude|Azimuth, Latitude} -> Function[
        Append[Cos[#[[2]]] * {Cos[#[[1]]], Sin[#[[1]]]}, Sin[#[[2]]]]],
      {PolarAngle, Longitude|Azimuth} -> Function[
        Append[Sin[#[[1]]] * {Cos[#[[2]]], Sin[#[[2]]]}, Cos[#[[1]]]]],
      {Longitude|Azimuth, PolarAngle} -> Function[
        Append[Sin[#[[2]]] * {Cos[#[[1]]], Sin[#[[1]]]}, Cos[#[[2]]]]],
      _ :> (
        Message[SphericalToCartesian::badopt, SphericalCoordinateStyle, "Style not recognized"];
        $Failed)}]},
  If[scFn === $Failed, $Failed, Map[scFn, dat]]];
SphericalToCartesian[dat : {Rule[{_, _}, _] ..}, OptionsPattern[]] := With[
  {scFn = Replace[
     OptionValue[SphericalCoordinateStyle],
     {{Latitude, Longitude|Azimuth} -> Function[
        Rule[
          Append[Cos[#[[1,1]]] * {Cos[#[[1,2]]], Sin[#[[1,2]]]}, Sin[#[[1,1]]]],
          #[[2]]]],
      {Longitude|Azimuth, Latitude} -> Function[
        Rule[
          Append[Cos[#[[1,2]]] * {Cos[#[[1,1]]], Sin[#[[1,1]]]}, Sin[#[[1,2]]]],
          #[[2]]]],
      {PolarAngle, Longitude|Azimuth} -> Function[
        Rule[
          Append[Sin[#[[1,1]]] * {Cos[#[[1,2]]], Sin[#[[1,2]]]}, Cos[#[[1,1]]]],
          #[[2]]]],
      {Longitude|Azimuth, PolarAngle} -> Function[
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

(* #CorticalSurface *******************************************************************************)
Options[CorticalSurface] = {
  SphericalCoordinateStyle -> {Longitude, Latitude}, 
  Filter -> None,
  Field -> None,
  Faces -> None,
  Polygons -> None};
CorticalSurface[data: {{_,_,_}..}, OptionsPattern[]] := Catch[
  With[
    {sym = Unique["surf"],
     field = Replace[
       OptionValue[Field],
       {l_List /; Length[l] == Length[data] :> l,
        None -> Table[None, {Length[data]}],
        _ :> (
          Message[CorticalSurface::badarg, Field, "must have same length as vertices"];
          Throw[$Failed])}],
     faces = Replace[
       OptionValue[Faces],
       {(l  : {{_Integer, _Integer, _Integer..}..} 
           /; And@@Map[(1 <= # <= Length[data])&, Flatten[l]]
         ) :> l,
        None -> None,
        _ :> (
          Message[CorticalSurface::badarg, Faces, "must consist of index lists"];
          Throw[$Failed])}],
     polygons = Replace[
       OptionValue[Polygons],
       {l_List /; TrueQ[Union[data] == Union[Join[data, Flatten[l, 2]]]] :> l,
        None -> None,
        _ :> (
          Message[CorticalSurface::badarg, Faces, "must consist of index lists"];
          Throw[$Failed])}],
     filtFn = Replace[
       OptionValue[Filter], 
       (False|None) :> (True&)]},
    sym /: SurfaceQ[sym] = True;
    sym /: Filter[sym] = OptionValue[Filter];
    With[
      {idcs = Select[Range[Length[data]], TrueQ[filtFn[data[[#]], field[[#]]]]&]},
      sym /: Vertices[sym] := TagSet[
        sym, 
        Vertices[sym],
        data[[idcs]]];
      sym /: Vertices[sym, toStyle_] := TagSet[
        sym, 
        Vertices[sym, toStyle],
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
      sym /: Faces[sym] := TagSet[
        sym, 
        Faces[sym],
        If[faces === None && polygons === None, 
          None,
          With[
            {keep = Normal[SparseArray[Thread[idcs -> 1], Length[data], 0]]},
            Select[
              If[faces =!= None,
                faces,
                polygons /. Dispatch[MapThread[Rule, {data, Range[Length[data]]}]]],
              Apply[And, Thread[Evaluate[Part[keep, #] == 1]]]&]]]];
      sym /: Polygons[sym] := TagSet[
        sym, 
        Polygons[sym],
        If[faces === None && polygons === None, 
          None,
          With[
            {keep = Dispatch[Append[Table[data[[i]] -> True, {i, idcs}], _ -> False]]},
            Select[
              If[polygons =!= None,
                polygons,
                faces /. Dispatch[MapThread[Rule, {Range[Length[data]], data}]]],
            Apply[And, Replace[#, keep, {1}]]&]]]];
      sym /: Polygons[sym, toStyle_] := TagSet[
        sym, 
        Polygons[sym, toStyle],
        Map[ConvertCoordinates[#, Cartesian -> toStyle]&, Polygons[sym]]];
      sym]]];
CorticalSurface[data: {Rule[{_,_,_},_]..}, opts:OptionsPattern[]] := CorticalSurface[
  data[[All,1]],
  opts,
  Field -> data[[All,2]]];
CorticalSurface[data: {{_,_}..}, opts:OptionsPattern[]] := With[
  {sc = OptionValue[SphericalCoordinateStyle]},
  CorticalSurface[
    ConvertCoordinates[data, sc -> Cartesian],
    opts]];
CorticalSurface[data: {Rule[{_,_},_]..}, opts:OptionsPattern[]] := CorticalSurface[
  data[[All, 1]],
  opts,
  Field -> data[[All, 2]]];

(* #SurfaceQ **************************************************************************************)
SurfaceQ[x_] = False;

(* #CorticalSurfaceFromVTK ************************************************************************)
Options[CorticalSurfaceFromVTK] = Options[CorticalSurface];
CorticalSurfaceFromVTK[filename_String, opts:OptionsPattern[]] := CorticalSurface[
  ReadVTK[filename],
  opts];

(* #MergeSurfaces *********************************************************************************)
Options[MergeSurfaces] = {
  Function -> Function[#2],
  Filter -> None};
MergeSurfaces[surfs:{_?SurfaceQ, _?SurfaceQ..}, OptionsPattern[]] := With[
  {n = Length[surfs],
   ideal = Range[Length[surfs]],
   mergeFn = OptionValue[Function],
   filt = OptionValue[Filter]},
  CorticalSurface[
    Map[
      Function[
        If[Length[Union[#[[All,3]]]] != Length[#],
          (Message[MergeSurfaces::dup]; #[[1,1]] -> Indeterminate),
          #[[1,1]] -> mergeFn[
            #[[1,1]],
            Replace[Range[n], Append[Map[(#[[3]] -> #[[2]])&, #], (_ -> None)], {1}]]]],
      GatherBy[
        Apply[
          Join,
          MapIndexed[
            Function[{surf, idx}, 
              MapThread[{#1, #2, idx[[1]]}&, {Vertices[surf], Field[surf]}]],
            surfs]],
        First]],
    Faces -> Faces[surfs[[1]]],
    Filter -> filt]];

(* #SurfaceProjection *****************************************************************************)
Options[SurfaceProjection] = {
  Center -> Automatic,
  ProjectionShear -> Automatic,
  Radius -> Automatic,
  OrientPoint -> Automatic,
  Duplicate -> None,
  Filter -> None};
SurfaceProjection[surf_, OptionsPattern[]] := Catch[
  With[
    {dflt = Replace[
       OptionValue[Duplicate],
       {(None|False) -> {
          OrientPoint -> ({0,1,0} -> 0),
          Center -> {0,0,1},
          Radius -> (Pi / 3.0),
          ProjectionShear -> {{1,0},{0,1}},
          Filter -> None},
        map_?MapQ :> {
          OrientPoint -> OrientPoint[map],
          Center -> Center[map],
          Radius -> Radius[map],
          ProjectionShear -> ProjectionShear[map],
          Filter -> Filter[map]},
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
            Message[SurfaceProjection::badarg, ProjectionShear, "must be of the form {{1, Sx}, {Sy, 1}}"];
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
       filt = OptionValue[Filter],
       sym = Unique["map"],
       scFn = Function[CartesianToSpherical[#, SphericalCoordinateStyle -> {Longitude, Latitude}]],
       iscFn = Function[SphericalToCartesian[#, SphericalCoordinateStyle -> {Longitude, Latitude}]]},
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
          With[
            {fn = Function[Dot[scFn[Dot[#, rotT]], opT]],
             ifn = Function[Dot[iscFn[Dot[#, iopT]], irotT]]},
            sym /: MapQ[sym] = True;
            sym /: Filter[sym] = filt;
            sym /: Surface[sym] = surf;
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
               s_?SurfaceQ :> MapThread[Rule, {fn[Vertices[s]], Field[s]}],
               _ :> (Message[SurfaceProjection::badfmt]; Indeterminate)}];
            sym /: InverseProjectionDispatch[sym] = Dispatch[
              {pat:{Rule[{_,_},_]..} :> MapThread[Rule, {ifn[pat[[All,1]]], pat[[All,2]]}],
               pat:{{_,_}..} :> ifn[pat],
               Rule[a:{_,_},b_] :> Rule[First[ifn[{a}]], b],
               a:{_?AtomQ,_?AtomQ} :> First[ifn[{a}]],
               s_?SurfaceQ :> MapThread[Rule, {ifn[Vertices[s]], Field[s]}],
               _ :> (Message[SurfaceProjection::badfmt]; Indeterminate)}];
            sym /: ProjectionTransform[sym] = Function[Replace[#, ProjectionDispatch[sym]]];
            sym /: InverseProjectionTransform[sym] = Function[
              Replace[#, InverseProjectionDispatch[sym]]];
            sym /: DomainIndices[sym] := TagSet[
              sym,
              DomainIndices[sym],
              With[
                {cosrad = Cos[rad],
                 normmu = Normalize[mu]},
                Flatten[
                  Position[
                    Normal[surf],
                    Rule[xyz_,_] /; Dot[Normalize[xyz], normmu] >= cosrad]]]];
            sym /: Domain[sym] := TagSet[
              sym,
              Domain[sym],
              Part[Normal[surf], DomainIndices[sym]]];
            sym /: Normal[sym] := TagSet[
              sym,
              Normal[sym],
              Replace[
                If[filt === None, Domain[sym], Select[Domain[sym], TrueQ[f[#[[1]], #[[2]]]]&]],
                ProjectionDispatch[sym]]];
            sym /: Vertices[sym] := TagSet[sym, Vertices[sym], Normal[sym][[All,1]]];
            sym /: Field[sym] := TagSet[sym, Field[sym], Normal[sym][[All,2]]];
            sym /: Faces[sym] := TagSet[
              sym,
              Faces[sym],
              With[
                {dom = Dispatch[
                   Append[
                     Map[(# -> True)&, DomainIndices[sym]],
                     _Integer -> False]],
                 tr = Dispatch[
                   MapThread[
                     Rule, 
                     {DomainIndices[sym], Range[Length@DomainIndices[sym]]}]]},
                Replace[
                  Select[Faces[surf], Apply[And, # /. dom]&],
                  tr,
                  {2}]]];
            sym /: Polygons[sym] := TagSet[
              sym,
              Polygons[sym],
              Map[Vertices[sym][[#]]&, Faces[sym]]];
            sym]]]]]];

(* #StructuredRescale * ***************************************************************************)
(* This is used be MapPlot and SurfacePlot when normalizing oddly structured  data sets in the color
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
   ColorFunctionScaling -> True}];
MapPlot[map_?MapQ, opts:OptionsPattern[]] := Graphics[
  With[
    {V = Normal[Vertices[map]],
     F = Normal[Faces[map]],
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

(* #MapQ ******************************************************************************************)
MapQ[x_] = False;

(* #WithFilter ************************************************************************************)
WithFilter[surf_?SurfaceQ, filt_] := CorticalSurface[Normal[surf], Filter -> filt];
WithFilter[map_?MapQ, filt_] := With[
  {newFilt = Replace[
     Filter[map],
     {None -> filt,
      f_ :> (And[f[#1,#2], filt[#1,#2]]&)}]},
  SurfaceProjection[Domain[map], Duplicate -> map, Filter -> newFilt]];

(* SurfaceSelect **********************************************************************************)
SurfaceSelect[surf_?SurfaceQ, filterFn_] := With[
  {V = Vertices[surf],
   Z = Field[surf],
   F = Faces[surf]},
  With[
    {pick = Map[filterFn, Normal[surf]]},
    With[
      {idcs = Flatten[Position[pick, True]]},
      With[
        {tr = Dispatch[MapThread[Rule, {idcs, Range[Length[idcs]]}]]},
        CorticalSurface[
          Pick[V, pick],
          Field -> Pick[Z, pick],
          Faces -> (Select[F, Apply[And, pick[[#]]]&] /. tr)]]]]];

(* SurfaceSelect **********************************************************************************)
SurfaceCases[surf_?SurfaceQ, form_] := With[
  {V = Vertices[surf],
   Z = Field[surf],
   F = Faces[surf]},
  With[
    {pick = Map[MatchQ[#, form]&, Normal[surf]]},
    With[
      {idcs = Flatten[Position[pick, True]]},
      With[
        {tr = Dispatch[MapThread[Rule, {idcs, Range[Length[idcs]]}]]},
        CorticalSurface[
          Pick[V, pick],
          Field -> Pick[Z, pick],
          Faces -> (Select[F, Apply[And, pick[[#]]]&] /. tr)]]]]];

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
  Options[Graphics3D],
  {ColorFunction -> Automatic,
   ColorFunctionScaling -> True}];
SurfacePlot[surf_?SurfaceQ, opts:OptionsPattern[]] := Graphics3D[
  With[
    {V = Normal[Vertices[surf]],
     F = Normal[Faces[surf]],
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
    {opts},
    Except[ColorFunction|ColorFunctionScaling]]];

(* #SurfaceMap ************************************************************************************)
SurfaceMap[fn_, surf_?SurfaceQ] := CorticalSurface[
  Map[fn, Normal[surf]],
  Faces -> Faces[surf]];
SurfaceMap[fn_, surfs_List /; And@@((ListQ[#] || SurfaceQ[#])& /@ surfs)] := CorticalSurface[
  Map[Apply[fn, #]&, Transpose[Normal/@surfs]],
  Faces -> Faces[surfs[[1]]]];

(* #SurfaceReplace ********************************************************************************)
SurfaceReplace[surf_?SurfaceQ, rules_] := CorticalSurface[
  Replace[Normal[surf], rules],
  Faces -> Faces[surf]];

(* #SurfaceResample *******************************************************************************)
Options[SurfaceResample] = Prepend[
  FilterRules[Options[CorticalSurface], Except[Field]],
  Method -> Nearest];
SurfaceResample[a_?SurfaceQ, b_?SurfaceQ, opts:OptionsPattern[]] := Catch[
  If[Length[Field[a]] == Length[Field[b]] && Vertices[a] == Vertices[b],
    WithField[b, Field[a]],
    CorticalSurface[
      Vertices[b],
      Sequence@@FilterRules[{opts}, Except[Method]],
      Faces -> Faces[b],
      Field -> With[
        {method = OptionValue[Method]},
        Switch[
          method,
          Nearest, With[
            {nearest = Nearest[Vertices[a] -> Automatic],
             F = Field[a]},
            Map[Function[F[[First[nearest[#]]]]], Vertices[b]]],
          (Interpolation|{Interpolation, ___}), With[
            {interp = Interpolation[
               MapThread[List, {Vertices[a], Field[a]}],
               Sequence@@If[ListQ[method], Rest[methods], {}]]},
            Map[interp, Vertices[b]]],
          (ListInterpolation|{ListInterpolation, ___}), With[
            {interp = ListInterpolation[
               MapThread[List, {Vertices[a], Field[a]}],
               Sequence@@If[ListQ[method], Rest[methods], {}]]},
            Map[interp, Vertices[b]]],
          _, (
            Message[
              SurfResample::badarg,
              Method,
              "must be Nearest/Interpolation/ListInterpolation"];
            Throw[$Failed])]]]]];
(* Alias for surface resample that self-hashes like surfaces *)
(*Unprotect[Rule];
Rule[a_?SurfaceQ, b_?SurfaceQ] := With[
  {res = SurfaceResample[a, b]},
  If[res =!= $Failed, 
    Unprotect[Rule];
    Set[Rule[a,b], res];
    Protect[Rule]];
  res];
Protect[Rule];*)

(* #WithField and Machinery for specifying Cortical Surfaces as field -> surface ******************)
WithField[s_?SurfaceQ, f_?FieldQ] := Rule[f, s];
Vertices[Rule[f_?FieldQ, s_?SurfaceQ]] := Vertices[s];
Field[Rule[f_?FieldQ, s_?SurfaceQ]] := Field[f];
Faces[Rule[f_?FieldQ, s_?SurfaceQ]] := Faces[s];
Polygons[Rule[f_?FieldQ, s_?SurfaceQ]] := Polygons[s];
SurfaceQ[Rule[f_?FieldQ, s_?SurfaceQ]] := True;

WithField[s_?SurfaceQ, f_List] := Rule[f, s];
Vertices[Rule[f_List, s_?SurfaceQ]] := Vertices[s];
Field[Rule[f_List, s_?SurfaceQ]] := f;
Faces[Rule[f_List, s_?SurfaceQ]] := Faces[s];
Polygons[Rule[f_List, s_?SurfaceQ]] := Polygons[s];
SurfaceQ[Rule[f_List, s_?SurfaceQ]] := True;

WithField[m_?MapQ, f_?FieldQ] := Rule[f, m];
Vertices[Rule[f_?FieldQ, m_?MapQ]] := Vertices[m];
Filter[Rule[f_?FieldQ, m_?MapQ]] := Filter[m];
Radius[Rule[f_?FieldQ, m_?MapQ]] := Radius[m];
Surface[Rule[f_?FieldQ, m_?MapQ]] := Surface[m];
OrientPoint[Rule[f_?FieldQ, m_?MapQ]] := m;
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

WithField[m_?MapQ, f_List] := Rule[f, m];
Vertices[Rule[f_List, m_?MapQ]] := Vertices[m];
Filter[Rule[f_List, m_?MapQ]] := Filter[m];
Radius[Rule[f_List, m_?MapQ]] := Radius[m];
Surface[Rule[f_List, m_?MapQ]] := Surface[m];
OrientPoint[Rule[f_List, m_?MapQ]] := m;
InverseProjectionDispatch[Rule[f_List, m_?MapQ]] := InverseProjectionDispatch[m];
InverseProjectionTransform[Rule[f_List, m_?MapQ]] := InverseProjectionTransform[m];
ProjectionDispatch[Rule[f_List, m_?MapQ]] := ProjectionDispatch[m];
ProjectionShear[Rule[f_List, m_?MapQ]] := ProjectionShear[m];
SurfaceRotation[Rule[f_List, m_?MapQ]] := SurfaceRotation[m];
ProjectionTransform[Rule[f_List, m_?MapQ]] := ProjectionTransform[m];
DomainIndices[Rule[f_List, m_?MapQ]] := DomainIndices[m];
Domain[Rule[f_List, m_?MapQ]] := Domain[m];
Field[Rule[f_List, m_?MapQ]] /; Length[f] == Length[Vertices[m]] := f;
Field[Rule[f_List, m_?MapQ]] /; Length[f] == Length[Vertices[Surface[m]]] := f[[DomainIndices[m]]];
Faces[Rule[f_List, m_?MapQ]] := Faces[m];
MapQ[Rule[f_List, m_?MapQ]] := True;

Unprotect[Rule];

Rule /: Center[Rule[f_?FieldQ, m_?MapQ]] := Center[m];
Rule /: Center[Rule[a_List, m_?MapQ]] := Center[m];

Rule /: Normal[Rule[f_?FieldQ, s_?SurfaceQ]] := MapThread[Rule, {Vertices[s], Field[f]}];
Rule /: Normal[Rule[f_?FieldQ, m_?MapQ]] := MapThread[
  Rule,
  {Vertices[m], Part[Field[f], DomainIndices[m]]}];
Rule /: Normal[Rule[a_List, s_?SurfaceQ]] := MapThread[Rule, {Vertices[s], a}];
Rule /: Normal[Rule[a_List, m_?MapQ]] /; Length[a] == Length[DomainIndices[m]] := MapThread[
  Rule,
  {Vertices[m], a}];
Rule /: Normal[Rule[a_List, m_?MapQ]] /; Length[a] == Length[Vertices[Surface[m]]] := MapThread[
  Rule,
  {Vertices[m], a[[DomainIndices[m]]]}];

Rule[f_?FieldQ, s_?SurfaceQ] /; Length[Field[f]] != Length[Field[s]] := (
  Message[WithField::incompat];
  $Failed);
Rule[f_?FieldQ, m_?MapQ] /; Length[Field[f]] != Length[Field[Surface[m]]] := (
  Message[WithField::incompat];
  $Failed);
Rule[a_List, b_?SurfaceQ] /; Length[a] != Length[Vertices[b]] := (
  Message[WithField::incompat];
  $Failed);
Rule[a_List, b_?MapQ] := (
  Message[WithField::incompat];
  $Failed
 ) /; (Length[a] != Length[Vertices[b]] && Length[a] != Length[Vertices[Surface[b]]]);

Protect[Rule];

ToField[dat_List] := (
  Unprotect[Field];
  Field[dat] = dat;
  Protect[Field];
  dat);

Protect[Azimuth, Cartesian, CartesianToSpherical, ConvertCoordinates,
        CorticalSurface, CorticalSurfaceFromVTK, Domain,
        DomainIndices, Duplicate, Faces, Field, FieldQ, Filter,
        InverseProjectionDispatch, InverseProjectionTransform,
        Latitude, Longitude, MapPlot, MapQ, MergeSurfaces,
        OrientMatrix, OrientPoint, PolarAngle, Polygons,
        ProjectionDispatch, ProjectionRotation, ProjectionShear,
        ProjectionTransform, Radius, ReadVTK,
        SphericalCoordinateStyle, SphericalToCartesian, Surface,
        SurfaceCases, SurfaceMap, SurfacePlot, SurfaceProjection,
        SurfaceQ, SurfaceReplace, SurfaceRotation, SurfaceResample,
        SurfaceSelect, ToField, Vertices, WithField, WithFilter];

End[];
EndPackage[];
