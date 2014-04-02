(* OccipitalPole.m
 *
 * Basic utility functions for dealing with cortical surface maps of the occipital pole in 
 * Mathematica.
 *
 * Copyright (C) 2013-2014 by Noah C. Benson.
 * This file is part of the MmaSurfer library, which is provided under the terms of the Eclipse
 * Public License, version 1.0. See the accompanying LICENSE file for more information.
 *)

(**************************************************************************************************)
BeginPackage["OccipitalPole`", {"CorticalSurface`", "FreeSurfer`", "ComputationalGeometry`"}];

Unprotect[FSAverageV1Mask, FSAverageSymV1Mask, FSAverageV1, FSAverageSymV1];
ClearAll[FSAverageV1Mask, FSAverageSymV1Mask, FSAverageV1, FSAverageSymV1];
FSAVerageMask::usage = "FSAverageV1[hemi] is a lazily-loaded label object of the points in the
 fsaverage V1, according to Hinds et al. (2008) NeuroImage. This data is loaded from the
 fsaverage subject, so that must be set appropriately (in the FreeSurfer package).";
FSAVerageSymMask::usage = "FSAverageSymV1 is a lazily-loaded label object of the points in the
 fsaverage_sym V1, according to Hinds et al. (2008) NeuroImage. This data is loaded from
 the fsaverage_sym subject, so that must be set appropriately (in the FreeSurfer package).";
FSAVerageV1::usage = "FSAverageV1[hemi] is a lazily-loaded list of vertex indices of those points
 in the fsaverage V1, according to Hinds et al. (2008) NeuroImage. This data is loaded from
 the fsaverage subject, so that must be set appropriately (in the FreeSurfer package). Note
 that the hull is determined for the surface projection in which the occipital pole is the
 center of the projection; indices are given, and the final and first indices are the
 same.";
FSAVerageSymV1::usage = "FSAverageSymV1 is a lazily-loaded list of vertex indices of those points
 in the fsaverage_sym V1, according to Hinds et al. (2008) NeuroImage. This data is loaded
 from the fsaverage_sym subject, so that must be set appropriately (in the FreeSurfer
 package). Note that the hull is determined for the surface projection in which the
 occipital pole is the center of the projection; indices are given, and the final and first
 indices are the same.";

Unprotect[FSAverageV1Hull, FSAverageSymV1Hull];
ClearAll[FSAverageV1Hull, FSAverageSymV1Hull];
FSAverageV1Hull::usage = "$FSAverageV1Hull[hemi] is the list of polygons in the convex hull of
 $FSAverageV1 (with the first polygon appended so that ListPlot or Lines will produce a
 correct polygon).";
FSAverageSymV1Hull::usage = "$FSAverageSymV1Hull is the list of polygons in the convex hull of
 $FSAverageSymV1 (with the first polygon appended so that ListPlot or Lines will produce
 a correct polygon).";

Unprotect[FSAverageOP, FSAveragePialOP, FSAverageInflatedOP, FSAverageSphereOP, 
          FSAverageSymOP, FSAverageSymPialOP, FSAverageSymInflatedOP, FSAverageSymSphereOP, 
          FSAverageSymRegisteredOP];
ClearAll[ FSAverageOP, FSAveragePialOP, FSAverageInflatedOP, FSAverageSphereOP, 
          FSAverageSymOP, FSAverageSymPialOP, FSAverageSymInflatedOP, FSAverageSymSphereOP, 
          FSAverageSymRegisteredOP];
FSAverageOP::usage = "FSAverageOP[hemi] is the (lazily-evaluated) index of the vertex in the
 occipital cortex (spherical hemisphere) that corresponds roughly to the tip of the
 occipital pole in the fsaverage.";
FSAveragePialOP::usage = "FSAveragePialOP[hemi] yields the vertex coordinate (not the index) of the
  occipital pole in the given hemisphere hemi of the fsaverage brain's pial surface.";
FSAverageInflatedOP::usage = "FSAverageInflatedOP[hemi] yields the vertex coordinate (not the index)
  of the occipital pole in the given hemisphere hemi of the fsaverage brain's inflated surface.";
FSAverageSphereOP::usage = "FSAverageSphereOP[hemi] yields the vertex coordinate (not the index) of
  the occipital pole in the given hemisphere hemi of the fsaverage brain's spherical surface.";
FSAverageSymOP::usage = "FSAverageSymOP is the (lazily-evaluated) index of the vertex in the
 occipital cortex that corresponds roughly to the tip of the occipital pole in the
 fsaverage_sym brain.";
FSAverageSymPialOP::usage = "FSAverageSymPialOP[hemi] yields the vertex coordinate (not the index)
 of the occipital pole in the given hemisphere hemi of the fsaverage brain's pial surface.";
FSAverageSymInflatedOP::usage = "FSAverageSymInflatedOP[hemi] yields the vertex coordinate (not the
 index) of the occipital pole in the given hemisphere hemi of the fsaverage brain's inflated
 surface.";
FSAverageSymSphereOP::usage = "FSAverageSymSphereOP[hemi] yields the vertex coordinate (not the
 index) of the occipital pole in the given hemisphere hemi of the fsaverage brain's spherical
 surface.";

Unprotect[SubjectOP, SubjectPialOP, SubjectInflatedOP, SubjectSphereOP, SubjectRegisteredOP];
ClearAll[ SubjectOP, SubjectPialOP, SubjectInflatedOP, SubjectSphereOP, SubjectRegisteredOP];
SubjectOP::usage = "SubjectOP[sub, hemi] yields the vertex index of the subject sub's occipital
 pole in the given hemisphere hemi.";
SubjectPialOP::usage = "SubjectPialOP[sub, hemi] yields the vertex coordinate (not the index) of the
 occipital pole in the given subject sub's given pial hemisphere hemi.";
SubjectInflatedOP::usage = "SubjectInflatedOP[sub, hemi] yields the vertex coordinate (not the
 index) of the occipital pole in the given subject sub's given inflated hemisphere hemi.";
SubjectSphereOP::usage = "SubjectSphereOP[sub, hemi] yields the vertex coordinate (not the index) of
 the occipital pole in the given subject sub's given spherical hemisphere hemi.";
SubjectRegisteredOP::usage = "SubjectPialOP[sub, hemi] yields the vertex coordinate (not the index)
 of the occipital pole in the given subject sub's given registered spherical hemisphere hemi.";


(**************************************************************************************************)
Begin["`Private`"];

FSAverageV1Mask[hemi:(LH|RH)] := With[
  {res = Check[
     Import[
       $FSAverage <> "/label/"<> If[hemi === LH, "lh", "rh"] <> ".v1.predict.label",
       "FreeSurferLabel",
       Max -> Length[Field[FSAverageCurvature[hemi]]]],
     $Failed]},
  If[res === $Failed,
    $Failed,
    Set[
      FSAverageV1Mask[hemi],
      res]]];
FSAverageSymV1Mask := With[
  {res = Check[
     Import[
       If[FileExistsQ[$FSAverageSym <> "/label/lh.v1.predict.label"],
         $FSAverageSym <> "/label/lh.v1.predict.label",
         $FSAverageSym <> "/label/lh.v1-predict.label"],
       "FreeSurferLabel",
       Max -> Length[Field[FSAverageSymCurvature]]],
     $Failed]},
  If[res === $Failed,
    $Failed,
    Set[
      FSAverageSymV1Mask,
      res]]];

FSAverageV1[hemi:(LH|RH)] := With[
  {idcs = Check[
     Flatten[Position[Normal[Field[FSAverageV1Mask[hemi]]], 1]],
     $Failed]},
  If[idcs === $Failed,
    $Failed,
    Set[
      FSAverageV1[hemi],
      idcs]]];
FSAverageSymV1 := With[
  {idcs = Check[
     Flatten[Position[Normal[Field[FSAverageSymV1Mask]], 1]],
     $Failed]},
  If[idcs === $Failed,
    $Failed,
    Set[
      FSAverageSymV1,
      idcs]]];

FSAverageOP[hemi:(LH|RH)] := With[
  {V = Check[Vertices[FSAverageInflatedSurface[hemi]], $Failed]},
  If[V === $Failed,
    $Failed,
    Set[
      FSAverageOP[hemi],
      First[Ordering[V[[All, 2]], 1]]]]];
FSAverageSymOP := With[
  {V = Check[Vertices[FSAverageSymInflatedSurface], $Failed]},
  If[V === $Failed,
    $Failed,
    Set[
      FSAverageSymOP,
      First[Ordering[V[[All, 2]], 1]]]]];

FSAverageV1Hull[hemi:(LH|RH)] := With[
  {res = Check[
     With[
       {Z = Normal[Field[FSAverageV1Mask[hemi]]],
        V = Vertices[FSAverageSphereSurface[hemi]]},
       With[
         {ids = Map[
            Rest[SortBy[#, Z[[#]] &]] &,
            Select[
              Faces[FSAverageSphereSurface[hemi]],
              Total[Z[[#]]] == 2 &]]},
         With[
           {disp = Dispatch[
              Map[
                #[[1, 1]] -> #[[All, 2]] &,
                GatherBy[
                  Flatten[
                    Map[
                      {#1[[1]] -> #1[[2]], #1[[2]] -> #1[[1]]} &,
                      ids]],
                  #[[1]] &]]],
            first = ids[[1, 2]],
            last = ids[[1, 1]]},
           NestWhileList[
             Function[{edge},
               With[
                 {from = edge[[1]],
                  opts = edge[[2]] /. disp},
                 {edge[[2]], If[opts[[1]] == from, opts[[2]], opts[[1]]]}]],
             First[ids],
             #[[2]] != last &]]]],
     $Failed]},
  If[res === $Failed, $Failed, Set[FSAverageV1Hull[hemi], res]]];
FSAverageSymV1Hull := With[
  {res = Check[
     With[
       {Z = Normal[Field[FSAverageSymV1Mask]],
        V = Vertices[FSAverageSymSphereSurface]},
       With[
         {ids = Map[
            Rest[SortBy[#, Z[[#]] &]] &,
            Select[
              Faces[FSAverageSymSphereSurface],
              Total[Z[[#]]] == 2 &]]},
         With[
           {disp = Dispatch[
              Map[
                #[[1, 1]] -> #[[All, 2]] &,
                GatherBy[
                  Flatten[
                    Map[
                      {#1[[1]] -> #1[[2]], #1[[2]] -> #1[[1]]} &,
                      ids]],
                  #[[1]] &]]],
            first = ids[[1, 2]],
            last = ids[[1, 1]]},
           NestWhileList[
             Function[{edge},
               With[
                 {from = edge[[1]],
                  opts = edge[[2]] /. disp},
                 {edge[[2]], If[opts[[1]] == from, opts[[2]], opts[[1]]]}]],
             First[ids],
             #[[2]] != last &]]]],
     $Failed]},
  If[res === $Failed, $Failed, Set[FSAverageSymV1Hull, res]]];


SubjectOP[sub_, hemi:(LH|RH)] := With[
  {V = Check[Vertices[SubjectInflatedSurface[sub, hemi]], $Failed]},
  If[V === $Failed,
    $Failed,
    Set[
      SubjectOP[sub, hemi],
      First[Ordering[V[[All, 2]], 1]]]]];
SubjectPialeOP[sub_, hemi:(LH|RH)] := With[
  {idx = SubjectOP[sub, hemi]},
  If[idx === $Failed, $Failed, Vertices[SubjectPialSurface[sub, hemi]][[idx]]]];
SubjectInflatedOP[sub_, hemi:(LH|RH)] := With[
  {idx = SubjectOP[sub, hemi]},
  If[idx === $Failed, $Failed, Vertices[SubjectInflatedSurface[sub, hemi]][[idx]]]];
SubjectSphereOP[sub_, hemi:(LH|RH)] := With[
  {idx = SubjectOP[sub, hemi]},
  If[idx === $Failed, $Failed, Vertices[SubjectSphereSurface[sub, hemi]][[idx]]]];
SubjectRegisteredOP[sub_, hemi:(LH|RH)] := With[
  {idx = SubjectOP[sub, hemi]},
  If[idx === $Failed, $Failed, Vertices[SubjectRegisteredSurface[sub, hemi]][[idx]]]];


End[];
EndPackage[];

