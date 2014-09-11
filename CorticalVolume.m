(* CorticalVolume.m
 *
 * Utility functions for dealing with (mostly FreeSurfer) cortical volume data in Mathematica.
 *
 * Copyright (C) 2013-2014 by Noah C. Benson.
 * This file is part of the MmaSurfer library, which is provided under the terms of the Eclipse
 * Public License, version 1.0. See the accompanying LICENSE file for more information.
 *)

(**************************************************************************************************)
BeginPackage["CorticalVolume`"];
Unprotect["CorticalVolume`*", "CorticalVolume`Private`*"];
ClearAll["CorticalVolume`*", "CorticalVolume`Private`*"];

VolumeShell::usage = "VolumeShell[vol, tag] yields the indices of the outer shell of voxels in the volume vol whose values match the pattern tag when they are in the volume.";

VolumeIndices::usage = "VolumeIndices[vol, tag] yields the indices of all voxels in the volume vol that match the given tag.";
VolumeMask::usage = "VolumeMask[indices] yields a volume mask in which voxels labeled in indices are given a value of 1 and other voxels are given a value of 0. An optional third argument, size, may give the dimensions (e.g., {256, 256, 256}); otherwise, the largest value in the indices is used. The options True and False may also be passed to specify what values should be given to voxels in or out of the mask, respectively.";

VolumeBoundaries::usage = "VolumeBoundaries[vol, tag] yields a list of {{xmin, xmax}, {ymin, ymax}, ...} such that everything outside of the given rectangular ranges fails to match tag. If no voxels match the given tag, then None is returned.";

SurfacePointToVoxelIndex::usage = "SurfacePointToVoxelIndex[p, voldims] gives a translation of the surface point p ({x, y, z}) to an index ({i, j, k}) such that rounding the index values will give an approximate position of the surface point in the volume with dims given in voldims. If not provided, then voldims is taken to ba {256, 256, 256}. The first argument may also be a list of points, in which case the equivalent list of indices is returned.
Note that this method is intended to work with typically processed FreeSurfer volumes and surfaces and is not designed for other programs or for volumes that are used as input to FreeSurfer.";


(**************************************************************************************************)
Begin["`Private`"];

(* #VolumeShell ***********************************************************************************)
With[
  {neis = Flatten[
     Outer[
       Plus,
       {{1, 0, 0}, {-1, 0, 0}, {0, 0, 0}},
       {{0, 1, 0}, {0, -1, 0}, {0, 0, 0}},
       {{0, 0, 1}, {0, 0, -1}, {0, 0, 0}},
       1, 1, 1],
     2]},
  VolumeShell[mask_List /; Length[Dimensions[mask]] > 2, tag_] := Union[
     Last@Reap[
       Scan[
         Function[{idx}, 
           Scan[
             Function[{neighborAt},
               With[
                 {neighborIdx = idx + neighborAt},
                 If[Not[MatchQ[Extract[mask, neighborIdx], tag]],
                   (Sow[neighborIdx]; Sow[idx])]]],
           neis]],
         Position[mask, tag, {3}, Heads -> False]],
       _,
       Apply[Sequence, #2]&]]];
Protect[VolumeShell];

(* #VolumeIndices *********************************************************************************)
VolumeIndices[vol_List /; Length[Dimensions[vol]] > 2, tag_] := Position[
  vol,
  tag,
  {Length[Dimensions[vol]]},
  Heads -> False];
VolumeIndices[vol_MGHQ, tag_] := Map[
  Function[{vol}, Position[vol, tag, {3}, Heads -> False]],
  Volumes[vol]];
Protect[VolumeIndices];

(* #VolumeBoundaries ******************************************************************************)
VolumeBoundaries[vol_List /; Length[Dimensions[vol]] > 2, tag_] := With[
  {match = VolumeIndices[vol, tag]},
  If[Length[match] == 0,
    None,
    Table[{Min[match[[All, k]]], Max[match[[All, k]]]}, {k, 1, Length[First[match]]}]]];

(* #VolumeMask ************************************************************************************)
Options[VolumeMask] = {True -> 1, False -> 0};
VolumeMask[idcs:{{_Integer, _Integer, _Integer}..}, 
           dims:{_Integer, _Integer, _Integer}, 
           OptionsPattern[]] := Check[
  With[
    {tval = OptionValue[True],
     fval = OptionValue[False]},
    SparseArray[
      Map[(# -> tval)&, idcs],
      dims,
      fval]],
  $Failed];
VolumeMask[idcs:{{_Integer, _Integer, _Integer}..}, OptionsPattern[]] := VolumeMask[
  idcs,
  Max/@Transpose[idcs],
  True -> OptionValue[True],
  False -> OptionValue[False]];
Protect[VolumeMask];

(* #SurfacePointToVoxelIndex **********************************************************************)
SurfacePointToVoxelIndex[p_List, voxdim_List] := Which[
  Length[p] == 0, p,
  ListQ[First[p]], Map[SurfacePointToVoxelIndex[#, voxdim]&, p],
  True, List[
    0.5*voxdim[[1]] - p[[1]] + 1.0, 
    0.5*voxdim[[2]] - p[[3]] + 1.0, 
    0.5*voxdim[[3]] + p[[2]] + 1.0]];
SurfacePointToVoxelIndex[p_List] := SurfacePointToVoxelIndex[p, {256.0, 256.0, 256.0}];
Protect[SurfacePointToVoxelIndex];

End[];
EndPackage[];

