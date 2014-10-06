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

CorticalVolume::usage = "CorticalVolume[data] yields a cortical volume object for the given matrix data. The data argument may be either an n x m x q matrix or a list of such identically sized matrices (in which case the elements of data are considered frames).
The following functions may be used on a cortical volume object: VolumeQ, VoxelToCoordinateMatrix, CoordinateToVoxelMatrix, OrientationMatrix, InverseOrientationMatrix, Frames, DataRange, Spacings, Normal.";
CorticalVolume::badarg = "Bad argument given to CorticalVolume: `1`";

VolumeQ::usage = "VolumeQ[x] yields True if and only if x is a cortical volume object; otherwise, yields False.";

VoxelToCoordinateMatrix::usage = "VoxelToCoordinateMatrix[vol, coordinate] yields the 3 x 4 matrix that can be used to translate voxel indices to coordinates.";

VoxelToCoordinate::usage = "VoxelToCoordinate[vol, coordinate] yields a translation of the given voxel index to an (x,y,z) coordinate, according to the volume vol's orientation matrix.";
VoxelToCoordinate::badarg = "Bad argument given to VoxelToCoordinate: `1`";

CoordinateToVoxelMatrix::usage = "CoordinateToVoxelMatrix[vol, coordinate] yields the 3 x 4 matrix that can be used to translate coordinates to voxel indices.";

CoordinateToVoxel::usage = "CoordinateToVoxel[vol, coordinate] yields a translation of the given (x,y,z) coordinate to a position in voxel space, according to the volume vol's orientation matrix.";
CoordinateToVoxel::badarg = "Bad argument given to CoordinateToVoxel: `1`";

OrientationMatrix::usage = "OrientationMatrix is an option to CortivalVolume that specifies the orientation of the given volume. The matrix should be a 3 x 4 matrix or Automatic (in which case the default matrix is {{-1,0,0,0},{0,0,1,0},{0,-1,0,0}}) in which the final column is the displacement.
OrientationMatrix[vol] yields the orientation matrix for the given cortical volume object vol.";

InverseOrientationMatrix::usage = "InverseOrientationMatrix[vol] yields the inverse of the orientation matrix for the given cortical volume object vol.";


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
SurfacePointToVoxelIndex[p_List, voxdim:{_NumericQ, _NumericQ, _NumericQ}] := Which[
  Length[p] == 0, p,
  ListQ[First[p]], Map[SurfacePointToVoxelIndex[#, voxdim]&, p],
  True, List[
    0.5*voxdim[[1]] - p[[1]] + 1.0, 
    0.5*voxdim[[2]] - p[[3]] + 1.0, 
    0.5*voxdim[[3]] + p[[2]] + 1.0]];
SurfacePointToVoxelIndex[p_List, vox_List /; ArrayQ[vox, 3]] := SurfacePointToVoxelIndex[
  p,
  Dimensions[vox]];
SurfacePointToVoxelIndex[p_List] := SurfacePointToVoxelIndex[p, {256.0, 256.0, 256.0}];
Protect[SurfacePointToVoxelIndex];

(* #CorticalVolume ********************************************************************************)
Options[CorticalVolume] = {
  Spacings -> {1.0, 1.0, 1.0},
  OrientationMatrix -> Automatic};
CorticalVolume[data_ /; ArrayQ[data, 3|4], OptionsPattern[]] := Catch[
  With[
    {spacing = Replace[
       OptionValue[Spacings],
       {l:{_?NumericQ, _?NumericQ, _?NumericQ} /; l == Re[l] && l == Abs[l] :> l,
        _ :> (
          Message[
            CorticalVolume::badarg,
            "Spacings must be a 3-element list of positive real numbers"];
          Throw[$Failed])}],
     matrix = Replace[
       OptionValue[OrientationMatrix],
       {m_List /; MatrixQ[m, (NumericQ[#] && Re[#] == #)&] && Dimensions[m] == {3,4} :> m,
        Automatic -> {{-1,0,0,0},{0,0,1,0},{0,-1,0,0}},
        _ :> (
          Message[
            CorticalVolume::badarg,
            "OrientationMatrix must be a 3 x 4 matrix of real numbers"];
          Throw[$Failed])}],
     sym = Unique["volume"],
     dims = Dimensions[data]},
    SetAttributes[Evaluate[sym], Temporary];
    sym /: VolumeQ[sym] = True;
    sym /: Spacings[sym] = spacing;
    sym /: OrientationMatrix[sym] = matrix;
    sym /: InverseOrientationMatrix[sym] = Most@Inverse[Append[matrix, {0,0,0,1}]];
    sym /: Normal[sym] = data;
    sym /: Dimensions[sym] = dims[[1;;3]];
    sym /: Frames[sym] = If[Length[dims] == 4, Last[dims], 1];
    sym /: DataRange[sym] = With[
      {r = dims * 0.5 * spacing},
      Transpose[{-r,r}]];
    sym /: VoxelToCoordinateMatrix[sym] = With[
      {mtx = MapThread[Most[#1]*#2 &, {matrix, spacing}]},
      MapThread[Append, {mtx, -Dot[mtx, 0.5*dims]}]];
    sym /: CoordinateToVoxelMatrix[sym] = Most[
      Inverse[Append[VoxelToCoordinateMatrix[sym], {0,0,0,1}]]];
    sym]];
Protect[CorticalVolume];

VolumeQ[v_] := False;
Protect[VolumeQ];

(* #VoxelToCoordinate *****************************************************************************)
VoxelToCoordinate[vol_?VolumeQ, idx:{_?NumericQ, _?NumericQ, _?NumericQ}] := With[
  {dims = Dimensions[vol],
   matrix = VoxelToCoordinateMatrix[vol]},
  If[And@@MapThread[TrueQ[0 < #1 <= #2]&, {idx, dims}],
    Dot[matrix, Append[(idx - 1.0) * (dims / (dims - 1.0)), 1.0]],
    (Message[VoxelToCoordinate::badarg, "index is not a valid voxel"];
     $Failed)]];
VoxelToCoordinate[vol_?VolumeQ, idcs:{{_?NumericQ, _?NumericQ, _?NumericQ}..}] := With[
  {dims = Dimensions[vol],
   matrix = VoxelToCoordinateMatrix[vol]},
  Map[
    Function[
      If[And@@MapThread[TrueQ[0 < #1 <= #2]&, {#, dims}],
        Dot[matrix, Append[(# - 1.0) * (dims / (dims - 1.0)), 1.0]],
        (Message[VoxelToCoordinate::badarg, "index is not a valid voxel"];
         $Failed)]],
    idcs]];
Protect[VoxelToCoordinate];

(* #CoordinateToVoxel *****************************************************************************)
CoordinateToVoxel[vol_?VolumeQ, coord:{_?NumericQ, _?NumericQ, _?NumericQ}] := With[
  {dims = Dimensions[vol]},
  With[
    {idx = Dot[CoordinateToVoxelMatrix[vol], Append[coord, 1.0]] * (dims - 1.0)/dims + 1.0},
    If[!And@@MapThread[TrueQ[0 < #1 <= #2]&, {idx, dims}],
      Message[CoordinateToVoxel::badarg, "coordinate is outside voxel range"]];
    idx]];
CoordinateToVoxel[vol_?VolumeQ, coords:{{_?NumericQ, _?NumericQ, _?NumericQ}..}] := With[
  {dims = Dimensions[vol],
   matrix = CoordinateToVoxelMatrix[vol]},
  Map[
    Function[
      With[
        {coord = Dot[matrix, Append[#,1.0]] * (dims - 1.0)/dims + 1.0},
        If[!And@@MapThread[TrueQ[0 < #1 <= #2]&, {coord, dims}],
          Message[CoordinateToVoxel::badarg, "coordinate is outside voxel range"]];
        coord]],
    coords]];
Protect[CoordinateToVoxel];

End[];
EndPackage[];

