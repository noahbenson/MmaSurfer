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

Unprotect[VolumeShell];
ClearAll[VolumeShell];
VolumeShell::usage = "VolumeShell[vol, tag] yields the indices of the outer shell of voxels in the
 volume vol whose values match the pattern tag when they are in the volume.";

Unprotect[VolumeIndices, VolumeMask];
ClearAll[VolumeIndices, VolumeMask];
VolumeIndices::usage = "VolumeIndices[vol, tag] yields the indices of all voxels in the volume vol
 that match the given tag.";
VolumeMask::usage = "VolumeMask[indices] yields a volume mask in which voxels labeled in indices
 are given a value of 1 and other voxels are given a value of 0. An optional third
 argument, size, may give the dimensions (e.g., {256, 256, 256}); otherwise, the largest
 value in the indices is used. The options True and False may also be passed to specify
 what values should be given to voxels in or out of the mask, respectively.";


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
     Reap[
       Scan[
        Function[{idx}, 
          Scan[
            Function[{neighborAt},
              With[
                {neighborIdx = idx + neighborAt},
                If[Not[MatchQ[Extract[mask, neighborIdx], tag]],
                  (Sow[neighborIdx]; Sow[idx])]]],
          neis]],
        Position[mask, tag, {3}]]
      ][[2, 1]]]];
Protect[VolumeShell];

(* #VolumeIndices *********************************************************************************)
VolumeIndices[vol_List /; Length[Dimensions[vol]] > 2, tag_] := Position[
  vol,
  tag,
  {3}];
Protect[VolumeIndices];

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

End[];
EndPackage[];

