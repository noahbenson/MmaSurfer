(* VisualCortex.m
 *
 * Functions and data structures for efficiently implementing retinotopy on the cortical surface,
 * including the Banded Double-Sech model of Schira, Tyler, Spehar, and Breakspear (2010), in
 * Mathematica.
 * References:
 *   [Retinotopic Templates]
 *   Benson NC, Butt OH, Brainard DH, Aguirre GK (2014) Correction of Distortion in Flattened 
 *     Representations of the Cortical Surface Allows Prediction of V1-V3 Functional Organization
 *     from Anatomy. PLoS Comput Biol 10(3): e1003538. doi: 10.1371/journal.pcbi.1003538
 *   [Schira Model of Retinotopy]
 *   Schira MM, Tyler CW, Spehar B, Breakspear M (2010) Modeling Magnification and Anisotropy in the
 *     Primate Foveal Confluence. PLoS Comput Biol 6(1): e1000651. doi: 10.1371/journal.pcbi.1000651
 *   [V1 Masks and Hulls]
 *   Hinds OP, Rajendran N, Polimeni JR, Augustinack JC, Wiggins G, Wald LL, Diana Rosas H,
 *     Potthast A, Schwartz EL, Fischl B (2008) Accurate prediction of V1 location from cortical
 *     folds in a surface coordinate system. Neuroimage. 2008 Feb 15;39(4):1585-99.
 *     doi: 10.1016/j.neuroimage.2007.10.033
 *
 * Copyright (C) 2013-2014 by Noah C. Benson.
 * This file is part of the MmaSurfer library, which is provided under the terms of the Eclipse
 * Public License, version 1.0. See the accompanying LICENSE file for more information.
 *)

(**************************************************************************************************)
BeginPackage["VisualCortex`", {"CorticalSurface`", "FreeSurfer`"}];
Unprotect["VisualCortex`*", "VisualCortex`Private`*"];
ClearAll[ "VisualCortex`*", "VisualCortex`Private`*"];

(* Eccentricity, PolarAngle, and VisualAreal Definitions ******************************************)
Eccentricity::usage = "Eccentricity is a key used by the visual cortex package to represent eccentricity, as measured in degrees of visual angle from the foveal confluence (center of the visual field).";
PolarAngle::usage = "PolarAngle is a key used by the visual cortex package to represent the polar angle, as measured in degrees of rotation about the foveal confluence (center of the visual field) from the upper to the lower vertical meridia.";
VisualArea::usage = "VisualArea is a key used by the retinotopy package to represent the visual area ID of a particular patch of cortex. See also VisualAreasData.";

(* Complex Representation of the Visual Field  ****************************************************)
VisualAngleToComplex::usage = "VisualAngleToComplex[polarAngle, eccentricity] yields an imaginary number that represents the visual field coordinate. This complex number is always of the form r * Exp[I * t] where t is between -Pi/2 and Pi/2.";
ComplexToVisualAngle::usage = "ComplexToVisualAngle[z] yields a {polarAngle, eccentricity} pair that represents the visual angle coordinates represented by the complex number z. This complex  number should be of the form r * Exp[I * t] where t is between -Pi/2 and Pi/2.";
ComplexToPolarAngle::usage = "ComplexToPolarAngle[z] yields the polar angle value that is represented by the complex number z. This complex number should be of the form r * Exp[I * t] where t is between -Pi/2 and Pi/2.";
ComplexToEccentricity::usage = "ComplexToEccentricity[z] yields the eccentricity value that is represented by the complex number z. This complex number should be of the form r * Exp[I * t] where t is between -Pi/2 and Pi/2.";

(* Visual Areas ***********************************************************************************)
VisualAreaQ::usage = "VisualAreaQ[area] yields true if and only if area is a valid visual area id.";
VisualAreaData::usage = "VisualAreaData[id] yields a list of data regarding the visual area whose ID is given. VisualAreaData[] yields a list of all known IDs and their data.";
VisualAreaName::usage = "VisualAreaName[id] yields the visual area name of the given visual area ID.";
VisualAreaSimplify::usage = "VisualAreaSimplify[id] yields a simplified version of id such that border IDs (e.g., V1/V2 border) are converted to the lower visual area, ventral and dorsal signatures are removed (all made dorsal), and anything outside of V1, V2, and V3 is converted to 0. Accordingly, this function always yields 0, 1, 2, or 3 for any valid visual area id.";
VisualAreaData::badarea = "An unknown visual area was given: `1`";

(* Retinotopic Templates **************************************************************************)
FSAverageSymPolarAngle::usage = "FSAverageSymPolarAngle yields the polar angle data for the fsaverage_sym subjectas a field on the fsaverage_sym hemisphere.";
FSAverageSymEccentricity::usage = "FSAverageSymEccentricity yields the eccentricity data for the  fsaverage_sym subject as a field on the fsaverage_sym hemisphere.";
FSAverageSymVisualArea::usage = "FSAverageSymVisualArea yields the visual area data for the  fsaverage_sym subject as a field on the fsaverage_sym hemisphere. Each vertex will either be labeled as 1, 2, 3 (for dorsal regions), or -1, -2, -3 (for ventral regions) or None.";
FSAverageSymRetinotopy::usage = "FSAverageSymRetinotopy yields the retinotopy field for the fsaverage_sym subject. The field for each vertex is {PA, E, A} where PA is the polar angle, E is the eccentricity, and A is the area.";

SubjectPolarAngle::usage = "SubjectPolarAngle[sub, hemi] yields the polar angle data for the subject sub and the hemisphere hemi as a field on the subject's registered fsaverage_sym hemisphere.";
SubjectEccentricity::usage = "SubjectEccentricity[sub, hemi] yields the eccentricity data for the  subject sub and the hemisphere hemi as a field on the subject's registered fsaverage_sym hemisphere.";
SubjectVisualArea::usage = "SubjectVisualArea[sub, hemi] yields the visual area data for the subject sub and the hemisphere hemi as a field on the subject's registered fsaverage_sym hemisphere. Each vertex will either be labeled as 1, 2, 3 (for dorsal regions), or -1, -2, -3 (for ventral regions) or None.";
SubjectRetinotopy::usage = "SubjectRetinotopy[sub, hemi] yields the retinotopy field for the given subject sub and the given hemisphere hemi. The field for each vertex is {PA, E, A} where PA is the polar angle, E is the eccentricity, and A is the area.";

(* V1 Masks and Hulls *****************************************************************************)
FSAverageV1Mask::usage = "FSAverageV1Mask[hemi] is a lazily-loaded label object of the points in the fsaverage V1, according to Hinds et al. (2008) NeuroImage. This data is loaded from the fsaverage subject, so that must be set appropriately (in the FreeSurfer package).";
FSAverageSymV1Mask::usage = "FSAverageSymV1Mask is a lazily-loaded label object of the points in the fsaverage_sym V1, according to Hinds et al. (2008) NeuroImage. This data is loaded from the fsaverage_sym subject, so that must be set appropriately (in the FreeSurfer package).";
FSAVerageV1::usage = "FSAverageV1[hemi] is a lazily-loaded list of vertex indices of those points in the fsaverage V1, according to Hinds et al. (2008) NeuroImage. This data is loaded from the fsaverage subject, so that must be set appropriately (in the FreeSurfer package). Note that the hull is determined for the surface projection in which the occipital pole is the center of the projection; indices are given, and the final and first indices are the same.";
FSAVerageSymV1::usage = "FSAverageSymV1 is a lazily-loaded list of vertex indices of those points in the fsaverage_sym V1, according to Hinds et al. (2008) NeuroImage. This data is loaded from the fsaverage_sym subject, so that must be set appropriately (in the FreeSurfer package). Note that the hull is determined for the surface projection in which the occipital pole is the center of the projection; indices are given, and the final and first indices are the same.";
FSAverageV1Hull::usage = "$FSAverageV1Hull[hemi] is the list of polygons in the convex hull of $FSAverageV1 (with the first polygon appended so that ListPlot or Lines will produce a correct polygon).";
FSAverageSymV1Hull::usage = "$FSAverageSymV1Hull is the list of polygons in the convex hull of $FSAverageSymV1 (with the first polygon appended so that ListPlot or Lines will produce a correct polygon).";

(* The Schira Model *******************************************************************************)
$DefaultSchiraA::usage = "The default value of the Schira model A parameter.";
$DefaultSchira\[CapitalLambda]::usage = "The default value of the Schira model lambda parameter.";
$DefaultSchira\[CapitalPsi]::usage = "The default rotation of the Schira model in radians.";
$DefaultSchiraV1Size::usage = "The default size of V1 in the Schira model.";
$DefaultSchiraV2Size::usage = "The default size of V2 in the Schira model.";
$DefaultSchiraV3Size::usage = "The default size of V3 in the Schira model.";
$DefaultSchiraHV4Size::usage = "The default size of the hV4 pseudo-area in the Schira model.";
$DefaultSchiraV3ASize::usage = "The default size of V3A-like pseudo-area in the Schira model.";
$DefaultSchiraFC::usage = "The default position of the foveal confluence in the Schira model.";
$DefaultSchiraScale::usage = "The default {x,y} scale of the Schira model.";
$DefaultSchiraShear::usage = "The default shear matrix of the Schira model.";
$SchiraParameters::usage = "$SchiraParameters yields a list of the default parameters for the Schira model. The list rules are delayed so any Block redefining the default Schira parameters will affect this list's values.";

A::usage="The A parameter of the Schira model.";
B::usage="The B parameter of the Schira model.";
\[CapitalPsi]::usage="The rotation parameter of the Schira model.";
V1Size::usage="The parameter of the Schira model that determines the size of V1.";
V2Size::usage="The parameter of the Schira model that determines the size of V2.";
V3Size::usage="The parameter of the Schira model that determines the size of V3.";
HV4Size::usage="The parameter of the Schira model that determines the size of the hV4 pseudo-area.";
V3ASize::usage="The parameter of the Schira model that determines the size of the V3A pseudo-area.";
FC::usage="The foveal confluence position parameter of the Schira model.";
Shear::usage="The shear matrix parameter of the Schira model.";

SchiraModel::usage = "SchiraModel[parameters...] yields a SchiraModelObject using the given
 parameters when possible and using the default Schira parameters otherwise.";
SchiraModel::badarg = "Bad argument(s); all arguments must be rules: `1`";
SchiraModel::badarea = "Bad area given to Schira function; areas must be 1, 2, 3, or 4.";
SchiraModelObject::usage = "A SchiraModelObject form stores the data for a Schira model.";
SchiraModelObject::badarg = "Unrecognized SchiraModelObject argument: `1`";
SchiraFunction::usage = "SchiraFunction[mdl] yields the forward tranformation function for the
 given Schira model mdl. This is equivalent to mdl[Function].";
SchiraFunction::badarg = "Bad argument given to Schira function: `1`";
SchiraInverse::usage = "SchiraFunction[mdl] yields the inverse tranformation function for the
 given Schira model mdl. This is equivalent to mdl[Inverse].";
SchiraInverse::badarg = "Bad argument given to Schira inverse function: `1`";

(**************************************************************************************************)
(**************************************************************************************************)
Begin["`Private`"];

(* Visual Angles and Visual Areas *****************************************************************)
PolarAngle = PolarAngle;
Eccentricity = Eccentricity;
VisualArea = VisualArea;
Protect[PolarAngle, Eccentricity, VisualArea];

VisualAngleToComplex[th_, r_] := r * Exp[I * (Pi/2 - th * Pi / 180)];
ComplexToVisualAngle[z_] := {180/Pi * (Pi/2 - Arg[z]), Abs[z]};
ComplexToPolarAngle[z_] := 180/Pi * (Pi/2 - Arg[z]);
ComplexToEccentricity[z_] := Abs[z];
SetAttribute[VisualAngleToComplex, Listable];
SetAttribute[ComplexToVisualAngle, Listable];
SetAttribute[ComplexToPolarAngle, Listable];
SetAttribute[ComplexToEccentricity, Listable];
Protect[VisualAngleToComplex, ComplexToVisualAngle, ComplexToPolarAngle, ComplexToEccentricity];

$VisualAreasData = {
  1  -> {"Name" -> "V1 Dorsal",  "Areas" -> {"V1"}, "Simple" -> 1, "Stream" -> "Dorsal"},
  -1 -> {"Name" -> "V1 Ventral", "Areas" -> {"V1"}, "Simple" -> 1, "Stream" -> "Ventral"},
  2  -> {"Name" -> "V2 Dorsal",  "Areas" -> {"V2"}, "Simple" -> 2, "Stream" -> "Dorsal"},
  -2 -> {"Name" -> "V2 Ventral", "Areas" -> {"V2"}, "Simple" -> 2, "Stream" -> "Ventral"},
  3  -> {"Name" -> "V3 Dorsal",  "Areas" -> {"V3"}, "Simple" -> 3, "Stream" -> "Dorsal"},
  -3 -> {"Name" -> "V3 Ventral", "Areas" -> {"V3"}, "Simple" -> 3, "Stream" -> "Ventral"},
  4  -> {"Name" -> "V3A", "Areas" -> {"V3A"}, "Simple" -> 0, "Stream" -> "Dorsal"},
  -4 -> {"Name" -> "HV4", "Areas" -> {"HV4"}, "Simple" -> 0, "Stream" -> "Ventral"},
  5  -> {"Name" -> "Unknown Dorsal Area",  "Areas" -> {}, "Simple" -> 0, "Stream" -> "Dorsal"},
  -5 -> {"Name" -> "Unknown Ventral Area", "Areas" -> {}, "Simple" -> 0, "Stream" -> "Ventral"},
  0  -> {"Name" -> "V1 Horizontal Meridian", "Areas" -> {"V1"}, "Simple" -> 1, "Stream" -> None},
  (3/2)  -> {"Name" -> "V1/V2 Lower Vertical Meridian", "Areas" -> {"V1", "V2"}, "Simple" -> 1, "Stream" -> "Dorsal"},
  (-3/2) -> {"Name" -> "V1/V2 Upper Vertical Meridian", "Areas" -> {"V1", "V2"}, "Simple" -> 1, "Stream" -> "Ventral"},
  (5/2)  -> {"Name" -> "V2/V3 Lower Vertical Meridian", "Areas" -> {"V2", "V3"}, "Simple" -> 2, "Stream" -> "Dorsal"},
  (-5/2) -> {"Name" -> "V2/V3 Upper Vertical Meridian", "Areas" -> {"V2", "V3"}, "Simple" -> 2, "Stream" -> "Ventral"},
  (7/2)  -> {"Name" -> "V3/V3A Lower Vertical Meridian", "Areas" -> {"V3", "V3A"}, "Simple" -> 3, "Stream" -> "Dorsal"},
  (-7/2) -> {"Name" -> "V3/HV4 Upper Vertical Meridian", "Areas" -> {"V3", "HV4"}, "Simple" -> 3, "Stream" -> "Ventral"},
  (9/2)  -> {"Name" -> "V3A Upper Vertical Meridian", "Areas" -> {"V3A"}, "Simple" -> 0, "Stream" -> "Dorsal"},
  (-9/2) -> {"Name" -> "HV4 Lower Vertical Meridian", "Areas" -> {"HV4"}, "Simple" -> 0, "Stream" -> "Ventral"}};
$VisualAreasDispatch = Dispatch[
  Append[
    $VisualAreasData,
    a_ :> (Message[VisualAreaData::badarea, a]; Indeterminate)]];

VisualAreaData[] := $VisualAreasData;
VisualAreaData[id_] := Replace[id, $VisualAreasDispatch];
VisualAreaName[id_] := Check[Replace["Name", Replace[id, $VisualAreasDispatch]], Indeterminate];
VisualAreaSimplify[id_] := Check[Replace["Simple", Replace[id, $VisualAreasDispatch]], Indeterminate];
SetAttribute[VisualAreaData, Listable];
SetAttribute[VisualAreaName, Listable];
SetAttribute[VisualAreaSimplify, Listable];
Protect[$VisualAreasData, $VisualAreasDispatch, VisualAreaData, VisualAreaName, VisualAreaSimplify];

(* Retinotopic Templates **************************************************************************)
FSAverageSymRetinotopy := With[
  {field = Check[
     MapThread[
       Function[{pa, e, a},
         Which[
           Abs[a] < 1 || Abs[a] > 3, {None, None, None},
           NumberQ[pa] && pa > 90.0, {pa, e, -Abs[a]},
           NumberQ[pa] && pa <= 90.0, {pa, e, Abs[a]},
           True, {None, None, None}]],
       Map[
         Function[{fl}, First[Surfaces[Import[fl, "MGH"]]]],
         {"https://cfn.upenn.edu/aguirreg/public/ES_template/mgh_files/angle-template.sym.mgh",
          "https://cfn.upenn.edu/aguirreg/public/ES_template/mgh_files/eccen-template.sym.mgh",
          "https://cfn.upenn.edu/aguirreg/public/ES_template/mgh_files/areas-template.sym.mgh"}]],
     $Failed]},
  If[filed === $Failed,
    $Failed,
    (Unprotect[FSAverageSymRetinotopy];
     Set[FSAverageSymRetinotopy, Rule[field, FSAverageSymSphereSurface]];
     Protect[FSAverageSymRetinotopy];
     FSAverageSymRetinotopy)]];
FSAverageSymPolarAngle := With[
  {retino = Check[FSAverageSymRetinotopy, $Failed]},
  If[retino === $Failed,
    $Failed,
    (Unprotect[FSAverageSymPolarAngle];
     Set[FSAverageSymPolarAngle, Rule[Field[retino][[All,1]], FSAverageSymSphereSurface]];
     Protect[FSAverageSymPolarAngle];
     FSAverageSymPolarAngle)]];
FSAverageSymEccentricity := With[
  {retino = Check[FSAverageSymRetinotopy, $Failed]},
  If[retino === $Failed,
    $Failed,
    (Unprotect[FSAverageSymEccentricity];
     Set[FSAverageSymEccentricity, Rule[Field[retino][[All,2]], FSAverageSymSphereSurface]];
     Protect[FSAverageSymEccentricity];
     FSAverageSymEccentricity)]];
FSAverageSymVisualArea := With[
  {retino = Check[FSAverageSymRetinotopy, $Failed]},
  If[retino === $Failed,
    $Failed,
    (Unprotect[FSAverageSymVisualArea];
     Set[FSAverageSymVisualArea, Rule[Field[retino][[All,3]], FSAverageSymSphereSurface]];
     Protect[FSAverageSymVisualArea];
     FSAverageSymVisualArea)]];

SubjectRetinotopy[sub_, hemi:(LH|RH)] := With[
  {retino = Check[
     SurfaceResample[
       FSAverageSymRetinotopy,
       SubjectSymSurface[sub, hemi]],
     $Failed]},
  If[retino === $Failed,
    $Failed,
    (Unprotect[SubjectRetinotopy];
     Set[SubjectRetinotopy[sub, hemi], retino];
     Protect[SubjectRetinotopy];
     retino)]];
SubjectPolarAngle[sub_, hemi:(LH|RH)] := With[
  {retino = Check[
     SurfaceResample[
       FSAverageSymPolarAngle,
       SubjectSymSurface[sub, hemi]],
     $Failed]},
  If[retino === $Failed,
    $Failed,
    (Unprotect[SubjectPolarAngle];
     Set[SubjectPolarAngle[sub, hemi], retino];
     Protect[SubjectPolarAngle];
     retino)]];
SubjectEccentricity[sub_, hemi:(LH|RH)] := With[
  {retino = Check[
    SurfaceResample[
      FSAverageSymEccentricity,
      SubjectSymSurface[sub, hemi]],
    $Failed]},
  If[retino === $Failed,
    $Failed,
    (Unprotect[SubjectEccentricity];
     Set[SubjectEccentricity[sub, hemi], retino];
     Protect[SubjectEccentricity];
     retino)]];
SubjectVisualArea[sub_, hemi:(LH|RH)] := With[
  {retino = Check[
     SurfaceResample[
       FSAverageSymVisualArea,
       SubjectSymSurface[sub, hemi]],
     $Failed]},
  If[retino === $Failed,
    $Failed,
    (Unprotect[SubjectVisualArea];
     Set[SubjectVisualArea[sub, hemi], retino];
     Protect[SubjectVisualArea];
     retino)]];

Protect[
  FSAverageSymRetinotopy, FSAverageSymPolarAngle, FSAverageSymEccentricity, FSAverageSymVisualArea,
  SubjectRetinotopy, SubjectPolarAngle, SubjectEccentricity, SubjectVisualArea];


(* V1 Masks and Hulls *****************************************************************************)
FSAverageV1Mask[hemi:(LH|RH)] := With[
  {res = Check[
     Import[
       $FSAverage <> "/label/"<> If[hemi === LH, "lh", "rh"] <> ".v1.predict.label",
       "FreeSurferLabel",
       Max -> Length[Field[FSAverageCurvature[hemi]]]],
     $Failed]},
  If[res === $Failed,
    $Failed,
    (Unprotect[FSAverageV1Mask];
     Set[FSAverageV1Mask[hemi], res];
     Protect[FSAverageV1Mask];
     res)]];
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
    (Unprotect[FSAverageSymV1Mask];
     Set[FSAverageSymV1Mask, res];
     Protect[FSAverageSymV1Mask];
     res)]];

FSAverageV1[hemi:(LH|RH)] := With[
  {idcs = Check[
     Flatten[Position[Normal[Field[FSAverageV1Mask[hemi]]], 1]],
     $Failed]},
  If[idcs === $Failed,
    $Failed,
    (Unprotect[FSAverageV1];
     Set[FSAverageV1[hemi], idcs];
     Protect[FSAverageV1];
     idcs)]];
FSAverageSymV1 := With[
  {idcs = Check[
     Flatten[Position[Normal[Field[FSAverageSymV1Mask]], 1]],
     $Failed]},
  If[idcs === $Failed,
    $Failed,
    (Unprotect[FSAverageSymV1];
     Set[FSAverageSymV1, idcs];
     Protect[FSAverageSymV1];
     idcs)]];

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
  If[res === $Failed,
    $Failed,
    (Unprotect[FSAverageV1Hull];
     Set[FSAverageV1Hull[hemi], res];
     Protect[FSAverageV1Hull];
     res)]];
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
  If[res === $Failed,
    $Failed,
    (Unprotect[FSAverageSymV1Hull];
     Set[FSAverageSymV1Hull, res];
     Protect[FSAverageSymV1Hull];
     res)]];

Protect[
  FSAverageV1, FSAverageV1Mask, FSAverageV1Hull, 
  FSAverageSymV1, FSAverageSymV1Mask, FSAverageSymV1Hull];


(* Default Schira Parameters **********************************************************************)
$DefaultSchiraA = 2.0;
$DefaultSchiraB = 60.0;
$DefaultSchira\[CapitalLambda] = 1.0;
$DefaultSchira\[CapitalPsi] = 0.30;
$DefaultSchiraV1Size = 0.99;
$DefaultSchiraV2Size = 0.79;
$DefaultSchiraV3Size = 0.38;
$DefaultSchiraHV4Size = 0.3;
$DefaultSchiraV3ASize = 0.3;
$DefaultSchiraFC = {-0.14, -0.1};
$DefaultSchiraScale = {0.35, -0.35};
$DefaultSchiraShear = {{1, 0}, {0, 1}};
$SchiraParameters = List[
   A :> $DefaultSchiraA,
   B :> $DefaultSchiraB,
   \[CapitalLambda] :> $DefaultSchira\[CapitalLambda],
   \[CapitalPsi] :> $DefaultSchira\[CapitalPsi],
   V1Size :> $DefaultSchiraV1Size,
   V2Size :> $DefaultSchiraV2Size,
   V3Size :> $DefaultSchiraV3Size,
   HV4Size :> $DefaultSchiraHV4Size,
   V3ASize :> $DefaultSchiraV3ASize,
   FC :> $DefaultSchiraFC,
   Scale :> $DefaultSchiraScale,
   Shear :> $DefaultSchiraShear];
Protect[
    $DefaultSchiraA, $DefaultSchiraB, $DefaultSchira\[CapitalLambda], $DefaultSchira\[CapitalPsi],
    $DefaultSchiraV1Size, $DefaultSchiraV2Size, $DefaultSchiraV3Size,
    $DefaultSchiraHV4Size, $DefaultSchiraV3ASize,
    $DefaultSchiraFC, $DefaultSchiraScale, $DefaultSchiraShear,
    $SchiraParameters];

(* Schira Parameters Keys *************************************************************************)
A = A; B = B; \[CapitalLambda] = \[CapitalLambda];
V1Size = V1Size; V2Size = V2Size; V3Size = V3Size; HV4Size = HV4Size; V3ASize = V3ASize; 
\[CapitalPsi] = \[CapitalPsi]; FC = FC; Shear = Shear;
Protect[A, B, \[CapitalLambda], 
       V1Size, V2Size, V3Size, V3ASize, HV4Size,
       \[CapitalPsi], FC, Shear];

(* The Schira Model Objects ***********************************************************************)

(* private variable used below *)
Unprotect[$SchiraParameterPositions];
ClearAll[$SchiraParameterPositions];
$SchiraParameterPositions = Dispatch[
  MapThread[
   Rule,
   {$SchiraParameters[[All, 1]], Range[Length[$SchiraParameters]]}]];
Protect[$SchiraParameterPositions];

(* These actually compile the low-level Schira calculation functions *)
Unprotect[CompileSchiraFunction, CompileSchiraInverse];
ClearAll[CompileSchiraFunction, CompileSchiraInverse];
CompileSchiraFunction[a_, b_, lambda_, psi_, shearMtx_, scale_, fc_, areas_] := Check[
  With[
    {v1b = N[areas[[1]]],
     v2b = N[areas[[2]]],
     v3b = N[areas[[3]]],
     hv4b = N[areas[[4]]],
     v3ab = N[areas[[5]]],
     dsech1 = 0.1821,
     dsech2 = 0.76,
     fcx0 = N[Log[(a + lambda) / (b + lambda)]],
     mtx = MapThread[
       Append,
       {Dot[
          RotationMatrix[N[psi]],
          N[shearMtx]],
        N[fc]}]},
    (* sanity checks should already be done: just compile the function *)
    Compile[
      {{z, _Complex}},
      (* Layerless Transorm:
         First, set zLayered to the layerless tranform for each of the 4 areas *)
      With[
        {zLayered = With[
          {zz = Arg[z]*2.0/Pi},
          With[
            {sgn = If[zz == 0, 1.0, Sign[zz]]},
            Abs[z] * Exp[I * {
              v1b*zz,
              sgn * (v1b + v2b*(1.0 - Abs[zz])),
              sgn * (v1b + v2b + v3b*Abs[zz]),
              v1b + v2b + v3b + hv4b * (1.0 - 0.5 * (zz + 1.0)),
              -(v1b + v2b + v3b + 0.5 * v3ab * (zz + 1.0))}]]]},
        (* Log-Polar Transform:
           Now, do the log-polar part of the transform *)
        With[
          {zLogPolar = Table[
             With[
               {argz = Arg[
                  If[Re[zz] >= 0, zz + lambda, zz + 2.0 * lambda * (1.0 - Abs[Arg[zz]]/Pi)]],
                absz = Abs[
                  If[Re[zz] >= 0, zz + lambda, zz + 2.0 * lambda * (1.0 - Abs[Arg[zz]]/Pi)]]},
               If[absz == 0,
                 Log[a/b] + 0.0 * I,
                 Log@Divide[
                   a + absz * Exp[I * argz * Sech[argz]^(dsech1 * Sech[dsech2*Log[absz/a]])],
                   b + absz * Exp[I * argz * Sech[argz]^(dsech1 * Sech[dsech2*Log[absz/b]])]]]],
             {zz, zLayered}]},
          (* We center the FC on zero to start, by subtracting fcx0, then we scale and shear, and,
             last, we push things back to the specified FC *)
          Flatten@Dot[
            {{1.0, I}},
            mtx,
            (* Note that we flip the z here so that the arrangement matchis the LH *)
            {Re[zLogPolar] - fcx0, -Im[zLogPolar], {1.0, 1.0, 1.0, 1.0, 1.0}}]]],
      RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
      Parallelization -> True,
      RuntimeAttributes -> {Listable}]],
  $Failed];

CompileSchiraInverse[a_, b_, lambda_, psi_, shearMtx_, scale_, fc_, areas_] := Check[
  With[
    {v1b = N[areas[[1]]],
     v2b = N[areas[[2]]],
     v3b = N[areas[[3]]],
     hv4b = N[areas[[4]]],
     v3ab = N[areas[[5]]],
     dsech1 = 0.1821,
     dsech2 = 0.76,
     tol = 0.000001,
     fcx0 = N[Log[(a + lambda) / (b + lambda)]],
     mtx = MapThread[
       Append,
       {Dot[
          RotationMatrix[N[psi]],
          N[shearMtx]],
        N[fc]}],
     sideFn = With[
       {u = If[ListQ[scale] && Times@@scale < 0,
          {Cos[psi + 0.5*Pi], Sin[psi + 0.5*Pi]},
          {Cos[psi - 0.5*Pi], Sin[psi - 0.5*Pi]}]},
       Function[Sign[Dot[{Re@#, Im@#} - fc, u]]]]},
    (* sanity checks should already be done: we compile the function next:
       this function is basically identical to the compiled function just above, but it does not
       perform the layered transform (ie, translating area V1-V4 into an imaginary number. This
       allows us to make a simple inverse function *)
    With[
      {forward = Compile[
         {{z, _Complex}},
         (* Log-Polar Transform: Do the log-polar part of the transform *)
         With[
           {ztr = With[
              {argz = Arg[
                 If[Re[z] >= 0, z + lambda, z + 2.0 * lambda * (1.0 - Abs[Arg[z]]/Pi)]],
               absz = Abs[
                 If[Re[z] >= 0, z + lambda, z + 2.0 * lambda * (1.0 - Abs[Arg[z]]/Pi)]]},
              If[absz == 0,
                Log[a/b] + 0.0 * I,
                Log@Divide[
                  a + absz * Exp[I * argz * Sech[argz]^(dsech1 * Sech[dsech2*Log[absz/a]])],
                  b + absz * Exp[I * argz * Sech[argz]^(dsech1 * Sech[dsech2*Log[absz/b]])]]]]},
           (* We center the FC on zero to start, by subtracting fcx0, then we scale and shear, 
              and, last, we push things back to the specified FC *)
           First@Dot[{{1.0, I}}, mtx, {Re[ztr] - fcx0, Im[ztr], 1.0}]],
         RuntimeOptions -> {"Speed", "EvaluateSymbolically" -> False},
         Parallelization -> True,
         RuntimeAttributes -> {Listable}]},
      (* Now, we create an inverse function for the forward function *)
      With[
        {inverse = Function[
             Conjugate@FindRoot[
               # == forward[w],
               {w, 0.0001 - 0.0001*I, 0.0001 + 0.0001*I}
              ][[1,2]]]},
        (* And wrap this inverse in a translator for the areas *)
        Function[{z},
          With[
            {w = Check[inverse[z], Undefined],
             side = sideFn[z]},
            If[!NumberQ[w],
              {0.0, 0},
              With[
                {argw = Arg[w],
                 absw = Abs[w]},
                Which[
                  Round[Abs@argw, tol] == 0, {absw + 0.0*I, 0},
                  Round[v1b - Abs@argw, tol] >= 0, {
                    absw*Exp[I * Pi/2 * argw / v1b], 
                    If[side == 0, 1, side] * If[Round[Abs@argw - v1b, tol] == 0, 3/2, 1]},
                  Round[v2b - (Abs@argw - v1b), tol] >= 0, {
                    absw*Exp[I * Pi/2 * Sign[argw] * (1.0 - (Abs[argw] - v1b)/v2b)], 
                    side * If[Round[Abs@argw - v1b - v2b, tol] == 0, 5/2, 2]},
                  Round[v3b - (Abs@argw - v1b - v2b), tol] >= 0, {
                    absw*Exp[I * Pi/2 * Sign[argw] * (Abs[argw] - v1b - v2b)/v3b],
                    side * If[Round[Abs@argw - v1b - v2b - v3b, tol] == 0, 7/2, 3]},
                  side < 0 && Round[hv4b - (Abs@argw - v1b - v2b - v3b), tol] >= 0, {
                    absw*Exp[I * Pi * (0.5 - (Abs@argw - v1b - v2b - v3b)/hv4b)],
                    If[Round[Abs@argw - v1b - v2b - v3b - hv4b, tol] == 0, -9/2, -4]},
                  side > 0 && Round[v3ab - (Abs@argw - v1b - v2b - v3b), tol] >= 0, {
                    absw*Exp[I * Pi * ((Abs@argw - v1b - v2b - v3b)/v3ab - 0.5)],
                    If[Round[Abs@argw - v1b - v2b - v3b - v3ab, tol] == 0, 9/2, 4]},
                  side < 0, {
                    abs*Exp[I * Pi * ((Abs@argw - v1b - v2b - v3b - hv4b)/(Pi - v1b - v2b - v3b - hv4b) - 0.5)],
                    -5},
                  side > 0, {
                    abs*Exp[I * Pi * (0.5 - (Abs@argw - v1b - v2b - v3b - v3ab)/(Pi - v1b - v2b - v3b - v3ab))],
                    5},
                  True, {-absw, Infinity}]]]],
          {Listable}]]]],
  $Failed];
Protect[CompileSchiraFunction, CompileSchiraInverse];


(* This is used to force preparation of the function and inverse below *)
Unprotect[SchiraModelObjectPrep];
ClearAll[SchiraModelObjectPrep];
SchiraModelObjectPrep[params_List] := With[
  {ff = Unique["fun"],
   if = Unique["inv"],
   a = params[[A /. $SchiraParameterPositions, 2]],
   b = params[[B /. $SchiraParameterPositions, 2]],
   lambda = params[[\[CapitalLambda] /. $SchiraParameterPositions, 2]],
   psi = params[[\[CapitalPsi] /. $SchiraParameterPositions, 2]] /. None -> 0,
   shearMtx = params[[Shear /. $SchiraParameterPositions, 2]] /. None -> {{1,0},{0,1}},
   scale = Replace[params[[Scale /. $SchiraParameterPositions, 2]], x_?NumericQ :> {x,x}],
   fc = Replace[params[[FC /. $SchiraParameterPositions, 2]], z_Complex :> {Re[x], Im[z]}],
   areas = params[[{V1Size, V2Size, V3Size, HV4Size, V3ASize} /. $SchiraParameterPositions, 2]]},
  (* sanity checking *)
  Which[
    !NumericQ[a] || a <= 0, Message[SchiraModelObject::badarg, "A must be numeric and > 0"],
    !NumericQ[b] || b <= 0, Message[SchiraModelObject::badarg, "B must be numeric and > 0"],
    !NumericQ[lambda] || lambda < 0, Message[
      SchiraModelObject::badarg,
      "\[CapitalLambda] must be numeric and >= 0"],
    !NumericQ[psi], Message[
      SchiraModelObject::badarg,
      "\[CapitalPsi] must be numeric and >= 0"],
    Dimensions[shearMtx] != {2,2} || shearMtx[[1,1]] != 1 || shearMtx[[2,2]] != 1, Message[
      SchiraModelObject::badarg,
      "Shear must be a 2 x 2 matrix with ones on the diagonal or None"],
    !NumericQ[shearMtx[[1,2]]] || !NumericQ[shearMtx[[2,1]]], Message[
      SchiraModelObject::badarg,
      "Shear matrix must have numeric off-diagonal elements"],
    !ListQ[scale] || Length[scale] != 2 || !NumericQ[scale[[1]]] || !NumericQ[scale[[2]]], Message[
      SchiraModelObject::badarg,
      "Scale must be a single numeric quantity or a pair of {x-scale, y-scale} numeric quantities"],
    !ListQ[fc] || Length[scale] != 2 || !NumericQ[fc[[1]]] || !NumericQ[fc[[2]]], Message[
      SchiraModelObject::badarg,
      "Foveal convluence (FC) must be a either a coordinate {x, y} with numeric quantities or a"
       <> " single complex number with numeric real and imaginary parts"],
    Not[And@@Map[NumericQ, areas]] || Not[And@@Map[(#>0)&, areas]], Message[
      SchiraModelObject::badarg,
      "V1Size, V2Size, V3Size, HV4Size, and V3ASize must all be numeric quantities > 0"],
    Total[Most[areas]] > Pi || Total[areas[[{1,2,3,5}]]] > Pi, Message[
      SchiraModelObject::badarg,
      "V1Size, V2Size, V3Size, and neither HV4Size nor V3ASize may sum to be >= Pi"]];
  (* Note that these are temporary variables *)
  SetAttributes[Evaluate[ff], Temporary];
  SetAttributes[Evaluate[if], Temporary];
  (* set these to auto-memoize themselves if requested *)
  ff := With[
    {f = Check[
       CompileSchiraFunction[a, b, lambda, psi, shearMtx, scale, fc, areas],
       $Failed]},
    If[f === $Failed, $Failed, (ff = f)]];
  if := With[
    {f = Check[
       CompileSchiraInverse[a, b, lambda, psi, shearMtx, scale, fc, areas],
       $Failed]},
    If[f === $Failed, $Failed, (if = f)]];
  (* And make a dispatch for this model *)
  SchiraModelObject[
    Dispatch[
      Join[
        params,
        {Function :> ff,
         Inverse :> if,
         All :> params,
         x_ :> Message[SchiraModelObject::badarg, x]}]]]];
Protect[SchiraModelObjectPrep];

SchiraModel[
   opts : Evaluate[
     Repeated[
      (Rule | RuleDelayed)[
       Apply[Alternatives, First /@ $SchiraParameters],
       _]]]
 ] := SchiraModelObjectPrep[
  ReplaceAll[
    $SchiraParameters,
    Map[
      Function[{rule},
        Rule[
          (First[rule] :> _),
          Head[rule][First[rule], Last[rule]]]],
        {opts}]]];
SchiraModel[
   SchiraModelObject[disp_],
   opts : Evaluate[
     Repeated[
      (Rule | RuleDelayed)[
       Apply[Alternatives, First /@ $SchiraParameters],
       _]]]
 ] := SchiraModelObjectPrep[
  ReplaceAll[
    Select[
      disp[[1]],
      And[Head[#[[1]]] =!= Pattern, MemberQ[$SchiraParameters[[All,1]], #[[1]]]]&],
    Map[
      Function[{rule},
        Rule[
          (Rule|DelayedRule)[First[rule], _],
          Head[rule][First[rule], Last[rule]]]],
        {opts}]]];
SchiraModel[] := SchiraModelObjectPrep[$SchiraParameters];
 
SchiraModelObject[disp_][x_] := Replace[x, disp];
SchiraFunction[SchiraModelObject[disp_]] := Replace[Function, disp];
SchiraInverse[SchiraModelObject[disp_]] := Replace[Inverse, disp];

Protect[SchiraModel, SchiraModelObject, SchiraFunction, SchiraInverse];

End[];
EndPackage[];

