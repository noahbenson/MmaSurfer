(* FreeSurfer.m
 *
 * Basic utility functions for dealing with FreeSurfer data in Mathematica.
 *
 * Copyright (C) 2013-2014 by Noah C. Benson.
 * This file is part of the MmaSurfer library, which is provided under the terms of the Eclipse
 * Public License, version 1.0. See the accompanying LICENSE file for more information.
 *)

(**************************************************************************************************)
BeginPackage["FreeSurfer`", {"CorticalSurface`", "CorticalVolume`"}];

Unprotect[Frames, ImageBufferType, DegreesOfFreedom, VOXToRASMatrix, Volumes, Surfaces, 
          OptionalData, Header, HeaderComment, ImportedData, LookupTable];
ClearAll[ Frames, ImageBufferType, DegreesOfFreedom, VOXToRASMatrix, Volumes, Surfaces, 
          OptionalData, Header, HeaderComment, ImportedData, LookupTable];
Frames::usage = "Frames[mgh] yields the number of frames in an mgh object.";
ImageBufferType::usage = "ImageBufferType[mgh] yields the image buffer type for the given mgh object.";
DegreesOfFreedom::usage = "DegreesOfFreedom[mgh] yields the number of degrees of freedom of the given
 mgh object.";
VOXToRASMatrix::usage = "VOXToRASMatrix[mgh] yields the vox2ras matrix of the given mgh object.";
Volumes::usage = "Volumes[mgh] yields the array of 3D images in the given mgh object.";
Surfaces::usage = "Surfaces[mgh] yields the surface data of the given mgh object.";
OptionalData::usage = "OptionalData[mgh] yields a list of the optional data encoded in the file
 from which the mgh object was read.  OptionalData[mgh, name] yields the value of the
 optional datum, if stored in the file, with the given name, or None if no such datum was
 stored in the file.";
Header::usage = "Header[obj] yields the header data from certain FreeSurfer objects that have been
 imported.";
HeaderComment::usage = "HeaderComment[object] yields the header comment associated with reading in
 some objects.";
ImportedData::usage = "ImportedData[object] yields the raw data that was imported for object, if
 applicable.";
LookupTable::usage = "LookupTable[object] yields the color/label lookup table for an annotation
 object.";
Protect[  Frames, ImageBufferType, DegreesOfFreedom, VOXToRASMatrix, Volumes, Surfaces, OptionalData,
          Header, HeaderComment, ImportedData];

Unprotect[ImportMGH, ImportSurface, ImportCurv, ImportWeights, ImportAnnotation, ImportLabel,
          AnnotationQ, MGHObjectFromData, MGHQ];
ClearAll[ ImportMGH, ImportSurface, ImportCurv, ImportWeights, ImportAnnotation, ImportLabel, 
          AnnotationQ, MGHObjectFromData, MGHQ];
ImportMGH::usage = "ImportMGH[filename, options...] is identical to Import[filename, \"MGH\",
 options].";
ImportSurface::usage = "ImportSurface[filename, options...] is identical to Import[filename,
 \"FreeSurferSurface\", options].";
ImportCurc::usage = "ImportCurv[filename, options...] is identical to Import[filename,
 \"FreeSurferCurv\", options].";
ImportWeights::usage = "ImportWeights[filename, options...] is identical to Import[filename,
 \"FreeSurferWeights\", options].";
ImportAnnotation::usage = "ImportAnnotation[filename, options...] is identical to Import[filename,
 \"FreeSurferAnnotation\", options].";
ImportLabel::usage = "ImportLabel[filename, options...] is identical to Import[filename, 
\"FreeSurferLabel\", options].";
ImportMGH::badfmt = "Bad format when reading `1`: `2`";
ImportMGH::version = "ImportMGH expects a version 1 file; this file is version `1`";
ImportMGH::wrongtype = "MGH object appears to be a `1` type, but `2` type requested.";
ImportMGH::nofile = "File `1` given to ImportMGH could not be found";
ImportSurface::badfmt = "Bad format when importing surface `1`: `2`";
ImportCurv::badfmt = "Bad format when importing curv `1`: `2`";
ImportWeights::badfmt = "Bad format when importing weights `1`: `2`";
ImportAnnotation::badfmt = "Bad format when importing annotation `1`: `2`";
ImportAnnotation::warning = "Warning from ImportAnnotation: `1`";
MGHObjectFromData::usage = "MGHObjectFromData[data] yields an MGH object from the data (ie obtained
 via Import[... {\"MGH\", \"Data\"}]).";
MGHObjectFromData::badfmt = "Bad data format given to MGHObjectFromData: `1`";
MGHQ::usage = "MGHQ[obj] yields true if and only if obj is an MGH object.";
AnnotationQ::usage = "AnnotationQ[obj] yields true if and only if obj is an Annotation object.";

Unprotect[ExportMGH, ExportSurface, ExportWeights, ExportCurv, ExportAnnotation, ExportLabel];
ClearAll[ ExportMGH, ExportSurface, ExportWeights, ExportCurv, ExportAnnotation, ExportLabel];
ExportMGH::usage = "ExportMGH[filename, data, options...] is equivalent to Export[filename, data,
 \"MGH\", options].";
ExportSurface::usage = "ExportSurface[filename, data, options...] is equivalent to Export[filename,
 data, \"FreeSurferSurface\", options].";
ExportWeights::usage = "ExportWeights[filename, data, options...] is equivalent to Export[filename,
 data, \"FreeSurferWeights\", options].";
ExportCurv::usage = "ExportCurv[filename, data, options...] is equivalent to Export[filename, data,
 \"FreeSurferCurv\", options].";
ExportAnnotation::usage = "ExportAnnotation[filename, data, options...] is equivalent to
 Export[filename, data, \"FreeSurferAnnotation\", options].";
ExportLabel::usage = "ExportLabel[filename, data, options...] is identical to Export[filename,
 data, \"FreeSurferLabel\", options].";
ExportLabel::badfmt = "ExportLabel given data not in the form of {vertexID -> {p, {R, A, S}} ...}";

Unprotect[$MGHHeaderSize, $MGHOptionalData];
ClearAll[ $MGHHeaderSize, $MGHOptionalData];
$MGHHeaderSize::usage = "$MGHHeaderSize is the number of bytes in the header of an MGH file.";
$MGHOptionalData::usage = "$MGHOptionalData is a list of the optional pieces of data that may be
 stored at the end of an MGH file.";

Unprotect[$FreeSurferHomes, $SubjectsDirs, $Subjects, $FreeSurferColorLUT];
ClearAll[ $FreeSurferHomes, $SubjectsDirs, $Subjects, $FreeSurferColorLUT];
$SubjectsDirs::usage = "$SubjectsDirs is a list of potential FreeSurfer subjects directories; it is
 obtained at runtime from the $SUBJECTS_DIR environment variable as well as by scanning
 common FreeSurfer paths.";
$Subjects::usage = "$Subjects is a list of potential FreeSurfer subjects (named by their
 directories). If you wish to add an arbitrary subject to FreeSurfer, use AddSubject.";
$FreeSurferHomes::usage = "$FreeSurferHomes is a list of potential home directories for FreeSurfer.
 Homes may be added with AddFreeSurferHome.";
$FreeSurferColorLUT::usage = "$FreeSurferColorLUT is a dispatch table that will replace either an
 integer or string volume label and yield the remaining data (i.e., a list of the string or
 integer label, whichever was not given, and a color.";

Unprotect[RH, LH, RHX];
ClearAll[RH, LH, RHX];
RH::usage = "RH is a keyword that represents the right hemisphere.";
LH::usage = "LH is a keyword that represents the left hemisphere.";
RHX::usage = "RHX is a keyword that represents the inverted right hemisphere.";
Protect[RH, LH, RHX];

Unprotect[AddFreeSurferHome, RemoveFreeSurferHome, AddSubjectsDir, RemoveSubjectsDir, AddSubject,
          RemoveSubject];
ClearAll[ AddFreeSurferHome, RemoveFreeSurferHome, AddSubjectsDir, RemoveSubjectsDir, AddSubject,
          RemoveSubject];
AddFreeSurferHome::usage = "AddFreeSurferHome[dir] adds the directory dir to the FreeSurfer home
 structures.";
RemoveFreeSurferHome::usage = "RemoveFreeSurferHome[dir] removes the directories matching the
 pattern dir from the FreeSurfer subject directories structures.";
AddSubjectsDir::usage = "AddFreeSurferHome[dir] adds the directory dir to the FreeSurfer home
 structures.";
RemoveSubjectsDir::usage = "RemoveFreeSurferHome[dir] removes the directories matching the pattern
 dir from the FreeSurfer subject directories structures.";
AddSubject::usage = "AddSubject[dir] adds the directory dir to the FreeSurfer subject structures.";
RemoveSubject::usage = "RemoveSurbject[dir] removes the directories matching the pattern dir from
 the FreeSurfer subject structures.";

Unprotect[Subject, SubjectQ];
ClearAll[Subject, SubjectQ];
Subject::usage = "Subject[directory] represents a FreeSurfer subject whose data is stored in the
 directory named by directory.";
Subject::notfound = "Subject directory `1` does not seem to exist";
Subject::baddata = "Subject directory `1` does not seem to contains FreeSurfer data: `2`";
SubjectQ::usage = "SubjectQ[directory] yields true if and only if directory is a string referring
 to a directory that contains a FreeSurfer subject.";

(* Volume Functions *)
Unprotect[SubjectSegments, SubjectSegments, SubjectBrain, SubjectWhiteMatter, SubjectFilledBrain,
          SubjectHemisphere];
ClearAll[ SubjectSegments, SubjectSegments, SubjectBrain, SubjectWhiteMatter, SubjectFilledBrain,
          SubjectHemisphere];
SubjectSegments::usage = "SubjectSegments[sub] yields an MGH object for subject sub whose values
 correspond to segments of the brain volume.";
SubjectSegment::usage = "SubjectSegment[sub, label] yeids a list of indices at which the subject's
 anatomical volume is labeled with label, which may be an integer or a string, either of
 which must be found in the FreeSurferColorLUT.";
SubjectSegment::badlbl = "Label given to SubjectSegment not found: `1`";
SubjectBrain::usage = "SubjectBrain[sub] yields the volume for the subject sub's normalized brain
 (after skull stripping).";
SubjectWhiteMatter::usage = "SubjectBrain[sub] yields the volume for the subject sub's white
 matter.";
SubjectFilledBrain::usage = "SubjectFilledBrain[sub] yields the volume for the subject sub's brain
 in which the right hemisphere white matter has values of Left and the left hemisphere has
 values of Right. Note that these can be exported and will translate to FreeSurfer's
 specified RH -> 127 and LH -> 255 values.";
SubjectHemisphere::usage = "SubjectHemisphere[sub, LH|RH] yields the volume for the subject sub's
 right or left hemisphere only (with values of 1 for the white matter and 0 elsewhere).";

(* Surface Functions *)
Unprotect[SubjectOriginalSurface, SubjectPialSurface, SubjectInflatedSurface, SubjectSphereSurface,
          SubjectRegisteredSurface, SubjectSymSurface];
ClearAll[ SubjectOriginalSurface, SubjectPialSurface, SubjectInflatedSurface, SubjectSphereSurface,
          SubjectRegisteredSurface, SubjectSymSurface];
SubjectOriginalSurface::usage = "SubjectOriginalSurface[sub, hemi] yields the original cortical
 surface tesselation for subject sub's specified hemishere.";
SubjectPialSurface::usage = "SubjectPialSurface[sub, hemi] yields the pial surface object for
 subject sub's specified hemishere.";
SubjectInflatedSurface::usage = "SubjectInflatedSurface[sub, hemi] yields the inflated surface
 object for subject sub's specified hemishere.";
SubjectSphereSurface::usage = "SubjectSphereSurface[sub, hemi] yields the spherical surface object
 for subject sub's specified hemishere; note that this is not the registration of the"
    " subject to the FSAverage hemisphere. For that, see SubjectRegisteredSurface..";
SubjectRegisteredSurface::usage = "SubjectRegisteredSurface[sub, hemi] yields the subject's
 cortical surface registered to the spherical fsaverage hemisphere for subject sub's
 specified hemishere.";
SubjectSymSurface::usage = "SubjectSymSurface[sub, hemi] yields the subject's
 cortical surface registered to the spherical fsaverage_sym hemisphere for subject sub's
 specified hemishere.";

(* Surface Field Functions *)
Unprotect[SubjectJacobian, SubjectCurvature, SubjectSulci, SubjectThickness, SubjectVertexArea,
          SubjectVertexAreaPial, SubjectVertexAreaWhite, SubjectVolume, SubjectParcellation, 
          SubjectParcellation2009, SubjectParcellation2005];
ClearAll[ SubjectJacobian, SubjectCurvature, SubjectSulci, SubjectThickness, SubjectVertexArea,
          SubjectVertexAreaPial, SubjectVertexAreaWhite, SubjectVolume, SubjectParcellation,
          SubjectParcellation2009, SubjectParcellation2005];
SubjectJacobian::usage = "SubjectJacobian[sub, hemi] yields the Jacobian field for the subejct
 sub's given hemisphere.";
SubjectCurvature::usage = "SubjectCurvature[sub, hemi] yields the curvature for subject sub's
 given hemisphere.";
SubjectSulci::usage = "SubjectSulci[sub, hemi] yields the sulcal depth for subject sub's given
 hemisphere.";
SubjectThickness::usage = "SubjectThickness[sub, hemi] yields the thickness of each vertex for
 subject sub's given hemisphere.";
SubjectVertexArea::usage = "SubjectVertexArea[sub, hemi] yields the area of each vertex for subject
 sub's given hemisphere (?h.mid.area).";
SubjectVertexAreaPial::usage = "SubjectVertexAreaPial[sub, hemi] yields the pial area of each
 vertex for subject sub's given hemisphere (?h.area.pial).";
SubjectVertexAreaWhite::usage = "SubjectVertexAreaWhite[sub, hemi] yields the white matter area of
 each vertex for subject sub's given hemisphere (?h.area.white).";
SubjectVertexVolume::usage = "SubjectVertexVolume[sub, hemi] yields the volume of each vertex for
 subject sub's given hemisphere";
SubjectParcellation::usage = "SubjectParcellation[sub, hemi] yields the 2009 cortical surface
 parcellation for subject sub's given hemisphere.";
SubjectParcellation2009::usage = "SubjectParcellation[sub, hemi] yields the 2009 cortical surface
 parcellation for subject sub's given hemisphere.";
SubjectParcellation2005::usage = "SubjectParcellation[sub, hemi] yields the 2005 cortical surface
 parcellation for subject sub's given hemisphere.";

Unprotect[FreeSurfer];
ClearAll[FreeSurfer];
FreeSurfer::nolabel = "The FreeSurferColorLUT.txt file could not be found. This file is normally in
 your FreeSurfer home directory. Try adding a FreeSurfer home directory (via
 AddFreeSurferHome[]) and re-evaluating.";
FreeSurfer::notfound = "Data not found: `1`";
Protect[FreeSurfer];

Unprotect[FSAverageOriginalSurface, FSAveragePialSurface, FSAverageInflatedSurface, 
          FSAverageSphereSurface, FSAverageCurvature, FSAverageSulci, FSAverageThickness,
          FSAverageRegisteredCurvature, FSAverageVertexArea, FSAverageVertexAreaPial,
          FSAverageVertexAreaWhite, FSAverageVolume, FSAverageParcellation, 
          FSAverageParcellation2009, FSAverageParcellation2005, $FSAverage];
ClearAll[ FSAverageOriginalSurface, FSAveragePialSurface, FSAverageInflatedSurface, 
          FSAverageSphereSurface, FSAverageCurvature, FSAverageSulci, FSAverageThickness,
          FSAverageRegisteredCurvature, FSAverageVertexArea, FSAverageVertexAreaPial, 
          FSAverageVertexAreaWhite, FSAverageVolume, FSAverageParcellation, 
          FSAverageParcellation2009, FSAverageParcellation2005, $FSAverage];
FSAverageOriginalSurface::usage = "FSAverageOriginalSurface[hemi] yields the fsaverage subject's
 original surface.";
FSAveragePialSurface::usage = "FSAveragePialSurface[hemi] yields the fsaverage subject's pial
 surface.";
FSAverageInflatedSurface::usage = "FSAverageInflatedSurface[hemi] yields the fsaverage subject's
 inflated surface.";
FSAverageSphereSurface::usage = "FSAverageSphereSurface[hemi] yields the fsaverage subject's
 spherical surface.";
FSAverageCurvature::usage = "FSAverageCurvature[hemi] yields the fsaverage subject's curvature.";
FSAverageSulci::usage = "FSAverageSulci[hemi] yields the fsaverage subject's sulcal depth.";
FSAverageThickness::usage = "FSAverageThickness[hemi] yields the fsaverage subject's cortical
 thickness.";
FSAverageRegisteredCurvature::usage = "FSAverageRegisteredCurvature[hemi] yields the fsaverage
 subject's average curvature.";
FSAverageVertexArea::usage = "FSAverageVertexArea[hemi] yields the fsaverage subject's vertex areas
 (mid.area).";
FSAverageVertexAreaPial::usage = "FSAverageVertexAreaPial[hemi] yields the fsaverage subject's
 vertex areas (area.pial)";
FSAverageVertexAreaWhite::usage = "FSAverageVertexAreaWhite[hemi] yields the fsaverage subject's
 white matter vertex areas (area.white).";
FSAverageVolume::usage = "FSAverageVolume[hemi] yields the fsaverage subject's vertex volumes.";
FSAverageParcellation::usage = "FSAverageParcellation[hemi] yields the fsaverage subject's region
 parcellations (2009).";
FSAverageParcellation2009::usage = "FSAverageParcellation2009[hemi] yields the fsaverage subject's
 region parcellations (2009).";
FSAverageParcellation2005::usage = "FSAverageParcellation2005[hemi] yields the fsaverage subject's
 region parcellations (2005).";
$FSAverage::usage = "$FSAverage is the subject directory for the fsaverage FreeSurfer subject. If
 you add your freesurfer directory (or it is auto-detected), this will be automatically
 discovered (see AddFreeSurferHome, AddFreeSurferSubjectsDir, and AddFreeSurferSubject).";

Unprotect[FSAverageSymOriginalSurface, FSAverageSymPialSurface, FSAverageSymInflatedSurface,
          FSAverageSymSphereSurface, FSAverageSymCurvature, FSAverageSymSulci, 
          FSAverageSymThickness, FSAverageSymRegisteredCurvature, FSAverageSymVertexArea, 
          FSAverageSymVertexAreaPial, FSAverageSymVertexAreaWhite, FSAverageSymVolume, 
          FSAverageSymParcellation, FSAverageSymParcellation2009, FSAverageSymParcellation2005,
          $FSAverageSym];
ClearAll[ FSAverageSymOriginalSurface, FSAverageSymPialSurface, FSAverageSymInflatedSurface, 
          FSAverageSymSphereSurface, FSAverageSymCurvature, FSAverageSymSulci, 
          FSAverageSymThickness, FSAverageSymRegisteredCurvature, FSAverageSymVertexArea, 
          FSAverageSymVertexAreaPial, FSAverageSymVertexAreaWhite, FSAverageSymVolume, 
          FSAverageSymParcellation, FSAverageSymParcellation2009, FSAverageSymParcellation2005, 
          $FSAverageSym];
FSAverageSymOriginalSurface::usage = "FSAverageSymOriginalSurface yields the fsaverage_sym
 subject's original surface.";
FSAverageSymPialSurface::usage = "FSAverageSymPialSurface yields the fsaverage_sym subject's pial
 surface.";
FSAverageSymInflatedSurface::usage = "FSAverageSymInflatedSurface yields the fsaverage_sym
 subject's inflated surface.";
FSAverageSymSphereSurface::usage = "FSAverageSymSphereSurface yields the fsaverage_sym subject's
 spherical surface.";
FSAverageSymCurvature::usage = "FSAverageSymCurvature yields the fsaverage_sym subject's
 curvature.";
FSAverageSymSulci::usage = "FSAverageSymSulci yields the fsaverage_sym subject's sulcal depth.";
FSAverageSymThickness::usage = "FSAverageSymThickness yields the fsaverage_sym subject's cortical
 thickness.";
FSAverageSymRegisteredCurvature::usage = "FSAverageSymRegisteredCurvature yields the fsaverage_sym
 subject's average curvature.";
FSAverageSymVertexArea::usage = "FSAverageSymVertexArea yields the fsaverage_sym subject's vertex
 areas (mid.area).";
FSAverageSymVertexAreaPial::usage = "FSAverageSymVertexAreaPial yields the fsaverage_sym subject's
 vertex areas (area.pial)";
FSAverageSymVertexAreaWhite::usage = "FSAverageSymVertexAreaWhite yields the fsaverage_sym
 subject's white matter vertex areas (area.white).";
FSAverageSymVolume::usage = "FSAverageSymVolume yields the fsaverage_sym subject's vertex volumes.";
FSAverageSymParcellation::usage = "FSAverageSymParcellation yields the fsaverage_sym subject's
 region parcellations (2009).";
FSAverageSymParcellation2009::usage = "FSAverageSymParcellation2009 yields the fsaverage_sym
 subject's region parcellations (2009).";
FSAverageSymParcellation2005::usage = "FSAverageSymParcellation2005 yields the fsaverage_sym
 subject's region parcellations (2005).";
$FSAverageSym::usage = "$FSAverageSym is the subject directory for the fsaverage_sym FreeSurfer
 subject. If you add your freesurfer directory (or it is auto-detected), this will be
 automatically discovered (see AddFreeSurferHome, AddFreeSurferSubjectsDir, and
 AddFreeSurferSubject).";


(**************************************************************************************************)
Begin["`Private`"];

$MGHHeaderSize = 284;
Protect[$MGHHeaderSize];

$MGHOptionalData = {{"TR", "Real32"},
                    {"FlipAngle", "Real32"},
                    {"TE", "Real32"},
                    {"TI", "Real32"},
                    {"FoV", "Real32"}};
Protect[$MGHOptionalData];

(* Translate between types and  *)
$MGHTypesToMMA = {
  0 -> "UnsignedInteger8",
  1 -> "Integer32",
  3 -> "Real32",
  4 -> "Integer16",
  t_ :> (
    Message[ImportMGH::badfmt, "header (type)", "Type ("<>ToString[t]<>") not recognized"];
    Throw[$Failed])};
$MMATypesToMGH = Append[
  Map[Rule[#[[2]],#[[1]]]&, Most[$MGHTypesToMMA]],
  t_ :> (
    Message[ExportMGH::badfmt, "header (type)", "Type ("<>ToString[t]<>") not recognized"];
    Throw[$Failed])];

MGHQ[obj_] := False;
Protect[MGHQ];

AnnotationQ[obj_] := False;
Protect[AnnotationQ];

Unprotect[$BytesPerType];
ClearAll[$BytesPerType];
$BytesPerType = {"Real32" -> 4, "UnsignedInteger8" -> 1, "Integer32" -> 4, "Integer16" -> 2};
Protect[$BytesPerType];

(* These are private; for use in registering FreeSurfer's import file formats. *)
Unprotect[ImportMGHHeader, ImportMGHFrames, ImportMGHFooter, ImportMGHData, ImportMGHObject];
ClearAll[ImportMGHHeader, ImportMGHFrames, ImportMGHFooter, ImportMGHData, ImportMGHObject];
ImportMGHHeader[stream_InputStream, opts___] := "Header" -> Catch[
  Block[
    {$ByteOrdering = 1},
    With[
      {version = BinaryRead[stream, "Integer32"],
       sane = TrueQ["SanityChecks" /. Append[{opts}, "SanityChecks" -> True]]},
      If[version != 1, Message[ImportMGH::version, version]];
      With[
        {width = BinaryRead[stream, "Integer32"],
         height = BinaryRead[stream, "Integer32"],
         depth = BinaryRead[stream, "Integer32"],
         nframes = BinaryRead[stream, "Integer32"],
         type = Replace[
           BinaryRead[stream, "Integer32"],
           $MGHTypesToMMA],
         dof = BinaryRead[stream, "Integer32"],
         goodRASFlag = BinaryRead[stream, "Integer16"]},
        With[
          {spacings = If[goodRASFlag == 1, BinaryReadList[stream, "Real32", 3], {1.0, 1.0, 1.0}],
           Xras = If[goodRASFlag == 1, BinaryReadList[stream, "Real32", 3], {-1.0, 0.0, 0.0}],
           Yras = If[goodRASFlag == 1, BinaryReadList[stream, "Real32", 3], {0.0, 0.0, -1.0}],
           Zras = If[goodRASFlag == 1, BinaryReadList[stream, "Real32", 3], {0.0, 1.0, 0.0}],
           Cras = If[goodRASFlag == 1, BinaryReadList[stream, "Real32", 3], {0.0, 0.0, 0.0}]},
          {Version -> version,
           Dimensions -> {width, height, depth},
           Frames -> nframes,
           ImageBufferType -> type,
           DegreesOfFreedom -> dof,
           Spacings -> spacings,
           VOXToRASMatrix -> Transpose[
             {Append[Xras, 0.0],
              Append[Yras, 0.0],
              Append[Zras, 0.0],
              Append[Cras, 1.0]}]}]]]]];
ImportMGHFrames[stream_InputStream, opts___] := "Frames" -> Catch[
  Block[
    {$ByteOrdering = 1},
    With[
      {header = Replace[
         "Header",
         Append[
           {opts}, 
           "Header" :> Replace[
             ("Header" /. ImportMGHHeader[stream, opts]),
             $Failed :> Throw[$Failed]]]],
       lopts = {opts}},
      With[
        {dims = Dimensions /. header,
         type = ImageBufferType /. header,
         nframes = Frames /. header},
        With[
          {volsz = Apply[Times, dims[[1;;3]]]},
          If[TrueQ["SanityChecks" /. Append[lopts, "SanityChecks" -> True]],
            (* make sure the number of things we're about to read is reasonable *)
            If[volsz < 1 || volsz > 10^9,
              Message[
                ImportMGH::badfmt,
                "import stream",
                StringJoin[
                  "size of volume (", ToString[volsz], 
                  ") seems unreasonable; run with option \"SanityChecks\" -> False to ignore this error"]];
              Throw[$Failed]];
            If[nframes < 1 || nframes > 10^5,
              Message[
                ImportMGH::badfmt,
                "import stream",
                StringJoin[
                  "number of frames (", ToString[nframes], 
                  ") seems unreasonable; run with option \"SanityChecks\" -> False to ignore this error"]];
              Throw[$Failed]]];
            SetStreamPosition[stream, $MGHHeaderSize];
            Table[
              Transpose[
                Partition[
                  Partition[
                    BinaryReadList[stream, type, volsz],
                    dims[[1]]],
                  dims[[2]]],
                {3,2,1}],
              {nframes}]]]]]];
ImportMGHFooter[stream_InputStream, opts___] := "OptionalData" -> Catch[
  Block[
    {$ByteOrdering = 1,
     lopts = {opts}},
    With[
      {header = Replace[
         "Header",
         Append[
           {opts}, 
           "Header" :> Replace[
             ("Header" /. ImportMGHHeader[stream, opts]),
             $Failed :> Throw[$Failed]]]]},
      With[
        {dims = Dimensions /. header,
         type = ImageBufferType /. header,
         nframes = Replace[Frames, header]},
        With[
          {volsz = Apply[Times, dims[[1;;3]]]},
          If[TrueQ["SanityChecks" /. Append[lopts, "SanityChecks" -> True]],
            (* make sure the number of things we're about to read is reasonable *)
            If[volsz < 1 || volsz > 10^9,
              Message[
                ImportMGH::badfmt,
                "import stream",
                StringJoin[
                  "size of volume (", ToString[volsz], 
                  ") seems unreasonable; run with option \"SanityChecks\" -> False to ignore this error"]];
              Throw[$Failed]];
            If[nframes < 1 || nframes > 10^5,
              Message[
                ImportMGH::badfmt,
                "import stream",
                StringJoin[
                  "number of frames (",
                  ToString[nframes],
                  ") seems unreasonable; run with option \"SanityChecks\" -> False to ignore this error"]];
              Throw[$Failed]]];
            SetStreamPosition[stream, $MGHHeaderSize + (volsz * nframes * (type /. $BytesPerType)) - 1];
            If[BinaryRead[stream, "Integer8"] === EndOfFile, Throw[None]];
            With[
              {data = Fold[
                 Function[
                   With[
                     {name = #2[[1]],
                      type = #2[[2]]},
                     Replace[
                       BinaryRead[stream, type],
                       {EndOfFile :> #1, x_ :> Append[#1, name -> x]}]]],
                 {},
                 $MGHOptionalData]},
              If[Length[data] < Length[$MGHOptionalData], 
                data,
                Replace[
                  BinaryReadList[stream, "UnsignedInteger8"],
                  {EndOfFile :> data,
                   chars_ :> Append[
                     data,
                     "tags" -> Cases[
                       Map[
                         Apply[StringJoin, FromCharacterCode /@ Cases[#, Except[0]]]&,
                         Split[chars, (#1 != 0 && #2 != 0)&]],
                       s_String /; StringLength[s] > 0]]}]]]]]]]];
ImportMGHData[stream_InputStream, opts___] := "Data" -> With[
  {header = "Header" /. ImportMGHHeader[stream, opts]},
  If[header === $Failed,
    $Failed,
    Replace[
      {"Header" -> header,
       ImportMGHFrames[stream, "Header" -> header, opts],
       ImportMGHFooter[stream, "Header" -> header, opts]},
      l_List /; Position[l, $Failed] != {} -> $Failed]]];
MGHObjectFromData[data_] := With[
  {sym = Unique["mgh"]},
  sym /: MGHQ[sym] = True;
  sym /: ImportedData[sym] = data;
  sym /: Header[sym] = Replace[
    ("Header" /. data),
    "Header" :> (
      Message[MGHObjectFromData::badfmt, "\"Header\" section of data required"];
      Throw[$Failed])];
  Scan[Function[sym /: #[[1]][sym] = #[[2]]], Header[sym]];
  sym /: OptionalData[sym] = Replace[("OptionalData" /. data), "OptionalData" -> None];
  If[OptionalData[sym] =!= None,
    Scan[Function[sym /: OptionalData[sym, #[[1]]] = #[[2]]], OptionalData[sym]]];
  With[
   {width = Dimensions[sym][[1]],
    height = Dimensions[sym][[2]],
    depth = Dimensions[sym][[3]],
    frames = Replace[
      ("Frames" /. data),
      Except[_List] :> (
        Message[MGHObjectFromData::badfmt, "\"Frames\" section of data required"];
        Throw[$Failed])]},
   sym /: Volumes[sym] = frames;
   sym /: Surfaces[sym] := TagSet[
     sym,
     Surfaces[sym],
     (If[Sort[{width, height, depth}][[1;;2]] != {1,1},
        Message[ImportMGH::wrongtype, "volume", "surface"]];
      Flatten /@ frames)];
   sym]];
ImportMGHObject[stream_InputStream, opts___] := Catch[
  With[
    {data = Replace[
       "Data" /. ImportMGHData[stream, opts],
        $Failed :> Throw[$Failed]]},
    MGHObjectFromData[data]]];
ImportMGH[filename_, opts___] := Import[filename, "MGH", opts];
Protect[ImportMGHHeader, ImportMGHFrames, ImportMGHFooter, ImportMGHData, ImportMGHObject, MGHObjectFromData, ImportMGH];
(* Register the importer *)
ImportExport`RegisterImport[
  "MGH",
  {"Header" :> ImportMGHHeader, 
   "Frames" :> ImportMGHFrames,
   "OptionalData" :> ImportMGHFooter,
   "Data" :> ImportMGHData,
   ImportMGHObject},
  "FunctionChannels" -> {"Streams"},
  "BinaryFormat" -> True];

(* Methods for setting an MGH object's data... *)
Unprotect[Volumes, Surfaces];
Volumes /: Set[Volumes[sym_], vol_] := Which[
  !MGHQ[sym], (
    Message[MGHObjectFromData::badfmt, "volumes symbol must be MGH object"];
    $Failed),
  !ListQ[vol] || Length[Dimensions[vol]] != 4, (
    Message[MGHObjectFromData::badfmt, "volumes must be a 4D list"];
    $Failed),
  True, (
    Evaluate[sym] /: Volumes[sym] = vol;
    Evaluate[sym] /: Surfaces[sym] := TagSet[
      Evaluate[sym],
      Surfaces[sym],
      (If[Dimensions[vol][[2]] <= 1 || Dimensions[vol][[3]] > 1 || Dimensions[vol][[4]] > 1,
         Message[ImportMGH::wrongtype, "volume", "surface"]];
       Flatten /@ vol)];
    vol)];
Surfaces /: Set[Surfaces[sym], vol_] := Which[
  !MGHQ[sym], (
    Message[MGHObjectFromData::badfmt, "volumes symbol must be MGH object"];
    $Failed),
  !ListQ[vol] || Length[Dimensions[vol]] != 2, (
    Message[MGHObjectFromData::badfmt, "surfaces must be a 2D list"];
    $Failed),
  True, (
    Evaluate[sym] /: Surfaces[sym] = vol;
    Evaluate[sym] /: Volumes[sym] := TagSet[
      Evaluate[sym],
      Volumes[sym],
      {{#}}& /@ vol];
    vol)];
Protect[Volumes, Surfaces];


(* Exporting MGH files ****************************************************************************)
ExportMGH[filename_, data_, opts___] := Block[
  {$ByteOrdering = 1},
  (* Make sure we have all the data we need... *)
  With[
    {dat = Which[
       MGHQ[data], 
       ImportedData[data],
       Length[Union[Cases[data, (Rule|RuleDelayed)[("Header"|"Frames"),_], {1}][[All, 1]]]] > 1,
       data,
       ListQ[data] && 3 <= Length[Dimensions[data]] <= 4,
       {"Header" -> {DegreesOfFreedom -> 0, Spacings -> {1.0,1.0,1.0}, 
                     VOXToRASMatrix -> {{-1., 0., 0.}, {0., 0., -1.}, {0., 1., 0.}, {0., 0., 0.}}},
        "Frames" -> If[Length[Dimensions[data]] == 3, {data}, data]},
       True,
       (Message[
          ExportMGH::badfmt,
          "Export data must be an MGH object or a list with \"Header\" and \"Frames\" rules"];
        $Failed)],
     outtype = Replace["OutputFormat", Append[Flatten[{opts}], "OutputFormat" -> "Real32"]]},
    If[dat === $Failed,
      $Failed,
      With[
        {fl = OpenWrite[filename, BinaryFormat -> True],
         header = "Header" /. dat,
         frames = If[MGHQ[data], Volumes[data], "Frames" /. dat],
         optdat = "OptionalData" /. dat},
        With[
          {res = Catch[
            If[fl === $Failed,
              Message[ExportMGH::nofile, filename];
              Throw[$Failed]];
            (* Write header... *)
            BinaryWrite[fl, 1, "Integer32"];
            BinaryWrite[fl, Rest[Dimensions[frames]], "Integer32"];
            BinaryWrite[fl, Length[frames], "Integer32"];
            BinaryWrite[fl, outtype /. $MMATypesToMGH, "Integer32"];
            BinaryWrite[fl, DegreesOfFreedom /. header, "Integer32"];
            BinaryWrite[fl, 1, "Integer16"];
            BinaryWrite[fl, Spacings /. header, "Real32"];
            BinaryWrite[fl, Flatten[Transpose[VOXToRASMatrix /. header][[All, 1;;3]]], "Real32"];
            BinaryWrite[fl, Table[0, {$MGHHeaderSize - StreamPosition[fl]}], "Integer8"];
            (* write frames... *)
            BinaryWrite[fl, Flatten[Map[Transpose[#,{3,2,1}]&, frames]], outtype];
            (* Optional data is not currently supported; zeros are written *)
            Scan[BinaryWrite[fl, 0, #[[2]]]&, $MGHOptionalData];
            True]},
          Close[fl];
          If[res === $Failed, $Failed, filename]]]]]];
Protect[ExportMGH];
ImportExport`RegisterExport["MGH", ExportMGH];

(* This just allows MGH objects to be treated as fields when appropriate *)
Unprotect[Field];
Field[m_?MGHQ] := First[Surfaces[m]];
Protect[Field];


(* Importing surface files ************************************************************************)
Unprotect[ImportSurfaceHeader, ImportSurfaceVertices, ImportSurfaceFaces, ImportSurfaceData, ImportSurfaceObject];
ClearAll[ImportSurfaceHeader, ImportSurfaceVertices, ImportSurfaceFaces, ImportSurfaceData, ImportSurfaceObject];
$SurfaceFileTrianglesID = -2;
ImportNLNLTerminatedString[stream_, opts___] := Catch[
  Apply[
    StringJoin,
    Map[
      FromCharacterCode,
      Reap[
        Module[
          {back2 = BinaryRead[stream, "Integer8"],
           back1 = BinaryRead[stream, "Integer8"]},
          While[back1 != 10 || back2 != 10, 
            Sow[back2];
            With[
              {c = BinaryRead[stream, "Integer8"]},
              If[c === EndOfFile,
                Message[ImportSurface::badfmt, "input stream", "string not terminated"];
                Throw[$Failed]];
              back2 = back1;
              back1 = c;]]];
       ][[2,1]]]]];
ImportSurfaceHeader[stream_, opts___] := "Header" -> Block[
  {$ByteOrdering = 1},
  Catch[
    With[
      {format = BinaryRead[stream, "Integer24"]},
      Which[
        format == $SurfaceFileTrianglesID,
        {"SurfaceFormat" -> "Triangle",
         "CreatorString" -> Replace[
           ImportNLNLTerminatedString[stream, opts],
           $Failed :> Throw[$Failed]],
         "VertexCount" -> BinaryRead[stream, "Integer32"],
         "FaceCount" -> BinaryRead[stream, "Integer32"]},
        True,
        (Message[
           ImportSurface::badfmt,
           "input stream",
           "Input type not recognized: "<>ToString[format]];
         Throw[$Failed])]]]];
ImportSurfaceVertices[stream_, opts__] := "Vertices" -> Block[
  {$ByteOrdering = 1},
  Catch[
    With[
      {header = Replace[
         "Header",
         Append[{opts}, "Header" :> ("Header" /. ImportSurfaceHeader[stream, opts])]]},
      If[header === $Failed, Throw[$Failed]];
      Partition[
        Replace[
          BinaryReadList[stream, "Real32", 3 * ("VertexCount" /. header)],
          {EndOfFile :> (
             Message[
               ImportSurface::badfmt,
               "input stream",
               "EndOfFile reached before end of vertices"];
             Throw[$Failed]),
           x_ /; Not[NumberQ[x]] :> (
             Message[
               ImportSurface::badfmt,
               "input stream",
               "NaN read from file while importing vertices"];
             Throw[$Failed])},
          {1}],
        3]]]];
ImportSurfaceFaces[stream_, opts__] := "Faces" -> Block[
  {$ByteOrdering = 1},
  Catch[
    With[
      {header = Replace[
         "Header",
         Append[{opts}, "Header" :> ("Header" /. ImportSurfaceHeader[stream, opts])]]},
      If[header === $Failed, Throw[$Failed]];
      If[!ListQ["Header" /. {opts}],
        Skip[stream, "Real32", 3 * ("VertexCount" /. header)]];
      With[
        {n = ("VertexCount" /. header)},
        Switch[("SurfaceFormat" /. header),
          "Triangle", Partition[
            Replace[
              BinaryReadList[stream, "Integer32", 3 * ("FaceCount" /. header)],
              {EndOfFile :> (
                 Message[
                   ImportSurface::badfmt,
                   "input stream",
                   "EndOfFile reached before end of vertices"];
                 Throw[$Failed]),
               x_Integer /; 0 <= x < n :> (x + 1),
               x_ :> (
                 Message[
                   ImportSurface::badfmt,
                   "input stream",
                   "out of range value ("<>ToString@x<>") read from file while importing vertices"];
                 Throw[$Failed])},
              {1}],
            3],
          _, (
            Message[
              ImportSurface::badfmt,
              "input stream",
              "Unsupported format: " <> ("SurfaceFormat" /. header)];
            Throw[$Failed])]]]]];
ImportSurfaceData[stream_, opts___] := "Data" -> Catch[
  With[
    {header = Replace[
       "Header",
       Append[{opts}, "Header" :> ("Header" /. ImportSurfaceHeader[stream, opts])]]},
    {"Header" -> header,
     "Vertices" -> Replace[
       "Vertices" /. ImportSurfaceVertices[stream, "Header" -> header, opts],
       $Failed :> Throw[$Failed]],
     "Faces" -> Replace[
       "Faces" /. ImportSurfaceFaces[stream, "Header" -> header, opts],
       $Failed :> Throw[$Failed]]}]];
ImportSurfaceObject[stream_, opts___] := Catch[
  With[
    {dat = Replace[
        "Data" /. ImportSurfaceData[stream, opts],
        {$Failed :> Throw[$Failed]}]},
    With[
      {sym = CorticalSurface[
         "Vertices" /. dat,
         Field -> (Field /. Append[{opts}, Field -> None]),
         Faces -> ("Faces" /. dat)]},
      If[sym === $Failed, Throw[$Failed]];
      sym /: HeaderComment[sym] = ("CreatorString" /. ("Header" /. dat));
      sym /: ImportedData[sym] = dat;
      sym]]];
ImportSurface[filename_String, opts___] := Import[filename, "FreeSurferSurface", opts];
Protect[ImportSurface, ImportSurfaceHeader, ImportSurfaceVertices, ImportSurfaceFaces, ImportSurfaceData, ImportSurfaceObject];
(* Register the importer *)
ImportExport`RegisterImport[
  "FreeSurferSurface",
  {"Header" :> ImportSurfaceHeader,
   "Vertices" :> ImportSurfaceVertices,
   "Faces" :> ImportSurfaceFaces,
   "Data" :> ImportSurfaceData,
   ImportSurfaceObject},
  "FunctionChannels" -> {"Streams"},
  "BinaryFormat" -> True];

(* Exporting Surface Format ***********************************************************************)
ExportSurface[filename_, data_, opts___] := Block[
  {$ByteOrdering = 1},
  (* Make sure we have all the data we need... *)
  With[
    {dat = Which[
       SurfaceQ[data] && ValueQ[ImportedData[data]], 
       ImportedData[data],
       3 == Length@Union[
         Cases[
           data,
           (Rule|RuleDelayed)[("Header"|"Vertices"|"Faces"),_], {1}
          ][[All, 1]]],
       dat,
       True,
       (Message[
          ExportSurface::badfmt,
          "Export data must be a surface object or a list with " <>
           "\"Header\", \"Faces\", and \"Vertices\" rules"];
        $Failed)],
     outtype = Replace["OutputFormat", Append[Flatten[{opts}], "OutputFormat" -> "Real32"]],
     createdString = Replace[
       "CreatedString",
       Append[Flatten[{opts}], "CreatedString" -> "Created by fmrilib for Mathematica"]],
     fl = OpenWrite[filename, BinaryFormat -> True]},
    If[fl === $Failed,
      $Failed,
      With[
        {header = "Header" /. dat,
         V = "Vertices" /. dat,
         F = "Face" /. dat},
        With[
          {res = Check[
             BinaryWrite[fl, $SurfaceFileTrianglesID, "Integer24"];
             BinaryWrite[
               fl,
               Join[ToCharacterCode /@ Characters[createdString], {10, 10}],
               "Integer8"];
             BinaryWrite[fl, Length[V], "Integer32"];
             BinaryWrite[fl, Length[F], "Integer32"];
             BinaryWrite[fl, Flatten[V], "Real32"];
             BinaryWrite[fl, Flatten[f], "Integer32"],
             $Failed]},
        Close[fl];
        If[res === $Failed, $Failed, filename]]]]]];
Protect[ExportSurface];
ImportExport`RegisterExport["FreeSurferSurface", ExporSurface];


(* Importing Weights Format ***********************************************************************)
Unprotect[ImportWeightsHeader, ImportWeightsField, ImportWeightsData, ImportWeightsObject];
ClearAll[ImportWeightsHeader, ImportWeightsField, ImportWeightsData, ImportWeightsObject];
ImportWeightsHeader[stream_, opts___] := "Header" -> Block[
  {$ByteOrdering = 1},
  {"Latency" -> BinaryRead[stream, "Integer16"],
   "Count" -> BinaryRead[stream, "Integer24"]}];
ImportWeightsField[stream_, opts___] := "Field" -> Block[
  {$ByteOrdering = 1},
  With[
    {header = Replace[
       "Header",
        Append[{opts}, "Header" :> ("Header" /. ImportWeightsHeader[stream, opts])]]},
    If[header === $Failed,
      $Failed,
      Replace[
        Range["Count" /. header],
        Append[
          Map[
            Function[#[[1]] -> #[[2]]],
            BinaryReadList[stream, {"Integer24", "Real32"}, "Count" /. header]],
          _ -> None],
        {1}]]]];
ImportWeightsData[stream_, opts___] := "Data" -> Catch[
  Block[
    {$ByeOrdering = 1},
    With[
      {header = Replace[
        Replace["Header", 
                Append[{opts}, "Header" :> ("Header" /. ImportWeightsHeader[stream, opts])]],
        $Failed :> Throw[$Failed]]},
      {"Header" -> header,
       "Field" -> Replace[
         "Field" /. ImportWeightsField[stream, "Header" -> header, opts],
         $Failed :> Throw[$Failed]]}]]];
ImportWeightsObject[stream_, opts___] := With[
  {data = "Data" /. ImportWeightsData[stream, opts]},
  If[data === $Failed,
    $Failed,
    With[
      {sym = Unique["weights"]},
      sym /: Field[sym] = ("Field" /. data);
      sym /: Header[sym] = ("Header" /. data);
      sym]]];;
ImportWeights[filename_String, opts___] := Import[filename, "FreeSurferWeights", opts];
Protect[ImportWeights, ImportWeightsHeader, ImportWeightsField, ImportWeightsData, ImportWeightsObject];
(* Register the importer *)
ImportExport`RegisterImport[
  "FreeSurferWeights",
  {"Header" :> ImportWeightsHeader,
   "Field" :> ImportWeightsField,
   "Data" :> ImportWeightsData,
   ImportWeightsObject},
  "FunctionChannels" -> {"Streams"},
  "BinaryFormat" -> True];

(* Exporting weights ******************************************************************************)
ExportWeights[filename_, data_, opts___] := Block[
  {$ByteOrdering = 1},
  (* Make sure we have all the data we need... *)
  With[
    {impdat = Which[
       FieldQ[data] && And@@(And[RuleQ[#], IntegerQ[#[[1]]], NumberQ[#[[2]]]]& /@ Field[data]), {
         "Header" -> {"Latency" -> ("Latency" /. Append[{opts}, "Latency" -> 0]),
                      "Count" -> Length[Field[data]]},
         "Field" -> Field[data]},
       ListQ[data] && And@@(And[RuleQ[#], IntegerQ[#[[1]]], NumberQ[#[[2]]]]& /@ data), {
         "Header" -> {"Latency" -> ("Latency" /. Append[{opts}, "Latency" -> 0]),
                      "Count" -> Length[data]},
         "Field" -> data},
       True, (
         Message[ExportWeights::badfmt, "output must be a field object or list of rules"];
         $Failed)]},
    If[impdat === $Failed,
      $Failed,
      With[
        {header = "Header" /. impdat,
         dat = "Field" /. impdat,
         fl = OpenWrite[filename, BinaryFormat -> True]},
        With[
          {res = Catch[
             If[fl === $Failed, Message[ExportWeights::nofile, filename]; Throw[$Failed]];
             BinaryWrite[fl, "Latency" /. header, "Integer16"];
             BinaryWrite[fl, "Count" /. header, "Integer16"];
             BinaryWrite[fl, (List@@#)& /@ dat, {"Integer24", "Real32"}];
             filename]},
          Close[fl];
          res]]]]];
Protect[ExportWeights];
ImportExport`RegisterExport["FreeSurferWeights", ExporWeights];

(* Importing curv files ***************************************************************************)
Unprotect[ImportCurv, ImportCurvObject, ImportCurvData, ImportCurvField, ImportCurvHeader];
ClearAll[ImportCurv, ImportCurvObject, ImportCurvData, ImportCurvField, ImportCurvHeader];
ImportCurvHeader[stream_, opts___] := "Header" -> Catch[
  Block[
    {$ByteOrdering = 1},
    With[
      {firstInt = BinaryRead[stream, "Integer24"]},
      Which[
        firstInt == -1,
        {"VertexCount" -> BinaryRead[stream, "Integer32"],
         "FaceCount" -> BinaryRead[stream, "Integer32"],
         "ValuesPerVertex" -> BinaryRead[stream, "Integer32"],
         "CurvType" -> "New"},
        IntegerQ[firstInt] && firstInt > 0, 
        {"VertexCount" -> firstInt,
         "FaceCount" -> BinaryRead[stream, "Integer24"],
         "ValuesPerVertex" -> 1,
         "CurvType" -> "Old"},
        True,
        (Message[ImportCurv::badfmt, "input stream", "Could not read header"];
         $Failed)]]]];
ImportCurvField[stream_, opts___] := "Field" -> Catch[
  Block[
    {$ByteOrdering = 1},
    With[
      {header = "Header" /. Append[
          {opts},
          "Header" :> Replace[
            "Header" /. ImportCurvHeader[stream, opts],
            $Failed :> Throw[$Failed]]]},
      Switch["CurvType" /. header,
        "Old", BinaryReadList[stream, "Integer16", "VertexCount" /. header] / 100.0,
        "New", With[
          {k = "ValuesPerVertex" /. header},
          If[k == 1,
            BinaryReadList[stream, "Real32", "VertexCount" /. header],
            BinaryReadList[stream, Table["Real32", {k}], "VertexCount" /. header]]]]]]];
ImportCurvData[stream_, opts___] := "Data" -> Catch[
  With[
    {header = Replace[
       "Header" /. ImportCurvHeader[stream, opts],
       $Failed :> Throw[$Failed]]},
    {"Header" -> header,
     "Field" -> Replace[
       "Field" /. ImportCurvField[stream, "Header" -> header, opts],
       $Failed :> Throw[$Failed]]}]];
ImportCurvObject[stream_, opts___] := Catch[
  With[
    {data = Replace[
       "Data" /. ImportCurvData[stream, opts],
       $Failed :> Throw[$Failed]],
     sym = Unique["curv"]},
    sym /: Field[sym] = "Field" /. data;
    sym /: Header[sym] = "Header" /. data;
    sym]];
ImportCurv[filename_String, opts___] := Import[filename, "FreeSurferCurv", opts];
Protect[ImportCurv, ImportCurvObject, ImportCurvData, ImportCurvField, ImportCurvHeader];
(* Register the importer *)
ImportExport`RegisterImport[
  "FreeSurferCurv",
  {"Header" :> ImportCurvHeader,
   "Field" :> ImportCurvField,
   "Data" :> ImportCurvData,
   ImportCurvObject},
  "FunctionChannels" -> {"Streams"},
  "BinaryFormat" -> True];

(* Exporting curv filess **************************************************************************)
ExportCurv[filename_, data_, opts___] := Block[
  {$ByteOrdering = 1},
  (* Make sure we have all the data we need... *)
  With[
    {impdat = Which[
       And[FieldQ[data], 
           Length[Dimensions[Field[data]]] == 2, 
           And@@(NumberQ /@ Flatten[Field[data]])], {
         "Header" -> {"VertexCount" -> Length[Field[data]],
                      "FaceCount" -> ("FaceCount" /. Append[{opts}, "FaceCount" -> 0]),
                      "ValuesPerVertex" -> Replace[
                        Dimensions[Field[data]], 
                        {{_,x_} :> x, {x_} :> 1}],
                      "CurvType" -> -1},
         "Field" -> Field[data]},
       And[ListQ[data], 
           Length[Dimensions[data]] == 2, 
           And@@(NumberQ /@ Flatten[data])], {
         "Header" -> {"VertexCount" -> Length[data],
                      "FaceCount" -> ("FaceCount" /. Append[{opts}, "FaceCount" -> 0]),
                      "ValuesPerVertex" -> Replace[
                        Dimensions[data], 
                        {{_,x_} :> x, {x_} :> 1}],
                      "CurvType" -> ("CurvType" /. Append[{opts}, "CurvType" -> "New"])},
         "Field" -> data},
       True, (
         Message[ExportCurv::badfmt, "output must be a field object or list or matrix"];
         $Failed)]},
    With[
      {header = "Header" /. impdat,
       dat = "Field" /. impdat},
      Which[
        impdat === $Failed, $Failed,
        ("CurvType" /. header) == "Old" && ("ValuesPerVertex" /. header) > 1, (
          Message[
            ExportCurv::badfmt, 
            "Cannot output old header type with more than 1 value per vertex"];
          $Failed),
        True, With[
           {fl = OpenWrite[filename, BinaryFormat -> True]},
          With[
            {res = Catch[
               If[fl === $Failed, Message[ExportCurv::nofile, filename]; Throw[$Failed]];
               If[("CurvType" /. header) == "Old",
                 (BinaryWrite[fl, Length[dat], "Integer24"];
                  BinaryWrite[fl, "FaceCount" /. header, "Integer24"];
                  BinaryWrite[fl, Round[100.0 * Flatten[dat]], "Integer16"]),
                 (BinaryWrite[fl, -1, "Integer24"];
                  BinaryWrite[fl, "VertexCount" /. header, "Integer32"];
                  BinaryWrite[fl, "FaceCount" /. header, "Integer32"];
                  BinaryWrite[fl, "ValuesPerVertex" /. header, "Integer32"];
                  BinaryWrite[fl, Flatten[dat], "Real32"])];
               BinaryWrite[fl, (List@@#)& /@ dat, {"Integer24", "Real32"}];
               filename]},
            Close[fl];
            res]]]]]];
Protect[ExportCurv];
ImportExport`RegisterExport["FreeSurferCurv", ExporCurv];

(* Importing Annotation files *********************************************************************)
Unprotect[ImportAnnotationData, ImportAnnotationObject];
ClearAll[ImportAnnotationData, ImportAnnotationObject];
ImportAnnotationData[stream_, opts___] := Block[
  {$ByteOrdering = 1},
  Catch[
    With[
      {n = Replace[
         BinaryRead[stream, "Integer32"],
         (x_Integer /; x < 1) :> (
           Message[ImportAnnotation::badfmt, "input stream", "number of vertices < 1"];
           Throw[$Failed])]},
      With[
        {vertices = Map[
           RGBColor[#[[1]]/255.0, #[[2]]/255.0, #[[3]]/255.0, 1.0 - #[[4]]/255.0]&,
           Reverse /@ SortBy[
             BinaryReadList[
               stream, 
               {"Integer32", 
                "UnsignedInteger8", 
                "UnsignedInteger8",
                "UnsignedInteger8",
                "UnsignedInteger8"}, 
               n],
             First
            ][[All, 2;;5]]]},
        With[
          {LUT = If[BinaryRead[stream, "Integer32"] == 0,
             ("LookupTable" /. Append[
                {opts},
                "LookupTable" :> (
                  Message[ImportAnnotation::warning, "Resorting to global label LookupTable"];
                  $FreeSurferColorLUT)]),
             With[
               {maxstruct = (
                  If[BinaryRead[stream, "Integer32"] >= 0,
                    Message[
                      ImportAnnotation::badfmt,
                      "input stream", "old annotation version not supported"];
                    Throw[$Failed]];
                  Replace[
                    BinaryRead[stream, "Integer32"],
                    (x_Integer /; x < 1) :> (
                      Message[ImportAnnotation::badfmt, "input stream", "max label < 1"];
                      Throw[$Failed])]),
                strlen = BinaryRead[stream, "Integer32"]},
               With[
                 {flname = StringJoin@@Most[BinaryReadList[stream, "Character8", strlen]],
                  entries = BinaryRead[stream, "Integer32"]},
                 Dispatch[
                   Flatten[
                     Table[
                       With[
                         {label = BinaryRead[stream, "Integer32"],
                          name = Apply[
                            StringJoin,
                            Most@BinaryReadList[stream, "Character8", BinaryRead[stream, "Integer32"]]],
                          clr = RGBColor@@{
                            BinaryRead[stream, "UnsignedInteger32"] / 255.0,
                            BinaryRead[stream, "UnsignedInteger32"] / 255.0,
                            BinaryRead[stream, "UnsignedInteger32"] / 255.0,
                            1.0 - BinaryRead[stream, "UnsignedInteger32"] / 255.0}},
                         {label -> {name, clr}, name -> {label, clr}, clr -> {label, name}}],
                       {entries}]]]]]]},
          If[LUT === $Failed,
            (Message[ImportAnnotation::badfmt, "input stream", "could not find a color LookupTable"];
             $Failed),
            {"Labels" -> Replace[
               (vertices /. LUT),
               {r_RGBColor :> (
                  Message[
                    ImportAnnotation::badfmt, 
                    "input stream", 
                    "some labels not found in annotation lookup table"];
                  r),
                {x_Integer, s_String} :> x},
               {1}],
             "LookupTable" -> LUT}]]]]]];
ImportAnnotationObject[stream_, opts___] := With[
  {data = ImportAnnotationData[stream, opts],
   sym = Unique["annot"]},
  If[data === $Failed,
    $Failed,
    (sym /: AnnotationQ[sym] = True;
     sym /: ImportedData[sym] = data;
     sym /: Field[sym] = ("Labels" /. data);
     sym /: LookupTable[sym] = ("LookupTable" /. Append[data, "LookupTable" -> None]);
     sym)]];
ImportAnnotation[filename_String, opts___] := Import[filename, "FreeSurferAnnotation", opts];
Protect[ImportAnnotation, ImportAnnotationData, ImportAnnotationObject];
(* Register the importer *)
ImportExport`RegisterImport[
  "FreeSurferAnnotation",
  {"Data" :> ImportAnnotationData,
   ImportAnnotationObject},
  "FunctionChannels" -> {"Streams"},
  "BinaryFormat" -> True];

(* Exporting Annotation filess ********************************************************************)
ExportAnnotation[filename_, data_, opts___] := Block[
  {$ByteOrdering = 1},
  (* Make sure we have all the data we need... *)
  With[
    {impdat = Which[
       And[FieldQ[data], 
           Length[Dimensions[Field[data]]] == 2, 
           And@@(NumberQ /@ Flatten[Field[data]])], {
         "Header" -> {"VertexCount" -> Length[Field[data]],
                      "FaceCount" -> ("FaceCount" /. Append[{opts}, "FaceCount" -> 0]),
                      "ValuesPerVertex" -> Replace[
                        Dimensions[Field[data]], 
                        {{_,x_} :> x, {x_} :> 1}],
                      "CurvType" -> -1},
         "Field" -> Field[data]},
       And[ListQ[data], 
           Length[Dimensions[data]] == 2, 
           And@@(NumberQ /@ Flatten[data])], {
         "Header" -> {"VertexCount" -> Length[data],
                      "FaceCount" -> ("FaceCount" /. Append[{opts}, "FaceCount" -> 0]),
                      "ValuesPerVertex" -> Replace[
                        Dimensions[data], 
                        {{_,x_} :> x, {x_} :> 1}],
                      "CurvType" -> ("CurvType" /. Append[{opts}, "CurvType" -> "New"])},
         "Field" -> data},
       True, (
         Message[ExportCurv::badfmt, "output must be a field object or list or matrix"];
         $Failed)]},
    With[
      {header = "Header" /. impdat,
       dat = "Field" /. impdat},
      Which[
        impdat === $Failed, $Failed,
        ("CurvType" /. header) == "Old" && ("ValuesPerVertex" /. header) > 1, (
          Message[
            ExportCurv::badfmt, 
            "Cannot output old header type with more than 1 value per vertex"];
          $Failed),
        True, With[
           {fl = OpenWrite[filename, BinaryFormat -> True]},
          With[
            {res = Catch[
               If[fl === $Failed, Message[ExportCurv::nofile, filename]; Throw[$Failed]];
               If[("CurvType" /. header) == "Old",
                 (BinaryWrite[fl, Length[dat], "Integer24"];
                  BinaryWrite[fl, "FaceCount" /. header, "Integer24"];
                  BinaryWrite[fl, Round[100.0 * Flatten[dat]], "Integer16"]),
                 (BinaryWrite[fl, -1, "Integer24"];
                  BinaryWrite[fl, "VertexCount" /. header, "Integer32"];
                  BinaryWrite[fl, "FaceCount" /. header, "Integer32"];
                  BinaryWrite[fl, "ValuesPerVertex" /. header, "Integer32"];
                  BinaryWrite[fl, Flatten[dat], "Real32"])];
               BinaryWrite[fl, (List@@#)& /@ dat, {"Integer24", "Real32"}];
               filename]},
            Close[fl];
            res]]]]]];
Protect[ExportAnnotation];;
ImportExport`RegisterExport["FreeSurferAnnotation", ExporAnnotation];

(* Importing Label Files **************************************************************************)
ImportLabelData[filename_, options___] := Check[
  "Data" -> Flatten[
    Last@Reap[
      Replace[
        Import[filename, "Table"],
        {u_Integer /; u > 0, r_, a_, s_, p_} :> Sow[(u + 1) -> {{r,a,s},p}],
        {1}]]],
  $Failed];
ImportLabel[filename_, options___] := Check[
  With[
    {dat = "Data" /. ImportLabelData[filename, options],
     sym = Unique["label"],
     opts = {options},
     inval = True /. Append[{options}, True -> 1]},
    sym /: ImportedData[sym] = dat;
    sym /: Field[sym] = SparseArray[
      Map[(#[[1]] -> inval)&, dat],
      {Max /. Append[opts, Max :> Max[dat[[All,1]]]]},
      False /. Append[opts, False -> 0]];
    sym],
  $Failed];
ImportExport`RegisterImport[
  "FreeSurferLabel",
  {"Data" :> ImportLabelData,
   ImportLabel}];

(* Exporting Label Files **************************************************************************)
ExportLabel[filename_, data, options___] := Check[
  Export[
    filename,
    Prepend[
      Replace[
        If[ListQ[data], data, ImportedData[data]],
        {Rule[id_Integer, {p_, {r_,a_,s}}] :> {id, r, a, s, p},
         _ :> Message[ExportLabel::badfmt]},
        {1}],
      {"#", "Auto-generated", "label", "file", "exported", "from", "Mathematica"}],
    "Table"],
  $Failed];


(* FreeSurfer Directory/filesystem data ***********************************************************)
$FreeSurferHomes = Union[
  Flatten[
    Last[
      Reap[
        Replace[
          {Environment["FREESURFER_HOME"],
           "/usr/local/freesurfer",
           "/opt/freesurfer",
           "/usr/freesurfer",
           "/Applications/freesurfer"},
          s_String /; DirectoryQ[s] && FileExistsQ[s <> "/FreeSurferEnv.sh"] :> Sow[s],
          {1}]]]]];
Protect[$FreeSurferHomes];

$SubjectsDirs = Union[
  Flatten[
    Last[
      Reap[
        Replace[
          Prepend[
            Map[(# <> "/subjects")&, $FreeSurferHomes],
            Environment["SUBJECTS_DIR"]],
          s_String /; DirectoryQ[s] :> Sow[s],
          {1}]]]]];
Protect[$SubjectsDirs];

$Subjects = Union[
  Flatten[
    Map[
      Function[
        Select[
          FileNames[# <> "/*"],
          DirectoryQ]],
      $SubjectsDirs]]];
Protect[$Subjects];

UpdateSubjectsDirs[] := (
  Unprotect[$SubjectsDirs];
  $SubjectsDirs = First/@Gather[
    Flatten[
      Last[
        Reap[
          Replace[
            Append[
              Join[
                Map[(# <> "/subjects")&, $FreeSurferHomes],
                $SubjectsDirs],
              Environment["SUBJECTS_DIR"]],
            s_String /; DirectoryQ[s] :> Sow[s],
            {1}]]]]];
  Protect[$SubjectsDirs]);
UpdateSubjects[] := (
  Unprotect[$Subjects];
  $Subjects = First/@Gather[
    Flatten[
      Join[
        Map[
          Function[
            Select[
              FileNames[# <> "/*"],
              DirectoryQ]],
          $SubjectsDirs],
        $Subjects]]];
  Protect[$Subjects]);


AddFreeSurferHome[s_String] := (
  Unprotect[$FreeSurferHomes];
  $FreeSurferHomes = Prepend[Complement[$FreeSurferHomes, {s}], s];
  UpdateSubjectsDirs[];
  UpdateSubjects[];
  Protect[$FreeSurferHomes]);
RemoveFreeSurferHome[s_] := (
  Unprotect[$FreeSurferHomes];
  $FreeSurferHomes = Complement[$FreeSurferHomes, Cases[$FreeSurferHomes, s]];
  Protect[$FreeSurferHomes]);
AddSubjectsDir[s_String] := (
  Unprotect[$SubjectsDirs];
  $SubjectsDirs = Prepend[Complement[$SubjectsDirs, {s}], s];
  UpdateSubjects[];
  Protect[$SubjectsDirs]);
RemoveSubjectsDir[s_] := (
  Unprotect[$SubjectsDirs];
  $SubjectsDirs = Complement[$SubjectsDirs, Select[$SubjectsDirs, StringMatchQ[#, s]&]];
  Unprotect[$Subjects];
  $Subjects = Select[
    $Subjects,
    Function[
      With[
        {dir = First[StringCases[#, RegularExpression["^(.+)/[^/]"] -> "$1"]]},
        StringMatchQ[dir, s]]]];
  Protect[$Subjets];
  Protect[$SubjectsDir]);
AddSubjects[s_String] := (
  Unprotect[$Subjects];
  $Subjects = Prepend[Complement[$Subjects, {s}], s];
  Protect[$Subjects]);
RemoveSubject[s_] := (
  Unprotect[$Subjects];
  $Subjects = Complement[$Subjects, Cases[$Subjects, s]];
  Protect[$Subjects]);
Protect[AddFreeSurferHome, RemoveFreeSurferHome, 
        AddSubjectsDir, RemoveSubjectsDir, 
        AddSubject, RemoveSubject];


(* Volume labels from the LUT (for use with aseg.mgz) *********************************************)
$FreeSurferColorLUT := With[
  {dat = Catch[
     Scan[
       Function[
         If[FileExistsQ[# <> "/FreeSurferColorLUT.txt"], 
           Throw[
             Select[
               Import[# <> "/FreeSurferColorLUT.txt", "Table"],
               Length[#] == 6 && NumberQ[#[[1]]] && StringQ[#[[2]]] &]]]],
       $FreeSurferHomes];
     $Failed]},
  If[ListQ[dat],
    Set[
      $FreeSurferColorLUT,
      Dispatch[
        Append[
          Flatten[
            Map[
              Function[
                {#[[1]] -> {#[[2]], RGBColor@@Append[#[[3;;5]] / 255.0, 1 - #[[6]]/255.0]}, 
                 #[[2]] -> {#[[1]], RGBColor@@Append[#[[3;;5]] / 255.0, 1 - #[[6]]/255.0]}}],
              dat]],
          _ -> Indeterminate]]],
    (Message[FreeSurfer::nolabel]; $Failed)]];

(* FreeSurfer's Volume data ***********************************************************************)
SubjectSegments[sub_String] := With[
  {dat = Check[
    If[FileExistsQ[sub <> "/mri/aseg.mgh"],
      Import[sub <> "/mri/aseg.mgh",  "MGH"],
      Import[sub <> "/mri/aseg.mgz", {"GZIP", "MGH"}]],
    $Failed]},
   If[dat === $Failed, 
     $Failed, 
     Set[SubjectSegments[sub], First@Volumes[dat]]]];
SubjectSegment[sub_String, label:(_String | _Integer)] := Check[
  With[
    {aseg = SubjectSegments[sub],
     lbl = label /. $FreeSurferColorLUT},
    Which[
      lbl === Indeterminate, Message[SubjectSegment::badlbl, label],
      IntegerQ[label], Position[aseg, label, {3}],
      True, Position[aseg, lbl[[1]], {3}]]],
  $Failed];
SubjectWhiteMatter[sub_String] := With[
  {dat = Check[
     If[FileExistsQ[sub <> "/mri/wm.mgh"],
      Import[sub <> "/mri/wm.mgh",  "MGH"],
      Import[sub <> "/mri/wm.mgz", {"GZIP", "MGH"}]],
    $Failed]},
   If[dat === $Failed, 
     $Failed, 
     Set[SubjectWhiteMatter[sub], First@Volumes[dat]]]];
SubjectBrain[sub_String] := With[
  {dat = Check[
     If[FileExistsQ[sub <> "/mri/brain.mgh"],
      Import[sub <> "/mri/brain.mgh",  "MGH"],
      Import[sub <> "/mri/brain.mgz", {"GZIP", "MGH"}]],
    $Failed]},
   If[dat === $Failed, 
     $Failed, 
     Set[SubjectBrain[sub], First@Volumes[dat]]]];
SubjectFilledBrain[sub_String] := With[
  {dat = Check[
     If[FileExistsQ[sub <> "/mri/brain.mgh"],
      Import[sub <> "/mri/filled.mgh",  "MGH"],
      Import[sub <> "/mri/filled.mgz", {"GZIP", "MGH"}]],
    $Failed]},
   If[dat === $Failed, 
     $Failed, 
     Set[
       SubjectFilledBrain[sub],
       Replace[First@Volumes[dat], {127 -> RH, 255 -> LH}, {3}]]]];
SubjectHemisphere[sub_String, hem:LH|RH] := Check[
  With[
    {dat = SubjectFilledBrain[sub]},
    If[dat =!= $Failed,
      Set[SubjectHemisphere[sub, hem], VolumeMask[VolumeIndices[dat, hem], Dimensions[dat]]]]],
  $Failed];



(* Subject Surface Data ***************************************************************************)
Unprotect[SubjectSimpleSurface, SubjectSimpleCurv, SubjectSimpleWeights, SubjectSimpleAnnot];
ClearAll[SubjectSimpleSurface, SubjectSimpleCurv, SubjectSimpleWeights, SubjectSimpleAnnot];
SubjectSimpleSurface[sub_String /; DirectoryQ[sub], hemi:(LH|RH|RHX), surf_String] := With[
  {dat = Check[
    With[
      {hemistr = Replace[hemi, {(LH|RHX) -> "lh", RH -> "rh"}],
       dirstr = sub <> If[hemi === RHX, "/xhemi/surf/", "/surf/"]},
      Import[dirstr <> hemistr <> "." <> surf, "FreeSurferSurface"]],
    $Failed]},
  If[dat === $Failed, 
    $Failed,
    Set[SubjectSimpleSurface[sub, hemi, surf], dat]]];
SubjectSimpleCurv[sub_String /; DirectoryQ[sub], hemi:(LH|RH|RHX), surf_String] := With[
  {dat = Check[
    With[
      {hemistr = Replace[hemi, {(LH|RHX) -> "lh", RH -> "rh"}],
       dirstr = sub <> If[hemi === RHX, "/xhemi/surf/", "/surf/"]},
      Import[dirstr <> hemistr <> "." <> surf, "FreeSurferCurv"]],
    $Failed]},
  If[dat === $Failed, 
    $Failed,
    Set[SubjectSimpleCurv[sub, hemi, surf], dat]]];
SubjectSimpleWeights[sub_String /; DirectoryQ[sub], hemi:(LH|RH|RHX), surf_String] := With[
  {dat = Check[
    With[
      {hemistr = Replace[hemi, {(LH|RHX) -> "lh", RH -> "rh"}],
       dirstr = sub <> If[hemi === RHX, "/xhemi/surf/", "/surf/"]},
      Import[dirstr <> hemistr <> "." <> surf, "FreeSurferWeights"]],
    $Failed]},
  If[dat === $Failed, 
    $Failed,
    Set[SubjectSimpleWeights[sub, hemi, surf], dat]]];
SubjectSimpleAnnot[sub_String /; DirectoryQ[sub], hemi:(LH|RH|RHX), label_String] := With[
  {dat = Check[
    With[
      {hemistr = Replace[hemi, {(LH|RHX) -> "lh", RH -> "rh"}],
       dirstr = sub <> If[hemi === RHX, "/xhemi/label/", "/label/"]},
      Import[dirstr <> hemistr <> "." <> label, "FreeSurferAnnotation"]],
    $Failed]},
  If[dat === $Failed, 
    $Failed,
    Set[SubjectSimpleAnnot[sub, hemi, surf], dat]]];

(* Subject specific surfaces **********************************************************************)
SubjectOriginalSurface[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleSurface[sub, hemi, "orig"],
  $Failed];
SubjectPialSurface[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleSurface[sub, hemi, "pial"],
  $Failed];
SubjectInflatedSurface[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleSurface[sub, hemi, "inflated"],
  $Failed];
SubjectSphereSurface[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleSurface[sub, hemi, "sphere"],
  $Failed];
SubjectRegisteredSurface[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleSurface[sub, hemi, "sphere.reg"],
  $Failed];
SubjectSymSurface[sub_String, hemi:(LH|RH)] := Check[
  SubjectSimpleSurface[sub, If[hemi === RH, RHX, hemi], "fsaverage_sym.sphere.reg"],
  $Failed];
SubjectFSAverage[sub_String, hemi:(LH|RH)] := Rule[
  SubjectRegisteredSurface[sub, hemi],
  $FSAverageSurface];
SubjectFSAverageSym[sub_String, hemi:(LH|RH)] := Rule[
  SubjectSymSurface[sub, hemi],
  $FSAverageSymSurface];
Protect[SubjectOriginalSurface, SubjectPialSurface, SubjectInflatedSurface, SubjectSphereSurface, 
        SubjectRegisteredSurface, SubjectSymSurface];


(* Data that can be merged with surfaces **********************************************************)
SubjectJacobian[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleCurv[sub, hemi, "jacobian_white"] -> SubjectRegisteredSurface[sub, hemi],
  $Failed];
SubjectCurvature[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleCurv[sub, hemi, "curv"] -> SubjectRegisteredSurface[sub, hemi],
  $Failed];
SubjectSulci[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleCurv[sub, hemi, "sulc"] -> SubjectRegisteredSurface[sub, hemi],
  $Failed];
SubjectThickness[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleCurv[sub, hemi, "thickness"] -> SubjectRegisteredSurface[sub, hemi],
  $Failed];
SubjectRegisteredCurvature[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleCurv[sub, hemi, "avg_curv"] -> SubjectRegisteredSurface[sub, hemi],
  $Failed];
SubjectVertexArea[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleCurv[sub, hemi, "mid.area"] -> SubjectRegisteredSurface[sub, hemi],
  $Failed];
SubjectVertexAreaPial[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleCurv[sub, hemi, "area.pial"] -> SubjectRegisteredSurface[sub, hemi],
  $Failed];
SubjectVertexAreaWhite[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleCurv[sub, hemi, "area.white"] -> SubjectRegisteredSurface[sub, hemi],
  $Failed];
SubjectVertexVolume[sub_String, hemi:(LH|RH|RHX)] := Check[
  SubjectSimpleCurv[sub, hemi, "volume"] -> SubjectRegisteredSurface[sub, hemi],
  $Failed];
SubjectParcellation[sub_String, hemi:(LH|RH|RHX)] := Check[
  With[
    {annot = SubjectSimpleAnnot[sub, hemi, "aparc.a2009s.annot"]},
    With[
      {surf = Rule[
        annot,
        SubjectRegisteredSurface[sub, hemi]]},
      LookupTable[surf] = LookupTable[annot];
      surf]],
  $Failed];
SubjectParcellation2009[sub_String, hemi:(LH|RH|RHX)] := SubjectParcellation[sub, hemi];;
SubjectParcellation2005[sub_String, hemi:(LH|RH|RHX)] := Check[
  With[
    {annot = SubjectSimpleAnnot[sub, hemi, "aparc.annot"]},
    With[
      {surf = Rule[
        annot,
        SubjectRegisteredSurface[sub, hemi]]},
      LookupTable[surf] = LookupTable[annot];
      surf]],
  $Failed];

Protect[SubjectJacobian, SubjectCurvature, SubjectVertexArea, SubjectVertexAreaPial, 
        SubjectVertexAreaWhite, SubjectVolume, SubjectParcellation, SubjectParcellation2009,
        SubjectParcellation2005, SubjectThickness];
    

(* FSAverage and FSAverageSym *********************************************************************)
$FSAverage := With[
  {possibles = Select[$Subjects, (Last[StringSplit[#, $PathnameSeparator]] == "fsaverage")&]},
  If[possibles == {},
    (Message[
       FreeSurfer::notfound,
       "no fsaverage subject found; you may beed to add freesurfer homes or subjects, or set" <> 
        " $FSAverage to the fsaverage subject directory manually."];
     $Failed),
    Set[$FSAverage, First[possibles]]]];
FSAverageOriginalSurface[hemi:(LH|RH)] := Check[
  SubjectSimpleSurface[$FSAverage, hemi, "orig"],
  $Failed];
FSAveragePialSurface[hemi:(LH|RH)] := Check[
  SubjectSimpleSurface[$FSAverage, hemi, "pial"],
  $Failed];
FSAverageInflatedSurface[hemi:(LH|RH)] := Check[
  SubjectSimpleSurface[$FSAverage, hemi, "inflated"],
  $Failed];
FSAverageSphereSurface[hemi:(LH|RH)] := Check[
  SubjectSimpleSurface[$FSAverage, hemi, "sphere"],
  $Failed];
FSAverageCurvature[hemi:(LH|RH)] := Check[
  SubjectSimpleCurv[$FSAverage, hemi, "curv"] -> FSAverageSphereSurface[hemi],
  $Failed];
FSAverageSulci[hemi:(LH|RH)] := Check[
  SubjectSimpleCurv[$FSAverage, hemi, "sulc"] -> FSAverageSphereSurface[hemi],
  $Failed];
FSAverageThickness[hemi:(LH|RH)] := Check[
  SubjectSimpleCurv[$FSAverage, hemi, "thickness"] -> FSAverageSphereSurface[hemi],
  $Failed];
FSAverageRegisteredCurvature[hemi:(LH|RH)] := Check[
  SubjectSimpleCurv[$FSAverage, hemi, "avg_curv"] -> FSAverageSphereSurface[hemi],
  $Failed];
FSAverageVertexArea[hemi:(LH|RH)] := Check[
  SubjectSimpleCurv[$FSAverage, hemi, "mid.area"] -> FSAverageSphereSurface[hemi],
  $Failed];
FSAverageVertexAreaPial[hemi:(LH|RH)] := Check[
  SubjectSimpleCurv[$FSAverage, hemi, "area.pial"] -> FSAverageSphereSurface[hemi],
  $Failed];
FSAverageVertexAreaWhite[hemi:(LH|RH)] := Check[
  SubjectSimpleCurv[$FSAverage, hemi, "area.white"] -> FSAverageSphereSurface[hemi],
  $Failed];
FSAverageVertexVolume[hemi:(LH|RH)] := Check[
  SubjectSimpleCurv[$FSAverage, hemi, "volume"] -> FSAverageSphereSurface[hemi],
  $Failed];
FSAverageParcellation[hemi:(LH|RH)] := Check[
  With[
    {annot = SubjectSimpleAnnot[$FSAverage, hemi, "aparc.a2009s.annot"]},
    With[
      {surf = Rule[
        annot,
        FSAverageSphereSurface[hemi]]},
      LookupTable[surf] = LookupTable[annot];
      surf]],
  $Failed];
FSAverageParcellation2009[hemi:(LH|RH)] := SubjectParcellation[$FSAverage, hemi];;
FSAverageParcellation2005[hemi:(LH|RH)] := Check[
  With[
    {annot = SubjectSimpleAnnot[$FSAverage, hemi, "aparc.annot"]},
    With[
      {surf = Rule[
        annot,
        FSAverageSphereSurface[hemi]]},
      LookupTable[surf] = LookupTable[annot];
      surf]],
  $Failed];
Protect[FSAverageOriginalSurface, FSAveragePialSurface, FSAverageInflatedSurface, FSAverageSphereSurface, FSAverageCurvature, FSAverageSulci, FSAverageThickness, FSAverageRegisteredCurvature, FSAverageVertexArea, FSAverageVertexAreaPial, FSAverageVertexAreaWhite, FSAverageVolume, FSAverageParcellation, FSAverageParcellation2009, FSAverageParcellation2005];

$FSAverageSym := With[
  {possibles = Select[$Subjects, (Last[StringSplit[#, $PathnameSeparator]] == "fsaverage_sym")&]},
  If[possibles == {},
    (Message[
       FreeSurfer::notfound,
       "no fsaverage_sym subject found; you may beed to add freesurfer homes or subjects, or" <> 
        " set $FSAverageSym to the fsaverage_sym subject directory manually."];
     $Failed),
    Set[$FSAverageSym, First[possibles]]]];
FSAverageSymOriginalSurface := Check[
  SubjectSimpleSurface[$FSAverageSym, LH, "orig"],
  $Failed];
FSAverageSymPialSurface := Check[
  SubjectSimpleSurface[$FSAverageSym, LH, "pial"],
  $Failed];
FSAverageSymInflatedSurface := Check[
  SubjectSimpleSurface[$FSAverageSym, LH, "inflated"],
  $Failed];
FSAverageSymSphereSurface := Check[
  SubjectSimpleSurface[$FSAverageSym, LH, "sphere.left_right"],
  $Failed];
FSAverageSymCurvature := Check[
  SubjectSimpleCurv[$FSAverageSym, LH, "curv"] -> FSAverageSymSphereSurface,
  $Failed];
FSAverageSymSulci := Check[
  SubjectSimpleCurv[$FSAverageSym, LH, "sulc"] -> FSAverageSymSphereSurface,
  $Failed];
FSAverageSymThickness := Check[
  SubjectSimpleCurv[$FSAverageSym, LH, "thickness"] -> FSAverageSymSphereSurface,
  $Failed];
FSAverageSymRegisteredCurvature := Check[
  SubjectSimpleCurv[$FSAverageSym, LH, "avg_curv"] -> FSAverageSymSphereSurface,
  $Failed];
FSAverageSymVertexArea := Check[
  SubjectSimpleCurv[$FSAverageSym, LH, "mid.area"] -> FSAverageSymSphereSurface,
  $Failed];
FSAverageSymVertexAreaPial := Check[
  SubjectSimpleCurv[$FSAverageSym, LH, "area.pial"] -> FSAverageSymSphereSurface,
  $Failed];
FSAverageSymVertexAreaWhite := Check[
  SubjectSimpleCurv[$FSAverageSym, LH, "area.white"] -> FSAverageSymSphereSurface,
  $Failed];
FSAverageSymVertexVolume := Check[
  SubjectSimpleCurv[$FSAverageSym, LH, "volume"] -> FSAverageSymSphereSurface,
  $Failed];
FSAverageSymParcellation := Check[
  With[
    {annot = SubjectSimpleAnnot[$FSAverageSym, LH, "aparc.a2009s.annot"]},
    With[
      {surf = Rule[
        annot,
        FSAverageSymSphereSurface]},
      LookupTable[surf] = LookupTable[annot];
      surf]],
  $Failed];
FSAverageSymParcellation2009[hemi:(LH|RH)] := SubjectParcellation[$FSAverageSym, hemi];;
FSAverageSymParcellation2005 := Check[
  With[
    {annot = SubjectSimpleAnnot[$FSAverageSym, LH, "aparc.annot"]},
    With[
      {surf = Rule[
        annot,
        FSAverageSymSphereSurface]},
      LookupTable[surf] = LookupTable[annot];
      surf]],
  $Failed];
Protect[FSAverageSymOriginalSurface, FSAverageSymPialSurface, FSAverageSymInflatedSurface, FSAverageSymSphereSurface, FSAverageSymCurvature, FSAverageSymSulci, FSAverageSymThickness, FSAverageSymRegisteredCurvature, FSAverageSymVertexArea, FSAverageSymVertexAreaPial, FSAverageSymVertexAreaWhite, FSAverageSymVolume, FSAverageSymParcellation, FSAverageSymParcellation2009, FSAverageSymParcellation2005];


End[];
EndPackage[];

