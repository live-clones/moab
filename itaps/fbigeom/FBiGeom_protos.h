#ifndef FBIGEOM_PROTOS_H
#define FBIGEOM_PROTOS_H

#include "moab/MOABConfig.h"

#if defined( MOAB_FC_FUNC_ )
#define ITAPS_FC_WRAPPER MOAB_FC_FUNC_
#elif defined( MOAB_FC_FUNC )
#define ITAPS_FC_WRAPPER MOAB_FC_FUNC
#else
#define ITAPS_FC_WRAPPER( name, NAME ) name
#endif

#define FBiGeom_getDescription  ITAPS_FC_WRAPPER( fbigeom_getdescription, FBIGEOM_GETDESCRIPTION )
#define FBiGeom_getErrorType    ITAPS_FC_WRAPPER( fbigeom_geterrortype, FBIGEOM_GETERRORTYPE )
#define FBiGeom_newGeom         ITAPS_FC_WRAPPER( fbigeom_newgeom, FBIGEOM_NEWGEOM )
#define FBiGeom_dtor            ITAPS_FC_WRAPPER( fbigeom_dtor, FBIGEOM_DTOR )
#define FBiGeom_load            ITAPS_FC_WRAPPER( fbigeom_load, FBIGEOM_LOAD )
#define FBiGeom_save            ITAPS_FC_WRAPPER( fbigeom_save, FBIGEOM_SAVE )
#define FBiGeom_getRootSet      ITAPS_FC_WRAPPER( fbigeom_getrootset, FBIGEOM_GETROOTSET )
#define FBiGeom_getBoundBox     ITAPS_FC_WRAPPER( fbigeom_getboundbox, FBIGEOM_GETBOUNDBOX )
#define FBiGeom_getEntities     ITAPS_FC_WRAPPER( fbigeom_getentities, FBIGEOM_GETENTITIES )
#define FBiGeom_getNumOfType    ITAPS_FC_WRAPPER( fbigeom_getnumoftype, FBIGEOM_GETNUMOFTYPE )
#define FBiGeom_getEntType      ITAPS_FC_WRAPPER( fbigeom_getenttype, FBIGEOM_GETENTTYPE )
#define FBiGeom_getArrType      ITAPS_FC_WRAPPER( fbigeom_getarrtype, FBIGEOM_GETARRTYPE )
#define FBiGeom_getEntAdj       ITAPS_FC_WRAPPER( fbigeom_getentadj, FBIGEOM_GETENTADJ )
#define FBiGeom_getArrAdj       ITAPS_FC_WRAPPER( fbigeom_getarradj, FBIGEOM_GETARRADJ )
#define FBiGeom_getEnt2ndAdj    ITAPS_FC_WRAPPER( fbigeom_getent2ndadj, FBIGEOM_GETENT2NDADJ )
#define FBiGeom_getArr2ndAdj    ITAPS_FC_WRAPPER( fbigeom_getarr2ndadj, FBIGEOM_GETARR2NDADJ )
#define FBiGeom_isEntAdj        ITAPS_FC_WRAPPER( fbigeom_isentadj, FBIGEOM_ISENTADJ )
#define FBiGeom_isArrAdj        ITAPS_FC_WRAPPER( fbigeom_isarradj, FBIGEOM_ISARRADJ )
#define FBiGeom_getTopoLevel    ITAPS_FC_WRAPPER( fbigeom_gettopolevel, FBIGEOM_GETTOPOLEVEL )
#define FBiGeom_getEntClosestPt ITAPS_FC_WRAPPER( fbigeom_getentclosestpt, FBIGEOM_GETENTCLOSESTPT )
#define FBiGeom_getEntClosestPtTrimmed \
    ITAPS_FC_WRAPPER( fbigeom_getentclosestpttrimmed, FBIGEOM_GETENTCLOSESTPTTRIMMED )
#define FBiGeom_getArrClosestPt    ITAPS_FC_WRAPPER( fbigeom_getarrclosestpt, FBIGEOM_GETARRCLOSESTPT )
#define FBiGeom_getEntNrmlXYZ      ITAPS_FC_WRAPPER( fbigeom_getentnrmlxyz, FBIGEOM_GETENTNRMLXYZ )
#define FBiGeom_getArrNrmlXYZ      ITAPS_FC_WRAPPER( fbigeom_getarrnrmlxyz, FBIGEOM_GETARRNRMLXYZ )
#define FBiGeom_getEntNrmlPlXYZ    ITAPS_FC_WRAPPER( fbigeom_getentnrmlplxyz, FBIGEOM_GETENTNRMLPLXYZ )
#define FBiGeom_getArrNrmlPlXYZ    ITAPS_FC_WRAPPER( fbigeom_getarrnrmlplxyz, FBIGEOM_GETARRNRMLPLXYZ )
#define FBiGeom_getEntTgntXYZ      ITAPS_FC_WRAPPER( fbigeom_getenttgntxyz, FBIGEOM_GETENTTGNTXYZ )
#define FBiGeom_getArrTgntXYZ      ITAPS_FC_WRAPPER( fbigeom_getarrtgntxyz, FBIGEOM_GETARRTGNTXYZ )
#define FBiGeom_getFcCvtrXYZ       ITAPS_FC_WRAPPER( fbigeom_getfccvtrxyz, FBIGEOM_GETFCCVTRXYZ )
#define FBiGeom_getEgCvtrXYZ       ITAPS_FC_WRAPPER( fbigeom_getegcvtrxyz, FBIGEOM_GETEGCVTRXYZ )
#define FBiGeom_getEntArrCvtrXYZ   ITAPS_FC_WRAPPER( fbigeom_getentarrcvtrxyz, FBIGEOM_GETENTARRCVTRXYZ )
#define FBiGeom_getEgEvalXYZ       ITAPS_FC_WRAPPER( fbigeom_getegevalxyz, FBIGEOM_GETEGEVALXYZ )
#define FBiGeom_getFcEvalXYZ       ITAPS_FC_WRAPPER( fbigeom_getfcevalxyz, FBIGEOM_GETFCEVALXYZ )
#define FBiGeom_getArrEgEvalXYZ    ITAPS_FC_WRAPPER( fbigeom_getarregevalxyz, FBIGEOM_GETARREGEVALXYZ )
#define FBiGeom_getArrFcEvalXYZ    ITAPS_FC_WRAPPER( fbigeom_getarrfcevalxyz, FBIGEOM_GETARRFCEVALXYZ )
#define FBiGeom_getEntBoundBox     ITAPS_FC_WRAPPER( fbigeom_getentboundbox, FBIGEOM_GETENTBOUNDBOX )
#define FBiGeom_getArrBoundBox     ITAPS_FC_WRAPPER( fbigeom_getarrboundbox, FBIGEOM_GETARRBOUNDBOX )
#define FBiGeom_getVtxCoord        ITAPS_FC_WRAPPER( fbigeom_getvtxcoord, FBIGEOM_GETVTXCOORD )
#define FBiGeom_getVtxArrCoords    ITAPS_FC_WRAPPER( fbigeom_getvtxarrcoords, FBIGEOM_GETVTXARRCOORDS )
#define FBiGeom_getPntRayIntsct    ITAPS_FC_WRAPPER( fbigeom_getpntrayintsct, FBIGEOM_GETPNTRAYINTSCT )
#define FBiGeom_getPntArrRayIntsct ITAPS_FC_WRAPPER( fbigeom_getpntarrrayintsct, FBIGEOM_GETPNTARRRAYINTSCT )
#define FBiGeom_getPntClsf         ITAPS_FC_WRAPPER( fbigeom_getpntclsf, FBIGEOM_GETPNTCLSF )
#define FBiGeom_getPntArrClsf      ITAPS_FC_WRAPPER( fbigeom_getpntarrclsf, FBIGEOM_GETPNTARRCLSF )
#define FBiGeom_getEntNrmlSense    ITAPS_FC_WRAPPER( fbigeom_getentnrmlsense, FBIGEOM_GETENTNRMLSENSE )
#define FBiGeom_getArrNrmlSense    ITAPS_FC_WRAPPER( fbigeom_getarrnrmlsense, FBIGEOM_GETARRNRMLSENSE )
#define FBiGeom_getEgFcSense       ITAPS_FC_WRAPPER( fbigeom_getegfcsense, FBIGEOM_GETEGFCSENSE )
#define FBiGeom_getEgFcArrSense    ITAPS_FC_WRAPPER( fbigeom_getegfcarrsense, FBIGEOM_GETEGFCARRSENSE )
#define FBiGeom_getEgVtxSense      ITAPS_FC_WRAPPER( fbigeom_getegvtxsense, FBIGEOM_GETEGVTXSENSE )
#define FBiGeom_getEgVtxArrSense   ITAPS_FC_WRAPPER( fbigeom_getegvtxarrsense, FBIGEOM_GETEGVTXARRSENSE )
#define FBiGeom_measure            ITAPS_FC_WRAPPER( fbigeom_measure, FBIGEOM_MEASURE )
#define FBiGeom_getFaceType        ITAPS_FC_WRAPPER( fbigeom_getfacetype, FBIGEOM_GETFACETYPE )
#define FBiGeom_getParametric      ITAPS_FC_WRAPPER( fbigeom_getparametric, FBIGEOM_GETPARAMETRIC )
#define FBiGeom_isEntParametric    ITAPS_FC_WRAPPER( fbigeom_isentparametric, FBIGEOM_ISENTPARAMETRIC )
#define FBiGeom_isArrParametric    ITAPS_FC_WRAPPER( fbigeom_isarrparametric, FBIGEOM_ISARRPARAMETRIC )
#define FBiGeom_getEntUVtoXYZ      ITAPS_FC_WRAPPER( fbigeom_getentuvtoxyz, FBIGEOM_GETENTUVTOXYZ )
#define FBiGeom_getArrUVtoXYZ      ITAPS_FC_WRAPPER( fbigeom_getarruvtoxyz, FBIGEOM_GETARRUVTOXYZ )
#define FBiGeom_getEntUtoXYZ       ITAPS_FC_WRAPPER( fbigeom_getentutoxyz, FBIGEOM_GETENTUTOXYZ )
#define FBiGeom_getArrUtoXYZ       ITAPS_FC_WRAPPER( fbigeom_getarrutoxyz, FBIGEOM_GETARRUTOXYZ )
#define FBiGeom_getEntXYZtoUV      ITAPS_FC_WRAPPER( fbigeom_getentxyztouv, FBIGEOM_GETENTXYZTOUV )
#define FBiGeom_getEntXYZtoU       ITAPS_FC_WRAPPER( fbigeom_getentxyztou, FBIGEOM_GETENTXYZTOU )
#define FBiGeom_getArrXYZtoUV      ITAPS_FC_WRAPPER( fbigeom_getarrxyztouv, FBIGEOM_GETARRXYZTOUV )
#define FBiGeom_getArrXYZtoU       ITAPS_FC_WRAPPER( fbigeom_getarrxyztou, FBIGEOM_GETARRXYZTOU )
#define FBiGeom_getEntXYZtoUVHint  ITAPS_FC_WRAPPER( fbigeom_getentxyztouvhint, FBIGEOM_GETENTXYZTOUVHINT )
#define FBiGeom_getArrXYZtoUVHint  ITAPS_FC_WRAPPER( fbigeom_getarrxyztouvhint, FBIGEOM_GETARRXYZTOUVHINT )
#define FBiGeom_getEntUVRange      ITAPS_FC_WRAPPER( fbigeom_getentuvrange, FBIGEOM_GETENTUVRANGE )
#define FBiGeom_getEntURange       ITAPS_FC_WRAPPER( fbigeom_getenturange, FBIGEOM_GETENTURANGE )
#define FBiGeom_getArrUVRange      ITAPS_FC_WRAPPER( fbigeom_getarruvrange, FBIGEOM_GETARRUVRANGE )
#define FBiGeom_getArrURange       ITAPS_FC_WRAPPER( fbigeom_getarrurange, FBIGEOM_GETARRURANGE )
#define FBiGeom_getEntUtoUV        ITAPS_FC_WRAPPER( fbigeom_getentutouv, FBIGEOM_GETENTUTOUV )
#define FBiGeom_getVtxToUV         ITAPS_FC_WRAPPER( fbigeom_getvtxtouv, FBIGEOM_GETVTXTOUV )
#define FBiGeom_getVtxToU          ITAPS_FC_WRAPPER( fbigeom_getvtxtou, FBIGEOM_GETVTXTOU )
#define FBiGeom_getArrUtoUV        ITAPS_FC_WRAPPER( fbigeom_getarrutouv, FBIGEOM_GETARRUTOUV )
#define FBiGeom_getVtxArrToUV      ITAPS_FC_WRAPPER( fbigeom_getvtxarrtouv, FBIGEOM_GETVTXARRTOUV )
#define FBiGeom_getVtxArrToU       ITAPS_FC_WRAPPER( fbigeom_getvtxarrtou, FBIGEOM_GETVTXARRTOU )
#define FBiGeom_getEntNrmlUV       ITAPS_FC_WRAPPER( fbigeom_getentnrmluv, FBIGEOM_GETENTNRMLUV )
#define FBiGeom_getArrNrmlUV       ITAPS_FC_WRAPPER( fbigeom_getarrnrmluv, FBIGEOM_GETARRNRMLUV )
#define FBiGeom_getEntTgntU        ITAPS_FC_WRAPPER( fbigeom_getenttgntu, FBIGEOM_GETENTTGNTU )
#define FBiGeom_getArrTgntU        ITAPS_FC_WRAPPER( fbigeom_getarrtgntu, FBIGEOM_GETARRTGNTU )
#define FBiGeom_getEnt1stDrvt      ITAPS_FC_WRAPPER( fbigeom_getent1stdrvt, FBIGEOM_GETENT1STDRVT )
#define FBiGeom_getArr1stDrvt      ITAPS_FC_WRAPPER( fbigeom_getarr1stdrvt, FBIGEOM_GETARR1STDRVT )
#define FBiGeom_getEnt2ndDrvt      ITAPS_FC_WRAPPER( fbigeom_getent2nddrvt, FBIGEOM_GETENT2NDDRVT )
#define FBiGeom_getArr2ndDrvt      ITAPS_FC_WRAPPER( fbigeom_getarr2nddrvt, FBIGEOM_GETARR2NDDRVT )
#define FBiGeom_getFcCvtrUV        ITAPS_FC_WRAPPER( fbigeom_getfccvtruv, FBIGEOM_GETFCCVTRUV )
#define FBiGeom_getFcArrCvtrUV     ITAPS_FC_WRAPPER( fbigeom_getfcarrcvtruv, FBIGEOM_GETFCARRCVTRUV )
#define FBiGeom_isEntPeriodic      ITAPS_FC_WRAPPER( fbigeom_isentperiodic, FBIGEOM_ISENTPERIODIC )
#define FBiGeom_isArrPeriodic      ITAPS_FC_WRAPPER( fbigeom_isarrperiodic, FBIGEOM_ISARRPERIODIC )
#define FBiGeom_isFcDegenerate     ITAPS_FC_WRAPPER( fbigeom_isfcdegenerate, FBIGEOM_ISFCDEGENERATE )
#define FBiGeom_isFcArrDegenerate  ITAPS_FC_WRAPPER( fbigeom_isfcarrdegenerate, FBIGEOM_ISFCARRDEGENERATE )
#define FBiGeom_getTolerance       ITAPS_FC_WRAPPER( fbigeom_gettolerance, FBIGEOM_GETTOLERANCE )
#define FBiGeom_getEntTolerance    ITAPS_FC_WRAPPER( fbigeom_getenttolerance, FBIGEOM_GETENTTOLERANCE )
#define FBiGeom_getArrTolerance    ITAPS_FC_WRAPPER( fbigeom_getarrtolerance, FBIGEOM_GETARRTOLERANCE )
#define FBiGeom_initEntIter        ITAPS_FC_WRAPPER( fbigeom_initentiter, FBIGEOM_INITENTITER )
#define FBiGeom_initEntArrIter     ITAPS_FC_WRAPPER( fbigeom_initentarriter, FBIGEOM_INITENTARRITER )
#define FBiGeom_getNextEntIter     ITAPS_FC_WRAPPER( fbigeom_getnextentiter, FBIGEOM_GETNEXTENTITER )
#define FBiGeom_getNextEntArrIter  ITAPS_FC_WRAPPER( fbigeom_getnextentarriter, FBIGEOM_GETNEXTENTARRITER )
#define FBiGeom_resetEntIter       ITAPS_FC_WRAPPER( fbigeom_resetentiter, FBIGEOM_RESETENTITER )
#define FBiGeom_resetEntArrIter    ITAPS_FC_WRAPPER( fbigeom_resetentarriter, FBIGEOM_RESETENTARRITER )
#define FBiGeom_endEntIter         ITAPS_FC_WRAPPER( fbigeom_endentiter, FBIGEOM_ENDENTITER )
#define FBiGeom_endEntArrIter      ITAPS_FC_WRAPPER( fbigeom_endentarriter, FBIGEOM_ENDENTARRITER )
#define FBiGeom_copyEnt            ITAPS_FC_WRAPPER( fbigeom_copyent, FBIGEOM_COPYENT )
#define FBiGeom_sweepEntAboutAxis  ITAPS_FC_WRAPPER( fbigeom_sweepentaboutaxis, FBIGEOM_SWEEPENTABOUTAXIS )
#define FBiGeom_deleteAll          ITAPS_FC_WRAPPER( fbigeom_deleteall, FBIGEOM_DELETEALL )
#define FBiGeom_deleteEnt          ITAPS_FC_WRAPPER( fbigeom_deleteent, FBIGEOM_DELETEENT )
#define FBiGeom_createSphere       ITAPS_FC_WRAPPER( fbigeom_createsphere, FBIGEOM_CREATESPHERE )
#define FBiGeom_createPrism        ITAPS_FC_WRAPPER( fbigeom_createprism, FBIGEOM_CREATEPRISM )
#define FBiGeom_createBrick        ITAPS_FC_WRAPPER( fbigeom_createbrick, FBIGEOM_CREATEBRICK )
#define FBiGeom_createCylinder     ITAPS_FC_WRAPPER( fbigeom_createcylinder, FBIGEOM_CREATECYLINDER )
#define FBiGeom_createCone         ITAPS_FC_WRAPPER( fbigeom_createcone, FBIGEOM_CREATECONE )
#define FBiGeom_createTorus        ITAPS_FC_WRAPPER( fbigeom_createtorus, FBIGEOM_CREATETORUS )
#define FBiGeom_moveEnt            ITAPS_FC_WRAPPER( fbigeom_moveent, FBIGEOM_MOVEENT )
#define FBiGeom_rotateEnt          ITAPS_FC_WRAPPER( fbigeom_rotateent, FBIGEOM_ROTATEENT )
#define FBiGeom_reflectEnt         ITAPS_FC_WRAPPER( fbigeom_reflectent, FBIGEOM_REFLECTENT )
#define FBiGeom_scaleEnt           ITAPS_FC_WRAPPER( fbigeom_scaleent, FBIGEOM_SCALEENT )
#define FBiGeom_uniteEnts          ITAPS_FC_WRAPPER( fbigeom_uniteents, FBIGEOM_UNITEENTS )
#define FBiGeom_subtractEnts       ITAPS_FC_WRAPPER( fbigeom_subtractents, FBIGEOM_SUBTRACTENTS )
#define FBiGeom_intersectEnts      ITAPS_FC_WRAPPER( fbigeom_intersectents, FBIGEOM_INTERSECTENTS )
#define FBiGeom_sectionEnt         ITAPS_FC_WRAPPER( fbigeom_sectionent, FBIGEOM_SECTIONENT )
#define FBiGeom_imprintEnts        ITAPS_FC_WRAPPER( fbigeom_imprintents, FBIGEOM_IMPRINTENTS )
#define FBiGeom_mergeEnts          ITAPS_FC_WRAPPER( fbigeom_mergeents, FBIGEOM_MERGEENTS )
#define FBiGeom_createEntSet       ITAPS_FC_WRAPPER( fbigeom_createentset, FBIGEOM_CREATEENTSET )
#define FBiGeom_destroyEntSet      ITAPS_FC_WRAPPER( fbigeom_destroyentset, FBIGEOM_DESTROYENTSET )
#define FBiGeom_isList             ITAPS_FC_WRAPPER( fbigeom_islist, FBIGEOM_ISLIST )
#define FBiGeom_getNumEntSets      ITAPS_FC_WRAPPER( fbigeom_getnumentsets, FBIGEOM_GETNUMENTSETS )
#define FBiGeom_getEntSets         ITAPS_FC_WRAPPER( fbigeom_getentsets, FBIGEOM_GETENTSETS )
#define FBiGeom_addEntToSet        ITAPS_FC_WRAPPER( fbigeom_addenttoset, FBIGEOM_ADDENTTOSET )
#define FBiGeom_rmvEntFromSet      ITAPS_FC_WRAPPER( fbigeom_rmventfromset, FBIGEOM_RMVENTFROMSET )
#define FBiGeom_addEntArrToSet     ITAPS_FC_WRAPPER( fbigeom_addentarrtoset, FBIGEOM_ADDENTARRTOSET )
#define FBiGeom_rmvEntArrFromSet   ITAPS_FC_WRAPPER( fbigeom_rmventarrfromset, FBIGEOM_RMVENTARRFROMSET )
#define FBiGeom_addEntSet          ITAPS_FC_WRAPPER( fbigeom_addentset, FBIGEOM_ADDENTSET )
#define FBiGeom_rmvEntSet          ITAPS_FC_WRAPPER( fbigeom_rmventset, FBIGEOM_RMVENTSET )
#define FBiGeom_isEntContained     ITAPS_FC_WRAPPER( fbigeom_isentcontained, FBIGEOM_ISENTCONTAINED )
#define FBiGeom_isEntArrContained  ITAPS_FC_WRAPPER( fbigeom_isentarrcontained, FBIGEOM_ISENTARRCONTAINED )
#define FBiGeom_isEntSetContained  ITAPS_FC_WRAPPER( fbigeom_isentsetcontained, FBIGEOM_ISENTSETCONTAINED )
#define FBiGeom_addPrntChld        ITAPS_FC_WRAPPER( fbigeom_addprntchld, FBIGEOM_ADDPRNTCHLD )
#define FBiGeom_rmvPrntChld        ITAPS_FC_WRAPPER( fbigeom_rmvprntchld, FBIGEOM_RMVPRNTCHLD )
#define FBiGeom_isChildOf          ITAPS_FC_WRAPPER( fbigeom_ischildof, FBIGEOM_ISCHILDOF )
#define FBiGeom_getNumChld         ITAPS_FC_WRAPPER( fbigeom_getnumchld, FBIGEOM_GETNUMCHLD )
#define FBiGeom_getNumPrnt         ITAPS_FC_WRAPPER( fbigeom_getnumprnt, FBIGEOM_GETNUMPRNT )
#define FBiGeom_getChldn           ITAPS_FC_WRAPPER( fbigeom_getchldn, FBIGEOM_GETCHLDN )
#define FBiGeom_getPrnts           ITAPS_FC_WRAPPER( fbigeom_getprnts, FBIGEOM_GETPRNTS )
#define FBiGeom_createTag          ITAPS_FC_WRAPPER( fbigeom_createtag, FBIGEOM_CREATETAG )
#define FBiGeom_destroyTag         ITAPS_FC_WRAPPER( fbigeom_destroytag, FBIGEOM_DESTROYTAG )
#define FBiGeom_getTagName         ITAPS_FC_WRAPPER( fbigeom_gettagname, FBIGEOM_GETTAGNAME )
#define FBiGeom_getTagSizeValues   ITAPS_FC_WRAPPER( fbigeom_gettagsizevalues, FBIGEOM_GETTAGSIZEVALUES )
#define FBiGeom_getTagSizeBytes    ITAPS_FC_WRAPPER( fbigeom_gettagsizebytes, FBIGEOM_GETTAGSIZEBYTES )
#define FBiGeom_getTagHandle       ITAPS_FC_WRAPPER( fbigeom_gettaghandle, FBIGEOM_GETTAGHANDLE )
#define FBiGeom_getTagType         ITAPS_FC_WRAPPER( fbigeom_gettagtype, FBIGEOM_GETTAGTYPE )
#define FBiGeom_setEntSetData      ITAPS_FC_WRAPPER( fbigeom_setentsetdata, FBIGEOM_SETENTSETDATA )
#define FBiGeom_setEntSetIntData   ITAPS_FC_WRAPPER( fbigeom_setentsetintdata, FBIGEOM_SETENTSETINTDATA )
#define FBiGeom_setEntSetDblData   ITAPS_FC_WRAPPER( fbigeom_setentsetdbldata, FBIGEOM_SETENTSETDBLDATA )
#define FBiGeom_setEntSetEHData    ITAPS_FC_WRAPPER( fbigeom_setentsetehdata, FBIGEOM_SETENTSETEHDATA )
#define FBiGeom_setEntSetESHData   ITAPS_FC_WRAPPER( fbigeom_setentseteshdata, FBIGEOM_SETENTSETESHDATA )
#define FBiGeom_getEntSetData      ITAPS_FC_WRAPPER( fbigeom_getentsetdata, FBIGEOM_GETENTSETDATA )
#define FBiGeom_getEntSetIntData   ITAPS_FC_WRAPPER( fbigeom_getentsetintdata, FBIGEOM_GETENTSETINTDATA )
#define FBiGeom_getEntSetDblData   ITAPS_FC_WRAPPER( fbigeom_getentsetdbldata, FBIGEOM_GETENTSETDBLDATA )
#define FBiGeom_getEntSetEHData    ITAPS_FC_WRAPPER( fbigeom_getentsetehdata, FBIGEOM_GETENTSETEHDATA )
#define FBiGeom_getEntSetESHData   ITAPS_FC_WRAPPER( fbigeom_getentseteshdata, FBIGEOM_GETENTSETESHDATA )
#define FBiGeom_getAllEntSetTags   ITAPS_FC_WRAPPER( fbigeom_getallentsettags, FBIGEOM_GETALLENTSETTAGS )
#define FBiGeom_rmvEntSetTag       ITAPS_FC_WRAPPER( fbigeom_rmventsettag, FBIGEOM_RMVENTSETTAG )
#define FBiGeom_getArrData         ITAPS_FC_WRAPPER( fbigeom_getarrdata, FBIGEOM_GETARRDATA )
#define FBiGeom_getIntArrData      ITAPS_FC_WRAPPER( fbigeom_getintarrdata, FBIGEOM_GETINTARRDATA )
#define FBiGeom_getDblArrData      ITAPS_FC_WRAPPER( fbigeom_getdblarrdata, FBIGEOM_GETDBLARRDATA )
#define FBiGeom_getEHArrData       ITAPS_FC_WRAPPER( fbigeom_geteharrdata, FBIGEOM_GETEHARRDATA )
#define FBiGeom_getESHArrData      ITAPS_FC_WRAPPER( fbigeom_getesharrdata, FBIGEOM_GETESHARRDATA )
#define FBiGeom_setArrData         ITAPS_FC_WRAPPER( fbigeom_setarrdata, FBIGEOM_SETARRDATA )
#define FBiGeom_setIntArrData      ITAPS_FC_WRAPPER( fbigeom_setintarrdata, FBIGEOM_SETINTARRDATA )
#define FBiGeom_setDblArrData      ITAPS_FC_WRAPPER( fbigeom_setdblarrdata, FBIGEOM_SETDBLARRDATA )
#define FBiGeom_setEHArrData       ITAPS_FC_WRAPPER( fbigeom_seteharrdata, FBIGEOM_SETEHARRDATA )
#define FBiGeom_setESHArrData      ITAPS_FC_WRAPPER( fbigeom_setesharrdata, FBIGEOM_SETESHARRDATA )
#define FBiGeom_rmvArrTag          ITAPS_FC_WRAPPER( fbigeom_rmvarrtag, FBIGEOM_RMVARRTAG )
#define FBiGeom_getData            ITAPS_FC_WRAPPER( fbigeom_getdata, FBIGEOM_GETDATA )
#define FBiGeom_getIntData         ITAPS_FC_WRAPPER( fbigeom_getintdata, FBIGEOM_GETINTDATA )
#define FBiGeom_getDblData         ITAPS_FC_WRAPPER( fbigeom_getdbldata, FBIGEOM_GETDBLDATA )
#define FBiGeom_getEHData          ITAPS_FC_WRAPPER( fbigeom_getehdata, FBIGEOM_GETEHDATA )
#define FBiGeom_getESHData         ITAPS_FC_WRAPPER( fbigeom_geteshdata, FBIGEOM_GETESHDATA )
#define FBiGeom_setData            ITAPS_FC_WRAPPER( fbigeom_setdata, FBIGEOM_SETDATA )
#define FBiGeom_setIntData         ITAPS_FC_WRAPPER( fbigeom_setintdata, FBIGEOM_SETINTDATA )
#define FBiGeom_setDblData         ITAPS_FC_WRAPPER( fbigeom_setdbldata, FBIGEOM_SETDBLDATA )
#define FBiGeom_setEHData          ITAPS_FC_WRAPPER( fbigeom_setehdata, FBIGEOM_SETEHDATA )
#define FBiGeom_setESHData         ITAPS_FC_WRAPPER( fbigeom_seteshdata, FBIGEOM_SETESHDATA )
#define FBiGeom_getAllTags         ITAPS_FC_WRAPPER( fbigeom_getalltags, FBIGEOM_GETALLTAGS )
#define FBiGeom_rmvTag             ITAPS_FC_WRAPPER( fbigeom_rmvtag, FBIGEOM_RMVTAG )
#define FBiGeom_subtract           ITAPS_FC_WRAPPER( fbigeom_subtract, FBIGEOM_SUBTRACT )
#define FBiGeom_intersect          ITAPS_FC_WRAPPER( fbigeom_intersect, FBIGEOM_INTERSECT )
#define FBiGeom_unite              ITAPS_FC_WRAPPER( fbigeom_unite, FBIGEOM_UNITE )
#define FBiGeom_getFacets          ITAPS_FC_WRAPPER( fbigeom_getfacets, FBIGEOM_GETFACETS )

#endif
