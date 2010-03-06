#ifndef IMESH_MOAB_HPP
#define IMESH_MOAB_HPP

#include "iMesh.h"
#include "MBForward.hpp"
#include <cstring>
#include <cstdlib>

/* map from MB's entity type to TSTT's entity topology */
extern const iMesh_EntityTopology tstt_topology_table[MBMAXTYPE+1];

/* map from MB's entity type to TSTT's entity type */
extern const iBase_EntityType tstt_type_table[MBMAXTYPE+1];

/* map to MB's entity type from TSTT's entity topology */
extern const MBEntityType mb_topology_table[MBMAXTYPE+1];

/* map from TSTT's tag types to MOAB's */
extern const MBDataType mb_data_type_table[4];

/* map from MOAB's tag types to tstt's */
extern const iBase_TagValueType tstt_data_type_table[MB_MAX_DATA_TYPE+1];

/* map from MOAB's MBErrorCode to tstt's */
extern "C" const iBase_ErrorType iBase_ERROR_MAP[MB_FAILURE+1];

/* Create ITAPS iterator */
iMesh_EntityIterator create_itaps_iterator( MBRange& swap_range,
                                            int array_size = 1 ); 

/* Define macro for quick reference to MBInterface instance */
static inline MBInterface* MBI_cast( iMesh_Instance i )  
  { return reinterpret_cast<MBInterface*>(i); }         
#define MBI MBI_cast(instance)

/* Most recently returned error code */
extern "C" iBase_Error iMesh_LAST_ERROR;

#define RETURN(a) do {iMesh_LAST_ERROR.error_type = *err = (a); \
                      iMesh_LAST_ERROR.description[0] = '\0'; \
                      return;} while(false)

#include "MBCore.hpp"

class MBiMesh : public MBCore
{
private:
  bool haveDeletedEntities;
public:
  MBiMesh();

  virtual ~MBiMesh();
  bool have_deleted_ents( bool reset ) {
    bool result = haveDeletedEntities;
    if (reset)
      haveDeletedEntities = false;
    return result;
  }

  virtual MBErrorCode delete_mesh();
  virtual MBErrorCode delete_entities( const MBEntityHandle*, const int );
  virtual MBErrorCode delete_entities( const MBRange& );
  int AdjTable[16];
};

#define MBimesh reinterpret_cast<MBiMesh*>(MBI)


static inline void
iMesh_processError( int code, const char* desc ) 
{
  std::strncpy( iMesh_LAST_ERROR.description, desc,
                sizeof(iMesh_LAST_ERROR.description) );
  iMesh_LAST_ERROR.error_type = (iBase_ErrorType)code;
}

#define ERROR(CODE,MSG) do { iMesh_setLastError( MBI, (CODE), (MSG) ); RETURN((CODE)); } while(false)
#define IBASE_ERROR(CODE,MSG) iMesh_processError( *err = (CODE), (MSG) )

static inline void iMesh_setLastError( MBInterface*, int code, const char* msg )
  { iMesh_processError( code, msg ); }  
static inline void iMesh_setLastError( MBInterface* mbi, MBErrorCode code, const char* msg )
  { 
    std::string message(msg);
    message += "  (MOAB Error Code: ";
    message += mbi->get_error_string(code);
    message += ")";
    iMesh_processError( iBase_ERROR_MAP[code], message.c_str() ); 
  }

#define CHKERR(CODE,MSG) \
  if (iMesh_isError((CODE))) ERROR((CODE),(MSG))

static inline bool iMesh_isError(int code)
  { return (iBase_SUCCESS != code); }
static inline bool iMesh_isError(MBErrorCode code)
  { return (MB_SUCCESS != code); }


// Check the array size, and allocate the array if necessary.
// Free the array upon leaving scope unless KEEP_ARRAY
// is invoked.
#define ALLOC_CHECK_ARRAY(array, this_size) \
  iMeshArrayManager array ## _manager ( reinterpret_cast<void**>(array), *(array ## _allocated), *(array ## _size), this_size, sizeof(**array), err ); \
  if (iBase_SUCCESS != *err) return

#define ALLOC_CHECK_TAG_ARRAY(array, this_size) \
  iMeshArrayManager array ## _manager ( reinterpret_cast<void**>(array), *(array ## _allocated), *(array ## _size), this_size, 1, err ); \
  if (iBase_SUCCESS != *err) return

#define KEEP_ARRAY(array) \
  array ## _manager .keep_array()

// Check the array size, and allocate the array if necessary.
// Do NOT free the array upon leaving scope.
#define ALLOC_CHECK_ARRAY_NOFAIL(array, this_size) \
  ALLOC_CHECK_ARRAY(array, this_size); KEEP_ARRAY(array)


// Implement RAII pattern for allocated arrays
class iMeshArrayManager
{
  void** arrayPtr;

public:


  iMeshArrayManager( void** array_ptr,
                     int& array_allocated_space,
                     int& array_size,
                     int count,
                     int val_size,
                     int* err ) : arrayPtr(0)
  {
    if (!*array_ptr) {
      *array_ptr = std::malloc(val_size * count);
      array_allocated_space = array_size = count;
      if (!*array_ptr) {
        IBASE_ERROR(iBase_MEMORY_ALLOCATION_FAILED, "Couldn't allocate array.");
      }
      arrayPtr = array_ptr;
    }
    else {
      array_size = count;
      if (array_allocated_space < count) {
        IBASE_ERROR(iBase_BAD_ARRAY_DIMENSION, 
          "Allocated array not large enough to hold returned contents.");
      }
    }
    RETURN(iBase_SUCCESS);
  }
  
  ~iMeshArrayManager() 
  {
    if (arrayPtr) {
      std::free(*arrayPtr);
      *arrayPtr = 0;
    }
  }
  
  void keep_array()
    { arrayPtr = 0; }
};

#endif // IMESH_MOAB_HPP
