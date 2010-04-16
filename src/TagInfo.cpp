#include "TagInfo.hpp"

namespace moab {

int TagInfo::size_from_data_type( DataType t )
{
  static const int sizes[] = { 1, 
                               sizeof(int), 
                               sizeof(double), 
                               1, 
                               sizeof(EntityHandle),
                               0 };  
   return sizes[t];
}


void TagInfo::invalidate()
{
  mTagName.clear();
  isValid = false;
  free( mMeshValue );
  free( mDefaultValue );
  mDefaultValue = mMeshValue = 0;
  mDefaultValueSize = mMeshValueSize = 0;
}

  
TagInfo::TagInfo( const TagInfo& copy )
  : mDefaultValue(0),
    mMeshValue(0),
    mDefaultValueSize(copy.mDefaultValueSize),
    mMeshValueSize(copy.mMeshValueSize),
    mDataSize(copy.mDataSize),
    dataType(copy.dataType),
    mTagName(copy.mTagName),
    isValid(copy.isValid)
{
  if (mDefaultValueSize) {
    MALLOC(mDefaultValue, mDefaultValueSize, void* );
    memcpy( mDefaultValue, copy.mDefaultValue, mDefaultValueSize );
  }
  if (mMeshValueSize) {
    MALLOC(mMeshValue, mMeshValueSize, void* );
    memcpy( mMeshValue, copy.mMeshValue, mMeshValueSize );
  }
}

TagInfo& TagInfo::operator=( const TagInfo& copy )
{
  if (copy.mDefaultValue) {
    if (mDefaultValueSize != copy.mDefaultValueSize)
      REALLOC(mDefaultValue, mDefaultValue, 
              std::min(mDefaultValueSize, copy.mDefaultValueSize), 
              copy.mDefaultValueSize, void*);
//      mDefaultValue = realloc( mDefaultValue,copy.mDefaultValueSize);
    mDefaultValueSize = copy.mDefaultValueSize;
    memcpy( mDefaultValue, copy.mDefaultValue, copy.mDefaultValueSize );
  }
  else if (mDefaultValue) {
    free( mDefaultValue );
    mDefaultValue = 0;
    mDefaultValueSize = 0;
  }

  if (copy.mMeshValue) {
    if (mMeshValueSize != copy.mMeshValueSize)
      REALLOC(mMeshValue, mMeshValue, 
              std::min(mMeshValueSize, copy.mMeshValueSize), 
              copy.mMeshValueSize, void*);
//      mMeshValue = realloc( mMeshValue,copy.mMeshValueSize);
    mMeshValueSize = copy.mMeshValueSize;
    memcpy( mMeshValue, copy.mMeshValue, copy.mMeshValueSize );
  }
  else if (mMeshValue) {
    free( mMeshValue );
    mMeshValue = 0;
    mMeshValueSize = 0;
  }
  
  mDataSize = copy.mDataSize;
  dataType = copy.dataType;
  mTagName = copy.mTagName;
  isValid = copy.isValid;
  return *this;
}

void TagInfo::set_mesh_value( const void* data, int size )
{
  if (mMeshValueSize != size) {
    REALLOC(mMeshValue, mMeshValue, 
            std::min(mMeshValueSize, size), size, void*);
    mMeshValueSize = size;
//    mMeshValue = realloc( mMeshValue, size );
  }
  memcpy( mMeshValue, data, size );
}

    //! remove mesh value
void TagInfo::remove_mesh_value() 
{
  free( mMeshValue );
  mMeshValue = 0;
  mMeshValueSize = 0;
}

  
    // Check that all lengths are valid multiples of the type size.
    // Returns true if all lengths are valid, false othersize.
bool TagInfo::check_valid_sizes( const int* sizes, int num_sizes ) const
{
  unsigned sum = 0;
  const unsigned size = size_from_data_type( get_data_type() );
  if (size == 1)
    return true;
  for (int i = 0; i < num_sizes; ++i)
    sum |= ((unsigned)sizes[i]) % size;
  return (sum == 0);
}

} // namespace moab

