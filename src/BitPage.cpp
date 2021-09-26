#include "BitPage.hpp"
#include "moab/Range.hpp"
#include <cstdlib>
#include <cstring>

namespace moab
{

void BitPage::search( unsigned char value, int offset, int count, int per_ent, Range& results,
                      EntityHandle start ) const
{
    const int end        = offset + count;
    Range::iterator hint = results.begin();
    while( offset != end )
    {
        if( get_bits( offset, per_ent ) == value ) hint = results.insert( hint, start );
        ++offset;
        ++start;
    }
}

BitPage::BitPage( int per_ent, unsigned char init_val )
{
    unsigned char mask = (unsigned char)( 1 << per_ent ) - 1;  // 2^per_ent - 1
    init_val &= (unsigned char)mask;
    switch( per_ent )
    {
        default:
            assert( false );
            abort();
            break;  // must be power of two

            // Note: fall through such that all bits in init_val are set, but with odd structure to avoid
	    // fall-through warnings
        case 1:
        case 2:
        case 4:
        case 8:
	  if (1 == per_ent) init_val |= (unsigned char)( init_val << 1 );
	  if (2 >= per_ent) init_val |= (unsigned char)( init_val << 2 );
	  if (4 >= per_ent) init_val |= (unsigned char)( init_val << 4 );
    }
    memset( byteArray, init_val, BitTag::PageSize );
}

}  // namespace moab
