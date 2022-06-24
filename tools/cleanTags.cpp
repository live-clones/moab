/*
 * cleanTags.cpp
 * this tool will remove a list of tags from a file, or keep only the tags from a given list
 * will use just the names of the tags
 *
 * example of usage:
 * ./mbcleantags -i input_file.h5m -o outfile.h5m -d list_tags_separated_by':' -k list_tags separated by ':'
 *
 */
#include "moab/MOABConfig.h"

#include "moab/ProgOptions.hpp"
#include "moab/Core.hpp"

using namespace moab;
using namespace std;

vector< string > split( const string& i_str, const string& i_delim )
{
    vector< string > result;

    size_t found      = i_str.find( i_delim );
    size_t startIndex = 0;

    while( found != string::npos )
    {
        result.push_back( string( i_str.begin() + startIndex, i_str.begin() + found ) );
        startIndex = found + i_delim.size();
        found      = i_str.find( i_delim, startIndex );
    }
    if( startIndex != i_str.size() ) result.push_back( string( i_str.begin() + startIndex, i_str.end() ) );
    return result;
}

int main( int argc, char* argv[] )
{

    ProgOptions opts;

    string inputfile, outputfile, deleteTags, keepTags;
    opts.addOpt< string >( "input,i", "input filename ", &inputfile );
    opts.addOpt< string >( "output,o", "output file", &outputfile );

    opts.addOpt< string >( "deleteTags,d", "delete tags ", &deleteTags );
    opts.addOpt< string >( "keepTags,k", "keep tags ", &keepTags );
    opts.parseCommandLine( argc, argv );

    Core core;
    Interface* mb = &core;
    ErrorCode rval;
    rval = mb->load_file( inputfile.c_str() );MB_CHK_ERR( rval );
    vector< Tag > existingTags;
    rval = mb->tag_get_tags( existingTags );MB_CHK_ERR( rval );
    vector< string > tagsToDelete;
    if( !keepTags.empty() )
    {
        vector< string > tagsToKeep = split( keepTags, string( ":" ) );
        for( size_t i = 0; i < existingTags.size(); i++ )
        {
            string tname;
            rval = mb->tag_get_name( existingTags[i], tname );MB_CHK_ERR( rval );
            bool deleteTag = false;
            for( size_t k = 0; k < tagsToKeep.size() && !deleteTag; k++ )
            {
                if( tname.compare( tagsToKeep[k] ) == 0 ) deleteTag = true;
            }
            if( !deleteTag ) tagsToDelete.push_back( tname );
        }
    }
    if( !deleteTags.empty() )
    {
        tagsToDelete = split( deleteTags, string( ":" ) );
    }
    for( size_t i = 0; i < tagsToDelete.size(); i++ )
    {
        Tag tag;
        rval = mb->tag_get_handle( tagsToDelete[i].c_str(), tag );
        if( rval == MB_SUCCESS && tag != NULL )
        {
            rval = mb->tag_delete( tag );MB_CHK_ERR( rval );
        }
    }
    cout << "write file " << outputfile << endl;
    rval = mb->write_file( outputfile.c_str() );MB_CHK_ERR( rval );

    return 0;
}
