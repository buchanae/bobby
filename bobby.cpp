#include "api/BamMultiReader.h"
#include <iostream>
#include <map>
#include <locale>
#include <stdint.h>
#include "api/BamWriter.h"

#define MIN_GAP 100
#define MAX_GAP 10000

using namespace std;
using namespace BamTools;


typedef int32_t refID_t;

struct key_t {
    refID_t refID;
    char mateID;
    bool rev;
    bool operator < ( const key_t &other ) const {
        return refID < other.refID 
               || ( refID == other.refID && mateID < other.mateID )
               || ( refID == other.refID && mateID == other.mateID && rev < other.rev );
    }
};

typedef multimap< key_t, BamAlignment > group_t;
typedef pair< group_t::iterator, group_t::iterator > group_range_t;


int gap_size( BamAlignment& a, BamAlignment& b ){

    if( a.Position < b.Position ){
        return b.Position - ( a.Position + a.Length );
    } else {
        return a.Position - ( b.Position + b.Length );
    }
}

void _output_valid( group_range_t a, group_range_t b ){

    group_t::iterator a_it;
    group_t::iterator b_it;

    for( a_it = a.first; a_it != a.second; ++a_it ){
        for( b_it = b.first; b_it != b.second; ++b_it ){
            
            int gap = gap_size( a_it->second, b_it->second );

            if( gap >= MIN_GAP && gap <= MAX_GAP ){
                // TOOD arrange complete paired BAM record and write
                cout << a_it->second.Name << endl;
            }
        }
    }
}

void output_valid( group_t& group, vector< refID_t >& refs ){

    key_t key_a;
    key_t key_b;
    group_range_t range_a;
    group_range_t range_b;

    vector< refID_t >::iterator refs_it;
    for( refs_it = refs.begin(); refs_it < refs.end(); refs_it++ ){

        key_a.refID = *refs_it;
        key_b.refID = *refs_it;

        key_a.mateID = '1';
        key_b.mateID = '2';

        key_a.rev = true;
        key_b.rev = false;

        _output_valid( group.equal_range( key_a ), group.equal_range( key_b ) );

        key_a.rev = false;
        key_b.rev = true;

        _output_valid( group.equal_range( key_a ), group.equal_range( key_b ) );
    }
}

int main( int argc, char * argv[] ){

    // TODO accept filename args
    // TODO accept min/max gap size args
    vector< string > inputFilenames;
    inputFilenames.push_back("alignments.sorted-by-name.picard.bam");
    inputFilenames.push_back("splats.sorted-by-name.picard.bam");

    BamMultiReader reader;
    if( !reader.Open( inputFilenames ) ){
        cerr << reader.GetErrorString();
        return 1;
    }

    string current;

    group_t group;

    vector< refID_t > refs;

    BamAlignment a;
    while( reader.GetNextAlignment( a ) ){

        if( current != a.Name.substr( 0, a.Name.length() - 11 ) ){

            output_valid( group, refs );

            group.clear();
            refs.clear();
        }
        
        refs.push_back( a.RefID );

        key_t key;
        key.refID = a.RefID;
        key.mateID = a.Name.at( a.Name.length() - 12 );
        key.rev = a.IsReverseStrand();

        group.insert( make_pair( key, a ) );

        current = a.Name.substr( 0, a.Name.length() - 11 );
    }
    output_valid( group, refs );
}
