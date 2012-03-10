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


void _output_valid( group_range_t range_a, group_range_t range_b ){

    group_t::iterator a_it;
    group_t::iterator b_it;

    for( a_it = range_a.first; a_it != range_a.second; ++a_it ){
        for( b_it = range_b.first; b_it != range_b.second; ++b_it ){
            
            BamAlignment a = a_it->second;
            BamAlignment b = b_it->second;

            // ensure 'a' is the most 5' alignment
            if( b.Position < a.Position ){
                swap( a, b );
            }

            int gap = b.Position - ( a.Position + a.Length );

            if( gap >= MIN_GAP && gap <= MAX_GAP ){

                BamAlignment o = new BamAlignment();

                o.RefID = a.RefID;
                o.MateRefID = b.RefID;

                o.MatePosition = b.Position;

                o.MapQuality = a.MapQuality;

                o.SetIsPaired( true );
                o.SetIsReverseStrand( a.IsReverseStrand() );
                o.SetIsMateReverseStrand( b.IsReverseStrand() );
                o.SetIsMapped( true );
                o.SetIsMateMapped( true );

                o.AddTag( "XQ", "Z", b.Name );
                o.AddTag( "R2", "Z", b.QueryBases );
                o.AddTag( "XM", "Z", b.CigarData ); // TODO CigarData to string

                vector< string > cigars;
                vector< string > flanks;
                vector< string > seqs;
                vector< int > positions;

                bool aIsSplat = a.HasTag( "XC" );
                bool bIsSplat = b.HasTag( "XC" );

                if( aIsSplat ){

                    string a_flanks;
                    a.GetTag( "XD", a_flanks );
                    flanks.push_back( a_flanks );

                    int pos;
                    a.GetTag( "XT", pos ); // TODO does bamtools parse out int automatically?
                    o.Position = pos;
                    positions.push_back( a.Position );

                    string seq;
                    a.GetTag( "XS", seq );
                    o.QueryBases = seq;
                    seqs.push_back( a.QueryBases );

                    string cigar;
                    a.GetTag( "XC", cigar );
                    o.CigarData = cigar; // TODO cigar string to CigarData
                    cigars.push_back( a.CigarData ); // TODO CigarData to string

                } else {

                    o.Position = a.Position;
                    o.CigarData = a.CigarData;
                    o.QueryBases = a.QueryBases;
                }

                if( bIsSplat ){

                    string b_flanks;
                    b.GetTag( "XD", b_flanks );
                    flanks.push_back( b_flanks );

                    int pos;
                    b.GetTag( "XT", pos ); // TODO does bamtools parse out int automatically?
                    positions.push_back( pos );

                    string seq;
                    b.GetTag( "XS", seq );
                    seqs.push_back( seq );

                    string cigar;
                    b.GetTag( "XC", cigar );
                    cigars.push_back( cigar );

                }

                // both are regular alignments
                } else {

                    o.AddTag( "R2", "Z", b.QueryBases );
                    o.CigarData = a.CigarData; // TODO
                    o.AddTag( "XM", "Z", b.CigarData ); // TODO CigarData to string
                }

                o.AddTag( "XC", cigars );
                o.AddTag( "XD", flanks );
                o.AddTag( "XS", seqs );
                o.AddTag( "XT", positions );
                
                // TODO pass BamWriter and write this alignment.
                // includes defining header information. how to do this given we're reading multiple files?
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
