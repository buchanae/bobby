#include "api/BamMultiReader.h"
#include <iostream>
#include <map>
#include <locale>
#include <stdint.h>
#include "api/BamWriter.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define MIN_GAP 100
#define MAX_GAP 10000

using namespace std;
using namespace BamTools;


typedef int32_t refID_t;

struct group_key_t {
    refID_t refID;
    char mateID;
    bool rev;
    bool operator < ( const group_key_t &other ) const {
        return refID < other.refID 
               || ( refID == other.refID && mateID < other.mateID )
               || ( refID == other.refID && mateID == other.mateID && rev < other.rev );
    }
};

typedef multimap< group_key_t, BamAlignment > group_t;
typedef pair< group_t::iterator, group_t::iterator > group_range_t;


string cigar_to_string( vector<CigarOp>& cd ){
    string out = "";
    vector<CigarOp>::iterator iter = cd.begin();
    for( ; iter != cd.end() ; ++iter ){
        char* len;
        sprintf( len, "%d", iter->Length );
        out.append( len );
        out.push_back( iter->Type );
    }
    return out;
}

vector<CigarOp> string_to_cigar( string& s ){
    vector<CigarOp> ret;
    string::iterator iter = s.begin();
    string cur;
    for( ; iter != s.end(); iter++ ){
        if( !isdigit(*iter) ){
            if( !cur.empty() ){
                CigarOp op;
                op.Type = *iter;
                op.Length = atoi( cur.c_str() );
                ret.push_back( op );
            }
            cur.clear();
        } else {
            cur.push_back( *iter );
        }
    }
    return ret;
}


void _output_valid( BamWriter& writer, group_range_t range_a, group_range_t range_b ){

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

                BamAlignment o;

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
                o.AddTag( "XM", "Z", cigar_to_string( b.CigarData ) );

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
                    o.CigarData = string_to_cigar( cigar );
                    cigars.push_back( cigar_to_string( a.CigarData ) );

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

                o.AddTag( "XC", cigars );
                o.AddTag( "XD", flanks );
                o.AddTag( "XS", seqs );
                o.AddTag( "XT", positions );
                
                writer.SaveAlignment( o );
            }
        }
    }
}

void output_valid( BamWriter& writer, group_t& group, vector< refID_t >& refs ){

    group_key_t key_a;
    group_key_t key_b;
    group_range_t range_a;
    group_range_t range_b;
    vector<refID_t>::iterator refs_it;

    for( refs_it = refs.begin(); refs_it < refs.end(); refs_it++ ){

        key_a.refID = *refs_it;
        key_b.refID = *refs_it;

        key_a.mateID = '1';
        key_b.mateID = '2';

        key_a.rev = true;
        key_b.rev = false;

        _output_valid( writer, group.equal_range( key_a ), group.equal_range( key_b ) );

        key_a.rev = false;
        key_b.rev = true;

        _output_valid( writer, group.equal_range( key_a ), group.equal_range( key_b ) );
    }
}

int main( int argc, char * argv[] ){

    string outputFilename = argv[1];

    // TODO validate args
    // TODO accept min/max gap size args

    vector<string> inputFilenames;
    for( int i = 2; i < argc; i++ ){
        inputFilenames.push_back(argv[i]);
    }

    BamMultiReader reader;
    if( !reader.Open( inputFilenames ) ){
        cerr << reader.GetErrorString();
        return 1;
    }

    BamWriter writer;
    if( !writer.Open( outputFilename, reader.GetHeader(), reader.GetReferenceData() ) ){
        cerr << writer.GetErrorString();
        return 1;
    }

    string current;
    group_t group;
    vector<refID_t> refs;
    BamAlignment a;

    while( reader.GetNextAlignment( a ) ){

        if( current != a.Name.substr( 0, a.Name.length() - 11 ) ){

            output_valid( writer, group, refs );

            group.clear();
            refs.clear();
        }
        
        refs.push_back( a.RefID );

        group_key_t key;
        key.refID = a.RefID;
        key.mateID = a.Name.at( a.Name.length() - 12 );
        key.rev = a.IsReverseStrand();

        group.insert( make_pair( key, a ) );

        current = a.Name.substr( 0, a.Name.length() - 11 );
    }
    output_valid( writer, group, refs );
}
