#ifndef _TYPES_H
#define _TYPES_H

#include <map>
#include <vector>
#include <string>
#include <sstream>

#include "api/BamAlignment.h"

using std::string;
using std::stringstream;
using std::vector;

struct group_key_t {
    int32_t refID;
    char mateID;
    bool rev;
    bool operator< ( const group_key_t &other ) const {
        return refID < other.refID 
               || ( refID == other.refID && mateID < other.mateID )
               || ( refID == other.refID && mateID == other.mateID && rev < other.rev );
    }
};

struct alignment_key_t {
    int32_t refID;
    char mateID;
    bool rev;
    int32_t Position;
    string CigarData;

    bool operator < ( const alignment_key_t &other ) const {
        return refID < other.refID
               || ( refID == other.refID && mateID < other.mateID )
               || ( refID == other.refID && mateID == other.mateID && rev < other.rev )
               || ( refID == other.refID && mateID == other.mateID && rev == other.rev &&                    Position < other.Position)
               || ( refID == other.refID && mateID == other.mateID && rev == other.rev &&                    Position == other.Position && 
                    CigarData.compare(other.CigarData) == -1);
    }
};

struct splat_key_t {
    int32_t ref;
    int a_start;
    int a_end;
    int b_start;
    int b_end;

    bool operator < ( const splat_key_t& other ) const {
        return ref < other.ref
               || ( ref == other.ref && a_start < other.a_start )
               || ( ref == other.ref && a_start == other.a_start && a_end < other.a_end )
               || ( ref == other.ref && a_start == other.a_start 
                    && a_end == other.a_end && b_start < other.b_start )
               || ( ref == other.ref && a_start == other.a_start && a_end == other.a_end 
                    && b_start == other.b_start && b_end < other.b_end );
    }
};

struct splat_t {
    string ref;
    string flanks;
    int a_start;
    int a_end;
    int b_start;
    int b_end;
    string seq;
    vector<std::string> readIDs;

    string str(void) {

        stringstream buffer;

        buffer << ref << "\t";
        buffer << flanks << "\t";
        buffer << a_end - a_start + 1 << "\t";
        buffer << b_end - b_start + 1 << "\t";
        buffer << b_start - a_end << "\t";
        buffer << a_start << "\t" << a_end << "\t";
        buffer << b_start << "\t" << b_end << "\t";
        buffer << seq << "\t";

        stringstream joinedReadIDs;
        vector<string>::iterator r_it = readIDs.begin();
        for( ; r_it != readIDs.end(); ++r_it ){

            if (r_it != readIDs.begin()) joinedReadIDs << ",";
            joinedReadIDs << *r_it;
        }
        
        buffer << readIDs.size() << "\t";
        buffer << joinedReadIDs.str() << "\t";
        buffer.flush();

        return buffer.str();
    }

    bool shouldMerge (splat_t& other) {
        return ref == other.ref && a_start == other.a_start && a_end == other.a_end 
               && b_start == other.b_start && b_end == other.b_end;
    }

    void merge (splat_t& other) {
        readIDs.insert( readIDs.end(), other.readIDs.begin(), other.readIDs.end() );
    }
};

typedef std::multimap< group_key_t, BamTools::BamAlignment > group_t;

typedef std::pair< group_t::iterator, group_t::iterator > group_range_t;

#endif
