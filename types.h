#ifndef _TYPES_H
#define _TYPES_H

#include <map>
#include <vector>
#include <string>
#include <sstream>

#include "api/BamAlignment.h"

#include "utils.h"

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

struct splat_pos_t {
    int a_start;
    int a_end;
    int b_start;
    int b_end;

    bool operator== ( const splat_pos_t& other ) const {
        return a_start == other.a_start && a_end == other.a_end 
               && b_start == other.b_start && b_end == other.b_end;
    }

};

struct splat_t {
    string ref;
    string flanks;
    splat_pos_t pos;
    string seq;
    vector<std::string> readIDs;

    string str(void) {

        stringstream buffer;

        buffer << ref << "\t";
        buffer << flanks << "\t";
        buffer << pos.a_end - pos.a_start + 1 << "\t";
        buffer << pos.b_end - pos.b_start + 1 << "\t";
        buffer << pos.b_start - pos.a_end << "\t";
        buffer << pos.a_start << "\t" << pos.a_end << "\t";
        buffer << pos.b_start << "\t" << pos.b_end << "\t";
        buffer << seq << "\t";
        buffer << joinString(',', readIDs);
        buffer.flush();

        return buffer.str();
    }

    bool shouldMerge (splat_t& other) {
        return ref == other.ref && pos == other.pos;
    }

    void merge (splat_t& other) {
        readIDs.insert( readIDs.end(), other.readIDs.begin(), other.readIDs.end() );
    }
};

typedef std::multimap< group_key_t, BamTools::BamAlignment > group_t;

typedef std::pair< group_t::iterator, group_t::iterator > group_range_t;

#endif
