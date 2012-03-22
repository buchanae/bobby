#include <sstream>

#include "api/BamAlignment.h"

#include "types.h"


splat_pos_t getSplatPosition(BamTools::BamAlignment& al) {

    splat_pos_t pos;
    pos.a_start = al.Position + 1;
    pos.a_end = pos.a_start + al.CigarData.at(0).Length - 1;
    pos.b_start = pos.a_end + al.CigarData.at(1).Length + 1;
    pos.b_end = pos.b_start + al.CigarData.at(2).Length - 1;

    return pos;
}

bool isSplat(BamTools::BamAlignment& al) {
    return al.HasTag("XD");
}

splat_t bam2splat(BamTools::BamAlignment& al, BamTools::RefVector& refData) {

    splat_t splat;

    splat.ref = refData.at(al.RefID).RefName;
    al.GetTag("XD", splat.flanks);

    splat.pos = getSplatPosition(al);
    splat.seq = al.QueryBases;

    stringstream readID (stringstream::in | stringstream::out);
    readID << (al.IsReverseStrand())? "-" : "+";
    readID << al.Name;
    splat.readIDs.push_back(readID.str());

    return splat;
}
