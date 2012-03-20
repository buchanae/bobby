#include <iostream>
#include <map>
#include <locale>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

using namespace std;
using namespace BamTools;

int MIN_GAP = 10;
int MAX_GAP = 10000;

BamWriter Alignments;
ofstream splat_fh("splats", ios::out | ios::trunc);

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

struct alignment_key_t {
    refID_t refID;
    char mateID;
    bool rev;
    int32_t Position;
    string CigarData;
    bool operator < ( const alignment_key_t &other ) const {
        return refID < other.refID
               || ( refID == other.refID && mateID < other.mateID )
               || ( refID == other.refID && mateID == other.mateID && rev < other.rev )
               || ( refID == other.refID && mateID == other.mateID && rev == other.rev && Position < other.Position)
               || ( refID == other.refID && mateID == other.mateID && rev == other.rev && Position == other.Position && CigarData.compare(other.CigarData) == -1);
    }
};

typedef multimap< group_key_t, BamAlignment > group_t;

typedef pair< group_t::iterator, group_t::iterator > group_range_t;

string cigar_to_string( vector<CigarOp>& );

void _output_valid( BamWriter&, group_range_t, group_range_t, map<alignment_key_t, BamAlignment>& );

void output_valid( BamWriter&, group_t&, map<refID_t,bool>&, RefVector&);

void parseID (string&, string&, char&);

string bam2splat(BamAlignment &, RefVector&);

string joinString(const char, vector <string>&);

int main( int argc, char * argv[] ){
    string outputFilename, finalAlignments;

    if (argc > 2) {
        outputFilename = argv[1];
        finalAlignments = "aligned.bam";
    } else {
        cerr << "Error - you must give an output file name followed by 1 or more bam files." << std::endl;
        return 1;
    }

    // TODO validate args
    // TODO accept min/max gap size args

    vector<string> inputFilenames;
    for (int i = 3; i <= argc; i++)
        inputFilenames.push_back(argv[i - 1]);

    BamMultiReader reader;
    reader.Open(inputFilenames);

    BamWriter writer;
    if(!writer.Open(outputFilename, reader.GetHeader(), reader.GetReferenceData())){
        cerr << writer.GetErrorString() << endl;
        return 1;
    }
    if (!Alignments.Open(finalAlignments, reader.GetHeader(), reader.GetReferenceData())) {
        cerr << Alignments.GetErrorString() << endl;
        return 1;
    }

    string current, prev;
    char mateID;
    group_t group;
    map<refID_t, bool> refs;
    RefVector refData = reader.GetReferenceData();

    BamAlignment a;
    while(reader.GetNextAlignment(a)){
        parseID(a.Name, current, mateID);

        if(current.compare(prev) && prev.size() > 0){
            output_valid(writer, group, refs, refData);
            group.clear();
            refs.clear();
        }
        
        if(refs.find(a.RefID) == refs.end()) refs.insert(pair <refID_t, bool> (a.RefID, true));

        group_key_t key;
        key.refID = a.RefID;
        key.mateID = mateID;
        key.rev = a.IsReverseStrand();

        group.insert( make_pair( key, a ) );

        prev = current;
    }
    output_valid(writer, group, refs, refData);

    Alignments.Close();
    splat_fh.close();
}

void parseID (string& id, string& groupID, char& mateID) {
    groupID.clear();

    int i = 0;
    while (i < id.size()) {
        if (id.at(i) != ' ') groupID += id.at(i);
        else {
            mateID = (i + 1 < id.size()) ? id.at(i + 1) : '0';
            i = id.size();
        }
        i++;
    }
}

string cigar_to_string( vector<CigarOp>& cd ){
    stringstream ss;
    string out = "";
    vector<CigarOp>::iterator iter = cd.begin();
    for( ; iter != cd.end() ; ++iter ){
        ss << iter->Length << iter->Type;
    }
    ss >> out;
    return out;
}

string bam2splat(BamAlignment& al, RefVector& refData) {
    stringstream buffer (stringstream::in | stringstream::out);

    string flanking, strand;
    int a_start, a_end, b_start, b_end;

    al.GetTag("XD", flanking);
    a_start = al.Position + 1;
    a_end = a_start + al.CigarData.at(0).Length - 1;
    b_start = a_end + al.CigarData.at(1).Length + 1;
    b_end = b_start + al.CigarData.at(2).Length - 1;
    strand = (al.IsReverseStrand())? "-" : "+";

    buffer << refData.at(al.RefID).RefName << "\t";
    buffer << al.CigarData.at(0).Length << "\t";
    buffer << al.CigarData.at(2).Length << "\t";
    buffer << al.CigarData.at(1).Length << "\t";
    buffer << a_start << "\t" << a_end << "\t";
    buffer << b_start << "\t" << b_end << "\t";
    buffer << al.QueryBases << "\t1\t" << strand << al.Name;
    buffer.flush();

    return buffer.str();
}

void _output_valid( BamWriter& writer, group_range_t range_a, group_range_t range_b , map<alignment_key_t, BamAlignment> &alignments){
    group_t::iterator a_it;
    group_t::iterator b_it;

    for( a_it = range_a.first; a_it != range_a.second; ++a_it ){
        for( b_it = range_b.first; b_it != range_b.second; ++b_it ){
            group_key_t a_key = a_it->first;
            BamAlignment a = a_it->second;

            group_key_t b_key = b_it->first;
            BamAlignment b = b_it->second;

            // ensure 'a' is the most 5' alignment
            if( b.Position < a.Position ) {
                swap(a, b);
                swap(a_key, b_key);  // Ensure that we can check that alignments are oriented correctly (5'--------->3' 3'<---------5')
            }

            int gap = 0;

            // Calculate paired-end gapped alignments using the CigarOp data - NOT the length of the query sequence
            int a_start = a.Position, a_end, a_length = 0, b_start = b.Position, b_end, b_length = 0;
            for (int i = 0; i < a.CigarData.size(); i++)
                a_length += a.CigarData.at(i).Length;
            a_end = a_start + a_length - 1;
            for (int i = 0; i < b.CigarData.size(); i++)
                b_length += b.CigarData.at(i).Length;
            b_end = b_start + b_length - 1;

            // calculate gap only if fragments do not overlap
            if (a_end < b_start) gap = (b_start - a_end + 1) - 2; // Calculate the correct size of the insert

            if(gap >= MIN_GAP && gap <= MAX_GAP && !a_key.rev && b_key.rev){
            //if(gap >= MIN_GAP && gap <= MAX_GAP){
                alignment_key_t a_t;
                a_t.refID = a.RefID;
                a_t.mateID = a_key.mateID;
                a_t.rev = a.IsReverseStrand();
                a_t.Position = a.Position;
                a_t.CigarData = cigar_to_string(a.CigarData);
                if (alignments.find(a_t) == alignments.end()) alignments.insert(pair <alignment_key_t, BamAlignment> (a_t, a));

                alignment_key_t b_t;
                b_t.refID = b.RefID;
                b_t.mateID = b_key.mateID;
                b_t.rev = b.IsReverseStrand();
                b_t.Position = b.Position;
                b_t.CigarData = cigar_to_string(b.CigarData);
                if (alignments.find(b_t) == alignments.end()) alignments.insert(pair <alignment_key_t, BamAlignment> (b_t, b));

                BamAlignment o;

                // Clone main alignment data
                o.Name = a.Name;
                o.AlignmentFlag = a.AlignmentFlag;
                o.RefID = a.RefID;
                o.Position = a.Position;
                o.MapQuality = 255; // Do not keep track of quality data...
                o.CigarData = a.CigarData;
                o.Length = a.Length;
                o.QueryBases = a.QueryBases;
                o.Qualities = "*"; // Again - do no keep track of quality data...
                o.Bin = a.Bin;

                // Add in mate pair info
                o.MateRefID = b.RefID;
                o.MatePosition = b.Position;
                o.InsertSize = gap;

                // Update Mapping information appropriately
                o.SetIsPaired(true);
                o.SetIsReverseStrand(a.IsReverseStrand());
                o.SetIsMateReverseStrand(b.IsReverseStrand());
                o.SetIsMapped(true);
                o.SetIsMateMapped(true);

                // Keep track of all data of mate pair
                o.AddTag("XQ", "Z", b.Name);
                o.AddTag("R2", "Z", b.QueryBases);
                o.AddTag("XM", "Z", cigar_to_string(b.CigarData));

                //Add Splat tag data
                if (a.HasTag("XD") || b.HasTag("XD")) {
                    vector< string > flanks;
                    string a_flanks;
                    if (a.HasTag("XD")) {
                        a.GetTag("XD", a_flanks);
                        flanks.push_back(a_flanks);
                    }
                    if (b.HasTag("XD")) {
                        b.GetTag("XD", a_flanks );
                        flanks.push_back(a_flanks);
                    }
                    o.AddTag("XD", "Z", joinString(',', flanks));
                }
 
                writer.SaveAlignment(o);
            }
        }
    }
}

string joinString(const char token, vector <string>& tokens) {
    string data;

    for (int i = 0; i < tokens.size(); i++) {
       data.append(tokens.at(i));
       if (i < tokens.size() - 1) data.append(&token);
    }

    return data;
}

void output_valid( BamWriter& writer, group_t& group, map<refID_t,bool>& refs, RefVector& refData){
    group_key_t key_a;
    group_key_t key_b;
    group_range_t range_a;
    group_range_t range_b;
    map<refID_t,bool>::iterator refs_it;
    map<alignment_key_t,BamAlignment> good_al;
    map<alignment_key_t,BamAlignment>::iterator al_it;

    for(refs_it = refs.begin(); refs_it != refs.end(); refs_it++ ){
        key_a.refID = refs_it->first;
        key_b.refID = refs_it->first;

        key_a.mateID = '1';
        key_b.mateID = '2';

        key_a.rev = true;
        key_b.rev = false;

        _output_valid( writer, group.equal_range( key_a ), group.equal_range( key_b ), good_al );

        key_a.rev = true;
        key_b.rev = true;

        //_output_valid( writer, group.equal_range( key_a ), group.equal_range( key_b ), good_al );

        key_a.rev = false;
        key_b.rev = true;

        _output_valid( writer, group.equal_range( key_a ), group.equal_range( key_b ), good_al );

        key_a.rev = false;
        key_b.rev = false;

        //_output_valid( writer, group.equal_range( key_a ), group.equal_range( key_b ), good_al );
    }

    for (al_it = good_al.begin(); al_it != good_al.end(); al_it++) {
        if (al_it->second.CigarData.size() == 3) {
            string splat = bam2splat(al_it->second, refData);
            splat_fh << splat << endl;
        } else {
            Alignments.SaveAlignment(al_it->second);
        }
    }
}
