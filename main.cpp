#include <iostream>
#include <map>
#include <locale>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "tclap/CmdLine.h"

#include "utils.h"
#include "types.h"
#include "splats.h"
#include "SplatPool.h"

#define VERSION "0.1"

#define DEFAULT_MIN_GAP 100
#define DEFAULT_MAX_GAP 10000

using namespace std;
using namespace BamTools;

int MIN_GAP;
int MAX_GAP;

BamWriter AlignmentsOut;
ofstream splatsOut;

// TODO defaulting to a null stream (e.g. /dev/null or some cross-platform version)
// could clean up the code
bool outputAlignments;
bool outputSplats;

string cigar_to_string( vector<CigarOp>& );

void _output_valid( BamWriter&, group_range_t, group_range_t, map<alignment_key_t, BamAlignment>& );

void output_valid( BamWriter&, SplatPool&, group_t&, map<int32_t,bool>&, RefVector&);

void parseID (string&, string&, char&);

int main( int argc, char * argv[] ){

    vector<string> inputFilenames;
    string combinedOutFilename, alignmentsOutFilename, splatsOutFilename;

    try {
        TCLAP::CmdLine cmd("Program description", ' ', VERSION);

        TCLAP::ValueArg<string> combinedOutputArg("o", "out", 
            "Combined output filename (BAM format)", true, "", "combined.bam");

        TCLAP::ValueArg<string> alignmentsOutputArg("a", "alignments", 
            "Alignments output file name (BAM format)", false, 
            "alignments.bam", "alignments.bam");

        TCLAP::ValueArg<string> splatsOutputArg("s", "splats", 
            "Splat output file name (Splat format)", false, 
            "splats.splat", "splats.splat");

        TCLAP::ValueArg<int> minInsertArg("n", "min-insert", 
            "Minimum insert size", false, DEFAULT_MIN_GAP, "min insert size");

        TCLAP::ValueArg<int> maxInsertArg("x", "max-insert", 
            "Maximum insert size", false, DEFAULT_MAX_GAP, "max insert size");

        TCLAP::UnlabeledMultiArg<string> inputArgs("", "Input files (BAM format)", 
            true, "input.bam");

        cmd.add(combinedOutputArg);
        cmd.add(alignmentsOutputArg);
        cmd.add(splatsOutputArg);
        cmd.add(minInsertArg);
        cmd.add(maxInsertArg);

        // must be added last
        cmd.add(inputArgs);

        cmd.parse(argc, argv);

        combinedOutFilename = combinedOutputArg.getValue();

        if (alignmentsOutputArg.isSet()) {
            outputAlignments = true;
            alignmentsOutFilename = alignmentsOutputArg.getValue();
        }

        if (splatsOutputArg.isSet()) {
            outputSplats = true;
            splatsOutFilename = splatsOutputArg.getValue();
        }

        MIN_GAP = minInsertArg.getValue();
        MAX_GAP = maxInsertArg.getValue();
        inputFilenames = inputArgs.getValue();

    } catch (TCLAP::ArgException &e) {
        cerr << "Error: " << e.error() << " " << e.argId() << endl;
    }

    BamMultiReader reader;
    reader.Open(inputFilenames);

    BamWriter writer;
    if(!writer.Open(combinedOutFilename, reader.GetHeader(), reader.GetReferenceData())){
        cerr << writer.GetErrorString() << endl;
        return 1;
    }

    if (outputAlignments) {
        if (!AlignmentsOut.Open(alignmentsOutFilename, reader.GetHeader(), reader.GetReferenceData())) {
            cerr << AlignmentsOut.GetErrorString() << endl;
            return 1;
        }
    }

    if (outputSplats) {
        splatsOut.open(splatsOutFilename.c_str(), ios::out | ios::trunc);
    }

    string current, prev;
    char mateID;
    group_t group;
    map<int32_t, bool> refs;
    RefVector refData = reader.GetReferenceData();

    // TODO
    SplatPool pool(refData);

    BamAlignment a;
    while(reader.GetNextAlignment(a)){
        parseID(a.Name, current, mateID);

        if(current.compare(prev) && prev.size() > 0){
            output_valid(writer, pool, group, refs, refData);
            group.clear();
            refs.clear();
        }
        
        if(refs.find(a.RefID) == refs.end()) refs.insert(pair <int32_t, bool> (a.RefID, true));

        group_key_t key;
        key.refID = a.RefID;
        key.mateID = mateID;
        key.rev = a.IsReverseStrand();

        group.insert( make_pair( key, a ) );

        prev = current;
    }
    output_valid(writer, pool, group, refs, refData);

    if (outputSplats) {

        SplatPool::Reader* reader = pool.reader();
        RefVector refs = reader->GetReferenceData();

        splat_t prev;
        BamAlignment al;
        while( reader->GetNextAlignment(al) ){

            splat_t splat = bam2splat(al, refs);

            if( prev.ref.empty() ) prev = splat;
            else {
                if( splat.shouldMerge( prev ) ) splat.merge( prev );
                else splatsOut << prev.str() << endl;

                prev = splat;
            }
        }
        // TODO what if prev is empty?
        splatsOut << prev.str() << endl;
        splatsOut.close();
    }

    if (outputAlignments) AlignmentsOut.Close();
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

                // create copies of BamAlignments a and b
                BamAlignment x(a);
                BamAlignment y(b);

                // TODO I don't like that we lose the original read ID
                String baseName;

                parseID(x.Name, baseName, String a);
                x.Name = baseName;

                parseID(y.Name, baseName, String a);
                y.Name = baseName;

                // Add in mate pair info
                x.MateRefID = b.RefID;
                x.MatePosition = b.Position;
                x.InsertSize = gap;

                y.MateRefID = a.RefID;
                y.MatePosition = a.Position;
                y.InsertSize = gap;

                // Update Mapping information appropriately
                x.SetIsPaired(true);
                x.SetIsReverseStrand(a.IsReverseStrand());
                x.SetIsMateReverseStrand(b.IsReverseStrand());
                x.SetIsMapped(true);
                x.SetIsMateMapped(true);
                x.SetIsSecondMate(true);
                x.SetIsProperPair(true);

                y.SetIsPaired(true);
                y.SetIsReverseStrand(b.IsReverseStrand());
                y.SetIsMateReverseStrand(a.IsReverseStrand());
                y.SetIsMapped(true);
                y.SetIsMateMapped(true);
                y.SetIsSecondMate(true);
                y.SetIsProperPair(true);

                //Add Splat tag data
                if (isSplat(a)){
                    string flanks;
                    a.GetTag("XD", flanks);
                    x.AddTag("XD", "Z", flanks);
                }

                if (isSplat(b)){
                    string flanks;
                    b.GetTag("XD", flanks);
                    y.AddTag("XD", "Z", flanks);
                }
 
                writer.SaveAlignment(x);
                writer.SaveAlignment(y);
            }
        }
    }
}

void output_valid( BamWriter& writer, SplatPool& pool, group_t& group, map<int32_t,bool>& refs, RefVector& refData){
    group_key_t key_a;
    group_key_t key_b;
    group_range_t range_a;
    group_range_t range_b;
    map<int32_t,bool>::iterator refs_it;
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

        key_a.rev = false;
        key_b.rev = true;

        _output_valid( writer, group.equal_range( key_a ), group.equal_range( key_b ), good_al );
    }

    if (outputAlignments || outputSplats) {

        for (al_it = good_al.begin(); al_it != good_al.end(); al_it++) {

            if (isSplat(al_it->second)) {
                if (outputSplats) pool.add(al_it->second);
            } else if (outputAlignments) AlignmentsOut.SaveAlignment(al_it->second);
        }
    }
}
