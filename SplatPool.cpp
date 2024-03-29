#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <map>
#include <vector>

#include "api/algorithms/Sort.h"
#include "api/BamAlignment.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"

#include "types.h"
#include "splats.h"
#include "SplatPool.h"

using std::vector;
using std::string;

using BamTools::BamAlignment;
using BamTools::BamMultiReader;
using BamTools::BamWriter;
using BamTools::RefVector;

SplatPool::SplatPool (const RefVector& r) {
    MAX_SIZE = DEFAULT_MAX_SIZE;
    TMP_DIR = DEFAULT_TMP_DIR;
    references = r;
}

SplatPool::SplatPool (int m, const RefVector& r) {
    MAX_SIZE = m;
    TMP_DIR = DEFAULT_TMP_DIR;
    references = r;
}

SplatPool::SplatPool (string tmp, const BamTools::RefVector& r ){
    MAX_SIZE = DEFAULT_MAX_SIZE;
    TMP_DIR = tmp;
    references = r;
}

SplatPool::SplatPool (int m, string tmp, const BamTools::RefVector& r ) {
    MAX_SIZE = m;
    TMP_DIR = tmp;
    references = r;
}

SplatPool::~SplatPool () {

    // remove temporary files
    vector<string>::iterator it = filenames.begin();
    for( ; it != filenames.end(); ++it ){
        remove(it->c_str());
    }
}

void SplatPool::flush (void) {

    std::sort(buffer.begin(), buffer.end(), BamTools::Algorithms::Sort::ByName());

    // open temp. file
    string path = TMP_DIR + "/SplatPool-XXXXXX";
    std::vector<char> dst_path(path.begin(), path.end());
    dst_path.push_back('\0');

    int fd = mkstemp(&dst_path[0]);
    if(fd != -1) {
        path.assign(dst_path.begin(), dst_path.end() - 1);
        close(fd);
    }

    BamWriter writer;
    if(!writer.Open(path, "@HD\tVN:1.0\tSO:queryname", references)){
        std::cerr << writer.GetErrorString() << std::endl;
        // TODO should really throw error here
    }
    filenames.push_back( path );

    // dump buffer to file
    vector<BamAlignment>::iterator it = buffer.begin();
    for( ; it != buffer.end(); ++it ){
        writer.SaveAlignment(*it);
    }

    // TODO does closing the writer allow it to be removed?
    //      need to save reference to file handle?
    writer.Close();
    buffer.clear();
}

void SplatPool::add (BamAlignment& splat) {

    if (!isSplat(splat)) return;

    splat.AddTag("XN", "Z", splat.Name);

    splat_pos_t pos = getSplatPosition(splat);

    stringstream key;
    key << splat.RefID << "-";
    key << pos.a_start << "-" << pos.a_end << "-";
    key << pos.b_start << "-" << pos.b_end;

    splat.Name = key.str();

    if ((int)buffer.size() >= MAX_SIZE) flush();
    buffer.push_back( splat );
}

SplatPool::Reader* SplatPool::reader(void){

    if (buffer.size() > 0) flush();

    Reader* r = new Reader();
    r->Open( filenames );
    return r;
}

bool SplatPool::Reader::GetNextAlignment(BamAlignment& alignment) {
    bool ret = BamMultiReader::GetNextAlignment(alignment);
    if (ret) {
        alignment.GetTag("XN", alignment.Name);
        alignment.RemoveTag("XN");
    }
    return ret;
}
