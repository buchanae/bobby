#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <map>
#include <vector>

#include "SplatPool.h"

using std::vector;
using std::string;

SplatPool::SplatPool () {
    MAX_SIZE = DEFAULT_MAX_SIZE;
}

SplatPool::SplatPool (int m) {
    MAX_SIZE = m;
}

SplatPool::~SplatPool () {

    // remove temporary files
    vector<string>::iterator it = tempFilenames.begin();
    for( ; it != tempFilenames.end(); ++it ){
        remove(it->c_str());
    }
}

void SplatPool::flush (void) {

    std::sort(buffer.begin(), buffer.end(), Sort::ByName);

    // open temp. file
    std::fstream tempfile;
    string path = open_temp("/tmp", tempfile);
    tempFilenames.push_back( path );

    // dump buffer to file
    vector<BamAlignment>::iterator it = buffer.begin();
    for( ; it != buffer.end(); ++it ){
        tempfile << it->str() << std::endl;
    }

    tempfile.close();
    buffer.clear();
}

void SplatPool::add (BamAlignment& splat) {

    if (buffer.size() >= MAX_SIZE) flush();
    // TODO hack alignment name to be used as splat sort key
    buffer.push_back( splat );
}

SplatPool::Reader* SplatPool::reader (void){
    return new SplatPool::Reader(tempFilenames);
}

SplatPool::Reader::Reader (vector<string>& filenames) {

    vector<string>::iterator it = filenames.begin();
    for ( ; it != filenames.end(); ++it ){

        std::fstream f( it->c_str() );
        updateHead( f );
    }
}

SplatPool::Reader::~Reader (){
    // TODO close file handles
}

void SplatPool::Reader::updateHead (std::fstream& f) {
    if (f.good()) {
        string line;
        getline( f, line );

        heads.insert( std::pair<BamAlignment, std::fstream*>( next, &f ) );
    }
}

bool SplatPool::Reader::hasNext (){
    return heads.empty();
}

BamAlignment SplatPool::Reader::read (){
    BamAlignment s = heads.begin()->first;
    std::fstream* f = heads.begin()->second;
    heads.erase( heads.begin() );
    updateHead( *f );
}

std::string open_temp(std::string path, std::fstream& f) {
    path += "/SplatPool-XXXXXX";
    std::vector<char> dst_path(path.begin(), path.end());
    dst_path.push_back('\0');

    int fd = mkstemp(&dst_path[0]);
    if(fd != -1) {
        path.assign(dst_path.begin(), dst_path.end() - 1);
        f.open(path.c_str(), 
               std::ios_base::trunc | std::ios_base::in | std::ios_base::out);
        close(fd);
    }
    return path;
}
