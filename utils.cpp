#include <string>
#include <vector>

#include "utils.h"

std::string joinString(const char token, std::vector<std::string>& tokens) {
    std::string data;

    for (int i = 0; i < tokens.size(); i++) {
       data.append(tokens.at(i));
       if (i < tokens.size() - 1) data.append(&token);
    }

    return data;
}
