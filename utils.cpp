string joinString(const char token, vector <string>& tokens) {
    string data;

    for (int i = 0; i < tokens.size(); i++) {
       data.append(tokens.at(i));
       if (i < tokens.size() - 1) data.append(&token);
    }

    return data;
}
