class CombinerAlignment : public BamAlignment {

    public:
        splat_key_t splatKey(void);
        bool isSplat(void);
        std::string cigarToString(void);
};
