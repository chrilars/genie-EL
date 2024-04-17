std::vector<size_t> preprocess(const std::vector<bool>& label) {
    std::vector<size_t> converted;
    for (size_t i = 0; i < label.size(); ++i) {
        if (label[i]) converted.push(i);
    }
    return converted;
}

template <class UBDD>
UBDD setToBdd(const UBDD &base,
              const std::vector<size_t> &set,
              const std::vector<size_t> &ubddVars,
              const size_t nvars) {
    UBDD set_bdd = base.zero();
    for (auto i = std::begin(set); i != std::end(set); ++i) {
        set_bdd |= elementToBdd<UBDD>(base, *i, ubddVars, nvars);
    }
    return set_bdd;
}
