def add_cassette_exon(splicegraph, new_exon, exons_pre, exons_aft)
    ### exon_pre contains the indices of preceding exons
    ### exon_aft contains the indices of successing exons
    
    splicegraph["exons"] = sp.r_[splicegraph["exons"], new_exon.T]

    splicegraph["graph"] = sp.r_[splicegraph["graph"], sp.zeros((splicegraph["graph"].shape[1],))]
    splicegraph["graph"] = sp.c_[splicegraph["graph"], sp.zeros((splicegraph["graph"].shape[0],))]

    splicegraph["graph"][exons_pre, -1] = 1
    splicegraph["graph"][exons_aft, -1] = 1
    splicegraph["graph"][-1, :] = splicegraph["graph"][:, -1].T

    if "term" in splicegraph:
        splicegrap["term"] = sp.c_[splicegrap["term"], [[0], [0]]]

    return splicegraph
