//
// Created by mrychlicki on 8/12/21.
//

#ifndef MASCOTSDS_SYLVANTESTHELPER_H
#define MASCOTSDS_SYLVANTESTHELPER_H
#include "ubdd/SylvanUBDD.hh"

using namespace fairsyn;

struct SylvanTestHelper {
    SylvanTestHelper() {
        init_sylvan();
    }

    ~SylvanTestHelper() {
        stop_sylvan();
    }

    void init_sylvan() {
        /* initiate the lace work-stealing framework and the sylvan parallel bdd library */
        int dqsize = 100000;
        lace_start(2, dqsize);
        // use at most 1 GB, nodes:cache ratio 2:1, initial size 1/32 of maximum
        sylvan::sylvan_set_limits(1 * 1024 * 1024 * 1024, 1, 5);
        sylvan::sylvan_init_package();
        sylvan::sylvan_init_mtbdd();
    }

    void stop_sylvan() {
        sylvan::sylvan_stats_report(stdout);
        sylvan::sylvan_quit();
        lace_stop();
        SylvanUBDD::nodes_map = std::map<size_t, SylvanUBDD *>();
        SylvanUBDD::size_ = 0;
    }
};

#endif//MASCOTSDS_SYLVANTESTHELPER_H
