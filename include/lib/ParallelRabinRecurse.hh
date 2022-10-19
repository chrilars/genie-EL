//
// Created by rychh on 23.09.22.
//

#pragma once
#include "ubdd/SylvanUBDD.hh"
#include "BaseFixedPoint.hh"
#include "BaseRabinAutomaton.hh"

namespace fairsyn {

    TASK_DECL_4(SylvanUBDD, RabinRecurseInit, BaseFixedPoint<SylvanUBDD> *, SylvanUBDD *, struct const_arg_recursive_rabin<SylvanUBDD> *,
                struct nconst_arg_recursive_rabin<SylvanUBDD> *)

#define RabinRecurseInit(fp, controller, arg_const, arg_nconst) CALL(RabinRecurseInit, (fp), (controller), (arg_const), (arg_nconst))

        SylvanUBDD ParallelRabinRecurse(BaseFixedPoint<SylvanUBDD> *rabin,
                                        SylvanUBDD controller,
                                        const_arg_recursive_rabin<SylvanUBDD> rrConst,
                                        nconst_arg_recursive_rabin<SylvanUBDD> rrVars) {

            return RUN(RabinRecurseInit, rabin, &controller, &rrConst, &rrVars);
        }

    TASK_DECL_5(SylvanUBDD, RabinRecurseForLoop, size_t, BaseFixedPoint<SylvanUBDD> *, SylvanUBDD *,
                struct const_arg_recursive_rabin<SylvanUBDD> *, struct nconst_arg_recursive_rabin<SylvanUBDD> *);
#define RabinRecurseForLoop(i, fp, controller, arg_const, arg_nconst) (CALL(RabinRecurseForLoop, (i), (fp), (controller), (arg_const), (arg_nconst)))

    TASK_IMPL_4(SylvanUBDD,
                RabinRecurseInit,
                BaseFixedPoint<SylvanUBDD> *, fp,
                SylvanUBDD *, controller,
                const_arg_recursive_rabin<SylvanUBDD> *, arg_const,
                nconst_arg_recursive_rabin<SylvanUBDD> *, arg_nconst) {
        /* initialize a vector for storing the winning domain computed in each thread */

        std::vector<SylvanUBDD> Y;
        std::vector<SylvanUBDD *> C;
        int dupa = arg_const->pairs.size();
        for (size_t i = 0; i < arg_const->pairs.size(); i++) {
            Y.push_back(fp->base_.zero());
            SylvanUBDD controller_copy = *controller;
            C.push_back(&controller_copy);
        }

        /* spawn as many additional worker threads as the number of remaining pairs minus 1 */
        for (size_t i = 1; i < arg_const->pairs.size(); i++)
            sylvan::bdd_refs_spawn(SPAWN(RabinRecurseForLoop, i, fp, C[i], arg_const, arg_nconst));

        /* the current worker thread is repsonsible for the first rabin pair */
        SylvanUBDD y = CALL(RabinRecurseForLoop, 0, fp, C[0], arg_const, arg_nconst);
        Y[0] = y;
        sylvan::bdd_refs_push(y.bdd_.GetBDD());

        /* synchronize all the additional threads */
        for (size_t i = arg_const->pairs.size() - 1; i >= 1; i--)
            Y[i] = SylvanUBDD(sylvan::BDD(sylvan::bdd_refs_sync(SYNC(RabinRecurseForLoop).bdd_.GetBDD())));// todo check if this work properly

        /* dereference the refs pushed during the call operation */
        sylvan::bdd_refs_pop(1);

        /* compute the union of all the winning domains and the controllers computed by the different worker threads */
        SylvanUBDD U = fp->base_.zero();
        for (size_t i = 0; i < arg_const->pairs.size(); i++) {
            U |= Y[i];
            SylvanUBDD N = *C[i] & (!(controller->existAbstract(fp->CubeNotState())));
            *controller |= N;
        }

        /* return the underlying MTBDD of U */
        return U;
    }

    ///* Sequential rabin fixpoint */

    TASK_IMPL_5(SylvanUBDD,
                RabinRecurseForLoop,
                const size_t, i,
                BaseFixedPoint<SylvanUBDD> *, fp,
                SylvanUBDD *, controller,
                struct const_arg_recursive_rabin<SylvanUBDD> *, arg_const,
                struct nconst_arg_recursive_rabin<SylvanUBDD> *, arg_nconst) {
        /* unpack the inputs */
        const bool accl_on = arg_const->accl_on;
        const size_t M = arg_const->M; /* the bound on the iteration count for memorizing the BDDs from the past iterations */
        const int depth = arg_const->depth;
        const std::vector<rabin_pair_<SylvanUBDD>> pairs = arg_const->pairs;
        const SylvanUBDD initial_seed = arg_const->initial_seed;
        SylvanUBDD seqR = arg_nconst->seqR; // & fp->tr_; // todo diff
        SylvanUBDD right = arg_nconst->right; // & fp->tr_; // todo diff
        /* the original scheme from piterman pnueli paper */
        auto hist_Y = arg_nconst->hist_Y;
        auto hist_X = arg_nconst->hist_X;
        /* the scheme used in the Mascot-SDS paper */

        /* create variables for remembering the current indices of the fixpoint variables and the indices of the rabin pairs */
        std::vector<size_t> *indexY = new std::vector<size_t>;
        *indexY = *(arg_nconst->indexY);
        std::vector<size_t> *indexX = new std::vector<size_t>;
        *indexX = *(arg_nconst->indexX);
        std::vector<size_t> *indexRP = new std::vector<size_t>;
        *indexRP = *(arg_nconst->indexRP);

        if (accl_on)
            indexRP->push_back(pairs[i].rabin_index_);

        SylvanUBDD G = pairs[i].G_;
        SylvanUBDD nR = pairs[i].nR_;
        std::vector<rabin_pair_<SylvanUBDD>> remPairs = pairs;
        remPairs.erase(remPairs.begin() + i);

        const int verbose = arg_const->verbose;
//        if (verbose >= 2) { // todo diff
//            fp->printTabs(3 * depth - 1);
//            std::cout << "\n Current rabin pair: {";
//            fp->print_bdd_info(G, fp->preVars_);
//            std::cout << "}, {";
//            fp->print_bdd_info(!nR & fp->nodes_, fp->preVars_);
//            std::cout << "}\n";
//        }
        if (verbose >= 2) {
            fp->printTabs(3 * depth - 1);
            std::cout << "Remaining pairs " << pairs.size() << "\n\n";
        }
        /* initialize a local copy for the controller */
        SylvanUBDD C = fp->base_.zero();

        /////////////////////////  INITIALISING BDDS in YY

        /* initialize the sets for the nu fixed point */
        SylvanUBDD Y = fp->base_.zero();
        SylvanUBDD YY;
        if (accl_on && fp->check_threshold(*indexY,M-1) && fp->check_threshold(*indexX,M-1)) {
            YY = (*hist_Y)[depth - 1]
            [fp->findRank(fp->RabinPairs_.size(), *indexRP)]
            [std::min((*indexY)[0], M - 1)]
            [fp->to_dec(M, *indexX)];
        } else {
            YY = initial_seed;
        }

        for (int j = 0;
             Y.existAbstract(fp->CubeNotState()) != YY.existAbstract(fp->CubeNotState()); j++) { // todo diff
            Y = YY;
            if (accl_on)
                indexY->push_back(j);
            fp->print_rabin_info(Y, "Y", verbose, j, depth);
            SylvanUBDD term1 = (right | (seqR & (nR & (G & fp->cpre(Y)))));

            /* reset the local copy of the controller to the most recently added state-input pairs */
            SylvanUBDD N = term1 & (!(controller->existAbstract(fp->CubeNotState()))); // todo diff
            C = *controller | N;

            /////////////////////////  INITIALISING BDDS in XX

            /* initialize the sets for the mu fixed point */
            SylvanUBDD X = fp->base_.one();
            SylvanUBDD XX;

            if (accl_on && fp->check_threshold(*indexY, M - 1) && fp->check_threshold(*indexX, M - 1)) {
                XX = (*hist_X)[depth - 1]
                [fp->findRank(fp->RabinPairs_.size(), *indexRP)]
                [std::min((*indexX)[0], M - 1)]
                [fp->to_dec(M, *indexY)];
            } else {
                XX = fp->base_.zero();
            }

            for (int k = 0; X.existAbstract(fp->CubeNotState()) !=
                            XX.existAbstract(fp->CubeNotState()); k++) { // todo diff
                X = XX;
                if (accl_on)
                    indexX->push_back(k);
                fp->print_rabin_info(X, "X", verbose, k, depth);
                SylvanUBDD term2 = (term1 | (seqR & (nR & fp->apre(Y, X))));
                /* add the recently added state-input pairs to the controller */
                N = term2 & (!(C.existAbstract(fp->CubeNotState()))); // todo diff
                C |= N;
                if (remPairs.empty()) {
                    XX = term2;
                } else {
                    const_arg_recursive_rabin<SylvanUBDD> arg_const_new = {
                            accl_on,
                            M, /* the bound on the iteration count for memorizing the BDDs from the past iterations */
                            depth + 1,
                            remPairs,
                            initial_seed,
                            verbose};
                    nconst_arg_recursive_rabin<SylvanUBDD> arg_nconst_new = {
                            (seqR & nR), //.existAbstract(fp->CubeNotState()), // todo diff
                            term2, //.existAbstract(fp->CubeNotState()), // todo diff
                            indexRP,
                            indexY,
                            indexX,
                            hist_Y,
                            hist_X};
                    XX = CALL(RabinRecurseInit, fp, &C, &arg_const_new, &arg_nconst_new);
                }
                if (accl_on)
                    indexX->pop_back();
            }

            YY = XX;
            /////////////////////////  UPDATING BDDS in XX

            if (accl_on) { // todo diff * fp->cubeOther_) * fp->tr_
                if (accl_on && fp->check_threshold(*indexY, M - 1) && fp->check_threshold(*indexX, M - 1)) {
                    if ((*hist_X)[depth - 1]
                        [fp->findRank(fp->RabinPairs_.size(), *indexRP)]
                        [std::min((*indexX)[0] + 1, M - 1)]
                        // [fp->to_dec(M, *indexY)] <= (XX.existAbstract(fp->CubeNotState()) * fp->tr_)) {
                        [fp->to_dec(M, *indexY)] <= (XX.existAbstract(fp->cubePost_))) {
                        (*hist_X)[depth - 1]
                        [fp->findRank(fp->RabinPairs_.size(), *indexRP)]
                        [std::min((*indexX)[0] + 1, M - 1)]
                        // [fp->to_dec(M, *indexY)] = (XX.existAbstract(fp->CubeNotState()) * fp->tr_);
                        [fp->to_dec(M, *indexY)] = (XX.existAbstract(fp->cubePost_));
                    }
                }
                indexY->pop_back();
            }
        }

        *controller = C;

        /////////////////////////  UPDATING BDDS in YY
        if (accl_on) { // todo diff * fp->cubeOther_) * fp->tr_
            if (accl_on && fp->check_threshold(*indexY, M - 1) && fp->check_threshold(*indexX, M - 1)) {

                // if ((YY.existAbstract(fp->CubeNotState()) * fp->tr_) <= (*hist_Y)[depth - 1]
                if ((YY.existAbstract(fp->cubePost_)) <= (*hist_Y)[depth - 1]
                [fp->findRank(fp->RabinPairs_.size(), *indexRP)]
                [std::min((*indexY)[0] + 1, M - 1)]
                [fp->to_dec(M, *indexX)]) {
                    (*hist_Y)[depth - 1]
                    [fp->findRank(fp->RabinPairs_.size(), *indexRP)]
                    [std::min((*indexY)[0] + 1, M - 1)]
                    // [fp->to_dec(M, *indexX)] = (YY.existAbstract(fp->CubeNotState()) * fp->tr_);
                    [fp->to_dec(M, *indexX)] = (YY.existAbstract(fp->cubePost_));
                }
            }
        }


        if (accl_on)
            indexRP->pop_back();
        // YY = YY.existAbstract(fp->CubeNotState()); // todo diff

        return YY;
    }

}
