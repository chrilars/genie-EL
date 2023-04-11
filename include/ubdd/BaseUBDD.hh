/** @file BaseUBDD.hh
 *  Contains the class BaseUBDD definition
 *
 *  @date 10.09.2021
 *  @author Mateusz Rychlicki
 */

#pragma once

#include "UBDDMintermIterator.hh"
#include <cstdio>
#include <string>
#include <vector>

namespace genie {
    /**
     * @brief An interface for Universal BDD
     */
    template <typename SubUBDD>// CRTP - Curiously recurring template pattern
    class BaseUBDD {
    public:
        /**
         * @brief Returns the %UBDD representing "False"
         * @return the %UBDD with "False" value
         */
        virtual SubUBDD zero() const = 0;

        /**
         * @brief Returns the %UBDD representing "True"
         * @return the %UBDD with "True" value
         */
        virtual SubUBDD one() const = 0;

        /**
         * @brief Returns a new %UBDD variable.
         * @details The new variable has an index equal to the largest previous
         * index plus 1.
         * @return new %UBDD variable
         */
        virtual SubUBDD var() const = 0;


        /**
         * @brief Returns a %UBDD variable with this index. If doesn't exist then create it.
         * @details If nodeSize() < index then creates ubdds for all indexes between nodeSize() and index.
         * @param index of variable.
         * @return %UBDD variable with this index
         */
        virtual SubUBDD var(size_t index) const = 0;

        /**
         * @brief Returns the number of %UBDD variables in existance.
         * @return the number of %UBDD variables in existance
         */
        virtual size_t nodeSize() const = 0;

        /**
         * @brief Gets the index of this node.
         * @return index of node
         */
        virtual size_t nodeIndex() const = 0;

        /**
         * @brief Returns the %UBDD representing a cube of variables, according to the given values.
         * @details The length of string must be equal to the number of variables in the cube.
         * For every ith char in string, if it is 0, the corresponding variable will appear in its
         * negative form, otherwise it will appear in its positive form.
         * @param variables the variables that will be in the cube in their positive or negative form
         * @param string a character array describing how the variables will appear in the result
         * @return cube of variables
         */
        virtual SubUBDD cube(std::vector<SubUBDD> variables, std::vector<int> values) const = 0;

        /**
         * @brief Returns the %UBDD representing a cube of variables, according to the given values.
         * @details The length of string must be equal to the number of variables in the cube.
         * For every ith char in string, if it is 0, the corresponding variable will appear in its
         * negative form, otherwise it will appear in its positive form.
         * @param variables the variables that will be in the cube in their positive or negative form
         * @param values string a character array describing how the variables will appear in the result
         * @return cube of variables
         */
        virtual SubUBDD cube(std::vector<SubUBDD> variables, std::vector<uint8_t> values) const = 0;

        /**
         * @brief Returns the %UBDD representing a cube of variables.
         * @details We assume that "values" contains only "1"
         * @param variables the variables that will be in the cube in their positive form
         * @return cube of variables
         */
        virtual SubUBDD cube(std::vector<SubUBDD> variables) const = 0;

        /**
         * @brief Substitute all variables in the array "from" by the corresponding variables in "to".
         * @details Variable with index from[i] became a variable with index to[i]
         * @param  from
         * @param  to
         * @return Permuted UBDD
         */
        virtual SubUBDD permute(const std::vector<size_t> &from, const std::vector<size_t> &to) const = 0;

        /**
         * @brief Compute the number of satisfying variable assignments (minterms), using the given number of variables.
         * @details Works like CountMinterm from cudd.
         * @param nvars number of variables
         * @return number of satisfying variable assignment
         */
        virtual double countMinterm(size_t nvars) const = 0;

        /**
         * @brief print a minterm of UBDD
         */
        virtual void printMinterm() const = 0;

        /**
         * @brief Computes Exists cube: f
         */
        virtual SubUBDD existAbstract(const SubUBDD &cube) const = 0;

        /**
         * @brief Computes Exists cube: f \and g
         */
        virtual SubUBDD andAbstract(const SubUBDD &g, const SubUBDD &cube) const = 0;

        /**
         * @brief Generate minterm iterator
         * @details Look at UBDDMintermIterator.hh and Polymorphism
         * @param ivars
         * @return iterator for ubbds miniterms
         */
        virtual UBDDMintermIterator *generateMintermIterator(std::vector<size_t> &ivars) const = 0;

        virtual SubUBDD complement(SubUBDD &symbolicSet,
                                   std::vector<size_t> nofGridPoints,
                                   std::vector<size_t> nofBddVars,
                                   std::vector<std::vector<size_t>> indBddVars,
                                   size_t dim) = 0;

        /**
         * @brief compute smallest radius of ball that contains L*cell
         */
        virtual SubUBDD computePolytope(const size_t p,
                                        const std::vector<double> &H,
                                        const std::vector<double> &h,
                                        int type,
                                        size_t dim,
                                        const std::vector<double> &eta,
                                        const std::vector<double> &z,
                                        const std::vector<double> &firstGridPoint,
                                        const std::vector<size_t> &nofBddVars,
                                        const std::vector<std::vector<size_t>> &indBddVars,
                                        const std::vector<size_t> &nofGridPoints) = 0;

        /**
         * @brief Save %UBDD to file.
         * @param file the file to which %UBDD will be written
         * @return 1 if success otherwise 0
         */
        virtual int save(FILE *file) = 0;

        /**
         * @brief Load %UBDD to file.
         * @details
         * If (newID) load(file,{},0).permute({1,2,3,4...n}, composeids});
         * newID=0 - the %UBDD variable ids (used for the symbolic set) are taken from file
         * newID=1 - the %UBDD variable ids (used for the symbolic set) are newly generated
         *
         * Cudd has bug when newID == 0.
         *
         * @param file
         * @param composeids
         * @param newID
         * @return Loaded %UBDD
         */
        virtual SubUBDD load(FILE *file, std::vector<int> composeids, int newID) = 0;

        /**
         * @brief Convert a %UBDD from a manager to another one.
         * This is used for paralleling in cudd.
         */
        virtual SubUBDD transfer(const SubUBDD &ubdd) const = 0;

        virtual bool isCoverEqual(const SubUBDD &other) const = 0;

        /**
         * @brief returns true if UBDD class is parallel safe
         */
        virtual bool isParallelSafe() const = 0;

        virtual SubUBDD &operator=(const SubUBDD &right) = 0;
        virtual bool operator==(const SubUBDD &other) const = 0;
        virtual bool operator!=(const SubUBDD &other) const = 0;
        virtual bool operator<=(const SubUBDD &other) const = 0;
        virtual bool operator>=(const SubUBDD &other) const = 0;
        virtual bool operator<(const SubUBDD &other) const = 0;
        virtual bool operator>(const SubUBDD &other) const = 0;

        /** @brief Not*/
        virtual SubUBDD operator!() const = 0;

        /** @brief Not*/
        virtual SubUBDD operator~() const = 0;

        /** @brief And*/
        virtual SubUBDD operator*(const SubUBDD &other) const = 0;
        virtual SubUBDD operator*=(const SubUBDD &other) = 0;

        /** @brief And*/
        virtual SubUBDD operator&(const SubUBDD &other) const = 0;
        virtual SubUBDD operator&=(const SubUBDD &other) = 0;

        /** @brief Or*/
        virtual SubUBDD operator+(const SubUBDD &other) const = 0;
        virtual SubUBDD operator+=(const SubUBDD &other) = 0;

        /** @brief Or*/
        virtual SubUBDD operator|(const SubUBDD &other) const = 0;
        virtual SubUBDD operator|=(const SubUBDD &other) = 0;

        /** @brief Xor*/
        virtual SubUBDD operator^(const SubUBDD &other) const = 0;
        virtual SubUBDD operator^=(const SubUBDD &other) = 0;

        /** @brief A and not B */
        virtual SubUBDD operator-(const SubUBDD &other) const = 0;
        virtual SubUBDD operator-=(const SubUBDD &other) = 0;
    };
}// namespace genie