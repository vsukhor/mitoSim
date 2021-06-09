/* =============================================================================
   Copyright (C) 2015 Valerii Sukhorukov & Michael Meyer-Hermann,
   Helmholtz Center for Infection Research (Braunschweig, Germany).
   All Rights Reserved.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

================================================================================
*/

/**
* @file core_transformer.h
* @brief Contains low-level network transformations.
* @details Low-level network transformations are implemented in
*          class CoreTransformer.
* @author Valerii Sukhorukov
*/

#ifndef MITOSIM_CORE_TRANSFORMER_H
#define MITOSIM_CORE_TRANSFORMER_H

#include <vector>

#include "utils/common/constants.h"
#include "utils/common/misc.h"
#include "utils/common/msgr.h"

#include "structure.h"

namespace mitosim {

/**
 * @brief The CoreTransformer class.
 * @details Encapsulates a elementary dynamics and updates the network for it.
 * Forms base for clases adding more specific tapes of dynamics.
 * @tparam Mt type of the network Edge.
 */
template<typename Mt>
class CoreTransformer
    : public Structure<Mt> {

protected:

    using Structure<Mt>::msgr;
    using Structure<Mt>::mt;
    using Structure<Mt>::mtnum;
    using Structure<Mt>::clnum;

public:

    using Msgr = utils::common::Msgr;
    using szt = utils::common::szt;
    using ulong = utils::common::ulong;

    /**
     * @brief Constructor,
     * @param msgr Output message processor.
     */
    explicit CoreTransformer(Msgr& msgr);

    /**
     * @brief Update the structure of diaconnected components produced by fusion.
     * @param w1 Index of the 1st participant segment.
     * @param w2 Index of the 2nd participant segment.
     */
    constexpr void update_mtcl_fuse(szt w1, szt w2) noexcept;

    /**
     * @brief Update the structure of diaconnected components produced by fusion.
     * @param c1 Index of the 1st participant component.
     * @param c2 Index of the 2nd participant component.
     */
    constexpr void update_cl_fuse(szt c1, szt c2) noexcept;

    /**
     * @brief Update the diaconnected component indexes.
     * @param cf Initial index.
     * @param ct Final index.
     */
    constexpr void update_cl(szt cf, szt ct) noexcept;

    /**
     * @brief Update the diaconnected component indexes.
     * @param c Initial index.
     */
    constexpr void update_gIndcl(szt c) noexcept;

protected:

    /**
     * @brief Perform fusion of two antiparallel oriented segments.
     * @details Antiparallel orientation is defined as
     * fusion of end 1 of the 1st segment to to end 1 of the 2nd segment, or
     * fusion of end 2 of the 1st segment to to end 2 of the 2nd segment.
     * @param end End index at the fusion position.
     * @param w1 Index of the 1st segment.
     * @param w2 Index of the 2nd segment.
     */
    auto fuse_antiparallel(szt end, szt w1, szt w2) noexcept -> std::array<szt,2>;

    /**
     * @brief Perform fusion of two parallel oriented segments.
     * @details Parallel orientation is defined as
     * fusion of end 1 of the 1st segment to to end 2 of the 2nd segment, or
     * fusion of end 2 of the 1st segment to to end 1 of the 2nd segment.
     * @param w1 Index of the 1st segment.
     * @param w2 Index of the 2nd segment.
     */
    auto fuse_parallel(szt w1, szt w2) noexcept -> std::array<szt,2>;

    /**
     * @brief Update network segment indexes.
     * @details The network indexes are updatted such that segment
     *          indexed as source will become target.
     * @param f Source segment index.
     * @param t Target segment index.
    */
    void rename_mito(szt f, szt t);

    /**
     * @brief Copy connection partners to a new segment.
     * @param f Source segment index.
     * @param ef Source segment end.
     * @param t Target segment index.
     * @param et Target segment end.
    */
    void copy_neigs(szt f, szt ef, szt t, szt et) noexcept;

    /**
     * @brief Update segment neighbours.
     * @param oldn Old neighbour segment index.
     * @param oend Old neighbour segment end.
     * @param n1 Neighbour initial index.
     * @param n2 Neighbour final index.
     * @param newn New neighbour segment index.
     * @param nend New neighbour segment end.
     * @param removefromneigs Flag if the neighbouring segments should be removed.
    */
    void update_neigs(szt oldn, szt oend,
                      szt n1, szt n2,
                      szt newn, szt nend,
                      bool removefromneigs
    ) noexcept;

};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<typename Mt>
CoreTransformer<Mt>::
CoreTransformer(
        Msgr& msgr
    )
    : Structure<Mt> {msgr}
{}


template<typename Mt>
void CoreTransformer<Mt>::
rename_mito( const szt f, const szt t )
{
    copy_neigs(f, 1, t, 1);
    copy_neigs(f, 2, t, 2);
    mt[t].g = std::move(mt[f].g);
    mt[t].set_cl(mt[f].get_cl());
}


template<typename Mt>
void CoreTransformer<Mt>::
copy_neigs(
    const szt f, const szt ef,
    const szt t, const szt et
) noexcept
{
    for (szt j=1; j<=mt[f].nn[ef]; j++) {
        mt[t].neig[et][j] = mt[f].neig[ef][j];
        mt[t].neen[et][j] = mt[f].neen[ef][j];
    }
    mt[t].nn[et] = mt[f].nn[ef];

    // Substitute f in f's neig's neigs for t:
    update_neigs(f, ef, 1, mt[f].nn[ef], t, et, false);
}


template<typename Mt>
void CoreTransformer<Mt>::
update_neigs( const szt oldn, const szt oend,
              const szt n1, const szt n2,
              const szt newn, const szt nend,
              const bool removefromneigs ) noexcept
{
    for (szt j=n1; j<=n2; j++) {
        const auto cn = mt[oldn].neig[oend][j];  // our neig currently processed
        const auto ce = mt[oldn].neen[oend][j];  // our neig currently processed

        szt i1 {};
        while (true) {
            ++i1;
            XASSERT(i1 <= mt[cn].nn[ce],
                    "A neig was not found: "+
                    std::to_string(oldn)+" "+
                    std::to_string(cn)+"\n");
            if (mt[cn].neig[ce][i1] == oldn &&
                mt[cn].neen[ce][i1] == oend)
                break;
        }
        if (removefromneigs) {
            // Remove oldn from the neig list of its j-th neig
            mt[cn].neig[ce][i1] = mt[cn].neig[ce][mt[cn].nn[ce]];
            mt[cn].neen[ce][i1] = mt[cn].neen[ce][mt[cn].nn[ce]--];
            // Remove the j-th neig from the oldn's list of neigs
            mt[oldn].neig[oend][j] = mt[oldn].neig[oend][mt[oldn].nn[oend]];
            mt[oldn].neen[oend][j] = mt[oldn].neen[oend][mt[oldn].nn[oend]--];
        }
        else {
            mt[cn].neig[ce][i1] = newn;
            mt[cn].neen[ce][i1] = nend;
        }
    }
}


template<typename Mt>
auto CoreTransformer<Mt>::
fuse_antiparallel(
    const szt end,
    const szt w1,
    const szt w2
) noexcept -> std::array<szt,2>
{
    [[maybe_unused]] const auto len1 = mt[w1].g.size();
    [[maybe_unused]] const auto len2 = mt[w2].g.size();
    const auto cl1 = mt[w1].get_cl();
    const auto cl2 = mt[w2].get_cl();

    if constexpr (verbose) {
        using utils::common::STR;
        msgr.print("Fusion11a:  ",
                   w1, "(of ", len1, ") with ",
                   w2, "(of ", len2, ") at end ", end);
        mt[w1].print(w1, "     before a: ");
        mt[w2].print(w2, "     before a: ");
    }
    XASSERT(w1 != w2,
            "Error during antiparallel fusion: w1 == w2: fuse_toLoop should be used instead.\n");
    XASSERT(!mt[w1].nn[end],
            "Error during antiparallel fusion: end of w1 is not free.\n");
    XASSERT(!mt[w2].nn[end],
            "Error during antiparallel fusion: end of w2 is not free.\n");

    const szt opend = (end==2) ? 1 : 2;
    if (end == 1)
        copy_neigs(w1, 2, w1, 1);    // copy w1's 1-end neigs to its 0-end
    copy_neigs(w2, opend, w1, 2);    // copy w2's 1-end neigs to w1's 1-end

    if (mt[w2].get_cl() != mt[w1].get_cl())
        update_mtcl_fuse(w1, w2);

    if (end == 1)
        // Edge::a takes values [1:g.size()];
        // for w1 reflect positions if 1-ends are joined;
        mt[w1].reflect_g();
    else
        // Edge::a takes values [1:g.size()];
        // for w2 reflect positions if 2-ends are joined;
        mt[w2].reflect_g();

    std::move(mt[w2].g.begin(), mt[w2].g.end(), std::back_inserter(mt[w1].g));
    mt[w2].g.clear();

    if (w2 != mtnum)
        rename_mito(mtnum, w2);
    mt.pop_back();
    mtnum--;

    update_gIndcl(cl1);
    if (cl1 != cl2)
        update_gIndcl(cl2);

    if constexpr (verbose) {
        if (w1 == mtnum+1) {
            mt[w2].print(w2, "       producing ");
            if (msgr.so) *msgr.so << std::endl;
            if (msgr.sl) *msgr.sl << std::endl;
        }
        else {
            mt[w1].print(w1, "       producing ");
            if (msgr.so) *msgr.so << std::endl;
            if (msgr.sl) *msgr.sl << std::endl;
        }
    }
    return {cl1, cl2};
}


template<typename Mt>
auto CoreTransformer<Mt>::
fuse_parallel(
    const szt w1,
    const szt w2
) noexcept -> std::array<szt,2>
{
    [[maybe_unused]] const auto len1 = mt[w1].g.size();
    [[maybe_unused]] const auto len2 = mt[w2].g.size();
    const auto cl1 = mt[w1].get_cl();
    const auto cl2 = mt[w2].get_cl();

    if constexpr (verbose) {
        using utils::common::STR;
        msgr.print("Fusion11p:  ",
                   w1, "(of ", len1, ") with ",
                   w2, "(of ", len2, ")");
        mt[w1].print(w1, "     before p: ");
        mt[w2].print(w2, "     before p: ");
    }
    XASSERT(w1 != w2,
            "Error during parallel fusion: w1 == w2: fuse_toLoop should be used instead.\n");
    XASSERT(!mt[w1].nn[1],
            "Error during parallel fusion: end 1 of w1 is not free.\n");
    XASSERT(!mt[w2].nn[2],
            "Error during parallel fusion: end 2 of w2 is not free.\n");

    copy_neigs(w2, 1, w1, 1);
    if (mt[w2].get_cl() != mt[w1].get_cl())
        update_mtcl_fuse(w1, w2);

    std::move(mt[w1].g.begin(), mt[w1].g.end(), std::back_inserter(mt[w2].g));
    mt[w1].g = std::move(mt[w2].g);

    if (w2 != mtnum)
        rename_mito(mtnum, w2);
    mt.pop_back();
    mtnum--;

    update_gIndcl(cl1);
    if (cl1 != cl2)
        update_gIndcl(cl2);

    if constexpr (verbose) {
        if (w1 == mtnum+1) {
            mt[w2].print(w2, "       producing ");
            if (msgr.so) *msgr.so << std::endl;
            if (msgr.sl) *msgr.sl << std::endl;
        }
        else {
            mt[w1].print(w1, "       producing ");
            if (msgr.so) *msgr.so << std::endl;
            if (msgr.sl) *msgr.sl << std::endl;
        }
    }

    return {cl1, cl2};
}


template<typename Mt> constexpr
void CoreTransformer<Mt>::
update_mtcl_fuse(
    const szt w1,
    const szt w2
) noexcept
{
    const auto w1cl = mt[w1].get_cl();
    const auto w2cl = mt[w2].get_cl();

    for (szt i=1; i<=mtnum; i++)
        if (mt[i].get_cl() == w2cl)
            mt[i].set_cl(w1cl);

    if (w2cl != clnum-1)
        for(szt i=1; i<=mtnum; i++)
            if(mt[i].get_cl() == clnum-1)
                mt[i].set_cl(w2cl);
    clnum--;
}


template<typename Mt> constexpr
void CoreTransformer<Mt>::
update_cl_fuse(
    const szt c1,
    const szt c2
) noexcept   // args by value
{
    update_cl(c2, c1);
    if (c2 != clnum - 1)
        update_cl(clnum-1, c2);
    clnum--;
}


template<typename Mt> constexpr
void CoreTransformer<Mt>::
update_cl(
    const szt cf,
    const szt ct
) noexcept        // args by value
{
    for (szt i=1; i<=mtnum; i++)
        if (mt[i].get_cl() == cf)
            mt[i].set_cl(ct);

    update_gIndcl(ct);
}


template<typename Mt> constexpr
void CoreTransformer<Mt>::
update_gIndcl( const szt cl ) noexcept
{
    szt indcl {};
    for (szt j=1; j<=mtnum; j++)
        if (mt[j].get_cl() == cl)
            indcl = mt[j].set_gCl(cl, indcl);
}

}  // namespace mitosim

#endif  // MITOSIM_CORE_TRANSFORMER_H
