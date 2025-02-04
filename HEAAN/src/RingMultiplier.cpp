/*
 * Copyright (c) by CryptoLab inc.
 * This program is licensed under a
 * Creative Commons Attribution-NonCommercial 3.0 Unported License.
 * You should have received a copy of the license along with this
 * work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
 */

#include "RingMultiplier.h"

#include <NTL/BasicThreadPool.h>
#include <NTL/tools.h>
#include <cmath>
#include <cstdlib>
#include <iterator>

using namespace std;
using namespace NTL;

namespace heaan
{

    RingMultiplier::RingMultiplier()
    {

        uint64_t primetest = (1ULL << pbnd) + 1;
        for (long i = 0; i < nprimes; ++i)
        {
            while (true)
            {
                primetest += M;
                if (primeTest(primetest))
                {
                    pVec[i] = primetest;
                    break;
                }
            }
        }

        for (long i = 0; i < nprimes; ++i)
        {
            red_ss_array[i] = _ntl_general_rem_one_struct_build(pVec[i]);
            pInvVec[i] = inv(pVec[i]);
            prVec[i] = (static_cast<unsigned __int128>(1) << kbar2) / pVec[i];
            uint64_t root = findMthRootOfUnity(M, pVec[i]); // pVec[i] 下的 M 阶单位根
            uint64_t rootinv = invMod(root, pVec[i]);		// root 在 pVec[i] 下的乘法逆元
            uint64_t NInv = invMod(N, pVec[i]);				// N 在 pVec[i] 下的乘法逆元
            mulMod(scaledNInv[i], NInv, (1ULL << 32), pVec[i]);
            mulMod(scaledNInv[i], scaledNInv[i], (1ULL << 32), pVec[i]);
            scaledRootPows[i] = new uint64_t[N]();
            scaledRootInvPows[i] = new uint64_t[N]();
            uint64_t power = 1;
            uint64_t powerInv = 1;
            for (long j = 0; j < N; ++j)
            {
                uint32_t jprime = bitReverse(static_cast<uint32_t>(j)) >> (32 - logN); // j 在 N-bits 下的位反转结果
                uint64_t rootpow = power;											   // root^j
                mulMod(scaledRootPows[i][jprime], rootpow, (1ULL << 32), pVec[i]);
                mulMod(scaledRootPows[i][jprime], scaledRootPows[i][jprime], (1ULL << 32), pVec[i]);
                uint64_t rootpowInv = powerInv;
                mulMod(scaledRootInvPows[i][jprime], rootpowInv, (1ULL << 32), pVec[i]);
                mulMod(scaledRootInvPows[i][jprime], scaledRootInvPows[i][jprime], (1ULL << 32), pVec[i]);
                mulMod(power, power, root, pVec[i]);
                mulMod(powerInv, powerInv, rootinv, pVec[i]);
            }
        }

        for (long i = 0; i < nprimes; ++i)
        {
            coeffpinv_array[i] = new mulmod_precon_t[i + 1];
            pProd[i] = (i == 0) ? to_ZZ((long)pVec[i]) : pProd[i - 1] * (long)pVec[i];
            pProdh[i] = pProd[i] / 2;
            pHat[i] = new ZZ[i + 1];
            pHatInvModp[i] = new uint64_t[i + 1];
            for (long j = 0; j < i + 1; ++j)
            {
                pHat[i][j] = ZZ(1);
                for (long k = 0; k < j; ++k)
                {
                    pHat[i][j] *= (long)pVec[k];
                }
                for (long k = j + 1; k < i + 1; ++k)
                {
                    pHat[i][j] *= (long)pVec[k];
                }
                pHatInvModp[i][j] = to_long(pHat[i][j] % (long)pVec[j]);
                pHatInvModp[i][j] = invMod(pHatInvModp[i][j], pVec[j]);
                coeffpinv_array[i][j] = PrepMulModPrecon(pHatInvModp[i][j], pVec[j]);
            }
        }
    }

    // 米勒-拉宾素性测试
    bool RingMultiplier::primeTest(uint64_t p)
    {
        if (p < 2)
            return false;
        if (p != 2 && p % 2 == 0)
            return false;
        uint64_t s = p - 1;
        while (s % 2 == 0)
        {
            s /= 2;
        }
        for (long i = 0; i < 200; i++)
        {
            uint64_t temp1 = rand();
            temp1 = (temp1 << 32) | rand();
            temp1 = temp1 % (p - 1) + 1;
            uint64_t temp2 = s;
            uint64_t mod = powMod(temp1, temp2, p);
            while (temp2 != p - 1 && mod != 1 && mod != p - 1)
            {
                mulMod(mod, mod, mod, p);
                temp2 *= 2;
            }
            if (mod != p - 1 && temp2 % 2 == 0)
                return false;
        }
        return true;
    }

    // 基于 CT 蝴蝶的 NTT
    void RingMultiplier::NTT(uint64_t *a, long index)
    {
        long t = N;
        long logt1 = logN + 1;
        uint64_t p = pVec[index];
        uint64_t pInv = pInvVec[index];
        for (long m = 1; m < N; m <<= 1)
        {
            t >>= 1;
            logt1 -= 1;
            for (long i = 0; i < m; i++)
            {
                long j1 = i << logt1;
                long j2 = j1 + t - 1;
                uint64_t W = scaledRootPows[index][m + i];
                for (long j = j1; j <= j2; j++)
                {
                    butt(a[j], a[j + t], W, p, pInv);
                }
            }
        }
    }

    // 基于 GS 蝴蝶的 NTT
    void RingMultiplier::INTT(uint64_t *a, long index)
    {
        uint64_t p = pVec[index];
        uint64_t pInv = pInvVec[index];
        long t = 1;
        for (long m = N; m > 1; m >>= 1)
        {
            long j1 = 0;
            long h = m >> 1;
            for (long i = 0; i < h; i++)
            {
                long j2 = j1 + t - 1;
                uint64_t W = scaledRootInvPows[index][h + i];
                for (long j = j1; j <= j2; j++)
                {
                    ibutt(a[j], a[j + t], W, p, pInv);
                }
                j1 += (t << 1);
            }
            t <<= 1;
        }

        uint64_t NScale = scaledNInv[index];
        for (long i = 0; i < N; i++)
        {
            idivN(a[i], NScale, p, pInv);
        }
    }

    //----------------------------------------------------------------------------------
    //   FFT
    //----------------------------------------------------------------------------------

    // 批量将多项式系数使用 CRT 拆分，然后转换为 NTT 形式。
    void RingMultiplier::CRT(uint64_t *rx, ZZ *x, const long np)
    {
        NTL_EXEC_RANGE(np, first, last);
        for (long i = first; i < last; ++i)
        {
            uint64_t *rxi = rx + (i << logN);
            uint64_t pi = pVec[i];
            _ntl_general_rem_one_struct *red_ss = red_ss_array[i];
            for (long n = 0; n < N; ++n)
            {
                // 转换为蒙哥马利形式
                rxi[n] = _ntl_general_rem_one_struct_apply(x[n].rep, pi, red_ss);
            }
            NTT(rxi, i);
        }
        NTL_EXEC_RANGE_END;
    }

    // 批量进行 NTT 形式的加法
    void RingMultiplier::addNTTAndEqual(uint64_t *ra, uint64_t *rb, const long np)
    {
        for (long i = 0; i < np; ++i)
        {
            uint64_t *rai = ra + (i << logN);
            uint64_t *rbi = rb + (i << logN);
            uint64_t pi = pVec[i];
            for (long n = 0; n < N; ++n)
            {
                rai[n] += rbi[n];
                if (rai[n] > pi)
                    rai[n] -= pi;
            }
        }
    }

    // 使用 CRT 还原
    void RingMultiplier::reconstruct(ZZ *x, uint64_t *rx, long np, const ZZ &q)
    {
        ZZ *pHatnp = pHat[np - 1];
        uint64_t *pHatInvModpnp = pHatInvModp[np - 1];
        mulmod_precon_t *coeffpinv_arraynp = coeffpinv_array[np - 1];
        ZZ &pProdnp = pProd[np - 1];
        ZZ &pProdhnp = pProdh[np - 1];
        NTL_EXEC_RANGE(N, first, last);
        for (long n = first; n < last; ++n)
        {
            ZZ &acc = x[n];
            QuickAccumBegin(acc, pProdnp.size());
            for (long i = 0; i < np; i++)
            {
                long p = pVec[i];
                long tt = pHatInvModpnp[i];
                mulmod_precon_t ttpinv = coeffpinv_arraynp[i];
                long s = MulModPrecon(rx[n + (i << logN)], tt, p, ttpinv);
                QuickAccumMulAdd(acc, pHatnp[i], s);
            }
            QuickAccumEnd(acc);
            rem(x[n], x[n], pProdnp);
            if (x[n] > pProdhnp)
                x[n] -= pProdnp;
            x[n] %= q;
        }
        NTL_EXEC_RANGE_END;
    }

    /*	输入：a(X)，多项式形式，存储于 a。
     *	输入：b(X)，多项式形式，存储于 b。
     *	输出：a(X)*b(X)，多项式形式，存储到 x。
     *	系数模除 mod。
     */
    void RingMultiplier::mult(ZZ *x, ZZ *a, ZZ *b, long np, const ZZ &mod)
    {
        uint64_t *ra = new uint64_t[np << logN]();
        uint64_t *rb = new uint64_t[np << logN]();
        uint64_t *rx = new uint64_t[np << logN]();

        NTL_EXEC_RANGE(np, first, last);
        for (long i = first; i < last; ++i)
        {
            uint64_t *rai = ra + (i << logN);
            uint64_t *rbi = rb + (i << logN);
            uint64_t *rxi = rx + (i << logN);
            uint64_t pi = pVec[i];
            uint64_t pri = prVec[i];
            _ntl_general_rem_one_struct *red_ss = red_ss_array[i];
            for (long n = 0; n < N; ++n)
            {
                rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
                rbi[n] = _ntl_general_rem_one_struct_apply(b[n].rep, pi, red_ss);
            }
            NTT(rai, i);
            NTT(rbi, i);
            for (long n = 0; n < N; ++n)
            {
                mulModBarrett(rxi[n], rai[n], rbi[n], pi, pri);
            }
            INTT(rxi, i);
        }
        NTL_EXEC_RANGE_END;

        reconstruct(x, rx, np, mod);

        delete[] ra;
        delete[] rb;
        delete[] rx;
    }

    /*	输入：a(X)，多项式形式，存储于 a。
     *	输入：b(X)，点值形式，存储于 rb。
     *	输出：a(X)*b(X)，多项式形式，存储到 x。
     *	系数模除 mod。
     */
    void RingMultiplier::multNTT(ZZ *x, ZZ *a, uint64_t *rb, long np, const ZZ &mod)
    {
        uint64_t *ra = new uint64_t[np << logN]();
        uint64_t *rx = new uint64_t[np << logN]();
        NTL_EXEC_RANGE(np, first, last);
        for (long i = first; i < last; ++i)
        {
            uint64_t *rai = ra + (i << logN);
            uint64_t *rbi = rb + (i << logN);
            uint64_t *rxi = rx + (i << logN);
            uint64_t pi = pVec[i];
            uint64_t pri = prVec[i];
            _ntl_general_rem_one_struct *red_ss = red_ss_array[i];
            for (long n = 0; n < N; ++n)
            {
                rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
            }
            NTT(rai, i);
            for (long n = 0; n < N; ++n)
            {
                mulModBarrett(rxi[n], rai[n], rbi[n], pi, pri);
            }
            INTT(rxi, i);
        }
        NTL_EXEC_RANGE_END;

        reconstruct(x, rx, np, mod);

        delete[] ra;
        delete[] rx;
    }

    /*	输入：a(X)，点值形式，存储于 ra。
     *	输入：b(X)，点值形式，存储于 rb。
     *	输出：a(X)*b(X)，多项式形式，存储到 x。
     *	系数模除 mod。
     */
    void RingMultiplier::multDNTT(ZZ *x, uint64_t *ra, uint64_t *rb, long np, const ZZ &mod)
    {
        uint64_t *rx = new uint64_t[np << logN]();

        NTL_EXEC_RANGE(np, first, last);
        for (long i = first; i < last; ++i)
        {
            uint64_t *rai = ra + (i << logN);
            uint64_t *rbi = rb + (i << logN);
            uint64_t *rxi = rx + (i << logN);
            uint64_t pi = pVec[i];
            uint64_t pri = prVec[i];
            for (long n = 0; n < N; ++n)
            {
                mulModBarrett(rxi[n], rai[n], rbi[n], pi, pri);
            }
            INTT(rxi, i);
        }
        NTL_EXEC_RANGE_END;

        reconstruct(x, rx, np, mod);

        delete[] rx;
    }

    /*	输入：a(X)，多项式形式，存储于 a。
     *	输入：b(X)，多项式形式，存储于 b。
     *	输出：a(X)*b(X)，多项式形式，存储到 a。
     *	系数模除 mod。
     */
    void RingMultiplier::multAndEqual(ZZ *a, ZZ *b, long np, const ZZ &mod)
    {
        uint64_t *ra = new uint64_t[np << logN]();
        uint64_t *rb = new uint64_t[np << logN]();

        NTL_EXEC_RANGE(np, first, last);
        for (long i = first; i < last; ++i)
        {
            uint64_t *rai = ra + (i << logN);
            uint64_t *rbi = rb + (i << logN);
            uint64_t pi = pVec[i];
            uint64_t pri = prVec[i];
            _ntl_general_rem_one_struct *red_ss = red_ss_array[i];
            for (long n = 0; n < N; ++n)
            {
                rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
                rbi[n] = _ntl_general_rem_one_struct_apply(b[n].rep, pi, red_ss);
            }
            NTT(rai, i);
            NTT(rbi, i);
            for (long n = 0; n < N; ++n)
            {
                mulModBarrett(rai[n], rai[n], rbi[n], pi, pri);
            }
            INTT(rai, i);
        }
        NTL_EXEC_RANGE_END;

        ZZ *pHatnp = pHat[np - 1];
        uint64_t *pHatInvModpnp = pHatInvModp[np - 1];

        reconstruct(a, ra, np, mod);

        delete[] ra;
        delete[] rb;
    }

    /*	输入：a(X)，多项式形式，存储于 a。
     *	输入：b(X)，点值形式，存储于 rb。
     *	输出：a(X)*b(X)，多项式形式，存储到 a。
     *	系数模除 mod。
     */
    void RingMultiplier::multNTTAndEqual(ZZ *a, uint64_t *rb, long np, const ZZ &mod)
    {
        uint64_t *ra = new uint64_t[np << logN]();

        NTL_EXEC_RANGE(np, first, last);
        for (long i = first; i < last; ++i)
        {
            uint64_t *rai = ra + (i << logN);
            uint64_t *rbi = rb + (i << logN);
            uint64_t pi = pVec[i];
            uint64_t pri = prVec[i];
            _ntl_general_rem_one_struct *red_ss = red_ss_array[i];
            for (long n = 0; n < N; ++n)
            {
                rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
            }
            NTT(rai, i);
            for (long n = 0; n < N; ++n)
            {
                mulModBarrett(rai[n], rai[n], rbi[n], pi, pri);
            }
            INTT(rai, i);
        }
        NTL_EXEC_RANGE_END;

        ZZ *pHatnp = pHat[np - 1];
        uint64_t *pHatInvModpnp = pHatInvModp[np - 1];

        reconstruct(a, ra, np, mod);

        delete[] ra;
    }

    /*	输入：a(X)，多项式形式，存储于 a。
     *	输出：a(X)^2，多项式形式，存储到 x。
     *	系数模除 mod。
     */
    void RingMultiplier::square(ZZ *x, ZZ *a, long np, const ZZ &mod)
    {
        uint64_t *ra = new uint64_t[np << logN]();
        uint64_t *rx = new uint64_t[np << logN]();

        NTL_EXEC_RANGE(np, first, last);
        for (long i = first; i < last; ++i)
        {
            uint64_t *rai = ra + (i << logN);
            uint64_t *rxi = rx + (i << logN);
            uint64_t pi = pVec[i];
            uint64_t pri = prVec[i];
            _ntl_general_rem_one_struct *red_ss = red_ss_array[i];
            for (long n = 0; n < N; ++n)
            {
                rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
            }
            NTT(rai, i);
            for (long n = 0; n < N; ++n)
            {
                mulModBarrett(rxi[n], rai[n], rai[n], pi, pri);
            }
            INTT(rxi, i);
        }
        NTL_EXEC_RANGE_END;

        ZZ *pHatnp = pHat[np - 1];
        uint64_t *pHatInvModpnp = pHatInvModp[np - 1];

        reconstruct(x, rx, np, mod);

        delete[] ra;
        delete[] rx;
    }

    /*	输入：a(X)，点值形式，存储于 ra。
     *	输出：a(X)^2，多项式形式，存储到 x。
     *	系数模除 mod。
     */
    void RingMultiplier::squareNTT(ZZ *x, uint64_t *ra, long np, const ZZ &mod)
    {
        uint64_t *rx = new uint64_t[np << logN]();

        NTL_EXEC_RANGE(np, first, last);
        for (long i = first; i < last; ++i)
        {
            uint64_t *rai = ra + (i << logN);
            uint64_t *rxi = rx + (i << logN);
            uint64_t pi = pVec[i];
            uint64_t pri = prVec[i];
            for (long n = 0; n < N; ++n)
            {
                mulModBarrett(rxi[n], rai[n], rai[n], pi, pri);
            }
            INTT(rxi, i);
        }
        NTL_EXEC_RANGE_END;

        reconstruct(x, rx, np, mod);

        delete[] rx;
    }

    /*	输入：a(X)，多项式形式，存储于 a。
     *	输出：a(X)^2，多项式形式，存储到 a。
     *	系数模除 mod。
     */
    void RingMultiplier::squareAndEqual(ZZ *a, long np, const ZZ &mod)
    {
        uint64_t *ra = new uint64_t[np << logN]();

        NTL_EXEC_RANGE(np, first, last);
        for (long i = first; i < last; ++i)
        {
            uint64_t *rai = ra + (i << logN);
            uint64_t pi = pVec[i];
            uint64_t pri = prVec[i];
            _ntl_general_rem_one_struct *red_ss = red_ss_array[i];
            for (long n = 0; n < N; ++n)
            {
                // 生成多项式 a(x) 在 pi 上的蒙哥马利形式
                rai[n] = _ntl_general_rem_one_struct_apply(a[n].rep, pi, red_ss);
            }
            NTT(rai, i);
            for (long n = 0; n < N; ++n)
            {
                // 为什么使用巴雷特模约，而不是蒙哥马利模乘？？？
                mulModBarrett(rai[n], rai[n], rai[n], pi, pri);
            }
            INTT(rai, i);
        }
        NTL_EXEC_RANGE_END;

        reconstruct(a, ra, np, mod);

        delete[] ra;
    }

    // 计算模乘 a * b % m (128-bits)
    void RingMultiplier::mulMod(uint64_t &r, uint64_t a, uint64_t b, uint64_t m)
    {
        unsigned __int128 mul = static_cast<unsigned __int128>(a) * b;
        mul %= static_cast<unsigned __int128>(m);
        r = static_cast<uint64_t>(mul);
    }

    // 基于巴雷特模乘计算 r = a * b % p
    void RingMultiplier::mulModBarrett(uint64_t &r, uint64_t a, uint64_t b, uint64_t p, uint64_t pr)
    {
        unsigned __int128 mul = static_cast<unsigned __int128>(a) * b;
        uint64_t abot = static_cast<uint64_t>(mul);
        uint64_t atop = static_cast<uint64_t>(mul >> 64);
        unsigned __int128 tmp = static_cast<unsigned __int128>(abot) * pr;
        tmp >>= 64;
        tmp += static_cast<unsigned __int128>(atop) * pr;
        tmp >>= kbar2 - 64;
        tmp *= p;
        tmp = mul - tmp;
        r = static_cast<uint64_t>(tmp);
        if (r >= p)
            r -= p;
    }

    /*	在蒙哥马利形式下计算：
     *	a = a + b * w
     *  b = a - b * w
     */
    void RingMultiplier::butt(uint64_t &a, uint64_t &b, uint64_t W, uint64_t p, uint64_t pInv)
    {
        unsigned __int128 U = static_cast<unsigned __int128>(b) * W;
        uint64_t U0 = static_cast<uint64_t>(U);
        uint64_t U1 = U >> 64;
        uint64_t Q = U0 * pInv;
        unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * p;
        uint64_t H = Hx >> 64;
        uint64_t V = U1 < H ? U1 + p - H : U1 - H;
        b = a < V ? a + p - V : a - V;
        a += V;
        if (a > p)
            a -= p;
    }

    /*	在蒙哥马利形式下计算：
     *	a = a + b
     *  b = (a - b) * w
     */
    void RingMultiplier::ibutt(uint64_t &a, uint64_t &b, uint64_t W, uint64_t p, uint64_t pInv)
    {
        uint64_t T = a < b ? a + p - b : a - b;
        a += b;
        if (a > p)
            a -= p;
        unsigned __int128 UU = static_cast<unsigned __int128>(T) * W;
        uint64_t U0 = static_cast<uint64_t>(UU);
        uint64_t U1 = UU >> 64;
        uint64_t Q = U0 * pInv;
        unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * p;
        uint64_t H = Hx >> 64;
        b = (U1 < H) ? U1 + p - H : U1 - H;
    }

    // 在蒙哥马利形式下计算 a * NScale
    void RingMultiplier::idivN(uint64_t &a, uint64_t NScale, uint64_t p, uint64_t pInv)
    {
        unsigned __int128 U = static_cast<unsigned __int128>(a) * NScale;
        uint64_t U0 = static_cast<uint64_t>(U);
        uint64_t U1 = U >> 64;
        uint64_t Q = U0 * pInv;
        unsigned __int128 Hx = static_cast<unsigned __int128>(Q) * p;
        uint64_t H = Hx >> 64;
        a = (U1 < H) ? U1 + p - H : U1 - H;
    }

    // 计算 x ^ (-1) (mod m)
    uint64_t RingMultiplier::invMod(uint64_t x, uint64_t m)
    {
        return powMod(x, m - 2, m);
    }

    // 计算 x ^ y % modulus
    uint64_t RingMultiplier::powMod(uint64_t x, uint64_t y, uint64_t modulus)
    {
        uint64_t res = 1;
        while (y > 0)
        {
            if (y & 1)
            {
                mulMod(res, res, x, modulus);
            }
            y = y >> 1;
            mulMod(x, x, x, modulus);
        }
        return res;
    }

    // 计算模 2^64 下的逆元
    uint64_t RingMultiplier::inv(uint64_t x)
    {
        return pow(x, static_cast<uint64_t>(-1));
    }

    // 计算模 2^64 下的快速幂
    uint64_t RingMultiplier::pow(uint64_t x, uint64_t y)
    {
        uint64_t res = 1;
        while (y > 0)
        {
            if (y & 1)
            {
                res *= x;
            }
            y = y >> 1;
            x *= x;
        }
        return res;
    }

    // 基于 32 位进行位反转
    uint32_t RingMultiplier::bitReverse(uint32_t x)
    {
        x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
        x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
        x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
        x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
        return ((x >> 16) | (x << 16));
    }

    // 因式分解
    void RingMultiplier::findPrimeFactors(vector<uint64_t> &s, uint64_t number)
    {
        while (number % 2 == 0)
        {
            s.push_back(2);
            number /= 2;
        }
        for (uint64_t i = 3; i < sqrt(number); i++)
        {
            while (number % i == 0)
            {
                s.push_back(i);
                number /= i;
            }
        }
        if (number > 2)
        {
            s.push_back(number);
        }
    }

    // 求原根
    uint64_t RingMultiplier::findPrimitiveRoot(uint64_t modulus)
    {
        vector<uint64_t> s;
        uint64_t phi = modulus - 1;
        findPrimeFactors(s, phi);
        for (uint64_t r = 2; r <= phi; r++)
        {
            bool flag = false;
            for (auto it = s.begin(); it != s.end(); it++)
            {
                if (powMod(r, phi / (*it), modulus) == 1)
                {
                    flag = true;
                    break;
                }
            }
            if (flag == false)
            {
                return r;
            }
        }
        return -1;
    }

    // 求 M 阶单位根
    uint64_t RingMultiplier::findMthRootOfUnity(uint64_t M, uint64_t mod)
    {
        uint64_t res;
        res = findPrimitiveRoot(mod);
        if ((mod - 1) % M == 0)
        {
            uint64_t factor = (mod - 1) / M;
            res = powMod(res, factor, mod);
            return res;
        }
        else
        {
            return -1;
        }
    }

} // namespace heaan
