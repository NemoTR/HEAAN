/*
 * Copyright (c) by CryptoLab inc.
 * This program is licensed under a
 * Creative Commons Attribution-NonCommercial 3.0 Unported License.
 * You should have received a copy of the license along with this
 * work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
 */
#ifndef HEAAN_RINGMULTIPLIER_H_
#define HEAAN_RINGMULTIPLIER_H_

#include <cstdint>
#include <vector>
#include <NTL/ZZ.h>
#include "Params.h"

namespace heaan
{

	class RingMultiplier
	{
	public:
		uint64_t *pVec = new uint64_t[nprimes];	   // 比 2^pbnd 稍大的一些随机素数
		uint64_t *prVec = new uint64_t[nprimes];   // 用于巴雷特模约的参数
		uint64_t *pInvVec = new uint64_t[nprimes]; // 在 2^64 下的乘法逆元
		uint64_t **scaledRootPows = new uint64_t *[nprimes]; // scaledRootPows[i][j] = Mthroot^j << 64 (mod pVec[i])
		uint64_t **scaledRootInvPows = new uint64_t *[nprimes]; // scaledRootInvPows[i][j] = MthrootInv[i]^j << 64 (mod pVec[i])
		uint64_t *scaledNInv = new uint64_t[nprimes]; // scaledNInv[i] = Ninv << 64 % pVec[i]
		_ntl_general_rem_one_struct *red_ss_array[nprimes]; // red_ss_array[i] = _ntl_general_rem_one_struct_build(pVec[i]);
		NTL::mulmod_precon_t *coeffpinv_array[nprimes];

		NTL::ZZ *pProd = new NTL::ZZ[nprimes]; // pProd[i] = pVec[0] * ... * pVec[i]
		NTL::ZZ *pProdh = new NTL::ZZ[nprimes]; // pProdh[i] = pProd[i] / 2
		NTL::ZZ **pHat = new NTL::ZZ *[nprimes]; // pHat[i][j] = pVec[0] * ... * pVec[j - 1] * pVec[j + 1] * ... * pVec[i]
		uint64_t **pHatInvModp = new uint64_t *[nprimes]; // pHatInvModp[i][j] = pHat[i][j] ^ (-1) (mod pVec[j])

		RingMultiplier();

		bool primeTest(uint64_t p);

		void NTT(uint64_t *a, long index);
		void INTT(uint64_t *a, long index);

		void CRT(uint64_t *rx, NTL::ZZ *x, const long np);

		void addNTTAndEqual(uint64_t *ra, uint64_t *rb, const long np);

		void reconstruct(NTL::ZZ *x, uint64_t *rx, long np, const NTL::ZZ &QQ);

		void mult(NTL::ZZ *x, NTL::ZZ *a, NTL::ZZ *b, long np, const NTL::ZZ &QQ);

		void multNTT(NTL::ZZ *x, NTL::ZZ *a, uint64_t *rb, long np, const NTL::ZZ &QQ);

		void multDNTT(NTL::ZZ *x, uint64_t *ra, uint64_t *rb, long np, const NTL::ZZ &QQ);

		void multAndEqual(NTL::ZZ *a, NTL::ZZ *b, long np, const NTL::ZZ &QQ);

		void multNTTAndEqual(NTL::ZZ *a, uint64_t *rb, long np, const NTL::ZZ &QQ);

		void square(NTL::ZZ *x, NTL::ZZ *a, long np, const NTL::ZZ &QQ);

		void squareNTT(NTL::ZZ *x, uint64_t *ra, long np, const NTL::ZZ &QQ);

		void squareAndEqual(NTL::ZZ *a, long np, const NTL::ZZ &QQ);

		void mulMod(uint64_t &r, uint64_t a, uint64_t b, uint64_t p);

		void mulModBarrett(uint64_t &r, uint64_t a, uint64_t b, uint64_t p, uint64_t pr);
		void butt(uint64_t &a, uint64_t &b, uint64_t W, uint64_t p, uint64_t pInv);
		void ibutt(uint64_t &a, uint64_t &b, uint64_t W, uint64_t p, uint64_t pInv);
		void idivN(uint64_t &a, uint64_t NScale, uint64_t p, uint64_t pInv);

		uint64_t invMod(uint64_t x, uint64_t p);

		uint64_t powMod(uint64_t x, uint64_t y, uint64_t p);

		uint64_t inv(uint64_t x);

		uint64_t pow(uint64_t x, uint64_t y);

		uint32_t bitReverse(uint32_t x);

		void findPrimeFactors(std::vector<uint64_t> &s, uint64_t number);

		uint64_t findPrimitiveRoot(uint64_t m);

		uint64_t findMthRootOfUnity(uint64_t M, uint64_t p);
	};

} // namespace heaan

#endif /* RINGMULTIPLIER_H_ */
