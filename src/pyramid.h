#pragma once

#include <vector>
#include <cstring>
#include <immintrin.h>
#include "threadpool.h"

class Pyramid {
public:
	class Level {
	public:
		int width, height;
		int pitch;
		int m128_pitch;
		size_t bytes;
		float* data;
		int x;
		int y;
		bool x_shift;
		bool y_shift;
		bool upper_x_shift;
		int upper_m128_pitch;
		std::vector<int> bands;
		__m128* GetLine(int y);
		void GetExpandedLine(__m128* temp, int y);
	};

private:
	std::vector<Level> levels;
	__m128** lines;
	bool shared;
	bool no_alloc;
	float* lut = NULL;
	int lut_bits = 0;
	bool lut_gamma = false;
	int out_max;
	Threadpool* threadpool;
	void set_lut(int bits, bool gamma);
	void CopyInterleavedThread_8bit(uint8_t* src_p, int step, int pitch, int sy, int ey);
	void CopyInterleavedThread_16bit(uint16_t* src_p, int step, int pitch, int sy, int ey);
	void CopyPlanarThread_8bit(uint8_t* src_p, int pitch, bool gamma, int sy, int ey);
	void CopyPlanarThread_16bit(uint16_t* src_p, int pitch, bool gamma, int sy, int ey);
	void CopyPlanarThread_32bit(__m128* src_p, int pitch, bool gamma, int sy, int ey);
	void Subsample_Squeeze(__m128* in, __m128* Out, int m128_pitch_in, int m128_pitch_out, __m128* mul);
	void ShrinkThread(__m128* line, __m128* hi, __m128* lo, int m128_pitch_hi, int m128_pitch_lo, int first_bad_line, int height_odd, int sy, int ey, const bool x_shift, const bool y_shift);

	void Squeeze(__m128* line, __m128* lo, int m128_pitch_lo, int m128_pitch_hi, __m128 final_mul, bool x_shift); // was __forceinline
	void LaplaceThreadWrapper(Level* upper_level, Level* lower_level, int sy, int ey);
	void LaplaceThread(Level* upper_level, Level* lower_evel, int sy, int ey, __m128* temp1, __m128* temp2, __m128* temp3);
	void FuseThread(__m128* a, __m128* b, __m128* m, int m128_pitch, int sy, int ey, bool pre, int black);
	void LaplaceCollapse(int n_levels, bool Collapse);
	size_t total_bytes = 0;

#define PLL(X) [](float* p, __m128 d, __m128 z, __m128 m) { return X; }
#define   N _mm_load_ps(p)
#define   G _mm_sqrt_ps(_mm_load_ps(p))
#define   D _mm_add_ps(_mm_load_ps(p), d)
#define   C _mm_min_ps(_mm_max_ps(_mm_load_ps(p), z), m)
#define  DG _mm_add_ps(_mm_sqrt_ps(_mm_load_ps(p)), d)
#define  CG _mm_min_ps(_mm_sqrt_ps(_mm_max_ps(_mm_load_ps(p), z)), m)
#define  CD _mm_min_ps(_mm_add_ps(_mm_max_ps(_mm_load_ps(p), z), d), m)
#define CDG _mm_min_ps(_mm_add_ps(_mm_sqrt_ps(_mm_max_ps(_mm_load_ps(p), z)), d), m)

#define _OP_ (void* _dst_p, F _loader, int pitch, int sy, int ey, int level, bool chroma)

	template <typename F>
	void OutPlanar8  _OP_;
	template <typename F>
	void OutPlanar16 _OP_;
	template <typename F>
	void OutPlanar32 _OP_;

	template <typename T, typename F>
	void OutInterleaved(T dst_p, F _loader, int pitch, int sy, int ey, int level, bool chroma, int step, int offset);

public:
	Pyramid(int width, int height, int _levels = 0, Pyramid* share = NULL);
	Pyramid(int width, int height, int _levels, int x, int y, bool no_alloc, Pyramid* share = NULL);
	~Pyramid();
	static int DefaultNumLevels(int width, int height) { return 8; /* (int)ceil(log2(max(width, height))); */ };
	void Copy(uint8_t* src_p, int step, int pitch, bool gamma, int bits);
	void Subsample(int sub_w, int sub_h, Pyramid* source);
	void Shrink();
	void Laplace() { LaplaceCollapse((int)levels.size(), false); };
	void Collapse() { LaplaceCollapse((int)levels.size(), true); };
	void Collapse(int n_levels) { LaplaceCollapse(n_levels, true); };
	float Average();
	void Add(float add, int _levels = 0x7fff);
	void MultiplyAndAdd(float add, float mul, int _levels = 0x7fff);
	void MultiplyAddClamp(float add, float mul, int level = 0x7fff);
	void Multiply(int level, float mul);
	void MultplyByPyramid(Pyramid* b);
	void Fuse(Pyramid* b, Pyramid* mask, bool pre, int black);
	void Fuse(Pyramid* b, float weight);
	void Denoise(int level, float power, bool gamma);
	void Blend(Pyramid* b);
	void BlurX(float radius, Pyramid* transpose);
	void BlurXThread(float radius, Pyramid* transpose, int sy, int ey);
	template <typename T>
	void Out(T dst_p, int pitch, bool gamma, bool dither, bool clamp, int level = 0, int step = 0, int offset = 0, bool chroma = false);
	int GetNLevels() { return (int)levels.size(); };
	int GetPitch(int level = 0) { return levels[level].pitch; };
	float* GetData(int level = 0) { return levels[level].data; };
	int GetWidth(int level = 0) { return levels[level].width; };
	int GetHeight(int level = 0) { return levels[level].height; };
	int GetX(int level = 0) { return levels[level].x; };
	int GetY(int level = 0) { return levels[level].y; };
	static void LaplaceLine2(__m128* hi, __m128* temp1, __m128* temp2, int m128_pitch);
	static void LaplaceLine3(__m128* hi, __m128* temp1, __m128* temp2, __m128* temp3, int m128_pitch);
	static void LaplaceExpand(__m128* hi, __m128* lo, int m128_pitch_hi, int m128_pitch_lo);
	static void LaplaceExpandShifted(__m128* hi, __m128* lo, int m128_pitch_hi, int m128_pitch_lo);
	size_t GetTotalBytes() { return total_bytes; }
	std::vector<Level>& GetLevels() { return levels; };
	Level& GetLevel(int level) { return levels[level]; };
#ifdef PNGER
	void Png(const char* filename);
#endif
#ifdef SIMPLE
	void shrink_simple();
	void laplace_simple(bool collapse);
#endif
};

/***********************************************************************
************************************************************************
* Output
************************************************************************
***********************************************************************/
#define LOAD _loader((float*)p_p++, dither_add, zeroes, maxes)

/***********************************************************************
* 8-bit planar
***********************************************************************/
template <typename F>
void Pyramid::OutPlanar8 _OP_{
	int x, y;
	uint8_t* dst_p = (uint8_t*)_dst_p;
	uint8_t black = chroma ? 0x80 : 0x00;

	__m128 zeroes = _mm_setzero_ps();
	__m128 maxes = _mm_set_ps1(255.0f);

	__m128i shuffle1 = _mm_set_epi32(0x80808080, 0x80808080, 0x80808080, 0x0c080400);
	__m128i shuffle2 = _mm_set_epi32(0x80808080, 0x80808080, 0x0c080400, 0x80808080);
	__m128i shuffle3 = _mm_set_epi32(0x80808080, 0x0c080400, 0x80808080, 0x80808080);
	__m128i shuffle4 = _mm_set_epi32(0x0c080400, 0x80808080, 0x80808080, 0x80808080);
	__m128i four_shuffle = _mm_set_epi32(0x80808080, 0x80808080, 0x80808080, 0x0c080400);
	__m128i pixels;
	__m128* p_p;
	__m128* p_pt = (__m128*)levels[level].data;

	__m128i* dst_pp_m;
	int* dst_pp_i;
	uint8_t* dst_pp_b;

	int m128_pitch = levels[level].pitch >> 2;

	int sixteens = levels[level].width >> 4;
	int fours = (levels[level].width >> 2) - (sixteens << 2);
	int singles = levels[level].width & 3;

	if (level) {
		if (sy == 0) sy++;
		if (ey == levels[level].height) ey--;
	}

	dst_p += (size_t)(sy - (level ? 1 : 0)) * pitch;

	p_pt += (size_t)m128_pitch * sy;

	__m128 dither_add;

	for (y = sy; y < ey; y++) {
		switch (y & 3) {
			case 0: dither_add = _mm_set_ps(0.4999f, 0.0f, 0.375f, -0.125f); break;
			case 1: dither_add = _mm_set_ps(-0.25f, 0.25f, -0.375f, 0.125f); break;
			case 2: dither_add = _mm_set_ps(0.3125f, -0.1875f, 0.4375f, -0.0625f); break;
			case 3: dither_add = _mm_set_ps(-0.4375f, 0.0625f, -0.3125f, 0.1875f); break;
		}

		dst_pp_m = (__m128i*)dst_p;
		p_p = p_pt;
		for (x = 0; x < sixteens; ++x) {
			pixels = _mm_shuffle_epi8(_mm_cvtps_epi32(LOAD), shuffle1);
			pixels = _mm_or_si128(pixels, _mm_shuffle_epi8(_mm_cvtps_epi32(LOAD), shuffle2));
			pixels = _mm_or_si128(pixels, _mm_shuffle_epi8(_mm_cvtps_epi32(LOAD), shuffle3));
			pixels = _mm_or_si128(pixels, _mm_shuffle_epi8(_mm_cvtps_epi32(LOAD), shuffle4));

			_mm_storeu_si128(dst_pp_m++, pixels);
		}

		dst_pp_i = (int*)dst_pp_m;
		for (x = 0; x < fours; ++x) {
			*dst_pp_i++ = _mm_cvtsi128_si32(_mm_shuffle_epi8(_mm_cvtps_epi32(LOAD), four_shuffle));
		}

		if (singles) {
			dst_pp_b = (uint8_t*)dst_pp_i;
			__m128i a;
			a = _mm_cvtps_epi32(LOAD);
			*dst_pp_b++ = (uint8_t)_mm_extract_epi8(a, 0);
			if (singles > 1) {
				*dst_pp_b++ = (uint8_t)_mm_extract_epi8(a, 4);
				if (singles == 3) *dst_pp_b++ = (uint8_t)_mm_extract_epi8(a, 8);
			}
		}

		if (level) {
			memcpy(dst_p, dst_p + 1, levels[level].width - 2);
			dst_p[levels[level].width - 1] = (uint8_t)black;
			dst_p[levels[level].width - 2] = (uint8_t)black;
		}

		p_pt += m128_pitch;
		dst_p += pitch;
	}
};

/***********************************************************************
* 16-bit planar
***********************************************************************/
template <typename F>
void Pyramid::OutPlanar16 _OP_{
	int x, y;
	uint16_t* dst_p = (uint16_t*)_dst_p;
	uint16_t black = chroma ? 0x8000 : 0x0000;

	__m128 zeroes = _mm_setzero_ps();
	__m128 maxes = _mm_set_ps1(65535.0f);

	__m128i shuffle1 = _mm_set_epi32(0x80808080, 0x80808080, 0x0d0c0908, 0x05040100);
	__m128i shuffle2 = _mm_set_epi32(0x0d0c0908, 0x05040100, 0x80808080, 0x80808080);
	__m128i pixels;
	__m128* p_p;
	__m128* p_pt = (__m128*)levels[level].data;

	__m128i* dst_pp_m;
	uint16_t* dst_pp_w;

	int m128_pitch = levels[level].pitch >> 2;

	int eights = levels[level].width >> 3;
	int four = levels[level].width & 4;
	int singles = levels[level].width & 3;

	if (level) {
		if (sy == 0) sy++;
		if (ey == levels[level].height) ey--;
	}

	dst_p += (sy - (level ? 1 : 0)) * pitch;

	p_pt += m128_pitch * sy;

	__m128 dither_add;

	for (y = sy; y < ey; y++) {
		switch (y & 3) {
			case 0: dither_add = _mm_set_ps(0.4999f, 0.0f, 0.375f, -0.125f); break;
			case 1: dither_add = _mm_set_ps(-0.25f, 0.25f, -0.375f, 0.125f); break;
			case 2: dither_add = _mm_set_ps(0.3125f, -0.1875f, 0.4375f, -0.0625f); break;
			case 3: dither_add = _mm_set_ps(-0.4375f, 0.0625f, -0.3125f, 0.1875f); break;
		}

		dst_pp_m = (__m128i*)dst_p;
		p_p = p_pt;
		for (x = 0; x < eights; ++x) {
			pixels = _mm_shuffle_epi8(_mm_cvtps_epi32(LOAD), shuffle1);
			pixels = _mm_or_si128(pixels, _mm_shuffle_epi8(_mm_cvtps_epi32(LOAD), shuffle2));

			_mm_storeu_si128(dst_pp_m++, pixels);
		}

		__m128i a;

		dst_pp_w = (uint16_t*)dst_pp_m;
		if (four) {
			a = _mm_cvtps_epi32(LOAD);
			*dst_pp_w++ = (uint16_t)_mm_extract_epi16(a, 0);
			*dst_pp_w++ = (uint16_t)_mm_extract_epi16(a, 2);
			*dst_pp_w++ = (uint16_t)_mm_extract_epi16(a, 4);
			*dst_pp_w++ = (uint16_t)_mm_extract_epi16(a, 6);
		}

		if (singles) {
			a = _mm_cvtps_epi32(LOAD);
			*dst_pp_w++ = (uint16_t)_mm_extract_epi16(a, 0);
			if (singles > 1) {
				*dst_pp_w++ = (uint16_t)_mm_extract_epi16(a, 2);
				if (singles == 3) *dst_pp_w++ = (uint16_t)_mm_extract_epi16(a, 4);
			}
		}

		if (level) {
			memcpy(dst_p, dst_p + 1, (levels[level].width - 2) << 1);
			dst_p[levels[level].width - 1] = black;
			dst_p[levels[level].width - 2] = black;
		}

		p_pt += m128_pitch;
		dst_p += pitch;
	}
}

/***********************************************************************
* 32-bit planar
***********************************************************************/
#undef LOAD
#define LOAD _loader((float*)&p_p[x], dither_add, zeroes, maxes)

template <typename F>
void Pyramid::OutPlanar32 _OP_{
	int x, y;
	__m128* dst_p = (__m128*)_dst_p;

	__m128 zeroes;
	__m128 maxes;
	__m128 dither_add = _mm_set_ps1(0.0f);

	pitch >>= 2; // number of floats to number of __m128s

	if (chroma) {
		zeroes = _mm_set_ps1(-0.5f);
		maxes = _mm_set_ps1(0.5f);
	} else {
		zeroes = _mm_set_ps1(0.0f);
		maxes = _mm_set_ps1(1.0f);
	}

	__m128* p_p = (__m128*)levels[level].data;

	float* dst_pp_f;

	int m128_pitch = levels[level].pitch >> 2;

	int fours = (levels[level].width + 3) >> 2;
	int wipes = (fours << 2) - levels[level].width;

	if (level) {
		if (sy == 0) sy++;
		if (ey == levels[level].height) ey--;
	}

	dst_p += (sy - (level ? 1 : 0)) * pitch;

	p_p += m128_pitch * sy;

	for (y = sy; y < ey; y++) {
		for (x = 0; x < fours; x++) {
			dst_p[x] = LOAD;
		}

		if (wipes) {
			dst_pp_f = (float*)&dst_p[x];
			int w = wipes;
			while (w) {
				*--dst_pp_f = 0;
				w--;
			}
		}

		if (level) {
			memcpy(dst_p, ((float*)dst_p) + 1, (levels[level].width - 2) << 2);
			((float*)dst_p)[levels[level].width - 1] = 0.0f;
			((float*)dst_p)[levels[level].width - 2] = 0.0f;
		}

		p_p += m128_pitch;
		dst_p += pitch;
	}
}

/***********************************************************************
* Interleaved
***********************************************************************/
template <typename T, typename F>
void Pyramid::OutInterleaved(T dst_p, F _loader, int pitch, int sy, int ey, int level, bool chroma, int step, int offset) {
	int x, y;

	__m128 zeroes = _mm_setzero_ps();
	__m128 maxes = _mm_set_ps1(sizeof(*dst_p) == 1 ? 255.0f : 65535.0f);

	if (level) {
		if (sy == 0) sy++;
		if (ey == levels[level].height) ey--;
	}

	dst_p += (sy - (level ? 1 : 0)) * pitch + offset;

	int m128_pitch = levels[level].pitch >> 2;
	__m128* p_p = (__m128*)levels[level].data + m128_pitch * sy;
	T dst_pp;

	__m128i a;
	int fours = levels[level].width >> 2;
	int singles = levels[level].width & 3;

	if (level) {
		singles--;
		if (singles < 0) {
			singles = 3;
			fours--;
		}
	}

	__m128 dither_add;

	// loop
	for (y = sy; y < ey; y++) {
		switch (y & 3) {
			case 0: dither_add = _mm_set_ps(0.4999f, 0.0f, 0.375f, -0.125f); break;
			case 1: dither_add = _mm_set_ps(-0.25f, 0.25f, -0.375f, 0.125f); break;
			case 2: dither_add = _mm_set_ps(0.3125f, -0.1875f, 0.4375f, -0.0625f); break;
			case 3: dither_add = _mm_set_ps(-0.4375f, 0.0625f, -0.3125f, 0.1875f); break;
		}

		dst_pp = dst_p;
		x = 0;
		if (level) {
			a = _mm_cvtps_epi32(LOAD);
			*dst_pp = _mm_extract_epi16(a, 2); dst_pp += step;
			*dst_pp = _mm_extract_epi16(a, 4); dst_pp += step;
			*dst_pp = _mm_extract_epi16(a, 6); dst_pp += step;
			x++;
		}

		for (; x < fours; x++) {
			a = _mm_cvtps_epi32(LOAD);
			*dst_pp = _mm_extract_epi16(a, 0); dst_pp += step;
			*dst_pp = _mm_extract_epi16(a, 2); dst_pp += step;
			*dst_pp = _mm_extract_epi16(a, 4); dst_pp += step;
			*dst_pp = _mm_extract_epi16(a, 6); dst_pp += step;
		}

		if (singles) {
			a = _mm_cvtps_epi32(LOAD);
			*dst_pp = _mm_extract_epi8(a, 0);
			if (singles > 1) {
				dst_pp += step; *dst_pp = _mm_extract_epi16(a, 2);
				if (singles == 3) {
					dst_pp += step; *dst_pp = _mm_extract_epi16(a, 2);
				}
			}
		}

		p_p += m128_pitch;
		dst_p += pitch;
	}
}

/***********************************************************************
* Out
***********************************************************************/
template <typename T>
void Pyramid::Out(T dst_p, int pitch, bool gamma, bool dither, bool clamp, int level, int step, int offset, bool chroma) {
	int bytes = sizeof(*dst_p);
	int eb = (int)(levels[level].bands.size() - 1);

	typedef typename std::conditional<sizeof(*dst_p) == 1, uint8_t*, uint16_t*>::type Type; // used to avoid generating a float* version of OutInterleaved which would cause warnings

	int s = (gamma ? 1 : 0) | (dither && bytes != 4 ? 2 : 0) | (clamp ? 4 : 0);

	if (step) { // interleaved
		for (int b = 0; b < eb; ++b) {
			switch (s) {
				case 0: threadpool->Queue([=] { OutInterleaved((Type)dst_p, PLL(N), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma, step, offset); }); break;
				case 1: threadpool->Queue([=] { OutInterleaved((Type)dst_p, PLL(G), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma, step, offset); }); break;
				case 2: threadpool->Queue([=] { OutInterleaved((Type)dst_p, PLL(D), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma, step, offset); }); break;
				case 3: threadpool->Queue([=] { OutInterleaved((Type)dst_p, PLL(DG), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma, step, offset); }); break;
				case 4: threadpool->Queue([=] { OutInterleaved((Type)dst_p, PLL(C), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma, step, offset); }); break;
				case 5: threadpool->Queue([=] { OutInterleaved((Type)dst_p, PLL(CG), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma, step, offset); }); break;
				case 6: threadpool->Queue([=] { OutInterleaved((Type)dst_p, PLL(CD), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma, step, offset); }); break;
				case 7: threadpool->Queue([=] { OutInterleaved((Type)dst_p, PLL(CDG), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma, step, offset); }); break;
			}
		}
	} else { // planar
		switch (bytes) {
			case 1: {
				for (int b = 0; b < eb; ++b) {
					switch (s) {
						case 0: threadpool->Queue([=] { OutPlanar8(dst_p, PLL(N), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 1: threadpool->Queue([=] { OutPlanar8(dst_p, PLL(G), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 2: threadpool->Queue([=] { OutPlanar8(dst_p, PLL(D), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 3: threadpool->Queue([=] { OutPlanar8(dst_p, PLL(DG), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 4: threadpool->Queue([=] { OutPlanar8(dst_p, PLL(C), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 5: threadpool->Queue([=] { OutPlanar8(dst_p, PLL(CG), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 6: threadpool->Queue([=] { OutPlanar8(dst_p, PLL(CD), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 7: threadpool->Queue([=] { OutPlanar8(dst_p, PLL(CDG), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
					}
				}
			} break;
			case 2: {
				for (int b = 0; b < eb; ++b) {
					switch (s) {
						case 0: threadpool->Queue([=] { OutPlanar16(dst_p, PLL(N), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 1: threadpool->Queue([=] { OutPlanar16(dst_p, PLL(G), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 2: threadpool->Queue([=] { OutPlanar16(dst_p, PLL(D), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 3: threadpool->Queue([=] { OutPlanar16(dst_p, PLL(DG), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 4: threadpool->Queue([=] { OutPlanar16(dst_p, PLL(C), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 5: threadpool->Queue([=] { OutPlanar16(dst_p, PLL(CG), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 6: threadpool->Queue([=] { OutPlanar16(dst_p, PLL(CD), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 7: threadpool->Queue([=] { OutPlanar16(dst_p, PLL(CDG), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
					}
				}
			} break;
			case 4: {
				for (int b = 0; b < eb; ++b) {
					switch (s) {
						case 0: threadpool->Queue([=] { OutPlanar32(dst_p, PLL(N), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 1: threadpool->Queue([=] { OutPlanar32(dst_p, PLL(G), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 4: threadpool->Queue([=] { OutPlanar32(dst_p, PLL(C), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
						case 5: threadpool->Queue([=] { OutPlanar32(dst_p, PLL(CG), pitch, levels[level].bands[b], levels[level].bands[b + 1], level, chroma); }); break;
					}
				}
			} break;
		}
	}

	threadpool->Wait();
}

#undef PLL
#undef   N
#undef   G
#undef   D
#undef   C
#undef  DG
#undef  CG
#undef  CD
#undef CDG
