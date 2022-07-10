/*
	Multiblend 2.0 (c) 2021 David Horman

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program. If not, see <https://www.gnu.org/licenses/>.

	The author can be contacted at davidhorman47@gmail.com
*/

#define NOMINMAX
#include <stdio.h>
#include <stdint.h>
#include <vector>
#include <algorithm>
#ifdef __APPLE__
#define memalign(a,b) malloc((b))
#else
#include <malloc.h>
#endif

#include "tiffio.h"
#include "jpeglib.h"

#ifndef _WIN32
#include <strings.h>
int _stricmp(const char* a, const char* b) { return strcasecmp(a, b); }
#define ZeroMemory(a,b) memset(a,0,b)
#define sprintf_s sprintf
#define sscanf_s sscanf
void* _aligned_malloc(size_t size, int boundary) { return memalign(boundary, size); }
void _aligned_free(void* a) { free(a); }
void fopen_s(FILE** f, const char* filename, const char* mode) { *f = fopen(filename, mode); }
#endif

int verbosity = 1;

#include "pnger.cpp"
#include "pyramid.cpp"
#include "functions.cpp"

#include "mapalloc.cpp"
#include "threadpool.cpp"
#include "geotiff.cpp"

class PyramidWithMasks : public Pyramid {
public:
	using Pyramid::Pyramid;
	std::vector<Flex*> masks;
};

enum class ImageType { MB_NONE, MB_TIFF, MB_JPEG, MB_PNG };

#include "image.cpp"

#ifdef _WIN32
FILE _iob[] = { *stdin, *stdout, *stderr };

extern "C" FILE * __cdecl __iob_func(void) {
	return _iob;
}
#pragma comment(lib, "legacy_stdio_definitions.lib")
// the above are required to support VS 2010 build of libjpeg-turbo 2.0.6
#pragma comment(lib, "tiff.lib")
#pragma comment(lib, "turbojpeg-static.lib")
#pragma comment(lib, "libpng16.lib")
#pragma comment(lib, "zlib.lib")
#endif

#define MASKVAL(X) (((X) & 0x7fffffffffffffff) | images[(X) & 0xffffffff]->mask_state)

int main(int argc, char* argv[]) {
// This is here because of a weird problem encountered during development with Visual Studio. It should never be triggered.
	if (verbosity != 1) {
		printf("bad compile?\n");
		exit(EXIT_FAILURE);
	}

	int i;
	Timer timer_all, timer;
	timer_all.Start();

	TIFFSetWarningHandler(NULL);

/***********************************************************************
* Variables
***********************************************************************/
	std::vector<Image*> images;
	int fixed_levels = 0;
	int add_levels = 0;

	int width = 0;
	int height = 0;

	bool no_mask = false;
	bool big_tiff = false;
	bool bgr = false;
	bool wideblend = false;
	bool reverse = false;
	bool timing = false;
	bool dither = true;
	bool gamma = false;
	bool all_threads = true;
	int wrap = 0;

	TIFF* tiff_file = NULL;
	FILE* jpeg_file = NULL;
	Pnger* png_file = NULL;
	ImageType output_type = ImageType::MB_NONE;
	int jpeg_quality = -1;
	int compression = -1;
	char* seamsave_filename = NULL;
	char* seamload_filename = NULL;
	char* xor_filename = NULL;
	char* output_filename = NULL;
	int output_bpp = 0;

	double images_time = 0;
	double copy_time = 0;
	double seam_time = 0;
	double shrink_mask_time = 0;
	double shrink_time = 0;
	double laplace_time = 0;
	double blend_time = 0;
	double collapse_time = 0;
	double wrap_time = 0;
	double out_time = 0;
	double write_time = 0;

/***********************************************************************
* Help
***********************************************************************/
	if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help") || !strcmp(argv[1], "/?")) {
		Output(1, "\n");
		Output(1, "Multiblend v2.0.0 (c) 2021 David Horman        http://horman.net/multiblend/\n");
		Output(1, "----------------------------------------------------------------------------\n");

		printf("Usage: multiblend [options] [-o OUTPUT] INPUT [X,Y] [INPUT] [X,Y] [INPUT]...\n");
		printf("\n");
		printf("Options:\n");
		printf("  --levels X / -l X      X: set number of blending levels to X\n");
		printf("                        -X: decrease number of blending levels by X\n");
		printf("                        +X: increase number of blending levels by X\n");
		printf("  --depth D / -d D       Override automatic output image depth (8 or 16)\n");
		printf("  --bgr                  Swap RGB order\n");
		printf("  --wideblend            Calculate number of levels based on output image size,\n");
		printf("                         rather than input image size\n");
		printf("  -w, --wrap=[mode]      Blend around images boundaries (NONE (default),\n");
		printf("                         HORIZONTAL, VERTICAL). When specified without a mode,\n");
		printf("                         defaults to HORIZONTAL.\n");
		printf("  --compression=X        Output file compression. For TIFF output, X may be:\n");
		printf("                         NONE (default), PACKBITS, or LZW\n");
		printf("                         For JPEG output, X is JPEG quality (0-100, default 75)\n");
		printf("                         For PNG output, X is PNG filter (0-9, default 3)\n");
		printf("  --cache-threshold=     Allocate memory beyond X bytes/[K]ilobytes/\n");
		printf("      X[K/M/G]           [M]egabytes/[G]igabytes to disk\n");
		printf("  --no-dither            Disable dithering\n");
		printf("  --tempdir <dir>        Specify temporary directory (default: system temp)\n");
		printf("  --save-seams <file>    Save seams to PNG file for external editing\n");
		printf("  --load-seams <file>    Load seams from PNG file\n");
		printf("  --no-output            Do not blend (for use with --save-seams)\n");
		printf("                         Must be specified as last option before input images\n");
		printf("  --bigtiff              BigTIFF output\n");
		printf("  --reverse              Reverse image priority (last=highest) for resolving\n");
		printf("                         indeterminate pixels\n");
		printf("  --quiet                Suppress output (except warnings)\n");
		printf("  --all-threads          Use all available CPU threads\n");
		printf("  [X,Y]                  Optional position adjustment for previous input image\n");
		exit(EXIT_SUCCESS);
	}

/***********************************************************************
************************************************************************
* Parse arguments
************************************************************************
***********************************************************************/
	std::vector<char*> my_argv;

	bool skip = false;

	for (i = 1; i < argc; ++i) {
		my_argv.push_back(argv[i]);

		if (!skip) {
			int c = 0;

			while (argv[i][c]) {
				if (argv[i][c] == '=') {
					argv[i][c++] = 0;
					if (argv[i][c]) {
						my_argv.push_back(&argv[i][c]);
					}
					break;
				}
				++c;
			}

			if (!strcmp(argv[i], "-o") || !strcmp(argv[i], "--output")) {
				skip = true;
			}
		}
	}

	if ((int)my_argv.size() < 3) die("Error: Not enough arguments (try -h for help)");

	for (i = 0; i < (int)my_argv.size(); ++i) {
		if (!strcmp(my_argv[i], "-d") || !strcmp(my_argv[i], "--d") || !strcmp(my_argv[i], "--depth") || !strcmp(my_argv[i], "--bpp")) {
			if (++i < (int)my_argv.size()) {
				output_bpp = atoi(my_argv[i]);
				if (output_bpp != 8 && output_bpp != 16) {
					die("Error: Invalid output depth specified");
				}
			} else {
				die("Error: Missing parameter value");
			}
		} 		else if (!strcmp(my_argv[i], "-l") || !strcmp(my_argv[i], "--levels")) {
			if (++i < (int)my_argv.size()) {
				int n;
				if (my_argv[i][0] == '+' || my_argv[i][0] == '-') {
					sscanf_s(my_argv[i], "%d%n", &add_levels, &n);
				} else {
					sscanf_s(my_argv[i], "%d%n", &fixed_levels, &n);
					if (fixed_levels == 0) fixed_levels = 1;
				}
				if (my_argv[i][n]) die("Error: Bad --levels parameter");
			} else {
				die("Error: Missing parameter value");
			}
		} 		else if (!strcmp(my_argv[i], "--wrap") || !strcmp(my_argv[i], "-w")) {
			if (i + 1 >= (int)my_argv.size()) {
				die("Error: Missing parameters");
			}
			if (!strcmp(my_argv[i + 1], "none") || !strcmp(my_argv[i + 1], "open")) ++i;
			else if (!strcmp(my_argv[i + 1], "horizontal") || !strcmp(my_argv[i + 1], "h")) { wrap = 1; ++i; } else if (!strcmp(my_argv[i + 1], "vertical") || !strcmp(my_argv[i + 1], "v")) { wrap = 2; ++i; } else if (!strcmp(my_argv[i + 1], "both") || !strcmp(my_argv[i + 1], "hv")) { wrap = 3; ++i; } else wrap = 1;
		} 		else if (!strcmp(my_argv[i], "--cache-threshold")) {
			if (i + 1 >= (int)my_argv.size()) {
				die("Error: Missing parameters");
			}
			++i;
			int shift = 0;
			int n = 0;
			size_t len = strlen(my_argv[i]);
			size_t threshold;
			sscanf_s(my_argv[i], "%zu%n", &threshold, &n);
			if (n != len) {
				if (n == len - 1) {
					switch (my_argv[i][len - 1]) {
						case 'k':
						case 'K': shift = 10; break;
						case 'm':
						case 'M': shift = 20; break;
						case 'g':
						case 'G': shift = 30; break;
						default: die("Error: Bad --cache-threshold parameter");
					}
					threshold <<= shift;
				} else {
					die("Error: Bad --cache-threshold parameter");
				}
			}
			MapAlloc::CacheThreshold(threshold);
		} 		else if (!strcmp(my_argv[i], "--nomask") || !strcmp(my_argv[i], "--no-mask")) no_mask = true;
		else if (!strcmp(my_argv[i], "--timing") || !strcmp(my_argv[i], "--timings")) timing = true;
		else if (!strcmp(my_argv[i], "--bigtiff"))   big_tiff = true;
		else if (!strcmp(my_argv[i], "--bgr"))       bgr = true;
		else if (!strcmp(my_argv[i], "--wideblend")) wideblend = true;
		else if (!strcmp(my_argv[i], "--reverse"))   reverse = true;
		else if (!strcmp(my_argv[i], "--gamma"))     gamma = true;
		else if (!strcmp(my_argv[i], "--no-dither") || !strcmp(my_argv[i], "--nodither")) dither = false;
		//		else if (!strcmp(my_argv[i], "--force"))     force_coverage = true;
		else if (!strncmp(my_argv[i], "-f", 2)) Output(0, "ignoring Enblend option -f\n");
		else if (!strcmp(my_argv[i], "-a")) Output(0, "ignoring Enblend option -a\n");
		else if (!strcmp(my_argv[i], "--no-ciecam")) Output(0, "ignoring Enblend option --no-ciecam\n");
		else if (!strcmp(my_argv[i], "--primary-seam-generator")) {
			Output(0, "ignoring Enblend option --primary-seam-generator\n");
			++i;
		}

		else if (!strcmp(my_argv[i], "--compression")) {
			if (++i < (int)my_argv.size()) {
				if (strcmp(my_argv[i], "0") == 0) jpeg_quality = 0;
				else if (atoi(my_argv[i]) > 0) jpeg_quality = atoi(my_argv[i]);
				else if (_stricmp(my_argv[i], "lzw") == 0) compression = COMPRESSION_LZW;
				else if (_stricmp(my_argv[i], "packbits") == 0) compression = COMPRESSION_PACKBITS;
				//				else if (_stricmp(my_argv[i], "deflate") == 0) compression = COMPRESSION_DEFLATE;
				else if (_stricmp(my_argv[i], "none") == 0) compression = COMPRESSION_NONE;
				else die("Error: Unknown compression codec %s", my_argv[i]);
			} else {
				die("Error: Missing parameter value");
			}
		} else if (!strcmp(my_argv[i], "-v") || !strcmp(my_argv[i], "--verbose")) ++verbosity;
		else if (!strcmp(my_argv[i], "-q") || !strcmp(my_argv[i], "--quiet")) --verbosity;
		else if ((!strcmp(my_argv[i], "--saveseams") || !strcmp(my_argv[i], "--save-seams")) && i < (int)my_argv.size() - 1) seamsave_filename = my_argv[++i];
		else if ((!strcmp(my_argv[i], "--loadseams") || !strcmp(my_argv[i], "--load-seams")) && i < (int)my_argv.size() - 1) seamload_filename = my_argv[++i];
		else if ((!strcmp(my_argv[i], "--savexor") || !strcmp(my_argv[i], "--save-xor")) && i < (int)my_argv.size() - 1) xor_filename = my_argv[++i];
		else if (!strcmp(my_argv[i], "--tempdir") || !strcmp(my_argv[i], "--tmpdir") && i < (int)my_argv.size() - 1) MapAlloc::SetTmpdir(my_argv[++i]);
		else if (!strcmp(my_argv[i], "--all-threads")) all_threads = true;
		else if (!strcmp(my_argv[i], "-o") || !strcmp(my_argv[i], "--output")) {
			if (++i < (int)my_argv.size()) {
				output_filename = my_argv[i];
				char* ext = strrchr(output_filename, '.');

				if (!ext) {
					die("Error: Unknown output filetype");
				}

				++ext;
				if (!(_stricmp(ext, "jpg") && _stricmp(ext, "jpeg"))) {
					output_type = ImageType::MB_JPEG;
					if (jpeg_quality == -1) jpeg_quality = 75;
				} else if (!(_stricmp(ext, "tif") && _stricmp(ext, "tiff"))) {
					output_type = ImageType::MB_TIFF;
				} else if (!_stricmp(ext, "png")) {
					output_type = ImageType::MB_PNG;
				} else {
					die("Error: Unknown file extension");
				}

				++i;
				break;
			}
		} else if (!strcmp(my_argv[i], "--no-output")) {
			++i;
			break;
		} else {
			die("Error: Unknown argument \"%s\"", my_argv[i]);
		}
	}

	if (compression != -1) {
		if (output_type != ImageType::MB_TIFF) {
			Output(0, "Warning: non-TIFF output; ignoring TIFF compression setting\n");
		}
	} else if (output_type == ImageType::MB_TIFF) {
		compression = COMPRESSION_LZW;
	}

	if (jpeg_quality != -1 && output_type != ImageType::MB_JPEG && output_type != ImageType::MB_PNG) {
		Output(0, "Warning: non-JPEG/PNG output; ignoring compression quality setting\n");
	}

	if ((jpeg_quality < -1 || jpeg_quality > 9) && output_type == ImageType::MB_PNG) {
		die("Error: Bad PNG compression quality setting\n");
	}

	if (output_type == ImageType::MB_NONE && !seamsave_filename) die("Error: No output file specified");
	if (seamload_filename && seamsave_filename) die("Error: Cannot load and save seams at the same time");
	if (wrap == 3) die("Error: Wrapping in both directions is not currently supported");

	if (!strcmp(my_argv[i], "--")) ++i;

/***********************************************************************
* Push remaining arguments to images vector
***********************************************************************/
	int x, y, n;

	while (i < (int)my_argv.size()) {
		if (images.size()) {
			n = 0;
			sscanf_s(my_argv[i], "%d,%d%n", &x, &y, &n);
			if (!my_argv[i][n]) {
				images.back()->xpos_add = x;
				images.back()->ypos_add = y;
				i++;
				continue;
			}
		}
		images.push_back(new Image(my_argv[i++]));
	}

	int n_images = (int)images.size();

	if (n_images == 0) die("Error: No input files specified");
	if (seamsave_filename && n_images > 256) { seamsave_filename = NULL; Output(0, "Warning: seam saving not possible with more than 256 images"); }
	if (seamload_filename && n_images > 256) { seamload_filename = NULL; Output(0, "Warning: seam loading not possible with more than 256 images"); }
	if (xor_filename && n_images > 255) { xor_filename = NULL; Output(0, "Warning: XOR map saving not possible with more than 255 images"); }

/***********************************************************************
* Print banner
***********************************************************************/
	Output(1, "\n");
	Output(1, "Multiblend v2.0.0 (c) 2021 David Horman        http://horman.net/multiblend/\n");
	Output(1, "----------------------------------------------------------------------------\n");

	Threadpool* threadpool = Threadpool::GetInstance(all_threads ? 2 : 0);

/***********************************************************************
************************************************************************
* Open output
************************************************************************
***********************************************************************/
	switch (output_type) {
		case ImageType::MB_TIFF: {
			if (!big_tiff) tiff_file = TIFFOpen(output_filename, "w"); else tiff_file = TIFFOpen(output_filename, "w8");
			if (!tiff_file) die("Error: Could not open output file");
		} break;
		case ImageType::MB_JPEG: {
			if (output_bpp == 16) die("Error: 16bpp output is incompatible with JPEG output");
			fopen_s(&jpeg_file, output_filename, "wb");
			if (!jpeg_file) die("Error: Could not open output file");
		} break;
		case ImageType::MB_PNG: {
			fopen_s(&jpeg_file, output_filename, "wb");
			if (!jpeg_file) die("Error: Could not open output file");
		} break;
	}

/***********************************************************************
************************************************************************
* Process images
************************************************************************
***********************************************************************/
	timer.Start();

/***********************************************************************
* Open images to get prelimary info
***********************************************************************/
	size_t untrimmed_bytes = 0;

	for (i = 0; i < n_images; ++i) {
		images[i]->Open();
		untrimmed_bytes = std::max(untrimmed_bytes, images[i]->untrimmed_bytes);
	}

/***********************************************************************
* Check paramters, display warnings
***********************************************************************/
	for (i = 1; i < n_images; ++i) {
		if (images[i]->tiff_xres != images[0]->tiff_xres || images[i]->tiff_yres != images[0]->tiff_yres) {
			Output(0, "Warning: TIFF resolution mismatch (%f %f/%f %f)\n", images[0]->tiff_xres, images[0]->tiff_yres, images[i]->tiff_xres, images[i]->tiff_yres);
		}
	}

	for (i = 0; i < n_images; ++i) {
		if (output_bpp == 0 && images[i]->bpp == 16) output_bpp = 16;
		if (images[i]->bpp != images[0]->bpp) {
			die("Error: mixture of 8bpp and 16bpp images detected (not currently handled)\n");
		}
	}

	if (output_bpp == 0) output_bpp = 8;
	else if (output_bpp == 16 && output_type == ImageType::MB_JPEG) {
		Output(0, "Warning: 8bpp output forced by JPEG output\n");
		output_bpp = 8;
	}

/***********************************************************************
* Allocate working space for reading/trimming/extraction
***********************************************************************/
	void* untrimmed_data = MapAlloc::Alloc(untrimmed_bytes);

/***********************************************************************
* Read/trim/extract
***********************************************************************/
	for (i = 0; i < n_images; ++i) {
		try {
			images[i]->Read(untrimmed_data, gamma);
		} catch (char* e) {
			printf("\n\n");
			printf("%s\n", e);
			exit(EXIT_FAILURE);
		}
	}

/***********************************************************************
* Clean up
***********************************************************************/
	MapAlloc::Free(untrimmed_data);

/***********************************************************************
* Tighten
***********************************************************************/
	int min_xpos = 0x7fffffff;
	int min_ypos = 0x7fffffff;
	width = 0;
	height = 0;

	for (i = 0; i < n_images; ++i) {
		min_xpos = std::min(min_xpos, images[i]->xpos);
		min_ypos = std::min(min_ypos, images[i]->ypos);
	}

	for (i = 0; i < n_images; ++i) {
		images[i]->xpos -= min_xpos;
		images[i]->ypos -= min_ypos;
		width = std::max(width, images[i]->xpos + images[i]->width);
		height = std::max(height, images[i]->ypos + images[i]->height);
	}

	images_time = timer.Read();

/***********************************************************************
* Determine number of levels
***********************************************************************/
	int blend_wh;
	int blend_levels;

	if (!fixed_levels) {
		if (!wideblend) {
			std::vector<int> widths;
			std::vector<int> heights;

			for (auto image : images) {
				widths.push_back(image->width);
				heights.push_back(image->height);
			}

			std::sort(widths.begin(), widths.end());
			std::sort(heights.begin(), heights.end());

			size_t halfway = (widths.size() - 1) >> 1;

			blend_wh = std::max(
				widths.size() & 1 ? widths[halfway] : (widths[halfway] + widths[halfway + 1] + 1) >> 1,
				heights.size() & 1 ? heights[halfway] : (heights[halfway] + heights[halfway + 1] + 1) >> 1
			);
		} else {
			blend_wh = (std::max)(width, height);
		}

		blend_levels = (int)floor(log2(blend_wh + 4.0f) - 1);
		if (wideblend) blend_levels++;
	} else {
		blend_levels = fixed_levels;
	}

	blend_levels += add_levels;

	if (n_images == 1) {
		blend_levels = 0;
		Output(1, "\n%d x %d, %d bpp\n\n", width, height, output_bpp);
	} else {
		Output(1, "\n%d x %d, %d levels, %d bpp\n\n", width, height, blend_levels, output_bpp);
	}
	
/***********************************************************************
************************************************************************
* Seaming
************************************************************************
***********************************************************************/
	timer.Start();

	Output(1, "Seaming");
	switch (((!!seamsave_filename) << 1) | !!xor_filename) {
		case 1: Output(1, " (saving XOR map)"); break;
		case 2: Output(1, " (saving seam map)"); break;
		case 3: Output(1, " (saving XOR and seam maps)"); break;
	}
	Output(1, "...\n");

	int min_count;
	int xor_count;
	int xor_image;
	uint64_t utemp;
	int stop;

	uint64_t best;
	uint64_t a, b, c, d;

#define DT_MAX 0x9000000000000000
	uint64_t* prev_line = NULL;
	uint64_t* this_line = NULL;
	bool last_pixel = false;
	bool arbitrary_seam = false;

	Flex* seam_flex = new Flex(width, height);
	int max_queue = 0;

/***********************************************************************
* Backward distance transform
***********************************************************************/
	int n_threads = std::max(2, threadpool->GetNThreads());
	uint64_t** thread_lines = new uint64_t*[n_threads];

	if (!seamload_filename) {
		std::mutex* flex_mutex_p = new std::mutex;
		std::condition_variable* flex_cond_p = new std::condition_variable;

		uint8_t** thread_comp_lines = new uint8_t*[n_threads];

		for (i = 0; i < n_threads; ++i) {
			thread_lines[i] = new uint64_t[width];
			thread_comp_lines[i] = new uint8_t[width];
		}

// set all image masks to bottom right
		for (i = 0; i < n_images; ++i) {
			images[i]->tiff_mask->End();
		}

		for (y = height - 1; y >= 0; --y) {
			int t = y % n_threads;
			this_line = thread_lines[t];
			uint8_t* comp = thread_comp_lines[t];

// set initial image mask states
			for (i = 0; i < n_images; ++i) {
				images[i]->mask_state = 0x8000000000000000;
				if (y >= images[i]->ypos && y < images[i]->ypos + images[i]->height) {
					images[i]->mask_count = width - (images[i]->xpos + images[i]->width);
					images[i]->mask_limit = images[i]->xpos;
				} else {
					images[i]->mask_count = width;
					images[i]->mask_limit = width;
				}
			}

			x = width - 1;

			{ // make sure the last compression thread to use this chunk of memory is finished
				std::unique_lock<std::mutex> mlock(*flex_mutex_p);
				flex_cond_p->wait(mlock, [=] { return seam_flex->y > (height - 1) - y - n_threads; });
			}

			while (x >= 0) {
				min_count = x + 1;
				xor_count = 0;

// update image mask states
				for (i = 0; i < n_images; ++i) {
					if (!images[i]->mask_count) {
						if (x >= images[i]->mask_limit) {
							utemp = images[i]->tiff_mask->ReadBackwards32();
							images[i]->mask_state = ((~utemp) << 32) & 0x8000000000000000;
							images[i]->mask_count = utemp & 0x7fffffff;
						} else {
							images[i]->mask_state = 0x8000000000000000;
							images[i]->mask_count = min_count;
						}
					}

					if (images[i]->mask_count < min_count) min_count = images[i]->mask_count;
					if (!images[i]->mask_state) { // mask_state is inverted
						++xor_count;
						xor_image = i;
					}
				}

				stop = x - min_count;

				if (xor_count == 1) {
					images[xor_image]->seam_present = true;
					while (x > stop) this_line[x--] = xor_image;
				} else {
					if (y == height - 1) { // bottom row
						if (x == width - 1) { // first pixel(s)
							while (x > stop) this_line[x--] = DT_MAX; // max
						} else {
							utemp = this_line[x + 1];
							utemp = MASKVAL(utemp);
							while (x > stop) {
								utemp += 0x300000000;
								this_line[x--] = utemp; // was min(temp, DT_MAX) but this is unlikely to happen
							}
						}
					} else { // other rows
						if (x == width - 1) { // first pixel(s)
							utemp = prev_line[x - 1] + 0x400000000;
							a = MASKVAL(utemp);

							utemp = prev_line[x] + 0x300000000;
							b = MASKVAL(utemp);

							d = a < b ? a : b;

							this_line[x--] = d;

							if (x == stop) {
								for (i = 0; i < n_images; ++i) {
									images[i]->mask_count -= min_count;
								}
								continue;
							}

							c = b + 0x100000000;
							b = a - 0x100000000;
							d += 0x300000000;
						} else {
							utemp = prev_line[x] + 0x300000000;
							b = MASKVAL(utemp);

							utemp = prev_line[x + 1] + 0x400000000;
							c = MASKVAL(utemp);

							utemp = this_line[x + 1] + 0x300000000;
							d = MASKVAL(utemp);
						}

						if (stop == -1) {
							stop = 0;
							last_pixel = true;
						}

						while (x > stop) {
							utemp = prev_line[x - 1] + 0x400000000;
							a = MASKVAL(utemp);

							if (a < d) d = a;
							if (b < d) d = b;
							if (c < d) d = c;
							
							this_line[x--] = d;

							c = b + 0x100000000;
							b = a - 0x100000000;
							d += 0x300000000;
						}

						if (last_pixel) {
							// d is the new "best" to compare against
							if (b < d) d = b;
							if (c < d) d = c;

							this_line[x--] = d;

							last_pixel = false;
						}
					}
				}

				for (i = 0; i < n_images; ++i) {
					images[i]->mask_count -= min_count;
				}
			}

			if (y) {
				threadpool->Queue([=] {
					int p = CompressSeamLine(this_line, comp, width);
					if (p > width) {
						printf("bad p: %d at line %d", p, y);
						exit(0);
					}

					{
						std::unique_lock<std::mutex> mlock(*flex_mutex_p);
						flex_cond_p->wait(mlock, [=] { return seam_flex->y == (height - 1) - y; });
						seam_flex->Copy(comp, p);
						seam_flex->NextLine();
					}
					flex_cond_p->notify_all();
				});
			}

			prev_line = this_line;
		} // end of row loop

		threadpool->Wait();

		for (i = 0; i < n_images; ++i) {
			if (!images[i]->seam_present) {
				Output(1, "Warning: %s is fully obscured by other images\n", images[i]->filename);
			}
		}

		for (i = 0; i < n_threads; ++i) {
			if (i >= 2) delete[] thread_lines[i];
			delete[] thread_comp_lines[i];
		}

		delete[] thread_comp_lines;
		delete flex_mutex_p;
		delete flex_cond_p;
	} else { // if seamload_filename:
		for (i = 0; i < n_images; ++i) {
			images[i]->tiff_mask->Start();
		}
	}

// create top level masks
	for (i = 0; i < n_images; ++i) {
		images[i]->masks.push_back(new Flex(width, height));
	}

	Pnger* xor_map = xor_filename ? new Pnger(xor_filename, "XOR map", width, height, PNG_COLOR_TYPE_PALETTE) : NULL;
	Pnger* seam_map = seamsave_filename ? new Pnger(seamsave_filename, "Seam map", width, height, PNG_COLOR_TYPE_PALETTE) : NULL;

/***********************************************************************
* Forward distance transform
***********************************************************************/
	int current_count = 0;
	int64 current_step;
	uint64_t dt_val;

	prev_line = thread_lines[1];

	uint64_t total_pixels = 0;
	uint64_t channel_totals[3] = { 0 };

	Flex full_mask(width, height);
	Flex xor_mask(width, height);

	bool alpha = false;

	for (y = 0; y < height; ++y) {
		for (i = 0; i < n_images; ++i) {
			images[i]->mask_state = 0x8000000000000000;
			if (y >= images[i]->ypos && y < images[i]->ypos + images[i]->height) {
				images[i]->mask_count = images[i]->xpos;
				images[i]->mask_limit = images[i]->xpos + images[i]->width;
			} else {
				images[i]->mask_count = width;
				images[i]->mask_limit = width;
			}
		}

		x = 0;
		int mc = 0;
		int prev_i = -1;
		int current_i = -1;
		int best_temp;

		while (x < width) {
			min_count = width - x;
			xor_count = 0;

			for (i = 0; i < n_images; ++i) {
				if (!images[i]->mask_count) {
					if (x < images[i]->mask_limit) {
						utemp = images[i]->tiff_mask->ReadForwards32();
						images[i]->mask_state = ((~utemp) << 32) & 0x8000000000000000;
						images[i]->mask_count = utemp & 0x7fffffff;
					} else {
						images[i]->mask_state = 0x8000000000000000;
						images[i]->mask_count = min_count;
					}
				}

				if (images[i]->mask_count < min_count) min_count = images[i]->mask_count;
				if (!images[i]->mask_state) {
					++xor_count;
					xor_image = i;
				}
			}

			stop = x + min_count;

			if (!xor_count) {
				alpha = true;
			}
			full_mask.MaskWrite(min_count, xor_count);
			xor_mask.MaskWrite(min_count, xor_count == 1);

			if (xor_count == 1) {
				if (xor_map) memset(&xor_map->line[x], xor_image, min_count);

				size_t p = (y - images[xor_image]->ypos) * images[xor_image]->width + (x - images[xor_image]->xpos);

				int total_count = min_count;
				total_pixels += total_count;
				if (gamma) {
					switch (images[xor_image]->bpp) {
						case 8: {
							uint16_t v;
							while (total_count--) {
								v = ((uint8_t*)images[xor_image]->channels[0]->data)[p];
								channel_totals[0] += v * v;
								v = ((uint8_t*)images[xor_image]->channels[1]->data)[p];
								channel_totals[1] += v * v;
								v = ((uint8_t*)images[xor_image]->channels[2]->data)[p];
								channel_totals[2] += v * v;
								++p;
							}
						} break;
						case 16: {
							uint32_t v;
							while (total_count--) {
								v = ((uint16_t*)images[xor_image]->channels[0]->data)[p];
								channel_totals[0] += v * v;
								v = ((uint16_t*)images[xor_image]->channels[1]->data)[p];
								channel_totals[1] += v * v;
								v = ((uint16_t*)images[xor_image]->channels[2]->data)[p];
								channel_totals[2] += v * v;
								++p;
							}
						} break;
					}
				} else {
					switch (images[xor_image]->bpp) {
						case 8: {
							while (total_count--) {
								channel_totals[0] += ((uint8_t*)images[xor_image]->channels[0]->data)[p];
								channel_totals[1] += ((uint8_t*)images[xor_image]->channels[1]->data)[p];
								channel_totals[2] += ((uint8_t*)images[xor_image]->channels[2]->data)[p];
								++p;
							}
						} break;
						case 16: {
							while (total_count--) {
								channel_totals[0] += ((uint16_t*)images[xor_image]->channels[0]->data)[p];
								channel_totals[1] += ((uint16_t*)images[xor_image]->channels[1]->data)[p];
								channel_totals[2] += ((uint16_t*)images[xor_image]->channels[2]->data)[p];
								++p;
							}
						} break;
					}
				}

				if (!seamload_filename) {
					RECORD(xor_image, min_count);
					while (x < stop) {
						this_line[x++] = xor_image;
					}
				} else {
					x = stop;
				}

				best = xor_image;
			} else {
				if (xor_map) memset(&xor_map->line[x], 0xff, min_count);

				if (!seamload_filename) {
					if (y == 0) {
						// top row
						while (x < stop) {
							best = this_line[x];

							if (x > 0) {
								utemp = this_line[x - 1] + 0x300000000;
								d = MASKVAL(utemp);

								if (d < best) best = d;
							}

							if (best & 0x8000000000000000 && xor_count) {
								arbitrary_seam = true;
								for (i = 0; i < n_images; ++i) {
									if (!images[i]->mask_state) {
										best = 0x8000000000000000 | i;
										if (!reverse) break;
									}
								}
							}

							best_temp = best & 0xffffffff;
							RECORD(best_temp, 1);
							this_line[x++] = best;
						}
					} else {
						// other rows
						if (x == 0) {
							SEAM_DT;
							best = dt_val;

							utemp = *prev_line + 0x300000000;
							b = MASKVAL(utemp);
							if (b < best) best = b;

							utemp = prev_line[1] + 0x400000000;
							c = MASKVAL(utemp);
							if (c < best) best = c;

							if (best & 0x8000000000000000 && xor_count) {
								arbitrary_seam = true;
								for (i = 0; i < n_images; ++i) {
									if (!images[i]->mask_state) {
										best = 0x8000000000000000 | i;
										if (!reverse) break;
									}
								}
							}

							best_temp = best & 0xffffffff;
							RECORD(best_temp, 1);
							this_line[x++] = best;

							if (x == stop) {
								for (i = 0; i < n_images; ++i) {
									images[i]->mask_count -= min_count;
								}
								continue;
							}

							a = b + 0x100000000;
							b = c - 0x100000000;
						} else {
							utemp = prev_line[x - 1] + 0x400000000;
							a = MASKVAL(utemp);

							utemp = prev_line[x] + 0x300000000;
							b = MASKVAL(utemp);
						}

						utemp = best + 0x300000000;
						d = MASKVAL(utemp);

						if (stop == width) {
							stop--;
							last_pixel = true;
						}

						while (x < stop) {
							utemp = prev_line[x + 1] + 0x400000000;
							c = MASKVAL(utemp);

							SEAM_DT;
							best = dt_val;

							if (a < best) best = a;
							if (b < best) best = b;
							if (c < best) best = c;
							if (d < best) best = d;

							if (best & 0x8000000000000000 && xor_count) {
								arbitrary_seam = true;
								for (i = 0; i < n_images; ++i) {
									if (!images[i]->mask_state) {
										best = 0x8000000000000000 | i;
										if (!reverse) break;
									}
								}
							}

							best_temp = best & 0xffffffff;
							RECORD(best_temp, 1);
							this_line[x++] = best; // best;

							a = b + 0x100000000;
							b = c - 0x100000000;
							d = best + 0x300000000;
						}

						if (last_pixel) {
							SEAM_DT;
							best = dt_val;

							if (a < best) best = a;
							if (b < best) best = b;
							if (d < best) best = d;

							if (best & 0x8000000000000000 && xor_count) {
								arbitrary_seam = true;
								for (i = 0; i < n_images; ++i) {
									if (!images[i]->mask_state) {
										best = 0x8000000000000000 | i;
										if (!reverse) break;
									}
								}
							}

							best_temp = best & 0xffffffff;
							RECORD(best_temp, 1);
							this_line[x++] = best; // best;

							last_pixel = false;
						}
					}
				} else { // if (seamload_filename)...
					x = stop;
				}
			}

			for (i = 0; i < n_images; ++i) {
				images[i]->mask_count -= min_count;
			}
		}
			
		if (!seamload_filename) {
			RECORD(-1, 0);

			for (i = 0; i < n_images; ++i) {
				images[i]->masks[0]->NextLine();
			}
		}

		full_mask.NextLine();
		xor_mask.NextLine();

		if (xor_map) xor_map->Write();
		if (seam_map) seam_map->Write();

		std::swap(this_line, prev_line);
	}

	if (!seamload_filename) {
		delete[] thread_lines[0];
		delete[] thread_lines[1];
		delete[] thread_lines;
	}

	delete xor_map;
	delete seam_map;

	if (!alpha || output_type == ImageType::MB_JPEG) no_mask = true;

/***********************************************************************
* Seam load
***********************************************************************/
	if (seamload_filename) {
		int png_depth, png_colour;
		png_uint_32 png_width, png_height;
		uint8_t sig[8];
		png_structp png_ptr;
		png_infop info_ptr;
		FILE* f;

		fopen_s(&f, seamload_filename, "rb");
		if (!f) die("Error: Couldn't open seam file");

		size_t r = fread(sig, 1, 8, f); // assignment suppresses g++ -Ofast warning
		if (!png_check_sig(sig, 8)) die("Error: Bad PNG signature");

		png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
		if (!png_ptr) die("Error: Seam PNG problem");
		info_ptr = png_create_info_struct(png_ptr);
		if (!info_ptr) die("Error: Seam PNG problem");

		png_init_io(png_ptr, f);
		png_set_sig_bytes(png_ptr, 8);
		png_read_info(png_ptr, info_ptr);
		png_get_IHDR(png_ptr, info_ptr, &png_width, &png_height, &png_depth, &png_colour, NULL, NULL, NULL);

		if (png_width != width || png_height != png_height) die("Error: Seam PNG dimensions don't match workspace");
		if (png_depth != 8 || png_colour != PNG_COLOR_TYPE_PALETTE) die("Error: Incorrect seam PNG format");

		png_bytep png_line = (png_bytep)malloc(width);

		for (y = 0; y < height; ++y) {
			png_read_row(png_ptr, png_line, NULL);

			int ms = 0;
			int mc = 0;
			int prev_i = -1;
			int current_i = -1;

			for (x = 0; x < width; ++x) {
				if (png_line[x] > n_images) die("Error: Bad pixel found in seam file: %d,%d", x, y);
				RECORD(png_line[x], 1);
			}

			RECORD(-1, 0);

			for (i = 0; i < n_images; ++i) {
				images[i]->masks[0]->NextLine();
			}

		}

		free(png_line);
	}

	seam_time = timer.Read();

/***********************************************************************
* No output?
***********************************************************************/
	void* output_channels[3] = { NULL, NULL, NULL };

	if (output_type != ImageType::MB_NONE) {

/***********************************************************************
* Shrink masks
***********************************************************************/
		Output(1, "Shrinking masks...\n");

		timer.Start();

		for (i = 0; i < n_images; ++i) {
			threadpool->Queue([=] {
				ShrinkMasks(images[i]->masks, blend_levels);
			});
		}
		threadpool->Wait();

		shrink_mask_time = timer.Read();

/***********************************************************************
* Create shared input pyramids
***********************************************************************/
// wrapping
		std::vector<PyramidWithMasks*> wrap_pyramids;
		int wrap_levels_h = 0;
		int wrap_levels_v = 0;

		if (wrap & 1) {
			wrap_levels_h = (int)floor(log2((width >> 1) + 4.0f) - 1);
			wrap_pyramids.push_back(new PyramidWithMasks(width >> 1, height, wrap_levels_h, 0, 0, true));
			wrap_pyramids.push_back(new PyramidWithMasks((width + 1) >> 1, height, wrap_levels_h, width >> 1, 0, true));
		}

		if (wrap & 2) {
			wrap_levels_v = (int)floor(log2((height >> 1) + 4.0f) - 1);
			wrap_pyramids.push_back(new PyramidWithMasks(width, height >> 1, wrap_levels_v, 0, 0, true));
			wrap_pyramids.push_back(new PyramidWithMasks(width, (height + 1) >> 1, wrap_levels_v, 0, height >> 1, true));
		}

// masks
		for (auto& py : wrap_pyramids) {
			threadpool->Queue([=] {
				py->masks.push_back(new Flex(width, height));
				for (int y = 0; y < height; ++y) {
					if (y < py->GetY() || y >= py->GetY() + py->GetHeight()) {
						py->masks[0]->Write32(0x80000000 | width);
					} else {
						if (py->GetX()) {
							py->masks[0]->Write32(0x80000000 | py->GetX());
							py->masks[0]->Write32(0xc0000000 | py->GetWidth());
						} else {
							py->masks[0]->Write32(0xc0000000 | py->GetWidth());
							if (py->GetWidth() != width) py->masks[0]->Write32(0x80000000 | (width - py->GetWidth()));
						}
					}
					py->masks[0]->NextLine();
				}

				ShrinkMasks(py->masks, py->GetWidth() == width ? wrap_levels_v : wrap_levels_h);
				});
		}

		threadpool->Wait();
// end wrapping

		int total_levels = std::max({ blend_levels, wrap_levels_h, wrap_levels_v, 1 });

		for (int i = 0; i < n_images; ++i) {
			images[i]->pyramid = new Pyramid(images[i]->width, images[i]->height, blend_levels, images[i]->xpos, images[i]->ypos, true);
		}

		for (int l = total_levels - 1; l >= 0; --l) {
			size_t max_bytes = 0;

			if (l < blend_levels) {
				for (auto& image : images) {
					max_bytes = std::max(max_bytes, image->pyramid->GetLevel(l).bytes);
				}
			}

			for (auto& py : wrap_pyramids) {
				if (l < py->GetNLevels()) max_bytes = std::max(max_bytes, py->GetLevel(l).bytes);
			}

			float* temp;

			try {
				temp = (float*)MapAlloc::Alloc(max_bytes);
			} catch (char* e) {
				printf("%s\n", e);
				exit(EXIT_FAILURE);
			}

			if (l < blend_levels) {
				for (auto& image : images) {
					image->pyramid->GetLevel(l).data = temp;
				}
			}

			for (auto& py : wrap_pyramids) {
				if (l < py->GetNLevels()) py->GetLevel(l).data = temp;
			}
		}

/***********************************************************************
* Create output pyramid
***********************************************************************/
		Pyramid* output_pyramid = NULL;

		output_pyramid = new Pyramid(width, height, total_levels, 0, 0, true);

		for (int l = total_levels - 1; l >= 0; --l) {
			float* temp;

			try {
				temp = (float*)MapAlloc::Alloc(output_pyramid->GetLevel(l).bytes);
			} catch (char* e) {
				printf("%s\n", e);
				exit(EXIT_FAILURE);
			}

			output_pyramid->GetLevel(l).data = temp;
		}

/***********************************************************************
* Blend
***********************************************************************/
		if (n_images == 1) {
			if (wrap) Output(1, "Wrapping...\n"); else Output(1, "Processing...\n");
		} else {
			if (wrap) Output(1, "Blending/wrapping...\n"); else Output(1, "Blending...\n");
		}

		for (int c = 0; c < 3; ++c) {
			if (n_images > 1) {
				for (i = 0; i < n_images; ++i) {
					timer.Start();

					images[i]->pyramid->Copy((uint8_t*)images[i]->channels[c]->data, 1, images[i]->width, gamma, images[i]->bpp);
					if (output_bpp != images[i]->bpp) images[i]->pyramid->Multiply(0, gamma ? (output_bpp == 8 ? 1.0f / 66049 : 66049) : (output_bpp == 8 ? 1.0f / 257 : 257));

					delete images[i]->channels[c];
					images[i]->channels[c] = NULL;

					copy_time += timer.Read();

					timer.Start();
					images[i]->pyramid->Shrink();
					shrink_time += timer.Read();

					timer.Start();
					images[i]->pyramid->Laplace();
					laplace_time += timer.Read();

		// blend into output pyramid...

					timer.Start();

					for (int l = 0; l < blend_levels; ++l) {
						auto in_level = images[i]->pyramid->GetLevel(l);
						auto out_level = output_pyramid->GetLevel(l);

						int x_offset = (in_level.x - out_level.x) >> l;
						int y_offset = (in_level.y - out_level.y) >> l;

						for (int b = 0; b < (int)out_level.bands.size() - 1; ++b) {
							int sy = out_level.bands[b];
							int ey = out_level.bands[b + 1];

							threadpool->Queue([=] {
								for (int y = sy; y < ey; ++y) {
									int in_line = y - y_offset;
									if (in_line < 0) in_line = 0; else if (in_line > in_level.height - 1) in_line = in_level.height - 1;
									float* input_p = in_level.data + (size_t)in_line * in_level.pitch;
									float* output_p = out_level.data + (size_t)y * out_level.pitch;

									CompositeLine(input_p, output_p, i, x_offset, in_level.width, out_level.width, out_level.pitch, images[i]->masks[l]->data, images[i]->masks[l]->rows[y]);
								}
								});
						}

						threadpool->Wait();
					}

					blend_time += timer.Read();
				}

				timer.Start();
				output_pyramid->Collapse(blend_levels);
				collapse_time += timer.Read();
			} else {
				timer.Start();

				output_pyramid->Copy((uint8_t*)images[0]->channels[c]->data, 1, images[0]->width, gamma, images[0]->bpp);
				if (output_bpp != images[0]->bpp) output_pyramid->Multiply(0, gamma ? (output_bpp == 8 ? 1.0f / 66049 : 66049) : (output_bpp == 8 ? 1.0f / 257 : 257));

				delete images[0]->channels[c];
				images[0]->channels[c] = NULL;

				copy_time += timer.Read();
			}

/***********************************************************************
* Wrapping
***********************************************************************/
			if (wrap) {
				timer.Start();

				int p = 0;

				for (int w = 1; w <= 2; ++w) {
					if (wrap & w) {

						if (w == 1) {
							SwapH(output_pyramid);
						} else {
							SwapV(output_pyramid);
						}

						int wrap_levels = (w == 1) ? wrap_levels_h : wrap_levels_v;
						for (int wp = 0; wp < 2; ++wp) {
							wrap_pyramids[p]->Copy((uint8_t*)(output_pyramid->GetData() + wrap_pyramids[p]->GetX() + wrap_pyramids[p]->GetY() * (int64)output_pyramid->GetPitch()), 1, output_pyramid->GetPitch(), false, 32);
							wrap_pyramids[p]->Shrink();
							wrap_pyramids[p]->Laplace();

							for (int l = 0; l < wrap_levels; ++l) {
								auto in_level = wrap_pyramids[p]->GetLevel(l);
								auto out_level = output_pyramid->GetLevel(l);

								int x_offset = (in_level.x - out_level.x) >> l;
								int y_offset = (in_level.y - out_level.y) >> l;

								for (int b = 0; b < (int)out_level.bands.size() - 1; ++b) {
									int sy = out_level.bands[b];
									int ey = out_level.bands[b + 1];

									threadpool->Queue([=] {
										for (int y = sy; y < ey; ++y) {
											int in_line = y - y_offset;
											if (in_line < 0) in_line = 0; else if (in_line > in_level.height - 1) in_line = in_level.height - 1;
											float* input_p = in_level.data + (size_t)in_line * in_level.pitch;
											float* output_p = out_level.data + (size_t)y * out_level.pitch;

											CompositeLine(input_p, output_p, wp + (l == 0), x_offset, in_level.width, out_level.width, out_level.pitch, wrap_pyramids[p]->masks[l]->data, wrap_pyramids[p]->masks[l]->rows[y]);
										}
										});
								}

								threadpool->Wait();
							}
							++p;
						}

						output_pyramid->Collapse(wrap_levels);

						if (w == 1) {
							UnswapH(output_pyramid);
						} else {
							UnswapV(output_pyramid);
						}
					} // if (wrap & w)
				} // w loop

				wrap_time += timer.Read();
			} // if (wrap)
// end wrapping

/***********************************************************************
* Offset correction
***********************************************************************/
			if (total_pixels) {
				double channel_total = 0; // must be a double
				float* data = output_pyramid->GetData();
				xor_mask.Start();

				for (y = 0; y < height; ++y) {
					x = 0;
					while (x < width) {
						uint32_t v = xor_mask.ReadForwards32();
						if (v & 0x80000000) {
							v = x + v & 0x7fffffff;
							while (x < (int)v) {
								channel_total += data[x++];
							}
						} else {
							x += v;
						}
					}

					data += output_pyramid->GetPitch();
				}

				float avg = (float)channel_totals[c] / total_pixels;
				if (output_bpp != images[0]->bpp) {
					switch (output_bpp) {
						case 8: avg /= 256; break;
						case 16: avg *= 256; break;
					}
				}
				float output_avg = (float)channel_total / total_pixels;
				output_pyramid->Add(avg - output_avg, 1);
			}

/***********************************************************************
* Output
***********************************************************************/
			timer.Start();

			try {
				output_channels[c] = MapAlloc::Alloc(((size_t)width * height) << (output_bpp >> 4));
			} catch (char* e) {
				printf("%s\n", e);
				exit(EXIT_FAILURE);
			}

			switch (output_bpp) {
				case 8: output_pyramid->Out((uint8_t*)output_channels[c], width, gamma, dither, true); break;
				case 16: output_pyramid->Out((uint16_t*)output_channels[c], width, gamma, dither, true); break;
			}

			out_time += timer.Read();
		}

/***********************************************************************
* Write
***********************************************************************/
#define ROWS_PER_STRIP 64

		Output(1, "Writing %s...\n", output_filename);

		timer.Start();

		struct jpeg_compress_struct cinfo;
		struct jpeg_error_mgr jerr;

		JSAMPARRAY scanlines = NULL;

		int spp = no_mask ? 3 : 4;

		int bytes_per_pixel = spp << (output_bpp >> 4);
		int bytes_per_row = bytes_per_pixel * width;

		int n_strips = (int)((height + ROWS_PER_STRIP - 1) / ROWS_PER_STRIP);
		int remaining = height;
		void* strip = malloc((ROWS_PER_STRIP * (int64)width) * bytes_per_pixel);
		void* oc_p[3] = { output_channels[0], output_channels[1], output_channels[2] };
		if (bgr) std::swap(oc_p[0], oc_p[2]);

		switch (output_type) {
			case ImageType::MB_TIFF: {
				TIFFSetField(tiff_file, TIFFTAG_IMAGEWIDTH, width);
				TIFFSetField(tiff_file, TIFFTAG_IMAGELENGTH, height);
				TIFFSetField(tiff_file, TIFFTAG_COMPRESSION, compression);
				TIFFSetField(tiff_file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
				TIFFSetField(tiff_file, TIFFTAG_ROWSPERSTRIP, ROWS_PER_STRIP);
				TIFFSetField(tiff_file, TIFFTAG_BITSPERSAMPLE, output_bpp);
				if (no_mask) {
					TIFFSetField(tiff_file, TIFFTAG_SAMPLESPERPIXEL, 3);
				} else {
					TIFFSetField(tiff_file, TIFFTAG_SAMPLESPERPIXEL, 4);
					uint16_t out[1] = { EXTRASAMPLE_UNASSALPHA };
					TIFFSetField(tiff_file, TIFFTAG_EXTRASAMPLES, 1, &out);
				}

				TIFFSetField(tiff_file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
				if (images[0]->tiff_xres != -1) { TIFFSetField(tiff_file, TIFFTAG_XRESOLUTION, images[0]->tiff_xres); TIFFSetField(tiff_file, TIFFTAG_XPOSITION, (float)(min_xpos / images[0]->tiff_xres)); }
				if (images[0]->tiff_yres != -1) { TIFFSetField(tiff_file, TIFFTAG_YRESOLUTION, images[0]->tiff_yres); TIFFSetField(tiff_file, TIFFTAG_YPOSITION, (float)(min_ypos / images[0]->tiff_yres)); }

				if (images[0]->geotiff.set) {
					// if we got a georeferenced input, store the geotags in the output
					GeoTIFFInfo info(images[0]->geotiff);
					info.XGeoRef = min_xpos * images[0]->geotiff.XCellRes;
					info.YGeoRef = -min_ypos * images[0]->geotiff.YCellRes;
					Output(1, "Output georef: UL: %f %f, pixel size: %f %f\n", info.XGeoRef, info.YGeoRef, info.XCellRes, info.YCellRes);
					geotiff_write(tiff_file, &info);
				}
			} break;
			case ImageType::MB_JPEG: {
				cinfo.err = jpeg_std_error(&jerr);
				jpeg_create_compress(&cinfo);
				jpeg_stdio_dest(&cinfo, jpeg_file);

				cinfo.image_width = width;
				cinfo.image_height = height;
				cinfo.input_components = 3;
				cinfo.in_color_space = JCS_RGB;

				jpeg_set_defaults(&cinfo);
				jpeg_set_quality(&cinfo, jpeg_quality, true);
				jpeg_start_compress(&cinfo, true);
			} break;
			case ImageType::MB_PNG: {
				png_file = new Pnger(output_filename, NULL, width, height, no_mask ? PNG_COLOR_TYPE_RGB : PNG_COLOR_TYPE_RGB_ALPHA, output_bpp, jpeg_file, jpeg_quality);
			} break;
		}

		if (output_type == ImageType::MB_PNG || output_type == ImageType::MB_JPEG) {
			scanlines = new JSAMPROW[ROWS_PER_STRIP];
			for (i = 0; i < ROWS_PER_STRIP; ++i) {
				scanlines[i] = (JSAMPROW) & ((uint8_t*)strip)[i * bytes_per_row];
			}
		}

		full_mask.Start();

		for (int s = 0; s < n_strips; ++s) {
			int strip_p = 0;
			int rows = std::min(remaining, ROWS_PER_STRIP);

			for (int strip_y = 0; strip_y < rows; ++strip_y) {
				x = 0;
				while (x < width) {
					uint32_t cur = full_mask.ReadForwards32();
					if (cur & 0x80000000) {
						int lim = x + (cur & 0x7fffffff);
						switch (output_bpp) {
							case 8: {
								while (x < lim) {
									((uint8_t*)strip)[strip_p++] = ((uint8_t*)(oc_p[0]))[x];
									((uint8_t*)strip)[strip_p++] = ((uint8_t*)(oc_p[1]))[x];
									((uint8_t*)strip)[strip_p++] = ((uint8_t*)(oc_p[2]))[x];
									if (!no_mask) ((uint8_t*)strip)[strip_p++] = 0xff;
									++x;
								}
							} break;
							case 16: {
								while (x < lim) {
									((uint16_t*)strip)[strip_p++] = ((uint16_t*)(oc_p[0]))[x];
									((uint16_t*)strip)[strip_p++] = ((uint16_t*)(oc_p[1]))[x];
									((uint16_t*)strip)[strip_p++] = ((uint16_t*)(oc_p[2]))[x];
									if (!no_mask) ((uint16_t*)strip)[strip_p++] = 0xffff;
									++x;
								}
							} break;
						}
					} else {
						size_t t = (size_t)cur * bytes_per_pixel;
						switch (output_bpp) {
							case 8: {
								ZeroMemory(&((uint8_t*)strip)[strip_p], t);
							} break;
							case 16: {
								ZeroMemory(&((uint16_t*)strip)[strip_p], t);
							} break;
						}
						strip_p += cur * spp;
						x += cur;
					}
				}

				switch (output_bpp) {
					case 8: {
						oc_p[0] = &((uint8_t*)(oc_p[0]))[width];
						oc_p[1] = &((uint8_t*)(oc_p[1]))[width];
						oc_p[2] = &((uint8_t*)(oc_p[2]))[width];
					} break;
					case 16: {
						oc_p[0] = &((uint16_t*)(oc_p[0]))[width];
						oc_p[1] = &((uint16_t*)(oc_p[1]))[width];
						oc_p[2] = &((uint16_t*)(oc_p[2]))[width];
					} break;
				}
			}

			switch (output_type) {
				case ImageType::MB_TIFF: {
					TIFFWriteEncodedStrip(tiff_file, s, strip, rows * (int64)bytes_per_row);
				} break;
				case ImageType::MB_JPEG: {
					jpeg_write_scanlines(&cinfo, scanlines, rows);
				} break;
				case ImageType::MB_PNG: {
					png_file->WriteRows(scanlines, rows);
				} break;
			}

			remaining -= ROWS_PER_STRIP;
		}

		switch (output_type) {
			case ImageType::MB_TIFF: {
				TIFFClose(tiff_file);
			} break;
			case ImageType::MB_JPEG: {
				jpeg_finish_compress(&cinfo);
				jpeg_destroy_compress(&cinfo);
				fclose(jpeg_file);
			} break;
		}

		write_time = timer.Read();
	}

/***********************************************************************
* Timing
***********************************************************************/
	if (timing) {
		printf("\n");
		printf("Images:   %.3fs\n", images_time);
		printf("Seaming:  %.3fs\n", seam_time);
		if (output_type != ImageType::MB_NONE) {
			printf("Masks:    %.3fs\n", shrink_mask_time);
			printf("Copy:     %.3fs\n", copy_time);
			printf("Shrink:   %.3fs\n", shrink_time);
			printf("Laplace:  %.3fs\n", laplace_time);
			printf("Blend:    %.3fs\n", blend_time);
			printf("Collapse: %.3fs\n", collapse_time);
			if (wrap) printf("Wrapping: %.3fs\n", wrap_time);
			printf("Output:   %.3fs\n", out_time);
			printf("Write:    %.3fs\n", write_time);
		}
	}

/***********************************************************************
* Clean up
***********************************************************************/
	if (timing) {
		if (output_type == ImageType::MB_NONE) {
			timer_all.Report("\nExecution complete. Total execution time");
		} else {
			timer_all.Report("\nBlend complete. Total execution time");
		}
	}

	return EXIT_SUCCESS;
}
