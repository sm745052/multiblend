int hist_red[256];
int hist_grn[256];
int hist_blu[256];

class Channel {
public:
	Channel(size_t _bytes) : bytes(_bytes) {
		data = MapAlloc::Alloc(bytes);
	};
	~Channel() { MapAlloc::Free(data); };
	void* data;
	size_t bytes;
	FILE* file = NULL;
};

class Image {
public:
	Image(char* _filename);
	~Image();
	char* filename;
	ImageType type;
	int width;
	int height;
	int xpos;
	int ypos;
	int xpos_add = 0;
	int ypos_add = 0;
	std::vector<Channel*> channels;
	Pyramid* pyramid = NULL;
	GeoTIFFInfo geotiff;
	int tiff_width;
	int tiff_height;
	int tiff_u_height;
	int rows_per_strip;
	int first_strip;
	int end_strip;
	uint16_t bpp;
	uint16_t spp;
	void Open();
	void Read(void* data, bool gamma);
//	size_t untrimmed_pixels;
	size_t untrimmed_bytes;
	Flex* tiff_mask;
	float tiff_xres, tiff_yres;
	uint64_t mask_state;
	int mask_count;
	int mask_limit;
	bool seam_present;
	std::vector<Flex*> masks;
	void MaskPng(int i);

private:
	TIFF* tiff;
	FILE* file;
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	png_structp png_ptr;
};

Image::Image(char* _filename) : filename(_filename) {
}

Image::~Image() {
	for (auto it = channels.begin(); it < channels.end(); ++it) delete (*it);
	for (auto it = masks.begin(); it < masks.end(); ++it) delete (*it);
	channels.clear();
	masks.clear();

	delete pyramid;
}

/***********************************************************************
************************************************************************
** Open
************************************************************************
***********************************************************************/
void Image::Open() {
	float tiff_xpos, tiff_ypos;
	uint16_t compression;

	char* ext = strrchr(filename, '.');
	if (!ext) {
		die("Could not identify file extension: %s", filename);
	}
	++ext;

	if (!_stricmp(ext, "tif") || !_stricmp(ext, "tiff")) {
		type = ImageType::MB_TIFF;
	} else if (!_stricmp(ext, "jpg") || !_stricmp(ext, "jpeg")) {
		type = ImageType::MB_JPEG;
	} else if (!_stricmp(ext, "png")) {
		type = ImageType::MB_PNG;
	} else {
		die("Unknown file extension: %s", filename);
	}

	switch (type) {
		case ImageType::MB_TIFF: {
			tiff = TIFFOpen(filename, "r");
			if (!tiff) die("Could not open %s", filename);

			if (!TIFFGetField(tiff, TIFFTAG_XPOSITION, &tiff_xpos)) tiff_xpos = -1;
			if (!TIFFGetField(tiff, TIFFTAG_YPOSITION, &tiff_ypos)) tiff_ypos = -1;
			if (!TIFFGetField(tiff, TIFFTAG_XRESOLUTION, &tiff_xres)) tiff_xres = -1;
			if (!TIFFGetField(tiff, TIFFTAG_YRESOLUTION, &tiff_yres)) tiff_yres = -1;
			TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &tiff_width);
			TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &tiff_height);
			TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bpp);
			TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &spp);
			TIFFGetField(tiff, TIFFTAG_ROWSPERSTRIP, &rows_per_strip);
			TIFFGetField(tiff, TIFFTAG_COMPRESSION, &compression);

			if (bpp != 8 && bpp != 16) {
				printf("Invalid bpp %d (%s)", bpp, filename);
				printf("%d, %d\n", tiff_width, tiff_height);
				TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &bpp);
				if (bpp != 8 && bpp != 16) die("Invalid bpp %d (%s)", bpp, filename);
			}
//			if (spp != 4) die("Images must be RGBA (%s)", filename);

			geotiff.set = false;

			if (tiff_xpos == -1 && tiff_ypos == -1) {
				// try to read geotiff tags
				if (geotiff_read(tiff, &geotiff)) {
					xpos = (int)(geotiff.XGeoRef / geotiff.XCellRes);
					ypos = (int)(geotiff.YGeoRef / geotiff.YCellRes);
				} else {
					xpos = 0;
					ypos = 0;
				}
			} else {
				if (tiff_xpos != -1 && tiff_xres > 0) xpos = (int)(tiff_xpos * tiff_xres + 0.5);
				if (tiff_ypos != -1 && tiff_yres > 0) ypos = (int)(tiff_ypos * tiff_yres + 0.5);
			}

			first_strip = 0;
			end_strip = TIFFNumberOfStrips(tiff);

			int s;
			tmsize_t temp;

			if (tiff_xpos <= 0 && tiff_ypos <= 0 && compression != 1 && spp == 4) {
				tmsize_t min_stripsize = 0xffffffff;
				int min_stripcount = 0;
				for (s = 0; s < (int)TIFFNumberOfStrips(tiff) - 1; ++s) {
					temp = TIFFRawStripSize(tiff, s);
					if (temp < min_stripsize) { min_stripsize = temp; min_stripcount = 1; } else if (temp == min_stripsize) min_stripcount++;
				}

				if (min_stripcount > 2) {
					first_strip = -1;
					for (s = 0; s < (int)TIFFNumberOfStrips(tiff); ++s) {
						temp = TIFFRawStripSize(tiff, s);
						if (temp != min_stripsize) {
							if (first_strip == -1) first_strip = s;
							end_strip = s + 1;
						}
					}
					if (first_strip == -1) first_strip = 0;
				}
			}

			if (first_strip || end_strip != TIFFNumberOfStrips(tiff)) { // double check that min strips are (probably) transparent
				tdata_t buf = _TIFFmalloc(TIFFScanlineSize(tiff));
				if (first_strip) TIFFReadScanline(tiff, buf, 0); else TIFFReadScanline(tiff, buf, tiff_u_height - 1);
				bool trans;
				switch (bpp) {
					case 8: trans = !(((uint32_t*)buf)[0] && 0xff000000); break;
					case 16: trans = !(((uint64_t*)buf)[0] && 0xffff000000000000); break;
				}
				if (!trans) {
					first_strip = 0;
					end_strip = TIFFNumberOfStrips(tiff);
				}
				_TIFFfree(buf);
			}

			ypos += first_strip * rows_per_strip;
			int rows_missing = TIFFNumberOfStrips(tiff) * rows_per_strip - tiff_height;

			tiff_u_height = (end_strip - first_strip) * rows_per_strip;
			if (end_strip == TIFFNumberOfStrips(tiff)) tiff_u_height -= rows_missing;
		} break;
		case ImageType::MB_JPEG: {
			fopen_s(&file, filename, "rb");
			if (!file) die("Could not open %s", filename);

			cinfo.err = jpeg_std_error(&jerr);
			jpeg_create_decompress(&cinfo);
			jpeg_stdio_src(&cinfo, file);
			jpeg_read_header(&cinfo, TRUE);
			jpeg_start_decompress(&cinfo);

			if (!cinfo.output_width || !cinfo.output_height) {
				die("Unknown JPEG format (%s)", filename);
			}

			if (cinfo.out_color_components != 3) {
				die("Unknown JPEG format (%s)", filename);
			}

			tiff_width = cinfo.output_width;
			tiff_height = tiff_u_height = cinfo.output_height;

			bpp = 8;
			spp = 3;

			xpos = ypos = 0;
			tiff_xpos = tiff_ypos = 0;
			tiff_xres = tiff_yres = 90;
		} break;
		case ImageType::MB_PNG: {
			fopen_s(&file, filename, "rb");
			if (!file) die("Could not open %s", filename);

			uint8_t sig[8];
			size_t r = fread(sig, 1, 8, file); // assignment suppresses g++ -Ofast warning
			if (!png_check_sig(sig, 8)) die("Bad PNG signature (%s)", filename);

			png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
			if (!png_ptr) die("Error: libpng problem");
			png_infop info_ptr = png_create_info_struct(png_ptr);
			if (!info_ptr) die("Error: libpng problem");

			png_init_io(png_ptr, file);
			png_set_sig_bytes(png_ptr, 8);
			png_read_info(png_ptr, info_ptr);

			int png_colour;
			uint32_t png_width, png_height;
			int _bpp;
			png_get_IHDR(png_ptr, info_ptr, &png_width, &png_height, &_bpp, &png_colour, NULL, NULL, NULL);
			bpp = _bpp;
			tiff_width = png_width;
			tiff_u_height = tiff_height = png_height;

			switch (png_colour) {
				case PNG_COLOR_TYPE_RGB: spp = 3; break;
				case PNG_COLOR_TYPE_RGBA: spp = 4; break;
				default: die("Bad PNG colour type (%s)", filename);
			}

			if (bpp != 8 && bpp != 16) {
				die("Bad bit depth (%s)", filename);
			}

			xpos = ypos = 0;
			tiff_xpos = tiff_ypos = 0;
			tiff_xres = tiff_yres = 90;
		} break;
	}

	xpos += xpos_add;
	ypos += ypos_add;

	size_t untrimmed_pixels = (size_t)tiff_u_height * tiff_width;
	untrimmed_bytes = (untrimmed_pixels * spp ) << (bpp >> 4);
}

/***********************************************************************
************************************************************************
** Read
************************************************************************
***********************************************************************/
void Image::Read(void* data, bool gamma) {
	Output(1, "Processing %s...", filename);

	switch (type) {
		case ImageType::MB_TIFF: {
			char* pointer = (char*)data;

			for (int s = first_strip; s < end_strip; s++) {
				auto strip_size = TIFFReadEncodedStrip(tiff, s, pointer, -1);
				pointer += strip_size;
			}
		} break;
		case ImageType::MB_JPEG: {
			uint8_t* pointer = (uint8_t*)data;

			while (cinfo.output_scanline < cinfo.output_height) {
				jpeg_read_scanlines(&cinfo, &pointer, 1);
				pointer += (tiff_width * spp) << (bpp >> 4);
			}
		} break;
		case ImageType::MB_PNG: {
			uint8_t* pointer = (uint8_t*)data;

			for (int y = 0; y < tiff_height; ++y) {
				png_read_row(png_ptr, pointer, NULL);
				pointer += (tiff_width * spp) << (bpp >> 4);
			}
		} break;
	}

/***********************************************************************
* Trim
***********************************************************************/
	int x, y;
	int top, left, bottom, right;

	if (spp == 4) {
		switch (bpp) {
			case 8: {
				uint32_t* p = (uint32_t*)data;
				p--;

				for (y = 0; y < tiff_u_height; ++y) {
					for (x = 0; x < tiff_width; ++x) {
						if (*++p >= 0xff000000) {
							left = x;
							top = y;
							y = tiff_u_height;
							break;
						}
					}
				}

				p = ((uint32_t*)data) + tiff_width * tiff_u_height;
				for (y = tiff_u_height - 1; y >= 0; --y) {
					for (x = tiff_width - 1; x >= 0; --x) {
						if (*--p >= 0xff000000) {
							right = x;
							bottom = y;
							y = -1;
							break;
						}
					}
				}

				if (right < left) std::swap(left, right);

				p = ((uint32_t*)data) + tiff_width * top;
				for (y = top; y <= bottom; ++y) {
					for (x = 0; x < left; ++x) {
						if (p[x] >= 0xff000000) {
							left = x;
							break;
						}
					}

					for (x = tiff_width - 1; x > right; --x) {
						if (p[x] >= 0xff000000) {
							right = x;
							break;
						}
					}

					p += tiff_width;
				}
			} break;
			case 16: {
				uint16_t* p = ((uint16_t*)data) + 3;

				for (y = 0; y < tiff_u_height; ++y) {
					for (x = 0; x < tiff_width; ++x) {
						if (p[x << 2] == 0xffff) {
							left = x;
							top = y;
							y = tiff_u_height;
							break;
						}
					}

					p += tiff_width << 2;
				}

				p = ((uint16_t*)data) + (tiff_width << 2) * (tiff_u_height - 1) + 3;
				for (y = tiff_u_height - 1; y >= 0; --y) {
					for (x = tiff_width - 1; x >= 0; --x) {
						if (p[x << 2] == 0xffff) {
							right = x;
							bottom = y;
							y = -1;
							break;
						}
					}

					p -= tiff_width << 2;
				}

				if (right < left) std::swap(left, right);

				p = ((uint16_t*)data) + (tiff_width << 2) * top + 3;
				for (y = top; y <= bottom; ++y) {
					for (x = 0; x < left; ++x) {
						if (p[x << 2] == 0xffff) {
							left = x;
							break;
						}
					}

					for (x = tiff_width - 1; x > right; --x) {
						if (p[x << 2] == 0xffff) {
							right = x;
							break;
						}
					}

					p += tiff_width << 2;
				}

			} break;
		}

		width = right + 1 - left;
		height = bottom + 1 - top;

		xpos += left;
		ypos += top;

/***********************************************************************
* Inpaint
***********************************************************************/
		int temp_copy;
		uint32_t a, b, c, d;
		uint32_t* this_line = NULL;
		uint32_t* prev_line = NULL;
		Threadpool* threadpool = Threadpool::GetInstance();

		tiff_mask = new Flex(width, height);
		Flex* dt = new Flex(width, height);
		int mc;

		int n_threads = (std::max)(2, threadpool->GetNThreads());

		uint32_t** thread_lines = new uint32_t * [n_threads];
		uint32_t** thread_comp_lines = new uint32_t * [n_threads];

		std::mutex* flex_mutex_p = new std::mutex;
		std::condition_variable* flex_cond_p = new std::condition_variable;

		for (int i = 0; i < n_threads; ++i) {
			thread_lines[i] = new uint32_t[width];
			thread_comp_lines[i] = new uint32_t[width];
		}

		uint32_t* bitmap32 = NULL;
		uint64_t* bitmap64 = NULL;
		if (bpp == 8) bitmap32 = ((uint32_t*)data) + top * tiff_width + left; else bitmap64 = ((uint64_t*)data) + top * tiff_width + left;

		for (y = 0; y < height; ++y) {
			int t = y % n_threads;
			this_line = thread_lines[t];
			uint32_t* comp = thread_comp_lines[t];
			bool first;

			x = 0;

			{ // make sure the last compression thread to use this chunk of memory is finished
				std::unique_lock<std::mutex> mlock(*flex_mutex_p);
				flex_cond_p->wait(mlock, [=] { return dt->y > y - n_threads; });
			}

			while (x < width) {
				mc = 0;
				first = true;

				while (x < width && (bpp == 8 ? (bitmap32[x] < 0xff000000) : (bitmap64[x] < 0xffff000000000000))) {
					if (y == 0) {
						if (x == 0) {
							this_line[0] = 0x80000000;
						} else {
							this_line[x] = this_line[x - 1] + 3;
							if (bpp == 8) {
								bitmap32[x] = bitmap32[x - 1];
							} else {
								bitmap64[x] = bitmap64[x - 1];
							}
						}
					} else {
						if (x == 0) {
							b = prev_line[x] + 3;
							c = prev_line[x + 1] + 4;

							if (b < c) {
								d = b;
								temp_copy = x - tiff_width;
							} else {
								d = c;
								temp_copy = x - tiff_width + 1;
							}
							d = b < c ? b : c;

							this_line[x] = d;
							if (bpp == 8) {
								bitmap32[x] = bitmap32[temp_copy];
							} else {
								bitmap64[x] = bitmap64[temp_copy];
							}

							a = b + 1;
							b = c - 1;
							d += 3;
						} else {
							if (first) {
								a = prev_line[x - 1] + 4;
								b = prev_line[x] + 3;
								d = this_line[x - 1] + 3;
								first = false;
							}

							if (x < width - 1) c = prev_line[x + 1] + 4; else c = 0xffffffff;
							temp_copy = x - 1;

							if (a < d) { d = a; temp_copy = x - tiff_width - 1; }
							if (b < d) { d = b; temp_copy = x - tiff_width; }
							if (c < d) { d = c; temp_copy = x - tiff_width + 1; }

							this_line[x] = d;
							if (bpp == 8) {
								bitmap32[x] = bitmap32[temp_copy];
							} else {
								bitmap64[x] = bitmap64[temp_copy];
							}

							a = b + 1;
							b = c - 1;
							d += 3;
						}
					}
					++x;
					++mc;
				}

				if (mc) {
					tiff_mask->Write32(mc);
					mc = 0;
				}

				switch (bpp) {
					case 8: {
						while (x < width && (bitmap32[x] >= 0xff000000)) {
							this_line[x++] = 0;
							++mc;
						}
					} break;
					case 16: {
						while (x < width && (bitmap64[x] >= 0xffff000000000000)) {
							this_line[x++] = 0;
							++mc;
						}
//						}
					} break;
				}

				if (mc) {
					tiff_mask->Write32(0x80000000 | mc);
					mc = 0;
				}
			}

			if (bpp == 8) bitmap32 += tiff_width; else bitmap64 += tiff_width;

			if (y < height - 1) {
				threadpool->Queue([=] {
					int p = CompressDTLine(this_line, (uint8_t*)comp, width);
					{
						std::unique_lock<std::mutex> mlock(*flex_mutex_p);
						flex_cond_p->wait(mlock, [=] { return dt->y == y; });
						dt->Copy((uint8_t*)comp, p);
						dt->NextLine();
					}
					flex_cond_p->notify_all();
				});
			}

			tiff_mask->NextLine();
			prev_line = this_line;
		}

		threadpool->Wait();

		// backward
		int current_count = 0;
		int current_step;
		uint32_t dt_val;

		uint32_t mask;

		prev_line = thread_lines[(y - 2) % n_threads];

		// first line
		x = width - 1;
		if (bpp == 8) bitmap32 -= tiff_width; else bitmap64 -= tiff_width;

		while (x >= 0) {
			mask = tiff_mask->ReadBackwards32();
			if (mask & 0x80000000) { // solid
				x -= mask & 0x7fffffff;
				d = 3;
			} else {
				if (x == width - 1) {
					d = this_line[x] + 3;
					--mask;
					--x;
				}
				while (mask) {
					uint32_t best = this_line[x];
					if (d < best) {
						this_line[x] = d;
						if (bpp == 8) bitmap32[x] = bitmap32[x + 1]; else bitmap64[x] = bitmap64[x + 1];
						d += 3;
					} else {
						d = best + 3;
					}
					--x;
					--mask;
				}
			}
		}

		std::swap(this_line, prev_line);

		// other lines
		for (y = height - 2; y >= 0; --y) {
			x = width - 1;
			c = d = 0x80000000;
			if (bpp == 8) bitmap32 -= tiff_width; else bitmap64 -= tiff_width;

			while (x >= 0) {
				mask = tiff_mask->ReadBackwards32();

				if (mask & 0x80000000) { // solid
					mc = mask & 0x7fffffff;
					x -= mc;
					memset(&this_line[x + 1], 0, mc << 2);
					d = 3;
				} else {
					b = prev_line[x] + 3;
					c = (x < width - 1) ? prev_line[x + 1] + 4 : 0x80000000;
					while (mask) {
						a = (x > 0) ? prev_line[x - 1] + 4 : 0x80000000;
						INPAINT_DT;
						int copy = 0;
						uint32_t best = dt_val;
						if (a < best) { best = a; copy = tiff_width + x - 1; }
						if (b < best) { best = b; copy = tiff_width + x; }
						if (c < best) { best = c; copy = tiff_width + x + 1; }
						if (d < best) { best = d; copy = x + 1; }

						if (copy) {
							if (bpp == 8) bitmap32[x] = bitmap32[copy]; else bitmap64[x] = bitmap64[copy];
						}

						this_line[x--] = best;
						c = b + 1;
						d = best + 3;
						b = a - 1;
						mask--;
					}
				}
			}

			std::swap(this_line, prev_line);
		}

		delete flex_mutex_p;
		delete flex_cond_p;

		for (int i = 0; i < n_threads; ++i) {
			delete[] thread_lines[i];
			delete[] thread_comp_lines[i];
		}

		delete[] thread_lines;
		delete[] thread_comp_lines;
	} else {
		width = tiff_width;
		height = tiff_height;

		tiff_mask = new Flex(width, height);

		for (y = 0; y < height; ++y) {
			tiff_mask->Write32(0x80000000 | width);
			tiff_mask->NextLine();
		}
	}

/***********************************************************************
* Extract channels
***********************************************************************/
	size_t channel_bytes = ((size_t)width * height) << (bpp >> 4);

	for (int c = 0; c < 3; ++c) {
		channels.push_back(new Channel(channel_bytes));
	}

	if (spp == 4) {
		switch (bpp) {
			case 8: {
				uint32_t* line = ((uint32_t*)data) + top * tiff_width + left;
				int p = 0;
				for (y = 0; y < height; ++y) {
					for (x = 0; x < width; ++x) {
						uint32_t pixel = line[x];
						((uint8_t*)channels[0]->data)[p + x] = pixel & 0xff;
						((uint8_t*)channels[1]->data)[p + x] = (pixel >> 8) & 0xff;
						((uint8_t*)channels[2]->data)[p + x] = (pixel >> 16) & 0xff;
					}
					p += width;
					line += tiff_width;
				}
			} break;
			case 16: {
				uint64_t* line = ((uint64_t*)data) + top * tiff_width + left;
				int p = 0;
				for (y = 0; y < height; ++y) {
					for (x = 0; x < width; ++x) {
						uint64_t pixel = line[x];
						((uint16_t*)channels[0]->data)[p + x] = pixel & 0xffff;
						((uint16_t*)channels[1]->data)[p + x] = (pixel >> 16) & 0xffff;
						((uint16_t*)channels[2]->data)[p + x] = (pixel >> 32) & 0xffff;
					}
					p += width;
					line += tiff_width;
				}
			} break;
		}
	} else {
		switch (bpp) {
			case 8: {
				uint8_t byte;
				uint8_t* bytes = (uint8_t*)data;
				size_t p = 0;
				for (y = 0; y < height; ++y) {
					for (x = 0; x < width; ++x) {
						byte = *bytes++;
						((uint8_t*)channels[0]->data)[p] = byte;
//						channel_totals[0] += gamma ? byte*byte : byte;

						byte = *bytes++;
						((uint8_t*)channels[1]->data)[p] = byte;
//						channel_totals[1] += gamma ? byte * byte : byte;

						byte = *bytes++;
						((uint8_t*)channels[2]->data)[p] = byte;
//						channel_totals[2] += gamma ? byte * byte : byte;

						++p;
					}
				}
			} break;
			case 16: {
				uint16_t word;
				uint16_t* words = (uint16_t*)data;
				size_t p = 0;
				for (y = 0; y < height; ++y) {
					for (x = 0; x < width; ++x) {
						word = *words++;
						((uint16_t*)channels[0]->data)[p] = word;
//						channel_totals[0] += gamma ? word * word : word;

						word = *words++;
						((uint16_t*)channels[1]->data)[p] = word;
//						channel_totals[1] += gamma ? word * word : word;

						word = *words++;
						((uint16_t*)channels[2]->data)[p] = word;
//						channel_totals[2] += gamma ? word * word : word;

						++p;
					}
				}
			} break;
		}

//		total_pixels += width * height;
	}

	Output(1, "\n");
}

/***********************************************************************
* Debugging
***********************************************************************/
void Image::MaskPng(int i) {
	char filename[256];
	sprintf_s(filename, "masks-%d.png", i);

	int width = masks[0]->width;
	int height = masks[0]->height; // +1 + masks[1]->height;

	size_t size = (size_t)width * height;
	uint8_t* temp = (uint8_t*)malloc(size);
	memset(temp, 32, size);

	int px = 0, py = 0;
	uint32_t cur;

	for (int l = 0; l < (int)masks.size(); ++l) {
		uint32_t* data = (uint32_t*)masks[l]->data;
		uint8_t* line = temp + py * masks[0]->width + px;
		for (int y = 0; y < masks[l]->height; ++y) {
			int x = 0;
			while (x < masks[l]->width) {
				float val;
				int count;

				cur = *data++;

				if (cur & 0x80000000) {
					count = cur & 0x00ffffff;
					if (cur & 0x20000000) val = *(float*)data++; else val = (float)((cur >> 30) & 1);
				} else {
					val = *((float*)&cur);
					count = 1;
				}

				int t = x + count;
				while (x < t) line[x++] = (uint8_t)(val * 255);
			}

			line += masks[0]->width;
		}
		if (l & 1) px += masks[l]->width + 1; else py += masks[l]->height + 1;
		break;
	}

	Pnger::Quick(filename, temp, width, height, width, PNG_COLOR_TYPE_GRAY);
}
