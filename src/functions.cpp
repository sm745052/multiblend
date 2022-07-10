#include <chrono>

/***********************************************************************
* Flexible data class
***********************************************************************/
class Flex {
public:
	Flex(int _width, int _height) : width(_width), height(_height) {
		rows = new int[height];
		NextLine();
	};

	void NextLine() {
		end_p = p;
		if (y < height - 1) rows[++y] = p;

		if (p + (width << 2) > size) {
			if (y == 0) {
				size = (std::max)(height, 16) << 4; // was << 2
				data = (uint8_t*)malloc(size);
			} else if (y < height) {
				int prev_size = size;
				int new_size1 = (p / y) * height + (width << 4);
				int new_size2 = size << 1;
				size = (std::max)(new_size1, new_size2);
				data = (uint8_t*)realloc(data, size);
			}
		}

		MaskFinalise();
		first = true;
	}

	void Shrink() {
		data = (uint8_t*)realloc(data, p);
		end_p = p;
	}

	void Write32(uint32_t w)      { *((uint32_t*)&data[p]) = w; p += 4; }
	void Write64(uint64_t w)      { *((uint64_t*)&data[p]) = w; p += 8; }

	void MaskWrite(int count, bool white) {
		if (first) {
			mask_count = count;
			mask_white = white;
			first = false;
		} else {
			if (white == mask_white) {
				mask_count += count;
			} else {
				Write32((mask_white << 31) | mask_count);
				mask_count = count;
				mask_white = white;
			}
		}
	}

	void MaskFinalise() {
		if (mask_count) Write32((mask_white << 31) | mask_count);
	}

	void IncrementLast32(int inc) { *((uint32_t*)&data[p - 4]) += inc; }

	uint8_t  ReadBackwards8 () { return data[--p]; }
	uint16_t ReadBackwards16() { return *((uint16_t*)&data[p -= 2]); }
	uint32_t ReadBackwards32() { return *((uint32_t*)&data[p -= 4]); }
	uint64_t ReadBackwards64() { return *((uint64_t*)&data[p -= 8]); }

	uint32_t ReadForwards32() {
		uint32_t out = *((uint32_t*)&data[p]);
		p += 4;
		return out;
	}

	void Copy(uint8_t* src, int len) {
		memcpy(&data[p], src, len);
		p += len;
	}

	void Start() {
		p = 0;
	}

	void End() {
		p = end_p;
	}

	~Flex() {
		free(data);
		delete[] rows;
	}

	uint8_t* data = NULL;
	int width;
	int height;
	int* rows;
	int y = -1;

private:
	int size = 0;
	int p = 0;
	int end_p = 0;
	int mask_count = 0;
	bool mask_white;
	bool first;
};

/***********************************************************************
* Output
***********************************************************************/
void Output(int level, const char* fmt, ...) {
	va_list args;

	if (level <= verbosity) {
		va_start(args, fmt);
		vprintf(fmt, args);
		va_end(args);
	}
	fflush(stdout);
}

/***********************************************************************
* Timer
***********************************************************************/
class Timer {
public:
	void Start() {
		start_time = std::chrono::high_resolution_clock::now();
	};

	double Read() {
		std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start_time;
		return elapsed.count();
	};

	void Report(const char* name) {
		printf("%s: %.3fs\n", name, Read());
	};

private:
	std::chrono::high_resolution_clock::time_point start_time;
};

/***********************************************************************
* Die
***********************************************************************/
void die(const char* error, ...) {
	va_list args;

	va_start(args,error);
	vprintf(error,args);
	va_end(args);
	printf("\n");

  exit(EXIT_FAILURE);
}

/***********************************************************************
* ShrinkMasks
***********************************************************************/
int Squish(uint32_t* in, uint32_t* out, int in_width, int out_width) {
	float current_val;
	uint32_t cur;
	float a, b, c, d, e;
	int last_int = -1;

	int in_p = 0;
	int out_p = 0;
	int read = 0;
	int wrote = 0;
	int in_count;
	int out_count;

	cur = in[in_p++];
	if (cur & 0x80000000) {
		in_count = cur & 0x00ffffff;
		if (cur & 0x20000000) {
			current_val = ((float*)in)[in_p++];
		} else {
			last_int = (cur >> 30) & 1;
		}
	} else {
		in_count = 1;
		current_val = *((float*)&cur);
	}

	d = -1; e = -1;

	// trivial: full width block
	if (in_count == in_width) {
		if (last_int >= 0) {
			*out = (cur & 0xff000000) | out_width;
		} else {
			out[out_p++] = (cur & 0xa0000000) | out_width;
			((float*)out)[out_p] = current_val;
		}
		return in_p;
	}

	out_count = (in_count + 1) >> 1;
	if (last_int >= 0) {
		out[out_p++] = (cur & 0xff000000) | out_count;
		current_val = (float)last_int;
	} else {
		if (out_count > 1) out[out_p++] = (cur & 0xa0000000) | out_count;
		((float*)out)[out_p++] = current_val;
	}
	wrote = out_count;

	a = b = c = current_val;
	read = in_count;
	in_count -= (in_count - 1) | 1;

	int c_read = in_p;
	int e_read;

	while (wrote < out_width) {
		// loop always starts needing to read two more pixels

		// read d
		if (!in_count && read < in_width) {
			cur = in[in_p++];
			if (cur & 0x80000000) { // float
				in_count = cur & 0x00ffffff;
				read += in_count;
				if (cur & 0x20000000) { // repeated float
					last_int = -1;
					current_val = *((float*)&in[in_p++]);
				} else {
					last_int = (cur >> 30) & 1;
					current_val = (float)last_int;
				}
			} else {
				current_val = *((float*)&cur);
				last_int = -1;
				in_count = 1;
				++read;
			}
		}

		if (read == in_width) in_count = 0x7fffffff;

		if (in_count >= 2) {
			d = e = current_val;
			e_read = in_p;
			in_count -= 2;
		} else {
			d = current_val;

			cur = in[in_p++];
			if (cur & 0x80000000) {
				in_count = cur & 0x00ffffff;
				read += in_count;
				if (cur & 0x20000000) {
					last_int = -1;
					current_val = ((float*)in)[in_p++];
				} else {
					last_int = (cur >> 30) & 1;
					current_val = (float)last_int;
				}
			} else {
				current_val = *((float*)&cur);
				last_int = -1;
				in_count = 1;
				++read;
			}

			e = current_val;
			e_read = in_p;
			--in_count;
		}

		// calculate
		((float*)out)[out_p++] = (b + d) * 0.25f + (a + e + 6 * c) * 0.0625f;
		wrote++;

		if (c_read == in_p && in_count >= 4) {
			out_count = in_count >> 1;
			if (out_count > (out_width - wrote)) out_count = out_width - wrote;
			if (last_int >= 0) {
				out[out_p++] = ((0x2 | last_int) << 30) | out_count;
			} else {
				out[out_p++] = 0xa0000000 | out_count;
				((float*)out)[out_p++] = e;
			}
			wrote += out_count;
			in_count -= out_count << 1;
		}
		a = c;
		b = d;
		c = e;
		c_read = e_read;
	}

	return in_p;
}

void ShrinkMasks(std::vector<Flex*>& masks, int n_levels) {
	int i;
	Flex flex_temp(masks[0]->width, masks[0]->height);
	uint32_t cur;
	uint32_t* lines[5];
	uint32_t* real_lines[5];
	int count[5];
	int pointer[5];
	float vals[5];
	int min_count;

	for (int i = 0; i < 5; ++i) real_lines[i] = new uint32_t[masks[0]->width];

	for (int l = 1; l < n_levels; ++l) {
		int in_width = masks[l - 1]->width;
		int out_width = (in_width + 6) >> 1;
		int out_height = (masks[l - 1]->height + 6) >> 1;
		masks.push_back(new Flex(out_width, out_height));

		int input_p = 0;

		// first line
		input_p += Squish((uint32_t*)&masks[l - 1]->data[input_p], real_lines[0], in_width, out_width) << 2;
		int lines_read = 1;

		int wrote = 0;
		pointer[0] = 0;
		while (wrote < out_width) {
			cur = real_lines[0][pointer[0]++];
			masks[l]->Write32(cur);
			if (!(cur & 0x80000000)) {
				++wrote;
			} else {
				if (cur & 0x20000000) {
					masks[l]->Write32(real_lines[0][pointer[0]++]);
				}
				wrote += cur & 0x00ffffff;
			}
		}

		// other lines
		lines[0] = lines[1] = lines[2] = real_lines[0];
		lines[3] = real_lines[3];
		lines[4] = real_lines[4];

		input_p += Squish((uint32_t*)&masks[l - 1]->data[input_p], lines[3], in_width, out_width) << 2;
		input_p += Squish((uint32_t*)&masks[l - 1]->data[input_p], lines[4], in_width, out_width) << 2;
		lines_read += 2;

		for (int y = 1; y < masks[l]->height; ++y) {
			masks[l]->NextLine();

			wrote = 0;
			pointer[0] = pointer[1] = pointer[2] = pointer[3] = pointer[4] = 0;
			count[0] = count[1] = count[2] = count[3] = count[4] = 0;
			int last_int;

			while (wrote < out_width) {
				min_count = 0x7fffffff;
				for (i = 0; i < 5; ++i) {
					if (!count[i]) {
						cur = lines[i][pointer[i]++];
						if (!(cur & 0x80000000)) {
							vals[i] = *((float*)&cur);
							count[i] = 1;
						} else {
							if (cur & 0x20000000) {
								vals[i] = ((float*)lines[i])[pointer[i]++];
							} else {
								last_int = (cur >> 30) & 1;
								vals[i] = (float)last_int;
							}
							count[i] = cur & 0x00ffffff;
						}
					}

					if (count[i] < min_count) min_count = count[i];
				}

				float val = (vals[0] + vals[4]) * 0.0625f + (vals[1] + vals[3]) * 0.25f + vals[2] * 0.375f;
				if (val == 0 || val == 1) {
					if (min_count == 1) {
						masks[l]->Write32(*((uint32_t*)&val));
					} else {
						int _int = (int)val;
						masks[l]->Write32(0x80000000 | (_int << 30) | min_count);
					}
				} else {
					if (min_count > 1) {
						masks[l]->Write32(0xa0000000 | min_count);
					}
					masks[l]->Write32(*((uint32_t*)&val));
				}

				wrote += min_count;
				for (i = 0; i < 5; ++i) count[i] -= min_count;
			}

			if (y == 1) {
				lines[0] = real_lines[1];
				lines[1] = real_lines[2];
			}

			uint32_t* temp = lines[0];
			lines[0] = lines[2];
			lines[2] = lines[4];
			lines[4] = lines[1];
			lines[1] = lines[3];
			lines[3] = temp;

			if (lines_read < masks[l - 1]->height) {
				input_p += Squish((uint32_t*)&masks[l - 1]->data[input_p], lines[3], in_width, out_width) << 2;
				++lines_read;
			} else {
				lines[3] = lines[2];
			}

			if (lines_read < masks[l - 1]->height) {
				input_p += Squish((uint32_t*)&masks[l - 1]->data[input_p], lines[4], in_width, out_width) << 2;
				++lines_read;
			} else {
				lines[4] = lines[3];
			}
		}
	}

	for (i = 0; i < 5; ++i) delete[] real_lines[i];
}

/***********************************************************************
* Composite line
***********************************************************************/
void CompositeLine(float* input_p, float* output_p, int i, int x_offset, int in_level_width, int out_level_width, int out_level_pitch, uint8_t* _mask, size_t mask_p) {
	int x = 0;
	int mask_count;
	float current_val;
	int last_int;

	mask_p >>= 2;
	uint32_t* mask = (uint32_t*)_mask;

	while (x < out_level_width) {
		uint32_t cur = mask[mask_p++];
		last_int = -1;
		if (cur & 0x80000000) {
			mask_count = cur & 0x0fffffff;
			if (cur & 0x20000000) {
				current_val = ((float*)mask)[mask_p++];
			} else {
				last_int = (cur >> 30) & 1;
			}
		} else {
			mask_count = 1;
			current_val = *((float*)&cur);
		}

		int lim = x + mask_count;

		if (last_int == 0) {
			if (i == 0) memset(&output_p[x], 0, mask_count << 2);
			x += mask_count;
		} else if (last_int == 1) {
			if (x < x_offset) {
				float f = input_p[0];
				while (x < lim && x < x_offset) {
					output_p[x++] = f;
				}
			}

			while (x < lim && x < x_offset + in_level_width) {
				output_p[x] = input_p[x - x_offset];
				++x;
			}

			if (x < lim) {
				float f = input_p[in_level_width - 1];
				while (x < lim) {
					output_p[x++] = f;
				}
			}
		} else {
			if (x < x_offset) {
				float f = input_p[0] * current_val;
				while (x < lim && x < x_offset) {
					if (i == 0) output_p[x++] = f; else output_p[x++] += f;
				}
			}

			while (x < lim && x < x_offset + in_level_width) {
				float f = input_p[x - x_offset] * current_val;
				if (i == 0) output_p[x++] = f; else output_p[x++] += f;
			}

			if (x < lim) {
				float f = input_p[in_level_width - 1] * current_val;
				while (x < lim) {
					if (i == 0) output_p[x++] = f; else output_p[x++] += f;
				}
			}
		}
	}

	if (x < out_level_pitch) {
		float f = output_p[x - 1];
		do {
			output_p[x++] = f;
		} while (x < out_level_pitch);
	}
}

/***********************************************************************
* Seam macro
***********************************************************************/
// Record macro
#define RECORD(I, C)\
	if ((I) != current_i) {\
		if (mc > 0) {\
			if (seam_map) memset(&seam_map->line[x - mc], current_i, mc);\
			for (i = 0; i < n_images; ++i) {\
				if (i == current_i) {\
					images[i]->masks[0]->Write32(0xc0000000 | mc);\
				} else if (i == prev_i || prev_i == -1) {\
					images[i]->masks[0]->Write32(0x80000000 | mc);\
				} else {\
					images[i]->masks[0]->IncrementLast32(mc);\
				}\
			}\
		}\
		prev_i = current_i;\
		mc = C;\
		current_i = I;\
	} else {\
		mc += C;\
	}
// end macro

/***********************************************************************
* DT reader and macros
***********************************************************************/
void ReadInpaintDT(Flex* flex, int& current_count, int& current_step, uint32_t& dt_val) {
	if (current_count) {
		--current_count;
		dt_val += current_step;
	} else {
		uint8_t _byte = flex->ReadBackwards8();
		if (_byte == 255) {
			dt_val = flex->ReadBackwards32();
			return;
		} else {
			current_step = ((_byte & 7) - 3);
			if (!(_byte & 0x80)) { // 0b0ccccsss
				current_count = _byte >> 3;
			} else if (!(_byte & 0x40)) { // 0b10ssssss
				current_step = (_byte & 0x3f);
				current_count = 0;
			} else if (!(_byte & 0x20)) { // 0b11000000
				current_count = flex->ReadBackwards8();
			} else if (!(_byte & 0x10)) { // 0b11100000
				current_count = flex->ReadBackwards16();
			} else { //if (!(_byte & 0x08)) {
				current_count = flex->ReadBackwards32();
			}
			dt_val += current_step;
		}
	}
}

void ReadSeamDT(Flex* flex, int& current_count, int64& current_step, uint64_t& dt_val) {
	if (current_count) {
		--current_count;
		dt_val += current_step;
	} else {
		uint8_t _byte = flex->ReadBackwards8();
		if (_byte == 255) {
			dt_val = flex->ReadBackwards64();
			return;
		} else {
			current_step = ((int64)(_byte & 7) - 3) << 32;
			if (!(_byte & 0x80)) { // 0b0ccccsss
				current_count = _byte >> 3;
			} else if (!(_byte & 0x40)) { // 0b10ssssss
				current_step = (int64)(_byte & 0x3f) << 32;
				current_count = 0;
			} else if (!(_byte & 0x20)) { // 0b11000000
				current_count = flex->ReadBackwards8();
			} else if (!(_byte & 0x10)) { // 0b11100000
				current_count = flex->ReadBackwards16();
			} else { //if (!(_byte & 0x08)) {
				current_count = flex->ReadBackwards32();
			}
			dt_val += current_step;
		}
	}
}


#define SEAM_DT    ReadSeamDT(seam_flex, current_count, current_step, dt_val);
#define INPAINT_DT ReadInpaintDT(    dt, current_count, current_step, dt_val);

/***********************************************************************
* Seam line compress
***********************************************************************/
void cwrite(int current_count, int current_step, uint8_t * output, int& p) {
	if (current_count) {
		if (--current_count < 16) {
			output[p++] = current_count << 3 | current_step;
		} else if (current_count < 256) {
			output[p++] = current_count;
			output[p++] = 0xc0 | current_step;
		} else if (current_count < 65536) {
			*((uint16_t*)&output[p]) = current_count;
			p += 2;
			output[p++] = 0xe0 | current_step;
		} else {
			*((uint32_t*)&output[p]) = current_count;
			p += 4;
			output[p++] = 0xf0 | current_step;
		}
	}
}

#define CWRITE cwrite(current_count, current_step, output, p);

int CompressDTLine(uint32_t* input, uint8_t* output, int width) {
	int current_step = -100;
	int current_count = 0;

	int step;
	int x = 0;
	int p = 0;
	uint32_t left_val;
	uint32_t right_val;

	while (!(left_val = input[x++]) && x < width);
	if (!left_val) return p;

	while (x < width) {
		while (!(right_val = input[x++]) && x < width);
		if (!right_val) break;

		if ((step = (int)(left_val - right_val) + 3) < 67 && step >= 0) {
			if (step <= 7) {
				if (step == current_step) {
					++current_count;
				} else {
					CWRITE;
					current_step = step;
					current_count = 1;
				}
			} else {
				CWRITE;
				output[p++] = 0x80 | (step - 3);
				current_step = -100;
				current_count = 0;
			}
		} else {
			CWRITE;
			*((uint32_t*)&output[p]) = left_val;
			p += 4;
			output[p++] = -1;
			current_step = -100;
			current_count = 0;
		}

		left_val = right_val;
	}

	CWRITE;
	*((uint32_t*)&output[p]) = left_val;
	p += 4;

	output[p++] = -1;

	return p;
}

int CompressSeamLine(uint64_t* input, uint8_t* output, int width) {
	int current_step = -100;
	int current_count = 0;

	int64 step;
	int x = width;
	int p = 0;
	uint64_t left_val;
	uint64_t right_val;
	
	while (!((right_val = input[--x]) & 0xffffffff00000000) && x > 0);
	if (!(right_val & 0xffffffff00000000)) return p;

	while (x > 0) {
		while (!((left_val = input[--x]) & 0xffffffff00000000) && x > 0);
		if (!(left_val & 0xffffffff00000000)) break;

		if (!((right_val ^ left_val) & 0xffffffff) && (step = ((int64)(right_val - left_val) >> 32) + 3) < 67 && step >= 0) { // was <= 7
			if (step <= 7) {
				if (step == current_step) {
					++current_count;
				} else {
					CWRITE;
					current_step = (int)step;
					current_count = 1;
				}
			} else {
				CWRITE;
				output[p++] = 0x80 | ((int)step - 3);
				current_step = -100;
				current_count = 0;
			}
		} else {
			CWRITE;
			*((uint64_t*)&output[p]) = right_val;
			p += 8;
			output[p++] = -1;
			current_step = -100;
			current_count = 0;
		}

		right_val = left_val;
	}

	CWRITE;
	*((uint64_t*)&output[p]) = right_val;
	p += 8;

	output[p++] = -1;

	return p;
}

/***********************************************************************
* Wrap juggling
***********************************************************************/
void SwapUnswapH(Pyramid* py, bool unswap) {
	int width = py->GetWidth();
	int height = py->GetHeight();

	int minor_bytes = (width >> 1) << 2;
	int major_bytes = ((width + 1) >> 1) << 2;
	float* data = py->GetData();
	uint8_t* temp = (uint8_t*)malloc(major_bytes);

	for (int y = 0; y < height; ++y) {
		if (unswap) {
			memcpy(temp, &((uint8_t*)data)[minor_bytes], major_bytes);
			memcpy(&((uint8_t*)data)[major_bytes], data, minor_bytes);
			memcpy(data, temp, major_bytes);
		} else {
			memcpy(temp, data, major_bytes);
			memcpy(data, &((uint8_t*)data)[major_bytes], minor_bytes);
			memcpy(&((uint8_t*)data)[minor_bytes], temp, major_bytes);
		}

		data += py->GetPitch();
	}

	free(temp);
}

void SwapH(Pyramid* py) {
	SwapUnswapH(py, false);
}

void UnswapH(Pyramid* py) {
	SwapUnswapH(py, true);
}

void SwapUnswapV(Pyramid* py, bool unswap) {
	int height = py->GetHeight();
	int byte_pitch = py->GetPitch() << 2;
	uint8_t* temp = (uint8_t*)malloc(byte_pitch);
	int half_height = height >> 1;

	if (height & 1) {
		uint8_t* temp2 = (uint8_t*)malloc(byte_pitch);
		if (unswap) {
			uint8_t* upper = (uint8_t*)(py->GetData() + ((int64)height >> 1) * py->GetPitch());
			uint8_t* lower = (uint8_t*)(py->GetData() + ((int64)height - 1) * py->GetPitch());

			memcpy(temp, upper, byte_pitch);

			for (int y = 0; y < half_height; ++y) {
				memcpy(temp2, lower, byte_pitch);
				memcpy(lower, temp, byte_pitch);
				upper -= byte_pitch;
				memcpy(temp, upper, byte_pitch);
				memcpy(upper, temp2, byte_pitch);
				lower -= byte_pitch;
			}

			memcpy(lower, temp, byte_pitch);
		} else {
			uint8_t* upper = (uint8_t*)py->GetData();
			uint8_t* lower = (uint8_t*)(py->GetData() + ((int64)height >> 1) * py->GetPitch());

			memcpy(temp, lower, byte_pitch);

			for (int y = 0; y < half_height; ++y) {
				memcpy(temp2, upper, byte_pitch);
				memcpy(upper, temp, byte_pitch);
				lower += byte_pitch;
				memcpy(temp, lower, byte_pitch);
				memcpy(lower, temp2, byte_pitch);
				upper += byte_pitch;
			}

			memcpy(upper, temp, byte_pitch);
		}

		free(temp2);
	} else {
		uint8_t* upper = (uint8_t*)py->GetData();
		uint8_t* lower = (uint8_t*)(py->GetData() + (int64)half_height * py->GetPitch());
		for (int y = 0; y < half_height; ++y) {
			memcpy(temp, upper, byte_pitch);
			memcpy(upper, lower, byte_pitch);
			memcpy(lower, temp, byte_pitch);
			upper += byte_pitch;
			lower += byte_pitch;
		}
	}

	free(temp);
}

void SwapV(Pyramid* py) {
	SwapUnswapV(py, false);
}

void UnswapV(Pyramid* py) {
	SwapUnswapV(py, true);
}
