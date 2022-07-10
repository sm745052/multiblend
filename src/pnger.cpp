#include "png.h"
#include <math.h>
#define PNGER

void Output(int level, const char* fmt, ...);

class Pnger {
public:
	Pnger(const char* filename, const char* name, int w, int _h, int type, int bpp = 8, FILE* _f = NULL, int compression = -1);
	~Pnger();
	bool Ready() { return !!f; };
	void WriteRows(uint8_t** rows, int num_rows);
	void Write();
	static void Quick(char* filename, uint8_t* data, int width, int height, int pitch, int type);
	uint8_t* line;

private:
	png_structp png_ptr;
	png_infop info_ptr;
	static png_color* palette;
	FILE* f;
	int y;
	int h;
};

png_color* Pnger::palette = NULL;

Pnger::Pnger(const char* filename, const char* name, int w, int _h, int type, int bpp, FILE* _f, int compression) {
	y = 0;
	h = _h;

	if (type == PNG_COLOR_TYPE_PALETTE && !palette) {
		palette = (png_color*)malloc(256 * sizeof(png_color));

		double base = 2;
		double rad;
		double r, g, b;

		for (int i = 0; i < 255; ++i) {
			rad = base;
			r = std::max(0.0, std::min(1.0, std::min(rad, 4 - rad)));
			rad += 2; if (rad >= 6) rad -= 6;
			g = std::max(0.0, std::min(1.0, std::min(rad, 4 - rad)));
			rad += 2; if (rad >= 6) rad -= 6;
			b = std::max(0.0, std::min(1.0, std::min(rad, 4 - rad)));
			base += 6 * 0.618033988749895;
			if (base >= 6) base -= 6;
			palette[i].red = (png_byte)(sqrt(r) * 255 + 0.5);
			palette[i].green = (png_byte)(sqrt(g) * 255 + 0.5);
			palette[i].blue = (png_byte)(sqrt(b) * 255 + 0.5);
		}

		palette[255].red = 0;
		palette[255].green = 0;
		palette[255].blue = 0;
	}

	if (!_f) {
		fopen_s(&f, filename, "wb");
		if (!f) {
			if (name) Output(0, "WARNING: Could not save %s\n", name);
			return;
		}
	} else {
		f = _f;
	}

	png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	if (!png_ptr) {
		if (name) Output(0, "WARNING: PNG create failed\n");
		fclose(f);
		remove(filename);
		f = NULL;
		return;
	}

	info_ptr = png_create_info_struct(png_ptr);
	if (!info_ptr) {
		if (name) Output(0, "WARNING: PNG create failed\n");
		png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
		fclose(f);
		remove(filename);
		f = NULL;
		return;
	}

	png_init_io(png_ptr, f);

	png_set_IHDR(png_ptr, info_ptr, w, h, bpp, type, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
	if (type == PNG_COLOR_TYPE_PALETTE) png_set_PLTE(png_ptr, info_ptr, palette, 256);

	png_write_info(png_ptr, info_ptr);
	png_set_compression_level(png_ptr, compression < 0 ? 3 : compression);
	if (bpp == 16) png_set_swap(png_ptr);

	line = name ? new uint8_t[(type == PNG_COLOR_TYPE_RGB_ALPHA ? (w << 2) : type == PNG_COLOR_TYPE_RGB ? w * 3 : w) << (bpp >> 4)] : NULL;
};

void Pnger::Write() {
	if (!f) return;

	png_write_row(png_ptr, (uint8_t*)line);

	if (++y == h) {
//		printf("png close\n");
		png_write_end(png_ptr, NULL);
		png_destroy_write_struct(&png_ptr, &info_ptr);
		fclose(f);
		f = NULL;
	}
}

void Pnger::WriteRows(uint8_t** rows, int num_rows) {
	if (!f) return;

	png_write_rows(png_ptr, rows, num_rows);

	if ((y += num_rows) == h) {
		png_write_end(png_ptr, NULL);
		png_destroy_write_struct(&png_ptr, &info_ptr);
		fclose(f);
		f = NULL;
	}
}

void Pnger::Quick(char* filename, uint8_t* data, int width, int height, int pitch, int type) {
	Pnger temp(filename, NULL, width, height, type);

	for (int y = 0; y < height; ++y) {
		temp.WriteRows(&data, 1);
		data += (type == PNG_COLOR_TYPE_RGB_ALPHA) ? (pitch << 2) : pitch;
	}
}

Pnger::~Pnger() {
	delete line;
	if (f) fclose(f);
}
