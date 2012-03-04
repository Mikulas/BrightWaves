/**
 * Bright Waves
 *
 * @Author Mikuláš Dítì
 * @Author Lubomír Grund
 * @License Original BSD, see LICENSE.txt
 */


#include <stdio.h>
#include <tchar.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <strstream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <boost/thread.hpp>



// System
#define MAX_THREADS 10
#define WIDTH 1024

// Config
#define DEPTH 1400
#define LAMBDA 25
#define AMPLIT 1
#define IND 1.33
#define STEP 1
#define SCALE 8
#define F 1
#define FI 0

#define PI 3.14159265
#define DER_STEP 1e7

// Static computations
#define PI2DLTS PI * 2 / (LAMBDA * SCALE) // 0.03141592

// Surfaces
#define S_PICASO 100
#define S_WOW 101
#define S_INTER 102
#define S_RAND2 103
#define S_RAND3 104
#define S_TRANS 105
#define S_STAND 106



// Globals
const int offset = WIDTH / 8;
const int boundary = WIDTH + offset;



// Forward declarations

/* Thread worker, computes all refractions at a time */
void worker(const double timeInMs);

/* Computes where surface at (x, y, t in PI2Ft) refracts into refractionsArray */
void refract(int* refractionsArray, const int x, const int y, const double PI2Ft);

/* Returns absolute water height at (x, y, t in PI2Ft) */
double surface(const double x, const double y, const double PI2Ft);

/* Saves refractions into a bitmap */
void writeBmp(int* refractionsArray, int brightnessFactor, const char* absoluteFilepath);

/* Saves configuration for future references */
void writeConfig(const char* absoluteFilepath);



int main(int argc, char *argv[])
{
	double time = 0;
	boost::thread workers[MAX_THREADS];
	int workerId = 0;
	
	std::cout << "Bright Waves\n############\n\n";
	std::cout << "computing...\n";

	writeConfig("run_1.txt");

	while (true) {
		workers[workerId] = boost::thread(worker, time);
		workerId++;
		time += 0.01;

		if (workerId == MAX_THREADS) {
			for (int i = 0; i < MAX_THREADS; ++i) {
				workers[i].join();
				std::cout << "\tworker[" << i << "]: done\n";
			}
			workerId = 0;
			std::cout << "started new batch of " << MAX_THREADS << " workers\n";
		}
	}

	return 0;
}

void worker(const double time)
{
	int *arref = (int*) malloc(WIDTH * WIDTH * sizeof(int));
	for (int i = 0; i < WIDTH * WIDTH; ++i) {
		arref[i] = 0;
	}

	const double PI2Ft = PI * 2 * F * time;
	for (int x = -offset; x < boundary; x += STEP) {
		for (int y = -offset; y < boundary; y += STEP) {
			refract(arref, x, y, PI2Ft);
		}
	}
	std::stringstream filename (std::stringstream::in | std::stringstream::out);
	filename << "run_1_" << std::setw(7) << std::setfill('0') << (time * 100) << ".bmp";
	writeBmp(arref, 8, filename.str().c_str());
	free(arref);
}

void writeBmp(int *arref, int brightness, const char* filename)
{
	FILE *f;
	unsigned char *img = NULL;
	int filesize = 54 + 3 * WIDTH * WIDTH;
	img = (unsigned char *) malloc(3 * WIDTH * WIDTH);
	memset(img, 0, sizeof(img));

	for (int x = 0; x < WIDTH; ++x) {
		for (int y = 0; y < WIDTH; ++y) {
			int v = arref[x * WIDTH + y] * brightness;
			if (v > 255)
				v = 255;

			img[(x + y * WIDTH) *3 + 2] = (unsigned char)(v);
			img[(x + y * WIDTH) *3 + 1] = (unsigned char)(v);
			img[(x + y * WIDTH) *3 + 0] = (unsigned char)(v);
		}
	}

	unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
	unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
	unsigned char bmppad[3] = {0,0,0};

	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(WIDTH);
	bmpinfoheader[5] = (unsigned char)(WIDTH >> 8);
	bmpinfoheader[6] = (unsigned char)(WIDTH >> 16);
	bmpinfoheader[7] = (unsigned char)(WIDTH >> 24);
	bmpinfoheader[8] = (unsigned char)(WIDTH);
	bmpinfoheader[9] = (unsigned char)(WIDTH >> 8);
	bmpinfoheader[10] = (unsigned char)(WIDTH >> 16);
	bmpinfoheader[11] = (unsigned char)(WIDTH >> 24);

	f = fopen(filename, "wb");
	fwrite(bmpfileheader,1,14,f);
	fwrite(bmpinfoheader,1,40,f);

	for (int i = 0; i < WIDTH; i++) {
		fwrite(img + (WIDTH * (WIDTH - i - 1) * 3), 3, WIDTH, f);
		fwrite(bmppad, 1, (4 - (WIDTH * 3) % 4) % 4, f);
	}
	free(img);
	fclose(f);
}

void refract(int *arref, const int x, const int y, const double PI2Ft)
{
	const double surf = surface(x, y, PI2Ft);
	const double dx = DER_STEP * (surface(x + (1/DER_STEP), y, PI2Ft) - surf);
    const double dy = DER_STEP * (surface(x, y + (1/DER_STEP), PI2Ft) - surf);

    const double gradh = sqrt(dx * dx + dy * dy);
	const double der_tmp = DEPTH * sqrt(1 / (IND * IND * (1 + gradh * gradh) - gradh * gradh));
    const int xxx = x + dx * der_tmp + 0.5; // intentional rounding
    const int yyy = y + dy * der_tmp + 0.5;
	
	if ((xxx > -1) && (xxx < WIDTH) && (yyy > -1) && (yyy < WIDTH)) {
		arref[xxx * WIDTH + yyy]++;
	}
}

double surface(const double x, const double y, const double PI2Ft)
{
	return AMPLIT * (
		sin(PI2DLTS * sqrt(x * x + y * y) + PI2Ft) +
		0.5 * sin(PI2DLTS * sqrt((200 - x) * (200 - x) + (200 - y) * (200 - y)) + PI2Ft + FI) + 
		0.5 * sin(PI2DLTS * x + PI2Ft + 2 * FI)
	) + DEPTH;
}

void writeConfig(const char* filename)
{
	std::ofstream out(filename); 
	out << "Width: \t" << WIDTH << "\n";
	out << "Depth: \t" << DEPTH << "\n";
	out << "Amplitude: \t" << AMPLIT << "\n";
	out << "Frequency: \t" << F << "\n";
	out << "FI: \t" << FI << "\n";
	out << "Lambda: \t" << LAMBDA << "\n";
	out << "IND: \t" << IND << "\n";
	out << "Step: \t" << STEP << "\n";
	out << "Derivate step: \t" << DER_STEP << "\n";
	out.close(); 
}

/*
// Implementation of other surface types

// Globals
double p = 0;
double q = 0;

double surface(int l, double x, double y)
{
	switch (l) {
		case 2:
			return (sin(PI2DLTS * x) + sin(PI2DLTS * y)) * AMPLIT + DEPTH;
		case 3:
			return sin(PI2DLTS * (x - y)) * AMPLIT + DEPTH;
		case 4:
			return (sin(PI2DLTS * sqrt(x * x + y * y)) + sin(PI2DLTS * sqrt(x * y))) * AMPLIT + DEPTH;
		case 5:
			return 3/2 * AMPLIT * (
				sin(PI2DLTS * sqrt((x - 102) * (x - 102) + (y - 78) * (y - 78))) + 
				sin(PI2DLTS * sqrt((x - 33) * (x - 33) + (y - 147) * (y - 147))) + 
				sin(PI2DLTS * sqrt((x - 168) * (x - 168) + (y - 27) * (y - 27)))
			) + DEPTH;
		case 6:
			return AMPLIT * (
				sin(PI2DLTS * sqrt((x - 102) * (x - 102) + (y - 78) * (y - 78)) + (2 * PI * F * t)) +
				sin(PI2DLTS * sqrt((x - 33) * (x - 33) + (y - 147) * (y - 147)) + (2 * PI * F * t))
			) + DEPTH;
		case 7:
			return AMPLIT * (
				1.5 * sin(PI2DLTS * sqrt((x - 102) * (x - 102) + (y - 78) * (y - 78)) + (2 * PI * F * t)) +
				0.5 * sin(PI2DLTS * sqrt((x - 33) * (x - 33) + (y - 147) * (y - 147)) + (2 * PI * F * t))
			) + DEPTH;
		case 8:
			return AMPLIT * (1.5 * sin(PI2DLTS * sqrt((x - 102) * (x - 102) + (y - 78) * (y - 78)))) + DEPTH;
		case 9:
			return 1/2 * AMPLIT * (
				sin(PI2DLTS * sqrt((x - 50) * (x - 50) + (y - 50) * (y - 50))) +
				sin(PI2DLTS * sqrt((x - 50) * (x - 50) + (y - 150) * (y - 150))) + 
				sin(PI2DLTS * sqrt((x - 150) * (x - 150) + (y - 50) * (y - 50))) + 
				sin(PI2DLTS * sqrt((x - 150) * (x - 150) + (y - 150) * (y - 150)))
			) + DEPTH;
		case 10:
			return AMPLIT * (sin(PI2DLTS * y / 2) + sin(PI2DLTS * (y - x * sqrt(3.0)) / (sqrt(2.0) * 2))) + DEPTH;
		case 11:
			return AMPLIT * (sin(PI2DLTS * y / 4) + sin(PI2DLTS * (y - x * sqrt(3.0)) / (4))) + DEPTH;
		case 12:
			return AMPLIT * (sin(PI2DLTS * y / 4) + sin(PI2DLTS * (y - x * sqrt(3.0)) / (sqrt(2.0) * 4))) + DEPTH;
		case 13:
			return AMPLIT * (sin(PI2DLTS * y / 2) + sin(PI2DLTS * (y - x * sqrt(3.0)) / (sqrt(2.0) * 2))) + DEPTH;
		case 14:
			return AMPLIT * (
				sin(PI2DLTS * y / 2) + 
				sin(PI2DLTS * (y - x * sqrt(3.0)) / (sqrt(2.0) * 2)) +
				sin(PI2DLTS * sqrt((x - 102) * (x - 102) + (y - 78) * (y - 78)))
			) + DEPTH;
		case 15:
			return AMPLIT * (
				sin(PI2DLTS * (sqrt((x - 150) * (x - 150) + (y - 100) * (y - 100))) / 0.9) + 
				sin(PI2DLTS * (sqrt((x - 50) * (x - 50) + (y - 100) * (y - 100))) / 1.5)
			) + DEPTH;
		case 16:
			if ((2 * (x - 300) * (x - 300) - 3 * (x - 300) + (y - 300) - 2) > 0)
				return DEPTH + AMPLIT;
			else
				return DEPTH - AMPLIT;
		case 17:
			return AMPLIT * (
				sin(PI2DLTS * (sqrt((x - 150) * (x - 150) + (y - 100) * (y - 100))) / 0.5) + 
				sin(PI2DLTS * (sqrt((x - 50) * (x - 50) + (y - 100) * (y - 100))) / 0.7)
			) + DEPTH;
		case 18:
			return AMPLIT * (
				sin(PI2DLTS * y / 2 + (2 * PI * F * t)) + 
				sin(PI2DLTS * (y - x * sqrt(3.0)) / (sqrt(2.0) * 2) + (2 * PI * F * t) + FI)
			) + DEPTH;

		case S_WOW:
			return AMPLIT * (
				sin(PI2DLTS * sqrt(x * x + y * y)) + 
				sin(PI2DLTS * x)
			) + DEPTH;
		case S_INTER:
			return AMPLIT * (
				sin(PI2DLTS * sqrt(x * x + y * y)) +
				sin(PI2DLTS * sqrt((200 - x) * (200 - x) + (200 - y) * (200 - y)))
			) + DEPTH;
		case S_PICASO:
			return AMPLIT * (
				sin(PI2DLTS * sqrt(x * x + y * y) + (2 * PI * F * t)) +
				0.5 * sin(PI2DLTS * sqrt((200 - x) * (200 - x) + (200 - y) * (200 - y)) + (2 * PI * F * t) + FI) + 
				0.5 * sin(PI2DLTS * x + (2 * PI * F * t) + 2 * FI)
			) + DEPTH;
		case S_RAND2:
			return AMPLIT * sin(PI2DLTS * ((x - p) * (x - p) + (y - q) * (y - q))) + DEPTH;
		case S_RAND3:
			return AMPLIT * sin(PI2DLTS * (sqrt((x - 50) * (x - 50) + (y - 122) * (y - 122)))) + DEPTH;
		case S_TRANS:
			return AMPLIT * (
				sin(2 * PI2DLTS * y + (2 * PI * F * t)) +
				sin(2 * PI2DLTS * (y - x * sqrt(3.0)) / sqrt(2.0) + (2 * PI * F * t) + FI)
			) / 4 + DEPTH;
		case S_STAND:
			return AMPLIT * (
				sin(2 * PI2DLTS * y) * sin((2 * PI * F * t)) + 
				sin(2 * PI2DLTS * (y - x * sqrt(3.0)) / sqrt(2.0)) * sin((2 * PI * F * t) + FI)
			) + DEPTH;
		default:
			return surface(S_PICASO, x, y);
	}
}*/
