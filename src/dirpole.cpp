#include <math.h>   // dirpole, directional dipole by T. Hachisuka and J. R. Frisvad
#include <stdlib.h> // originally smallpt, a path tracer by Kevin Beason, 2008
#include <stdio.h>  // Usage: ./dirpole 100000 && xv image.ppm
#define PI 3.14159265358979 //          ^^^^^^:number of samples per pixel

// linear congruential PRNG
inline unsigned int next_rand(unsigned int current) {return current * 3125u + 49u;}
const double rand_range = 4294967296.0;

// Halton sequence with reverse permutation
const int primes[20]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71};
inline int rev(const int i, const int p) {return (i == 0) ?  i : p - i;}
inline double hal(const int b, int j) {
	const int p = primes[b]; double h = 0.0, f = 1.0 / (double)p, fct = f;
	while (j > 0) {h += rev(j % p, p) * fct; j /= p; fct *= f;} return h;
}

struct Vec {union {struct{double x, y, z;}; struct{double r, g, b;};}; // vector: position, also color (r,g,b)
	Vec(double x_ = 0, double y_ = 0, double z_ = 0) {x = x_; y = y_; z = z_;}
	inline double& operator[](const int i) {return *(&x + i);}
	inline const double& operator[](const int i) const {return *(&x + i);}
	inline Vec operator-() const {return Vec(-x, -y, -z);}
	inline Vec operator+(const Vec& b) const {return Vec(x + b.x, y + b.y, z + b.z);}
	inline Vec operator-(const Vec& b) const {return Vec(x - b.x, y - b.y, z - b.z);}
	inline Vec operator*(const Vec& b) const {return Vec(x * b.x, y * b.y, z * b.z);}
	inline Vec operator/(const Vec& b) const {return Vec(x / b.x, y / b.y, z / b.z);}
	inline Vec operator+(double b) const {return Vec(x + b, y + b, z + b);}
	inline Vec operator-(double b) const {return Vec(x - b, y - b, z - b);}
	inline Vec operator*(double b) const {return Vec(x * b, y * b, z * b);}
	inline Vec operator/(double b) const {return Vec(x / b, y / b, z / b);}
	inline friend Vec operator+(double b, const Vec& v) {return Vec(b + v.x, b + v.y, b + v.z);}
	inline friend Vec operator-(double b, const Vec& v) {return Vec(b - v.x, b - v.y, b - v.z);}
	inline friend Vec operator*(double b, const Vec& v) {return Vec(b * v.x, b * v.y, b * v.z);}
	inline friend Vec operator/(double b, const Vec& v) {return Vec(b / v.x, b / v.y, b / v.z);}
	inline double len() const {return sqrt(x * x + y * y + z * z);}
	inline Vec normalized() const {return (*this) / this->len();}
	inline friend Vec sqrt(const Vec& b) {return Vec(sqrt(b.x), sqrt(b.y), sqrt(b.z));}
	inline friend Vec exp(const Vec& b) {return Vec(exp(b.x), exp(b.y), exp(b.z));}
	inline double dot(const Vec& b) const {return x * b.x + y * b.y + z * b.z;}
	inline Vec operator%(const Vec& b) const {return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);}
};

struct Ray {Vec o, d; Ray(){}; Ray(Vec o_, Vec d_) : o(o_), d(d_) {}};

const struct Sphere {double rad; Vec p;
	Sphere(double r_,Vec p_) : rad(r_),p(p_){}
	inline double intersect(const Ray& r) const {
		// ray-sphere intersection returns the distance
		const Vec op = p - r.o;
		double t, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
		if (det < 0) {
			return 1e20;
		}
		else {
			det = sqrt(det);
		}
		return (t = b - det) > 1e-4 ? t : ((t = b + det) > 1e-4 ? t : 1e20);
	}
} sph = Sphere(25, Vec(50, 35, 60));

// find the closest intersection
inline bool intersect(const Ray& r,double& t) {
	double d, inf = 1e20;
	t = inf;
	d = sph.intersect(r);
	if (d < t) t = d; 
	return t < inf;
}

// uniformly sample a point on the sphere
const double samplePDF = 1.0 / (4.0 * PI * sph.rad * sph.rad);
inline void sample(const int i, const Vec& offset, Vec& sp, Vec& sn) {
	double r0 = hal(0, i) + offset.x;
	double r1 = hal(1, i) + offset.y;
	if (r0 > 1.0) r0 -= 1.0;
	if (r1 > 1.0) r1 -= 1.0;
	const double p = 2.0 * PI * r0, t = 2.0 * acos(sqrt(1.0 - r1));
	const Vec d = Vec(cos(p) * sin(t), cos(t), sin(p) * sin(t));
	sp = sph.p + d * sph.rad;
	sn = d;
}

// index of refraction 
const double eta = 1.3;

// potato
const Vec sigma_s = Vec(0.68, 0.70, 0.55);
const Vec sigma_a = Vec(0.0024, 0.0090, 0.12);
const Vec g = Vec(0.0, 0.0, 0.0);

// white grapefruit
//const Vec sigma_s = Vec(0.3513, 0.3669, 0.5237);
//const Vec sigma_a = Vec(0.3609, 0.3800, 0.5632) - sigma_s;
//const Vec g = Vec(0.548, 0.545, 0.565);

inline double min3(const double x, const double y, const double z) {
	const double r = x < y ? x : y; return r < z ? r : z;
}

inline double C1(const double n) {
	double r;
	if (n > 1.0) {
		r = -9.23372 + n * (22.2272 + n * (-20.9292 + n * (10.2291 + n * (-2.54396 + 0.254913 * n))));
	} else {
		r = 0.919317 + n * (-3.4793 + n * (6.75335 + n *  (-7.80989 + n *(4.98554 - 1.36881 * n))));
	}
	return r / 2.0;
}
inline double C2(const double n) {
	double r = -1641.1 + n * (1213.67 + n * (-568.556 + n * (164.798 + n * (-27.0181 + 1.91826 * n))));
	r += (((135.926 / n) - 656.175) / n + 1376.53) / n;
	return r / 3.0;
}

// setup some constants
const Vec sigma_t = sigma_s + sigma_a;
const Vec sigma_sp = sigma_s * (1.0 - g);
const Vec sigma_tp = sigma_sp + sigma_a;
const Vec albedo_p = sigma_sp / sigma_tp;
const Vec D = 1.0 / (3.0 * sigma_tp);
const Vec sigma_tr = sqrt(sigma_a / D);
const Vec de = 2.131 * D / sqrt(albedo_p);
const double Cp_norm = 1.0 / (1.0 - 2.0 * C1(1.0 / eta));
const double Cp = (1.0 - 2.0 * C1(eta)) / 4.0;
const double Ce = (1.0 - 3.0 * C2(eta)) / 2.0;
const double A = (1.0 - Ce) / (2.0 * Cp);
const double min_sigma_tr = min3(sigma_tr.r, sigma_tr.g, sigma_tr.b);

// directional dipole
// --------------------------------
inline double Sp_d(const Vec& x, const Vec& w, const double& r, const Vec& n, const int j) {
	// evaluate the profile
	const double s_tr_r = sigma_tr[j] * r;
	const double s_tr_r_one = 1.0 + s_tr_r;
	const double x_dot_w = x.dot(w);
	const double r_sqr = r * r;

	const double t0 = Cp_norm * (1.0 / (4.0 * PI * PI)) * exp(-s_tr_r) / (r * r_sqr);
	const double t1 = r_sqr / D[j] + 3.0 * s_tr_r_one * x_dot_w;
	const double t2 = 3.0 * D[j] * s_tr_r_one * w.dot(n);
	const double t3 = (s_tr_r_one + 3.0 * D[j] * (3.0 * s_tr_r_one + s_tr_r * s_tr_r) / r_sqr * x_dot_w) * x.dot(n);

	return t0 * (Cp * t1 - Ce * (t2 - t3));
}
inline double bssrdf(const Vec& xi, const Vec& ni, const Vec& wi, const Vec& xo, const Vec& no, const Vec& wo, const int j) {
	// distance
	const Vec xoxi = xo - xi;
	const double r = xoxi.len();

	// modified normal
	const Vec ni_s = (xoxi.normalized()) % ((ni % xoxi).normalized());

	// directions of ray sources
	const double nnt = 1.0 / eta, ddn = -wi.dot(ni);
	const Vec wr = (wi * -nnt - ni * (ddn * nnt + sqrt(1.0 - nnt * nnt * (1.0 - ddn * ddn)))).normalized();
	const Vec wv = wr - ni_s * (2.0 * wr.dot(ni_s));

	// distance to real sources
	const double cos_beta = -sqrt((r * r - xoxi.dot(wr) * xoxi.dot(wr)) / (r * r + de[j] * de[j]));
	double dr;
	const double mu0 = -no.dot(wr);
	if (mu0 > 0.0) {
		dr = sqrt((D[j] * mu0) * ((D[j] * mu0) - de[j] * cos_beta * 2.0) + r * r);
	} else {
		dr = sqrt(1.0 / (3.0 * sigma_t[j] * 3.0 * sigma_t[j]) + r * r);
	}

	// distance to virtual source
	const Vec xoxv = xo - (xi + ni_s * (2.0 * A * de[j]));
	const double dv = xoxv.len();

	// BSSRDF
	const double result = Sp_d(xoxi, wr, dr, no, j) - Sp_d(xoxv, wv, dv, no, j);

	// clamping to zero
	return (result < 0.0) ? 0.0 : result;
}
// --------------------------------

inline double fresnel(const double cos_theta, const double eta) {
	const double sin_theta_t_sqr = 1.0 / (eta * eta) * (1.0 - cos_theta * cos_theta);
	if (sin_theta_t_sqr >= 1.0) return 1.0;
	const double cos_theta_t = sqrt(1.0 - sin_theta_t_sqr);
	const double r_s = (cos_theta - eta * cos_theta_t) / (cos_theta + eta * cos_theta_t);
	const double r_p = (eta * cos_theta - cos_theta_t) / (eta * cos_theta + cos_theta_t);
	return (r_s * r_s + r_p * r_p) * 0.5;
}

Vec trace(const Ray& r, const int index, const int num_samps) {
	// ray-sphere intersection
	double t;
	if (!intersect(r, t)) return Vec();

	// compute the intersection data
	const Vec x = r.o + r.d * t, n = (x - sph.p).normalized(), w = -r.d;

	// compute Fresnel transmittance at point of emergence
	const double T21 = 1.0 - fresnel(w.dot(n), eta);

	// integration of the BSSRDF over the surface
	unsigned int xi = next_rand(index);
	Vec result;
	for (int i = 0; i < num_samps; i++) {
		// generate a sample
		Vec sp, sn, sw;
		sample(i, Vec(hal(2, index), hal(3, index), 0.0), sp, sn);
		sw = Vec(1, 1, -0.5).normalized();

		// directional light source
		const double radiance = 8.0*PI;
		const double cos_wi_n = sn.dot(sw);
		if (cos_wi_n > 0.0) {
			// Russian roulette (note that accept_prob can be any value in (0, 1))
			const double accept_prob = exp(-(sp - x).len() * min_sigma_tr);
			if ((xi / rand_range) < accept_prob) {
				const double T12 = 1.0 - fresnel(cos_wi_n, eta);
				const double Li_cos = T12 * radiance * cos_wi_n / (samplePDF * accept_prob);

				for (int j = 0; j < 3; j++) result[j] += bssrdf(sp, sn, sw, x, n, w, j) * Li_cos;

				// reciprocal evaulation with the reciprocity hack
				//for (int j = 0; j < 3; j++) result[j] += 0.5 * (bssrdf(sp, sn, sw, x, n, w, j) + bssrdf(x, n, w, sp, sn, sw, j)) * Li_cos;
			}
			xi = next_rand(xi);
		}
	}
	return T21 * result / (double)num_samps;
}

int main(int argc, char* argv[]) {
	const int w = 512, h = 512, samps = (argc == 2) ? atoi(argv[1]) : 1 << 13;

	// main rendering process
	const Ray cam(Vec(50, 45, 290), Vec(0, -0.042612, -1).normalized());
	const Vec cx = Vec(w * 0.25 / h), cy = (cx % cam.d).normalized() * 0.25;
	Vec *c = new Vec[w * h];
	#pragma omp parallel for schedule(dynamic, 1)
	for (int y = 0; y < h; y++) {
		fprintf(stderr, "\rRendering %5.2f%%", 100.0 * y / (h - 1));
		for (int x = 0; x < w; x++) {
			const Vec d = cx * ((x + 0.5) / w - 0.5) + cy * (-(y + 0.5) / h + 0.5) + cam.d;
			const int i = x + y * w;
			c[i] = trace(Ray(cam.o + d * 140, d.normalized()), i, samps);
		}
	}

	// save the HDR image
	float *cf = new float[w * h * 3];
	for (int i = 0; i < w * h; i++) {
		for (int j = 0; j < 3; j++) cf[i * 3 + j] = c[i][j];
	}
	FILE *f = fopen("image.pfm", "wb"); 
	fprintf(f, "PF\n%d %d\n%6.6f\n", w, h, -1.0);
	fwrite(cf, w * h * 3, sizeof(float), f);
}

int render(){
    return 0;
}
