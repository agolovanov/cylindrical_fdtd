#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <functional>
#include <H5Cpp.h>
using namespace std;

template <typename T>
class array3d {
private:
    int n1, n2, n3;
    std::unique_ptr<T[]> data;

public:
    array3d() : n1(0), n2(0), n3(0) {}

    array3d(int n1, int n2, int n3) : n1(n1), n2(n2), n3(n3) {
        data = std::unique_ptr<T[]>(new T[n1 * n2 * n3]);
    }

    inline T & operator()(int i, int j, int k) const {
        return data[n2 * n3 * i + n3 * j + k];
    }

    int get_n1() const { return n1; }
    int get_n2() const { return n2; }
    int get_n3() const { return n3; }

    array3d<T> & operator=(const array3d<T> & other) {
        n1 = other.n1;
        n2 = other.n2;
        n3 = other.n3;

        const int size = n1 * n2 * n3;
        data = std::unique_ptr<T[]>(new T[size]);

        for (int i = 0; i < size; i++) {
            data[i] = other.data[i];
        }

        return *this;
    }
};

struct vector3d {
    double r, phi, x;
};

struct cvector3d {
    double x, y, z;
};

struct complex {
    double re = 0, im = 0;

    complex() = default;
    complex(double re, double im) : re(re), im(im) {};
};

template <typename T>
class array2d {
private:
    int n1, n2;
    std::unique_ptr<T[]> data;

public:
    array2d() : n1(0), n2(0) {}

    array2d(int n1, int n2) : n1(n1), n2(n2) {
        data = std::unique_ptr<T[]>(new T[n1 * n2]);
    }

    inline T & operator()(int i, int j) const {
        return data[n2 * i + j];
    }

    int get_n1() const { return n1; }
    int get_n2() const { return n2; }

    array3d<T> & operator=(const array3d<T> & other) {
        n1 = other.n1;
        n2 = other.n2;

        const int size = n1 * n2;
        data = std::unique_ptr<T[]>(new T[size]);

        for (int i = 0; i < size; i++) {
            data[i] = other.data[i];
        }

        return *this;
    }
};

array3d<vector3d> e;
array3d<vector3d> b;
vector<double> r, phi;
vector<double> ex0;
vector<complex> ephi0, br0;

constexpr double pi = 3.14159265358979323846;
int nx, nr, nphi;
double dx, dr, dphi, dt;


void advance_e(double dt) {
    // advancing axis values
    for (int k = 0; k < nx; k++) {
        int km = (k > 0 ? k-1 : nx-1);
        //ex0[k] = 0;
        //ephi0[k].re = 0;
        //ephi0[k].im = 0;
        double tmp1 = 0;
        complex tmp2;
        tmp2.re = 0;
        tmp2.im = 0;
        for (int j = 0; j < nphi; j++) {
            tmp1 += b(0, j, k).phi;
            tmp2.re += b(0, j, k).x * cos(phi[j] + 0.5 * dphi);
            tmp2.im -= b(0, j, k).x * sin(phi[j] + 0.5 * dphi);
        }
        ex0[k] += dt * 2 * dphi / pi / dr * tmp1;
        ephi0[k].re += - dt * 2 * dphi / pi / dr * tmp2.re;
        ephi0[k].im += - dt * 2 * dphi / pi / dr * tmp2.im;
        ephi0[k].re += dt * (br0[k].re - br0[km].re) / dx;
        ephi0[k].im += dt * (br0[k].im - br0[km].im) / dx;
    }

    int i = 0;
    for (int j = 0; j < nphi; j++) {
        for (int k = 0; k < nx; k++) {
            int jm = (j > 0 ? j-1 : nphi-1);
            int km = (k > 0 ? k-1 : nx-1);
            e(i,j,k).r += dt * ((b(i,j,k).x - b(i,jm,k).x) / (r[i] + 0.5 * dr) / dphi - (b(i,j,k).phi - b(i,j,km).phi) / dx);
            e(i,j,k).phi = ephi0[k].re * cos(phi[j] + 0.5 * dphi) - ephi0[k].im * sin(phi[j] + 0.5 * dphi);
            e(i,j,k).x = ex0[k];
        }
    }
    for (i = 1; i < nr-1; i++) {
        for (int j = 0; j < nphi; j++) {
            for (int k = 0; k < nx; k++) {
                int jm = (j > 0 ? j-1 : nphi-1);
                int km = (k > 0 ? k-1 : nx-1);
                e(i,j,k).r += dt * ((b(i,j,k).x - b(i,jm,k).x) / (r[i] + 0.5 * dr) / dphi - (b(i,j,k).phi - b(i,j,km).phi) / dx);
                e(i,j,k).phi += dt * ((b(i,j,k).r - b(i,j,km).r) / dx - (b(i,j,k).x - b(i-1,j,k).x) / dr);
                e(i,j,k).x += dt * ((b(i,j,k).phi - b(i-1,j,k).phi) / dr + 0.5 * (b(i,j,k).phi + b(i-1,j,k).phi) / r[i] - (b(i,j,k).r - b(i,jm,k).r) / (r[i] * dphi));
            }
        }
    }
}

void advance_b(double dt) {
    for (int k = 0; k < nx; k++) {
        int kp = (k < nx-1 ? k+1 : 0);
        //br0[k].re = 0;
        //br0[k].im = 0;
        complex tmp;
        tmp.re = 0;
        tmp.im = 0;
        for (int j = 0; j < nphi; j++) {
            tmp.re += e(1, j, k).x * sin(phi[j]);
            tmp.im += e(1, j, k).x * cos(phi[j]);
        }
        br0[k].re -= dt * dphi / pi / dr * tmp.re;
        br0[k].im -= dt * dphi / pi / dr * tmp.im;
        br0[k].re += dt * (ephi0[kp].re - ephi0[k].re) / dx;
        br0[k].im += dt * (ephi0[kp].im - ephi0[k].im) / dx;
    }

    int i = 0;
    for (int j = 0; j < nphi; j++) {
        for (int k = 0; k < nx; k++) {
            int jp = (j < nphi-1 ? j+1 : 0);
            int kp = (k < nx-1 ? k+1 : 0);
            b(i,j,k).r = br0[k].re * cos(phi[j] + 0.5 * dphi) - br0[k].im * sin(phi[j] + 0.5 * dphi);
            b(i,j,k).phi -= dt * ((e(i,j,kp).r - e(i,j,k).r) / dx - (e(i+1,j,k).x - e(i,j,k).x) / dr);
            b(i,j,k).x -= dt * ((e(i+1,j,k).phi - e(i,j,k).phi) / dr + 0.5 * (e(i+1,j,k).phi + e(i,j,k).phi) / (r[i] + 0.5 * dr) - (e(i,jp,k).r - e(i,j,k).r) / (r[i] + 0.5 * dr) / dphi);
        }
    }
    for (i = 1; i < nr-1; i++) {
        for (int j = 0; j < nphi; j++) {
            for (int k = 0; k < nx; k++) {
                int jp = (j < nphi-1 ? j+1 : 0);
                int kp = (k < nx-1 ? k+1 : 0);
                b(i,j,k).r -= dt * ((e(i,jp,k).x - e(i,j,k).x) / r[i] / dphi - (e(i,j,kp).phi - e(i,j,k).phi) / dx);
                b(i,j,k).phi -= dt * ((e(i,j,kp).r - e(i,j,k).r) / dx - (e(i+1,j,k).x - e(i,j,k).x) / dr);
                b(i,j,k).x -= dt * ((e(i+1,j,k).phi - e(i,j,k).phi) / dr + 0.5 * (e(i+1,j,k).phi + e(i,j,k).phi) / (r[i] + 0.5 * dr) - (e(i,jp,k).r - e(i,j,k).r) / (r[i] + 0.5 * dr) / dphi);
            }
        }
    }
}

cvector3d cylindrical_to_cartesian(const vector3d v) {
    cvector3d res;
    res.x = v.x;
    res.y = v.r * cos(v.phi);
    res.z = v.r * sin(v.phi);
    return res;
}

void initialize_fields(function<double(cvector3d)> ex, function<double(cvector3d)> ey, function<double(cvector3d)> ez,
                       function<double(cvector3d)> bx, function<double(cvector3d)> by, function<double(cvector3d)> bz) {

    cout << "Initialize 1D fields" << endl;
    int i = 0;
    for (int k = 0; k < nx; k++) {
        ex0[k] = ex({0, 0, (k + 0.5) * dx});

        ephi0[k].re = ez({0, 0, k * dx});
        ephi0[k].im = ey({0, 0, k * dx});

        br0[k].re = by({0, 0, (k + 0.5) * dx});
        br0[k].im = - bz({0, 0, (k + 0.5) * dx});
    }

    cout << "Initialize 3D fields" << endl;
    for (i = 0; i < nr-1; i++) {
        for (int j = 0; j < nphi; j++) {
            for (int k = 0; k < nx; k++) {
                //cout << "Init (" << i << ", " << j << ", " << k << ")" << endl;
                cvector3d point = cylindrical_to_cartesian({r[i], phi[j], (k + 0.5) * dx});
                e(i,j,k).x = ex(point);

                double point_phi = phi[j];
                point = cylindrical_to_cartesian({r[i] + 0.5 * dr, point_phi, k * dx});
                e(i,j,k).r = ey(point) * cos(point_phi) + ez(point) * sin(point_phi);
                //cout << "Er (" << r[i] + 0.5 * dr << ", " << point_phi << ", " << k * dx << ")->(" 
                //       << point.x << ", " << point.y << ", " << point.z << ") = " << e(i,j,k).r << endl;
                

                point_phi = phi[j] + 0.5 * dphi;
                point = cylindrical_to_cartesian({r[i], point_phi, k * dx});
                e(i,j,k).phi = - ey(point) * sin(point_phi) + ez(point) * cos(point_phi);
                //cout << "Ephi (" << r[i] << ", " << point_phi << ", " << k * dx << ")->(" 
                //        << point.x << ", " << point.y << ", " << point.z << ") = " << e(i,j,k).phi << endl;

                point = cylindrical_to_cartesian({r[i] + 0.5 * dr, point_phi, k * dx});
                b(i,j,k).x = bx(point);

                point = cylindrical_to_cartesian({r[i], point_phi, (k + 0.5) * dx});
                b(i,j,k).r = by(point) * cos(point_phi) + bz(point) * sin(point_phi);

                point_phi = phi[j];
                point = cylindrical_to_cartesian({r[i] + 0.5 * dr, point_phi, (k + 0.5) * dx});
                b(i,j,k).phi = - by(point) * sin(point_phi) + bz(point) * cos(point_phi);
            }
        }
    }
    cout << "Finish initialize fields" << endl;
}

void write_array(array3d<double> & field, std::string name, H5::H5File file) {
    const hsize_t dims[3] {field.get_n1(), field.get_n2(), field.get_n3()};
    H5::DataSpace dataspace(3, dims);

    H5::DataSet dataset = file.createDataSet(name, H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&(field(0,0,0)), H5::PredType::NATIVE_DOUBLE);
}

void output_fields(int number) {
    string filename = "Fields" + to_string(number) + ".h5";
    H5::H5File fields_file(filename, H5F_ACC_TRUNC);

    cout << "Writing at " << to_string(number) << endl;

    array3d<double> er(nr, nphi, nx);
    array3d<double> ephi(nr, nphi, nx);
    array3d<double> ex(nr, nphi, nx);
    array3d<double> br(nr, nphi, nx);
    array3d<double> bphi(nr, nphi, nx);
    array3d<double> bx(nr, nphi, nx);

    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nphi; j++) {
            for (int k = 0; k < nx; k++) {
                er(i,j,k) = e(i,j,k).r;
                ephi(i,j,k) = e(i,j,k).phi;
                ex(i,j,k) = e(i,j,k).x;
                br(i,j,k) = b(i,j,k).r;
                bphi(i,j,k) = b(i,j,k).phi;
                bx(i,j,k) = b(i,j,k).x;
            }
        }
    }

    write_array(er, "Er", fields_file);
    write_array(ephi, "Ephi", fields_file);
    write_array(ex, "Ex", fields_file);
    write_array(bx, "Bx", fields_file);
    write_array(bphi, "Bphi", fields_file);
    write_array(br, "Br", fields_file);
}

int main() {
    double lx = 10;
    double lr = 10;

    nx = 200;
    // nx = 1;
    nr = 50;
    // nr = 10;
    nphi = 30;
    // nphi = 10;

    dx = lx / nx;
    dr = lr / nr;
    dphi = 2 * pi / nphi;

    double tmax = 10;
    double dt = 0.01;
    double output_dt = 0.5;

    r.resize(nr);
    for (int i = 0; i < nr; i++) {
        r[i] = i * dr;
    }
    phi.resize(nphi);
    for (int j = 0; j < nphi; j++) {
        phi[j] = j * dphi;
    }

    e = array3d<vector3d>{nr, nphi, nx};
    b = array3d<vector3d>{nr, nphi, nx};
    ex0.resize(nx);
    ephi0.resize(nx);
    br0.resize(nx);

    int current_iteration = 0;
    int output_iteration = 0;

    double x0 = lx / 2 - 1;
    double y0 = -2;
    double angle = 30 * pi / 180;
    double sigma = 3;

    auto empty_func = [] (cvector3d r) -> double {return 0.0;};
    auto ex_func = [sigma, x0, y0, angle] (cvector3d r) -> double {return -sin(angle) * exp(- (r.x - x0) * (r.x - x0) / sigma /sigma - (r.y - y0) * (r.y - y0) / sigma / sigma - r.z * r.z / sigma / sigma) * cos(2 * pi * (cos(angle) * r.x + sin(angle) * r.y));};
    auto ey_func = [sigma, x0, y0, angle] (cvector3d r) -> double {return cos(angle) * exp(- (r.x - x0) * (r.x - x0) / sigma /sigma - (r.y - y0) * (r.y - y0) / sigma / sigma - r.z * r.z / sigma / sigma) * cos(2 * pi * (cos(angle) * r.x + sin(angle) * r.y));};

    initialize_fields(ex_func, ey_func, empty_func, empty_func, empty_func, ey_func);
    cout << "After init" << endl;

    while (current_iteration * dt <= tmax) {
        if (current_iteration * dt >= (output_iteration * output_dt)) {
            output_fields(output_iteration);

            output_iteration++;
        }
        advance_b(0.5 * dt);
        advance_e(dt);
        advance_b(0.5 * dt);

        current_iteration++;
    }

	return 0;
}
