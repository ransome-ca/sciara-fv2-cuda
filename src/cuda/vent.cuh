#ifndef vent_h
#define vent_h

#include <stdio.h>
#include <vector>
using std::vector;

#define EMISSION_RATE_FILE_ERROR     0
#define EMISSION_RATE_FILE_OK        1

//---------------------------------------------------------------------------
class TEmissionRate {
public:
	int vent_id() const {
		return _vent_id;
	}
	void set_vent_id(int id) {
		_vent_id = id;
	}
	double& operator[](int i) {
		return _emission_rate[i];
	}
	vector<double>& emission_rate() {
		return _emission_rate;
	}
	const vector<double>& emission_rate() const {
		return _emission_rate;
	}
	unsigned int size() {
		return _emission_rate.size();
	}
	void print() {
		printf("vent_id = %d\nemission_rate\n", _vent_id);
		for (unsigned int i = 0; i < _emission_rate.size(); i++)
			printf("%f\n", _emission_rate[i]);
	}
private:
	int _vent_id;
	vector<double> _emission_rate;
};
//---------------------------------------------------------------------------


class TVent {
	public:

		__host__ __device__
		inline const auto& x() const {
			return this->__x;
		}

		__host__ __device__
		inline const auto& y() const {
			return this->__y;
		}

		__host__ __device__
		inline const auto& vent_id() const {
			return this->__vent_id;
		}

		__host__ __device__
		inline void set_x(double x) {
			__x = x;
		}

		__host__ __device__
		inline void set_y(double y) {
			__y = y;
		}

		__host__ __device__
		inline void set_vent_id(int id) {
			__vent_id = id;
		}


		__host__
		bool setEmissionRate(vector<TEmissionRate> er_vec, int id) {
			
			bool found = false;

			for(const auto& i : er_vec) {
				if(i.vent_id() == id) {

					cudaMallocManaged(&__emission_rate, i.emission_rate().size() * sizeof(double));

					size_t k = 0;
					for(const auto& j : i.emission_rate()) {
						__emission_rate[k++] = j;
					}

					__emission_rate_size = i.emission_rate().size();

					found = true;
				}
			}

			return found;

		}

		__host__ __device__
		const auto& operator[](int i) const {
			return __emission_rate[i];
		}

		__host__ __device__
		size_t size() const {
			return __emission_rate_size;
		}

		__host__ __device__
		inline double thickness(double sim_elapsed_time, double Pt, unsigned int emission_time, double Pac) const {

			size_t i = (sim_elapsed_time / emission_time);

			if((sim_elapsed_time / emission_time) >= size()) {
				return 0;
			} else {
				return __emission_rate[i] / Pac * Pt;
			}
		
		}


		__host__
		inline void print() const {
			
		}


	private:
		size_t __x;
		size_t __y;
		size_t __vent_id;
		size_t __emission_rate_size;
		double* __emission_rate;

};


//---------------------------------------------------------------------------
/*extern unsigned int emission_time;
 extern vector<TEmissionRate> emission_rate;
 extern vector<TVent> vent;*/

void initVents(int* Mv, int lx, int ly, vector<TVent>& vent);
void addVent(int x, int y, int vent_id, vector<TVent>& vent);
void removeLastVent(vector<TVent>& vent);
int loadEmissionRates(FILE *f, unsigned int& emission_time, vector<TEmissionRate>& er_vec, vector<TVent> &vent);
int loadOneEmissionRates(FILE *f, unsigned int vent_id, vector<TEmissionRate>& er_vec);
int defineVents(const vector<TEmissionRate>& emission_rate, vector<TVent>& vent);

void rebuildVentsMatrix(int* Mv, int lx, int ly, vector<TVent>& vent);
void saveEmissionRates(FILE *f, unsigned int emission_time, vector<TEmissionRate>& er_vec);
void printEmissionRates(unsigned int emission_time, vector<TEmissionRate>& er_vec);
//---------------------------------------------------------------------------

#endif
