/*                                                                      
 * GPL3 License                                                                         
 *                                                                      
 *                                                                      
 * Copyright (c) 2022 Ransome CA                              
 *                                                                      
 * This file is part of SCIARA-fv2-CUDA.                                          
 *                                                                      
 * SCIARA-fv2-CUDA is free software: you can redistribute it and/or modify        
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation, either version 3 of the License, or    
 * (at your option) any later version.                                  
 *                                                                      
 * SCIARA-fv2-CUDA is distributed in the hope that it will be useful,             
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        
 * GNU General Public License for more details.                         
 *                                                                      
 * You should have received a copy of the GNU General Public License    
 * along with SCIARA-fv2-CUDA.  If not, see <http://www.gnu.org/licenses/>.       
 */       

#include "Sciara.cuh"
#include "cal2DBuffer.cuh"


const int _Xi[] = {0, -1,  0,  0,  1, -1,  1,  1, -1}; // Xj: Moore neighborhood row coordinates (see below)
const int _Xj[] = {0,  0, -1,  1,  0, -1, -1,  1,  1}; // Xj: Moore neighborhood col coordinates (see below)


template <typename T>
static void __malloc(T** ptr, size_t size)
{
  *ptr = (T*) malloc(size);
  if (*ptr == NULL)
  {
    printf("Error: malloc failed\n");
    exit(1);
  }
}

void allocateSubstates(Sciara *sciara) {

  __malloc(&sciara->substates->Sz,       sciara->domain->rows * sciara->domain->cols * sizeof(double));
  __malloc(&sciara->substates->Sz_next,  sciara->domain->rows * sciara->domain->cols * sizeof(double));
  __malloc(&sciara->substates->Sh,       sciara->domain->rows * sciara->domain->cols * sizeof(double));
  __malloc(&sciara->substates->Sh_next,  sciara->domain->rows * sciara->domain->cols * sizeof(double));
  __malloc(&sciara->substates->ST,       sciara->domain->rows * sciara->domain->cols * sizeof(double));
  __malloc(&sciara->substates->ST_next,  sciara->domain->rows * sciara->domain->cols * sizeof(double));
  __malloc(&sciara->substates->Mf,       sciara->domain->rows * sciara->domain->cols * sizeof(double) * NUMBER_OF_OUTFLOWS);
  __malloc(&sciara->substates->Mb,       sciara->domain->rows * sciara->domain->cols * sizeof(bool));
  __malloc(&sciara->substates->Mhs,      sciara->domain->rows * sciara->domain->cols * sizeof(double));

  memset(sciara->substates->Sz,       0, sciara->domain->rows * sciara->domain->cols * sizeof(double));
  memset(sciara->substates->Sz_next,  0, sciara->domain->rows * sciara->domain->cols * sizeof(double));
  memset(sciara->substates->Sh,       0, sciara->domain->rows * sciara->domain->cols * sizeof(double));
  memset(sciara->substates->Sh_next,  0, sciara->domain->rows * sciara->domain->cols * sizeof(double));
  memset(sciara->substates->ST,       0, sciara->domain->rows * sciara->domain->cols * sizeof(double));
  memset(sciara->substates->ST_next,  0, sciara->domain->rows * sciara->domain->cols * sizeof(double));
  memset(sciara->substates->Mf,       0, sciara->domain->rows * sciara->domain->cols * sizeof(double) * NUMBER_OF_OUTFLOWS);
  memset(sciara->substates->Mb,       0, sciara->domain->rows * sciara->domain->cols * sizeof(bool));
  memset(sciara->substates->Mhs,      0, sciara->domain->rows * sciara->domain->cols * sizeof(double));

}

void deallocateSubstates(Sciara *sciara) {

  if(sciara->substates->Sz != NULL) {
    cudaFree(sciara->substates->Sz);
  }

  if(sciara->substates->Sz_next != NULL) {
    cudaFree(sciara->substates->Sz_next);
  }

  if(sciara->substates->Sh != NULL) {
    cudaFree(sciara->substates->Sh);
  }

  if(sciara->substates->Sh_next != NULL) {
    cudaFree(sciara->substates->Sh_next);
  }

  if(sciara->substates->ST != NULL) {
    cudaFree(sciara->substates->ST);
  }

  if(sciara->substates->ST_next != NULL) {
    cudaFree(sciara->substates->ST_next);
  }

  if(sciara->substates->Mf != NULL) {
    cudaFree(sciara->substates->Mf);
  }

  if(sciara->substates->Mb != NULL) {
    cudaFree(sciara->substates->Mb);
  }

  if(sciara->substates->Mhs != NULL) {
    cudaFree(sciara->substates->Mhs);
  }

}


void evaluatePowerLawParams(double PTvent, double PTsol, double value_sol, double value_vent, double &k1, double &k2) {
    k2 = (log10(value_vent) - log10(value_sol) ) / (PTvent - PTsol);
    k1 = (log10(value_sol) - k2 * (PTsol));
}

void simulationInitialize(Sciara* sciara) {

    unsigned int maximum_number_of_emissions = 0;

    sciara->simulation->step = 0;
    sciara->simulation->elapsed_time = 0;

    for (unsigned int i = 0; i < sciara->simulation->emission_rate.size(); i++) {

        if (maximum_number_of_emissions < sciara->simulation->emission_rate[i].size()) {
            maximum_number_of_emissions = sciara->simulation->emission_rate[i].size();
        }

    }

    sciara->simulation->effusion_duration = sciara->simulation->emission_time * maximum_number_of_emissions;
    sciara->simulation->total_emitted_lava = 0;

    MakeBorder(sciara);

    evaluatePowerLawParams(
        sciara->parameters->PTvent, 
        sciara->parameters->PTsol, 
        sciara->parameters->Pr_Tsol,  
        sciara->parameters->Pr_Tvent,  
        sciara->parameters->a, 
        sciara->parameters->b);
    evaluatePowerLawParams(
        sciara->parameters->PTvent,
        sciara->parameters->PTsol,
        sciara->parameters->Phc_Tsol,
        sciara->parameters->Phc_Tvent,
        sciara->parameters->c,
        sciara->parameters->d);

}


void init(Sciara*& sciara) {

    sciara = new Sciara;
    sciara->domain = new Domain;

    sciara->X = new NeighsRelativeCoords;
    sciara->X->Xi = new int[MOORE_NEIGHBORS];
    sciara->X->Xj = new int[MOORE_NEIGHBORS];

    for (size_t n = 0; n < MOORE_NEIGHBORS; n++) {
        sciara->X->Xi[n] = _Xi[n];
        sciara->X->Xj[n] = _Xj[n];
    }

    sciara->substates = new Substates;
    sciara->parameters = new Parameters;
    sciara->simulation = new Simulation;

}



void finalize(Sciara*& sciara) {

    deallocateSubstates(sciara);
    delete sciara->domain;
    delete sciara->X->Xi;
    delete sciara->X->Xj;
    delete sciara->X;
    delete sciara->substates;
    delete sciara->parameters;
    delete sciara->simulation;
    delete sciara;
    sciara = NULL;

}


void MakeBorder(Sciara *sciara) {
	int j, i;

	//prima riga
	i = 0;
	for (j = 0; j < sciara->domain->cols; j++)
		if (calGetMatrixElement(sciara->substates->Sz, sciara->domain->cols, i, j) >= 0)
			calSetMatrixElement(sciara->substates->Mb, sciara->domain->cols, i, j, true);

	//ultima riga
	i = sciara->domain->rows - 1;
	for (j = 0; j < sciara->domain->cols; j++)
		if (calGetMatrixElement(sciara->substates->Sz, sciara->domain->cols, i, j) >= 0)
			calSetMatrixElement(sciara->substates->Mb, sciara->domain->cols, i, j, true);

	//prima colonna
	j = 0;
	for (i = 0; i < sciara->domain->rows; i++)
		if (calGetMatrixElement(sciara->substates->Sz, sciara->domain->cols, i, j) >= 0)
			calSetMatrixElement(sciara->substates->Mb, sciara->domain->cols, i, j, true);
  
	//ultima colonna
	j = sciara->domain->cols - 1;
	for (i = 0; i < sciara->domain->rows; i++)
		if (calGetMatrixElement(sciara->substates->Sz, sciara->domain->cols, i, j) >= 0)
			calSetMatrixElement(sciara->substates->Mb, sciara->domain->cols, i, j, true);
	
	//il resto
	for (int i = 1; i < sciara->domain->rows - 1; i++)
		for (int j = 1; j < sciara->domain->cols - 1; j++)
			if (calGetMatrixElement(sciara->substates->Sz, sciara->domain->cols, i, j) >= 0) {
				for (int k = 1; k < MOORE_NEIGHBORS; k++)
					if (calGetMatrixElement(sciara->substates->Sz, sciara->domain->cols, i+sciara->X->Xi[k], j+sciara->X->Xj[k]) < 0)
          {
			      calSetMatrixElement(sciara->substates->Mb, sciara->domain->cols, i, j, true);
						break;
					}
			}
}
