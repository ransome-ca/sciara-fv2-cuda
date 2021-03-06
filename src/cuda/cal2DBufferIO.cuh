/*
 * Copyright (c) 2016 OpenCALTeam (https://github.com/OpenCALTeam),
 * Telesio Research Group,
 * Department of Mathematics and Computer Science,
 * University of Calabria, Italy.
 *
 * This file is part of OpenCAL (Open Computing Abstraction Layer).
 *
 * OpenCAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * OpenCAL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with OpenCAL. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef cal2DBufferIO_h
#define cal2DBufferIO_h

//#include <OpenCAL/calCommon.h>
#include <stdio.h>



/*! \brief Loads a byte matrix from file.
*/
void calfLoadMatrix2Db(bool* M, int rows, int columns, FILE* f);

/*! \brief Loads an int matrix from file.
*/
void calfLoadMatrix2Di(int* M, int rows, int columns, FILE* f);

/*! \brief Loads a real (floating point) matrix from file.
*/
void calfLoadMatrix2Dr(double* M, int rows, int columns, FILE* f);



/*! \brief Loads a byte matrix from file.
*/
bool calLoadMatrix2Db(bool* M, int rows, int columns, char* path);

/*! \brief Loads an int matrix from file.
*/
bool calLoadMatrix2Di(int* M, int rows, int columns, char* path);

/*! \brief Loads a real (floating point) matrix from file.
*/
bool calLoadMatrix2Dr(double* M, int rows, int columns, char* path);



/*! \brief Saves a byte matrix to file.
*/
void calfSaveMatrix2Db(bool* M, int rows, int columns, FILE* f);

/*! \brief Saves an int matrix to file.
*/
void calfSaveMatrix2Di(int* M, int rows, int columns, FILE* f);

/*! \brief Saves a real (floating point) matrix to file.
*/
void calfSaveMatrix2Dr(double* M, int rows, int columns, FILE* f);



/*! \brief Saves a byte matrix to file.
*/
bool calSaveMatrix2Db(bool* M, int rows, int columns, char* path);

/*! \brief Saves a int matrix to file.
*/
bool calSaveMatrix2Di(int* M, int rows, int columns, char* path);

/*! \brief Saves a real (floating point) matrix to file.
*/
bool calSaveMatrix2Dr(double* M, int rows, int columns, char* path);



#endif
