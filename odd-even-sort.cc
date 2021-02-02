#include <iostream>
#include <mpi.h>
#include <math.h>
#include <iomanip>
#include <boost/sort/spreadsort/spreadsort.hpp>
using namespace std;

float calMaxVal(float data[], int length)
{
	float temp = data[0];
	for(int i = 0;i < length;i++)
	{
		if(data[i] > temp)
		{
			temp = data[i];
		}
	}
	return temp;
}

int findPartner(int phase, int myRank, int size, int fsize)
{
	int partnerRank;

	if(phase % 2 == 0)
	{
		partnerRank = (myRank % 2 == 0)? myRank + 1: myRank - 1;
	}
	else
	{
		partnerRank = (myRank % 2 == 0)? myRank - 1: myRank + 1;
	}
	if(fsize < size)
		return (partnerRank < 0 || partnerRank >= fsize )? -1: partnerRank;
	else
		return (partnerRank < 0 || partnerRank >= size )? -1: partnerRank;	
}

void localSort(float *data, int length)
{
	boost::sort::spreadsort::spreadsort(data, data + length);	
}


void oddEvenSort(int phase, int myRank, int partnerRank, float *myData, int dataLength)
{
	float *mergedData;
	float *recvBuffer;
	//static double commTime;
	//double commStart, commEnd;
	MPI_Status status;
	if(myRank > partnerRank)
	{
		//commStart = MPI_Wtime();
		MPI_Send(myData, dataLength, MPI_FLOAT, partnerRank, 0, MPI_COMM_WORLD);
		//commEnd = MPI_Wtime();
		//commTime += commEnd - commStart;
	}
	else
	{
		mergedData = new float[dataLength * 2];
		recvBuffer = new float[dataLength];
		
		//commStart = MPI_Wtime();
		MPI_Recv(recvBuffer, dataLength, MPI_FLOAT, partnerRank, 0 ,MPI_COMM_WORLD, &status);
		//commEnd = MPI_Wtime();
		//commTime += commEnd - commStart;

		int i = 0, j = 0, k = 0;
		for(;i < dataLength && j < dataLength;)
		{
			if(myData[i] < recvBuffer[j])
			{
				mergedData[k] = myData[i];
				i++;
				k++;
			}
			
			else
			{
				mergedData[k] = recvBuffer[j];
				j++;
				k++;
			}
		
		}
		if(i == dataLength)
		{
			for(;j < dataLength;)
			{
				mergedData[k] = recvBuffer[j];
				j++;
				k++;
			}
		}
		if(j == dataLength)
		{
			for(;i < dataLength;)
			{
				mergedData[k] = myData[i];
				i++;
				k++;
			}
		}
		delete []recvBuffer;
	}

	if(myRank < partnerRank)
	{
		for(int i=0; i < dataLength; i++)
		{
			myData[i] = mergedData[i];
		}
		float *dividePoint;
	    dividePoint = mergedData + dataLength;

		//commStart = MPI_Wtime();
		MPI_Send(dividePoint, dataLength, MPI_FLOAT, partnerRank, 99, MPI_COMM_WORLD);
		//commEnd = MPI_Wtime();
		//commTime += commEnd - commStart;

		delete []mergedData;
	}
	else
	{
		//commStart = MPI_Wtime();
		MPI_Recv(myData, dataLength, MPI_FLOAT, partnerRank, 99, MPI_COMM_WORLD, &status);
		//commEnd = MPI_Wtime();
		//commTime += commEnd - commStart;
	}
	//std::cout << "Comm: " << myRank << " : " << commTime << std::endl;
}

int main(int argc, char** argv) 
{
	MPI_Init(&argc,&argv);
	double start, end, ioStart, ioEnd, ioTime;
	start = MPI_Wtime();
	int rank, size, total, numOfRead, numOfComplement, numOfWrite;
	bool needComplement = false;
	MPI_Offset fsize;
	MPI_Status status;
	total = atoi(argv[1]);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_File f_in, f_out;
	MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_RDONLY, MPI_INFO_NULL, &f_in);
	MPI_File_open(MPI_COMM_WORLD, argv[3], MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &f_out);
	//MPI_File_get_size(f_in, &fsize);
	fsize = total; //= sizeof(float);
	long stride = (long)ceil((long double)fsize / (long double)size);
	float *data = new float[stride];
	float maxVal;
	float *maxVals;
	if(fsize % size != 0 && fsize > size)
		needComplement = true;	
	if(rank == size - 1 && needComplement)
	{
		numOfRead = fsize % stride;
		numOfComplement = stride - numOfRead;
	}
	else
	{
		numOfRead = stride;
	}
	if(fsize < size)
	{
		if(rank >= fsize)
		{
			numOfRead = 0;
		}
	}
	//ioStart = MPI_Wtime();
	MPI_File_read_at(f_in, sizeof(float) * rank * stride, data, numOfRead, MPI_FLOAT, MPI_STATUS_IGNORE);
	//ioEnd = MPI_Wtime();
	//ioTime += ioEnd - ioStart;

	maxVal = calMaxVal(data, numOfRead);

	if(needComplement)
	{
		maxVals = new float[size];
		MPI_Gather(&maxVal, 1, MPI_FLOAT, maxVals, 1, MPI_FLOAT, size - 1 , MPI_COMM_WORLD);
	}


	if(rank == size - 1 && needComplement)
	{
		maxVal = maxVals[0];
		for(int i = 0; i < size; i++)
		{
			if(maxVal < maxVals[i])
			{
				maxVal = maxVals[i];
			}		
		}

		for(int i = fsize % stride; i < stride; i++)
		{
			data[i] = maxVal;
		}

		delete []maxVals;
	}

	// Processes do the local sort respectively.
	if(rank < fsize)
		localSort(data, stride);
	
	// Processes do the odd-even sort"
	for(int i = 0; i < size; i++)
	{
		int myPartnerRank = findPartner(i, rank, size, fsize);
		if(myPartnerRank != -1 && rank < fsize)
		{	
			oddEvenSort(i, rank, myPartnerRank, data, stride);	
		}
	}
	
	numOfWrite = numOfRead;

	//ioStart = MPI_Wtime();
	MPI_File_write_at(f_out, sizeof(float) * rank * stride, data, numOfWrite, MPI_FLOAT, MPI_STATUS_IGNORE);
	//ioEnd = MPI_Wtime();
	//ioTime += ioEnd - ioStart;

	//std::cout << "io: " <<"rank: " << rank << " : " << ioTime << std::endl;
	end = MPI_Wtime();
	std::cout << "total: " <<"rank: " << rank << " : " << end - start << std::endl;

	delete []data;
	MPI_File_close(&f_in);
	MPI_File_close(&f_out);
	MPI_Finalize();
}
