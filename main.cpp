#include <iostream>
#include <algorithm>
#include <iomanip> 
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <ios>
#include <cstring>
#include <stdexcept>

using namespace std;

const int maximum = INT16_MAX;
const float minValue = -2.0; // a normalized vector ranges from -1 to 1.
const float epsilon = 0.000001f; //for float point comparison


const string INPUT_FILE("input2.txt");
const string OUTPUT_FILE("output.txt");

/**
 * Method that is created to solely print the memo look-up table.
 * 
 * \param memo 2D array containing the memoized data
 * \param m the number of rows in the 2D look-up array.
 * \param n the number of columns in the 2D look-up array.
 * \param precision the precision of each element of in the look-up array.
 */
void PrintMemo(float** memo, int m, int n, int precision) {
	for (int i = 0; i < m + 1; i++) {
		for (int j = 0; j < n + 1; j++) {
			cout << fixed << setprecision(precision) << setw(precision * 2) << memo[i][j] << " | ";
		}
		cout << endl;
		cout.unsetf(ios_base::floatfield);
	}
	cout << endl << endl;
}

/**
 * Method that compares two floating points and checks to see if they are the same.
 * 
 * \param a floating point 1
 * \param b floating point 2
 * \return Returns true if the 2 floats have the same value.
 */
bool FloatCompare(float a, float b) {
	return fabs(a - b) < epsilon;
}

/**
 * Method that frees the memory that was dynamically assigned to a 2D matrix.
 * 
 * \param memo 2D matrix that was dynamically created
 * \param m number of rows in the 2D matrix
 * \param n number of columns in the 2D matrix
 */
void FreeMemo(float** memo, int m, int n) {
	for (int i = 0; i < m; i++)
		free(memo[i]);
	free(memo);
}

/**
 * Method that calculates the best possible subsequence between two input sequences against a given target
 * sequence via recursion ONLY. NOTE: The sequences are essentially normalized vectors. 
 * 
 * \param p Normalized input sequence.
 * \param q Normalized input sequence.
 * \param t Normalized target sequence.
 * \param m Number of elements in the p sequence.
 * \param n Number of elements in the q sequence.
 * \return The highest alignment score.
 */
float AlignmentScore_Recursive(float* p, float* q, float* t, int m, int n) {
	int i = m - 1; //index used for p sequence
	int j = n - 1; //index used for q sequence
	int k = m + n - 1; //index used for the t (target) sequence.

	if (m == 0 && n == 1) {
		return q[j] * t[k]; //Base Case
	}
	else if (m == 1 && n == 0) {
		return p[i] * t[k]; //Base Case
	}
	else if (m == 0 && n > 1) {
		//Look at the right of the recursion tree
		return AlignmentScore_Recursive(p, q, t, m, n - 1) + (t[k]*q[j]); 
	}
	else if (n == 0 && m > 1) {
		//Look at the left of the recursion tree
		return AlignmentScore_Recursive(p, q, t, m - 1, n) + (t[k] * p[i]); 
	}
	else {
		//Compare both the left and right side of the same level of the tree
		return max(AlignmentScore_Recursive(p, q, t, m - 1, n) + (t[k] * p[i]) , 
			AlignmentScore_Recursive(p, q, t, m, n - 1) + (t[k] * q[j]));
	}
}

/**
 * Method that dynamically allocates memory for a 2D array and initializes every element in said 2D matrix.
 * 
 * \param m number of rows in a 2D matrix
 * \param n number of columns in a 2D matrix
 * \param minvalue value assigned to every element in the 2D matrix
 * \return 
 */

float** InitializeMemo(int m, int n, float minvalue) {
	float** memo;
	memo = (float**)malloc((m + 1) * sizeof(float*));
	for (int i = 0; i < m + 1; i++) {
		memo[i] = (float*)malloc((n + 1) * sizeof(float));
	}

	for (int i = 0; i < m + 1; i++) {
		for (int j = 0; j < n + 1; j++) {
			memo[i][j] = minvalue;
		}
	}
	return memo;
}

/**
 * Method calculate the best possible subsequence between two input sequences against a given target
 * sequence via recursion dynamic programming. NOTE: The sequences are essentially normalized vectors. 
 * 
 * \param p Normalized input sequence.
 * \param q Normalized input sequence.
 * \param t Normalized target sequence.
 * \param m Number of elements in the p sequence.
 * \param n Number of elements in the q sequence.
 * \param memo 2D look up table that stores the best values of each recursion path.
 * \return The highest alignment score.
 */
float AlignmentScore_Memoized(float* p, float* q, float* t, int m, int n, float** memo) {

	int i = m - 1; //index used for p sequence
	int j = n - 1; //index used for q sequence
	int k = m + n - 1; //index used for the t (target) sequence.


	if (FloatCompare(memo[m][n],minValue)) {
		if (m == 0 && n == 1) {
			memo[m][n] = q[j] * t[k]; //Base Case
		}
		else if (m == 1 && n == 0) {

			memo[m][n] = p[i] * t[k]; //Base Case
		}
		else if (m == 0 && n > 1) {
			//Look at the right of the recursion tree ONLY.
			memo[m][n] = AlignmentScore_Memoized(p, q, t, m, n - 1, memo) + (t[k] * q[j]); 
		}
		else if (n == 0 && m > 1) {
			//Look at the left of the recursion tree ONLY.
			memo[m][n] = AlignmentScore_Memoized(p, q, t, m - 1, n, memo) + (t[k] * p[i]); 
		}
		else {
			//Compare both the left and right side of the same level of the tree
			memo[m][n] = max(AlignmentScore_Memoized(p, q, t, m - 1, n, memo) + (t[k] * p[i]),
				AlignmentScore_Memoized(p, q, t, m, n - 1, memo) + (t[k] * q[j]));
		}
	}
		return memo[m][n];
}

/**
 * Method calculate the best possible subsequence between two input sequences against a given target
 * sequence via ITERATIVE dynamic programming. NOTE: The sequences are essentially normalized vectors.
 * \param p Normalized input sequence.
 * \param q Normalized input sequence.
 * \param t Normalized target sequence.
 * \param m Number of elements in the p sequence.
 * \param n Number of elements in the q sequence.
 * \param memo 2D look up table that stores the best values of each recursion path.
 * \return The highest alignment score.
 */
float AlignmentScore_Iterative(float* p, float* q, float* t, int m, int n, float** memo) {
	//initialize the first element of the array 
	memo[0][0] = 0;

	//Look at the top and left child of the parent and assign the parent
	//max value of the children. 
	for (int i = 1; i <= m; i++) {
		//Base Case 1: Initialize the first row of the 2D look up table
		memo[i][0] = memo[i - 1][0] + t[i - 1] * p[i - 1];
		for (int j = 1; j <= n; j++) {
			if (i == 1) {
				//Base Case 2: Initialize the first column of the 2D look up table
				memo[0][j] = memo[0][j - 1] + t[j - 1] * q[j - 1];
			}
			memo[i][j] = max(memo[i - 1][j] + (t[i + j - 1] * p[i - 1]), 
				memo[i][j - 1] + (t[i + j - 1] * q[j - 1]));
		}
	}

	return memo[m][n];

}

/**
 * Method that accepts a generated memo table and outputs a combination subsequences of the two given input sequences.
 * 
 * \param memo 2D lookup table that contains all the possible paths the program took to find the highest alignment score
 * \param p Normalized input sequence.
 * \param q Normalized input sequence.
 * \param t Normalized target sequence.
 * \param m Number of elements in the p sequence.
 * \param n Number of elements in the q sequence.
 * \return The combination subsequence of the two input sequences with the highest alignment score. 
 */
float* GetBestSequence(float** memo, float* p, float* q, float* t, int m, int n) {
	int i = m, j = n, index = m + n - 1;
	float* best = (float*)malloc((m + n) * sizeof(float));

	while (index >= 0) {
		if (i == 0 && j > 0) {
			best[index--] = q[--j];
		}

		else if (j == 0 && i > 0) {
			best[index--] = p[--i];
		}

		else if (FloatCompare(memo[i][j] - p[i - 1] * t[index], memo[i - 1][j])) {
			best[index--] = p[--i];
		}

		else if (FloatCompare(memo[i][j] - q[j - 1] * t[index], memo[i][j - 1])) {
			best[index--] = q[--j];
		}
		
		else {
			throw new invalid_argument("Error in best sequence operation!");
		}
	}
	return best;
}




int main(int argc, char** argv) {

	if (argc > 2) {
		cout << "Too many inputs. Try again." << endl;
		exit(EXIT_FAILURE);
	}

	ifstream in(INPUT_FILE); //load input file into the stream
	ofstream out(OUTPUT_FILE); //load output file into the stream
	string line; //string that takes one line of the file
	istringstream iss; //buffer stream

	int m, n;
	float *seq1, *seq2, *target;
	float** memo;

	if (!in.is_open()) {
		cout << "Error opening " << INPUT_FILE << endl;
		return 1;
	}

	if (!out.is_open()) {
		cout << "Error opening " << OUTPUT_FILE << endl;
		return 1;
	}



	//***********LOAD INPUT FILE************************//
	getline(in, line);
	iss.str(line);
	iss >> m >> n;
	iss.clear();

	seq1 = (float*)malloc(m * sizeof(float));
	seq2 = (float*)malloc(n * sizeof(float));
	target = (float*)malloc((m + n) * sizeof(float));
	getline(in, line);
	iss.str(line);
	for (int i = 0; i < m; i++) {
		iss >> seq1[i];
	}
	iss.clear(); //clear buffer


	getline(in, line);
	iss.str(line);
	for (int i = 0; i < n; i++) {
		iss >> seq2[i];
	}
	iss.clear(); //clear buffer

	getline(in, line);
	iss.str(line);
	for (int i = 0; i < (m + n); i++) {
		iss >> target[i];
	}
	iss.clear(); //clear buffer
	in.close();
	//*******************************************************//
	if (argc == 1 || (argc == 2 && (strcmp(argv[1], "-m") == 0))) {
		memo = InitializeMemo(m, n, minValue);
		out << AlignmentScore_Memoized(seq1, seq2, target, m, n, memo) << endl;
		float* bestseq2 = GetBestSequence(memo, seq1, seq2, target, m, n);
		for (int i = 0; i < (m + n); i++) {
			out << bestseq2[i] << " ";
		}
		out << endl;
		out.close();
		free(bestseq2);
		FreeMemo(memo, m, n);
		cout << "Sequence obtained via. memoized recursion dynamic programming!" << endl;
	}

	else if (argc == 2 && (strcmp(argv[1], "-i") == 0)) {
		memo = InitializeMemo(m, n, minValue);
		out <<  AlignmentScore_Iterative(seq1, seq2, target, m, n, memo) << endl;
		float* bestseq3 = GetBestSequence(memo, seq1, seq2, target, m, n);
		for (int i = 0; i < (m + n); i++) {
			out << bestseq3[i] << " ";
		}
		out << endl;
		out.close();
		free(bestseq3);
		FreeMemo(memo, m, n);
		cout << "Sequence obtained via. iterative dynamic programming!" << endl;
	}
	
	else if (argc == 2 && (strcmp(argv[1], "-r") == 0)) {
		cout << "Best alignment score via pure recursion is: " << AlignmentScore_Recursive(seq1, seq2, target, m, n) << endl;
	}
	else {
		cout << "Entered incorrect options." << endl;
	}

	free(seq1);
	free(seq2);
	free(target);
	return 0;
}
