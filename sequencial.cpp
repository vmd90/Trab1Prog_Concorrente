#include <time.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <math.h>

void JacobiRichardson(float **A, int N, float *b, float erro, int *iteracoes, float *solucao);
bool LerArquivo(const std::string& arq, float **&A, float *&b, int &linha, int &N, float &erro, int &iteracoes);
float CalcularLinha(float **A, int N, float *solucao, int linha);

int main(int argc, char *argv[])
{
    if(argc < 2) {
        printf("Numero de argumentos insuficiente\n");
        return 1;
    }
	float **A; /* matriz */
	float *b; /* vetor */
	int N; /* tamanho do vetor b */
	int linha;
	float erro;
	int iteracoes; /* numero de iteracoes */
	float *solucao; /* vetor da solucao */
	int i;
	std::string arq = argv[1]; // nome do arquivo
	std::string arq_saida = std::string("out_") + arq; // nome do arquivo de saida
	std::ofstream out(arq_saida.c_str(), std::ios::app);

	if( !LerArquivo(arq, A, b, linha, N, erro, iteracoes) ) {
		printf("Erro ao ler arquivo.\n");
		return 1;
	}
	if (!out.is_open()) {
		printf("Erro ao criar arquivo.\n");
		return 1;
	}

	solucao = new float[N];

	time_t t1 = time(NULL);
	JacobiRichardson(A, N, b, erro, &iteracoes, solucao);
	time_t t2 = time(NULL);
	out << "\nTempo: " << t2 - t1 << "s\n";

	float r = CalcularLinha(A, N, solucao, linha);

	out << "-------------------------------\nIterations: "<< iteracoes << "\n";
	out << "RowTest: " << linha << " => [" << r << "] =? " << b[linha] << "\n-------------------------------\n";
	out.close();

	delete[] b;
	delete[] solucao;
	for(i = 0; i < N; ++i)
		delete[] A[i];
	delete[] A;
	return 0;
}

// criterio das linhas
int Linhas(float **A, int N)
{
	float max = 0.0f;
	float soma = 0.0f;
	int i;
	int j;
	for (i = 0; i < N; ++i)
	{
		soma = 0.0f;
		for (j = 0; j < N; ++j)
		{
			if (i != j)	{
				soma += fabs(A[i][j]); /* soma para cada linha */
			}
		}
		soma /= fabs(A[i][i]);
		if (max < soma){
			max = soma;
		}
	}

	if (max < 1)
		return 1;
	else
		return 0;
}

// criterio das colunas
int Colunas(float **A, int N)
{
	float max = 0.0f;
	float soma = 0.0f;
	int i;
	int j;
	for (i = 0; i < N; ++i)
	{
		soma = 0.0f;
		for (j = 0; j < N; ++j)
		{
			if (i != j) {
				soma += fabs(A[j][i]); /* soma para cada coluna */
			}
		}
		soma /= fabs(A[i][i]);
		if (max < soma){  /* recebe a maior soma */
			max = soma;
		}
	}

	if (max < 1)
		return 1;
	else
		return 0;
}

// Leitura de arquivo matriz
bool LerArquivo(const std::string& arq, float **&A, float *&b, int &linha, int &N, float &erro, int &iteracoes)
{
	std::ifstream file(arq.c_str());
	if(!file.is_open())
		return false;
	file >> N;
	file >> linha;
	file >> erro;
	file >> iteracoes;
	// Lendo a matriz
	A = new float*[N];
	for(int i = 0; i < N; ++i)
		A[i] = new float[N];

	for(int i = 0; i < N; ++i)
		for(int j = 0; j < N; ++j)
			file >> A[i][j];

	// Lendo o vetor
	b = new float[N];
	for(int i = 0; i < N; ++i)
		file >> b[i];
	file.close();
	return true;
}

// metodo Jacobi-Richardson
void JacobiRichardson(float **A, int N, float *b, float erro, int *iteracoes, float *solucao)
{
	float *x, *y;
	int contador = 0;
	float numerador;
	float denominador;
	int i,j;

	x = new float[N];
	y = new float[N];
	/* x comeca com 0 */
	for(i = 0; i < N; i++)
		x[i] = 0.0f;

	if( Linhas(A, N) || Colunas(A, N) ) {
		do
		{
			contador++;
			numerador = 0.0f;
			denominador = 0.0f;
			for(i = 0; i < N; ++i)
			{
				y[i] = 0; /* x(k+1) */
				for(j = 0; j < N; ++j)
				{
					if(i != j) {
						y[i] += (A[i][j] * x[j]);
					}
				}
				y[i] = (1/A[i][i]) * (b[i] - y[i]);
				if(numerador < fabs(y[i] - x[i]))
					numerador = fabs(y[i] - x[i]);
				if(denominador < fabs(y[i]))
					denominador = y[i];
			}
			for(i = 0; i < N; ++i) {
				x[i] = y[i];
			}

		} while( (numerador/denominador) > erro && contador < *iteracoes );

		/* vetor de saida */
		for(i = 0; i < N; ++i) {
			solucao[i] = x[i];
		}
		*iteracoes = contador;
	}
	else {
		printf("Criterio de linhas ou colunas nao satisfeito!\n");
	}
	delete[] x;
	delete[] y;
}

// Calculo da linha da matriz
float CalcularLinha(float **A, int N, float *solucao, int linha)
{
	float r = 0;
	for (int i = 0; i < N; ++i)
	{
		r += A[linha][i] * solucao[i];
	}
	return r;
}
