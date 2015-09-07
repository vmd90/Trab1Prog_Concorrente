#include <time.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include <pthread.h>

//void JacobiRichardson(float **A, int N, float *b, float erro, int *iteracoes, float *solucao);
void *JacobiRichardson(void *arg);
bool LerArquivo(const std::string& arq, float **&A, float *&b, int &linha, int &N, float &erro, int &iteracoes);
float CalcularLinha(float **A, int N, float *solucao, int linha);

struct Dados
{
	float **A; /* matriz */
	float *b; /* vetor */
	int N; /* tamanho do vetor b */
	float erro;
	int iteracoes; /* numero de iteracoes */
	float *solucao; /* vetor da solucao */
};


int main(int argc, char *argv[])
{
    if(argc < 2) {
        printf("Numero de argumentos insuficiente\n");
        return 1;
    }

	std::string arq = argv[1]; // nome do arquivo
	std::string arq_saida = std::string("out_") + arq; // nome do arquivo de saida
	std::ofstream out(arq_saida.c_str(), std::ios::app);

	float **A; /* matriz */
	float *b; /* vetor */
	int N; /* tamanho do vetor b */
	int linha;
	float erro;
	int iteracoes;
	float *solucao; /* vetor global da solucao */

	if( !LerArquivo(arq, A, b, linha, N, erro, iteracoes) ) {
		printf("Erro ao ler arquivo.\n");
		return 1;
	}
	if (!out.is_open()) {
		printf("Erro ao criar arquivo.\n");
		goto end;
	}

	solucao = new float[N];

	// Pthreads -----------------------------------------------------------------------------------
	{
		int num_threads = 4;
		pthread_t threads[num_threads];
        Dados d[num_threads];

		time_t t1 = time(NULL); // inicio do processamento
		for (int k = 0, i = 0, j = 0; k < num_threads && i < N && j < N; ++k)
		{
			d[k].N = N/2;
			d[k].A = new float*[d[k].N];
			for (int m = 0; m < d[k].N; ++m)
				d[k].A[m] = new float[d[k].N];
            for (int m = 0; m < d[k].N; ++m)
                for(int n = 0; n < d[k].N; ++n)
                    d[k].A[m][n] = A[i+m][j+n];

			if(k % 2 == 0) {
				d[k].b = &b[0];
			}
			else {
				d[k].b = &b[d[k].N];
			}
			d[k].solucao = new float[d[k].N];
			d[k].erro = erro;
			d[k].iteracoes = iteracoes;

			if( pthread_create(&threads[k], NULL, JacobiRichardson, (void*)&d[k]) ) {
				printf("Erro ao criar threads.\n");
				goto end;
			}
			j += N/2;
			if (j >= N) {
				j = 0;
				i += N/2;
			}
		}
		// esperando threads terminar
		for (int i = 0; i < num_threads; ++i)
		{
			pthread_join(threads[i], NULL);
		}
		//----------------------------------------------------------------------------------------------------------
		time_t t2 = time(NULL);
		out << "\nTempo: " << t2 - t1 << "s\n";

        memset(solucao, 0, sizeof(float[N]));
		for(int k = 0; k < num_threads; ++k)
        {
            int i;
            if(k % 2 == 0)
                i = 0;
            else
                i = d[k].N;
            for(int j = 0; i < N && j < d[k].N; ++i, ++j)
                solucao[i] += d[k].solucao[j];
        }

		float r = CalcularLinha(A, N, solucao, linha);

		out << "-------------------------------\nIterations: "<< iteracoes << "\n";
		out << "RowTest: " << linha << " => [" << r << "] =? " << b[linha] << "\n-------------------------------\n";
		out.close();

		// Liberando memoria
		for(int k = 0; k < num_threads; ++k)
        {
            for(int i = 0; i < d[k].N; ++i)
                delete[] d[k].A[i];
            delete[] d[k].A;
        }
	}

	end:
	delete[] b;
	delete[] solucao;
	for(int i = 0; i < N; ++i)
		delete[] A[i];
	delete[] A;
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
void *JacobiRichardson(void *arg)
{
	Dados d = *((Dados*)arg);
	float **A = d.A;
	float *b = d.b;
	int N = d.N;
	int iteracoes = d.iteracoes;
	float erro = d.erro;
	float *solucao = d.solucao;

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
    } while( (numerador/denominador) > erro && contador < iteracoes );

    /* vetor de saida */
    for(i = 0; i < N; ++i) {
        solucao[i] = x[i];
    }

    iteracoes = contador;
	delete[] x;
	delete[] y;
	pthread_exit(NULL);
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
