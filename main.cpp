#include <iostream>
#include <queue>
#include <random>
#include <functional>
#include <math.h>
#include <time.h>
#include "Event.h"

#define PROGRESS_BAR // Mostra progresso para cada rodada de simulação
//#define LOG // Coloca alguns outputs a mais para acompanhar valores
#define EMPLACE // Melhora o uso de memória (Comentado durante desenvolvimento para melhor uso da IDE)

// Pequeno contorno à diferença entre push e emplace
#ifdef EMPLACE
#define PUSH(type) emplace
#define ENDPUSH
#else
#define PUSH(type) push(type
#define ENDPUSH )
#endif

using namespace std;

/*PARÂMETROS*/
// Parâmetros da simulação
int SAMPLES = 100;                 // Número de amostras por rodada
int SIMULATIONS = 5;                   // Número de rodadas de simulação
double percentil90(int s){ // Valores retirados da tabela t-Student
	switch  (s){
		case 1:  return 6.314;
		case 2:  return 2.92;
		case 3:  return 2.353;
		case 4:  return 2.132;
		case 5:  return 2.015;
		case 6:  return 1.943;
		case 7:  return 1.895;
		case 8:  return 1.860;
		case 9:  return 1.833;
		case 10: return 1.812;
		case 11: return 1.796;
		case 12: return 1.782;
		case 13: return 1.771;
		case 14: return 1.761;
		case 15: return 1.753;
		case 16: return 1.746;
		case 17: return 1.740;
		case 18: return 1.734;
		case 19: return 1.729;
		case 20: return 1.725;
		case 21: return 1.721;
		case 22: return 1.717;
		case 23: return 1.714;
		case 24: return 1.711;
		case 25: return 1.708;
		case 26: return 1.706;
		case 27: return 1.703;
		case 28: return 1.701;
		case 29: return 1.1699;
		default: return 1.645;
	}
};
bool PREEMPTION = false;               // Interrupção dos pacotes de dados em caso de chegada de pacote de voz
double UTILIZATION_1 = 0.1;            // ρ1 (utilização da fila de dados)
constexpr double SERVER_SPEED = 2e6;   // Velocidade do servidor (2Mb/segundo)
int TRANSIENT_SAMPLE_NUMBER = 50000;  // Amostras colocadas no período transiente
enum class PrintMode { NORMAL, JSON, TRANSIENT, CONFIDENCE_INTERVAL };
PrintMode printMode = PrintMode::NORMAL;
struct SimulationRound {
	int n1Packages, n2Packages, n2Intervals;                    // Controle do número de pacotes da rodada
	double startTime, totalTime,                                // Controle da duração da rodada
			totalDataTimeInQueue, totalVoiceTimeInQueue,        // Tempo*Fregueses na fila durante a rodada
			T1Acc, X1Acc, T2Acc, jitterAcc, jitterAccSqr,        // Acumulado das estatísticas dos fregueses que chegaram durante a rodada
			T1, X1, Nq1, T2, Nq2, JitterMean, JitterVariance;    // Estatísticas finais da rodada
};

// Variáveis globais de Debug
int FIXED_ACTIVE_PERIOD = 0; // Fixar qt. pacotes no periodo ativo
double FIXED_SILENCE_PERIOD = -1; // Fixar periodo de silencio
double debug_tempo_entre_chegadas = 0;
double debug_max_time = 0;
int *debug_activePeriodLength;

double genRandUnitary() {
	return ((double) rand() / RAND_MAX);
}

/*Data*/
// Data channel random variable generators
int genDataPackageSize() {
	double x = genRandUnitary();
	if (x < .3)
		return 64;
	else if (x < .4)
		return 512;
	else if (x < .7)
		return 1500;
	else
		// TODO corrigir probabilidade das pontas
		return static_cast<int>(genRandUnitary() * 1436 + 64);
}

auto genDataServiceTime = []() {
	return genDataPackageSize() * 8 / SERVER_SPEED;
};
#define DATA_TIME_OF_SERVICE genDataServiceTime()

double DATA_ARRIVAL_RATE = UTILIZATION_1 / (755 * 8 / SERVER_SPEED); // λ1 = ρ1/E[X1] = ρ1/(E[L]bytes*8/(2Mb/s))
double genDataArrivalTime() {
	return -log(genRandUnitary()) / (DATA_ARRIVAL_RATE);
}

#define DATA_ARRIVAL_TIME genDataArrivalTime()

/*Voice*/
int VOICE_CHANNELS = 30;
double VOICE_ARRIVAL_TIME = .016; // Tempo até o próximo pacote de voz durante o período ativo em segundos
int VOICE_PACKAGE_SIZE_IN_BITS = 512;
double VOICE_TIME_OF_SERVICE = VOICE_PACKAGE_SIZE_IN_BITS / SERVER_SPEED; // Tempo de transmissão do pacote de voz em segundos
const int MEAN_N_VOICE_PACKAGE = 22;
const double MEAN_SILENCE_PERIOD_DURATION = .65; // In seconds
double X2 = VOICE_TIME_OF_SERVICE;

// Voice channel random variable generators
auto genEndOfActivePeriod = []() {
	return genRandUnitary() < (1.0 / MEAN_N_VOICE_PACKAGE);
};
auto genSilencePeriod = []() {
    if (FIXED_SILENCE_PERIOD != -1) return FIXED_SILENCE_PERIOD;
	return -log(genRandUnitary()) / (1.0 / MEAN_SILENCE_PERIOD_DURATION);
};
#define VOICE_SILENCE_TIME genSilencePeriod()

// Variáveis globais
Packet *server;
double *channelsLastDeparture; // controle de última saída dos canais de voz (usado no cálculo do jitter)


/*FUNÇÕES*/
// Configurações iniciais

/**
 * Calcula o tempo da primeira chegada de cada evento e coloca na heap de eventos
 * @param arrivals A heap de eventos
 */
void setup(priority_queue<Event> &arrivals) {
	arrivals.PUSH(Event)(DATA_ARRIVAL_TIME, EventType::DATA, new Packet(-1, DATA_TIME_OF_SERVICE))ENDPUSH;
    channelsLastDeparture = new double[VOICE_CHANNELS];
    debug_activePeriodLength = new int[VOICE_CHANNELS];

	for (int i = 0; i < VOICE_CHANNELS; ++i) {
		arrivals.PUSH(Event)(VOICE_SILENCE_TIME, EventType::VOICE, new Packet(-1, i, VOICE_TIME_OF_SERVICE))ENDPUSH;
		channelsLastDeparture[i] = -1;
		debug_activePeriodLength[i] = -1;
	}
}

/**
 * Zera todos os parâmetros da rodada de simulação
 * @param s
 */
void setupSimulationStats(SimulationRound &s) {
	s.X1Acc = s.T1Acc = s.T2Acc = s.jitterAcc = s.jitterAccSqr =
	s.T1 = s.X1 = s.Nq1 = s.T2 = s.Nq2 = s.JitterMean = s.JitterVariance =
	s.startTime = s.totalTime = s.totalDataTimeInQueue = s.totalVoiceTimeInQueue =
	s.n1Packages = s.n2Packages = s.n2Intervals = 0;
}

// Colocar eventos na heap
/**
 * Coloca o evento na fila de eventos esperados com o tempo adequado e dá o valor para a estatísca de tempo de espera
 * @param event Evento a ser servido (tipo VOICE ou DATA)
 * @param events_queue A heap de eventos
 * @param currentTime Momento em que o evento foi servido
 */
Packet *serveEvent(Packet *packet, priority_queue<Event> &events_queue, double currentTime) {
#ifdef LOG
	cout << "> Gerado novo serviço de dados terminando em: " << currentTime + packet->serviceTime << endl;
#endif
	events_queue.PUSH(Event)(currentTime + packet->serviceTime, EventType::SERVER, packet)ENDPUSH;
	return packet;
}

/**
 * Serve o próximo evento em espera
 * @param events_queue Heap de eventos
 * @param voice Fila de eventos de voz (maior prioridade)
 * @param data Fila de eventos de dados (menor prioridade)
 * @param currentTime Momento em que o evento será servido
 * @return O tipo de evento colocado no servidor (EMPTY caso nenhum evento tenha sido colocado)
 */
Packet *serveNextEvent(priority_queue<Event> &events_queue, queue<Packet *> &voice, queue<Packet *> &data, double currentTime) {
	if (!voice.empty()) {
		Packet *p = voice.front();
		voice.pop();
		return serveEvent(p, events_queue, currentTime);
	} else if (!data.empty()) {
		return serveEvent(data.front(), events_queue, currentTime);
	}
	return nullptr;
}

// Cálculo de estatísticas
/**
 * Incrementa o acumulado de pessoas na fila ao longo da rodada s
 * @param Nq2
 * @param Nq1
 * @param lastTime
 * @param currentTime
 * @param s
 */
void registerAreaStatistics(unsigned long Nq2, unsigned long Nq1, double lastTime, double currentTime, SimulationRound &s) {
	double t = currentTime - lastTime;
	if (server != nullptr && server->type == PacketType::DATA) {
		Nq1--;
	}
	s.totalDataTimeInQueue += Nq1 * t;
	s.totalVoiceTimeInQueue += Nq2 * t;
}

/**
 * Insere as estatísticas do pacote na rodada de simulação equivalente e deleta o objeto para evitar vazamento de memória
 * @param packet
 * @param rounds
 */
void countPacketIntoStatistics(Packet *packet, SimulationRound &s) {
	if (packet->simulation != -1) {
		switch (packet->type) {
			case PacketType::DATA:
				s.T1Acc += packet->totalTime;
				s.n1Packages++;
				s.X1Acc += packet->serviceTime + packet->property.wastedTime;
				break;
			case PacketType::VOICE:
				s.T2Acc += packet->totalTime;
				s.n2Packages++;
				// Jitter incrementado do lado de fora
				break;
		}
	}
	delete packet;
}

/**
 * Inclui o intervalo do canal no acumulado da rodada s
 * @param s
 * @param channel
 * @param currentTime
 */
void incrementJitter(SimulationRound &s, int channel, double currentTime) {
	double lastDeparture = channelsLastDeparture[channel];
	if (lastDeparture != -1) {
		double lastInterval = currentTime - lastDeparture;
		s.jitterAcc += lastInterval;
		s.jitterAccSqr += lastInterval * lastInterval;
		s.n2Intervals += 1;
	}
}

/**
 * Transforma os acumulados da rodada s em estatísticas
 * @param s
 */
void calculateRoundStatistics(SimulationRound &s) {
	s.T1 = s.n1Packages == 0 ? 0 : s.T1Acc / s.n1Packages;
	s.Nq1 = s.totalDataTimeInQueue == 0 ? 0 : s.totalDataTimeInQueue / s.totalTime;
	s.X1 = s.n1Packages == 0 ? 0 : s.X1Acc / s.n1Packages;

	s.T2 = s.n2Packages == 0 ? 0 :s.T2Acc / s.n2Packages;
	s.Nq2 = s.totalVoiceTimeInQueue == 0 ? 0 : (s.totalVoiceTimeInQueue / (s.totalTime * VOICE_CHANNELS));
	s.JitterMean = s.jitterAcc == 0 ? 0 : s.jitterAcc / s.n2Intervals;
	s.JitterVariance = (s.jitterAcc == 0 || s.jitterAccSqr == 0 || s.n2Intervals < 2) ? 0 :
					   (s.jitterAccSqr - s.jitterAcc * s.jitterAcc / s.n2Intervals) / (s.n2Intervals - 1);
}

// Outputs
void printStats(SimulationRound &s) {
#ifdef LOG
	cout << "Tempo médio entre chegadas: " << (debug_tempo_entre_chegadas / SAMPLES) << endl;
	cout << "lambda_1: " << n1Packages / debug_max_time << " pacotes dados/seg" << endl;
	cout << "lambda_2: " << n2Packages / debug_max_time << " pacotes voz/seg" << endl;
	cout << "Max Time: " << debug_max_time << endl;
#endif
	if (printMode >= PrintMode::JSON) {
		cout << R"({"t1":)" << s.T1 << R"(,"w1":)" << s.T1 - s.X1 << R"(,"x1":)" << s.X1 << R"(,"nq1":)" << s.Nq1 <<
			 R"(,"t2":)" << s.T2 << R"(,"w2":)" << s.T2 - X2 << R"(,"x2":)" << X2 << R"(,"nq2":)" << s.Nq2 <<
			 R"(,"jitterE":)" << s.JitterMean << R"(,"jitterV":)" << s.JitterVariance << R"(})";
	} else {
		cout << "E[T1]: " << s.T1 << ", E[W1]: " << s.T1 - s.X1 << ", E[X1]: " << s.X1 << ", E[Nq1]: " << s.Nq1 << endl;
		cout << "E[T2]: " << s.T2 << ", E[W2]: " << s.T2 - X2 <<
			 ", E[X2]: " << X2 << ", E[Nq2]: " << s.Nq2 << ", E[Δ]: " << s.JitterMean << ", V(Δ):" << s.JitterVariance << endl;
	}
}

void printConfidenceInterval(SimulationRound rounds[], int simulations = SIMULATIONS) {
	if (simulations < 2)
		return;
	SimulationRound &f = rounds[0];
	double  T1  = f.T1 / simulations,
			X1  = f.X1 / simulations,
			Nq1  = f.Nq1 / simulations,
			T2  = f.T2 / simulations,
			Nq2  = f.Nq2 / simulations,
			JitterMean  = f.JitterMean / simulations,
			JitterVariance  = f.JitterVariance / simulations;
	double  W1 = T1 - X1,
			W2 = T2 - X2,
			VarianceT1 = 0,
			VarianceW1 = 0,
			VarianceX1 = 0,
			VarianceNq1 = 0,
			VarianceT2 = 0,
			VarianceW2 = 0,
			VarianceNq2 = 0,
			VarianceJitterMean = 0,
			VarianceJitterVariance = 0;
	if (printMode != PrintMode::CONFIDENCE_INTERVAL){
		if (printMode >= PrintMode::JSON) {
			cout << R"(,"average":)";
		} else {
			cout << "=================================================" << endl;
			cout << "=================================================" << endl;
			cout << "=================================================" << endl;
			cout << "Resultado médio final" << endl;
		}
		SimulationRound temp{};
		temp.T1 = T1;
		temp.X1 = X1;
		temp.Nq1 = Nq1;
		temp.T2 = T2;
		temp.Nq2 = Nq2;
		temp.JitterMean = JitterMean;
		temp.JitterVariance = JitterVariance;
		printStats(temp);
	}
	for (int i = 1; i <= simulations; ++i) {
		SimulationRound &s = rounds[i];
		VarianceT1 += (s.T1 - T1) * (s.T1 - T1);
		VarianceW1 += (s.T1 - s.X1 - W1) * (s.T1 - s.X1 - W1);
		VarianceX1 += (s.X1 - X1) * (s.X1 - X1);
		VarianceNq1 += (s.Nq1 - Nq1) * (s.Nq1 - Nq1);
		VarianceT2 += (s.T2 - T2) * (s.T2 - T2);
		VarianceW2 += (s.T2 - X2 - W2) * (s.T2 - X2 - W2);
		VarianceNq2 += (s.Nq2 - Nq2) * (s.Nq2 - Nq2);
		VarianceJitterMean += (s.JitterMean - JitterMean) * (s.JitterMean - JitterMean);
		VarianceJitterVariance += (s.JitterVariance - JitterVariance) * (s.JitterVariance - JitterVariance);
	}
	VarianceT1 /= (simulations - 1);
	VarianceW1 /= (simulations - 1);
	VarianceX1 /= (simulations - 1);
	VarianceNq1 /= (simulations - 1);
	VarianceT2 /= (simulations - 1);
	VarianceW2 /= (simulations - 1);
	VarianceNq2 /= (simulations - 1);
	VarianceJitterMean /= (simulations - 1);
	VarianceJitterVariance /= (simulations - 1);

	VarianceT1 = sqrt(VarianceT1);
	VarianceW1 = sqrt(VarianceW1);
	VarianceX1 = sqrt(VarianceX1);
	VarianceNq1 = sqrt(VarianceNq1);
	VarianceT2 = sqrt(VarianceT2);
	VarianceW2 = sqrt(VarianceW2);
	VarianceNq2 = sqrt(VarianceNq2);
	VarianceJitterMean = sqrt(VarianceJitterMean);
	VarianceJitterVariance = sqrt(VarianceJitterVariance);

	double p = percentil90(simulations) / sqrt(simulations);
	VarianceT1 *= p;
	VarianceW1 *= p;
	VarianceX1 *= p;
	VarianceNq1 *= p;
	VarianceT2 *= p;
	VarianceW2 *= p;
	VarianceNq2 *= p;
	VarianceJitterMean *= p;
	VarianceJitterVariance *= p;
	if (printMode >= PrintMode::JSON) {
		cout << ",\"confidenceInterval\":{" <<
			 R"("t1":{"lesser":)" << T1 - VarianceT1 << ",\"mean\":" << T1 << ",\"upper\":" << T1 + VarianceT1 << ",\"percentage\":" << (T1 ==0?0:(VarianceT1 * 200 / T1))<< "}," <<
			 R"("w1":{"lesser":)" << W1 - VarianceW1 << ",\"mean\":" << W1 << ",\"upper\":" << W1 + VarianceW1 << ",\"percentage\":" << (W1 ==0?0:(VarianceW1 * 200 / W1)) << "}," <<
			 R"("x1":{"lesser":)" << X1 - VarianceX1 << ",\"mean\":" << X1 << ",\"upper\":" << X1 + VarianceX1 << ",\"percentage\":" << (X1 ==0?0:(VarianceX1 * 200 / X1)) << "}," <<
			 R"("nq1":{"lesser":)" << Nq1 - VarianceNq1 << ",\"mean\":" << Nq1 << ",\"upper\":" << Nq1 + VarianceNq1 << ",\"percentage\":" << (Nq1 ==0?0:(VarianceNq1 * 200 / Nq1)) << "}," <<
			 R"("t2":{"lesser":)" << T2 - VarianceT2 << ",\"mean\":" << T2 << ",\"upper\":" << T2 + VarianceT2 << ",\"percentage\":" << (T2 ==0?0:(VarianceT2 * 200 / T2)) << "}," <<
			 R"("w2":{"lesser":)" << W2 - VarianceW2 << ",\"mean\":" << W2 << ",\"upper\":" << W2 + VarianceW2 << ",\"percentage\":" << (W2 ==0?0:(VarianceW2 * 200 / W2)) << "}," <<
			 R"("nq2":{"lesser":)" << Nq2 - VarianceNq2 << ",\"mean\":" << Nq2 << ",\"upper\":" << Nq2 + VarianceNq2 << ",\"percentage\":" << (Nq2 ==0?0:(VarianceNq2 * 200 / Nq2 ))<< "}," <<
			 R"("jitterMean":{"lesser":)" << JitterMean - VarianceJitterMean << ",\"mean\":" << JitterMean << ",\"upper\":" << JitterMean + VarianceJitterMean <<
			 ",\"percentage\":" << (JitterMean ==0?0:(VarianceJitterMean * 200 / JitterMean)) << "}," <<
			 R"("jitterVariance":{"lesser":)" << JitterVariance - VarianceJitterVariance << ",\"mean\":" << JitterVariance << ",\"upper\":" << JitterVariance + VarianceJitterVariance <<
			 ",\"percentage\":" << (JitterVariance ==0?0:(VarianceJitterVariance * 200 / JitterVariance)) << "}}";
	} else {
		cout << "T1:\t\t\t\t" << T1 - VarianceT1 << "\t" << T1 << "\t" << T1 + VarianceT1 << "\t" << VarianceT1 * 200 / T1 << "%" << endl;
		cout << "W1:\t\t\t\t" << W1 - VarianceW1 << "\t" << W1 << "\t" << W1 + VarianceW1 << "\t" << VarianceW1 * 200 / W1 << "%" << endl;
		cout << "X1:\t\t\t\t" << X1 - VarianceX1 << "\t" << X1 << "\t" << X1 + VarianceX1 << "\t" << VarianceX1 * 200 / X1 << "%" << endl;
		cout << "Nq1:\t\t\t" << Nq1 - VarianceNq1 << "\t" << Nq1 << "\t" << Nq1 + VarianceNq1 << "\t" << VarianceNq1 * 200 / Nq1 << "%" << endl;
		cout << "T2:\t\t\t\t" << T2 - VarianceT2 << "\t" << T2 << "\t" << T2 + VarianceT2 << "\t" << VarianceT2 * 200 / T2 << "%" << endl;
		cout << "W2:\t\t\t\t" << W2 - VarianceW2 << "\t" << W2 << "\t" << W2 + VarianceW2 << "\t" << VarianceW2 * 200 / W2 << "%" << endl;
		cout << "Nq2:\t\t\t" << Nq2 - VarianceNq2 << "\t" << Nq2 << "\t" << Nq2 + VarianceNq2 << "\t" << VarianceNq2 * 200 / Nq2 << "%" << endl;
		cout << "JitterMean:\t\t" << JitterMean - VarianceJitterMean << "\t" << JitterMean << "\t" << JitterMean + VarianceJitterMean << "\t" << (JitterMean ==0?0:(VarianceJitterMean * 200 / JitterMean)) << "%" << endl;
		cout << "JitterVariance:\t" << JitterVariance - VarianceJitterVariance << "\t" << JitterVariance << "\t" << JitterVariance + VarianceJitterVariance << "\t" << (JitterVariance ==0?0:(VarianceJitterVariance * 200 / JitterVariance))<< "%" << endl;
	}
}

/**
 * Imprime em stdout o valor das estatísticas globais: E[T1], E[W1], E[X1], E[Nq1], E[T2], E[W2], E[X2], E[Nq2], E[Δ] e V(Δ)
 */
void printStats(SimulationRound rounds[], int r) {
#ifdef LOG
	cout << "Tempo médio entre chegadas: " << (debug_tempo_entre_chegadas / SAMPLES) << endl;
	cout << "lambda_1: " << n1Packages / debug_max_time << " pacotes dados/seg" << endl;
	cout << "lambda_2: " << n2Packages / debug_max_time << " pacotes voz/seg" << endl;
	cout << "Max Time: " << debug_max_time << endl;
#endif
	SimulationRound &s = rounds[r];
	if (printMode >= PrintMode::JSON) {
		cout << R"({"t1":)" << s.T1 << R"(,"w1":)" << s.T1 - s.X1 << R"(,"x1":)" << s.X1 << R"(,"nq1":)" << s.Nq1 <<
			 R"(,"t2":)" << s.T2 << R"(,"w2":)" << s.T2 - X2 << R"(,"x2":)" << X2 << R"(,"nq2":)" << s.Nq2 <<
			 R"(,"jitterE":)" << s.JitterMean << R"(,"jitterV":)" << s.JitterVariance;
		if (printMode == PrintMode::CONFIDENCE_INTERVAL) printConfidenceInterval(rounds, r);
		cout << R"(})";
	} else {
		cout << "E[T1]: " << s.T1 << ", E[W1]: " << s.T1 - s.X1 << ", E[X1]: " << s.X1 << ", E[Nq1]: " << s.Nq1 << endl;
		cout << "E[T2]: " << s.T2 << ", E[W2]: " << s.T2 - X2 <<
			 ", E[X2]: " << X2 << ", E[Nq2]: " << s.Nq2 << ", E[Δ]: " << s.JitterMean << ", V(Δ):" << s.JitterVariance << endl;
	}
}

/**
-  * Imprime texto de ajuda para mostrar as opções a serem passadas para o programa
-  */
void printHelp() {
	cout << "Utilização:" << endl;
	cout << "_exe_ [OPTIONS]" << endl;
	cout << "\t-h\t\tMostra essa ajuda" << endl;
	cout << "\t-a int\tNúmero de amostras" << endl;
	cout << "\t-t int\tNúmero de amostras do período transiente" << endl;
	cout << "\t-r int\tNúmero de rodadas de simulação" << endl;
	cout << "\t-v int\tNúmero de canais de voz" << endl;
	cout << "\t-s int\tNúmero pacotes por canal de voz" << endl;
	cout << "\t-u double\tUtilização da fila 1" << endl;
	cout << "\t-i\tFila de dados pode ser interrompida" << endl;
	cout << "\t-pj\tImprime no formato JSON" << endl;
	cout << "\t-pt\tImprime no formato JSON as médias sendo cálculadas a durante o período transiente" << endl;
	cout << "\t-pi\tImprime no formato JSON com os valores atualizados de I.C. a cada rodada da simulação" << endl;
}

void runSimulationRound(
		SimulationRound rounds[], int s, int samples,
		priority_queue<Event> &arrivals, queue<Packet *> &voice, queue<Packet *> &data,
		double &lastTime, int &interruptedDataPackages) {
	double t;
	for (int i = 0; i < samples; ++i) {
	    // Barra de progresso da simulação
		#ifdef PROGRESS_BAR
		if (printMode == PrintMode::NORMAL && i % (samples / 20) == 0) cout << "Rodando simulação " << s << ": " << (i * 100 / samples) << "%" << '\r' << flush;
		#endif
		Event arrival = arrivals.top();
		arrivals.pop();

		// Atualiza estatísticas
		registerAreaStatistics(voice.size(), data.size(), lastTime, arrival.time, rounds[s]);

		#ifdef LOG
		cout << "!! Time: " << arrival.time << (arrival.type == EventType::VOICE ? " (Voice)" : " (Server)") << endl;
			cout << "channelsLastDeparture[0]: " << channelsLastDeparture[0] << endl;
			cout << "Instante " << arrival.time << "s" << endl;
			cout << "Filas:" << endl << "> Dados: " << (data.size() == 0 ? 0 : data.size() - 1) << endl << "> Voz: " << (voice.size() == 0 ? 0 : voice.size() - 1) << endl << endl;
		#endif

		switch (arrival.type) {
		    // Chegada de um pacote de dados
			case EventType::DATA:
				// 'Total Time' se inicia com o negativo do instante em que o pacote entra no sistema
                // e posteriormente é somado ao instante em que ele sai, gerando o tempo total que ele passou no sistema
				arrival.packet->totalTime = -arrival.time;

				// Coloca pacote na fila
				data.push(arrival.packet);

				// Gera a próxima chegada de pacote de dados na fila de eventos
				arrivals.PUSH(Event)(arrival.time + DATA_ARRIVAL_TIME, EventType::DATA, new Packet(s, DATA_TIME_OF_SERVICE))ENDPUSH;

				if (server == nullptr) {
				    // Se servidor vazio, gera o seu evento de serviço
					server = serveEvent(arrival.packet, arrivals, arrival.time);
				}
				break;
            // Chegada de um pacote de voz
            case EventType::VOICE:
				// Coloca a próxima chegada do canal na heap de eventos
				arrival.packet->totalTime = -arrival.time;
				t = arrival.time + VOICE_ARRIVAL_TIME;

				// Fixar o período ativo (útil para testes de correção)
				if (FIXED_ACTIVE_PERIOD) {
                    if (debug_activePeriodLength[0] < FIXED_ACTIVE_PERIOD) {
                        debug_activePeriodLength[0]++;
                        arrival.packet->property.channel.lastVoicePacket = false;
                    } else {
                        t += VOICE_SILENCE_TIME;
                        arrival.packet->property.channel.lastVoicePacket = true;
                        debug_activePeriodLength[0] = 0;
                    }
				}
				else {
				// Se o período ativo deveria terminar, marca o pacote como o último do período ativo e
                // insere um periodo de silencio exponencial antes da chegada do próximo pacote
                    if (genEndOfActivePeriod()) {
                        t += VOICE_SILENCE_TIME;
                        arrival.packet->property.channel.lastVoicePacket = true;
                    }
                }

                // Insere próximo pacote de voz
				arrivals.PUSH(Event)(t, EventType::VOICE, new Packet(s, arrival.packet->property.channel.number, VOICE_TIME_OF_SERVICE))ENDPUSH;

				if (server == nullptr) {
				    // Se servidor vazio, entra no servidor e gera o evento de sua saída do servidor
					serveEvent(arrival.packet, arrivals, arrival.time);
					server = arrival.packet;
				} else if (PREEMPTION && server->type == PacketType::DATA) {
                    // Se servidor ocupado com dados e preempção ativada, entra no servidor e calcula o tempo desperdiçado no servidor do pacote de dados
                    // O evento de término de serviço do pacote de dados ainda ocorrerá, porém ele verá que 'interruptedDataPackages'
                    // é maior do que 0, e portanto não irá 'efetivar' o seu envio
                    interruptedDataPackages++;
					data.front()->property.wastedTime -= arrival.time;
					// Gera-se o seu evento de serviço
					server = serveEvent(arrival.packet, arrivals, arrival.time);
				} else {
				    // Servidor ocupado com voz, portanto insere o pacote de voz na fila de espera
					voice.push(arrival.packet);
				}
				break;
            // Término do serviço de um pacote (dados ou voz)
            case EventType::SERVER:
				i--; //Saída do simulador não conta como amostra para a contagem
				SimulationRound &r = rounds[arrival.packet->simulation];
                switch (arrival.packet->type) {
                    // Pacote de dados é transmitido
                    case PacketType::DATA:
                        if (interruptedDataPackages) {
                            // O pacote foi interrompido, portanto na verdade não terminou seu serviço
                            // Ele continua na sua posição na fila e decrementa-se o número de pacotes de dados interrompidos
                            // Seu tempo desperdiçado de serviço é atualizado
                            interruptedDataPackages--;
                            data.front()->property.wastedTime += arrival.time;
                        } else {
                            // Tira pacote de dados da fila, atualiza seu tempo total passado no sistema, atualiza estatísticas
                            data.pop();
                            arrival.packet->totalTime += arrival.time;
                            countPacketIntoStatistics(arrival.packet, r);
                            // Gera-se próximo serviço
                            server = serveNextEvent(arrivals, voice, data, arrival.time);
                        }
                        break;
                    // Pacote de voz é transmitido
                    case PacketType::VOICE:
					    // Calcula estatísticas de jitter
						incrementJitter(r, arrival.packet->property.channel.number, arrival.time);
						// Atualiza controle do tempo de saída do último pacote de voz do canal
						if (!arrival.packet->property.channel.lastVoicePacket) {
							channelsLastDeparture[arrival.packet->property.channel.number] = arrival.time;
						} else { ;
							channelsLastDeparture[arrival.packet->property.channel.number] = -1;
						}
						// Atualiza tempo total passado no sistema, atualiza estatisticas
						arrival.packet->totalTime += arrival.time;
						countPacketIntoStatistics(arrival.packet, r);
						// Gera-se o próximo evento de serviço
						server = serveNextEvent(arrivals, voice, data, arrival.time);
						break;
				}
				break;
		}

		// Atualiza o último instante da simulação
		lastTime = arrival.time;
		if (printMode == PrintMode::TRANSIENT){
			static int c = i, m = samples/200;
			if (c!=i && !(i%m)){
				c = i;
				SimulationRound &a = rounds[0];
				a.totalTime = arrival.time;
				calculateRoundStatistics(a);
				static bool first = true;
				if (first) first = false;
				else cout << "," ;
				printStats (a);
			}
		}
	}
}

void printParameters() {
	cout << "{\"samples\":" << SAMPLES;
}

int main(int argc, char *argv[]) {
	// Pegando opções na linha de comando
	for (int p = 1; p < argc; ++p) {
		string option = argv[p];
		if (option[0] != '-') {
			cout << "Opção " << argv[p] << " inválida" << endl;
			continue;
		}
		switch (option[1]) {
			case 'a': //Número de amostras
				SAMPLES = stoi(argv[++p]);
				break;
			case 'i': // Com ou sem preempção
				PREEMPTION = true;
				break;
			case 'r': // Número de rodadas de simulação
				SIMULATIONS = stoi(argv[++p]);
				break;
			case 't': // Número de amostras do período transiente
				TRANSIENT_SAMPLE_NUMBER = stoi(argv[++p]);
				break;
			case 'u': // Utilização da fila de dados
				UTILIZATION_1 = stod(argv[++p]);
				DATA_ARRIVAL_RATE = UTILIZATION_1 / (755 * 8 / SERVER_SPEED); // Î»1 = Ï1/E[X1] = Ï1/(E[L]bytes*8/(2Mb/s))
				break;
			case 'p':
				switch (option[2]){
					case 'i': // Imprimir os resultados como JSON para medir a largura do intervalo de confiança
						printMode = PrintMode::CONFIDENCE_INTERVAL;
						break;
					case 'j': // Imprimir os resultados como JSON
						printMode = PrintMode::JSON;
						break;
					case 't': // Imprimir os resultados como JSON para medir o período transiente
						printMode = PrintMode::TRANSIENT;
						break;
					default: // Imprimir os resultados de forma amigável para o usuário
						printMode = PrintMode::NORMAL;
						break;
				}
				break;
			case 'h':
				printHelp();
				return 0;
            case 'q': // Número de rodadas de simulação
                SIMULATIONS = stoi(argv[++p]);
                break;
            case 'v': // Numero de canais de voz
                VOICE_CHANNELS = stoi(argv[++p]);
                break;
            case 's': // Numero fixo de pacotes de voz por periodo ativo
                FIXED_ACTIVE_PERIOD = stoi(argv[++p]);
                break;
		    case 'w': // Mudar 'gap' entre pacotes de voz (original: 16ms)
                VOICE_ARRIVAL_TIME = stod(argv[++p]);
                break;
            case 'y': // Duraçao fixa do periodo de silencio
                FIXED_SILENCE_PERIOD = stod(argv[++p]);
                break;
		    case 'z': // Outro tamanho para o pacote de voz
                VOICE_PACKAGE_SIZE_IN_BITS = stoi(argv[++p]);
                VOICE_TIME_OF_SERVICE = VOICE_PACKAGE_SIZE_IN_BITS / SERVER_SPEED;
                X2 = VOICE_TIME_OF_SERVICE;
                break;
			default:
				cout << "Opção " << argv[p] << " inválida" << endl;
				continue;
		}
	}

	// Filas e variáveis de controle
	priority_queue<Event> arrivals; // Estrutura para organizar as chegadas
	queue<Packet *> data, voice; // Filas de data e voz
	int interruptedDataPackages = 0; // Pacotes de dados interrompidos (usado para invalidar os respectivos pacotes quando a antiga saída aparecer na lista)
	SimulationRound rounds[SIMULATIONS + 1]; // Estatísticas das várias rodadas de simulação.
	double lastTime = 0;

	// Configurações iniciais
	srand(static_cast<unsigned int>(time(nullptr)));
	setup(arrivals);
	for (int i = 0; i <= SIMULATIONS; ++i) {
		setupSimulationStats(rounds[i]);
	}

	// Período transiente
	if (printMode >= PrintMode::JSON){
		printParameters();
		cout << R"(,"simulations":[)";
	}
	runSimulationRound(rounds, 0, TRANSIENT_SAMPLE_NUMBER, arrivals, voice, data, lastTime, interruptedDataPackages);

	// Executando rodadas de simulação
	for (int s = 1; s <= SIMULATIONS; ++s) {
		rounds[s].startTime = lastTime;
		runSimulationRound(rounds, s, SAMPLES, arrivals, voice, data, lastTime, interruptedDataPackages);
		rounds[s].totalTime = lastTime - rounds[s].startTime;
	}

	// Rodadas extras para finalizar os fregueses da última rodada (necessário?)
	for (int s = 1; s <= SIMULATIONS; ++s) {
		runSimulationRound(rounds, 0, SAMPLES, arrivals, voice, data, lastTime, interruptedDataPackages);
	}

	if (printMode != PrintMode::TRANSIENT){
		// Cálculo e output das estatísticas de cada rodada
		SimulationRound &f = rounds[0];
		setupSimulationStats(f);
		for (int i = 1; i <= SIMULATIONS; ++i) {
			SimulationRound &s = rounds[i];
			calculateRoundStatistics(s);
			f.T1 += s.T1;
			f.X1 += s.X1;
			f.Nq1 += s.Nq1;
			f.T2 += s.T2;
			f.Nq2 += s.Nq2;
			f.JitterMean += s.JitterMean;
			f.JitterVariance += s.JitterVariance;
			if (printMode >= PrintMode::JSON) {
				if (i > 1) cout << ",";
			} else {
				cout << "Simulação " << i << endl;
			}
			printStats(rounds, i);
		}
		if (printMode >= PrintMode::JSON) cout << "]";

		// Cálculo e output do intervalo de confiança
		printConfidenceInterval(rounds);
		if (printMode >= PrintMode::JSON ) cout << "}" << endl;
	} else {
		cout << "]}" << endl;
	}
	return 0;
}
