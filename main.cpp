#include <iostream>
#include <queue>
#include <random>
#include <functional>
#include <math.h>
#include <time.h>
#include "Event.h"

#define PROGRESS_BAR // Mostra progresso para cada rodada de simulação
//#define LOG // Coloca alguns outputs a mais para acompanhar valores
//#define EMPLACE // Melhora o uso de memória (Comentado durante desenvolvimento para melhor uso da IDE)

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
int SAMPLES = 10000;                 // Número de amostras por rodada
int SIMULATIONS = 50;                   // Número de rodadas de simulação
const double percentil90 = 1.65;
bool PREEMPTION = false;               // Interrupção dos pacotes de dados em caso de chegada de pacote de voz
double UTILIZATION_1 = 0.1;            // ρ1 (utilização da fila de dados)
constexpr double SERVER_SPEED = 2e6;   // Velocidade do servidor (2Mb/segundo)
int TRANSIENT_SAMPLE_NUMBER = 10000;  // Amostras colocadas no período transiente

struct SimulationRound {
	int n1Packages, n2Packages, n2Intervals;                    // Controle do número de pacotes da rodada
	double startTime, totalTime,                                // Controle da duração da rodada
			totalDataTimeInQueue, totalVoiceTimeInQueue,        // Tempo*Fregueses na fila durante a rodada
			T1Acc, X1Acc, T2Acc, jitterAcc, jitterAccSqr,        // Acumulado das estatísticas dos fregueses que chegaram durante a rodada
			T1, X1, Nq1, T2, Nq2, JitterMean, JitterVariance;    // Estatísticas finais da rodada
};

double genRandUnitary() {
	return ((double) rand() / RAND_MAX);
}

/*Data*/
// Data channel random variable generators
int genDataPackageSize() {
//    return 755;
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
const int VOICE_CHANNELS = 30;
const double VOICE_ARRIVAL_TIME = .016; // Tempo até o próximo pacote de voz durante o período ativo em segundos
const int VOICE_PACKAGE_SIZE_IN_BITS = 512;
constexpr double VOICE_TIME_OF_SERVICE = VOICE_PACKAGE_SIZE_IN_BITS / SERVER_SPEED; // Tempo de transmissão do pacote de voz em segundos
const int MEAN_N_VOICE_PACKAGE = 22;
const double MEAN_SILENCE_PERIOD_DURATION = .65; // In seconds
const double X2 = VOICE_TIME_OF_SERVICE;

// Voice channel random variable generators
auto genEndOfActivePeriod = []() {
	return genRandUnitary() < (1.0 / MEAN_N_VOICE_PACKAGE);
};
auto genSilencePeriod = []() {
	return -log(genRandUnitary()) / (1.0 / MEAN_SILENCE_PERIOD_DURATION);
};
#define VOICE_SILENCE_TIME genSilencePeriod()

// Variáveis globais
Packet *server;
double aux_tempo_entre_chegadas = 0;
double max_time = 0;
double channelsLastDeparture[VOICE_CHANNELS];
int activePeriodLength[VOICE_CHANNELS];
bool JSON = false;

/*FUNÇÕES*/
// Configurações iniciais

/**
 * Calcula o tempo da primeira chegada de cada evento e coloca na heap de eventos
 * @param arrivals A heap de eventos
 */
void setup(priority_queue<Event> &arrivals) {
	arrivals.PUSH(Event)(DATA_ARRIVAL_TIME, EventType::DATA, new Packet(-1, DATA_TIME_OF_SERVICE))ENDPUSH;
	for (int i = 0; i < VOICE_CHANNELS; ++i) {
		arrivals.PUSH(Event)(VOICE_SILENCE_TIME, EventType::VOICE, new Packet(-1, i, VOICE_TIME_OF_SERVICE))ENDPUSH;
		channelsLastDeparture[i] = -1;
		activePeriodLength[i] = -1;
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
	if (server != nullptr && server->type == PackageType::DATA) {
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
			case PackageType::DATA:
				s.T1Acc += packet->totalTime;
				s.n1Packages++;
				s.X1Acc += packet->serviceTime + packet->property.wastedTime;
				break;
			case PackageType::VOICE:
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

	s.T2 = s.T2Acc / s.n2Packages;
	s.Nq2 = s.totalVoiceTimeInQueue == 0 ? 0 : (s.totalVoiceTimeInQueue / (s.totalTime * VOICE_CHANNELS));
	s.JitterMean = s.jitterAcc == 0 ? 0 : s.jitterAcc / s.n2Intervals;
	s.JitterVariance = (s.jitterAcc == 0 || s.jitterAccSqr == 0 || s.n2Intervals == 0) ? 0 :
					   (s.jitterAccSqr - s.jitterAcc * s.jitterAcc / s.n2Intervals) / (s.n2Intervals - 1);
}

// Outputs
/**
 * Imprime em stdout o valor das estatísticas globais: E[T1], E[W1], E[X1], E[Nq1], E[T2], E[W2], E[X2], E[Nq2], E[Δ] e V(Δ)
 */
void printStats(SimulationRound &s) {
#ifdef LOG
	cout << "Tempo médio entre chegadas: " << (aux_tempo_entre_chegadas / SAMPLES) << endl;
	cout << "lambda_1: " << n1Packages / max_time << " pacotes dados/seg" << endl;
	cout << "lambda_2: " << n2Packages / max_time << " pacotes voz/seg" << endl;
	cout << "Max Time: " << max_time << endl;
#endif
	if (JSON) {
		cout << R"({"t1":)" << s.T1 << R"(,"w1":)" << s.T1 - s.X1 << R"(,"x1":)" << s.X1 << R"(,"Nq1":)" << s.Nq1 <<
			 R"(,"t2":)" << s.T2 << R"(,"w2":)" << s.T2 - X2 << R"(,"x2":)" << X2 << R"(,"Nq2":)" << s.Nq2 <<
			 R"(,"jitterE":)" << s.JitterMean << R"(,"jitterV":)" << s.JitterVariance << R"(})";
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
	cout << "\t-s int\tNúmero de rodadas de simulação" << endl;
	cout << "\t-u double\tUtilização da fila 1" << endl;
	cout << "\t-p\tFila de dados pode ser interrompida ()" << endl;
}

void runSimulationRound(
		SimulationRound rounds[], int s, int samples,
		priority_queue<Event> &arrivals, queue<Packet *> &voice, queue<Packet *> &data,
		double &lastTime, int &interruptedDataPackages) {
	double t;
	for (int i = 0; i < samples; ++i) {
		#ifdef PROGRESS_BAR
		if (i % (samples / 20) == 0) cout << "Rodando simulação " << s << ": " << (i * 100 / samples) << "%" << '\r' << flush;
		#endif
		Event arrival = arrivals.top();
		arrivals.pop();
		registerAreaStatistics(voice.size(), data.size(), lastTime, arrival.time, rounds[s]);

		#ifdef LOG
		cout << "!! Time: " << arrival.time << (arrival.type == EventType::VOICE ? " (Voice)" : " (Server)") << endl;
			cout << "channelsLastDeparture[0]: " << channelsLastDeparture[0] << endl;
			cout << "Instante " << arrival.time << "s" << endl;
			cout << "Filas:" << endl << "> Dados: " << (data.size() == 0 ? 0 : data.size() - 1) << endl << "> Voz: " << (voice.size() == 0 ? 0 : voice.size() - 1) << endl << endl;
		#endif

		switch (arrival.type) {
			case EventType::DATA:
				// Configura o tempo de serviço do pacote e coloca o mesmo na fila
				arrival.packet->totalTime = -arrival.time;
				data.push(arrival.packet);

				// Coloca a próxima chegada de pacote de dados na fila de eventos
				arrivals.PUSH(Event)(arrival.time + DATA_ARRIVAL_TIME, EventType::DATA, new Packet(s, DATA_TIME_OF_SERVICE))ENDPUSH;

				if (server == nullptr) {
					server = serveEvent(arrival.packet, arrivals, arrival.time);
				}
				break;
			case EventType::VOICE:
				// Coloca a próxima chegada do canal na heap de eventos
				arrival.packet->totalTime = -arrival.time;
				t = arrival.time + VOICE_ARRIVAL_TIME;
				if (genEndOfActivePeriod()) {
					t += VOICE_SILENCE_TIME;
					arrival.packet->property.channel.lastVoicePackage = true;
				}

				// DEBUG Para trabalhar com tamanho fixo de pacotes no período ativo
//					if (activePeriodLength[0] < 4) {
//						activePeriodLength[0]++;
//						arrival.packet->property.channel.lastVoicePackage = false;
//					} else {
//						t += VOICE_SILENCE_TIME;
//						arrival.packet->property.channel.lastVoicePackage = true;
//						activePeriodLength[0] = 0;
//					}

				arrivals.PUSH(Event)(t, EventType::VOICE, new Packet(s, arrival.packet->property.channel.number, VOICE_TIME_OF_SERVICE))ENDPUSH;
				if (server == nullptr) {
					serveEvent(arrival.packet, arrivals, arrival.time);
					server = arrival.packet;
				} else if (PREEMPTION && server->type == PackageType::DATA) {
					interruptedDataPackages++;
					data.front()->property.wastedTime -= arrival.time;
					server = serveEvent(arrival.packet, arrivals, arrival.time);
				} else {
					voice.push(arrival.packet);
				}
				break;
			case EventType::SERVER:
				i--; //Saída do simulador não conta como amostra para a contagem
				SimulationRound &r = rounds[arrival.packet->simulation];
				switch (arrival.packet->type) {
					case PackageType::DATA:
						if (interruptedDataPackages) {
							interruptedDataPackages--;
							data.front()->property.wastedTime += arrival.time;
						} else {
							data.pop();
							arrival.packet->totalTime += arrival.time;
							countPacketIntoStatistics(arrival.packet, r);
							server = serveNextEvent(arrivals, voice, data, arrival.time);
						}
						break;
					case PackageType::VOICE:
						incrementJitter(r, arrival.packet->property.channel.number, arrival.time);
						if (!arrival.packet->property.channel.lastVoicePackage) {
							channelsLastDeparture[arrival.packet->property.channel.number] = arrival.time;
						} else { ;
							channelsLastDeparture[arrival.packet->property.channel.number] = -1;
						}
						arrival.packet->totalTime += arrival.time;
						countPacketIntoStatistics(arrival.packet, r);
						server = serveNextEvent(arrivals, voice, data, arrival.time);
						break;
				}
				break;
		}
		lastTime = arrival.time;
	}
}

void printConfidenceInterval(SimulationRound rounds[]) {
	SimulationRound &f = rounds[0];
	double W1 = f.T1 - f.X1,
			W2 = f.T2 - X2,
			VarianceT1 = 0,
			VarianceW1 = 0,
			VarianceX1 = 0,
			VarianceNq1 = 0,
			VarianceT2 = 0,
			VarianceW2 = 0,
			VarianceNq2 = 0,
			VarianceJitterMean = 0,
			VarianceJitterVariance = 0;
	for (int i = 1; i <= SIMULATIONS; ++i) {
		SimulationRound &s = rounds[i];
		VarianceT1 += (s.T1 - f.T1) * (s.T1 - f.T1);
		VarianceW1 += (s.T1 - s.X1 - W1) * (s.T1 - s.X1 - W1);
		VarianceX1 += (s.X1 - f.X1) * (s.X1 - f.X1);
		VarianceNq1 += (s.Nq1 - f.Nq1) * (s.Nq1 - f.Nq1);
		VarianceT2 += (s.T2 - f.T2) * (s.T2 - f.T2);
		VarianceW2 += (s.T2 - X2 - W2) * (s.T2 - X2 - W2);
		VarianceNq2 += (s.Nq2 - f.Nq2) * (s.Nq2 - f.Nq2);
		VarianceJitterMean += (s.JitterMean - f.JitterMean) * (s.JitterMean - f.JitterMean);
		VarianceJitterVariance += (s.JitterVariance - f.JitterVariance) * (s.JitterVariance - f.JitterVariance);
	}
	VarianceT1 /= (SIMULATIONS - 1);
	VarianceW1 /= (SIMULATIONS - 1);
	VarianceX1 /= (SIMULATIONS - 1);
	VarianceNq1 /= (SIMULATIONS - 1);
	VarianceT2 /= (SIMULATIONS - 1);
	VarianceW2 /= (SIMULATIONS - 1);
	VarianceNq2 /= (SIMULATIONS - 1);
	VarianceJitterMean /= (SIMULATIONS - 1);
	VarianceJitterVariance /= (SIMULATIONS - 1);

	VarianceT1 = sqrt(VarianceT1);
	VarianceW1 = sqrt(VarianceW1);
	VarianceX1 = sqrt(VarianceX1);
	VarianceNq1 = sqrt(VarianceNq1);
	VarianceT2 = sqrt(VarianceT2);
	VarianceW2 = sqrt(VarianceW2);
	VarianceNq2 = sqrt(VarianceNq2);
	VarianceJitterMean = sqrt(VarianceJitterMean);
	VarianceJitterVariance = sqrt(VarianceJitterVariance);

	VarianceT1 *= percentil90 / sqrt(SIMULATIONS);
	VarianceW1 *= percentil90 / sqrt(SIMULATIONS);
	VarianceX1 *= percentil90 / sqrt(SIMULATIONS);
	VarianceNq1 *= percentil90 / sqrt(SIMULATIONS);
	VarianceT2 *= percentil90 / sqrt(SIMULATIONS);
	VarianceW2 *= percentil90 / sqrt(SIMULATIONS);
	VarianceNq2 *= percentil90 / sqrt(SIMULATIONS);
	VarianceJitterMean *= percentil90 / sqrt(SIMULATIONS);
	VarianceJitterVariance *= percentil90 / sqrt(SIMULATIONS);
	if (JSON) {
		cout << ",\"confidenceInterval\":{" <<
			R"("T1":{"lesser":)" << f.T1 - VarianceT1 << ",\"mean\":" << f.T1 << ",\"upper\":" << f.T1 + VarianceT1 << ",\"percentage\":" << VarianceT1 * 200 / f.T1 << "}," <<
			R"("W1":{"lesser":)" << W1 - VarianceW1 << ",\"mean\":" << W1 << ",\"upper\":" << W1 + VarianceW1 << ",\"percentage\":" << VarianceW1 * 200 / W1 << "}," <<
			R"("X1":{"lesser":)" << f.X1 - VarianceX1 << ",\"mean\":" << f.X1 << ",\"upper\":" << f.X1 + VarianceX1 << ",\"percentage\":" << VarianceX1 * 200 / f.X1 << "}," <<
			R"("Nq1":{"lesser":)" << f.Nq1 - VarianceNq1 << ",\"mean\":" << f.Nq1 << ",\"upper\":" << f.Nq1 + VarianceNq1 << ",\"percentage\":" << VarianceNq1 * 200 / f.Nq1 << "}," <<
			R"("T2":{"lesser":)" << f.T2 - VarianceT2 << ",\"mean\":" << f.T2 << ",\"upper\":" << f.T2 + VarianceT2 << ",\"percentage\":" << VarianceT2 * 200 / f.T2 << "}," <<
			R"("W2":{"lesser":)" << W2 - VarianceW2 << ",\"mean\":" << W2 << ",\"upper\":" << W2 + VarianceW2 << ",\"percentage\":" << VarianceW2 * 200 / W2 << "}," <<
			R"("Nq2":{"lesser":)" << f.Nq2 - VarianceNq2 << ",\"mean\":" << f.Nq2 << ",\"upper\":" << f.Nq2 + VarianceNq2 << ",\"percentage\":" << VarianceNq2 * 200 / f.Nq2 << "}," <<
			R"("JitterMean":{"lesser":)" << f.JitterMean - VarianceJitterMean << ",\"mean\":" << f.JitterMean << ",\"upper\":" << f.JitterMean + VarianceJitterMean <<
																									",\"percentage\":" << VarianceJitterMean * 200 / f.JitterMean << "}," <<
			R"("JitterVariance":{"lesser":)" << f.JitterVariance - VarianceJitterVariance << ",\"mean\":" << f.JitterVariance << ",\"upper\":" << f.JitterVariance + VarianceJitterVariance <<
																									",\"percentage\":" << VarianceJitterVariance * 200 / f.JitterVariance << "}}";
	} else {
		cout << "T1:\t\t\t\t" << f.T1 - VarianceT1 << "\t" << f.T1 << "\t" << f.T1 + VarianceT1 << "\t" << VarianceT1 * 200 / f.T1 << "%" << endl;
		cout << "W1:\t\t\t\t" << W1 - VarianceW1 << "\t" << W1 << "\t" << W1 + VarianceW1 << "\t" << VarianceW1 * 200 / W1 << "%" << endl;
		cout << "X1:\t\t\t\t" << f.X1 - VarianceX1 << "\t" << f.X1 << "\t" << f.X1 + VarianceX1 << "\t" << VarianceX1 * 200 / f.X1 << "%" << endl;
		cout << "Nq1:\t\t\t" << f.Nq1 - VarianceNq1 << "\t" << f.Nq1 << "\t" << f.Nq1 + VarianceNq1 << "\t" << VarianceNq1 * 200 / f.Nq1 << "%" << endl;
		cout << "T2:\t\t\t\t" << f.T2 - VarianceT2 << "\t" << f.T2 << "\t" << f.T2 + VarianceT2 << "\t" << VarianceT2 * 200 / f.T2 << "%" << endl;
		cout << "W2:\t\t\t\t" << W2 - VarianceW2 << "\t" << W2 << "\t" << W2 + VarianceW2 << "\t" << VarianceW2 * 200 / W2 << "%" << endl;
		cout << "Nq2:\t\t\t" << f.Nq2 - VarianceNq2 << "\t" << f.Nq2 << "\t" << f.Nq2 + VarianceNq2 << "\t" << VarianceNq2 * 200 / f.Nq2 << "%" << endl;
		cout << "JitterMean:\t\t" << f.JitterMean - VarianceJitterMean << "\t" << f.JitterMean << "\t" << f.JitterMean + VarianceJitterMean << "\t" << VarianceJitterMean * 200 / f.JitterMean << "%" << endl;
		cout << "JitterVariance:\t" << f.JitterVariance - VarianceJitterVariance << "\t" << f.JitterVariance << "\t" << f.JitterVariance + VarianceJitterVariance << "\t" << VarianceJitterVariance * 200 / f.JitterVariance << "%" << endl;
	}
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
			case 'p': // Com ou sem preempção
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
			case 'j':
				JSON = true;
				break;
			case 'h':
				printHelp();
				return 0;
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

	// Cálculo e output das estatísticas de cada rodada
	SimulationRound &f = rounds[0];
	setupSimulationStats(f);
	if (JSON) cout << R"({"simulations":[)";
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
		if (JSON) {
			if (i > 1) cout << ",";
		} else {
			cout << "Simulação " << i << endl;
		}
		printStats(s);
	}
	if (JSON) cout << "],";

	// Cálculo e output do intervalo de confiança
	f.T1 /= SIMULATIONS;
	f.X1 /= SIMULATIONS;
	f.Nq1 /= SIMULATIONS;
	f.T2 /= SIMULATIONS;
	f.Nq2 /= SIMULATIONS;
	f.JitterMean /= SIMULATIONS;
	f.JitterVariance /= SIMULATIONS;
	if (JSON) {
		cout << R"("average":)";
	} else {
		cout << "=================================================" << endl;
		cout << "=================================================" << endl;
		cout << "=================================================" << endl;
		cout << "Resultado médio final" << endl;
	}
	printStats(f);
	printConfidenceInterval(rounds);
	if (JSON) cout << "}" << endl;
	return 0;
}
