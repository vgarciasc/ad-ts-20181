#include <iostream>
#include <queue>
#include <random>
#include <functional>
#include <math.h>
#include <time.h>
#include "Event.h"

//#define PROGRESS_BAR
//#define LOG

using namespace std;

//PARÂMETROS
// Parâmetros da simulação
int SAMPLES = 100000;
int SIMULATIONS = 5;
bool PREEMPTION = false;
double UTILIZATION_1 = 0.1; //ρ1
constexpr double SERVER_SPEED = 2e6; //2Mb/segundo
enum SimulationEvent {
	DATA_ENTER_QUEUE, DATA_ENTER_SERVER, DATA_LEAVE_SERVER, DATA_INTERRUPTED, VOICE_ENTER_QUEUE, VOICE_ENTER_SERVER, VOICE_LEAVE_SERVER
};
Packet *server;
//EventType serverOccupied = EventType::EMPTY; // Tipo de evento presente no servidor
int TRANSIENT_SAMPLE_NUMBER = 0;

double genRandUnitary() {
    return ((double) rand() / RAND_MAX);
}

//Data//
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
	 	return genRandUnitary() * 1436 + 64;
}

auto genDataServiceTime = []() {
	return genDataPackageSize() * 8 / SERVER_SPEED;
};
#define DATA_TIME_OF_SERVICE genDataServiceTime()

double DATA_ARRIVAL_RATE = UTILIZATION_1 / (755 * 8 / SERVER_SPEED); // λ1 = ρ1/E[X1] = ρ1/(E[L]bytes*8/(2Mb/s))
double genDataArrivalTime(){
    return - log(genRandUnitary()) / (DATA_ARRIVAL_RATE);
}
#define DATA_ARRIVAL_TIME genDataArrivalTime()

//Voice//
const int VOICE_CHANNELS = 2;
const double VOICE_ARRIVAL_TIME = .016; // Tempo até o próximo pacote de voz durante o período ativo em segundos
const int VOICE_PACKAGE_SIZE_IN_BITS = 512; //bits
constexpr double VOICE_TIME_OF_SERVICE = VOICE_PACKAGE_SIZE_IN_BITS / SERVER_SPEED; // Tempo de transmissão do pacote de voz em segundos
const int MEAN_N_VOICE_PACKAGE = 22;
const double MEAN_SILENCE_PERIOD_DURATION = .65; // In seconds

// Voice channel random variable generators
auto genEndOfActivePeriod = []() {
    return genRandUnitary() < (1.0 / MEAN_N_VOICE_PACKAGE);
};
auto genSilencePeriod = []() {
    return - log(genRandUnitary()) / (1.0 / MEAN_SILENCE_PERIOD_DURATION);
};
#define VOICE_SILENCE_TIME genSilencePeriod()

// Códigos de erro
const int INVALID_SERVICE_TYPE = 1;
const int INVALID_EVENT_ARRIVAL = 2;

// Estatísticas
const double X2 = VOICE_TIME_OF_SERVICE;
double W1, X1, Nq1, W2, Nq2, JitterMean, JitterVariance;
double channelsLastDeparture[VOICE_CHANNELS];
int activePeriodLength[VOICE_CHANNELS];

double aux_tempo_entre_chegadas = 0;

double max_time = 0;

struct SimulationRound {
	int n1Packages, n2Packages, n2Intervals;
	double startTime;
	double totalDataTime, totalVoiceTime, totalX1, totalTime, jitterAcc, jitterAccSqr;
	double W1, X1, Nq1, W2, Nq2, JitterMean, JitterVariance;
	double channelsLastDeparture[VOICE_CHANNELS];
};

//#define EMPLACE
#ifdef EMPLACE
#define PUSH(type) emplace(
#define ENDPUSH
#else
#define PUSH(type) push(type
#define ENDPUSH )
#endif

//FUNÇÕES
/**
 * Calcula o tempo da primeira chegada de cada evento e coloca na heap de eventos
 * @param arrivals A heap de eventos
 */
void setup(priority_queue<Event> &arrivals) {
	arrivals.PUSH(Event)(DATA_ARRIVAL_TIME, EventType::DATA, new Packet(-1, DATA_TIME_OF_SERVICE))ENDPUSH;
	for (int i = 0; i < VOICE_CHANNELS; ++i) {
		arrivals.PUSH(Event)(VOICE_SILENCE_TIME, EventType::VOICE, new Packet(-1, i, VOICE_TIME_OF_SERVICE))ENDPUSH;
	}
	for (double &i : channelsLastDeparture) {
		i = -1;
	}
	for (int &i : activePeriodLength) {
		i = 0;
	}
}

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

int n1Packages, n2Packages, n2Intervals;
double totalDataTime, totalVoiceTime, totalX1, totalTime, jitterAcc, jitterAccSqr;

double voiceSilenceTotalTime = 0;
int voiceSilenceTotalTimeGenerated = 0;

void registerAreaStatistics(unsigned long Nq2, unsigned long Nq1, double lastTime, double currentTime) {
	// Add nq1, nq2, w1, w2 and x1
	double t = currentTime - lastTime;
	if (server != nullptr && server->type == PackageType::DATA) {
	    totalX1 += t;
		Nq1--;
	}
	totalDataTime += Nq1 * t;
	totalVoiceTime += Nq2 * t;
}

void calculateAreaStatistics(SimulationRound &s) {
	double roundDuration = totalTime - s.startTime;

    // Estatísticas dos dados
	cout << "Total Data Time: " << totalDataTime << ", Round Duration: " << roundDuration << endl;
	s.Nq1 = totalDataTime == 0 ? 0 : totalDataTime / roundDuration;
	s.W1 = totalDataTime == 0 ? 0 : totalDataTime / n1Packages;
//	cout << "Total W1: " << totalDataTime << ", n1Packages: " << n1Packages << endl;
	cout << "Total X1: " << totalX1 << ", n1Packages: " << n1Packages << endl;

	s.X1 = totalX1 == 0 ? 0 : totalX1 / n1Packages;
	// Estatísticas dos pacotes de voz

	cout << "Total Voice Time: " << totalVoiceTime << ", Total Time: " << totalTime << endl;
	cout << "n2Packages: " << n2Packages << endl;
	s.Nq2 = totalVoiceTime == 0 ? 0 : (totalVoiceTime / (roundDuration * VOICE_CHANNELS));
	s.W2 = totalVoiceTime == 0 ? 0 : totalVoiceTime / n2Packages;
	//X2 constante
	cout << "Jitter Acc: " << jitterAcc << ", n2Intervals: " << n2Intervals << endl;
	s.JitterMean = jitterAcc == 0 ? 0 : jitterAcc / n2Intervals;
	s.JitterVariance = (jitterAcc == 0 || jitterAccSqr == 0 || n2Intervals == 0) ? 0 :
					   (jitterAccSqr * n2Intervals - jitterAcc * jitterAcc) / (n2Intervals * (n2Intervals - 1));
}

void incrementJitter(int channel, double currentTime) {
	double lastDeparture = channelsLastDeparture[channel];
	double lastInterval = currentTime - lastDeparture;

	if (lastDeparture != -1) {
		jitterAcc += lastInterval;
		jitterAccSqr += lastInterval * lastInterval;
		n2Intervals += 1;
	}
}

/**
 * Coloca o valor de todas as variáveis usadas para as estatísticas em zero para evitar
 */
void resetStats() {
	JitterMean = JitterVariance = n1Packages = n2Packages = n2Intervals = 0;
	jitterAcc = jitterAccSqr = 0;
	totalTime = totalDataTime = totalVoiceTime = totalX1 = 0.0;
}

void setupSimulationStats(SimulationRound &s) {
	s.JitterMean = s.JitterVariance = s.n1Packages = s.n2Packages = s.n2Intervals = 0;
	s.jitterAcc = s.jitterAccSqr = 0;
	s.totalTime = s.totalDataTime = s.totalVoiceTime = 0.0;
	s.startTime = 0.0;
}

/**
 * Imprime em stdout o valor das estatísticas globais: E[T1], E[W1], E[X1], E[Nq1], E[T2], E[W2], E[X2], E[Nq2], E[Δ] e V(Δ)
 */
void printStats(SimulationRound &s) {
#ifdef LOG
	cout << "Tempo médio entre chegadas: " << (aux_tempo_entre_chegadas / SAMPLES) << endl;
	cout << "Voice Silence Media: " << (voiceSilenceTotalTime / voiceSilenceTotalTimeGenerated) << endl;
	cout << "lambda_1: " << n1Packages / max_time << " pacotes dados/seg" << endl;
	cout << "lambda_2: " << n2Packages / max_time << " pacotes voz/seg" << endl;
	cout << "Max Time: " << max_time << endl;
#endif
    cout << "lambda_1: " << n1Packages / max_time << " pacotes dados/seg" << endl;
    cout << "lambda_2: " << n2Packages / max_time << " pacotes voz/seg" << endl;
	cout << "E[T1]: " << s.W1 + s.X1 << ", E[W1]: " << s.W1 << ", E[X1]: " << s.X1 << ", E[Nq1]: " << s.Nq1 << endl;
	cout << "E[T2]: " << s.W2 + X2 << ", E[W2]: " << s.W2 <<
		 ", E[X2]: " << X2 << ", E[Nq2]: " << s.Nq2 << ", E[Δ]; " << s.JitterMean << ", V(Δ):" << s.JitterVariance << endl;
}

/**
-  * Imprime texto de ajuda para mostrar as opções a serem passadas para o programa
-  */
void printHelp() {
	cout << "Utilização:" << endl;
	cout << "_exe_ [OPTIONS]" << endl;
	cout << "\t-h\t\tMostra essa ajuda" << endl;
	cout << "\t-a int\tNúmero de amostras" << endl;
	cout << "\t-s int\tNúmero de rodadas de simulação" << endl;
	cout << "\t-u double\tUtilização da fila 1" << endl;
	cout << "\t-p\tFila de dados pode ser interrompida ()" << endl;
}

/**
 * Insere as estatísticas do pacote na rodada de simulação equivalente e deleta o objeto para evitar vazamento de memória
 * @param packet
 */
// TODO
void countPacketIntoStatistics(Packet *packet, SimulationRound rounds[]) {
	if (packet->simulation > -1) {
		SimulationRound s = rounds[packet->simulation];
		switch (packet->type) {
			case PackageType::DATA:
				break;
			case PackageType::VOICE:
				break;
		}
	}
	delete packet;
}

int main(int argc, char *argv[]) {
    srand(time(0));

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
    		case 'u': // Utilização da fila de dados
    			UTILIZATION_1 = stod(argv[++p]);
    			DATA_ARRIVAL_RATE = UTILIZATION_1 / (755 * 8 / SERVER_SPEED); // Î»1 = Ï1/E[X1] = Ï1/(E[L]bytes*8/(2Mb/s))
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
	SimulationRound rounds[SIMULATIONS]; // Estatísticas das várias rodadas de simulação.
    double lastTime, lastArrivalTime = 0.0;

	setup(arrivals);

	for (int i = 0; i < SIMULATIONS; ++i) {
		setupSimulationStats(rounds[i]);
	}

	for (int i = 0; i < TRANSIENT_SAMPLE_NUMBER; ++i) {
		// Executa o período transiente da simulação
	}

	for (int s = 0; s < SIMULATIONS; ++s) {
		double t;
        rounds[s].startTime = totalTime;

        resetStats();

		for (int i = 0; i < SAMPLES; ++i) {
#ifdef PROGRESS_BAR
			if (i % (SAMPLES / 20) == 0) cout << (i * 100 / SAMPLES) << "%" << 'r' << flush;
#endif
//		getchar();
			Event arrival = arrivals.top();
			arrivals.pop();

			registerAreaStatistics(voice.size(), data.size(), lastTime, arrival.time);
			max_time = arrival.time;
			lastTime = arrival.time;
			totalTime = arrival.time;

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

					aux_tempo_entre_chegadas += arrival.time - lastArrivalTime;
					lastArrivalTime = arrival.time;
#ifdef LOG
				cout << "Tempo entre chegadas: " << aux_tempo_entre_chegadas << " ~ " << arrival.time - lastTime << endl;
				cout << "> Gerado novo pacote de dados chegando em: " << arrival.time  << " + DATA_ARRIVAL_TIME" << endl;
#endif
					// Coloca a próxima chegada de pacote de dados na fila de eventos
					arrivals.PUSH(Event)(arrival.time + DATA_ARRIVAL_TIME, EventType::DATA, new Packet(s, DATA_TIME_OF_SERVICE))ENDPUSH;

					if (server == nullptr) {
						server = serveEvent(arrival.packet, arrivals, arrival.time);
					}
					break;
				case EventType::VOICE:
//					cout << "Chegada de Pacote de Voz (" << (activePeriodLength[arrival.stats->channel] + 1) << "º do periodo ativo)" << endl;
					// Coloca a próxima chegada do canal na heap de eventos
					t = arrival.time + VOICE_ARRIVAL_TIME;
//					if (genEndOfActivePeriod()) {
						t += VOICE_SILENCE_TIME;
						arrival.packet->property.channel.lastVoicePackage = true;
//					}

					// DEBUG Para trabalhar com tamanho fixo de pacotes no período ativo
//					if (activePeriodLength[0] < 4) {
//						activePeriodLength[0]++;
//						arrival.packet->property.channel.lastVoicePackage = false;
//					} else {
//						t += VOICE_SILENCE_TIME;
//						arrival.packet->property.channel.lastVoicePackage = true;
//						activePeriodLength[0] = 0;
//					}

					n2Packages++;
					arrivals.PUSH(Event)(t, EventType::VOICE, new Packet(s, arrival.packet->property.channel.number, VOICE_TIME_OF_SERVICE))ENDPUSH;
					if (server == nullptr) {
						serveEvent(arrival.packet, arrivals, arrival.time);
						server = arrival.packet;
					} else if (PREEMPTION && server->type == PackageType::DATA) {
						interruptedDataPackages++;
						data.front()->totalTime = -arrival.time;
						server = serveEvent(arrival.packet, arrivals, arrival.time);
					} else {
						voice.push(arrival.packet);
					}
					break;
				case EventType::SERVER:
					i--; //Saída do simulador não conta como amostra para a contagem
					switch (arrival.packet->type) {
						case PackageType::DATA:
							if (interruptedDataPackages) {
								interruptedDataPackages--;
							} else {
								n1Packages += 1;
								data.pop();
								arrival.packet->totalTime += arrival.time;
								countPacketIntoStatistics(arrival.packet, rounds);
								server = serveNextEvent(arrivals, voice, data, arrival.time);
							}
							break;
						case PackageType::VOICE:
							incrementJitter(arrival.packet->property.channel.number, arrival.time);
							if (!arrival.packet->property.channel.lastVoicePackage) {
								channelsLastDeparture[arrival.packet->property.channel.number] = arrival.time;
							} else {
;								channelsLastDeparture[arrival.packet->property.channel.number] = -1;
							}
							arrival.packet->totalTime += arrival.time;
							countPacketIntoStatistics(arrival.packet, rounds);
							server = serveNextEvent(arrivals, voice, data, arrival.time);
							break;
					}
					break;
			}
		}
        calculateAreaStatistics(rounds[s]);
	}

	for (int i = 0; i < SIMULATIONS; ++i) {
		SimulationRound s = rounds[i];
		printStats(s);
	}

	return 0;
}
