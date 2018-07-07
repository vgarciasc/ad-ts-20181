#include <iostream>
#include <queue>
#include <random>
#include <functional>
#include "Event.h"

using namespace std;

//PARÂMETROS
// Parâmetros da simulação
const int SAMPLES = 1;
const int SIMULATIONS = 1;
const bool PREEMPION = true;
const double UTILIZATION_1 = 0.1; //ρ1
constexpr double SERVER_SPEED = 2e6; //2Mb/segundo
enum SimulationEvent {
	DATA_ENTER_QUEUE, DATA_ENTER_SERVER, DATA_LEAVE_SERVER, DATA_INTERRUPTED, VOICE_ENTER_QUEUE, VOICE_ENTER_SERVER, VOICE_LEAVE_SERVER
};
EventType serverOccupied = EventType::EMPTY; // Tipo de evento presente no servidor

//Data//
// Data channel random variable generators
int genDataPackageSize() {
	static auto engine = default_random_engine{1};
	static uniform_real_distribution dist{0.0, 1.0};
	double x = dist(engine);

	if (x < .3)
		return 64;
	else if (x < .4)
		return 512;
	else if (x < .7)
		return 1500;
	else
		// TODO corrigir probabilidade das pontas
		return (x - .7) * 1436 / .3 + 64;
}

auto genDataServiceTime = []() { return genDataPackageSize() * 8 / SERVER_SPEED; };
#define DATA_TIME_OF_SERVICE genDataServiceTime()

const double DATA_ARRIVAL_RATE = UTILIZATION_1 / (755 * 8 / SERVER_SPEED); // λ1 = ρ1/E[X1] = ρ1/(E[L]bytes*8/(2Mb/s))
auto genDataArrivalTime = bind(exponential_distribution{DATA_ARRIVAL_RATE}, default_random_engine{4});
#define DATA_ARRIVAL_TIME genDataArrivalTime()

//Voice//
const int VOICE_CHANNELS = 30;
const double VOICE_ARRIVAL_TIME = .016; // Tempo até o próximo pacote de voz durante o período ativo em segundos
const int VOICE_PACKAGE_SIZE_IN_BITS = 512; //bits
constexpr double VOICE_TIME_OF_SERVICE = VOICE_PACKAGE_SIZE_IN_BITS / SERVER_SPEED; // Tempo de transmissão do pacote de voz em segundos
const int MEAN_N_VOICE_PACKAGE = 22;
const double MEAN_SILENCE_PERIOD_DURATION = .65; // In seconds

// Voice channel random variable generators
auto genEndOfActivePeriod = bind(bernoulli_distribution{1.0 / MEAN_N_VOICE_PACKAGE}, default_random_engine{2}); // Returns true when the voice channel enters a silence period after a voice package arrival
auto genSilencePeriod = bind(exponential_distribution{1.0 / MEAN_SILENCE_PERIOD_DURATION}, default_random_engine{3});
#define VOICE_SILENCE_TIME genSilencePeriod()

// Códigos de erro
const int INVALID_SERVICE_TYPE = 1;
const int INVALID_EVENT_ARRIVAL = 2;

// Estatísticas
const double X2 = VOICE_TIME_OF_SERVICE;
double W1, X1, Nq1, W2, Nq2, JitterMean, JitterVariance;
double channelsLastDeparture[VOICE_CHANNELS];

//FUNÇÕES
/**
 * Calcula o tempo da primeira chegada de cada evento e coloca na heap de eventos
 * @param arrivals A heap de eventos
 */
void setup(priority_queue<Event> &arrivals) {
	arrivals.emplace(DATA_ARRIVAL_TIME, EventType::DATA);
	for (int i = 0; i < VOICE_CHANNELS; ++i) {
		arrivals.emplace(VOICE_SILENCE_TIME, EventType::VOICE, new Stats(i, VOICE_TIME_OF_SERVICE));
	}
}

/**
 * Coloca o evento na fila de eventos esperados com o tempo adequado e dá o valor para a estatísca de tempo de espera
 * @param event Evento a ser servido
 * @param events_queue A heap de eventos
 * @param currentTime Momento em que o evento foi servido
 */
void serveEvent(Event event, priority_queue<Event> &events_queue, double currentTime) {
	event.stats->waitTime += currentTime - event.stats->enterQueueTime;
	events_queue.emplace(currentTime + event.stats->serviceTime, EventType::SERVER, event.stats);
}

/**
 * Serve o próximo evento em espera
 * @param events_queue Heap de eventos
 * @param voice Fila de eventos de voz (maior prioridade)
 * @param data Fila de eventos de dados (menor prioridade)
 * @param currentTime Momento em que o evento será servido
 * @return O tipo de evento colocado no servidor (EMPTY caso nenhum evento tenha sido colocado)
 */
EventType serveNextEvent(priority_queue<Event> &events_queue, queue<Event> &voice, queue<Event> &data, double currentTime) {
	if (!voice.empty()) {
		serveEvent(voice.front(), events_queue, currentTime);
		voice.pop();
		return EventType::VOICE;
	} else if (!data.empty()) {
		serveEvent(data.front(), events_queue, currentTime);
		return EventType::DATA;
	} else {
		return EventType::EMPTY;
	}
}

int n1Packages, n2Packages, n2Intervals;
double totalDataTime, totalVoiceTime, totalX1, totalTime, jitterAcc, jitterAccSqr;

void registerAreaStatistics(unsigned long Nq2, unsigned long Nq1, double lastTime, double currentTime) {
	// Add nq1, nq2, w1, w2 and x1
	double t = currentTime - lastTime;
	totalDataTime += Nq1 * t;
	totalVoiceTime += Nq2 * t;
	if (serverOccupied == EventType::DATA) {
		totalX1 += t;
	}
}

void calculateAreaStatistics() {
	// Estatísticas dos dados
	Nq1 = totalDataTime / totalTime;
	W1 = totalDataTime / n1Packages;
	X1 = totalX1 / n1Packages;
	// Estatísticas dos pacotes de voz
	Nq2 = totalVoiceTime / totalTime;
	W2 = totalVoiceTime / n2Packages;
	JitterMean = jitterAcc / n2Intervals;
	JitterVariance = (jitterAccSqr * n2Intervals - jitterAcc * jitterAcc) / (n2Intervals * (n2Intervals - 1));
	//X2 constante
}

void incrementJitter(int channel, double currentTime) {
	double lastDeparture = channelsLastDeparture[channel];
	double lastInterval = lastDeparture - currentTime;

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
	serverOccupied = EventType::EMPTY;
	JitterMean = JitterVariance = n1Packages = n2Packages = n2Intervals = jitterAcc = jitterAccSqr = 0;
	totalTime = totalDataTime = totalVoiceTime = 0.0;

	for (int i = 0; i < VOICE_CHANNELS; i++) {
		channelsLastDeparture[i] = -1;
	}
}

/**
 * Imprime em stdout o valor das estatísticas globais: E[T1], E[W1], E[X1], E[Nq1], E[T2], E[W2], E[X2], E[Nq2], E[Δ] e V(Δ)
 */
void printStats() {
	cout << "E[T1]: " << W1 + X1 << ", E[W1]: " << W1 << ", E[X1]: " << X1 << ", E[Nq1]: " << Nq1 << ", E[T2]: " << W2 + X2 << ", E[W2]: " << W2 <<
		 ", E[X2]: " << X2 << ", E[Nq2]: " << Nq2 << ", E[Δ]; " << JitterMean << ", V(Δ):" << JitterVariance << endl;
}

int main(int argc, char *argv[]) {
	// Filas e variáveis de controle
	priority_queue<Event> arrivals; // Estrutura para organizar as chegadas
	queue<Event> data, voice; // Filas de data e voz
	int interruptedDataPackages = 0; // Pacotes de dados interrompidos (usado para invalidar os respectivos pacotes quando a antiga saída aparecer na lista)

	setup(arrivals);

	for (int j = 0; j < SIMULATIONS; ++j) {
		resetStats();
		double lastTime = 0;
		for (int i = 0; i < SAMPLES; ++i) {
			Event arrival = arrivals.top();
			arrivals.pop();
			double t;
			registerAreaStatistics(voice.size(), serverOccupied == EventType::DATA ? data.size() - 1 : data.size(), lastTime, arrival.time);
			switch (arrival.type) {
				case EventType::DATA:
					// Configura o tempo de serviço do pacote e coloca o mesmo na fila
					arrival.stats->serviceTime = DATA_TIME_OF_SERVICE;
					arrival.stats->enterQueueTime = arrival.time;
					data.push(arrival);

					// Coloca a próxima chegada de pacote de dados na fila de eventos
					arrivals.emplace(arrival.time + DATA_ARRIVAL_TIME, EventType::DATA);

					if (serverOccupied == EventType::EMPTY) {
						serveEvent(arrival, arrivals, arrival.time);
						serverOccupied = EventType::DATA;
					}
					break;
				case EventType::VOICE:
					// Coloca a próxima chegada do canal na heap de eventos
					t = arrival.time + VOICE_ARRIVAL_TIME;
					if (genEndOfActivePeriod()) {
						t += VOICE_SILENCE_TIME;
						arrival.stats->lastVoicePackage = true;
					}
					arrivals.emplace(t, EventType::VOICE, new Stats(arrival.stats->channel, VOICE_TIME_OF_SERVICE));

					if (serverOccupied == EventType::EMPTY) {
						serveEvent(arrival, arrivals, arrival.time);
						serverOccupied = EventType::VOICE;
					} else if (PREEMPION && serverOccupied == EventType::DATA) {
						interruptedDataPackages++;
						data.front().stats->enterQueueTime = arrival.time; // Necessário para contar o tempo que o pacote passou na fila
						serveEvent(arrival, arrivals, arrival.time);
						serverOccupied = EventType::VOICE;
					} else {
						voice.push(arrival);
					}
					break;
				case EventType::SERVER:
					switch (serverOccupied) {
						case EventType::VOICE :
							incrementJitter(arrival.stats->channel, arrival.time);
							if (!arrival.stats->lastVoicePackage) {
								channelsLastDeparture[channel] = arrival.time;
							} else {
								channelsLastDeparture[channel] = -1;
							}
							break;
						case EventType::DATA :
							if (interruptedDataPackages) interruptedDataPackages--;
							else data.pop();
							break;
						default:
							return INVALID_SERVICE_TYPE; //ERRO!!!
					}
					delete arrival.stats; // Controle de memória
					serverOccupied = serveNextEvent(arrivals, voice, data, arrival.time);
					break;
				default:
					return INVALID_EVENT_ARRIVAL;
			}
		}
		calculateAreaStatistics();
		printStats();
	}

	return 0;
}
