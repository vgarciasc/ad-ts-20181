#include <iostream>
#include <queue>
#include <random>
#include <functional>
#include "Event.h"

using namespace std;

//PARÂMETROS
// Parâmetros da simulação
const int AMOSTRAS = 1;
const bool PREEMPION = false;
constexpr double SERVER_SPEED = 2e6;//(1 >> 21); //bits/segundo
enum SimulationEvent {
	DATA_ENTER_QUEUE, DATA_ENTER_SERVER, DATA_LEAVE_SERVER, DATA_INTERRUPTED, VOICE_ENTER_QUEUE, VOICE_ENTER_SERVER, VOICE_LEAVE_SERVER
};

//Data//
const int DATA_ARRIVAL_TIME = 500; // TODO V.A. exp. lambda1 para tempo até a próxima chegada de dados em milisegundos
constexpr int DATA_PACKAGE_SIZE_IN_BYTES = 512;// TODO V.A. com densidade f(x)=p1*u0(x-64)+p2*u0(x-512)+p3*u0(x-1500)+(p/1436)[u-1(x-64)–u-1(x-1500)] com p1=30%, p2=10%, p3=30%, p=1-p1-p2-p3=30%
constexpr int DATA_PACKAGE_SIZE_IN_BITS = DATA_PACKAGE_SIZE_IN_BYTES * 8;
constexpr double DATA_TIME_OF_SERVICE = DATA_PACKAGE_SIZE_IN_BITS / SERVER_SPEED;

//Voice//
const int VOICE_CHANNELS = 30;
const int VOICE_SILENCE_TIME = 650; // TODO V.A. exp. com média 650 para fim do período de silêncio do canal de voz em segundos
const int VOICE_ARRIVAL_TIME = 16; // Tempo até o próximo pacote de voz durante o período ativo em segundos
constexpr int VOICE_PACKAGE_SIZE_IN_BITS = 512; //bits
constexpr double VOICE_TIME_OF_SERVICE = VOICE_PACKAGE_SIZE_IN_BITS / SERVER_SPEED; // Tempo de transmissão do pacote de voz em segundos
const int MEAN_N_VOICE_PACKAGE = 22;
const double MEAN_SILENCE_PERIOD_DURATION = .65; // In seconds

// Códigos de erro
const int INVALID_SERVICE_TYPE = 1;
const int INVALID_EVENT_ARRIVAL = 2;

// Geradores de V.A.
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
// Voice channel random variable generators
auto genEndOfActivePeriod = bind(bernoulli_distribution{1.0 / MEAN_N_VOICE_PACKAGE}, default_random_engine{2}); // Returns true when the voice channel enters a silence period after a voice package arrival
auto genSilencePeriod = bind(exponential_distribution{1.0 / MEAN_SILENCE_PERIOD_DURATION}, default_random_engine{3});

// Estatísticas
double W1, X1, Nq1, W2, X2, Nq2, JitterMean, JitterVariance;

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
	// TODO confirmar
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
		return EventType::VOICE;
	} else if (!data.empty()) {
		serveEvent(data.front(), events_queue, currentTime);
		return EventType::DATA;
	} else {
		return EventType::EMPTY;
	}
}

/**
 * Registra os dados passados nas respectivas estatísticas
 * @param event
 * @param time
 * @param voiceChannel
 */
void registerStatistics(SimulationEvent event, double time, int voiceChannel = 0) {
	switch (event) {
		case DATA_ENTER_QUEUE:
			break;
		case DATA_ENTER_SERVER:
			break;
		case DATA_LEAVE_SERVER:
			break;
		case DATA_INTERRUPTED:
			break;
		case VOICE_ENTER_QUEUE:
			break;
		case VOICE_ENTER_SERVER:
			break;
		case VOICE_LEAVE_SERVER:
			break;
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
	EventType serverOccupied = EventType::EMPTY; // Tipo de evento presente no servidor
	int interruptedDataPackages = 0; // Pacotes de dados interrompidos (usado para invalidar os respectivos pacotes quando a antiga saída aparecer na lista)

	setup(arrivals);

	for (int i = 0; i < AMOSTRAS; ++i) {
		Event arrival = arrivals.top();
		arrivals.pop();
		double t;
		switch (arrival.type) {
			case EventType::DATA:
				// Configura o tempo de serviço do pacote e coloca o mesmo na fila
				arrival.stats->serviceTime = genDataServiceTime();
				arrival.stats->enterQueueTime = arrival.time;
				data.push(arrival);

				arrivals.emplace(arrival.time + DATA_ARRIVAL_TIME, EventType::DATA);

				if (serverOccupied == EventType::EMPTY) {
					serveEvent(arrival, arrivals, arrival.time);
					serverOccupied = EventType::DATA;
				}
				break;
			case EventType::VOICE:
				// Configura o tempo de entrada na fila
				arrival.stats->enterQueueTime = arrival.time;

				// Coloca a próxima chegada do canal na heap de eventos
				t = arrival.time + VOICE_ARRIVAL_TIME;
				if (genEndOfActivePeriod()) {
					t += VOICE_SILENCE_TIME;
					//TODO finalizar jitter
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
					// TODO registrar estatísticas
					case EventType::VOICE :
//						voice.pop();
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

	return 0;
}
