#include <iostream>
#include <queue>
#include <random>
#include <functional>
#include "Event.h"
//#include <vld.h>

using namespace std;

constexpr double SERVER_SPEED = 2e6;//(1 >> 21); //bits/segundo
enum SimulationEvent {
	DATA_ENTER_QUEUE, DATA_ENTER_SERVER, DATA_LEAVE_SERVER, DATA_INTERRUPTED, VOICE_ENTER_QUEUE, VOICE_ENTER_SERVER, VOICE_LEAVE_SERVER
};

const int DATA_ARRIVAL_TIME = 500; // TODO V.A. exp. lambda1 para tempo até a próxima chegada de dados em milisegundos
constexpr int DATA_PACKAGE_SIZE_IN_BYTES = 512;// TODO V.A. com densidade f(x)=p1*u0(x-64)+p2*u0(x-512)+p3*u0(x-1500)+(p/1436)[u-1(x-64)–u-1(x-1500)] com p1=30%, p2=10%, p3=30%, p=1-p1-p2-p3=30%
constexpr int DATA_PACKAGE_SIZE_IN_BITS = DATA_PACKAGE_SIZE_IN_BYTES * 8;
constexpr double DATA_TIME_OF_SERVICE = DATA_PACKAGE_SIZE_IN_BITS / SERVER_SPEED;

const int VOICE_CHANNELS = 30;
const int VOICE_SILENCE_TIME = 650; // TODO V.A. exp. com média 650 para fim do período de silêncio do canal de voz em segundos
const int VOICE_ARRIVAL_TIME = 16; // Tempo até o próximo pacote de voz durante o período ativo em segundos
constexpr int VOICE_PACKAGE_SIZE_IN_BITS = 512; //bits
constexpr double VOICE_TIME_OF_SERVICE = VOICE_PACKAGE_SIZE_IN_BITS / SERVER_SPEED; // Tempo de transmissão do pacote de voz em segundos
const int MEAN_N_VOICE_PACKAGE = 22;
const double MEAN_SILENCE_PERIOD_DURATION = .65; // In seconds

// Parâmetros da simulação
const int AMOSTRAS = 1;
const bool PREEMPION = false;

// Códigos de erro
const int INVALID_SERVICE_TYPE = 1;
const int INVALID_EVENT_ARRIVAL = 2;

// Data channel random variable generators
// TODO refatorar para uma função normal
auto genDataPackageSize = []()->int {
    static auto engine = default_random_engine{1};
    static uniform_real_distribution dist {0.0, 1.0};
    double x = dist(engine);

    if(x < .3)
        return 64;
    else if(x <.4)
        return 512;
    else if(x < .7)
        return 1500;
    else
        // TODO corrigir probabilidade das pontas
        return (x - .7)*1436/.3 + 64;
};

auto genDataServiceTime = []() { return genDataPackageSize()*8/SERVER_SPEED; };
// Voice channel random variable generators
auto genEndOfActivePeriod = bind(bernoulli_distribution{1.0/MEAN_N_VOICE_PACKAGE}, default_random_engine{2}); // Returns true when the voice channel enters a silence period after a voice package arrival
auto genSilencePeriod = bind(exponential_distribution{1.0/MEAN_SILENCE_PERIOD_DURATION}, default_random_engine{3});
/**
 * Coloca o evento na fila de eventos esperados com o tempo adequado
 * @param event
 * @param events_queue
 * @param currentTime
 */
void serveEvent(Event event, priority_queue<Event> &events_queue, double currentTime) {
	events_queue.emplace(currentTime + event.serviceTime, EventType::SERVER);
}

/**
 * Serve o próximo evento em espera
 * @param events_queue
 * @param voice
 * @param data
 * @param currentTime
 * @return O tipo de evento colocado no servidor (EMPTY caso nenhum evento tenha sido colocado)
 */
EventType serveNextEvent(priority_queue<Event> &events_queue, queue<Event> &voice, queue<Event> &data, double currentTime) {
	// TODO separar os tempos para coletar as estatísticas
	// TODO pegar tempo de serviço associado ao evento
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

int main(int argc, char *argv[]) {
	priority_queue<Event> arrivals; // Estrutura para organizar as chegadas
	queue<Event> data, voice; // Filas de data e voz
	EventType serverOccupied = EventType::EMPTY; // Tipo de evento presente no servidor
	int interruptedDataPackages = 0; // Pacotes de dados interrompidos (usado para invalidar os respectivos pacotes quando a antiga saída aparecer na lista)

	// Calcula o tempo de chegada de cada evento
	arrivals.emplace(DATA_ARRIVAL_TIME, DATA_TIME_OF_SERVICE, EventType::DATA);
	for (int i = 0; i < VOICE_CHANNELS; ++i) {
		arrivals.emplace(VOICE_SILENCE_TIME, VOICE_TIME_OF_SERVICE, EventType::VOICE, i);
	}
/*	Event server(INTMAX_MAX, SERVER);
	for (int i = 0; i < AMOSTRAS; ++i) {
		Event arrival = arrivals.top();
		arrivals.pop();
		while (server < arrival) {
			// Freguês sai do servidor e puxa o próximo da fila
			if (!voice.empty()) {
				double start = server.time;
				server = voice.front();
				voice.pop();
				server.time = start + VOICE_TIME_OF_SERVICE; //TODO separar os tempos para coletar as estatísticas
			} else if (!data.empty()) {
				double start = server.time;
				server = data.front();
				data.pop();
				server.time = start + DATA_TIME_OF_SERVICE; //TODO separar os tempos para coletar as estatísticas
			} else {
				server = Event(INTMAX_MAX, SERVER);
			}
		}

		// Adiciona a nova chegada na fila
		switch (arrival.type) {
			case DATA:
				data.push(arrival);
				arrivals.push(Event(arrival.time + DATA_ARRIVAL_TIME, DATA));
				break;
			case VOICE:
				voice.push(arrival);
				if (voiceContinue)
					arrivals.push(Event(arrival.time + VOICE_ARRIVAL_TIME, VOICE));
				else
					arrivals.push(Event(arrival.time + VOICE_ARRIVAL_TIME + VOICE_SILENCE_TIME, VOICE));
				break;
			default:
				break;
		}
	}*/

	for (int i = 0; i < AMOSTRAS; ++i) {
		Event arrival = arrivals.top();
		arrivals.pop();
		double t;
		switch (arrival.type) {
			case EventType::DATA:
				data.push(arrival);
				arrivals.emplace(arrival.time + DATA_ARRIVAL_TIME, DATA_TIME_OF_SERVICE, EventType::DATA);

				if (serverOccupied == EventType::EMPTY) {
					serveEvent(arrival, arrivals, arrival.time);
					serverOccupied = EventType::DATA;
				}
				break;
			case EventType::VOICE:
				arrival.serviceTime = VOICE_TIME_OF_SERVICE;
//				voice.push(arrival);
				t = arrival.time + VOICE_ARRIVAL_TIME;
				if (genEndOfActivePeriod())
					t += VOICE_SILENCE_TIME;
				arrivals.emplace(t, VOICE_TIME_OF_SERVICE, EventType::VOICE, arrival.channel);

				if (serverOccupied == EventType::EMPTY) {
					serveEvent(arrival, arrivals, arrival.time);
					serverOccupied = EventType::VOICE;
				} else if (PREEMPION && serverOccupied == EventType::DATA) {
					interruptedDataPackages++;
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
				serverOccupied = serveNextEvent(arrivals, voice, data, arrival.time);
				break;
			default:
				return INVALID_EVENT_ARRIVAL;
		}
	}

	return 0;
}
