#include <iostream>
#include <queue>
#include "Event.h"
//#include <vld.h>

using namespace std;

constexpr double SERVER_SPEED = 2e6;//(1 >> 21); //bits/segundo

const int DATA_ARRIVAL_TIME = 500; // TODO V.A. exp. lambda1 para tempo até a próxima chegada de dados em milisegundos
constexpr int DATA_PACKAGE_SIZE_IN_BYTES = 512;// TODO V.A. com densidade f(x)=p1*u0(x-64)+p2*u0(x-512)+p3*u0(x-1500)+(p/1436)[u-1(x-64)–u-1(x-1500)] com p1=30%, p2=10%, p3=30%, p=1-p1-p2-p3=30%
constexpr int DATA_PACKAGE_SIZE_IN_BITS = DATA_PACKAGE_SIZE_IN_BYTES * 8;
constexpr double DATA_TIME_OF_SERVICE = DATA_PACKAGE_SIZE_IN_BITS / SERVER_SPEED;

const int VOICE_SILENCE_TIME = 650; // TODO V.A. exp. com média 650 para fim do período de silêncio do canal de voz em segundos
const int VOICE_ARRIVAL_TIME = 16; // Tempo até o próximo pacote de voz durante o período ativo em segundos
constexpr int VOICE_PACKAGE_SIZE_IN_BITS = 512; //bits
constexpr double VOICE_TIME_OF_SERVICE = VOICE_PACKAGE_SIZE_IN_BITS / SERVER_SPEED; // Tempo de transmissão do pacote de voz em segundos
const double voiceContinue = 1.0 / 22; // Probabilidade de continuar o período ativo do canal de voz

// Parâmetros da simulação
const int AMOSTRAS = 1;
const bool PREEMPION = false;

// Códigos de erro
const int INVALID_SERVICE_TYPE = 1;
const int INVALID_EVENT_ARRIVAL = 2;

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

int main(int argc, char *argv[]) {
	priority_queue<Event> arrivals; // Estrutura para organizar as chegadas
	queue<Event> data, voice; // Filas de data e voz
	EventType serverOccupied = EventType::EMPTY; // Tipo de evento presente no servidor
	int interruptedDataPackages = 0; // Pacotes de dados interrompidos (usado para invalidar os respectivos pacotes quando a antiga saída aparecer na lista)

	// Calcula o tempo de chegada de cada evento
	arrivals.emplace(DATA_ARRIVAL_TIME, EventType::DATA);
	for (int i = 0; i < 30; ++i) {
		arrivals.emplace(VOICE_SILENCE_TIME, EventType::VOICE);
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
		switch (arrival.type) {
			case EventType::DATA:
				arrival.serviceTime = DATA_TIME_OF_SERVICE;
				data.push(arrival);
				arrivals.emplace(arrival.time + DATA_ARRIVAL_TIME, EventType::DATA);

				if (serverOccupied == EventType::EMPTY) {
					serveEvent(arrival, arrivals, arrival.time);
					serverOccupied = EventType::DATA;
				}
				break;
			case EventType::VOICE:
				arrival.serviceTime = VOICE_TIME_OF_SERVICE;
//				voice.push(arrival);
				if (voiceContinue)
					arrivals.emplace(arrival.time + VOICE_ARRIVAL_TIME, EventType::VOICE);
				else
					arrivals.emplace(arrival.time + VOICE_ARRIVAL_TIME + VOICE_SILENCE_TIME, EventType::VOICE);

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
