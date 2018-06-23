#include <iostream>
#include <queue>
#include "Event.h"
//#include <vld.h>

using namespace std;

constexpr unsigned int SERVER_SPEED = (1 >> 21)/1000; //Bytes/milisegundo

const int TD = 500; // TODO V.A. exp. lambda1 para tempo até a próxima chegada de dados em milisegundos
constexpr int DATA_PACKAGE_SIZE = 512;// TODO V.A. com densidade f(x) = p1*u0(x-64) + p2*u0(x-512) + p3*u0(x-1500) + (p/1436)[u-1(x-64) – u-1(x-1500)] com p1=30%, p2=10%, p3 = 30%, p = 1 - p1 - p2 - p3 = 30%.
constexpr long long int DATA_TIME_OF_SERVICE = DATA_PACKAGE_SIZE/SERVER_SPEED;

const int TVS = 650; // TODO V.A. exp. com média 650 para fim do período de silêncio do canal de voz em milisegundos
const int TVA = 16; // Tempo até o próximo pacote de voz durante o período ativo em milisegundos
constexpr int VOICE_PACKAGE_SIZE = 64; //bytes
constexpr long long int VOICE_TIME_OF_SERVICE = VOICE_PACKAGE_SIZE/SERVER_SPEED; // Tempo de transmissão do pacote de voz em milisegundos
const double voiceContinue = 1.0 / 22; // Probabilidade de continuar o período ativo do canal de voz

const int AMOSTRAS = 1000;

int main(int argc, char *argv[]) {
	priority_queue<Event> arrivals; // Estrutura para organizar as chegadas
	queue<Event> data, voice; // Filas de data e voz

	// Calcula o tempo de chegada de cada evento
	arrivals.push(Event(TD, DATA));
	for (int i = 0; i < 30; ++i) {
		arrivals.push(Event(TVS, VOICE));
	}
	Event server(INTMAX_MAX, SERVER);
	for (int i = 0; i < AMOSTRAS; ++i) {
		Event arrival = arrivals.top();
		arrivals.pop();
		while (server < arrival) {
			// Freguês sai do servidor e puxa o próximo da fila
			if (!voice.empty()) {
				long long int start = server.time;
				server = voice.front();
				voice.pop();
				server.time = start + VOICE_TIME_OF_SERVICE; //TODO separar os tempos para coletar as estatísticas
			} else if (!data.empty()) {
				long long int start = server.time;
				server = data.front();
				data.pop();
				server.time = start + DATA_TIME_OF_SERVICE; //TODO separar os tempos para coletar as estatísticas
			} else {
				server = NULL;
			}
		}

		// Adiciona a nova chegada na fila
		switch (arrival.type) {
			case DATA:
				data.push(arrival);
				arrivals.push(Event(arrival.time + TD, DATA));
				break;
			case VOICE:
				voice.push(arrival);
				if (voiceContinue)
					arrivals.push(Event(arrival.time + TVA, VOICE));
				else
					arrivals.push(Event(arrival.time + TVA + TVS, VOICE));
				break;
			default:
				break;
		}
	}
	return 0;
}